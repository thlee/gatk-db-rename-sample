#!/usr/bin/env python3
"""
gatk_db_rename_sample.py - Rename samples in GATK GenomicsDB

GATK GenomicsDB does not natively support renaming samples after import.
This tool directly modifies the callset.json metadata files inside
GenomicsDB workspaces to rename samples without rebuilding the database.

Usage:
    # List samples
    gatk_db_rename_sample.py list /path/to/db
    gatk_db_rename_sample.py list /path/to/db --count

    # Rename a single sample
    gatk_db_rename_sample.py rename /path/to/db --old SAMPLE_A --new SAMPLE_B

    # Batch rename from a mapping file (TSV: old_name<TAB>new_name)
    gatk_db_rename_sample.py rename /path/to/db --map rename_map.tsv

    # Preview changes without modifying
    gatk_db_rename_sample.py rename /path/to/db --map rename_map.tsv --dry-run

    # Validate consistency across interval workspaces
    gatk_db_rename_sample.py validate /path/to/db

    # Restore from backup after a mistake
    gatk_db_rename_sample.py restore /path/to/db

How it works:
    GenomicsDB stores sample metadata in callset.json within each workspace.
    Each entry contains sample_name, row_idx, idx_in_file, and stream_name.
    This tool modifies only the sample_name field while keeping row_idx
    intact, preserving consistency with the underlying TileDB array.

    The given path is searched recursively for GenomicsDB workspaces
    (directories containing callset.json). This means ALL workspaces found
    under the path will be affected — including unrelated databases if
    they happen to be in subdirectories. Always use --dry-run first to
    verify which workspaces will be modified.

Tips:
    To effectively REMOVE samples from GenomicsDB output, rename them with
    a distinctive prefix (e.g. '__REMOVED__'), then filter them out when
    generating VCFs:

        # 1. Rename samples to mark for removal
        gatk_db_rename_sample.py rename /path/to/db \\
            --map remove_map.tsv   # e.g. "SAMPLE_A  __REMOVED__SAMPLE_A"

        # 2. After GenotypeGVCFs, filter out marked samples
        bcftools view -S ^exclude_list.txt input.vcf.gz -Oz -o output.vcf.gz

        # exclude_list.txt contains one __REMOVED__* name per line

    Entries cannot be deleted from callset.json because TileDB still holds
    data at those row indices; missing entries cause GATK to fail.

Disclaimer:
    USE AT YOUR OWN RISK. This tool modifies internal metadata files of
    GATK GenomicsDB. It is not endorsed by or affiliated with the Broad
    Institute or the GATK development team. The authors are not responsible
    for any data loss or corruption. Always use --dry-run first and verify
    results after applying modifications.

Author: TH Lee (concept, direction, review), Claude/Anthropic (code)
License: MIT
"""

import argparse
import json
import shutil
import sys
from collections import Counter
from datetime import datetime
from glob import glob
from pathlib import Path


__version__ = "0.1.0"


# ============================================================
# Core functions
# ============================================================

def find_workspaces(db_path):
    """Find all GenomicsDB workspace directories under the given path.

    A workspace is identified by the presence of callset.json.
    Supports both a single workspace and a parent directory containing
    multiple interval workspaces.
    """
    db_path = Path(db_path)

    if (db_path / "callset.json").exists():
        return [db_path]

    workspaces = sorted(p.parent for p in db_path.rglob("callset.json"))

    if not workspaces:
        print(f"Error: No GenomicsDB workspaces found under {db_path}", file=sys.stderr)
        print("  A workspace must contain a callset.json file.", file=sys.stderr)
        sys.exit(1)

    return workspaces


def load_callset(workspace):
    """Load callset.json from a workspace."""
    path = Path(workspace) / "callset.json"
    try:
        with open(path) as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        print(f"Error: malformed callset.json in {workspace}: {e}", file=sys.stderr)
        sys.exit(1)
    if "callsets" not in data:
        print(f"Error: callset.json in {workspace} missing 'callsets' key", file=sys.stderr)
        sys.exit(1)
    return data


def save_callset(workspace, data, backup=True):
    """Save callset.json, optionally creating a timestamped backup."""
    path = Path(workspace) / "callset.json"

    if backup:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        shutil.copy2(path, path.with_suffix(f".json.bak.{timestamp}"))

    with open(path, "w") as f:
        json.dump(data, f, indent=2)
        f.write("\n")


def get_sample_names(data):
    """Extract sample names from callset data."""
    return [entry["sample_name"] for entry in data["callsets"]]


def load_rename_map(map_file):
    """Load rename mapping from a TSV file (old_name<TAB>new_name).

    Lines starting with '#' are comments. Blank lines are skipped.
    Accepts both tab and space as delimiters.
    """
    rename_map = {}
    with open(map_file) as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            if len(parts) != 2:
                print(f"Warning: skipping malformed line {lineno}: {line!r}",
                      file=sys.stderr)
                continue
            old_name, new_name = parts
            if old_name in rename_map:
                print(f"Warning: duplicate mapping for '{old_name}' at line {lineno}",
                      file=sys.stderr)
            rename_map[old_name] = new_name
    return rename_map


# ============================================================
# Validation
# ============================================================

def validate_rename(current_names, rename_map):
    """Validate rename mapping. Returns (warnings, errors)."""
    warnings = []
    errors = []

    current_set = set(current_names)

    # Mappings not found in DB
    missing = set(rename_map.keys()) - current_set
    if missing:
        examples = ", ".join(sorted(missing)[:5])
        suffix = f" ... and {len(missing)-5} more" if len(missing) > 5 else ""
        warnings.append(f"{len(missing)} name(s) in mapping not found in DB: {examples}{suffix}")

    # Check for duplicates after rename
    new_names = [rename_map.get(n, n) for n in current_names]
    dups = {k: v for k, v in Counter(new_names).items() if v > 1}
    if dups:
        examples = ", ".join(f"'{k}' (x{v})" for k, v in sorted(dups.items())[:5])
        suffix = f" ... and {len(dups)-5} more" if len(dups) > 5 else ""
        errors.append(f"Rename would create {len(dups)} duplicate name(s): {examples}{suffix}")

    # New names colliding with existing unrenamed samples
    staying = current_set - set(rename_map.keys())
    collisions = staying & set(rename_map.values())
    if collisions:
        errors.append(
            f"New names collide with existing unrenamed samples: "
            f"{', '.join(sorted(collisions)[:5])}"
        )

    return warnings, errors


# ============================================================
# Subcommands
# ============================================================

def cmd_list(args):
    """List samples in a GenomicsDB."""
    workspaces = find_workspaces(args.db)
    data = load_callset(workspaces[0])
    names = get_sample_names(data)

    if args.count:
        print(len(names))
        return

    for name in sorted(names) if args.sort else names:
        print(name)

    if args.verbose:
        print(f"\nTotal: {len(names)} samples in {len(workspaces)} workspace(s)",
              file=sys.stderr)


def cmd_rename(args):
    """Rename samples in a GenomicsDB."""
    if args.map:
        rename_map = load_rename_map(args.map)
    elif args.old and args.new:
        rename_map = {args.old: args.new}
    elif args.old:
        print("Error: --old requires --new", file=sys.stderr)
        sys.exit(1)
    elif args.new:
        print("Error: --new requires --old", file=sys.stderr)
        sys.exit(1)
    else:
        print("Error: provide either --map FILE or --old NAME --new NAME",
              file=sys.stderr)
        sys.exit(1)

    if not rename_map:
        print("Error: no valid rename mappings found", file=sys.stderr)
        sys.exit(1)

    workspaces = find_workspaces(args.db)

    # First pass: scan all workspaces
    ws_info = []
    for ws in workspaces:
        data = load_callset(ws)
        names = get_sample_names(data)
        applicable = sum(1 for n in names if n in rename_map)
        if applicable > 0:
            warnings, errors = validate_rename(names, rename_map)
        else:
            warnings, errors = [], []
        ws_info.append((ws, data, len(names), applicable, warnings, errors))

    # Display summary table
    name_w = max((len(ws.name) for ws, *_ in ws_info), default=9)
    name_w = max(name_w, 9)

    print(f"Found {len(workspaces)} workspace(s), {len(rename_map)} mapping(s)\n")
    print(f"  {'#':>4}  {'Workspace':<{name_w}}  {'Samples':>7}  {'Renames':>7}  Status")
    print(f"  {'':->4}  {'':->{name_w}}  {'':->7}  {'':->7}  ------")

    for i, (ws, data, n_samples, applicable, warnings, errors) in enumerate(ws_info):
        if applicable == 0:
            status = "-"
        elif errors:
            status = f"{len(errors)} error(s)"
        elif warnings:
            status = f"{len(warnings)} warning(s)"
        else:
            status = "OK"
        print(f"  {i+1:>4}  {ws.name:<{name_w}}  {n_samples:>7}  {applicable:>7}  {status}")

    # Show unique warnings/errors
    all_warnings = sorted(set(w for *_, warnings, _ in ws_info for w in warnings))
    all_errors = sorted(set(e for *_, errors in ws_info for e in errors))
    if all_warnings or all_errors:
        print()
        for w in all_warnings:
            print(f"  Warning: {w}")
        for e in all_errors:
            print(f"  ERROR: {e}")

    # Dry-run: show rename preview from first applicable workspace
    if args.dry_run:
        for ws, data, _, applicable, _, _ in ws_info:
            if applicable > 0:
                print(f"\n[DRY-RUN] Renames that would be applied:")
                for entry in data["callsets"]:
                    old = entry["sample_name"]
                    if old in rename_map:
                        print(f"  {old}  ->  {rename_map[old]}")
                break
        print(f"\nNo changes made (dry-run mode).")
        return

    # Second pass: per-workspace confirmation and rename
    modified = 0
    skipped = 0
    for i, (ws, data, _, applicable, _, errors) in enumerate(ws_info):
        if applicable == 0:
            continue

        if errors and not args.force:
            skipped += 1
            continue

        if not args.yes:
            try:
                answer = input(f"\n  [{i+1}] {ws.name} — rename {applicable} sample(s)? [y/N] ")
            except (EOFError, KeyboardInterrupt):
                print("\nAborted.")
                sys.exit(1)
            if answer.lower() not in ("y", "yes"):
                skipped += 1
                continue

        for entry in data["callsets"]:
            if entry["sample_name"] in rename_map:
                entry["sample_name"] = rename_map[entry["sample_name"]]
        save_callset(ws, data, backup=not args.no_backup)
        modified += 1
        if args.verbose:
            print(f"  [{i+1}] {ws.name}: {applicable} renamed")

    if modified:
        print(f"\nDone. Modified {modified} workspace(s).")
        if skipped:
            print(f"  Skipped: {skipped}")
        if not args.no_backup:
            print(f"Backups saved as callset.json.bak.*")
    else:
        print(f"\nNo workspaces were modified.")


def cmd_validate(args):
    """Validate GenomicsDB sample consistency across all interval workspaces."""
    workspaces = find_workspaces(args.db)
    print(f"Validating {len(workspaces)} workspace(s)...")

    issues = []
    reference_names = None
    reference_ws = None

    reference_mapping = None

    for ws in workspaces:
        data = load_callset(ws)
        names = get_sample_names(data)

        # Duplicates within workspace
        dups = {k: v for k, v in Counter(names).items() if v > 1}
        if dups:
            issues.append(f"{ws.name}: duplicate sample names: {dups}")

        # Duplicate row_idx within workspace
        row_indices = [entry["row_idx"] for entry in data["callsets"]]
        dup_rows = {k: v for k, v in Counter(row_indices).items() if v > 1}
        if dup_rows:
            issues.append(f"{ws.name}: duplicate row_idx values: {dup_rows}")

        # Consistency across workspaces
        name_tuple = tuple(sorted(names))
        name_to_row = {entry["sample_name"]: entry["row_idx"] for entry in data["callsets"]}

        if reference_names is None:
            reference_names = name_tuple
            reference_ws = ws.name
            reference_mapping = name_to_row
        else:
            if name_tuple != reference_names:
                ref_set, cur_set = set(reference_names), set(names)
                issues.append(
                    f"{ws.name}: differs from {reference_ws} "
                    f"(only in ref: {len(ref_set - cur_set)}, only here: {len(cur_set - ref_set)})"
                )

            # row_idx consistency across workspaces
            mismatched = [
                n for n in set(name_to_row) & set(reference_mapping)
                if name_to_row[n] != reference_mapping[n]
            ]
            if mismatched:
                examples = ", ".join(sorted(mismatched)[:5])
                suffix = f" ... and {len(mismatched)-5} more" if len(mismatched) > 5 else ""
                issues.append(
                    f"{ws.name}: row_idx mismatch vs {reference_ws} "
                    f"for {len(mismatched)} sample(s): {examples}{suffix}"
                )

    print(f"Samples: {len(reference_names) if reference_names else 0}")

    if issues:
        print(f"\nFound {len(issues)} issue(s):")
        for issue in issues:
            print(f"  - {issue}")
        sys.exit(1)
    else:
        print("All checks passed.")


def cmd_restore(args):
    """Restore callset.json from the most recent backup."""
    workspaces = find_workspaces(args.db)
    restored = 0

    for ws in workspaces:
        backups = sorted(glob(str(ws / "callset.json.bak.*")))
        if not backups:
            continue
        shutil.copy2(backups[-1], ws / "callset.json")
        restored += 1
        if args.verbose:
            print(f"  {ws.name}: restored from {Path(backups[-1]).name}")

    total = len(workspaces)
    if restored:
        print(f"Restored {restored} workspace(s) from backup.")
        if restored < total:
            print(f"  Warning: {total - restored} workspace(s) had no backup and were not restored.",
                  file=sys.stderr)
    else:
        print("No backups found.")


# ============================================================
# CLI
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        prog="gatk_db_rename_sample",
        description="Rename samples in GATK GenomicsDB without rebuilding.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  %(prog)s list /path/to/db
  %(prog)s list /path/to/db --count
  %(prog)s rename /path/to/db --old OLD_NAME --new NEW_NAME
  %(prog)s rename /path/to/db --map rename_map.tsv
  %(prog)s rename /path/to/db --map rename_map.tsv --dry-run
  %(prog)s validate /path/to/db
  %(prog)s restore /path/to/db

tips:
  To effectively REMOVE samples, rename them with a prefix like
  '__REMOVED__', then filter them out of VCF output using bcftools:

    bcftools view -S ^exclude_list.txt input.vcf.gz -Oz -o output.vcf.gz
        """,
    )
    parser.add_argument("--version", action="version",
                        version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command", required=True)

    # --- list ---
    p_list = subparsers.add_parser("list", help="List samples in the database")
    p_list.add_argument("db", help="Path to search for GenomicsDB workspaces (recursive)")
    p_list.add_argument("--count", action="store_true", help="Print only the count")
    p_list.add_argument("--sort", action="store_true", help="Sort alphabetically")
    p_list.add_argument("-v", "--verbose", action="store_true")

    # --- rename ---
    p_rename = subparsers.add_parser("rename", help="Rename samples")
    p_rename.add_argument("db", help="Path to search for GenomicsDB workspaces (recursive)")
    group = p_rename.add_mutually_exclusive_group()
    group.add_argument("--map", metavar="FILE",
                       help="TSV mapping file (old_name<TAB>new_name)")
    group.add_argument("--old", metavar="NAME", help="Old sample name (with --new)")
    p_rename.add_argument("--new", metavar="NAME", help="New sample name (with --old)")
    p_rename.add_argument("--dry-run", action="store_true",
                          help="Preview changes without modifying")
    p_rename.add_argument("--no-backup", action="store_true",
                          help="Skip backup creation")
    p_rename.add_argument("--force", action="store_true",
                          help="Proceed despite validation errors")
    p_rename.add_argument("-y", "--yes", action="store_true",
                          help="Skip confirmation prompt")
    p_rename.add_argument("-v", "--verbose", action="store_true")

    # --- validate ---
    p_val = subparsers.add_parser("validate", help="Validate sample consistency")
    p_val.add_argument("db", help="Path to search for GenomicsDB workspaces (recursive)")
    p_val.add_argument("-v", "--verbose", action="store_true")

    # --- restore ---
    p_res = subparsers.add_parser("restore", help="Restore callset.json from backup")
    p_res.add_argument("db", help="Path to search for GenomicsDB workspaces (recursive)")
    p_res.add_argument("-v", "--verbose", action="store_true")

    args = parser.parse_args()

    commands = {
        "list": cmd_list,
        "rename": cmd_rename,
        "validate": cmd_validate,
        "restore": cmd_restore,
    }
    commands[args.command](args)


if __name__ == "__main__":
    main()
