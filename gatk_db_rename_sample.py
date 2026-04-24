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
    GenomicsDB stores sample metadata in callset.json within each interval
    workspace. Each entry contains sample_name, row_idx, idx_in_file, and
    stream_name. This tool modifies only the sample_name field while keeping
    row_idx intact, preserving consistency with the underlying TileDB array.

    The db path can be a single workspace (containing callset.json) or a
    parent directory containing multiple interval workspaces.

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
    with open(Path(workspace) / "callset.json") as f:
        return json.load(f)


def save_callset(workspace, data, backup=True):
    """Save callset.json, optionally creating a timestamped backup."""
    path = Path(workspace) / "callset.json"

    if backup:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
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
    else:
        print("Error: provide either --map FILE or --old NAME --new NAME",
              file=sys.stderr)
        sys.exit(1)

    if not rename_map:
        print("Error: no valid rename mappings found", file=sys.stderr)
        sys.exit(1)

    workspaces = find_workspaces(args.db)
    data = load_callset(workspaces[0])
    current_names = get_sample_names(data)
    applicable = sum(1 for n in current_names if n in rename_map)
    warnings, errors = validate_rename(current_names, rename_map)

    print(f"Workspaces:  {len(workspaces)}")
    print(f"Samples:     {len(current_names)}")
    print(f"Mappings:    {len(rename_map)} (applicable: {applicable})")

    for w in warnings:
        print(f"  Warning: {w}")
    for e in errors:
        print(f"  ERROR: {e}")

    if errors and not args.force:
        print("\nAborted. Use --force to override errors.", file=sys.stderr)
        sys.exit(1)

    if args.dry_run:
        print(f"\n[DRY-RUN] Changes that would be applied:")
        for entry in data["callsets"]:
            old = entry["sample_name"]
            if old in rename_map:
                print(f"  {old}  ->  {rename_map[old]}")
        print(f"\nNo changes made (dry-run mode).")
        return

    print(f"\nApplying renames to {len(workspaces)} workspace(s)...")
    for i, ws in enumerate(workspaces):
        data = load_callset(ws)
        renamed = 0
        for entry in data["callsets"]:
            if entry["sample_name"] in rename_map:
                entry["sample_name"] = rename_map[entry["sample_name"]]
                renamed += 1
        save_callset(ws, data, backup=not args.no_backup)
        if args.verbose:
            print(f"  [{i+1}/{len(workspaces)}] {ws.name}: {renamed} renamed")

    print(f"Done. Renamed {applicable} sample(s) across {len(workspaces)} workspace(s).")
    if not args.no_backup:
        print(f"Backups saved as callset.json.bak.*")


def cmd_validate(args):
    """Validate GenomicsDB sample consistency across all interval workspaces."""
    workspaces = find_workspaces(args.db)
    print(f"Validating {len(workspaces)} workspace(s)...")

    issues = []
    reference_names = None
    reference_ws = None

    for ws in workspaces:
        data = load_callset(ws)
        names = get_sample_names(data)

        # Duplicates within workspace
        dups = {k: v for k, v in Counter(names).items() if v > 1}
        if dups:
            issues.append(f"{ws.name}: duplicate sample names: {dups}")

        # Consistency across workspaces
        name_tuple = tuple(sorted(names))
        if reference_names is None:
            reference_names = name_tuple
            reference_ws = ws.name
        elif name_tuple != reference_names:
            ref_set, cur_set = set(reference_names), set(names)
            issues.append(
                f"{ws.name}: differs from {reference_ws} "
                f"(only in ref: {len(ref_set - cur_set)}, only here: {len(cur_set - ref_set)})"
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

    if restored:
        print(f"Restored {restored} workspace(s) from backup.")
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
    p_list.add_argument("db", help="GenomicsDB path (workspace or parent dir)")
    p_list.add_argument("--count", action="store_true", help="Print only the count")
    p_list.add_argument("--sort", action="store_true", help="Sort alphabetically")
    p_list.add_argument("-v", "--verbose", action="store_true")

    # --- rename ---
    p_rename = subparsers.add_parser("rename", help="Rename samples")
    p_rename.add_argument("db", help="GenomicsDB path (workspace or parent dir)")
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
    p_rename.add_argument("-v", "--verbose", action="store_true")

    # --- validate ---
    p_val = subparsers.add_parser("validate", help="Validate sample consistency")
    p_val.add_argument("db", help="GenomicsDB path (workspace or parent dir)")
    p_val.add_argument("-v", "--verbose", action="store_true")

    # --- restore ---
    p_res = subparsers.add_parser("restore", help="Restore callset.json from backup")
    p_res.add_argument("db", help="GenomicsDB path (workspace or parent dir)")
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
