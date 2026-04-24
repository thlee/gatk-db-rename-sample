# gatk-db-rename-sample

Rename samples in [GATK GenomicsDB](https://gatk.broadinstitute.org/hc/en-us/articles/360035891051-GenomicsDB) without rebuilding the database.

## Problem

GATK GenomicsDB does not support renaming or removing samples after import.
The only official workaround is to rebuild the entire database from scratch
— which can take days for large cohorts.

## Solution

This tool directly modifies the `callset.json` metadata inside GenomicsDB
workspaces. It changes the `sample_name` field while keeping `row_idx`
intact, preserving consistency with the underlying TileDB array.

- **No database rebuild required**
- Handles single renames or batch renames from a mapping file
- Validates for duplicate names before applying changes
- Automatic backup with one-command restore
- Pure Python, no dependencies beyond the standard library

## Requirements

- Python 3.6+

## Quick Start

```bash
# List samples
./gatk_db_rename_sample.py list /path/to/db

# Rename a single sample
./gatk_db_rename_sample.py rename /path/to/db --old SAMPLE_A --new SAMPLE_B

# Batch rename from a mapping file
./gatk_db_rename_sample.py rename /path/to/db --map rename_map.tsv

# Preview changes first (recommended)
./gatk_db_rename_sample.py rename /path/to/db --map rename_map.tsv --dry-run
```

## Commands

### `list` — Show samples in the database

```bash
./gatk_db_rename_sample.py list /path/to/db           # full list
./gatk_db_rename_sample.py list /path/to/db --sort     # sorted
./gatk_db_rename_sample.py list /path/to/db --count    # count only
```

### `rename` — Rename samples

Single rename:

```bash
./gatk_db_rename_sample.py rename /path/to/db --old OLD_NAME --new NEW_NAME
```

Batch rename using a mapping file (TSV: `old_name<TAB>new_name`):

```
# rename_map.tsv
OLD_SAMPLE_1    NEW_SAMPLE_1
OLD_SAMPLE_2    NEW_SAMPLE_2
```

```bash
./gatk_db_rename_sample.py rename /path/to/db --map rename_map.tsv
```

Options:

| Flag | Description |
|------|-------------|
| `--dry-run` | Preview changes without modifying files |
| `--no-backup` | Skip automatic backup |
| `--force` | Proceed despite validation errors (e.g., duplicates) |
| `-v, --verbose` | Show per-workspace progress |

### `validate` — Check consistency across workspaces

Verifies that all interval workspaces have the same sample list and no
duplicate names:

```bash
./gatk_db_rename_sample.py validate /path/to/db
```

### `restore` — Undo changes

Restores `callset.json` from the most recent backup in each workspace:

```bash
./gatk_db_rename_sample.py restore /path/to/db
```

## DB Path

The db path can be either:

- A **single workspace** directory (containing `callset.json`)
- A **parent directory** containing multiple interval workspaces

For example, if your GenomicsDB was created with interval-based workspaces:

```
db/
├── my_db.Chr01-00000001-05000000/
│   ├── callset.json
│   ├── vidmap.json
│   └── ...
├── my_db.Chr01-05000001-10000000/
│   ├── callset.json
│   └── ...
└── ...
```

You can point the tool at the `db/` directory, and it will find and update
all workspaces automatically.

## How It Works

Each GenomicsDB interval workspace contains a `callset.json` file:

```json
{
  "callsets": [
    {"sample_name": "SAMPLE_A", "row_idx": 0, "idx_in_file": 0, "stream_name": "..."},
    {"sample_name": "SAMPLE_B", "row_idx": 1, "idx_in_file": 0, "stream_name": "..."}
  ]
}
```

This tool modifies only the `sample_name` field. The `row_idx` field maps to
the physical data location in the TileDB array and **must not be changed** —
altering it will cause GATK to fail.

## Tips

### Effectively removing samples

GenomicsDB entries cannot be deleted because TileDB retains data at those row
indices. Removing an entry from `callset.json` causes GATK to fail with:

```
GenomicsDB JNI Error: No sample/CallSet name specified for TileDB row N
```

Instead, rename the samples with a distinctive prefix, then filter them out
from the VCF output:

```bash
# 1. Create a mapping file to mark samples for removal
cat > remove_map.tsv << 'EOF'
SAMPLE_TO_REMOVE_1    __REMOVED__SAMPLE_TO_REMOVE_1
SAMPLE_TO_REMOVE_2    __REMOVED__SAMPLE_TO_REMOVE_2
EOF

# 2. Apply the rename
./gatk_db_rename_sample.py rename /path/to/db --map remove_map.tsv

# 3. Run GenotypeGVCFs as usual (output will contain __REMOVED__* samples)
gatk GenotypeGVCFs -V gendb:///path/to/workspace -R ref.fa -O output.vcf

# 4. Filter out the marked samples
cat > exclude.txt << 'EOF'
__REMOVED__SAMPLE_TO_REMOVE_1
__REMOVED__SAMPLE_TO_REMOVE_2
EOF
bcftools view -S ^exclude.txt output.vcf -Oz -o final.vcf.gz
```

## Backups

By default, a timestamped backup is created before each modification:

```
callset.json.bak.20260424_103500
```

The backup is stored alongside the original `callset.json` in each workspace
directory, so you can always find and inspect it next to the file it protects.

Use `restore` to revert all workspaces to their most recent backup, or
`--no-backup` to skip backup creation.

## Disclaimer

**USE AT YOUR OWN RISK.** This tool modifies internal metadata files of GATK
GenomicsDB. It is not endorsed by or affiliated with the Broad Institute or
the GATK development team.

While this tool has been tested on real-world datasets, there is no guarantee
that it will work correctly with all versions of GATK or GenomicsDB. The
authors are not responsible for any data loss or corruption resulting from
the use of this tool.

**Always use `--dry-run` first** to preview changes, and verify your results
after applying modifications. Backups are created automatically, but you
should maintain your own backups of critical data.

## Acknowledgments

This tool was conceived and directed by TH Lee after encountering sample
naming issues in a large-scale genotyping project. The code was
written by [Claude](https://claude.ai/) (Anthropic) under TH Lee's
direction and review.

## License

MIT
