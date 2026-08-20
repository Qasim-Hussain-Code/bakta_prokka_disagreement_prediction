#!/usr/bin/env bash
# Download the Bakta database and record exactly which release was used.
#
# The database release determines which calls Bakta makes. It is part of the
# result, not metadata, so its version, DOI and on-disk checksum are recorded
# here and quoted in results/metrics/02_annotation_versions.json.
#
# TYPE IS 'light', NOT 'full'. This machine has ~13 GB free and the full
# database does not fit. This is a result-affecting constraint, recorded as
# such rather than buried in an install note.

set -euo pipefail
source "$(dirname "$0")/../scripts/lib_env.sh"
cd "$BPDP_ROOT"

DB_TYPE="light"
DB_PARENT="$BPDP_ROOT/db"

mkdir -p "$DB_PARENT" provenance

if [[ -f "$BAKTA_DB/version.json" ]]; then
  # Already downloaded. The free-space precondition below guards the download
  # only -- applying it here would refuse to record provenance for a database
  # that is already sitting on disk, which is exactly the state the download
  # leaves the machine in once the space has been spent.
  echo "database already present at $BAKTA_DB — not re-downloading"
else
  need_gb=8
  avail_gb=$(df -BG --output=avail "$BPDP_ROOT" | tail -1 | tr -dc '0-9')
  if (( avail_gb < need_gb )); then
    echo "FATAL: ${avail_gb}G free, need at least ${need_gb}G for the ${DB_TYPE} database." >&2
    exit 1
  fi

  echo "=== available database versions ==="
  run_bakta_db list | tee provenance/bakta_db_available.txt

  echo "=== downloading db-${DB_TYPE} into $DB_PARENT ==="
  run_bakta_db download --output "$DB_PARENT" --type "$DB_TYPE"
fi

if [[ ! -f "$BAKTA_DB/version.json" ]]; then
  echo "FATAL: expected $BAKTA_DB/version.json after download; it is absent." >&2
  exit 1
fi

echo "=== recording database provenance ==="
{
  echo "recorded: $(date -Iseconds)"
  echo "db_path: $BAKTA_DB"
  echo "db_type: $DB_TYPE"
  echo "db_type_reason: full database does not fit in available disk on this machine"
  echo "version_json: $(tr -d '\n ' < "$BAKTA_DB/version.json")"
  echo "db_size_bytes: $(du -sb "$BAKTA_DB" | cut -f1)"
  echo "db_file_count: $(find "$BAKTA_DB" -type f | wc -l)"
} > provenance/bakta_db_version.txt
cat provenance/bakta_db_version.txt
