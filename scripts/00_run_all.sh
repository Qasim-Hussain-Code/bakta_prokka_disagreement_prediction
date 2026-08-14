#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
for s in scripts/[0-9][0-9]_*; do
  [[ "$s" == *00_run_all.sh ]] && continue
  echo "=== $s ==="
  case "$s" in
    *.py) python3 "$s" ;;
    *.sh) bash "$s" ;;
  esac
done
