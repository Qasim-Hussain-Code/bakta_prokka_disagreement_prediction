#!/usr/bin/env bash
# Full pipeline, from an empty checkout.
#
# Part one (01-11) is globbed in numeric order. Part two is delegated to
# run_content_only.sh rather than globbed, because its dependency order is not
# its numeric order -- see that file.
set -euo pipefail
cd "$(dirname "$0")/.."

for s in scripts/0[1-9]_* scripts/1[01]_*; do
  echo "=== $s ==="
  case "$s" in
    *.py) python3 "$s" ;;
    *.sh) bash "$s" ;;
  esac
done

echo "=== scripts/run_content_only.sh ==="
bash scripts/run_content_only.sh
