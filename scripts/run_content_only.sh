#!/usr/bin/env bash
# Re-run part two (the annotation-content experiment) on the annotation output
# already on disk. No annotation, no database, no downloads.
#
# Deliberately NOT named 00_run_content.sh: 00_run_all.sh globs
# scripts/[0-9][0-9]_* and would pick that name up, running part two twice and
# doing it before the annotation exists on a fresh machine.
#
# 00_run_all.sh remains the from-scratch entry point and picks up steps 12-21
# through the same glob, in order, after annotation.
set -euo pipefail
cd "$(dirname "$0")/.."

for s in scripts/1[2-9]_*.py scripts/2[0-1]_*.py; do
  echo "=== $s ==="
  python3 "$s"
done
