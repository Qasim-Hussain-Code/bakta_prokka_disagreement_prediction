#!/usr/bin/env bash
# Re-run part two (the annotation-content experiment) on the annotation output
# already on disk. No annotation, no database, no downloads.
#
# The order is explicit, not a glob, because it is not the numeric order: the
# mechanism check (22) computes numbers the summary (20) quotes, so it has to
# run first. A glob would silently produce a summary with that block missing.
#
# Deliberately NOT named 00_run_content.sh -- 00_run_all.sh globs
# scripts/[0-9][0-9]_* and would pick that name up and run part two twice.
set -euo pipefail
cd "$(dirname "$0")/.."

for s in 12_content_cohort 13_name_rules 14_content_features \
         15_content_baselines 16_content_cv 17_content_final \
         18_content_importance 19_content_circularity \
         22_content_mechanism 20_content_summary 21_content_figures; do
  echo "=== scripts/${s}.py ==="
  python3 "scripts/${s}.py"
done
