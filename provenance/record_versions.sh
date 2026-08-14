#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
{
  echo "recorded: $(date -Iseconds)"
  echo "bakta:  $(bakta --version 2>&1 | head -1)"
  echo "prokka: $(prokka --version 2>&1 | head -1)"
  echo "python: $(python3 --version)"
  echo "sklearn: $(python3 -c 'import sklearn;print(sklearn.__version__)' 2>/dev/null || echo absent)"
} > tool_versions.txt
cat tool_versions.txt
