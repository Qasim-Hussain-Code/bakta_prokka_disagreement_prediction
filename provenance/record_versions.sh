#!/usr/bin/env bash
# Record the tool versions actually used, resolved the same way the pipeline
# resolves them. Calling bare `bakta`/`prokka` here would record "command not
# found" on any machine where they live in a cached env rather than on PATH,
# which is exactly the machine this pipeline is built to run on.
set -euo pipefail
source "$(dirname "$0")/../scripts/lib_env.sh"
cd "$BPDP_ROOT/provenance"
{
  echo "recorded: $(date -Iseconds)"
  echo "bakta:  $(run_bakta --version 2>&1 | head -1)"
  echo "bakta_bin: $BAKTA_BIN"
  echo "prokka: $(run_prokka --version 2>&1 | head -1)"
  echo "prokka_bin: $PROKKA_BIN"
  echo "bakta_db: $BAKTA_DB"
  echo "python: $(python3 --version)"
  echo "sklearn: $(python3 -c 'import sklearn;print(sklearn.__version__)' 2>/dev/null || echo absent)"
  echo "numpy: $(python3 -c 'import numpy;print(numpy.__version__)' 2>/dev/null || echo absent)"
  echo "matplotlib: $(python3 -c 'import matplotlib;print(matplotlib.__version__)' 2>/dev/null || echo absent)"
} > tool_versions.txt
cat tool_versions.txt
