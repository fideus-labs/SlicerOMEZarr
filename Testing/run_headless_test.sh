#!/usr/bin/env bash
# Run the OMEZarr module self-test in a headless Slicer (Xvfb).
# Usage: Testing/run_headless_test.sh [/path/to/Slicer]
# Set OMEZARR_TEST_REMOTE=1 to also exercise the remote (HTTPS) store test.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SLICER="${1:-${SLICER:-Slicer}}"
LAUNCHER=()
if command -v xvfb-run >/dev/null 2>&1 && [ -z "${OMEZARR_NO_XVFB:-}" ]; then
  LAUNCHER=(xvfb-run -a -s "-screen 0 1280x1024x24")
fi
# Git Bash on Windows: hand native paths to Slicer.
if command -v cygpath >/dev/null 2>&1; then
  HERE="$(cygpath -m "$HERE")"
fi
exec ${LAUNCHER[@]+"${LAUNCHER[@]}"} "$SLICER" --no-splash --testing \
  --additional-module-paths "$HERE/../OMEZarr" \
  --python-script "$HERE/run_module_test.py"
