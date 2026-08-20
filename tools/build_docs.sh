#!/usr/bin/env bash
# Thin wrapper around tools/build_docs.py, which holds the actual logic.
# Kept so CI and anyone on a shell can keep calling the same command.
#
#   ./tools/build_docs.sh            # regenerate pages + mkdocs build
#   ./tools/build_docs.sh serve      # regenerate pages + mkdocs serve
#
# On Windows without bash, call the Python script directly:
#   python tools\build_docs.py
set -euo pipefail
exec python3 "$(dirname "$0")/build_docs.py" "$@"
