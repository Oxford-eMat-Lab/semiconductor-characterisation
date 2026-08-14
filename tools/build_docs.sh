#!/usr/bin/env bash
# Regenerates the notebook-derived pages of the docs site, then builds it.
#
#   ./tools/build_docs.sh            # regenerate pages + mkdocs build
#   ./tools/build_docs.sh serve      # regenerate pages + mkdocs serve
#
# docs/index.md is hand-written; docs/eqe.md and docs/tlm.md are generated
# from the notebooks and should not be edited directly - edit the notebook
# (for EQE, edit EQE/build_notebook.py and re-run it) and run this again.
set -euo pipefail
cd "$(dirname "$0")/.."

REPO_URL="${REPO_URL:-https://github.com/YOUR_GH_USERNAME/semicon_characterisation}"

python3 tools/nb2md.py EQE/eqe_analysis.ipynb docs/eqe.md \
  --assets docs/assets/nb --assets-url assets/nb \
  --rewrite './figures/=../assets/' \
  --rewrite './eqe_helper.py='"$REPO_URL"'/blob/main/EQE/eqe_helper.py' \
  --repo-url "$REPO_URL"

python3 tools/nb2md.py TLM/contact_res_v2.ipynb docs/tlm.md \
  --assets docs/assets/nb --assets-url assets/nb \
  --rewrite './tlm2.jpg=../assets/tlm2.jpg' \
  --rewrite './tlm3.jpg=../assets/tlm3.jpg' \
  --repo-url "$REPO_URL"

if [ "${1:-build}" = "serve" ]; then
  exec mkdocs serve
else
  exec mkdocs build --strict
fi
