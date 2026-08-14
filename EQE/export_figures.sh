#!/usr/bin/env bash
# Renders figures.pptx to the JPEG files used by eqe_analysis.ipynb and the
# docs site. Requires LibreOffice (soffice) and pdftoppm (poppler-utils).
#
# Workflow:
#   node make_figures.js     # rebuild figures.pptx from source
#   ./export_figures.sh      # re-export the JPEGs
#
# Edit the slides in PowerPoint instead if you prefer - just re-run this
# script afterwards to refresh the JPEGs.
set -euo pipefail
cd "$(dirname "$0")"

# Slide order in figures.pptx -> output filename
NAMES=(
  fig_cell_structure
  fig_optical_losses
  fig_collection_efficiency
  fig_albsf_vs_perc
  fig_measurement_setup
)

mkdir -p figures
rm -f figures.pdf _slide-*.jpg

soffice --headless --convert-to pdf figures.pptx >/dev/null 2>&1
pdftoppm -jpeg -r 150 figures.pdf _slide

for i in "${!NAMES[@]}"; do
  n=$((i + 1))
  src="_slide-${n}.jpg"
  [ -f "$src" ] || { echo "missing $src"; exit 1; }
  # trim the uniform white margin around each rendered slide
  convert "$src" -trim +repage -bordercolor white -border 12 \
    "figures/${NAMES[$i]}.jpg"
  echo "figures/${NAMES[$i]}.jpg"
done

rm -f _slide-*.jpg figures.pdf

# keep the docs site copies in sync
if [ -d ../docs/assets ]; then
  cp figures/*.jpg ../docs/assets/
  echo "copied to ../docs/assets/"
fi
