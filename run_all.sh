#!/usr/bin/env bash
# ============================================================================
# Reproduce every figure and table in the Data Descriptor from the GEO files.
#
#   bash run_all.sh            # full run
#   bash run_all.sh 10 11 12   # only the stages whose number matches a prefix
#
# Requires: python (numpy, pandas, openpyxl, matplotlib, python-docx, PyMuPDF,
#           Pillow, scipy) and R (limma, jsonlite). See environment.yml.
# Run from the repository root.
# ============================================================================
set -euo pipefail

PY="${PYTHON:-python}"
RS="${RSCRIPT:-Rscript}"
LOG="${LOG:-run_all.log}"
: > "$LOG"

# stage scripts, in dependency order
STAGES=(
  "00_prepare_inputs.py"
  "01_prepare_inputs_large.py"
  "02_prepare_inputs_wide.py"
  "03_probe_map_mrna.py"
  "04_probe_map_lncrna_others.py"
  "10_limma.R"
  "11_normalisation_sensitivity.R"
  "12_detection_flags.py"
  "13_pctmod_detected.R"
  "20_pca_standard.R"
  "21_isg_detection.py"
  "30_export_figdata.R"
  "40_figures.py"
  "50_table4_data.py"
  "60_build_tables.py"
  "61_build_table_docx.py"
  "70_check_numbers.py"
)

selected=()
if [ "$#" -gt 0 ]; then
  for s in "${STAGES[@]}"; do
    for want in "$@"; do
      case "$s" in "$want"*) selected+=("$s");; esac
    done
  done
else
  selected=("${STAGES[@]}")
fi

echo "repository : $(pwd)"
echo "python     : $($PY -c 'import sys; print(sys.version.split()[0])')"
echo "R          : $($RS -e 'cat(R.version.string)' 2>/dev/null || echo 'not found')"
echo "log        : $LOG"
echo

fail=0
for s in "${selected[@]}"; do
  printf '>>> %-32s ' "$s"
  start=$(date +%s)
  case "$s" in
    *.py) run=("$PY" "scripts/$s");;
    *.R)  run=("$RS" "scripts/$s");;
  esac
  if "${run[@]}" >>"$LOG" 2>&1; then
    printf 'ok  (%ss)\n' "$(( $(date +%s) - start ))"
  else
    printf 'FAILED  (see %s)\n' "$LOG"
    tail -25 "$LOG"
    fail=1
    break
  fi
done

echo
if [ "$fail" -eq 0 ]; then
  echo "all stages completed. Outputs:"
  echo "  results/  figdata/  figures/  tables/"
else
  echo "stopped at the first failing stage. Full log: $LOG"
fi
exit "$fail"
