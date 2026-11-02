#!/usr/bin/env bash
# Rerun every CDI-dependent analysis on the corrected rasters.
# usage: bash run_corrected_pipeline.sh ROOT DATA VALIDATION_CSV
#   ROOT            output folder of cdi_corrected_run.py
#   DATA            folder with counties/ (Census cartographic boundaries), usdm_county_dsci.csv,
#                   nd_yields.tsv (NASS county yields) and climdiv/ (nClimDiv division files)
#   VALIDATION_CSV  annual reference table (USDM, nClimDiv and state yield anomalies by year)
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$1"; DATA="$2"; VALID="$3"
PY="python -P"
export PYTHONIOENCODING=utf-8 PYTHONWARNINGS=ignore

# Tables 1-3, Figures 7-11
$PY "$HERE/analyze_cdi_corrected.py" "$ROOT"

# validation table: same references, CDI columns from the corrected rasters
mkdir -p "$ROOT/validation"
$PY - "$VALID" "$ROOT" <<'EOF'
import sys, pandas as pd
V = pd.read_csv(sys.argv[1]).set_index("year")
S = pd.read_csv(sys.argv[2] + "/analysis/CDI_summary_statistics.csv").set_index("year")
V["cdi_mean"] = S.mean_cdi.round(4); V["cdi_median"] = S.median_cdi.round(4)
V["cdi_pct_ge_q75"] = S.high_cdi_pixels_pct_global_q75.round(2)
V.to_csv(sys.argv[2] + "/validation/CDI_validation_dataset.csv")
EOF

$PY "$HERE/analyze.py" "$ROOT" "$DATA/climdiv" > "$ROOT/analysis/benchmark.json"
$PY "$HERE/county_benchmark.py" "$ROOT" "$DATA" "$ROOT/county_output" > "$ROOT/analysis/county.txt"
$PY "$HERE/variance_decomposition.py" "$ROOT" > "$ROOT/analysis/variance.txt"
$PY "$HERE/fig6.py" "$ROOT" "$ROOT/analysis/Figure6_corrected.png" > "$ROOT/analysis/fig6.txt"
$PY "$HERE/checks_corrected.py" "$ROOT" > "$ROOT/analysis/checks.json"
echo PIPELINE_DONE
