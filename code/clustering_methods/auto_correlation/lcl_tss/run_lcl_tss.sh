#!/bin/bash
# Run the LCL TSS pausing-index autocorrelation stages.
#
# This wrapper exists only to set the library environment: rtracklayer will not load
# unless a newer libcurl precedes the system one on LD_LIBRARY_PATH (the system R build
# otherwise fails with "libcurl.so: undefined symbol: curl_easy_option_next"), and
# reticulate needs RETICULATE_PYTHON pointed at an interpreter that has numpy/scipy.
#
# Usage:
#   ./run_lcl_tss.sh            # all stages in order
#   ./run_lcl_tss.sh 02 03      # only the named stages
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export LD_LIBRARY_PATH="/project/spott/cshan/envs/dimelo/lib:/software/gcc-12.2.0-el8-x86_64/lib64:/software/openblas-0.3.29-el8-x86_64/lib:/software/glpk-5.0-el8-x86_64/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export R_LIBS_USER="/project/spott/cshan/Rlibs/x86_64-pc-linux-gnu-library/4.4"
export RETICULATE_PYTHON="/project/spott/cshan/envs/Jupyter-notebook/bin/python"
# keep the sourced Python single-threaded, as the parent pipeline does
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

RSCRIPT="${LCL_TSS_RSCRIPT:-/software/R-4.4.1-el8-x86_64/bin/Rscript}"

declare -A STAGE=(
  [00]=00_pause_index.R
  [01]=01_select_tss_genes.R
  [02]=02_tss_autocorrelation.R
  [03]=03_plot_tss_autocorrelation.R
)
ORDER=(00 01 02 03)
TARGETS=("$@")
[ ${#TARGETS[@]} -eq 0 ] && TARGETS=("${ORDER[@]}")

for s in "${TARGETS[@]}"; do
  script="${STAGE[$s]:-}"
  [ -n "$script" ] || { echo "unknown stage '$s' (use 00 01 02 03)" >&2; exit 2; }
  echo "=== stage $s: $script ==="
  "$RSCRIPT" --vanilla "$HERE/$script"
done
echo "=== done ==="
