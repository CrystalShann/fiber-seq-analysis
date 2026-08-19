#!/bin/bash

# Submit one sbatch job per extra TSS access window. Each job runs
# run_m6a_pausing_quartiles.R with DO_POSITION_SCAN=0, so only the per-read
# accessibility table for that window is produced. The existing 100 bp
# accessibility output is not touched.
#
# Usage:
#   bash submit_extra_aw_runs.sh             # default: 250 500 1000
#   bash submit_extra_aw_runs.sh 1000        # just 1000
#   bash submit_extra_aw_runs.sh 250 500     # subset
#
# Per-window memory (override with env var MEM_<window>=<size>):
#   250  -> 300G
#   500  -> 300G
#   1000 -> 600G   (300G OOM'd at chr1 chunk 1 in job 49377331)

set -euo pipefail

JOB_SCRIPT="$(dirname "$0")/run_m6a_pausing_quartiles_extra_aw.sh"

if [[ $# -ge 1 ]]; then
  WINDOWS=("$@")
else
  WINDOWS=(250 500 1000)
fi

mem_for() {
  case "$1" in
    250)  echo "${MEM_250:-300G}"  ;;
    500)  echo "${MEM_500:-300G}"  ;;
    1000) echo "${MEM_1000:-600G}" ;;
    *)    echo "${MEM_DEFAULT:-300G}" ;;
  esac
}

for AW in "${WINDOWS[@]}"; do
  MEM="$(mem_for "${AW}")"
  echo "Submitting access window ${AW} bp (mem=${MEM})..."
  sbatch \
    --job-name="m6a_aw${AW}" \
    --mem="${MEM}" \
    --export=ALL,TSS_ACCESS_WINDOW=${AW} \
    "${JOB_SCRIPT}"
done
