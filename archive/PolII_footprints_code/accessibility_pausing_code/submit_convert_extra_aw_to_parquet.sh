#!/bin/bash

# Convert the per-window accessibility tsv.gz outputs (aw250, aw500, aw1000)
# produced by submit_extra_aw_runs.sh into parquet. Run after the R jobs have
# completed. The existing aw100 (.parquet) file is not touched.

#SBATCH --job-name=acc_to_parquet_aw
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/acc_to_parquet_aw_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/acc_to_parquet_aw_%j.err

set -euo pipefail

PYTHON=/project/spott/cshan/envs/Jupyter-notebook/bin/python
SCRIPT=/project/spott/cshan/fiber-seq/code/PolII_footprints_code/accessibility_pausing_code/convert_accessibility_to_parquet.py
TABLE_DIR=/project/spott/cshan/fiber-seq/results/PolII/m6a_pausing_quartiles/AL10_bc2178_19130_1kb_bin10_modthresh0.9_all/tables
SAMPLE=AL10_bc2178_19130

for AW in 250 500 1000; do
  src="${TABLE_DIR}/${SAMPLE}_tss_read_accessibility_aw${AW}bp.tsv.gz"
  dst="${TABLE_DIR}/${SAMPLE}_tss_read_accessibility_aw${AW}bp.parquet"
  if [[ ! -f "${src}" ]]; then
    echo "[skip] ${src} not found" >&2
    continue
  fi
  echo "Converting aw=${AW} bp ..."
  "${PYTHON}" "${SCRIPT}" "${src}" "${dst}"
done
