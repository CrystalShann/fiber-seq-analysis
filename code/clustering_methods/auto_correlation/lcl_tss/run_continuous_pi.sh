#!/bin/bash
#SBATCH --job-name=lcl_tss_cpi
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/lcl_tss_cpi_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/lcl_tss_cpi_%j.err

#
# Usage:
#   ./run_continuous_pi.sh 04       # gene selection, seconds
#   ./run_continuous_pi.sh submit   # one sbatch job: 04, then 05 (ACF of all genes + tables)
#   ./run_continuous_pi.sh 05       # stage 05 in the current shell
#   ./run_continuous_pi.sh 05v      # the 12 pilot genes vs the pilot tables, after 05
#   ./run_continuous_pi.sh 06       # PI regression on the stage 05 gene table, seconds
# Stage 05 is resumable: a rerun only computes genes not yet in outputs/acf_all_genes.rds.
set -euo pipefail

# absolute, not $BASH_SOURCE: under sbatch the script runs from the SLURM spool copy
HERE=/project/spott/cshan/fiber-seq/code/clustering_methods/auto_correlation/lcl_tss

export LD_LIBRARY_PATH="/project/spott/cshan/envs/dimelo/lib:/software/gcc-12.2.0-el8-x86_64/lib64:/software/openblas-0.3.29-el8-x86_64/lib:/software/glpk-5.0-el8-x86_64/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export R_LIBS_USER="/project/spott/cshan/Rlibs/x86_64-pc-linux-gnu-library/4.4"
export RETICULATE_PYTHON="/project/spott/cshan/envs/Jupyter-notebook/bin/python"
# keep the sourced Python single-threaded, as the parent pipeline does
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

RSCRIPT="${LCL_TSS_RSCRIPT:-/software/R-4.4.1-el8-x86_64/bin/Rscript}"

[ $# -gt 0 ] || { echo "usage: $0 <stage> [<stage> ...]  (04 submit 05 05v 06)" >&2; exit 2; }

for s in "$@"; do
  echo "=== stage $s ==="
  case "$s" in
    04)     "$RSCRIPT" --vanilla "$HERE/04_select_all_genes.R" ;;
    05)     "$RSCRIPT" --vanilla "$HERE/05_tss_acf_all_genes.R" ;;
    05v)    "$RSCRIPT" --vanilla "$HERE/05_tss_acf_all_genes.R" validate ;;
    06)     "$RSCRIPT" --vanilla "$HERE/06_pi_regression.R" ;;
    submit) sbatch "$HERE/run_continuous_pi.sh" 04 05 ;;
    *)      echo "unknown stage '$s' (use 04 submit 05 05v 06)" >&2; exit 2 ;;
  esac
done
echo "=== done ==="
