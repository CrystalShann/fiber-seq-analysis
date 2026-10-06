#!/bin/bash -l
#SBATCH --job-name=cluster_fibers_acf_pooled_capped
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=160G
#SBATCH --time=30:00:00
#SBATCH --chdir=/project/spott/cshan/fiber-seq
#SBATCH --output=/project/spott/cshan/fiber-seq/code/enhancer/logs/enhancer_acf_pooled_capped_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/code/enhancer/logs/enhancer_acf_pooled_capped_%j.err

# Submit: sbatch code/enhancer/shell_scripts/run_enhancer_acf.sh
# Submit after the shared-fiber preparation job completes, or with afterok.
# ENHANCER_SHARED_FIBERS optionally overrides the prepared cohort RDS path.
# Both classes and all timepoints enter one ACF fit using that shared cohort;
# there is no extraction or resampling here. Constant fibers stay unclustered.
# Writes tables/enhancer_acf_results_pooled_capped.rds and automatic
# acf_*_pooled_capped TSV/PDF outputs. No notebook execution is needed.
set -euo pipefail
PROJECT=${ENHANCER_PROJECT_ROOT:-/project/spott/cshan/fiber-seq}
module load R/4.4.1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export OPENBLAS_NUM_THREADS=$OMP_NUM_THREADS MKL_NUM_THREADS=$OMP_NUM_THREADS
export NUMBA_NUM_THREADS=$OMP_NUM_THREADS R_DATATABLE_NUM_THREADS=$OMP_NUM_THREADS
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
Rscript --vanilla "$PROJECT/code/enhancer/acf/run_enhancer_acf.R"
