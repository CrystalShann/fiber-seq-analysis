#!/bin/bash -l
#SBATCH --job-name=enhancer_man
#SBATCH --account=pi-spott
#SBATCH --partition=bigmem
#SBATCH --array=1-2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=700G
#SBATCH --time=30:00:00
#SBATCH --chdir=/project/spott/cshan/fiber-seq
#SBATCH --output=/project/spott/cshan/fiber-seq/code/enhancer/enhancer_man_%A_%a.out
#SBATCH --error=/project/spott/cshan/fiber-seq/code/enhancer/enhancer_man_%A_%a.err

# Submit: sbatch code/enhancer/run_enhancer_manhattan.sh
# Task 1: active; task 2: inactive. Each reads its saved 5,000 sampled enhancers.
# Each pools all four timepoints and writes enhancer_manhattan_<class>_results.rds.
# ENHANCER_MANHATTAN_NEIGHBORS=approximate (default) or exact (quadratic time).
# Uniformly sample up to 10,000 fibers per class across all four timepoints.
# Use k=100 neighbors with raw 1-bp distances, without PCA or binning.
# Save selected m6A/footprint data and heatmap, composition and footprint PDFs.
set -euo pipefail
PROJECT=${ENHANCER_PROJECT_ROOT:-/project/spott/cshan/fiber-seq}
case "${SLURM_ARRAY_TASK_ID:-}" in
  1) export ENHANCER_MANHATTAN_CLASS=active ;;
  2) export ENHANCER_MANHATTAN_CLASS=inactive ;;
  *) echo "Expected SLURM_ARRAY_TASK_ID=1 (active) or 2 (inactive)" >&2; exit 1 ;;
esac
export ENHANCER_MANHATTAN_MAX_FIBERS=10000
export ENHANCER_MANHATTAN_K=100
module load R/4.4.1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export OPENBLAS_NUM_THREADS=$OMP_NUM_THREADS MKL_NUM_THREADS=$OMP_NUM_THREADS
export NUMBA_NUM_THREADS=$OMP_NUM_THREADS R_DATATABLE_NUM_THREADS=$OMP_NUM_THREADS
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
Rscript --vanilla "$PROJECT/code/enhancer/run_enhancer_manhattan.R"
