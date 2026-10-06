#!/bin/bash -l
#SBATCH --job-name=cluster_fibers_manhattan_pooled_capped
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=160G
#SBATCH --time=30:00:00
#SBATCH --chdir=/project/spott/cshan/fiber-seq
#SBATCH --output=/project/spott/cshan/fiber-seq/code/enhancer/logs/enhancer_man_pooled_capped_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/code/enhancer/logs/enhancer_man_pooled_capped_%j.err

# Submit after the shared-fiber preparation job, with --dependency=afterok:<job_id>.
# One job clusters active + inactive fibers and all four timepoints together.
# Reads ENHANCER_SHARED_FIBERS, or tables/enhancer_shared_fibers_pooled_capped.rds.
# Writes enhancer_manhattan_results_pooled_capped.rds; sampling is never repeated here.
# ENHANCER_MANHATTAN_NEIGHBORS=approximate (default) or exact (quadratic time).
# Use k=100 neighbors with raw 1-bp distances, without PCA or binning.
# Save contribution diagnostics and *_pooled_capped PDFs; display all selected fibers
# in one footprint figure with cluster blocks and cluster/timepoint/class color bars.
set -euo pipefail
PROJECT=${ENHANCER_PROJECT_ROOT:-/project/spott/cshan/fiber-seq}
export ENHANCER_MANHATTAN_K=${ENHANCER_MANHATTAN_K:-100}
module load R/4.4.1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export OPENBLAS_NUM_THREADS=$OMP_NUM_THREADS MKL_NUM_THREADS=$OMP_NUM_THREADS
export NUMBA_NUM_THREADS=$OMP_NUM_THREADS R_DATATABLE_NUM_THREADS=$OMP_NUM_THREADS
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
Rscript --vanilla "$PROJECT/code/enhancer/manhattan/run_enhancer_manhattan.R"
