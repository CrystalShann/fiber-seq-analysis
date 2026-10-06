#!/bin/bash -l
#SBATCH --job-name=enhancer_shared_sample
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=160G
#SBATCH --time=30:00:00
#SBATCH --chdir=/project/spott/cshan/fiber-seq
#SBATCH --output=/project/spott/cshan/fiber-seq/code/enhancer/logs/enhancer_shared_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/code/enhancer/logs/enhancer_shared_%j.err

# Submit through submit_enhancer_pooled_jobs.sh to schedule both dependent methods.
# Sampling defaults/environment overrides live in run_enhancer_shared_sampling.R.
set -euo pipefail
PROJECT=${ENHANCER_PROJECT_ROOT:-/project/spott/cshan/fiber-seq}
module load R/4.4.1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export OPENBLAS_NUM_THREADS=$OMP_NUM_THREADS MKL_NUM_THREADS=$OMP_NUM_THREADS
export NUMBA_NUM_THREADS=$OMP_NUM_THREADS R_DATATABLE_NUM_THREADS=$OMP_NUM_THREADS
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
Rscript --vanilla "$PROJECT/code/enhancer/shared_functions/run_enhancer_shared_sampling.R"
