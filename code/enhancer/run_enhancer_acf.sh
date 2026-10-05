#!/bin/bash -l
#SBATCH --job-name=enhancer_acf
#SBATCH --account=pi-spott
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=300G
#SBATCH --time=30:00:00
#SBATCH --chdir=/project/spott/cshan/fiber-seq
#SBATCH --output=/project/spott/cshan/fiber-seq/code/enhancer/enhancer_acf_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/code/enhancer/enhancer_acf_%j.err

# Submit: sbatch code/enhancer/run_enhancer_acf.sh
# Reads the saved sampled enhancers and writes tables/enhancer_acf_results.rds.
# Settings are in run_enhancer_acf.R. Both classes use all four timepoints;
# classes run sequentially to limit peak memory. No notebook execution is needed.
set -euo pipefail
PROJECT=${ENHANCER_PROJECT_ROOT:-/project/spott/cshan/fiber-seq}
module load R/4.4.1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export OPENBLAS_NUM_THREADS=$OMP_NUM_THREADS MKL_NUM_THREADS=$OMP_NUM_THREADS
export NUMBA_NUM_THREADS=$OMP_NUM_THREADS R_DATATABLE_NUM_THREADS=$OMP_NUM_THREADS
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
Rscript --vanilla "$PROJECT/code/enhancer/run_enhancer_acf.R"
