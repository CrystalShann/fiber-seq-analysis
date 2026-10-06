#!/bin/bash -l
#SBATCH --job-name=count_footprints_and_background
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=30:00:00
#SBATCH --chdir=/project/spott/cshan/fiber-seq
#SBATCH --output=/project/spott/cshan/fiber-seq/code/enhancer/logs/enhancer_footprints_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/code/enhancer/logs/enhancer_footprints_%j.err

# The Rmd submission chunk writes the input snapshot passed as the first argument.
set -euo pipefail
PROJECT=${ENHANCER_PROJECT_ROOT:-/project/spott/cshan/fiber-seq}
CONFIG=${1:?Supply the input RDS written by the Rmd submission chunk}
module load R/4.4.1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export OPENBLAS_NUM_THREADS=$OMP_NUM_THREADS MKL_NUM_THREADS=$OMP_NUM_THREADS
export R_DATATABLE_NUM_THREADS=$OMP_NUM_THREADS
Rscript --vanilla "$PROJECT/code/enhancer/co_occ_functions/run_enhancer_footprints.R" "$CONFIG"
