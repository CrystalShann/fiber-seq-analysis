#!/bin/bash -l
#SBATCH --job-name=enhancer_active_inactive_cobinding
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=30:00:00
#SBATCH --chdir=/project/spott/cshan/fiber-seq
#SBATCH --output=/project/spott/cshan/fiber-seq/code/enhancer/logs/enhancer_site_pairs_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/code/enhancer/logs/enhancer_site_pairs_%j.err
set -euo pipefail
module load R/4.4.1
export OPENBLAS_NUM_THREADS=1
export R_DATATABLE_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
PROJECT=${ENHANCER_PROJECT_ROOT:-/project/spott/cshan/fiber-seq}
Rscript --vanilla "$PROJECT/code/enhancer/co_occ_functions/run_enhancer_site_pairs.R"
