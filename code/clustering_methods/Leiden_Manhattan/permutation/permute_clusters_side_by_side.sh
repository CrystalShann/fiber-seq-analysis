#!/bin/bash
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --job-name=permute_side_by_side
#SBATCH --output=/project/spott/cshan/fiber-seq/LCL_project/Leiden_manhattan/parameters_test/logs/permutation_run.log
set -euo pipefail
export PERMUTE_PLOT_WORKERS="${SLURM_CPUS_PER_TASK:-1}"
bash /project/spott/cshan/fiber-seq/code/Pol2/PolII/run_gene_examples.sh \
  /project/spott/cshan/fiber-seq/code/clustering_methods/Leiden_Manhattan/permutation/permute_clusters.R
