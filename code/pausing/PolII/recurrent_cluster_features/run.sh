#!/bin/bash
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/recurrent_cluster_features_%j.log
set -euo pipefail
cd /project/spott/cshan/fiber-seq
bash code/pausing/PolII/run_gene_examples.sh code/pausing/PolII/recurrent_cluster_features/01_summarize_recurrent_patterns.R
bash code/pausing/PolII/run_gene_examples.sh code/pausing/PolII/recurrent_cluster_features/03_cluster_gene_summary.R
