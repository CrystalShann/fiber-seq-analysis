#!/bin/bash
#SBATCH --job-name=acc_to_parquet
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/acc_to_parquet_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/acc_to_parquet_%j.err

set -euo pipefail

/project/spott/cshan/envs/Jupyter-notebook/bin/python \
  /project/spott/cshan/fiber-seq/code/PolII_footprints_code/accessibility_pausing_code/convert_accessibility_to_parquet.py
