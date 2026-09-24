#!/bin/bash
#SBATCH --job-name=all_nuc_PI
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/all_nuc_PI_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/all_nuc_PI_%j.err

set -euo pipefail

module load R

SCRIPT="/project/spott/cshan/fiber-seq/code/PolII_footprints_code/combined_footprints_code/plot_all_nucleosomes_pausing_index.R"

Rscript "${SCRIPT}" \
  --sample="${SAMPLE:-AL10_bc2178_19130}" \
  --nuc-source="${NUC_SOURCE:-CAGE}" \
  --feature="${FEATURE:-Nucleosome}" \
  --window-bp="${WINDOW_BP:-2000}" \
  --bin-size="${BIN_SIZE:-10}" \
  --force-rebuild="${FORCE_REBUILD:-false}"
