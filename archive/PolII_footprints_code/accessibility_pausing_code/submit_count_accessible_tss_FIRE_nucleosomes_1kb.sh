#!/bin/bash
#SBATCH --job-name=tss_FIRE_nuc
#SBATCH --account=pi-spott
#SBATCH --partition=bigmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=160G
#SBATCH --time=24:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/tss_FIRE_nuc_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/tss_FIRE_nuc_%j.err

set -euo pipefail

module load R

R_SCRIPT="/project/spott/cshan/fiber-seq/code/PolII_footprints_code/accessibility_pausing_code/count_accessible_tss_FIRE_nucleosomes_1kb.R"
Rscript "$R_SCRIPT"
