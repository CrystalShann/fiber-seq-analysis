#!/bin/bash
#SBATCH --job-name=fire_codep
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/fire_codep_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/fire_codep_%j.err

# Full-genome FIRE peak codependency across the LPS timecourse.
# Requires intersect_msp_fire.sh to have finished for all four timepoints
# (LPS_0 first, since it builds the peak universe the others intersect against).
#
# Usage:
#   sbatch run_fire_codependency.sh
#   sbatch run_fire_codependency.sh --chrom chr21   # extra args pass through

set -uo pipefail

PYTHON=/project/spott/cshan/envs/Jupyter-notebook/bin/python3

"$PYTHON" /project/spott/cshan/fiber-seq/code/FIRE_MSP/fire_codependency.py \
    --shuffle-control "$@"
