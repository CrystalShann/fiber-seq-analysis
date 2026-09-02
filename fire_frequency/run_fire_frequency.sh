#!/bin/bash
#SBATCH --job-name=fire_freq
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/fire_freq_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/fire_freq_%j.err

# Per-region FIRE frequency across the LPS timecourse.
# Requires 01_make_cre_universe.sh and 02_read_spans.sh (co-accessibility) to have
# finished: uses their fire_peaks_union.bed and <s>.read_spans.bed.gz.
#
# Usage:
#   sbatch run_fire_frequency.sh
#   sbatch run_fire_frequency.sh --chrom chr21     # single-chromosome test

set -uo pipefail

PYTHON=/project/spott/cshan/envs/Jupyter-notebook/bin/python3

"$PYTHON" /project/spott/cshan/fiber-seq/code/fire_frequency/01_fire_frequency.py "$@"
