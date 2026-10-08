#!/bin/bash
#SBATCH --job-name=feature_count
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=8:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/feature_count_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/feature_count_%j.err

# Per-region FIRE element and FiberHMM footprint counts across the LPS timecourse,
# on the same spanning reads as 03_fire_frequency.py (joins one-to-one on
# region_id + timepoint). Requires 01_make_fire_universe.sh and 02_read_spans.sh
# to have finished: uses their fire_peaks_union.bed and <s>.read_spans.bed.gz.
#
# Usage:
#   sbatch run_feature_count.sh
#   sbatch run_feature_count.sh --chrom chr21     # single-chromosome test

set -uo pipefail

PYTHON=/project/spott/cshan/envs/Jupyter-notebook/bin/python3

"$PYTHON" /project/spott/cshan/fiber-seq/code/fire_frequency/03b_feature_count.py "$@"
