#!/bin/bash
#SBATCH --job-name=ft_dist
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=8:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/ft_dist_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/ft_dist_%j.err

# Per-read FiberHMM TF footprint counts and per-region footprint-count distributions
# (ft_dist_macrophage.py) on the same regions, reads and footprints as
# 03b_feature_count.py, validated against its feature_count_long.tsv.gz. Requires
# 01_make_fire_universe.sh, 02_read_spans.sh and run_feature_count.sh to have finished.
#
# Usage:
#   sbatch ft_dist_macrophage.sh                    # genome-wide, into fire_frequency/
#   sbatch ft_dist_macrophage.sh --chrom chr21      # test, into fire_frequency/ft_dist_test_chr21/
#   sbatch ft_dist_macrophage.sh --no-per-read      # distributions only
# Any other ft_dist_macrophage.py option is passed through (see --help).

set -euo pipefail

PYTHON=/project/spott/cshan/envs/Jupyter-notebook/bin/python3
SCRIPT=/project/spott/cshan/fiber-seq/code/fire_frequency/ft_dist_macrophage.py
LOG_DIR=/project/spott/cshan/fiber-seq/results/logs

# SLURM opens the log files before this runs, so LOG_DIR must already exist for sbatch;
# the python script creates its output directory
mkdir -p "$LOG_DIR"

echo "Job        : ${SLURM_JOB_ID:-interactive}"
echo "Arguments  : $*"
echo "Started    : $(date)"
"$PYTHON" "$SCRIPT" "$@"
echo "Finished   : $(date)"
