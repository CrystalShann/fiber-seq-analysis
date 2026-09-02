#!/bin/bash
#SBATCH --job-name=m6a_binned
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=8:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/m6a_binned_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/m6a_binned_%j.err

# Average m6A level in 10 bp bins across the union FIRE regions.
# Requires the ft extract m6a_by_chr files (FiberHMM/extract) and
# fire_peaks_union.bed (co-accessibility 01) to exist.
#
# Usage:
#   sbatch run_m6a_binned.sh
#   sbatch run_m6a_binned.sh --chrom chr22     # single-chromosome test

set -uo pipefail

PYTHON=/project/spott/cshan/envs/Jupyter-notebook/bin/python3

"$PYTHON" /project/spott/cshan/fiber-seq/code/fire_frequency/04_m6a_binned.py "$@"
