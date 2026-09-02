#!/bin/bash
#SBATCH --job-name=entropy
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/entropy_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/entropy_%j.err

# Accessibility entropy, pooled LPS timepoints: windows + FIRE top-10 (R),
# per-read m6A fractions (python), 4-bin entropy tables + region selection
# (R), then Leiden-cluster-based entropy at the 8 early-response promoters
# (R: pooled Manhattan-kNN Leiden clustering, cluster-proportion entropy per
# timepoint and within-cluster positional entropy, with plots). Extra
# arguments go to 02_read_fractions.py, e.g.
#   sbatch run_entropy.sh --chrom chr21     # single-chromosome test
# (a test run overwrites the tables -- rerun genome-wide afterwards; the
# Leiden step reads the FiberHMM extracts directly and is unaffected)

set -uo pipefail

module load R/4.4.1

PYTHON=/project/spott/cshan/envs/Jupyter-notebook/bin/python3
CODE=/project/spott/cshan/fiber-seq/code/clustering_TSS_methods/entropy

Rscript "$CODE/01_entropy_windows.R" || exit 1
"$PYTHON" "$CODE/02_read_fractions.py" "$@" || exit 1
Rscript "$CODE/03_entropy_violin.R" || exit 1
Rscript "$CODE/04_leiden_entropy.R" || exit 1
