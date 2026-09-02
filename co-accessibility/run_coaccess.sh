#!/bin/bash
#SBATCH --job-name=coaccess
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/coaccess_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/coaccess_%j.err

# cCRE co-accessibility across the LPS timecourse.
# Requires 01_make_cre_universe.sh and 02_read_spans.sh to have finished.
#
# Usage:
#   sbatch run_coaccess.sh
#   sbatch run_coaccess.sh --max-dist 10000        # Kevin's single-sample range
#   sbatch run_coaccess.sh --read-rule contain     # require fibers to span the cCRE
#   sbatch run_coaccess.sh --kevin-compat          # drop zero-accessibility pairs

set -uo pipefail

PYTHON=/project/spott/cshan/envs/Jupyter-notebook/bin/python3

"$PYTHON" /project/spott/cshan/fiber-seq/code/co-accessibility/03_coaccess_cres.py "$@"
