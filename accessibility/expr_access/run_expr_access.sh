#!/bin/bash
#SBATCH --job-name=expr_access
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/expr_access_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/expr_access_%j.err

# TSS m6A metaprofiles by expression bin: expression bins (R), profiles
# (python), plots (R). Extra arguments go to 02_tss_m6a_profiles.py, e.g.
#   ./run_expr_access.sh --chrom chr21     # single-chromosome test

set -uo pipefail


module load R/4.4.1
module load python

PYTHON=/software/python-anaconda-2022.05-el8-x86_64/bin/python
CODE=/project/spott/cshan/fiber-seq/code/accessibility/expr_access

Rscript "$CODE/01_expression_bins.R" || exit 1
"$PYTHON" "$CODE/02_tss_m6a_profiles.py" "$@" || exit 1
Rscript "$CODE/03_plot_expr_access.R"
