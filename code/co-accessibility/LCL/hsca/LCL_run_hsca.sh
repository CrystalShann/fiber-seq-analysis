#!/bin/bash
#SBATCH --job-name=LCL_hsca
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/LCL_hsca_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/LCL_hsca_%j.err

# Haplotype / allele skew tables for LCL_hsca.Rmd, written to LCL_project/co-accessibility/hsca/outputs/
#   sbatch LCL_run_hsca.sh                                   # genome-wide (autosomes)
#   sbatch LCL_run_hsca.sh --chrom chr21 --out-dir <dir>     # smoke test
# 1. LCL_hsca_intersect.py: FIRE hap-differences peaks x paired cCREs -> LCL_hsca_cre_status.tsv.gz
# 2. LCL_hsca_run.py: H1/H2 and REF/ALT counts and tests -> the other five tables
# Needs coaccess/LCL_coaccess_pairs.tsv.gz and LCL_pair_sample_counts.tsv.gz from LCL_run_coaccess.sh.
set -euo pipefail
cd /project/spott/cshan/fiber-seq/code/co-accessibility/LCL/hsca
PYTHON=/project/spott/cshan/envs/Jupyter-notebook/bin/python3
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

# the intersect step is genome-wide and only takes --out-dir
out_dir=()
args=("$@")
for i in "${!args[@]}"; do
    [[ ${args[$i]} == --out-dir ]] && out_dir=(--out-dir "${args[$((i + 1))]}")
done
"$PYTHON" LCL_hsca_intersect.py "${out_dir[@]}"
"$PYTHON" LCL_hsca_run.py "$@"
