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

# REF/ALT skew tables for LCL_hsca.Rmd, written to LCL_project/co-accessibility/hsca/outputs/
#   sbatch LCL_run_hsca.sh                                   # genome-wide (autosomes)
#   sbatch LCL_run_hsca.sh --chrom chr21 --out-dir <dir>     # smoke test
# 1. LCL_hsca_pool.py: 27-sample pooled/noXCI comparisons from existing per-sample counts
# 2. LCL_hsca_run.py: REF/ALT counts and tests -> four allele tables
# Needs CRE_co-accessibility_FIRE/coaccess/LCL_coaccess_pairs.tsv.gz and LCL_pair_sample_counts.tsv.gz from ../CRE_co-accessibility_FIRE/LCL_run_coaccess.sh.
set -euo pipefail
cd /project/spott/cshan/fiber-seq/code/co-accessibility/LCL/hsca
PYTHON=/project/spott/cshan/envs/Jupyter-notebook/bin/python3
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1

# Pooling and allele tests share the same root, cohort, and output directory.
# --chrom and --workers apply only to the main HSCA tests.
shared_args=()
args=("$@")
for ((i=0; i<${#args[@]}; i++)); do
    case "${args[$i]}" in
        --root|--sample-metatable|--out-dir)
            shared_args+=("${args[$i]}" "${args[$((i + 1))]}")
            i=$((i + 1)) ;;
        --root=*|--sample-metatable=*|--out-dir=*) shared_args+=("${args[$i]}") ;;
    esac
done
"$PYTHON" LCL_hsca_pool.py "${shared_args[@]}"
"$PYTHON" LCL_hsca_run.py "$@"
