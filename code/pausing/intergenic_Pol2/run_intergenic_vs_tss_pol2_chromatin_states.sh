#!/bin/bash
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=6:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/%x_%A_%a.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/%x_%A_%a.err

# Intergenic vs TSS Pol II chromatin states on caslake.
#
# Usage (from a login node):
#   bash run_intergenic_vs_tss_pol2_chromatin_states.sh
# submits: states -> m6a array (chr1-22, X, Y) -> summarize, each after the previous succeeds.
# Single steps: bash run_intergenic_vs_tss_pol2_chromatin_states.sh states|m6a|summarize

set -euo pipefail
SELF=/project/spott/cshan/fiber-seq/code/pausing/intergenic_Pol2/run_intergenic_vs_tss_pol2_chromatin_states.sh
SCRIPT=/project/spott/cshan/fiber-seq/code/pausing/intergenic_Pol2/intergenic_vs_tss_pol2_chromatin_states.r
CHROMS=(chr{1..22} chrX chrY)

if [ -z "${SLURM_JOB_ID:-}" ]; then
    step=${1:-all}
    if [ "$step" = all ]; then
        s=$(sbatch --parsable --job-name=pol2_states "$SELF" states)
        a=$(sbatch --parsable --job-name=pol2_m6a --array=1-${#CHROMS[@]} --dependency=afterok:$s "$SELF" m6a)
        z=$(sbatch --parsable --job-name=pol2_m6a_sum --dependency=afterok:$a "$SELF" summarize)
        echo "states $s -> m6a $a -> summarize $z"
    elif [ "$step" = m6a ]; then
        sbatch --job-name=pol2_m6a --array=1-${#CHROMS[@]} "$SELF" m6a
    else
        sbatch --job-name=pol2_$step "$SELF" "$step"
    fi
    exit 0
fi

module load R/4.4.1
case "$1" in
    states)    Rscript "$SCRIPT" states ;;
    m6a)       Rscript "$SCRIPT" m6a "${CHROMS[$((SLURM_ARRAY_TASK_ID - 1))]}" ;;
    summarize) Rscript "$SCRIPT" summarize ;;
    *)         echo "unknown step: $1" >&2; exit 1 ;;
esac
