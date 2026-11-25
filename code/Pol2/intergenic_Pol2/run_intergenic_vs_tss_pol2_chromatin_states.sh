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
#   bash run_intergenic_vs_tss_pol2_chromatin_states.sh rank
# submits: fire_fp array (31 samples) -> fire_fp_summary -> states -> summarize -> plots
# (reuses tables/m6a; the m6A region selection does not depend on the ranks).
# Single steps: bash run_intergenic_vs_tss_pol2_chromatin_states.sh
#   fire_fp|fire_fp_summary|states|m6a|summarize|plots

set -euo pipefail
SELF=/project/spott/cshan/fiber-seq/code/Pol2/intergenic_Pol2/run_intergenic_vs_tss_pol2_chromatin_states.sh
DIR=/project/spott/cshan/fiber-seq/code/Pol2/intergenic_Pol2
SCRIPT=$DIR/intergenic_vs_tss_pol2_chromatin_states.r
PYTHON=/project/spott/cshan/envs/Jupyter-notebook/bin/python3
MANIFEST=/project/spott/cshan/fiber-seq/LCL_project/co-accessibility/CRE_co-accessibility_FIRE/universe/sample_manifest.tsv
CHROMS=(chr{1..22} chrX chrY)
mapfile -t SAMPLES < <(tail -n +2 "$MANIFEST" | cut -f1)

if [ -z "${SLURM_JOB_ID:-}" ]; then
    step=${1:-all}
    if [ "$step" = all ]; then
        s=$(sbatch --parsable --job-name=pol2_states "$SELF" states)
        a=$(sbatch --parsable --job-name=pol2_m6a --array=1-${#CHROMS[@]} --dependency=afterok:$s "$SELF" m6a)
        z=$(sbatch --parsable --job-name=pol2_m6a_sum --dependency=afterok:$a "$SELF" summarize)
        echo "states $s -> m6a $a -> summarize $z"
    elif [ "$step" = rank ]; then
        f=$(sbatch --parsable --job-name=pol2_fire_fp --array=1-${#SAMPLES[@]} "$SELF" fire_fp)
        g=$(sbatch --parsable --job-name=pol2_fire_fp_sum --dependency=afterok:$f "$SELF" fire_fp_summary)
        s=$(sbatch --parsable --job-name=pol2_states --dependency=afterok:$g "$SELF" states)
        z=$(sbatch --parsable --job-name=pol2_m6a_sum --dependency=afterok:$s "$SELF" summarize)
        p=$(sbatch --parsable --job-name=pol2_plots --mem=32G --dependency=afterok:$z "$SELF" plots)
        echo "fire_fp $f -> fire_fp_summary $g -> states $s -> summarize $z -> plots $p"
    elif [ "$step" = m6a ]; then
        sbatch --job-name=pol2_m6a --array=1-${#CHROMS[@]} "$SELF" m6a
    elif [ "$step" = fire_fp ]; then
        sbatch --job-name=pol2_fire_fp --array=1-${#SAMPLES[@]} "$SELF" fire_fp
    else
        sbatch --job-name=pol2_$step "$SELF" "$step"
    fi
    exit 0
fi

module load R/4.4.1
case "$1" in
    fire_fp)         "$PYTHON" "$DIR/pol2_consensus_fire_footprints.py" --sample "${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}" ;;
    fire_fp_summary) "$PYTHON" "$DIR/pol2_consensus_fire_footprints.py" --summarize ;;
    states)          Rscript "$SCRIPT" states ;;
    m6a)             Rscript "$SCRIPT" m6a "${CHROMS[$((SLURM_ARRAY_TASK_ID - 1))]}" ;;
    summarize)       Rscript "$SCRIPT" summarize ;;
    plots)           Rscript --vanilla "$DIR/plot_intergenic_pol2_examples.R" ;;
    *)               echo "unknown step: $1" >&2; exit 1 ;;
esac
