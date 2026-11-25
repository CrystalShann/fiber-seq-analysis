#!/bin/bash
#SBATCH --job-name=GM12878_PINTS
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/LCL_project/divergent_tss_enhancer/GM12878_GROcap_PINTS_outputs/pints_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/LCL_project/divergent_tss_enhancer/GM12878_GROcap_PINTS_outputs/pints_%j.err
set -euo pipefail

PROJECT=/project/spott/cshan/fiber-seq
CODE="$PROJECT/code/divergent_tss_enhancer"
OUT="$PROJECT/LCL_project/divergent_tss_enhancer/GM12878_GROcap_PINTS_outputs"
PINTS=/project/spott/cshan/tools/PINTS
PYTHON="$PINTS/.venv/bin/python"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONHASHSEED=0 PYTHONDONTWRITEBYTECODE=1
# Keep all results directly in OUT; do not create metadata, comparisons, or run subfolders.
mkdir -p "$OUT"
PINTS_TMP=$(mktemp -d "${TMPDIR:-/tmp}/gm12878_pints.XXXXXX")
trap 'rm -rf -- "$PINTS_TMP"' EXIT
export NUMBA_CACHE_DIR="$PINTS_TMP/numba_cache"
mkdir -p "$NUMBA_CACHE_DIR"
test "$(git -C "$PINTS" rev-parse HEAD)" = 5196ded843428b81d46e38ce48e3d53eb5795c63
"$PYTHON" "$CODE/get_divergent_tss_pints.py" validate
PLUS=$(awk -F '\t' '$2 == "+" {print $7}' "$CODE/inputs.tsv")
MINUS=$(awk -F '\t' '$2 == "-" {print $7}' "$CODE/inputs.tsv")
command=("$PINTS/.venv/bin/pints_caller"
    --bw-pl "$PLUS" --bw-mn "$MINUS"
    --save-to "$OUT" --file-prefix GM12878_GROcap_GRCh38
    --thread "${SLURM_CPUS_PER_TASK:-1}"
    --fdr-target 0.1 --close-threshold 300 --summit-dist-min 1
    --chromosome-start-with chr --dont-check-updates)
"${command[@]}"
"$PYTHON" "$CODE/get_divergent_tss_pints.py" summarize \
    --native "$OUT/GM12878_GROcap_GRCh38_1_divergent_peaks.bed" \
    --output "$OUT" --plus-bw "$PLUS" --minus-bw "$MINUS"
printf 'Completed: %s\n' "$OUT"
