#!/usr/bin/env bash
#SBATCH --job-name=tss_acf
#SBATCH --account=pi-spott
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=300G
#SBATCH --time=30:00:00
#SBATCH --array=0-2
#SBATCH --output=/project/spott/cshan/fiber-seq/macrophage_project/auto_correlation/tss/logs/slurm_%A_%a.out
#SBATCH --error=/project/spott/cshan/fiber-seq/macrophage_project/auto_correlation/tss/logs/slurm_%A_%a.err
# Array task -> TSS window: 0 = 2000_tss, 1 = left_1000_tss, 2 = right_1000_tss.
# Run one window only with e.g. `sbatch --array=2 06_run_tss_autocorrelation.sh`.

set -euo pipefail
PROJECT=/project/spott/cshan/fiber-seq
CODE="$PROJECT/code/clustering_methods/auto_correlation/tss"
OUT="$PROJECT/macrophage_project/auto_correlation/tss"
PYTHON=${TSS_PYTHON:-/project/spott/cshan/envs/Jupyter-notebook/bin/python}
RSCRIPT=${TSS_RSCRIPT:-/software/R-4.4.1-el8-x86_64/bin/Rscript}
BINS="$PROJECT/macrophage_project/expr_access/tables/tss_expression_bins.tsv"
CANONICAL=/project/spott/cshan/annotations/gencodev46_Ensembl_canonical_TSS.bed
FT_ROOT="$PROJECT/macrophage_project/FiberHMM/extract/ft_result_dir"
PER_BIN=2500
SEED=0
MAX_LAG=""   # empty = all nonnegative lags of the window
N_PCS=50
N_NEIGHBORS=10
RESOLUTION=0.4
CHROMS=()
BINS_EXPLICIT=0
WINDOWS=(2000_tss left_1000_tss right_1000_tss)
OUT_EXPLICIT=""
TEST=0
while (($#)); do
  case "$1" in
    --test) TEST=1; PER_BIN=100; CHROMS=(chr22); shift ;;
    --out-dir) OUT_EXPLICIT="$2"; shift 2 ;;
    --bins-tsv) BINS="$2"; BINS_EXPLICIT=1; shift 2 ;;
    --canonical-bed) CANONICAL="$2"; shift 2 ;;
    --ft-root) FT_ROOT="$2"; shift 2 ;;
    --per-bin) PER_BIN="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --max-lag) MAX_LAG="$2"; shift 2 ;;
    --n-pcs) N_PCS="$2"; shift 2 ;;
    --n-neighbors) N_NEIGHBORS="$2"; shift 2 ;;
    --resolution) RESOLUTION="$2"; shift 2 ;;
    --chrom) IFS=, read -r -a CHROMS <<< "$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
  echo "Submit this pipeline with sbatch; heavy computation requires a SLURM allocation." >&2
  exit 2
fi
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" || -z "${WINDOWS[$SLURM_ARRAY_TASK_ID]:-}" ]]; then
  echo "SLURM_ARRAY_TASK_ID must be 0 (2000_tss), 1 (left_1000_tss) or 2 (right_1000_tss)" >&2
  exit 2
fi
WINDOW=${WINDOWS[$SLURM_ARRAY_TASK_ID]}
# Genomic offsets from the TSS (end exclusive); each window has its own output folder.
case "$WINDOW" in
  2000_tss) WINDOW_START=-1000; WINDOW_END=1000 ;;
  left_1000_tss) WINDOW_START=-1000; WINDOW_END=-100 ;;
  right_1000_tss) WINDOW_START=100; WINDOW_END=1000 ;;
esac
OUT="${OUT_EXPLICIT:-$OUT/$WINDOW}"
if ((TEST)); then OUT="$OUT/test"; fi
mkdir -p "$OUT"
OUT=$(cd -- "$OUT" && pwd -P)
# Lock the directory inode, without creating a persistent lock file.
exec 9<"$OUT"
flock -n 9 || { echo "Another pipeline is using $OUT" >&2; exit 1; }
# Respect an older runner still holding its existing lock; never create one.
if [[ -f "$OUT/.pipeline.lock" ]]; then
  exec 8<"$OUT/.pipeline.lock"
  flock -n 8 || { echo "An older pipeline is still using $OUT" >&2; exit 1; }
fi

# Worker scripts share an ephemeral workspace. Only final tables and plots
# are published after all computation and plotting checks have passed.
TSS_SCRATCH_ROOT=${SLURM_TMPDIR:-${TMPDIR:-/tmp}}
TSS_WORK_DIR=$(mktemp -d "${TSS_SCRATCH_ROOT%/}/tss_acf.${SLURM_JOB_ID}.XXXXXX")
TSS_PUBLISH_TEMP=""
cleanup() {
  local status=$?
  local cleanup_failed=0
  trap - EXIT
  if [[ -n "$TSS_PUBLISH_TEMP" ]] && ! rm -f -- "$TSS_PUBLISH_TEMP"; then cleanup_failed=1; fi
  if [[ -d "$TSS_WORK_DIR" ]] && ! rm -rf -- "$TSS_WORK_DIR"; then cleanup_failed=1; fi
  if ((cleanup_failed)); then
    echo "Could not completely remove temporary files for $TSS_WORK_DIR" >&2
    if ((status == 0)); then status=1; fi
  fi
  exit "$status"
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
trap 'echo "[$(date -Is)] FAILED at line $LINENO" >&2' ERR
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export NUMBA_NUM_THREADS=1 R_DATATABLE_NUM_THREADS=1 PYTHONHASHSEED="$SEED"
export PYTHONDONTWRITEBYTECODE=1
# Package caches and other temporary files share the same disposable workspace.
export NUMBA_CACHE_DIR="$TSS_WORK_DIR/intermediate/numba_cache"
export MPLCONFIGDIR="$TSS_WORK_DIR/intermediate/matplotlib"
mkdir -p "$TSS_WORK_DIR/tmp"
export TMPDIR="$TSS_WORK_DIR/tmp"
export R_LIBS_USER="${R_LIBS_USER:-/project/spott/cshan/Rlibs/x86_64-pc-linux-gnu-library/4.4}"
export LD_LIBRARY_PATH="/software/openblas-0.3.29-el8-x86_64/lib:/software/glpk-5.0-el8-x86_64/lib:${LD_LIBRARY_PATH:-}"
echo "[$(date -Is)] Start job $SLURM_JOB_ID task $SLURM_ARRAY_TASK_ID on $(hostname); output=$OUT; per_bin=$PER_BIN; window=$WINDOW [TSS${WINDOW_START}, TSS+${WINDOW_END})"
echo "[$(date -Is)] Temporary workspace: $TSS_WORK_DIR (removed on exit)"
"$PYTHON" -B -c 'import sys,numpy,pandas,pysam,scanpy,igraph,leidenalg; print(sys.version); print("scanpy",scanpy.__version__,"numpy",numpy.__version__,"pysam",pysam.__version__)'
"$RSCRIPT" --vanilla -e 'cat(R.version.string,"\n")'
if [[ ! -f "$BINS" ]]; then
  if ((BINS_EXPLICIT)); then echo "Explicit expression table is missing: $BINS" >&2; exit 1; fi
  echo "[$(date -Is)] Rebuilding absent expression table using unchanged original R code"
  "$RSCRIPT" --vanilla "$CODE/00_prepare_expression_bins.R" --out-dir "$TSS_WORK_DIR/inputs"
  BINS="$TSS_WORK_DIR/inputs/tss_expression_bins.tsv"
fi
EXTRA=()
if ((${#CHROMS[@]})); then EXTRA=(--chrom "${CHROMS[@]}"); fi
LAG_ARGS=()
if [[ -n "$MAX_LAG" ]]; then LAG_ARGS=(--max-lag "$MAX_LAG"); fi
echo "[$(date -Is)] 01: survey, deduplicate, balanced sample, then binary matrix"
"$PYTHON" -B "$CODE/01_sample_tss_molecules.py" --out-dir "$TSS_WORK_DIR" --bins-tsv "$BINS" \
  --canonical-bed "$CANONICAL" --ft-root "$FT_ROOT" --per-bin "$PER_BIN" --seed "$SEED" \
  --window-start "$WINDOW_START" --window-end "$WINDOW_END" "${EXTRA[@]}"
echo "[$(date -Is)] 02: unchanged parent ACF"
"$PYTHON" -B "$CODE/02_compute_tss_autocorrelations.py" --out-dir "$TSS_WORK_DIR" "${LAG_ARGS[@]}"
echo "[$(date -Is)] 03: joint parent PCA / correlation kNN / Leiden / UMAP"
"$PYTHON" -B "$CODE/03_leiden_cluster_autocorrelations.py" --out-dir "$TSS_WORK_DIR" \
  --n-pcs "$N_PCS" --n-neighbors "$N_NEIGHBORS" --resolution "$RESOLUTION" --seed "$SEED"
echo "[$(date -Is)] 04: plotting tables and descriptive summaries"
"$PYTHON" -B "$CODE/04_prepare_plot_tables.py" --out-dir "$TSS_WORK_DIR"
echo "[$(date -Is)] 05: final R plots"
"$RSCRIPT" --vanilla "$CODE/05_plot_tss_autocorrelation.R" --out-dir "$TSS_WORK_DIR"

echo "[$(date -Is)] Publishing final tables and plots"
shopt -s nullglob
TABLE_FILES=("$TSS_WORK_DIR"/tables/*.tsv "$TSS_WORK_DIR"/tables/*.tsv.gz)
PLOT_FILES=("$TSS_WORK_DIR"/plots/*.pdf "$TSS_WORK_DIR"/plots/*.png)
if ((${#TABLE_FILES[@]} < 10 || ${#PLOT_FILES[@]} < 11)); then
  echo "Expected at least 10 final tables and 11 plots; publishing aborted" >&2
  exit 1
fi
for source in "${TABLE_FILES[@]}" "${PLOT_FILES[@]}"; do
  [[ -s "$source" ]] || { echo "Empty result: $source" >&2; exit 1; }
done
publish_file() {
  local source="$1" destination="$2" name
  name=$(basename -- "$source")
  mkdir -p "$destination"
  TSS_PUBLISH_TEMP=$(mktemp "$destination/.${name}.${SLURM_JOB_ID}.XXXXXX")
  cp --preserve=mode -- "$source" "$TSS_PUBLISH_TEMP"
  mv -f -- "$TSS_PUBLISH_TEMP" "$destination/$name"
  TSS_PUBLISH_TEMP=""
}
for source in "${TABLE_FILES[@]}"; do publish_file "$source" "$OUT/tables"; done
for source in "${PLOT_FILES[@]}"; do publish_file "$source" "$OUT/plots"; done
if [[ -f "$TSS_WORK_DIR/inputs/tss_expression_bins.tsv" ]]; then
  publish_file "$TSS_WORK_DIR/inputs/tss_expression_bins.tsv" "$OUT/tables"
fi
echo "[$(date -Is)] COMPLETE: $OUT"
