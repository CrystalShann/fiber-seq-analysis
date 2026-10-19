#!/bin/bash
#SBATCH --job-name=LCL_tf_occupancy
#SBATCH --account=pi-spott
#SBATCH --partition=bigmem
#SBATCH --qos=bigmem
#SBATCH --cpus-per-task=2
#SBATCH --mem=400G
#SBATCH --time=24:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/LCL_tf_%A_%a.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/LCL_tf_%A_%a.err
set -euo pipefail
SCRIPT_DIR=${LCL_TF_CODE:-/project/spott/cshan/fiber-seq/code/co-accessibility/LCL/tf_co-occ}
PYTHON=${LCL_TF_PYTHON:-/project/spott/cshan/envs/Jupyter-notebook/bin/python3}
ROOT=${LCL_COACCESS_ROOT:-/project/spott/cshan/fiber-seq/LCL_project/co-accessibility}
OUT=${LCL_TF_OUT:-$ROOT/tf_co-occ}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
stage=${1:-help}
shift || true
case "$stage" in
  prepare|combine)
    exec "$PYTHON" "$SCRIPT_DIR/LCL_tf_cooccupancy.py" "$stage" --root "$ROOT" --out "$OUT" "$@" ;;
  count)
    exec "$PYTHON" "$SCRIPT_DIR/LCL_tf_cooccupancy.py" count --root "$ROOT" --out "$OUT" --sample-index "${SLURM_ARRAY_TASK_ID:-1}" "$@" ;;
  analyze|report)
    module load R/4.4.1
    module load pandoc/2.17.1.1
    Rscript - "$OUT" "$SCRIPT_DIR/LCL_tf_cooccupancy.Rmd" "$stage" <<'RSCRIPT'
main <- function() {
  args <- commandArgs(TRUE)
  out <- normalizePath(args[1],mustWork=TRUE)
  source <- normalizePath(args[2],mustWork=TRUE)
  run_models <- identical(args[3],'analyze')
  temp <- tempfile('lcl_tf_report_');dir.create(temp)
  on.exit(unlink(temp,recursive=TRUE),add=TRUE)
  local_source <- file.path(temp,basename(source))
  stopifnot(file.copy(source,local_source))
  rmarkdown::render(local_source,output_format='pdf_document',
    output_options=list(latex_engine='xelatex',keep_tex=FALSE),
    params=list(out=out,run_models=run_models),output_dir=temp,intermediates_dir=temp,
    clean=TRUE,envir=new.env(),quiet=TRUE)
  pdf <- file.path(temp,'LCL_tf_cooccupancy.pdf')
  stopifnot(file.exists(pdf),file.copy(pdf,file.path(out,basename(pdf)),overwrite=TRUE))
  # knitr may clean its figure directory during render. The report always embeds
  # vector plots; copy standalone figures when retained by the rendering backend.
  figures <- setdiff(list.files(temp,pattern='\\.pdf$',recursive=TRUE,full.names=TRUE),pdf)
  if(length(figures)) {
    dir.create(file.path(out,'plots'),showWarnings=FALSE)
    stopifnot(all(file.copy(figures,file.path(out,'plots',basename(figures)),overwrite=TRUE)))
  }
}
main()
RSCRIPT
    ;;
  submit)
    mkdir -p /project/spott/cshan/fiber-seq/results/logs
    prep=$(sbatch --parsable --partition=caslake --qos=caslake --mem=48G "$SCRIPT_DIR/LCL_tf_cooccupancy.sh" prepare)
    fimo=$(sbatch --parsable --dependency="afterok:$prep" "$SCRIPT_DIR/LCL_tf_fimo_scan.sh")
    counts=$(sbatch --parsable --partition=caslake --qos=caslake --mem=48G --dependency="afterok:$fimo" --array=1-31%6 "$SCRIPT_DIR/LCL_tf_cooccupancy.sh" count)
    combine=$(sbatch --parsable --partition=caslake --qos=caslake --mem=48G --dependency="afterok:$counts" "$SCRIPT_DIR/LCL_tf_cooccupancy.sh" combine)
    analysis=$(sbatch --parsable --dependency="afterok:$combine" "$SCRIPT_DIR/LCL_tf_cooccupancy.sh" analyze)
    printf 'prepare=%s fimo=%s counts=%s combine=%s analysis=%s\n' "$prep" "$fimo" "$counts" "$combine" "$analysis" ;;
  *) echo 'Usage: bash LCL_tf_cooccupancy.sh {prepare|count|combine|analyze|report|submit} [arguments]' >&2; exit 2 ;;
esac
