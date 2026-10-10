#!/bin/bash
#SBATCH --job-name=LCL_coaccess_plots
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --export=ALL,LC_ALL=C,LANG=C
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/LCL_coaccess_plots_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/LCL_coaccess_plots_%j.err

# Submit with afterok dependency on LCL_run_coaccess.sh.
set -euo pipefail
LCL_SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
if [[ ! -f "$LCL_SCRIPT_DIR/LCL_co-access.Rmd" ]]; then
    LCL_SCRIPT_DIR=/project/spott/cshan/fiber-seq/code/co-accessibility/LCL
fi
export LCL_SCRIPT_DIR
export LCL_COACCESS_ROOT=${LCL_COACCESS_ROOT:-/project/spott/cshan/fiber-seq/LCL_project/co-accessibility}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 TZ=America/Chicago LC_ALL=C LANG=C
export RSTUDIO_PANDOC=/software/pandoc-2.17.1.1-el8-x86_64/bin
export LD_LIBRARY_PATH=/software/gcc-13.2.0-el8-x86_64/lib64:/software/openblas-0.3.29-el8-x86_64/lib:${LD_LIBRARY_PATH:-}
exec /software/R-4.4.1-el8-x86_64/bin/Rscript --vanilla -e '
root <- Sys.getenv("LCL_COACCESS_ROOT")
out <- file.path(root, "coaccess", "plots")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
settings <- list(root = root)
if (nzchar(Sys.getenv("LCL_SAMPLE_METATABLE"))) settings$sample_sheet <- Sys.getenv("LCL_SAMPLE_METATABLE")
rmarkdown::render(file.path(Sys.getenv("LCL_SCRIPT_DIR"), "LCL_co-access.Rmd"),
  output_dir = out, intermediates_dir = out, knit_root_dir = out, params = settings)
'
