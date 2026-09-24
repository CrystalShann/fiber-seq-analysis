#!/bin/bash
#SBATCH --job-name=tss_autocorrelation
#SBATCH --account=pi-spott
#SBATCH --partition=bigmem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --time=30:00:00
#SBATCH --mem=400G
#SBATCH --chdir=/project/spott/cshan/fiber-seq/code/clustering_methods/auto_correlation/tss_test
#SBATCH --output=/project/spott/cshan/fiber-seq/LCL_project/auto_correlation/tss_test/tss_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/LCL_project/auto_correlation/tss_test/tss_%j.err


# Usage: sbatch run_tss.slurm.sh [run_tss.py arguments]


set -euo pipefail

TSS_WORKFLOW_DIR=/project/spott/cshan/fiber-seq/code/clustering_methods/auto_correlation/tss_test

TSS_PYTHON=${TSS_PYTHON:-/project/spott/cshan/envs/Jupyter-notebook/bin/python}
export AUTOCOR_RSCRIPT=${AUTOCOR_RSCRIPT:-/software/R-4.4.1-el8-x86_64/bin/Rscript}
export R_LIBS_USER=${R_LIBS_USER:-/home/cshan/R/x86_64-pc-linux-gnu-library/4.4}
export LD_LIBRARY_PATH="/software/openblas-0.3.29-el8-x86_64/lib:/software/glpk-5.0-el8-x86_64/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONUNBUFFERED=1

if [[ ! -x "$TSS_PYTHON" || ! -x "$AUTOCOR_RSCRIPT" ]]; then
    printf 'Missing executable: Python=%s; Rscript=%s\n' "$TSS_PYTHON" "$AUTOCOR_RSCRIPT" >&2
    exit 1
fi

printf 'Job: %s; host: %s; UTC start: %s\n' "${SLURM_JOB_ID:-interactive}" "$(hostname)" "$(date -u +%FT%TZ)"
printf 'Command: '
printf '%q ' "$TSS_PYTHON" -u -B "$TSS_WORKFLOW_DIR/run_tss.py" "$@"
printf '\n'

exec "$TSS_PYTHON" -u -B "$TSS_WORKFLOW_DIR/run_tss.py" "$@"
