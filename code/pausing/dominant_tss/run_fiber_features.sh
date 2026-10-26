#!/bin/bash
#SBATCH --job-name=dominant_fibers
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=1-23%4
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/dominant_fibers_%A_%a.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/dominant_fibers_%A_%a.err
# sbatch run_fiber_features.sh [options]; or bash run_fiber_features.sh --chrom chr22 ...
set -euo pipefail
export LD_LIBRARY_PATH="/project/spott/cshan/envs/dimelo/lib:/software/gcc-12.2.0-el8-x86_64/lib64:/software/openblas-0.3.29-el8-x86_64/lib:/software/glpk-5.0-el8-x86_64/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export R_LIBS_USER="/project/spott/cshan/Rlibs/x86_64-pc-linux-gnu-library/4.4"
chrom_args=()
if [[ -n ${SLURM_ARRAY_TASK_ID:-} ]]; then
    [[ $SLURM_ARRAY_TASK_ID -ge 1 && $SLURM_ARRAY_TASK_ID -le 23 ]]
    chrom="chr${SLURM_ARRAY_TASK_ID}"
    if [[ $SLURM_ARRAY_TASK_ID -eq 23 ]]; then chrom=chrX; fi
    chrom_args=(--chrom "$chrom")
fi
/software/R-4.4.1-el8-x86_64/bin/Rscript --vanilla \
    /project/spott/cshan/fiber-seq/code/pausing/dominant_tss/03_fiber_features.R \
    "${chrom_args[@]}" "$@"
