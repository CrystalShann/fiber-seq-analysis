#!/bin/bash
#SBATCH --job-name=pol2_clusters
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/pol2_clusters_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/pol2_clusters_%j.err
set -euo pipefail
export LD_LIBRARY_PATH="/project/spott/cshan/envs/dimelo/lib:/software/gcc-12.2.0-el8-x86_64/lib64:/software/openblas-0.3.29-el8-x86_64/lib:/software/glpk-5.0-el8-x86_64/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export R_LIBS_USER=/project/spott/cshan/Rlibs/x86_64-pc-linux-gnu-library/4.4
export RSTUDIO_PANDOC=/software/pandoc-2.17.1.1-el8-x86_64/bin
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
script_dir=/project/spott/cshan/fiber-seq/code/pausing/clustering_Pol2
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  # Extraction only. Chromosome chunks are pooled into one subsequent clustering.
  if (( SLURM_ARRAY_TASK_ID == 23 )); then
    chrom=chrX
  else
    chrom="chr${SLURM_ARRAY_TASK_ID}"
  fi
  exec /software/R-4.4.1-el8-x86_64/bin/Rscript --vanilla "$script_dir/run_pol2_clustering.R" \
    --stage features --chrom "$chrom" "$@"
fi
exec /software/R-4.4.1-el8-x86_64/bin/Rscript --vanilla "$script_dir/run_pol2_clustering.R" "$@"
