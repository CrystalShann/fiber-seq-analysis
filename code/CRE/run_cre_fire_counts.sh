#!/bin/bash
# Count missing/outdated chromosomes and rebuild the combined FIRE-frequency tables.
# Usage: sbatch code/CRE/run_cre_fire_counts.sh [--test]
# Create results/logs before submitting.
#SBATCH --job-name=cre_fire_counts
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/cre_fire_counts_%j.log

set -euo pipefail
cd /project/spott/cshan/fiber-seq
export LC_ALL=C
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export LD_LIBRARY_PATH="/project/spott/cshan/envs/dimelo/lib:/software/gcc-12.2.0-el8-x86_64/lib64:/software/openblas-0.3.29-el8-x86_64/lib:/software/glpk-5.0-el8-x86_64/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export R_LIBS_USER=/project/spott/cshan/Rlibs/x86_64-pc-linux-gnu-library/4.4

exec /software/R-4.4.1-el8-x86_64/bin/Rscript --vanilla \
  code/CRE/count_cre_fire_frequency.R "$@"
