#!/bin/bash
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/nucleosome_entropy_%j.log
# Run locally: bash code/entropy/nucleosome/run_entropy.sh
# Slurm: sbatch --account=pi-spott --partition=caslake --cpus-per-task=4 \
#   --mem=8G --time=01:00:00 code/entropy/nucleosome/run_entropy.sh
# Optional trailing gene names select a subset; ENTROPY_OUT sets its output folder.
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH="/software/gcc-12.2.0-el8-x86_64/lib64:/software/openblas-0.3.29-el8-x86_64/lib:/software/glpk-5.0-el8-x86_64/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export R_LIBS_USER="/project/spott/cshan/Rlibs/x86_64-pc-linux-gnu-library/4.4"
/software/R-4.4.1-el8-x86_64/bin/Rscript --vanilla \
  /project/spott/cshan/fiber-seq/code/entropy/nucleosome/nucleosome_entropy.R "$@"
