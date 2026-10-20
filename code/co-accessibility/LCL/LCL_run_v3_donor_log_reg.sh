#!/bin/bash
#SBATCH --job-name=LCL_v3_donor
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/LCL_v3_donor_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/LCL_v3_donor_%j.err

# Fits the donor-adjusted models for LCL_co-access_v3_donor_log_reg.rmd by running
# LCL_fit_v3_donor_log_reg.R; knitting the Rmd only reads the tables this job writes.
#   sbatch LCL_run_v3_donor_log_reg.sh chr21   # smoke test
#   sbatch LCL_run_v3_donor_log_reg.sh         # genome-wide
# Per-chromosome fits are checkpointed under coaccess/donor_adjusted/summary/;
# a rerun with the same R, inputs and code reuses them instead of refitting.
set -euo pipefail
/software/R-4.4.1-el8-x86_64/bin/Rscript \
  /project/spott/cshan/fiber-seq/code/co-accessibility/LCL/LCL_fit_v3_donor_log_reg.R \
  "$SLURM_CPUS_PER_TASK" "$@"
