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

# Runs the code of LCL_co-access_v3_donor_log_reg.rmd without knitting.
#   sbatch LCL_run_v3_donor_log_reg.sh chr21   # smoke test
#   sbatch LCL_run_v3_donor_log_reg.sh         # genome-wide
# Per-chromosome fits are checkpointed under coaccess/donor_adjusted/summary/;
# a later knit with the same R, inputs and code reuses them instead of refitting.
set -euo pipefail
RMD=/project/spott/cshan/fiber-seq/code/co-accessibility/LCL/LCL_co-access_v3_donor_log_reg.rmd
RSCRIPT=/software/R-4.4.1-el8-x86_64/bin/Rscript
"$RSCRIPT" -e '
args <- commandArgs(trailingOnly = TRUE)
code <- parse(knitr::purl(args[1], output = tempfile(fileext = ".R"), documentation = 0, quiet = TRUE))
stopifnot(identical(code[[1]][[2]], as.name("params")))
eval(code[[1]])
params$cores <- as.integer(args[2])
if (length(args) == 3) params$chrom <- args[3]
for (expr in code[-1]) eval(expr)
' "$RMD" "$SLURM_CPUS_PER_TASK" "$@"
