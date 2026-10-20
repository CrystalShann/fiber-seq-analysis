#!/bin/bash
#SBATCH --job-name=LCL_render
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/LCL_render_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/LCL_render_%j.err

# Renders one notebook in this folder to /project/spott/cshan/fiber-seq/docs/<name>.html.
#   sbatch LCL_render_rmd.sh LCL_co-access.Rmd                       # after counts + LCL_coaccess_noXCI.py
#   sbatch LCL_render_rmd.sh LCL_co-access_v2_perm_methy_prop.Rmd    # after LCL_co-access.Rmd (reads its top-ten tables)
#   sbatch LCL_render_rmd.sh LCL_co-access_v3_donor_log_reg.rmd      # after LCL_run_v3_donor_log_reg.sh
set -euo pipefail
module load R/4.4.1
module load pandoc/2.17.1.1
Rscript -e '
args <- commandArgs(trailingOnly = TRUE)
tmp <- tempfile("render_"); dir.create(tmp)
rmarkdown::render(args[1], output_dir = args[2], intermediates_dir = tmp, envir = new.env(), quiet = TRUE)
' "/project/spott/cshan/fiber-seq/code/co-accessibility/LCL/$1" /project/spott/cshan/fiber-seq/docs
