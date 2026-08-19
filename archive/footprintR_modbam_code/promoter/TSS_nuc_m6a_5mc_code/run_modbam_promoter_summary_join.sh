#!/bin/bash
#SBATCH --job-name=modbam_join
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --account=pi-spott
#SBATCH --partition=bigmem
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/modbam_join_%j.err

export LANG=C
export LC_ALL=C

module load R/4.4.1
Rscript /project/spott/cshan/fiber-seq/code/footprintR_modbam_code/run_modbam_promoter_summary_join.R
