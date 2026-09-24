#!/bin/bash
#SBATCH --job-name=modbam_m6a
#SBATCH --time=12:00:00
#SBATCH --mem=400G
#SBATCH --cpus-per-task=3
#SBATCH --account=pi-spott
#SBATCH --partition=bigmem
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/modbam_m6a_%j.err


export LANG=C
export LC_ALL=C

module load R/4.4.1
Rscript /project/spott/cshan/fiber-seq/code/footprintR_modbam_code/promoter/TSS_nuc_m6a_5mc_code/run_modbam_promoter_summary_m6A.R