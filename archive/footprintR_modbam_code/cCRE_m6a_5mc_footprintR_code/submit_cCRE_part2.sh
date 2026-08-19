#!/bin/bash
#SBATCH --job-name=cCRE_part2
#SBATCH --time=24:00:00
#SBATCH --mem=400G
#SBATCH --cpus-per-task=3
#SBATCH --account=pi-spott
#SBATCH --partition=bigmem
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/cCRE_part2_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/cCRE_part2_%j.err

export LANG=C
export LC_ALL=C

export SAMPLE="${SAMPLE:-AL10_bc2178_19130}"
export CHROM="${CHROM:-chr1}"
export WINDOW_BP="${WINDOW_BP:-1000}"
export BIN_SIZE="${BIN_SIZE:-5}"
export MOD_PROB_THRESHOLD="${MOD_PROB_THRESHOLD:-0.9}"
export CCRE_CHUNK_SIZE="${CCRE_CHUNK_SIZE:-500}"

module load R/4.4.1

Rscript /project/spott/cshan/fiber-seq/code/footprintR_modbam_code/cCRE_m6a_5mc_code/cCRE_modbam_summary_part2.R
