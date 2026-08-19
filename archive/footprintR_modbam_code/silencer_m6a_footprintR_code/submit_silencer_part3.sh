#!/bin/bash
#SBATCH --job-name=silencer_part3
#SBATCH --time=24:00:00
#SBATCH --mem=180G
#SBATCH --cpus-per-task=2
#SBATCH --account=pi-spott
#SBATCH --partition=bigmem
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/silencer_part3_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/silencer_part3_%j.err

export LANG=C
export LC_ALL=C

export SAMPLE="${SAMPLE:-AL10_bc2178_19130}"
export CHROM="${CHROM:-chr1}"
export WINDOW_BP="${WINDOW_BP:-300}"
export SILENCER_CLASS="${SILENCER_CLASS:-ALL}"
export OUTPUT_SCOPE="${OUTPUT_SCOPE:-sample_chr}"

module load R/4.4.1

if ! command -v Rscript >/dev/null 2>&1; then
	echo "ERROR: Rscript not found after loading R/4.4.1. Check module environment." >&2
	exit 127
fi

Rscript /project/spott/cshan/fiber-seq/code/footprintR_modbam_code/silencer_m6a_footprintR_code/silencer_modbam_summary_part3.R