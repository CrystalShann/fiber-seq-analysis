#!/bin/bash
#SBATCH --job-name=fire_features_allchr
#SBATCH --time=12:00:00
#SBATCH --mem=180G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/fire_features_allchr_%j.err
#SBATCH --out=/project/spott/cshan/fiber-seq/results/logs/fire_features_allchr_%j.out


export LANG=C
export LC_ALL=C

module load R/4.4.1

SAMPLE="${SAMPLE:-AL10_bc2178_19130}"
MODE="${MODE:-combine}"
CHROM="${CHROM:-}"

if [[ "${MODE}" == "chunk" && -z "${CHROM}" ]]; then
  echo "MODE=chunk requires CHROM to be set" >&2
  exit 1
fi

if [[ "${MODE}" == "chunk" ]]; then
  Rscript /project/spott/cshan/fiber-seq/code/all_chr_aggregated_FIRE_features.r sample="${SAMPLE}" mode="${MODE}" chrom="${CHROM}"
else
  Rscript /project/spott/cshan/fiber-seq/code/all_chr_aggregated_FIRE_features.r sample="${SAMPLE}" mode="${MODE}"
fi

Rscript /project/spott/cshan/fiber-seq/code/all_chr_aggregated_FIRE_features.r sample=AL10_bc2178_19130 mode=combine

# Rscript /project/spott/cshan/fiber-seq/code/all_chr_aggregated_FIRE_features.r sample=AL10_bc2178_19130 mode=chunk chrom=chr1
# sbatch /project/spott/cshan/fiber-seq/code/all_chr_aggregated_FIRE_features.sh
