#!/bin/bash
#SBATCH --job-name=perread_cage_allchr
#SBATCH --time=12:00:00
#SBATCH --mem=180G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/perread_cage_allchr_%j.err
#SBATCH --out=/project/spott/cshan/fiber-seq/results/logs/perread_cage_allchr_%j.out

set -euo pipefail

export LANG=C
export LC_ALL=C

module load R/4.4.1

SAMPLE="${SAMPLE:-AL10_bc2178_19130}"
MODE="${MODE:-combine}"     # chunk | combine | full
CHROM="${CHROM:-}"          # required when MODE=chunk

if [[ "${MODE}" == "chunk" && -z "${CHROM}" ]]; then
  echo "MODE=chunk requires CHROM to be set" >&2
  exit 1
fi

if [[ "${MODE}" == "chunk" ]]; then
  Rscript /project/spott/cshan/fiber-seq/code/CAGE_TSS_per_read_all_chr/all_chr_per_read_CAGE_TSS_nucleosome_meta_profile.r \
    sample="${SAMPLE}" mode="${MODE}" chrom="${CHROM}"
else
  Rscript /project/spott/cshan/fiber-seq/code/CAGE_TSS_per_read_all_chr/all_chr_per_read_CAGE_TSS_nucleosome_meta_profile.r \
    sample="${SAMPLE}" mode="${MODE}"
fi

# Examples:
# SAMPLE=AL10_bc2178_19130 MODE=chunk CHROM=chr1 sbatch /project/spott/cshan/fiber-seq/code/CAGE_TSS_per_read_all_chr/all_chr_per_read_CAGE_TSS_nucleosome_meta_profile.sh
# SAMPLE=AL10_bc2178_19130 MODE=combine sbatch /project/spott/cshan/fiber-seq/code/CAGE_TSS_per_read_all_chr/all_chr_per_read_CAGE_TSS_nucleosome_meta_profile.sh
