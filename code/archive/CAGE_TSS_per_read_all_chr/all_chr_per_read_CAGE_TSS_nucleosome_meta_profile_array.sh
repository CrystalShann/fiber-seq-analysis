#!/bin/bash
#SBATCH --job-name=perread_cage_chr
#SBATCH --time=12:00:00
#SBATCH --mem=180G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --array=1-22
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/perread_cage_chr%A_%a.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/perread_cage_chr%A_%a.err

set -euo pipefail

export LANG=C
export LC_ALL=C

module load R/4.4.1

SAMPLE="${SAMPLE:-AL10_bc2178_19130}"
CHROM="chr${SLURM_ARRAY_TASK_ID}"

Rscript /project/spott/cshan/fiber-seq/all_chr_per_read_CAGE_TSS_nucleosome_meta_profile.r \
  sample="${SAMPLE}" mode=chunk chrom="${CHROM}"

# Submit combine after array completes, e.g.:
# jid=$(sbatch --parsable /project/spott/cshan/fiber-seq/code/all_chr_per_read_CAGE_TSS_nucleosome_meta_profile_array.sh)
# sbatch --dependency=afterok:${jid} /project/spott/cshan/fiber-seq/code/all_chr_per_read_CAGE_TSS_nucleosome_meta_profile.sh
