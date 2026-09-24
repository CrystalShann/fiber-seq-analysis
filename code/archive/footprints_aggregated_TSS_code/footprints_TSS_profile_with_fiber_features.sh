#!/bin/bash
#SBATCH --job-name=ft_FIRE_tss
#SBATCH --time=24:00:00
#SBATCH --mem=120G
#SBATCH --cpus-per-task=3
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/ft_FIRE_tss%A_%a.err

export LANG=C
export LC_ALL=C

module load R/4.4.1
Rscript /project/spott/cshan/fiber-seq/code/footprints_TSS_profiles_with_fiber_features.r sample=AL10_bc2178_19130 chrom=chr1 run_all_chr=false