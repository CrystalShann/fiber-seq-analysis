#!/bin/bash
#SBATCH --job-name=tss_features
#SBATCH --time=12:00:00
#SBATCH --mem=1500G
#SBATCH --cpus-per-task=3
#SBATCH --account=pi-spott
#SBATCH --partition=bigmem
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/tss_features%A_%a.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/tss_features%A_%a.err

# Fix locale - use C which is always available
export LANG=C
export LC_ALL=C

module load R/4.4.1
Rscript /project/spott/cshan/fiber-seq/code/expr_stratified_TSS_profiles_with_fiber_features.r sample=AL10_bc2178_19130 chrom=all