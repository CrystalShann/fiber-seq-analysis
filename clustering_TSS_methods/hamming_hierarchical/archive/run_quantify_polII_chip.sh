#!/bin/bash
#SBATCH --job-name=quantify_polII_chip
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/quantify_polII_chip_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/quantify_polII_chip_%j.err

# Lifts the GSE32324 THP-1 RNAPII ChIP read coordinates hg19 -> hg38 and
# counts overlaps with the promoter/pause/body windows from
# promoter_regions_gencode_v49.tsv (see quantify_polII_chip.R).
#
#   sbatch run_quantify_polII_chip.sh
#
# Requires make_promoter_table.R to have been run first.

set -uo pipefail

module load R/4.4.1

Rscript /project/spott/cshan/fiber-seq/code/clustering_TSS_methods/hamming_hierarchical/quantify_polII_chip.R
