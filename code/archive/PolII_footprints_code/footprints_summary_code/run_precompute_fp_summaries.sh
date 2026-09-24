#!/bin/bash
#SBATCH --job-name=fp_summaries
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=06:00:00
#SBATCH --array=1-24
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/fp_summaries_%A_%a.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/fp_summaries_%A_%a.err

module load python/miniforge-25.3.0
source activate /project/spott/cshan/envs/Jupyter-notebook

# Map SLURM array index (1-24) to chromosome name
CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
        chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 \
        chr21 chr22 chrX chrY)
CHROM=${CHROMS[$((SLURM_ARRAY_TASK_ID - 1))]}

SCRIPT="/project/spott/cshan/fiber-seq/code/PolII_footprints_code/precompute_fp_summaries.py"
OUT_DIR="/project/spott/cshan/fiber-seq/results/PolII/footprint_summaries"

echo "Array task : ${SLURM_ARRAY_TASK_ID}"
echo "Chromosome : ${CHROM}"
echo "Output dir : ${OUT_DIR}"
echo "Started    : $(date)"

python3 "${SCRIPT}" \
    --chroms "${CHROM}" \
    --out-dir "${OUT_DIR}"

echo "Finished   : $(date)"
