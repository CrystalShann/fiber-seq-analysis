#!/bin/bash
#SBATCH --job-name=fp_summary_bed
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=04:00:00
#SBATCH --array=1-24
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/fp_summary_bed_%A_%a.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/fp_summary_bed_%A_%a.err

module load python/miniforge-25.3.0
module load htslib
source activate /project/spott/cshan/envs/Jupyter-notebook

CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
        chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 \
        chr21 chr22 chrX chrY)
CHROM=${CHROMS[$((SLURM_ARRAY_TASK_ID - 1))]}

SCRIPT="/project/spott/cshan/fiber-seq/code/PolII_footprints_code/footprints_summary_code/index_fp_summaries_as_bed.py"
INPUT_DIR="/project/spott/cshan/fiber-seq/results/PolII/footprint_summaries"
OUT_DIR="/project/spott/cshan/fiber-seq/results/PolII/footprint_summary_beds"
TMP_DIR="/tmp/${USER}/fp_summary_beds_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

mkdir -p "${TMP_DIR}"

echo "Array task : ${SLURM_ARRAY_TASK_ID}"
echo "Chromosome : ${CHROM}"
echo "Input dir  : ${INPUT_DIR}"
echo "Output dir : ${OUT_DIR}"
echo "Temp dir   : ${TMP_DIR}"
echo "Started    : $(date)"

python3 "${SCRIPT}" \
    --input-dir "${INPUT_DIR}" \
    --out-dir "${OUT_DIR}" \
    --tmp-dir "${TMP_DIR}" \
    --chroms "${CHROM}" \
    --classes all \
    --force

rm -rf "${TMP_DIR}"

echo "Finished   : $(date)"
