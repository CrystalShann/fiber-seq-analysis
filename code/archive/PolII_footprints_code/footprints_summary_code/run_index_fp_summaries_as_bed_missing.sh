#!/bin/bash
#SBATCH --job-name=fp_summary_bed_missing
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=04:00:00
#SBATCH --array=1-9
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/fp_summary_bed_missing_%A_%a.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/fp_summary_bed_missing_%A_%a.err

module load python/miniforge-25.3.0
module load htslib
source activate /project/spott/cshan/envs/Jupyter-notebook

# (class, chrom) pairs that are missing from
# /project/spott/cshan/fiber-seq/results/PolII/footprint_summary_beds/
PAIRS=(
    "FIRE_nucleosome chr5"
    "FIRE_nucleosome chrX"
    "unknown chr5"
    "unknown chr6"
    "unknown chr10"
    "unknown chr12"
    "unknown chr13"
    "unknown chr14"
    "unknown chrX"
)
PAIR=${PAIRS[$((SLURM_ARRAY_TASK_ID - 1))]}
CLASS=${PAIR%% *}
CHROM=${PAIR##* }

SCRIPT="/project/spott/cshan/fiber-seq/code/PolII_footprints_code/footprints_summary_code/index_fp_summaries_as_bed.py"
INPUT_DIR="/project/spott/cshan/fiber-seq/results/PolII/footprint_summaries"
OUT_DIR="/project/spott/cshan/fiber-seq/results/PolII/footprint_summary_beds"
TMP_DIR="/tmp/${USER}/fp_summary_beds_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

mkdir -p "${TMP_DIR}"

echo "Array task : ${SLURM_ARRAY_TASK_ID}"
echo "Class      : ${CLASS}"
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
    --classes "${CLASS}" \
    --force

rm -rf "${TMP_DIR}"

echo "Finished   : $(date)"
