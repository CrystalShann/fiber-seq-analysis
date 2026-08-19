#!/bin/bash
#SBATCH --job-name=extract_region
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=1-36
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/extract_region_%A_%a.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/extract_region_%A_%a.err

# Drives extract_region_result_macrophage.sh over the 9 early-response promoter
# regions x 4 LPS timepoints (36 array tasks).
#
#   sbatch run_extract_regions.sh
#
# Task index -> (region, sample):
#   region index = (task - 1) / 4
#   sample index = (task - 1) % 4

set -uo pipefail

SCRIPT="/project/spott/cshan/fiber-seq/code/topic_model/read_region_extraction/extract_region_result_macrophage.sh"
# original:
# OUT_ROOT="/project/spott/cshan/fiber-seq/macrophage_project/topic_model/region_plots/early_repsonse_genes"
OUT_ROOT="/project/spott/cshan/fiber-seq/macrophage_project/topic_model/region_plots/early_repsonse_genes/promoters"

# Same 9 regions as the topic model (code/topic_model/1_m6a_promoter_topic_modelling.Rmd).
# "<gene> <chr> <start> <end>"
# original (explicit genomic intervals):
# REGIONS=(
#     "IL1B      chr2   112836323  112837385"
#     "TNF       chr6   31575265   31575824"
#     "IL1A      chr2   112793233  112794969"
#     "PTGS2     chr1   186679985  186681163"
#     "IL8       chr4   73740006   73741332"
#     "CXCL2     chr4   74098863   74099651"
#     "IL6       chr7   22726376   22727094"
#     "RELB      chr19  45001011   45002417"
#     "rs2836882 chr21  39094468   39094712"
# )
# TSS +/- 1000 bp, matching the TSS-centered regions in
# 1_m6a_promoter_topic_modelling.Rmd (dominant THP-1 CAGE peak: IL1B/TNF/IL8;
# GENCODE v49 TSS fallback: IL1A/PTGS2/CXCL2/IL6/RELB). rs2836882 is a SNP
# locus with no TSS - it keeps its explicit interval.
REGIONS=(
    "IL1B      chr2   112835779  112837779"
    "TNF       chr6   31574568   31576568"
    "IL1A      chr2   112783493  112785493"
    "PTGS2     chr1   186679922  186681922"
    "IL8       chr4   73739569   73741569"
    "CXCL2     chr4   74098196   74100196"
    "IL6       chr7   22724884   22726884"
    "RELB      chr19  45000461   45002461"
    "rs2836882 chr21  39094468   39094712"
)

SAMPLES=(LPS_0 LPS_5 LPS_10 LPS_15)

t=$((SLURM_ARRAY_TASK_ID - 1))
region_idx=$((t / ${#SAMPLES[@]}))
sample_idx=$((t % ${#SAMPLES[@]}))

read -r gene chr start end <<< "${REGIONS[$region_idx]}"
sample_name=${SAMPLES[$sample_idx]}

region="${chr}:${start}-${end}"
# original:
# outname="${gene}_${chr}_${start}_${end}"
outname="${gene}"

echo "Array task : ${SLURM_ARRAY_TASK_ID}"
echo "Gene       : ${gene}"
echo "Region     : ${region}"
echo "Sample     : ${sample_name}"
echo "Out root   : ${OUT_ROOT}"
echo "Started    : $(date)"

bash "${SCRIPT}" "${sample_name}" "${region}" "${outname}" "${OUT_ROOT}"
rc=$?

echo "Finished   : $(date) (exit ${rc})"
exit ${rc}
