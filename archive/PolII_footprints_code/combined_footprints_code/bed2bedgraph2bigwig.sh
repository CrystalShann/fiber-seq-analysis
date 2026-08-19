#!/bin/bash
#SBATCH --job-name=bedgraph_bw
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/bedgraph_bw_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/bedgraph_bw_%j.err

module load python/miniforge-25.3.0
source activate /project/spott/cshan/envs/bedtools

INPUT_DIR="/project/spott/cshan/fiber-seq/results/PolII/FIRE_nuc_by_chr_combined_sample/parsed_bed_files/sorted_indexed"
BEDGRAPH_DIR="/project/spott/cshan/fiber-seq/results/PolII/FIRE_nuc_by_chr_combined_sample/bedgraph_files"
BIGWIG_DIR="/project/spott/cshan/fiber-seq/results/PolII/FIRE_nuc_by_chr_combined_sample/bigwig_files"
CHROM_SIZES="/project/spott/cshan/annotations/hg38.chrom.sizes"
BEDGRAPH_TO_BIGWIG="/project/spott/cshan/tools/bedGraphToBigWig"

mkdir -p "${BEDGRAPH_DIR}" "${BIGWIG_DIR}"

mapfile -t BED_FILES < <(
    find "${INPUT_DIR}" -maxdepth 1 -type f -name "all_samples_nuc_chr*_parsed_sorted.bed.gz" \
        ! -name "*_sorted_parsed_sorted.bed.gz" \
        -size +1k | sort -V
)

if [[ ${#BED_FILES[@]} -eq 0 ]]; then
    echo "No BED files found in ${INPUT_DIR}"
    exit 1
fi

echo "Found ${#BED_FILES[@]} BED files"

for BED_FILE in "${BED_FILES[@]}"; do
    BASE_NAME="$(basename "${BED_FILE}" .bed.gz)"
    BEDGRAPH_FILE="${BEDGRAPH_DIR}/${BASE_NAME}.bedGraph"
    BIGWIG_FILE="${BIGWIG_DIR}/${BASE_NAME}.bw"

    echo "Input BED: ${BED_FILE}"
    echo "Output bedGraph: ${BEDGRAPH_FILE}"
    echo "Output bigWig: ${BIGWIG_FILE}"

    bedtools genomecov \
        -i <(zcat "${BED_FILE}") \
        -g "${CHROM_SIZES}" \
        -bg \
        > "${BEDGRAPH_FILE}"

    "${BEDGRAPH_TO_BIGWIG}" "${BEDGRAPH_FILE}" "${CHROM_SIZES}" "${BIGWIG_FILE}"

    echo "Finished ${BASE_NAME}"
done
