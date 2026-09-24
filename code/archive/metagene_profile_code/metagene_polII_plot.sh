#!/bin/bash
#SBATCH --job-name=deeptools_2
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --time=12:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/deeptools_PolII_2_%x_%j.out

set -euo pipefail

# -----------------------------
# 0) Environment
# -----------------------------
module load python/miniforge-25.3.0
source activate /project/spott/cshan/envs/deeptools
unset MPLBACKEND
export MPLBACKEND=Agg
echo "Using MPLBACKEND=${MPLBACKEND}"

BIGWIGMERGE="/project/spott/cshan/tools/bigWigMerge"
BEDGRAPHTOBIGWIG="/project/spott/cshan/tools/bedGraphToBigWig"

# -----------------------------
# 1) Inputs
# -----------------------------
FP_DIR="/project/spott/1_Shared_projects/LCL_Fiber_seq/FiberHMM/merged/combined/joint_trained_tracks"
CHROM_SIZES="/project/spott/cshan/annotations/hg38.chrom.sizes"

SAMPLE="AL10_bc2178_19130"

DNASE_BW="/project/spott/cshan/annotations/DNase.ENCFF743ULW.bigWig"
GROCAP_BW="/project/spott/cshan/annotations/GRO_cap_ENCFF799JUK.bigWig"
FIRE_BW_DIR="/project/spott/1_Shared_projects/LCL_Fiber_seq/FIRE/results/${SAMPLE}/trackHub-v0.1/bw"
FIRE_ALL_BW="${FIRE_BW_DIR}/all.fire.coverage.bw"
FIRE_LINKER_BW="${FIRE_BW_DIR}/all.linker.coverage.bw"
FIRE_NUC_BW="${FIRE_BW_DIR}/all.nucleosome.coverage.bw"

FP_SIZES=("10-30" "30-60" "60-80" "140-160")

TSS_NAMES=("gencode_v49_TSS" "fantom5_CAGE")

# fantom5: convert BED9 gz -> uncompressed BED6 (deeptools sortMatrix cannot read gzipped BED)
FANTOM5_GZ="/project/spott/cshan/annotations/fantom5/fantom5.hg38.LCL.consensus.CAGE_peaks.bed.gz"
FANTOM5_BED6="/project/spott/cshan/annotations/fantom5/fantom5.hg38.LCL.consensus.CAGE_peaks.bed6.bed"
if [[ ! -s "${FANTOM5_BED6}" ]]; then
  echo "Converting FANTOM5 BED9 gz -> BED6 uncompressed"
  zcat "${FANTOM5_GZ}" | awk 'BEGIN{OFS="\t"} !/^#/ {print $1,$2,$3,$4,$5,$6}' > "${FANTOM5_BED6}"
  echo "  written: ${FANTOM5_BED6} ($(wc -l < "${FANTOM5_BED6}") lines)"
fi

TSS_FILES=(
  "/project/spott/cshan/annotations/TSS.gencode.v49.bed"
  "${FANTOM5_BED6}"
)

UPSTREAM=2000
DOWNSTREAM=2000
BIN_SIZE=20
THREADS="${SLURM_CPUS_PER_TASK:-4}"

OUT_DIR="/project/spott/cshan/fiber-seq/results/PolII/deeptools_v2"
TMP_DIR="${OUT_DIR}/tmp"
mkdir -p "${OUT_DIR}" "${TMP_DIR}"

# -----------------------------
# 2) Merge per-chromosome footprint bigWigs
# -----------------------------
MERGED_FP_BWS=()

for size in "${FP_SIZES[@]}"; do
  list_file="${TMP_DIR}/fp_${size}.list"
  find "${FP_DIR}" -maxdepth 2 -type f -name "combined_chr*_${size}bp_fp_cov.bw" | sort > "${list_file}"

  if [[ ! -s "${list_file}" ]]; then
    echo "ERROR: no footprint bigWigs found for ${size}" >&2
    exit 1
  fi

  merged_bg="${TMP_DIR}/combined_${size}bp_fp_cov.bedGraph"
  merged_bg_sorted="${TMP_DIR}/combined_${size}bp_fp_cov.sorted.bedGraph"
  merged_bw="${OUT_DIR}/combined_hg38_${size}bp_fp_cov.bw"

  if [[ -s "${merged_bw}" ]]; then
    echo "Merged bigWig exists; skipping rebuild for ${size}: ${merged_bw}"
  else
    echo "Merging footprint bigWigs for size ${size}"
    "${BIGWIGMERGE}" $(tr '\n' ' ' < "${list_file}") "${merged_bg}"
    LC_ALL=C sort -k1,1 -k2,2n "${merged_bg}" > "${merged_bg_sorted}"
    "${BEDGRAPHTOBIGWIG}" "${merged_bg_sorted}" "${CHROM_SIZES}" "${merged_bw}"
  fi

  MERGED_FP_BWS+=("${merged_bw}")
done

SCORES=(
  "${MERGED_FP_BWS[0]}"
  "${MERGED_FP_BWS[1]}"
  "${MERGED_FP_BWS[2]}"
  "${MERGED_FP_BWS[3]}"
  "${DNASE_BW}"
  "${GROCAP_BW}"
  "${FIRE_ALL_BW}"
  "${FIRE_LINKER_BW}"
  "${FIRE_NUC_BW}"
)

LABELS=(
  "FP_10-30bp"
  "FP_30-60bp"
  "FP_60-80bp"
  "FP_140-160bp"
  "DNase"
  "GROcap"
  "FIRE_all"
  "FIRE_linker"
  "FIRE_nucleosome"
)

# -----------------------------
# 3) computeMatrix + heatmap + combined profile for each TSS set
# -----------------------------
for i in "${!TSS_FILES[@]}"; do
  NAME="${TSS_NAMES[$i]}"
  TSS="${TSS_FILES[$i]}"

  MATRIX_GZ="${OUT_DIR}/${NAME}.${SAMPLE}.matrix.gz"
  MATRIX_TAB="${OUT_DIR}/${NAME}.${SAMPLE}.matrix.tab"
  SORTED_REGIONS="${OUT_DIR}/${NAME}.${SAMPLE}.sorted_regions.bed"
  PROFILE_TAB="${OUT_DIR}/${NAME}.${SAMPLE}.profile_data.tab"

  echo "Processing ${NAME}"

  if [[ -s "${MATRIX_GZ}" ]]; then
    echo "Matrix exists; skipping computeMatrix for ${NAME}: ${MATRIX_GZ}"
  else
    MPLBACKEND=Agg computeMatrix reference-point \
      --referencePoint TSS \
      -R "${TSS}" \
      -S "${SCORES[@]}" \
      --samplesLabel "${LABELS[@]}" \
      -b "${UPSTREAM}" \
      -a "${DOWNSTREAM}" \
      --binSize "${BIN_SIZE}" \
      --missingDataAsZero \
      --skipZeros \
      --smartLabels \
      -p "${THREADS}" \
      -o "${MATRIX_GZ}" \
      --outFileNameMatrix "${MATRIX_TAB}" \
      --outFileSortedRegions "${SORTED_REGIONS}"
  fi

  # Heatmap with summary profile on top
  # --zMin 0 ... (one per sample) sets each track's color scale minimum independently;
  # omitting --zMax lets each track auto-scale its maximum independently (no shared scale)
  MPLBACKEND=Agg plotHeatmap \
    -m "${MATRIX_GZ}" \
    --samplesLabel "${LABELS[@]}" \
    --sortRegions descend \
    --sortUsing mean \
    --whatToShow 'plot, heatmap and colorbar' \
    --plotType lines \
    --averageTypeSummaryPlot mean \
    --colorMap RdYlBu \
    --zMin 0 0 0 0 0 0 0 0 0 \
    --plotTitle "${NAME} (${SAMPLE}): TSS metagene heatmap" \
    --xAxisLabel "Distance to TSS (bp)" \
    --refPointLabel "TSS" \
    --yAxisLabel "Mean signal" \
    --dpi 300 \
    -out "${OUT_DIR}/${NAME}.${SAMPLE}.heatmap.pdf"

  # Combined profile: all tracks overlaid on the same panel (--perGroup)
  MPLBACKEND=Agg plotProfile \
    -m "${MATRIX_GZ}" \
    --samplesLabel "${LABELS[@]}" \
    --plotType lines \
    --averageType mean \
    --perGroup \
    --plotTitle "${NAME} (${SAMPLE}): TSS metagene profile" \
    --legendLocation upper-right \
    --yAxisLabel "Mean signal" \
    --refPointLabel "TSS" \
    --dpi 300 \
    --outFileNameData "${PROFILE_TAB}" \
    -out "${OUT_DIR}/${NAME}.${SAMPLE}.profile.pdf"
done

echo "Done. Outputs: ${OUT_DIR}"