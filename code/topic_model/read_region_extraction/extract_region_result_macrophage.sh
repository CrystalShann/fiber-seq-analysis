#!/bin/bash

# Extract per-region, per-sample results for read-level region plots (macrophage LPS series).
#
# Produces the file set that Kevin Luo's get_region_combined_results() (plots.R)
# reads out of a `parsed` directory, using the v2 (fibertools, non-dimelo) contract:
#
#   parsed/extracted.m6a.bed.gz      tabix slice of the per-chr `ft extract --m6a` bed
#   parsed/extracted.cpg.bed.gz      tabix slice of the per-chr `ft extract --cpg` bed
#   parsed/region.fiberhmm_tf.bed         per-chr FiberHMM recalled_tf bed (TF footprints)
#   parsed/region.fiberhmm_footprint.bed  per-chr FiberHMM recalled_footprint bed (nucleosomes)
#   parsed/fire_elements.bed         tabix slice of the FIRE per-read elements
#   parsed/fire_peaks.bed            tabix slice of the FIRE peak calls (may be empty)
#   parsed/fire.bed                  ft add-nucleosomes | ft fire | ft fire --extract
#
# Only fire.bed needs the BAM; everything else is a tabix slice of an existing file.
# No padding: tabix returns whole BED12 read records, so read spans stay full length
# and the plotting code clips positions to the region itself.
#
# Usage:
#   extract_region_result_macrophage.sh <sample_name> <region> <outname> <out_root>
# original example:
#   extract_region_result_macrophage.sh LPS_0 chr4:74098863-74099651 \
#       CXCL2_chr4_74098863_74099651 \
#       /project/spott/cshan/fiber-seq/macrophage_project/topic_model/region_plots/early_repsonse_genes
# Example (TSS +/- 1000 bp windows; outname = gene, one folder per gene under promoters/):
#   extract_region_result_macrophage.sh LPS_0 chr4:74098196-74100196 \
#       CXCL2 \
#       /project/spott/cshan/fiber-seq/macrophage_project/topic_model/region_plots/early_repsonse_genes/promoters

set -uo pipefail

sample_name=$1
region=$2
outname=$3
out_root=$4

## Fixed paths
FT_RESULT_DIR="/project/spott/cshan/fiber-seq/macrophage_project/FiberHMM/extract/ft_result_dir"
FIBERHMM_TF_DIR="/project/spott/cshan/fiber-seq/macrophage_project/FiberHMM/extract/firehmm_tf"
FIBERHMM_FP_DIR="/project/spott/cshan/fiber-seq/macrophage_project/FiberHMM/extract/firehmm_footprint"
FIRE_ROOT="/project/spott/lizarraga/pacbio_analysis/macrophage_project/merged_hifi_bams/FIRE"
BAM_DIR="/project/spott/lizarraga/pacbio_analysis/macrophage_project/merged_hifi_bams/FiberHMM"
fire_ver="v0.1"

## Binaries: absolute paths so this works under sbatch without conda activation
SAMTOOLS="/project/spott/cshan/envs/dimelo/bin/samtools"
TABIX="/project/spott/cshan/envs/dimelo/bin/tabix"
BGZIP="/project/spott/cshan/envs/dimelo/bin/bgzip"
FT="/project/spott/cshan/envs/fire-env/bin/ft"

die() { echo "ERROR: $*" >&2; exit 1; }

## Check region format
if [[ ! $region =~ ^chr[0-9XYM]+:[0-9]+-[0-9]+$ ]]; then
    # original:
    # die "Region must be chr:start-end (e.g. chr4:74098863-74099651), got '$region'"
    die "Region must be chr:start-end (e.g. chr4:74098196-74100196), got '$region'"
fi
region_chr=$(echo "$region" | cut -d: -f1)

## The FiberHMM per-chr split uses labels without the underscore (LPS_0 -> LPS0)
fh_label=${sample_name/_/}

## Input files
m6a_file="${FT_RESULT_DIR}/${sample_name}/extracted_results/m6a_by_chr/${sample_name}.ft_extracted_m6a.${region_chr}.bed.gz"
cpg_file="${FT_RESULT_DIR}/${sample_name}/extracted_results/cpg_by_chr/${sample_name}.ft_extracted_cpg.${region_chr}.bed.gz"
tf_file="${FIBERHMM_TF_DIR}/${fh_label}/${fh_label}_hmm_extracted_tf_${region_chr}.bed.gz"
footprint_file="${FIBERHMM_FP_DIR}/${fh_label}/${fh_label}_hmm_extracted_footprint_${region_chr}.bed.gz"
fire_elements_file="${FIRE_ROOT}/${sample_name}/additional-outputs-${fire_ver}/fire-peaks/${sample_name}-${fire_ver}-fire-elements.bed.gz"
fire_peaks_file="${FIRE_ROOT}/${sample_name}/${sample_name}-fire-${fire_ver}-peaks.bed.gz"
bam_file="${BAM_DIR}/${sample_name}.recalled.bam"

for f in "$m6a_file" "$cpg_file" "$tf_file" "$footprint_file" \
         "$fire_elements_file" "$fire_peaks_file" "$bam_file"; do
    [ -s "$f" ] || die "missing input: $f"
done

## Output layout
sample_dir="${out_root}/${outname}/${sample_name}"
parsed_dir="${sample_dir}/parsed"
mkdir -p "$parsed_dir" || die "cannot create $parsed_dir"

echo "sample_name=$sample_name"
echo "region=$region"
echo "outname=$outname"
echo "parsed_dir=$parsed_dir"

## 1. m6A / CpG: tabix slice, re-bgzip, re-index.
##    read_ft_extracted_region() uses the .tbi when present.
echo "Slicing m6A calls..."
"$TABIX" "$m6a_file" "$region" | "$BGZIP" -c > "${parsed_dir}/extracted.m6a.bed.gz" \
    || die "m6a slice failed"
"$TABIX" -f -p bed "${parsed_dir}/extracted.m6a.bed.gz" || die "m6a index failed"

echo "Slicing CpG calls..."
"$TABIX" "$cpg_file" "$region" | "$BGZIP" -c > "${parsed_dir}/extracted.cpg.bed.gz" \
    || die "cpg slice failed"
"$TABIX" -f -p bed "${parsed_dir}/extracted.cpg.bed.gz" || die "cpg index failed"

## 2. FiberHMM calls. `tf` holds the TF footprints; `footprint` holds every protected
##    segment, which is where the nucleosome-sized calls live - fibertools v0.8.2 does
##    not emit nucleosomes in `ft fire --extract`, so they have to come from here.
echo "Slicing FiberHMM TF footprints..."
"$TABIX" "$tf_file" "$region" > "${parsed_dir}/region.fiberhmm_tf.bed" || die "tf slice failed"

echo "Slicing FiberHMM footprints (nucleosomes)..."
"$TABIX" "$footprint_file" "$region" > "${parsed_dir}/region.fiberhmm_footprint.bed" \
    || die "footprint slice failed"

## 3. FIRE per-read elements and peak calls (plain text; read with read.table).
echo "Slicing FIRE elements and peaks..."
"$TABIX" "$fire_elements_file" "$region" > "${parsed_dir}/fire_elements.bed" || die "fire elements slice failed"
"$TABIX" "$fire_peaks_file" "$region" > "${parsed_dir}/fire_peaks.bed" || die "fire peaks slice failed"

## 4. fire.bed: the only step that needs the BAM.
##    Follows extract_region_result_v2.sh - no modkit update-tags (that is the dimelo
##    path) and no --min-msp/--skip-no-m6a filters, which would drop reads we need.
echo "Subsetting BAM..."
region_bam="${sample_dir}/region.bam"
"$SAMTOOLS" view -b "$bam_file" "$region" > "$region_bam" || die "samtools view failed"
"$SAMTOOLS" index "$region_bam" || die "samtools index failed"

echo "Computing FIRE over the region BAM..."
"$FT" add-nucleosomes "$region_bam" \
    | "$FT" fire - \
    | "$FT" fire --extract - > "${parsed_dir}/fire.bed" || die "ft fire pipeline failed"

## 5. Verify.
##    fire_peaks.bed may legitimately be empty (no FIRE peak in the window) - the R
##    code try()s it and simply drops the FIRE-peak panel. Everything else must have rows.
echo "Verifying outputs..."
for f in "${parsed_dir}/extracted.m6a.bed.gz" \
         "${parsed_dir}/extracted.m6a.bed.gz.tbi" \
         "${parsed_dir}/extracted.cpg.bed.gz" \
         "${parsed_dir}/extracted.cpg.bed.gz.tbi" \
         "${parsed_dir}/region.fiberhmm_tf.bed" \
         "${parsed_dir}/region.fiberhmm_footprint.bed" \
         "${parsed_dir}/fire_elements.bed" \
         "${parsed_dir}/fire.bed"; do
    [ -s "$f" ] || die "empty or missing output: $f"
done
[ -f "${parsed_dir}/fire_peaks.bed" ] || die "missing output: ${parsed_dir}/fire_peaks.bed"

n_tf_cols=$(head -1 "${parsed_dir}/region.fiberhmm_tf.bed" | awk -F'\t' '{print NF}')
[ "$n_tf_cols" -eq 15 ] || die "region.fiberhmm_tf.bed has $n_tf_cols columns, expected 15"

n_fp_cols=$(head -1 "${parsed_dir}/region.fiberhmm_footprint.bed" | awk -F'\t' '{print NF}')
[ "$n_fp_cols" -eq 13 ] || die "region.fiberhmm_footprint.bed has $n_fp_cols columns, expected 13"

n_fire_cols=$(head -1 "${parsed_dir}/fire.bed" | awk -F'\t' '{print NF}')
[ "$n_fire_cols" -eq 11 ] || die "fire.bed has $n_fire_cols columns, expected 11"

n_elem_cols=$(head -1 "${parsed_dir}/fire_elements.bed" | awk -F'\t' '{print NF}')
[ "$n_elem_cols" -eq 11 ] || die "fire_elements.bed has $n_elem_cols columns, expected 11"

if [ -s "${parsed_dir}/fire_peaks.bed" ]; then
    n_peak_cols=$(head -1 "${parsed_dir}/fire_peaks.bed" | awk -F'\t' '{print NF}')
    [ "$n_peak_cols" -eq 29 ] || die "fire_peaks.bed has $n_peak_cols columns, expected 29"
else
    echo "NOTE: no FIRE peak in $region for $sample_name (fire_peaks.bed is empty)."
fi

## Read-level sanity: how many of the m6A reads survived into fire.bed
n_m6a_reads=$(zcat "${parsed_dir}/extracted.m6a.bed.gz" | cut -f4 | sort -u | wc -l)
n_fire_reads=$(cut -f4 "${parsed_dir}/fire.bed" | sort -u | wc -l)
echo "reads in m6A slice: $n_m6a_reads ; reads in fire.bed: $n_fire_reads"

echo "Done: $parsed_dir"
