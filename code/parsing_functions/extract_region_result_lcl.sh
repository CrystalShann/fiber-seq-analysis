#!/bin/bash

# Extract per-region, per-sample results for read-level region plots (LCL samples).
#
# Sibling of extract_region_result_macrophage.sh: it writes the same parsed/ file
# set, plus the ft nucleosomes, so load_region_results() reads both datasets
# (use nucleosome_source = "ft" for this one):
#
#   parsed/extracted.m6a.bed.gz      tabix slice of the per-chr `ft extract --m6a` bed
#   parsed/extracted.cpg.bed.gz      tabix slice of the per-chr `ft extract --cpg` bed
#   parsed/extracted.nuc.bed.gz      tabix slice of the per-chr `ft extract --nuc` bed (nucleosomes)
#   parsed/region.fiberhmm_tf.bed    tabix slice of the FiberHMM v2 recalled_tf bed (TF footprints)
#   parsed/fire_elements.bed         tabix slice of the FIRE per-read elements
#   parsed/fire_peaks.bed            tabix slice of the FIRE peak calls (may be empty)
#   parsed/fire.bed                  tabix slice of <sample>.fire_all.bed.gz (`ft fire --extract --all`)
#
# Every output is a tabix slice of an existing genome-wide file: no BAM is read and
# no ft command runs. fire.bed holds the segments that overlap the region rather than
# whole reads, so coverage inside the region is exact and read ends outside it stop
# at the last overlapping segment. Its nucleosome rows are the same calls as
# extracted.nuc.bed.gz; load_region_results(nucleosome_source = "ft") drops them.
#
# Usage:
#   extract_region_result_lcl.sh <sample_name> <region> <outname> <out_root>
# Example:
#   extract_region_result_lcl.sh AL10_bc2178_19130 chr21:36384224-36390964 \
#       CHAF1B /project/spott/cshan/fiber-seq/LCL_project/region_plots
# Called from R by extract_region_results() / plot_region_example()
# (code/parsing_functions/).

set -uo pipefail

sample_name=$1
region=$2
outname=$3
out_root=$4

## Fixed paths
FIRE_ROOT="/project/spott/1_Shared_projects/LCL_Fiber_seq/FIRE/results"
FIBERHMM_V2_DIR="/project/spott/1_Shared_projects/LCL_Fiber_seq/FiberHMM/FiberHMM_v2/results/results_v1_trained_model"
fire_ver="v0.1"

## Binaries: absolute paths so this works under sbatch without conda activation
TABIX="/project/spott/cshan/envs/dimelo/bin/tabix"
BGZIP="/project/spott/cshan/envs/dimelo/bin/bgzip"

die() { echo "ERROR: $*" >&2; exit 1; }

## Check region format
if [[ ! $region =~ ^chr[0-9XYM]+:[0-9]+-[0-9]+$ ]]; then
    die "Region must be chr:start-end (e.g. chr21:36384224-36390964), got '$region'"
fi
region_chr=$(echo "$region" | cut -d: -f1)

## Input files
ft_dir="${FIRE_ROOT}/${sample_name}/extracted_results"
m6a_file="${ft_dir}/m6a_by_chr/${sample_name}.ft_extracted_m6a.${region_chr}.bed.gz"
cpg_file="${ft_dir}/cpg_by_chr/${sample_name}.ft_extracted_cpg.${region_chr}.bed.gz"
nuc_file="${ft_dir}/nuc_by_chr/${sample_name}.ft_extracted_nuc.${region_chr}.bed.gz"
fire_all_file="${ft_dir}/${sample_name}.fire_all.bed.gz"
tf_file="${FIBERHMM_V2_DIR}/${sample_name}/${sample_name}.recalled_tf.bed.gz"
fire_elements_file="${FIRE_ROOT}/${sample_name}/additional-outputs-${fire_ver}/fire-peaks/${sample_name}-${fire_ver}-fire-elements.bed.gz"
fire_peaks_file="${FIRE_ROOT}/${sample_name}/${sample_name}-fire-${fire_ver}-peaks.bed.gz"

for f in "$m6a_file" "$cpg_file" "$nuc_file" "$fire_all_file" "$tf_file" \
         "$fire_elements_file" "$fire_peaks_file"; do
    [ -s "$f" ] || die "missing input: $f"
    [ -s "${f}.tbi" ] || die "missing tabix index: ${f}.tbi"
done

## Output layout
sample_dir="${out_root}/${outname}/${sample_name}"
parsed_dir="${sample_dir}/parsed"
mkdir -p "$parsed_dir" || die "cannot create $parsed_dir"

echo "sample_name=$sample_name"
echo "region=$region"
echo "outname=$outname"
echo "parsed_dir=$parsed_dir"

## 1. m6A / CpG / nucleosomes: tabix slice, re-bgzip, re-index.
for kind in m6a cpg nuc; do
    in_var="${kind}_file"
    echo "Slicing ${kind} calls..."
    "$TABIX" "${!in_var}" "$region" | "$BGZIP" -c > "${parsed_dir}/extracted.${kind}.bed.gz" \
        || die "${kind} slice failed"
    "$TABIX" -f -p bed "${parsed_dir}/extracted.${kind}.bed.gz" || die "${kind} index failed"
done

## 2. FiberHMM v2 TF footprints (BED15, the same format as the macrophage `tf` calls).
echo "Slicing FiberHMM v2 TF footprints..."
"$TABIX" "$tf_file" "$region" > "${parsed_dir}/region.fiberhmm_tf.bed" || die "tf slice failed"

## 3. FIRE per-read elements and peak calls (plain text; read with read.table).
echo "Slicing FIRE elements and peaks..."
"$TABIX" "$fire_elements_file" "$region" > "${parsed_dir}/fire_elements.bed" || die "fire elements slice failed"
"$TABIX" "$fire_peaks_file" "$region" > "${parsed_dir}/fire_peaks.bed" || die "fire peaks slice failed"

## 4. fire.bed: FIRE / linker / nucleosome segments from `ft fire --extract --all`.
echo "Slicing the FIRE segmentation..."
"$TABIX" "$fire_all_file" "$region" > "${parsed_dir}/fire.bed" || die "fire_all slice failed"

## 5. Verify.
##    The m6A, nucleosome and FIRE segmentation slices must have rows. TF footprints,
##    FIRE elements and peaks may legitimately be empty in a small window of one
##    sample; load_region_results() reads an empty TF file as no footprints.
echo "Verifying outputs..."
for kind in m6a nuc; do
    [ -n "$(zcat "${parsed_dir}/extracted.${kind}.bed.gz" | head -c 1)" ] \
        || die "no ${kind} rows in $region for $sample_name"
    [ -s "${parsed_dir}/extracted.${kind}.bed.gz.tbi" ] || die "missing index: extracted.${kind}.bed.gz.tbi"
done
[ -s "${parsed_dir}/fire.bed" ] || die "empty output: ${parsed_dir}/fire.bed"

n_fire_cols=$(head -1 "${parsed_dir}/fire.bed" | awk -F'\t' '{print NF}')
[ "$n_fire_cols" -eq 11 ] || die "fire.bed has $n_fire_cols columns, expected 11"

check_cols() {   # <file> <expected columns> <label>
    if [ -s "$1" ]; then
        n=$(head -1 "$1" | awk -F'\t' '{print NF}')
        [ "$n" -eq "$2" ] || die "$3 has $n columns, expected $2"
    else
        echo "NOTE: no $3 rows in $region for $sample_name ($(basename "$1") is empty)."
    fi
}
check_cols "${parsed_dir}/region.fiberhmm_tf.bed" 15 "FiberHMM TF footprint"
check_cols "${parsed_dir}/fire_elements.bed" 11 "FIRE element"
check_cols "${parsed_dir}/fire_peaks.bed" 29 "FIRE peak"

## Read-level sanity: how many of the m6A reads appear in fire.bed
n_m6a_reads=$(zcat "${parsed_dir}/extracted.m6a.bed.gz" | cut -f4 | sort -u | wc -l)
n_fire_reads=$(cut -f4 "${parsed_dir}/fire.bed" | sort -u | wc -l)
echo "reads in m6A slice: $n_m6a_reads ; reads in fire.bed: $n_fire_reads"

echo "Done: $parsed_dir"
