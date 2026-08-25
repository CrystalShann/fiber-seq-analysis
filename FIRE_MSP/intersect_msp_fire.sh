#!/bin/bash
#SBATCH --job-name=fire_msp_intersect
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --array=1-4
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/fire_msp_intersect_%A_%a.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/fire_msp_intersect_%A_%a.err

# Intersect per-read MSPs with FIRE peaks (macrophage LPS series)

# Outputs, per sample, in ${OUT_ROOT}/<sample>/:
#   <s>.FIRE_peaks_pass_coverage.autosomes.bed.gz   BED6 peaks; score = fire_cov/cov
#   <s>.MSP_BED6_single_<LEN>.bed.gz                BED6, one row per MSP (read name in col4)
#   <s>.MSP_<LEN>_FIRE_intersect.vsLPS0peaks.bed.gz peak coords + read name, one row per
#                                                   (MSP, peak) pair at >=50% overlap
#   <s>.FIRE_read_intersect.vsLPS0peaks.bed.gz      BED12 reads fully containing a peak
#
# Usage:
#   bash   intersect_msp_fire.sh LPS_0      # LPS_0 first, builds the peak universe
#   sbatch intersect_msp_fire.sh            # then the remaining timepoints as an array

set -uo pipefail

# define samples
SAMPLES=(LPS_0 LPS_5 LPS_10 LPS_15)
if [ -n "${1:-}" ]; then
    sample_name=$1
else
    sample_name=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}
fi

LEN=150   # min single-MSP length


FIRE_ROOT="/project/spott/lizarraga/pacbio_analysis/macrophage_project/merged_hifi_bams/FIRE"
OUT_ROOT="/project/spott/cshan/fiber-seq/macrophage_project/FIRE_MSP"
SORT_TMP="/scratch/midway3/cshan"
REF="/project/spott/reference/human/GRCh38/hg38.fa"

FT=/project/spott/cshan/envs/fire-env/bin/ft
SAMTOOLS=/project/spott/cshan/envs/dimelo/bin/samtools
BGZIP=/project/spott/cshan/envs/dimelo/bin/bgzip
TABIX=/project/spott/cshan/envs/dimelo/bin/tabix
BEDTOOLS=/project/spott/cshan/envs/bedtools/bin/bedtools
BED12TOBED6=/project/spott/cshan/envs/bedtools/bin/bed12ToBed6

ncore=8

die() { echo "ERROR: $*" >&2; exit 1; }

cram="${FIRE_ROOT}/${sample_name}/${sample_name}-fire-v0.1-filtered.cram"
peaks_in="${FIRE_ROOT}/${sample_name}/${sample_name}-fire-v0.1-peaks.bed.gz"
[ -s "$cram" ]     || die "FIRE CRAM not found: $cram"
[ -s "$peaks_in" ] || die "FIRE peaks not found: $peaks_in"

out_dir="${OUT_ROOT}/${sample_name}"
mkdir -p "$out_dir" || die "cannot create output dir"
if ! mkdir -p "$SORT_TMP" 2>/dev/null; then
    SORT_TMP="${SLURM_TMPDIR:-/tmp}"
    echo "WARNING: /scratch unavailable on this node, using sort tmp ${SORT_TMP}"
fi

peaks="${out_dir}/${sample_name}.FIRE_peaks_pass_coverage.autosomes.bed.gz"
peaks_lps0="${OUT_ROOT}/LPS_0/LPS_0.FIRE_peaks_pass_coverage.autosomes.bed.gz"
msp12="${out_dir}/${sample_name}.MSP_${LEN}.bed.gz"
msp6="${out_dir}/${sample_name}.MSP_BED6_single_${LEN}.bed.gz"
msp_x_fire="${out_dir}/${sample_name}.MSP_${LEN}_FIRE_intersect.vsLPS0peaks.bed.gz"
read_x_fire="${out_dir}/${sample_name}.FIRE_read_intersect.vsLPS0peaks.bed.gz"

echo "Sample     : ${sample_name}"
echo "CRAM       : ${cram}"
echo "Peaks      : ${peaks_in}"
echo "Output     : ${out_dir}"
echo "Started    : $(date)"


###########################
# filter FIRE peaks to calculate FIRE actuation --> BED6
###########################
## 1. FIRE peaks -> BED6: pass_coverage only, autosomes only, score = fire_coverage/coverage.
echo "Filtering FIRE peaks..."

# skip header
# in the last column of the FIRE peaks bed = pass_coverage, only keep the columns that pass the coverage
  # pass coverage = true
# only keep pass coverage peaks on autosomes
# create a bed6 file: calculate FIRE actuation = fraction of spanning fibers classified has FIRE
## score = FIRE coverage / coverage
## Coverage = The total number of individual DNA sequencing fibers (reads) spanning that genomic peak or window
## FIRE coverage = The number of individual fibers at that locus that are specifically predicted
## to contain an active regulatory element
  # 1. chr
  # 2. start
  # 3. end
  # 4. name = .
  # 5. score 
  # 6. strand = .
# sort, index and compress the bed6 file

zcat "$peaks_in" \
    | awk -F'\t' -v OFS='\t' '
        /^#/ { next }
        $NF == "true" && $1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ { print $1, $2, $3, ".", $7/$6, "." }' \
    | LC_ALL=C sort -k1,1 -k2,2n -k3,3n -T "$SORT_TMP" \
    | "$BGZIP" > "$peaks" || die "peak filtering failed"
"$TABIX" -f -p bed "$peaks" || die "tabix failed: $peaks"
n_peaks=$(zcat "$peaks" | wc -l)
[ "$n_peaks" -gt 0 ] || die "no peaks survived filtering"
echo "  ${n_peaks} peaks pass"


###########################
# Extract MSP > 150 bp from FIRE CRAM
###########################
## 2. Per-read MSPs from the FIRE CRAM (as/al tags already present), then one row per MSP


# extract MSP from FIRE CRAM using fibertools extract -x MSP length > 150 bp, this creates a bed12 file
## this contains a read as a row, and MSP are blocks
# convert BED12 block features to BED6 using bed12ToBed6
# now each row is a MSP block; the awk length check (> 100 bp) only drops malformed
# blocks - the operative MSP length filter is len(msp) > LEN at ft extract above
# sort, index and compress the file


if [ -s "$msp12" ]; then
    echo "Reusing existing ft extract output: $(basename "$msp12")"
else
    echo "ft extract --msp from the CRAM..."
    "$FT" extract -r -t ${ncore} -x "len(msp)>${LEN}" "$cram" --msp "$msp12" \
        || die "ft extract --msp failed"
    [ -s "$msp12" ] || die "empty output: $msp12"
fi

echo "Converting BED12 -> single-MSP BED6..."

"$BED12TOBED6" -i "$msp12" \
    | awk -F'\t' -v OFS='\t' '($3 - $2) > 100' \
    | LC_ALL=C sort -k1,1 -k2,2n -k3,3n -T "$SORT_TMP" \
    | "$BGZIP" > "$msp6" || die "bed12ToBed6 pipeline failed"
"$TABIX" -f -p bed "$msp6" || die "tabix failed: $msp6"
n_msp=$(zcat "$msp6" | wc -l)
[ "$n_msp" -gt 0 ] || die "no MSPs > ${LEN} bp"
echo "  ${n_msp} single MSPs"


###########################
# Find MSP that overlaps FIRE peaks at LPS_0 with at least 50% of either the MSP or the peak = numerator
###########################


## 3. MSPs vs FIREpeaks, >=50% overlap in either direction.
##    Output rows are peak coords + the MSP's read name (cols 4-6 of the MSP BED6),


# intersect MSP BED6 with the LPS_0 peak universe; -f 0.5 -F 0.5 -e = keep a match
# when the overlap covers >=50% of the MSP OR >=50% of the peak (either is enough)
  ## include all columns of the MSP BED6 and the peak in the output
# rearrange and only keep: peak_chr    peak_start    peak_end    read_name    score    strand


[ -s "$peaks_lps0" ] || die "LPS_0 peak universe not found (run LPS_0 first): $peaks_lps0"

echo "Intersecting MSPs with LPS_0 FIRE peaks..."
"$BEDTOOLS" intersect -wb -f 0.5 -F 0.5 -e -a "$msp6" -b "$peaks_lps0" \
    | awk -F'\t' -v OFS='\t' '{ print $7, $8, $9, $4, $5, $6 }' \
    | LC_ALL=C sort -k1,1 -k2,2n -k3,3n -T "$SORT_TMP" \
    | "$BGZIP" > "$msp_x_fire" || die "MSP x peak intersect failed"
"$TABIX" -f -p bed "$msp_x_fire" || die "tabix failed: $msp_x_fire"
echo "  $(zcat "$msp_x_fire" | wc -l) MSP-peak pairs"

###########################
# Find all fibers that fully span the FIRE peak = denominator
###########################

## 4. Reads fully containing a peak for per-peak fiber denominators

# drop secondary and supplementary alignments using -F 0x900, so only primary
# alignments count (secondary/supplementary would fabricate distant co-occurrences)
  ## -u requests uncompressed BAM output that can be piped into BEDtools
# intersect the reads from the FIRE CRAM with the LPS_0 peak universe
  ## LPS_0.FIRE_peaks_pass_coverage.autosomes.bed.gz
  ## -F 1 requires 100% of the FIRE peak to lie inside the read alignment


echo "Intersecting reads with LPS_0 FIRE peaks..."

"$SAMTOOLS" view -@ ${ncore} -u -F 0x900 -T "$REF" "$cram" \
    | "$BEDTOOLS" intersect -bed -F 1 -abam stdin -b "$peaks_lps0" \
    | "$BGZIP" > "$read_x_fire" || die "read x peak intersect failed"
echo "  $(zcat "$read_x_fire" | wc -l) read-peak pairs"

echo "Done -> ${out_dir}"
echo "Disk: $(du -sh "$out_dir" | cut -f1)"
echo "Finished   : $(date)"
