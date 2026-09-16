#!/bin/bash
# Build the shared FIRE peak universe for the FIRE-frequency analysis.
#
# Usage:  bash 01_make_fire_universe.sh

set -uo pipefail

SAMPLES=(LPS_0 LPS_5 LPS_10 LPS_15)
FIRE_ROOT="/project/spott/lizarraga/pacbio_analysis/macrophage_project/merged_hifi_bams/FIRE"
OUT_ROOT="/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency"
SORT_TMP="/scratch/midway3/cshan"

BEDTOOLS=/project/spott/cshan/envs/bedtools/bin/bedtools

# chr1-22, X, Y. chrM is dropped because there are no FIRE calls there.
CHROM_RE='^chr([1-9]|1[0-9]|2[0-2]|X|Y)$'

die() { echo "ERROR: $*" >&2; exit 1; }

mkdir -p "$SORT_TMP" 2>/dev/null || SORT_TMP="${SLURM_TMPDIR:-/tmp}"

uni="${OUT_ROOT}/universe"
mkdir -p "$uni" || die "cannot create $uni"

peaks_union="${uni}/fire_peaks_union.bed"
echo "Output     : ${uni}"
echo "Started    : $(date)"


###########################
# 1. FIRE peak union across timepoints
###########################
# The peaks bed carries a '#chrom' header and 29 columns. Columns 1-3 are the FULL
# peak interval; columns 4-5 are a narrower core interval around the local maximum,
# which is NOT what Kevin selects (he takes peak_start/peak_end). Use 1-3.
#
# No coverage filter is applied when building the shared peak universe.

echo "Pooling FIRE peaks across timepoints..."
: > "${uni}/.peaks_tmp"
for s in "${SAMPLES[@]}"; do
    pk="${FIRE_ROOT}/${s}/${s}-fire-v0.1-peaks.bed.gz"
    [ -s "$pk" ] || die "FIRE peaks not found: $pk"
    n=$(zcat "$pk" | awk -v re="$CHROM_RE" 'BEGIN{OFS="\t"} $1 ~ /^#/ {next} $1 ~ re {print $1,$2,$3}' \
        | tee -a "${uni}/.peaks_tmp" | wc -l)
    echo "  ${s}: ${n} peaks"
done
LC_ALL=C sort -k1,1 -k2,2n -T "$SORT_TMP" "${uni}/.peaks_tmp" \
    | "$BEDTOOLS" merge -i - > "$peaks_union" || die "peak merge failed"
rm -f "${uni}/.peaks_tmp"
echo "  union: $(wc -l < "$peaks_union") merged peak intervals"

echo "Done -> ${uni}"
echo "Finished   : $(date)"
