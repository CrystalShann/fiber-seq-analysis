#!/bin/bash
# Build the shared cCRE universe and gene windows for the co-accessibility analysis.
#
# Follows coaccess_fire_CREs_combined_samples_around_genes.R:
#   ENCODE cCREs  ->  keep those overlapping a FIRE peak  ->  assign to TSS +/- 10 kb
#   gene windows  ->  pairs are enumerated within a gene window
#
# Differences from Kevin, all forced or flagged:
#   * peaks are pooled across the four LPS timepoints (he pools 17 LCL samples), so
#     one universe serves every timepoint and pair identity is comparable. Peak calls
#     drift between timepoints; a per-timepoint universe would make boundary drift
#     look like biology.
#   * gene windows are keyed on gene_id, not gene_name. In the GENCODE v46 canonical
#     TSS bed 488 gene names are duplicated over 2,103 rows (Y_RNA alone appears 756
#     times); keying on the name fuses those loci into one pseudo-gene carrying ~1,034
#     cCREs, which alone would generate ~534k spurious pairs. Kevin's mapgen annotation
#     hits the same problem and he works around it with which.max(n_CREs); gene_id is
#     the correct fix.
#
# Usage:  bash 01_make_cre_universe.sh          (a few minutes)

set -uo pipefail

SAMPLES=(LPS_0 LPS_5 LPS_10 LPS_15)
WINDOW=10000          # gene window half-width around the TSS (Kevin: 10 kb)

FIRE_ROOT="/project/spott/lizarraga/pacbio_analysis/macrophage_project/merged_hifi_bams/FIRE"
OUT_ROOT="/project/spott/cshan/fiber-seq/macrophage_project/co-accessibility"
CRE_BED="/project/spott/cshan/annotations/GRCh38-cCREs.bed"
TSS_BED="/project/spott/cshan/annotations/gencodev46_Ensembl_canonical_TSS.bed"
SORT_TMP="/scratch/midway3/cshan"

BEDTOOLS=/project/spott/cshan/envs/bedtools/bin/bedtools
BGZIP=/project/spott/cshan/envs/dimelo/bin/bgzip
TABIX=/project/spott/cshan/envs/dimelo/bin/tabix

# chr1-22, X, Y.  chrM is dropped: no cCREs and no FIRE calls there.
CHROM_RE='^chr([1-9]|1[0-9]|2[0-2]|X|Y)$'

die() { echo "ERROR: $*" >&2; exit 1; }

[ -s "$CRE_BED" ] || die "cCRE bed not found: $CRE_BED"
[ -s "$TSS_BED" ] || die "canonical TSS bed not found (run make_gencode_v46_all_tss.sh): $TSS_BED"
mkdir -p "$SORT_TMP" 2>/dev/null || SORT_TMP="${SLURM_TMPDIR:-/tmp}"

uni="${OUT_ROOT}/universe"
mkdir -p "$uni" || die "cannot create $uni"

peaks_union="${uni}/fire_peaks_union.bed"
cre_universe="${uni}/cre_universe.bed"
windows="${uni}/gene_windows.bed"
cre_gene_map="${uni}/cre_gene_map.tsv.gz"

echo "cCREs      : ${CRE_BED}"
echo "TSS        : ${TSS_BED}"
echo "Window     : TSS +/- ${WINDOW} bp"
echo "Output     : ${uni}"
echo "Started    : $(date)"


###########################
# 1. FIRE peak union across timepoints
###########################
# The peaks bed carries a '#chrom' header and 29 columns. Columns 1-3 are the FULL
# peak interval; columns 4-5 are a narrower core interval around the local maximum,
# which is NOT what Kevin selects (he takes peak_start/peak_end). Use 1-3.
#
# Kevin applies no pass_coverage filter when building the peak set he screens cCREs
# against, so none is applied here either.

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


###########################
# 2. cCREs overlapping a FIRE peak  ->  the element universe
###########################
# subsetByOverlaps(CREs.gr, fire_peaks.gr) = any overlap >= 1 bp, keeping the cCRE
# coordinates (never the peak's).  CRE_ID = accession1.accession2.CRE_label, exactly
# Kevin's paste0(accession1,".",accession2,".",CRE_label).
#
# Output BED5: chrom, start, end, CRE_ID, CRE_label

echo "Selecting cCREs that overlap a FIRE peak..."
awk -v re="$CHROM_RE" 'BEGIN{OFS="\t"} $1 ~ re {print $1,$2,$3,$4,$5,$6}' "$CRE_BED" \
    | LC_ALL=C sort -k1,1 -k2,2n -T "$SORT_TMP" \
    | "$BEDTOOLS" intersect -a stdin -b "$peaks_union" -u \
    | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4"."$5"."$6,$6}' \
    | LC_ALL=C sort -k1,1 -k2,2n -T "$SORT_TMP" > "$cre_universe" || die "cCRE selection failed"
n_cre=$(wc -l < "$cre_universe")
[ "$n_cre" -gt 0 ] || die "no cCREs survived the FIRE peak filter"
echo "  ${n_cre} cCREs in the universe"
echo "  by class:"
cut -f5 "$cre_universe" | sort | uniq -c | sort -rn | sed 's/^/    /'


###########################
# 3. Gene windows: TSS +/- WINDOW
###########################
# The TSS bed holds 20 bp intervals centred on each canonical TSS, so TSS = start+10.
# col4 = gene_id;transcript_id;gene_name;transcript_type;tags
# Output BED7: chrom, start, end, gene_id, gene_name, transcript_type, strand

echo "Building gene windows..."
awk -F'\t' -v OFS='\t' -v W="$WINDOW" -v re="$CHROM_RE" '
    $1 ~ re {
        tss = $2 + 10
        s = tss - W; if (s < 0) s = 0
        split($4, a, ";")
        print $1, s, tss + W, a[1], a[3], a[4], $6
    }' "$TSS_BED" \
    | LC_ALL=C sort -k1,1 -k2,2n -T "$SORT_TMP" > "$windows" || die "gene window build failed"
echo "  $(wc -l < "$windows") gene windows"


###########################
# 4. cCRE -> gene window membership
###########################
# join_overlap_inner(gene_windows.gr, CREs_with_fire.gr): any overlap, so a cCRE in
# two overlapping windows belongs to both.

echo "Assigning cCREs to gene windows..."
{
  echo -e "CRE_ID\tchrom\tstart\tend\tCRE_label\tgene_id\tgene_name\ttranscript_type"
  "$BEDTOOLS" intersect -a "$cre_universe" -b "$windows" -wa -wb \
    | awk 'BEGIN{OFS="\t"} {print $4,$1,$2,$3,$5,$9,$10,$11}'
} | "$BGZIP" > "$cre_gene_map" || die "cCRE/gene intersect failed"

n_map=$(( $(zcat "$cre_gene_map" | wc -l) - 1 ))
[ "$n_map" -gt 0 ] || die "no cCREs fall inside a gene window"
n_genes=$(zcat "$cre_gene_map" | tail -n +2 | cut -f6 | sort -u | wc -l)
n_multi=$(zcat "$cre_gene_map" | tail -n +2 | cut -f6 | sort | uniq -c | awk '$1>1' | wc -l)
echo "  ${n_map} (cCRE, gene) memberships across ${n_genes} genes"
echo "  ${n_multi} genes carry more than one cCRE (these are the ones that yield pairs)"


###########################
# 5. Per-timepoint peak membership flags
###########################
# Whether each universe cCRE overlaps THAT timepoint's own peaks. Carried through to
# the results so a per-timepoint-peak-restricted view stays available as a filter,
# without changing the shared universe.

echo "Flagging per-timepoint peak membership..."
flags="${uni}/cre_in_timepoint_peaks.tsv.gz"
tmpdir=$(mktemp -d "${SORT_TMP}/creflag.XXXXXX") || die "mktemp failed"
cut -f4 "$cre_universe" > "${tmpdir}/ids"
for s in "${SAMPLES[@]}"; do
    zcat "${FIRE_ROOT}/${s}/${s}-fire-v0.1-peaks.bed.gz" \
        | awk -v re="$CHROM_RE" 'BEGIN{OFS="\t"} $1 ~ /^#/ {next} $1 ~ re {print $1,$2,$3}' \
        | LC_ALL=C sort -k1,1 -k2,2n -T "$SORT_TMP" \
        | "$BEDTOOLS" intersect -a "$cre_universe" -b stdin -c \
        | awk '{print ($6 > 0) ? 1 : 0}' > "${tmpdir}/${s}"
done
{
  printf "CRE_ID"; for s in "${SAMPLES[@]}"; do printf "\tin_%s_peaks" "$s"; done; printf "\n"
  paste "${tmpdir}/ids" $(for s in "${SAMPLES[@]}"; do echo "${tmpdir}/${s}"; done)
} | "$BGZIP" > "$flags" || die "peak flag build failed"
rm -rf "$tmpdir"
echo "  wrote $(basename "$flags")"

echo "Done -> ${uni}"
echo "Finished   : $(date)"
