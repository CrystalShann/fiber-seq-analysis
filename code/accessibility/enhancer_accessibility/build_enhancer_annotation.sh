#!/usr/bin/env bash
# Build the THP-1 enhancer anchor set used by the two notebooks in this directory.
#
# Definition (per the analysis plan):
#   anchor = the FANTOM5 enhancer "consensus TSS" (thickStart/thickEnd of the BED12,
#            i.e. the inferred mid position between the divergent CAGE initiation
#            blocks) that falls inside an ENCODE cCRE enhancer-like element.
#   class  = the ENCODE cCRE class (dELS / pELS), carried through from SCREEN.
#
# CAGE expression is deliberately NOT used to filter or rank (per request).
# THP-1 DNase overlap is added as an annotation column only, not a filter, so it
# can be applied later without rebuilding.
#
# Output: thp1_enhancer_anchors.bed  (BED6+, coordinate sorted, bgzip+tabix)
#   1 chrom
#   2 start        1-bp consensus TSS (0-based)
#   3 end
#   4 enhancer_id  FANTOM5 id (chrom:start-end of the enhancer region)
#   5 f5_score     FANTOM5 BED score
#   6 strand       always "." (FANTOM5 enhancers are unstranded)
#   7 encode_class dELS or pELS
#   8 ccre_id      ENCODE cCRE accession
#   9 enh_start    FANTOM5 enhancer region start
#  10 enh_end      FANTOM5 enhancer region end
#  11 dnase        1 if the anchor falls in a THP-1 ENCODE DNase peak (either file), else 0
set -euo pipefail

BEDTOOLS=/project/spott/cshan/envs/bedtools/bin/bedtools
BGZIP=/project/spott/cshan/envs/bedtools/bin/bgzip
TABIX=/project/spott/cshan/envs/bedtools/bin/tabix

ANN=/project/spott/cshan/fiber-seq/macrophage_project/annotations
PUB="${ANN}/enhancer_public_data"
OUT="${ANN}/enhancer_anchors"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$OUT"
cd "$OUT"

F5="${PUB}/fantom5_enhancer_hg38/F5.hg38.enhancers.bed.gz"
CCRE="${PUB}/GRCh38-cCREs.bed"
DNASE1="${PUB}/encode_thp1_dnase_hg38/ENCFF359IAQ.bed.gz"
DNASE2="${PUB}/encode_thp1_dnase_hg38/ENCFF891DFB.bed.gz"

# 1. FANTOM5 consensus TSS: thickStart/thickEnd of each enhancer, keeping the
#    enclosing enhancer region so it can be drawn later.
zcat "$F5" \
  | awk 'BEGIN{OFS="\t"} {print $1, $7, $8, $4, $5, ".", $2, $3}' \
  | LC_ALL=C sort -k1,1 -k2,2n > "${TMP}/f5_tss.bed"
echo "FANTOM5 enhancers            : $(wc -l < "${TMP}/f5_tss.bed")"

# 2. ENCODE cCRE enhancer-like elements only (distal + proximal ELS).
awk 'BEGIN{OFS="\t"} $6=="dELS" || $6=="pELS" {print $1, $2, $3, $5, $6}' "$CCRE" \
  | LC_ALL=C sort -k1,1 -k2,2n > "${TMP}/ccre_els.bed"
echo "ENCODE cCRE ELS              : $(wc -l < "${TMP}/ccre_els.bed")"

# 3. Keep consensus TSSs landing inside an ELS. -wa -wb, then collapse the rare
#    multi-cCRE hits by keeping the first (cCREs of the same class rarely overlap).
"$BEDTOOLS" intersect -a "${TMP}/f5_tss.bed" -b "${TMP}/ccre_els.bed" -wa -wb -sorted \
  | awk 'BEGIN{OFS="\t"} !seen[$1":"$2]++ {print $1, $2, $3, $4, $5, $6, $13, $12, $7, $8}' \
  > "${TMP}/anchors_noDNase.bed"
echo "F5 consensus TSS in cCRE ELS : $(wc -l < "${TMP}/anchors_noDNase.bed")"

# 4. THP-1 DNase peaks (union of the two ENCODE files) -> annotation column.
zcat "$DNASE1" "$DNASE2" \
  | awk 'BEGIN{OFS="\t"} !/^#/ {print $1, $2, $3}' \
  | LC_ALL=C sort -k1,1 -k2,2n \
  | "$BEDTOOLS" merge -i - > "${TMP}/dnase_union.bed"
echo "THP-1 DNase peaks (merged)   : $(wc -l < "${TMP}/dnase_union.bed")"

"$BEDTOOLS" intersect -a "${TMP}/anchors_noDNase.bed" -b "${TMP}/dnase_union.bed" -c -sorted \
  | awk 'BEGIN{OFS="\t"} {$11 = ($11 > 0 ? 1 : 0); print}' \
  > thp1_enhancer_anchors.bed

"$BGZIP" -f -c thp1_enhancer_anchors.bed > thp1_enhancer_anchors.bed.gz
"$TABIX" -f -p bed thp1_enhancer_anchors.bed.gz

echo
echo "Wrote ${OUT}/thp1_enhancer_anchors.bed(.gz)"
echo "  anchors            : $(wc -l < thp1_enhancer_anchors.bed)"
echo "  by ENCODE class    :"
cut -f7 thp1_enhancer_anchors.bed | sort | uniq -c | sed 's/^/    /'
echo "  in THP-1 DNase peak: $(awk '$11==1' thp1_enhancer_anchors.bed | wc -l)"
