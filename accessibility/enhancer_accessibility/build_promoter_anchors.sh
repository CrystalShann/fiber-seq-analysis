#!/usr/bin/env bash
# Build the THP-1 promoter anchor set used by enhancer_promoter_coupling.ipynb.
#
# Anchor rule (per request): use the THP-1 FANTOM5 dominant CAGE peak for a gene;
# where a gene has no CAGE peak, fall back to its GENCODE TSS.
# Restricted to protein-coding genes, matching the Tullius analysis.
#
# Output: thp1_promoter_anchors.bed  (BED6+, coordinate sorted, bgzip+tabix)
#   1 chrom
#   2 start       1-bp TSS (0-based)
#   3 end
#   4 gene_id     GENCODE v49 versioned id
#   5 gene_name
#   6 strand      GENCODE gene strand
#   7 tss_source  "cage" or "gencode"
set -euo pipefail

BGZIP=/project/spott/cshan/envs/bedtools/bin/bgzip
TABIX=/project/spott/cshan/envs/bedtools/bin/tabix

ANN=/project/spott/cshan/fiber-seq/macrophage_project/annotations
OUT="${ANN}/enhancer_anchors"
GTF=/project/spott/cshan/annotations/gencode.v49.primary_assembly.annotation.gtf
CAGE="${ANN}/fantom5_thp1/THP1.consensus.CAGE_peaks.dominantPerGene.bed.gz"
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
mkdir -p "$OUT"
cd "$OUT"

# 1. Protein-coding gene TSSs from GENCODE (TSS = gene start on +, gene end on -).
awk 'BEGIN{OFS="\t"}
     !/^#/ && $3=="gene" && /gene_type "protein_coding"/ {
        match($0, /gene_id "[^"]+"/);   gid  = substr($0, RSTART+9,  RLENGTH-10)
        match($0, /gene_name "[^"]+"/); gname= substr($0, RSTART+11, RLENGTH-12)
        tss = ($7=="+") ? $4-1 : $5-1
        print $1, tss, tss+1, gid, gname, $7
     }' "$GTF" | LC_ALL=C sort -k4,4 > "${TMP}/gencode_pc.bed"
echo "GENCODE protein-coding genes : $(wc -l < "${TMP}/gencode_pc.bed")"

# 2. THP-1 dominant CAGE peak per gene, keyed by gene_id.
zcat "$CAGE" | awk 'BEGIN{OFS="\t"} {print $4, $1, $2, $3}' \
  | LC_ALL=C sort -k1,1 > "${TMP}/cage_by_gene.tsv"
echo "THP-1 dominant CAGE peaks    : $(wc -l < "${TMP}/cage_by_gene.tsv")"

# 3. CAGE where available, GENCODE otherwise. The CAGE peak must sit on the same
#    chromosome as the gene (guards against a stray nearest-gene assignment).
LC_ALL=C join -a1 -1 4 -2 1 -t $'\t' -o '1.1,1.2,1.3,1.4,1.5,1.6,2.2,2.3' \
    "${TMP}/gencode_pc.bed" "${TMP}/cage_by_gene.tsv" \
  | awk 'BEGIN{OFS="\t"}
         {
            if ($7 != "" && $7 == $1) { start = $8; src = "cage" }
            else                      { start = $2; src = "gencode" }
            print $1, start, start+1, $4, $5, $6, src
         }' \
  | LC_ALL=C sort -k1,1 -k2,2n > thp1_promoter_anchors.bed

"$BGZIP" -f -c thp1_promoter_anchors.bed > thp1_promoter_anchors.bed.gz
"$TABIX" -f -p bed thp1_promoter_anchors.bed.gz

echo
echo "Wrote ${OUT}/thp1_promoter_anchors.bed(.gz)"
echo "  anchors        : $(wc -l < thp1_promoter_anchors.bed)"
echo "  by TSS source  :"
cut -f7 thp1_promoter_anchors.bed | sort | uniq -c | sed 's/^/    /'
