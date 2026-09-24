#!/usr/bin/env bash
# Build a THP-1 CAGE consensus TSS annotation from FANTOM5 CTSS files.
#
# Steps:
#   1. Download the three THP-1 hCAGE CTSS BED files (hg38) from FANTOM5.
#   2. Consensus = CTSS positions present in ALL THREE replicates, tag counts summed.
#   3. Assign the nearest gencode gene to each CTSS (bedtools closest),
#      keeping the CAGE strand (not the gencode strand).
#
# Output (mirrors the LCL reference layout, BED9):
#   col1-3  chrom / start / end   (single-bp CTSS)
#   col4    gene_id               (nearest gencode gene)
#   col5    summed tag count      (across the 3 reps)
#   col6    strand                (from CAGE)
#   col7-8  thickStart/thickEnd   (= the CTSS point)
#   col9    itemRgb (0,0,0)
#
# Run from anywhere; all paths are absolute.
set -euo pipefail

BEDTOOLS=/project/spott/cshan/envs/bedtools/bin/bedtools
GENCODE_TSS=/project/spott/cshan/annotations/TSS.gencode.v49.bed
OUTDIR=/project/spott/cshan/fiber-seq/macrophage_project/annotations/fantom5_thp1
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# htslib provides bgzip / tabix
. /software/modules/init/bash
module load htslib/1.16

# ---------------------------------------------------------------------------
# 1. Download the three THP-1 CTSS replicates
# ---------------------------------------------------------------------------
BASE="https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/basic/human.cell_line.hCAGE"

# associative array uses named keys
declare -A REPS=(
  [THP1_fresh]="acute%2520myeloid%2520leukemia%2520%2528FAB%2520M5%2529%2520cell%2520line%253aTHP-1%2520%2528fresh%2529.CNhs10722.10399-106A3.hg38.nobarcode.ctss.bed.gz"
  [THP1_revived]="acute%2520myeloid%2520leukemia%2520%2528FAB%2520M5%2529%2520cell%2520line%253aTHP-1%2520%2528revived%2529.CNhs10723.10400-106A4.hg38.nobarcode.ctss.bed.gz"
  [THP1_thawed]="acute%2520myeloid%2520leukemia%2520%2528FAB%2520M5%2529%2520cell%2520line%253aTHP-1%2520%2528thawed%2529.CNhs10724.10405-106A9.hg38.nobarcode.ctss.bed.gz"
)

for rep in "${!REPS[@]}"; do
  # -s is true when the file exists and has a size greater than zero
  if [[ ! -s "${rep}.ctss.bed.gz" ]]; then
    echo "Downloading ${rep} ..."
    # -O Writes the downloaded content to the specified local filename
    wget -q -O "${rep}.ctss.bed.gz" "${BASE}/${REPS[$rep]}"
  fi
done

# CTSS BED columns: chrom start end  name(chr:s..e,strand)  tag_count  strand
# ---------------------------------------------------------------------------
# 2. Consensus: keep CTSS present in all three reps, sum tag counts
# ---------------------------------------------------------------------------
# Key each CTSS by chrom:start:end:strand; sum col5 tags and count how many
# replicates contain it; keep those seen in all three.

# decompress gzip files and writes the output
# awk 'BEGIN{OFS="\t"} snippet initializes an awk script by setting the Output Field Separator
# (OFS) to a tab character (\t) before any file lines are processed
zcat THP1_fresh.ctss.bed.gz THP1_revived.ctss.bed.gz THP1_thawed.ctss.bed.gz \
  | awk 'BEGIN{OFS="\t"}
         {
           key=$1":"$2":"$3":"$6
           tags[key]+=$5 # for each CTDD key, it adds tag counts from column 5
           n[key]++
           chrom[key]=$1; s[key]=$2; e[key]=$3; str[key]=$6
         }
         END{
           for (k in n)
             if (n[k]==3) # loops through all the keys, and keep CTSS occurs exactly three times
               print chrom[k], s[k], e[k], ".", tags[k], str[k]
         }' \
  | sort -k1,1 -k2,2n \
  > THP1.consensus.CTSS.bed6

# save the consensus CTSS file
echo "Consensus CTSS (present in all 3 reps): $(wc -l < THP1.consensus.CTSS.bed6)"

# ---------------------------------------------------------------------------
# 3. Assign nearest gencode gene (keep CAGE strand)
# ---------------------------------------------------------------------------
# bedtools closest needs both inputs sorted the same way
sort -k1,1 -k2,2n "$GENCODE_TSS" > gencode_tss.sorted.bed

# closest gene per CTSS; -t first breaks ties deterministically.
# Output cols after closest: <ctss 1-6> <gene 7-12>; gene_id is col 10.

# For each CTSS, bedtools finds the nearest GENCODE TSS
# -t first chooses the first tied feature encountered in the sorted -b file
"$BEDTOOLS" closest -a THP1.consensus.CTSS.bed6 -b gencode_tss.sorted.bed -t first \
  | awk 'BEGIN{OFS="\t"}
         $10!="." {
           # BED9: gene_id in col4, summed tags col5, CAGE strand col6, thick=CTSS
           print $1, $2, $3, $10, $5, $6, $2, $3, "0,0,0"
         }' \
  | sort -k1,1 -k2,2n \
  > THP1.consensus.CAGE_peaks.withGene.bed

bgzip -f THP1.consensus.CAGE_peaks.withGene.bed
tabix -p bed THP1.consensus.CAGE_peaks.withGene.bed.gz

echo "Wrote $OUTDIR/THP1.consensus.CAGE_peaks.withGene.bed.gz"
echo "  $(zcat THP1.consensus.CAGE_peaks.withGene.bed.gz | wc -l) TSS with gene assignment"

# ---------------------------------------------------------------------------
# 4. Dominant CTSS per gene (one anchor per gene = the highest-tag CTSS)
# ---------------------------------------------------------------------------
# Per-CTSS windows are far too dense (median gap ~1bp) for genome-wide
# metaprofiles, so the analysis anchors on one dominant TSS per gene: the CTSS
# with the largest summed tag count for each gene_id.

# sort by gene ID in column 4. within each gene, the CTSS with the largest tag count appears first
# descending order sort CTSS score
zcat THP1.consensus.CAGE_peaks.withGene.bed.gz \
  | sort -k4,4 -k5,5nr \
  | awk 'BEGIN{OFS="\t"} !seen[$4]++' \
  | sort -k1,1 -k2,2n \
  > THP1.consensus.CAGE_peaks.dominantPerGene.bed

bgzip -f THP1.consensus.CAGE_peaks.dominantPerGene.bed
tabix -p bed THP1.consensus.CAGE_peaks.dominantPerGene.bed.gz

echo "Wrote $OUTDIR/THP1.consensus.CAGE_peaks.dominantPerGene.bed.gz"
echo "  $(zcat THP1.consensus.CAGE_peaks.dominantPerGene.bed.gz | wc -l) genes (dominant CTSS)"
