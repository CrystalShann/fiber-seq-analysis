#!/usr/bin/env bash
# LiftOver the public THP-1 ATAC (GSE201351) and H3K27ac ChIP (GSE201352) tracks
# from hg19 -> hg38 with rtracklayer (see liftover_tracks.R for the per-file logic).
#
# For every *.bw and *_peaks.narrowPeak.gz under enhancer_public_data/{atac,chip}
# it writes an hg38 copy into the sibling {atac,chip}_hg38/ dirs. narrowPeaks are
# additionally coordinate-sorted, bgzip'd and tabix-indexed to match the rest of
# the annotation pipeline. Already-present outputs are skipped, so re-runs resume.
#
# Requires: the UCSC chain file, fetched once with
#   curl -sO https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
#   gunzip hg19ToHg38.over.chain.gz




set -euo pipefail

RSCRIPT=/software/R-4.4.1-el8-x86_64/bin/Rscript
BGZIP=/project/spott/cshan/envs/bedtools/bin/bgzip
TABIX=/project/spott/cshan/envs/bedtools/bin/tabix

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LIFT="${HERE}/liftover_tracks.R"

ANN=/project/spott/cshan/fiber-seq/macrophage_project/annotations
PUB="${ANN}/enhancer_public_data"
CHAIN=/project/spott/cshan/annotations/hg19ToHg38.over.chain
SIZES=/project/spott/cshan/annotations/hg38.chrom.sizes

[[ -f "$CHAIN" ]] || { echo "ERROR: chain file not found: $CHAIN" >&2; exit 1; }

for assay in atac chip; do
  indir="${PUB}/${assay}"
  outdir="${PUB}/${assay}_hg38"
  mkdir -p "$outdir"
  echo "=== ${assay}: ${indir} -> ${outdir} ==="

  # bigWig signal tracks
  for bw in "${indir}"/*.bw; do
    [[ -e "$bw" ]] || continue
    out="${outdir}/$(basename "$bw")"
    if [[ -s "$out" ]]; then echo "  skip (exists): $(basename "$out")"; continue; fi
    "$RSCRIPT" "$LIFT" "$CHAIN" "$SIZES" "$bw" "$out"
  done

  # narrowPeak calls
  for np in "${indir}"/*_peaks.narrowPeak.gz; do
    [[ -e "$np" ]] || continue
    base="$(basename "${np%.gz}")"                          # ..._peaks.narrowPeak
    out="${outdir}/${base}"
    if [[ -s "${out}.gz" ]]; then echo "  skip (exists): ${base}.gz"; continue; fi
    "$RSCRIPT" "$LIFT" "$CHAIN" "$SIZES" "$np" "$out"
    LC_ALL=C sort -k1,1 -k2,2n "$out" > "${out}.sorted" && mv "${out}.sorted" "$out"
    "$BGZIP" -f -c "$out" > "${out}.gz"
    "$TABIX" -f -p bed "${out}.gz"
  done
done

echo
echo "Done. hg38 tracks written to ${PUB}/atac_hg38 and ${PUB}/chip_hg38"
