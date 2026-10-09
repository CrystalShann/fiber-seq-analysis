#!/bin/bash
# Genome-wide cCRE universe from FIRE peaks across all 31 LCL samples.
# Run before LCL_read_spans.sh and LCL_run_coaccess.sh.
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PYTHON=/project/spott/cshan/envs/Jupyter-notebook/bin/python3
BEDTOOLS=/project/spott/cshan/envs/bedtools/bin/bedtools
BGZIP=/project/spott/cshan/envs/dimelo/bin/bgzip
OUT_ROOT=${LCL_COACCESS_ROOT:-/project/spott/cshan/fiber-seq/LCL_project/co-accessibility}
CRE_BED=/project/spott/cshan/annotations/GRCh38-cCREs.bed
TSS_BED=/project/spott/cshan/annotations/gencodev46_Ensembl_canonical_TSS.bed
WINDOW=10000
CHROM_RE='^chr([1-9]|1[0-9]|2[0-2]|X|Y)$'

[[ -s "$CRE_BED" && -s "$TSS_BED" ]] || { echo "Missing annotations; see ../make_gencode_v46_all_tss.sh" >&2; exit 1; }
"$PYTHON" "$SCRIPT_DIR/LCL_coaccess_cres.py" --root "$OUT_ROOT" --prepare-manifest
uni="$OUT_ROOT/universe"
tmpdir=$(mktemp -d "${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/lcl-universe.XXXXXX")
trap 'rm -rf -- "$tmpdir"' EXIT

# Peak selection retains macrophage's >=1bp rule. The 50% rule applies later to
# per-read FIRE ELEMENTS, not to this population-level peak universe.
: > "$tmpdir/peaks.bed"
while IFS=$'\t' read -r sample_name cram_path peaks_path elements_path spans_path; do
    echo "FIRE peaks: $sample_name"
    gzip -cd -- "$peaks_path" | awk -v re="$CHROM_RE" 'BEGIN{OFS="\t"} $1 ~ re {print $1,$2,$3}' \
        | LC_ALL=C sort -k1,1 -k2,2n -T "$tmpdir" > "$tmpdir/$sample_name.peaks.bed"
    cat "$tmpdir/$sample_name.peaks.bed" >> "$tmpdir/peaks.bed"
done < "$uni/sample_inputs.tsv"
LC_ALL=C sort -k1,1 -k2,2n -T "$tmpdir" "$tmpdir/peaks.bed" \
    | "$BEDTOOLS" merge -i stdin > "$tmpdir/fire_peaks_union.bed"

awk -v re="$CHROM_RE" 'BEGIN{OFS="\t"} $1 ~ re {print $1,$2,$3,$4,$5,$6}' "$CRE_BED" \
    | "$BEDTOOLS" intersect -a stdin -b "$tmpdir/fire_peaks_union.bed" -u \
    | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4"."$5"."$6,$6}' \
    | LC_ALL=C sort -k1,1 -k2,2n -T "$tmpdir" > "$tmpdir/cre_universe.bed"
[[ -s "$tmpdir/cre_universe.bed" ]] || { echo "No cCRE overlaps LCL FIRE peaks" >&2; exit 1; }

awk -F'\t' -v OFS='\t' -v W="$WINDOW" -v re="$CHROM_RE" '
    $1 ~ re {s=$2-W; if(s<0)s=0; split($4,a,";"); print $1,s,$2+W,a[1],a[3],a[4],$6}' "$TSS_BED" \
    | LC_ALL=C sort -k1,1 -k2,2n -T "$tmpdir" > "$tmpdir/gene_windows.bed"
{
    printf 'CRE_ID\tchrom\tstart\tend\tCRE_label\tgene_id\tgene_name\ttranscript_type\n'
    "$BEDTOOLS" intersect -a "$tmpdir/cre_universe.bed" -b "$tmpdir/gene_windows.bed" -wa -wb \
        | awk 'BEGIN{OFS="\t"} {print $4,$1,$2,$3,$5,$9,$10,$11}'
} | "$BGZIP" -c > "$tmpdir/cre_gene_map.tsv.gz"

# Sample peak membership is metadata, never a per-sample testing restriction.
cut -f4 "$tmpdir/cre_universe.bed" > "$tmpdir/ids"
flag_files=("$tmpdir/ids")
printf 'CRE_ID' > "$tmpdir/flags_header"
while IFS=$'\t' read -r sample_name cram_path peaks_path elements_path spans_path; do
    "$BEDTOOLS" intersect -a "$tmpdir/cre_universe.bed" -b "$tmpdir/$sample_name.peaks.bed" -c \
        | awk '{print ($6>0)?1:0}' > "$tmpdir/$sample_name.flags"
    flag_files+=("$tmpdir/$sample_name.flags")
    printf '\tin_%s_peaks' "$sample_name" >> "$tmpdir/flags_header"
done < "$uni/sample_inputs.tsv"
printf '\n' >> "$tmpdir/flags_header"
{ cat "$tmpdir/flags_header"; paste "${flag_files[@]}"; } | "$BGZIP" -c > "$tmpdir/cre_in_sample_peaks.tsv.gz"

for name in fire_peaks_union.bed cre_universe.bed gene_windows.bed cre_gene_map.tsv.gz cre_in_sample_peaks.tsv.gz; do
    cp -- "$tmpdir/$name" "$uni/$name.partial"
    mv -- "$uni/$name.partial" "$uni/$name"
done
printf 'Completed all31 LCL universe\n' > "$uni/universe.complete"
echo "LCL universe ready: $uni"
wc -l "$uni/cre_universe.bed" "$uni/gene_windows.bed"
