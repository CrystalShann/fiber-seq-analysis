#!/bin/bash
#SBATCH --job-name=cre_fire_freq
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/cre_fire_freq_%A_%a.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/cre_fire_freq_%A_%a.err

# Per-read FIRE frequency at GM12878 active cCREs, Low-DNase cCREs and no-DNase control
# windows in the LCL Fiber-seq samples (GRCh38).
#   bash   cre_fire_freq_dnase.sh regions               Step 1 (once)
#   sbatch --array=1-23 cre_fire_freq_dnase.sh chrom    Step 2, one chromosome per task (chr1-22, chrX)
#   bash   cre_fire_freq_dnase.sh chrom chr21 [SAMPLE]  Step 2 for one chromosome (optionally one sample)
#   bash   cre_fire_freq_dnase.sh test                  regions + Step 2 for chr21 x first sample
# Steps 3-4 (per-region counts, combine, summary) are in cre_fire_freq_dnase.Rmd (params test: true -> test_chr21/).


# measures how often reads carry a FIRE element (a Fiber-seq call for open chromatin) at three kinds of regions
# The regions are active GM12878 cCREs, Low-DNase cCREs, and matched control windows with no DNase signal

module load R/4.4.1
set -euo pipefail

OUT_ROOT=${CRE_ACCESS_ROOT:-/project/spott/cshan/fiber-seq/LCL_project/CRE_access}
SAMPLE_SHEET=/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv
CCRE=/project/spott/cshan/annotations/LCL_public_data/SCREEN_CRE_20261009/GM12878_core_cCREs.bed
DNASE_PEAKS=/project/spott/cshan/annotations/LCL_public_data/DNase_ENCFF073ORT.bed
DNASE_BW=/project/spott/cshan/annotations/LCL_public_data/DNase.ENCFF743ULW.bigWig
REF=/project/spott/reference/human/GRCh38.p14/hg38.fa
FAI=$REF.fai
SAMTOOLS=/project/spott/cshan/envs/dimelo/bin/samtools
BEDTOOLS=/project/spott/cshan/envs/bedtools/bin/bedtools
BGZIP=/project/spott/cshan/envs/dimelo/bin/bgzip
TABIX=/project/spott/cshan/envs/dimelo/bin/tabix
CHROMS=(chr{1..22} chrX)
SETS=(active_cre low_dnase_cre no_dnase_cre_region)
SEED=42
REGION_DIR=$OUT_ROOT/regions
PER_READ_DIR=$OUT_ROOT/per_read
threads=${SLURM_CPUS_PER_TASK:-4}
export LC_ALL=C
tmpdir=$(mktemp -d "${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/cre-access.XXXXXX")
trap 'rm -rf -- "$tmpdir"' EXIT

# sample_name <TAB> FIRE-filtered CRAM <TAB> FIRE elements, from the sheet's sample_name and fire_dir columns
manifest() {
    awk -F, -v OFS='\t' '
        NR == 1 {for (i = 1; i <= NF; i++) col[$i] = i; next}
        {s = $col["sample_name"]; d = $col["fire_dir"]
         print s, d "/" s "-fire-v0.1-filtered.cram", d "/additional-outputs-v0.1/fire-peaks/" s "-v0.1-fire-elements.bed.gz"}
    ' "$SAMPLE_SHEET"
}

# Step 1: region files
make_regions() {
    local chrom_re='^chr([1-9]|1[0-9]|2[0-2]|X)$' n_cand sample cram fire
    mkdir -p "$REGION_DIR"
    awk -F'\t' '$10 == "Low-DNase"' "$CCRE" > "$REGION_DIR/low_dnase_cre.bed"
    awk -F'\t' '$10 != "Low-DNase"' "$CCRE" > "$REGION_DIR/active_cre.bed"

    # Controls: cCRE intervals (any class) moved to random positions on the same chromosome, so
    # widths follow the cCRE width distribution; never placed on any cCRE, DNase peak, assembly gap
    # (N run) or FIRE unreliable-coverage region of any sample (collapsed repeats).
    awk -v OFS='\t' -v re="$chrom_re" '$1 ~ re {print $1, $2}' "$FAI" > "$tmpdir/genome.txt"
    "$SAMTOOLS" faidx "$REF" "${CHROMS[@]}" | awk -v OFS='\t' '
        /^>/ {c = substr($1, 2); p = 0; next}
        /[Nn]/ {line = $0; off = 0
                while (match(line, /[Nn]+/)) {s = p + off + RSTART - 1; print c, s, s + RLENGTH
                                              off += RSTART - 1 + RLENGTH; line = substr(line, RSTART + RLENGTH)}}
        {p += length($0)}' | "$BEDTOOLS" merge > "$tmpdir/gaps.bed"
    while IFS=$'\t' read -r -u 3 sample cram fire; do
        gzip -cd -- "${cram%/*}/additional-outputs-v0.1/coverage/unreliable-coverage-regions.bed.gz" | cut -f1-3
    done 3< <(manifest) | sort -k1,1 -k2,2n | "$BEDTOOLS" merge > "$tmpdir/unreliable.bed"
    awk '{b += $3 - $2} END {printf "Excluded assembly gaps: %.1f Mb\n", b / 1e6}' "$tmpdir/gaps.bed"
    awk -v re="$chrom_re" '$1 ~ re {b += $3 - $2} END {printf "Excluded unreliable coverage (union of samples): %.1f Mb\n", b / 1e6}' "$tmpdir/unreliable.bed"
    cut -f1-3 "$CCRE" "$DNASE_PEAKS" "$tmpdir/gaps.bed" "$tmpdir/unreliable.bed" \
        | sort -k1,1 -k2,2n | "$BEDTOOLS" merge > "$tmpdir/exclude.bed"
    n_cand=$(awk -F'\t' -v re="$chrom_re" '$1 ~ re' "$REGION_DIR/active_cre.bed" | wc -l)
    awk -F'\t' -v OFS='\t' -v re="$chrom_re" '$1 ~ re {print $1, $2, $3}' "$CCRE" \
        | "$BEDTOOLS" sample -i stdin -n "$n_cand" -seed "$SEED" \
        | "$BEDTOOLS" shuffle -i stdin -g "$tmpdir/genome.txt" -chrom -excl "$tmpdir/exclude.bed" \
            -noOverlapping -seed "$SEED" \
        | sort -k1,1 -k2,2n \
        | "$BEDTOOLS" intersect -v -a stdin -b "$DNASE_PEAKS" \
        | "$BEDTOOLS" intersect -v -a stdin -b "$CCRE" \
        | awk -v OFS='\t' '{print $1, $2, $3, sprintf("ctrl_%06d", NR)}' > "$tmpdir/candidates.bed"

    # Signal filter: keep windows whose mean DNase signal is <= the genome-wide mean (column 5 = window mean).
    Rscript - "$DNASE_BW" "$tmpdir/candidates.bed" "$tmpdir/genome.txt" "$REGION_DIR/no_dnase_cre_region.bed" <<'EOF' | tee "$REGION_DIR/no_dnase_cre_region.log"
suppressPackageStartupMessages({library(rtracklayer); library(data.table)})
a <- commandArgs(trailingOnly = TRUE)
bw <- BigWigFile(a[1])
cand <- fread(a[2], header = FALSE, col.names = c("chrom", "start", "end", "id"))
genome <- fread(a[3], header = FALSE, col.names = c("chrom", "len"))
# The bigWig tiles every base (zeros included), so bp-weighted means are true window means.
chrom_mean <- sapply(summary(bw, which = GRanges(genome$chrom, IRanges(1, genome$len)),
                             size = 1L, type = "mean", defaultValue = 0), score)
threshold <- sum(chrom_mean * genome$len) / sum(genome$len)
win <- GRanges(cand$chrom, IRanges(cand$start + 1, cand$end))
sig <- import(bw, which = reduce(win))
hit <- findOverlaps(win, sig)
q <- queryHits(hit); s <- subjectHits(hit)
bp <- pmin(end(win)[q], end(sig)[s]) - pmax(start(win)[q], start(sig)[s]) + 1
sums <- data.table(q, v = bp * score(sig)[s])[, .(v = sum(v)), by = q]
cand[, dnase_mean := 0]
cand[sums$q, dnase_mean := sums$v / width(win)[sums$q]]
keep <- cand[dnase_mean <= threshold]
cat(sprintf("Control candidates (no cCRE / DNase-peak overlap): %d\n", nrow(cand)))
cat(sprintf("DNase signal threshold = genome-wide mean over chr1-22,X: %.5f\n", threshold))
cat(sprintf("Passed signal filter: %d (%.1f%%), of which zero signal: %d\n",
            nrow(keep), 100 * nrow(keep) / nrow(cand), keep[dnase_mean == 0, .N]))
fwrite(keep[, .(chrom, start, end, id, sprintf("%.5f", dnase_mean))], a[4], sep = "\t", col.names = FALSE)
EOF
    echo "Controls overlapping a DNase peak (must be 0): $("$BEDTOOLS" intersect -u -a "$REGION_DIR/no_dnase_cre_region.bed" -b "$DNASE_PEAKS" | wc -l)"
    echo "Controls overlapping a cCRE (must be 0):       $("$BEDTOOLS" intersect -u -a "$REGION_DIR/no_dnase_cre_region.bed" -b "$CCRE" | wc -l)"
    echo "Controls overlapping an assembly gap (must be 0):      $("$BEDTOOLS" intersect -u -a "$REGION_DIR/no_dnase_cre_region.bed" -b "$tmpdir/gaps.bed" | wc -l)"
    echo "Controls overlapping unreliable coverage (must be 0):  $("$BEDTOOLS" intersect -u -a "$REGION_DIR/no_dnase_cre_region.bed" -b "$tmpdir/unreliable.bed" | wc -l)"
    wc -l "$REGION_DIR"/*.bed
}

# Step 2: per read, its reference alignment span (primary records of the FIRE-filtered CRAM) and its
# FIRE elements, both with the read name in column 4, kept only for reads overlapping any region.
per_read() {
    local sample=$1 cram=$2 fire=$3 chrom=$4
    local d=$PER_READ_DIR/$sample
    local reads=$d/$sample.$chrom.reads.bed.gz fires=$d/$sample.$chrom.fire.bed.gz
    if [[ -s $reads && -s $fires && $reads -nt $REGION_DIR/no_dnase_cre_region.bed ]]; then return 0; fi
    mkdir -p "$d"
    "$SAMTOOLS" view -@ "$threads" -u -F 2304 "$cram" "$chrom" \
        | "$BEDTOOLS" bamtobed -i stdin | cut -f1-4 \
        | "$BEDTOOLS" intersect -sorted -u -a stdin -b "$tmpdir/any_region.bed" \
        | "$BGZIP" -c > "$reads.partial"
    "$TABIX" "$fire" "$chrom" | cut -f1-4 \
        | "$BEDTOOLS" intersect -sorted -u -a stdin -b "$tmpdir/any_region.bed" \
        | awk -F'\t' 'NR == FNR {keep[$4]; next} $4 in keep' <("$BGZIP" -cd "$reads.partial") - \
        | "$BGZIP" -c > "$fires.partial"
    mv -- "$reads.partial" "$reads"
    mv -- "$fires.partial" "$fires"
}

# Step 2 for one chromosome: per-read BEDs for every sample (Step 3 counts are in cre_fire_freq_dnase.Rmd).
run_chrom() {
    local chrom=$1 only=${2:-} set sample cram fire
    for set in "${SETS[@]}"; do
        awk -F'\t' -v c="$chrom" '$1 == c' "$REGION_DIR/$set.bed" > "$tmpdir/$set.bed"
        cut -f1-4 "$tmpdir/$set.bed" > "$tmpdir/$set.bed4"
    done
    cat "$tmpdir"/*.bed4 | cut -f1-3 | sort -k1,1 -k2,2n | "$BEDTOOLS" merge > "$tmpdir/any_region.bed"
    : > "$tmpdir/samples.txt"
    while IFS=$'\t' read -r -u 3 sample cram fire; do
        [[ -z $only || $sample == "$only" ]] || continue
        echo "[$chrom] per-read BEDs: $sample"
        per_read "$sample" "$cram" "$fire" "$chrom"
        echo "$sample" >> "$tmpdir/samples.txt"
    done 3< <(manifest)
    [[ -s $tmpdir/samples.txt ]] || { echo "No sample named '$only' in $SAMPLE_SHEET" >&2; exit 1; }
}

case ${1:-} in
    regions) make_regions ;;
    chrom)   run_chrom "${2:-${CHROMS[$(( ${SLURM_ARRAY_TASK_ID:?Give a chromosome or submit --array=1-23} - 1 ))]}}" "${3:-}" ;;
    test)    [[ -s $REGION_DIR/no_dnase_cre_region.bed ]] || make_regions
             run_chrom chr21 "$(manifest | awk 'NR == 1 {print $1}')" ;;
    *)       grep '^#   ' "$0" >&2; exit 1 ;;
esac
