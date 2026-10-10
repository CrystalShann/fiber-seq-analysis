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
#   sbatch --array=1-23 cre_fire_freq_dnase.sh chrom    Steps 2-3, one chromosome per task (chr1-22, chrX)
#   bash   cre_fire_freq_dnase.sh chrom chr21 [SAMPLE]  Steps 2-3 for one chromosome (optionally one sample)
#   bash   cre_fire_freq_dnase.sh combine [CHROM ...]   by_chrom tables -> fire_freq.<set>.tsv, then summarize.R
#   bash   cre_fire_freq_dnase.sh test                  regions + chr21 x first sample + combine, in test_chr21/
module load R/4.4.1
set -euo pipefail

CODE_DIR=/project/spott/cshan/fiber-seq/code/CRE
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
RUN_DIR=$OUT_ROOT
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

# Steps 2-3 for one chromosome: per-read BEDs for every sample, then per-region counts in R.
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
    mkdir -p "$RUN_DIR/by_chrom"

    # Step 3: n_reads = reads whose span covers the whole region; n_fire = those reads with >= 1 FIRE
    # element overlapping the region; each read counted once per region.
    Rscript - "$chrom" "$tmpdir" "$PER_READ_DIR" "$RUN_DIR/by_chrom" "$BEDTOOLS" <<'EOF'
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
chrom <- a[1]; tmp <- a[2]; per_read <- a[3]; out_dir <- a[4]; bedtools <- a[5]
samples <- readLines(file.path(tmp, "samples.txt"))

# unique (region id, read name) pairs; -a is BED4, so the read name is column 8
pairs <- function(regions, b, opt = "") {
  x <- fread(cmd = sprintf("%s intersect -sorted %s -wa -wb -a %s -b %s | cut -f4,8", bedtools, opt, regions, b),
             header = FALSE, sep = "\t", colClasses = "character")
  if (!nrow(x)) return(data.table(id = character(), read = character()))
  unique(setnames(x, c("id", "read")))
}

for (set in c("active_cre", "low_dnase_cre", "no_dnase_cre_region")) {
  bed4 <- file.path(tmp, paste0(set, ".bed4"))
  reg <- fread(file.path(tmp, paste0(set, ".bed")), header = FALSE, sep = "\t")
  out <- data.table(chrom = reg$V1, start = reg$V2, end = reg$V3, id = reg$V4)
  if (ncol(reg) >= 10) out[, class := reg$V10]
  for (s in samples) {
    f <- function(kind) file.path(per_read, s, sprintf("%s.%s.%s.bed.gz", s, chrom, kind))
    span <- pairs(bed4, f("reads"), "-f 1.0")
    fire <- span[pairs(bed4, f("fire")), on = .(id, read), nomatch = NULL]
    nr <- paste0("n_reads_", s); nf <- paste0("n_fire_", s)
    out[, (nr) := 0L][, (nf) := 0L]
    out[span[, .N, by = id], on = "id", (nr) := i.N]
    out[fire[, .N, by = id], on = "id", (nf) := i.N]
  }
  nr <- paste0("n_reads_", samples); nf <- paste0("n_fire_", samples)
  out[, n_reads := as.integer(rowSums(.SD)), .SDcols = nr]
  out[, n_fire := as.integer(rowSums(.SD)), .SDcols = nf]
  out[, fire_freq := fifelse(n_reads > 0, n_fire / n_reads, NA_real_)]
  setcolorder(out, c(intersect(c("chrom", "start", "end", "id", "class"), names(out)),
                     "n_reads", "n_fire", "fire_freq", as.vector(rbind(nr, nf))))
  bad <- sum(out$n_fire > out$n_reads) + sum(mapply(function(r, f) sum(out[[f]] > out[[r]]), nr, nf))
  cat(sprintf("[%s] %-20s regions %7d | n_reads > 0: %7d | median n_reads %5.0f | n_fire > n_reads: %d\n",
              chrom, set, nrow(out), sum(out$n_reads > 0), median(out$n_reads), bad))
  fwrite(out, file.path(out_dir, sprintf("%s.%s.tsv", set, chrom)), sep = "\t", na = "NA", quote = FALSE)
}
EOF
}

# Concatenate chromosome tables into fire_freq.<set>.tsv, then Step 4 (summarize.R).
combine() {
    local set chrom files
    for set in "${SETS[@]}"; do
        files=()
        for chrom in "$@"; do
            [[ -s $RUN_DIR/by_chrom/$set.$chrom.tsv ]] || { echo "Missing $RUN_DIR/by_chrom/$set.$chrom.tsv" >&2; exit 1; }
            files+=("$RUN_DIR/by_chrom/$set.$chrom.tsv")
        done
        [[ $(head -qn1 "${files[@]}" | sort -u | wc -l) -eq 1 ]] \
            || { echo "$set: chromosome tables have different sample columns" >&2; exit 1; }
        { head -n1 "${files[0]}"; tail -qn +2 "${files[@]}"; } > "$RUN_DIR/fire_freq.$set.tsv"
    done
    echo "Combined $# chromosome(s) into $RUN_DIR/fire_freq.<set>.tsv"
    Rscript "$CODE_DIR/summarize.R" "$RUN_DIR"
}

case ${1:-} in
    regions) make_regions ;;
    chrom)   run_chrom "${2:-${CHROMS[$(( ${SLURM_ARRAY_TASK_ID:?Give a chromosome or submit --array=1-23} - 1 ))]}}" "${3:-}" ;;
    combine) shift; if (( $# )); then combine "$@"; else combine "${CHROMS[@]}"; fi ;;
    test)    RUN_DIR=$OUT_ROOT/test_chr21
             [[ -s $REGION_DIR/no_dnase_cre_region.bed ]] || make_regions
             run_chrom chr21 "$(manifest | awk 'NR == 1 {print $1}')"
             combine chr21 ;;
    *)       grep '^#   ' "$0" >&2; exit 1 ;;
esac
