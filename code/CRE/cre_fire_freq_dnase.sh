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

# Per-read FIRE frequency at GM12878 active cCREs and Low-DNase cCREs in the LCL Fiber-seq
# samples (GRCh38).
#   bash   cre_fire_freq_dnase.sh regions               Step 1 (once)
#   sbatch --array=1-24 cre_fire_freq_dnase.sh chrom    Step 2, one chromosome per task (chr1-22, chrX, chrY)
#   bash   cre_fire_freq_dnase.sh chrom chr21 [SAMPLE]  Step 2 for one chromosome (optionally one sample)
#   bash   cre_fire_freq_dnase.sh test                  regions + Step 2 for chr21 x first sample
# Step 3 (per-region counts and combine): count_cre_fire_frequency.R [--test].
# Step 4 (read existing counts, summaries and plots): cre_fire_freq_dnase.Rmd (params test: true -> test_chr21/).


# measures how often reads carry a FIRE element (a Fiber-seq call for open chromatin) at two kinds of regions
# The regions are active GM12878 cCREs and Low-DNase cCREs

set -euo pipefail

OUT_ROOT=${CRE_ACCESS_ROOT:-/project/spott/cshan/fiber-seq/LCL_project/CRE_access}
SAMPLE_SHEET=/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv
CCRE=/project/spott/cshan/annotations/LCL_public_data/SCREEN_CRE_20261009/GM12878_core_cCREs.bed
SAMTOOLS=/project/spott/cshan/envs/dimelo/bin/samtools
BEDTOOLS=/project/spott/cshan/envs/bedtools/bin/bedtools
BGZIP=/project/spott/cshan/envs/dimelo/bin/bgzip
TABIX=/project/spott/cshan/envs/dimelo/bin/tabix
FIBERHMM=/project/spott/1_Shared_projects/LCL_Fiber_seq/FiberHMM/FiberHMM_v2/results/results_v1_trained_model
CHROMS=(chr{1..22} chrX chrY)
SETS=(active_cre low_dnase_cre)
REGION_DIR=$OUT_ROOT/regions
PER_READ_DIR=$OUT_ROOT/per_read
threads=${SLURM_CPUS_PER_TASK:-4}
export LC_ALL=C
tmpdir=$(mktemp -d "${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/cre-access.XXXXXX")
trap 'rm -rf -- "$tmpdir"' EXIT

# sample_name <TAB> FIRE-filtered CRAM <TAB> FIRE elements <TAB> FiberHMM v2 calls, from the sheet's
# sample_name and fire_dir columns
manifest() {
    awk -F, -v OFS='\t' -v hmm="$FIBERHMM" '
        NR == 1 {for (i = 1; i <= NF; i++) col[$i] = i; next}
        {s = $col["sample_name"]; d = $col["fire_dir"]
         print s, d "/" s "-fire-v0.1-filtered.cram", d "/additional-outputs-v0.1/fire-peaks/" s "-v0.1-fire-elements.bed.gz",
               hmm "/" s "/" s ".recalled_tf.bed.gz"}
    ' "$SAMPLE_SHEET"
}

# Step 1: region files
make_regions() {
    mkdir -p "$REGION_DIR"
    awk -F'\t' '$10 == "Low-DNase"' "$CCRE" > "$REGION_DIR/low_dnase_cre.bed"
    awk -F'\t' '$10 != "Low-DNase"' "$CCRE" > "$REGION_DIR/active_cre.bed"
    wc -l "$REGION_DIR"/*.bed
}

# Step 2: per read, its aligned blocks (one line per block of the primary records of the FIRE-filtered CRAM,
# split at deletions and N gaps) and its FIRE elements, both with the read name in column 4, kept only
# for reads with a FiberHMM v2 record on the chromosome (as n_covered in cre_summary/count_cre_footprints.R),
# for blocks overlapping any region and those blocks' reads.
per_read() {
    local sample=$1 cram=$2 fire=$3 tf=$4 chrom=$5
    local d=$PER_READ_DIR/$sample
    local reads=$d/$sample.$chrom.reads.bed.gz fires=$d/$sample.$chrom.fire.bed.gz
    if [[ -s $reads && -s $fires && $reads -nt $REGION_DIR/active_cre.bed ]]; then return 0; fi
    mkdir -p "$d"
    "$TABIX" "$tf" "$chrom" | cut -f4 > "$tmpdir/fiberhmm_reads.txt"
    "$SAMTOOLS" view -@ "$threads" -u -F 2304 "$cram" "$chrom" \
        | "$BEDTOOLS" bamtobed -splitD -i stdin | cut -f1-4 \
        | awk -F'\t' 'NR == FNR {keep[$1]; next} $4 in keep' "$tmpdir/fiberhmm_reads.txt" - \
        | sort -T "$tmpdir" -k2,2n \
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
    local chrom=$1 only=${2:-} set sample cram fire tf
    for set in "${SETS[@]}"; do
        awk -F'\t' -v c="$chrom" '$1 == c' "$REGION_DIR/$set.bed" > "$tmpdir/$set.bed"
        cut -f1-4 "$tmpdir/$set.bed" > "$tmpdir/$set.bed4"
    done
    cat "$tmpdir"/*.bed4 | cut -f1-3 | sort -k1,1 -k2,2n | "$BEDTOOLS" merge > "$tmpdir/any_region.bed"
    : > "$tmpdir/samples.txt"
    while IFS=$'\t' read -r -u 3 sample cram fire tf; do
        [[ -z $only || $sample == "$only" ]] || continue
        echo "[$chrom] per-read BEDs: $sample"
        per_read "$sample" "$cram" "$fire" "$tf" "$chrom"
        echo "$sample" >> "$tmpdir/samples.txt"
    done 3< <(manifest)
    [[ -s $tmpdir/samples.txt ]] || { echo "No sample named '$only' in $SAMPLE_SHEET" >&2; exit 1; }
}

case ${1:-} in
    regions) make_regions ;;
    chrom)   run_chrom "${2:-${CHROMS[$(( ${SLURM_ARRAY_TASK_ID:?Give a chromosome or submit --array=1-24} - 1 ))]}}" "${3:-}" ;;
    test)    [[ -s $REGION_DIR/active_cre.bed ]] || make_regions
             run_chrom chr21 "$(manifest | awk 'NR == 1 {print $1}')" ;;
    *)       grep '^#   ' "$0" >&2; exit 1 ;;
esac
