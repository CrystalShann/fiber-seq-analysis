#!/bin/bash
#SBATCH --job-name=LCL_coaccess_spans
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --array=1-31%6
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/LCL_coaccess_spans_%A_%a.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/LCL_coaccess_spans_%A_%a.err

# One indexed genome-wide BED6 per sample; no donor filtering, MAPQ cutoff, or
# reconstruction from fire_all segments. Same primary alignment rule as macrophages.
set -euo pipefail
OUT_ROOT=${LCL_COACCESS_ROOT:-/project/spott/cshan/fiber-seq/LCL_project/co-accessibility}
INPUTS="$OUT_ROOT/universe/sample_inputs.tsv"
REF=/project/spott/reference/human/GRCh38.p14/hg38.fa
SAMTOOLS=/project/spott/cshan/envs/dimelo/bin/samtools
BEDTOOLS=/project/spott/cshan/envs/bedtools/bin/bedtools
BGZIP=/project/spott/cshan/envs/dimelo/bin/bgzip
TABIX=/project/spott/cshan/envs/dimelo/bin/tabix
[[ -s "$INPUTS" && -s "$OUT_ROOT/universe/universe.complete" ]] || { echo "Run LCL_fire_universe.sh first" >&2; exit 1; }
[[ $(wc -l < "$INPUTS") -eq 31 ]] || { echo "Expected all31 samples" >&2; exit 1; }

if [[ -n ${1:-} ]]; then
    record=$(awk -F'\t' -v sample="$1" '$1==sample' "$INPUTS")
else
    task=${SLURM_ARRAY_TASK_ID:?Supply sample name or submit the 31-task array}
    [[ "$task" -ge 1 && "$task" -le 31 ]] || { echo "Invalid array task" >&2; exit 1; }
    record=$(sed -n "${task}p" "$INPUTS")
fi
[[ -n "$record" ]] || { echo "Sample absent from the all31 manifest" >&2; exit 1; }
IFS=$'\t' read -r sample_name cram_path peaks_path elements_path spans_path <<< "$record"
[[ -s "$cram_path" && -s "$cram_path.crai" && -s "$REF" ]] || { echo "Missing CRAM/reference" >&2; exit 1; }
mkdir -p -- "$(dirname -- "$spans_path")"

if [[ -s "$spans_path" && -s "$spans_path.tbi" && -s "$spans_path.source.tsv" ]]; then
    if [[ $(head -n 1 "$spans_path.source.tsv") == "$cram_path" && "$spans_path" -nt "$cram_path" ]]; then
        "$TABIX" -l "$spans_path" > /dev/null
        echo "Reusing primary read spans for $sample_name"
        exit 0
    fi
fi
tmpdir=$(mktemp -d "${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/lcl-spans.XXXXXX")
partial="$spans_path.partial.${SLURM_JOB_ID:-$$}.gz"
trap 'rm -rf -- "$tmpdir"; rm -f -- "$partial" "$partial.tbi"' EXIT
threads=${SLURM_CPUS_PER_TASK:-8}
export LC_ALL=C
echo "Extracting genome-wide primary spans for $sample_name"
for chrom in chr{1..22} chrX chrY; do
    "$SAMTOOLS" view -T "$REF" -F 0x900 -@ "$threads" -b "$cram_path" "$chrom" \
        | "$BEDTOOLS" bamtobed -i stdin
done | "$BGZIP" -@ "$threads" -c > "$partial"
"$TABIX" -p bed "$partial"
# Counting uses read IDs within samples; duplicate primary names must not be fused.
gzip -cd -- "$partial" | cut -f4 | sort -S 2G -T "$tmpdir" | uniq -d > "$tmpdir/duplicate_ids"
[[ ! -s "$tmpdir/duplicate_ids" ]] || { echo "Duplicate primary read IDs in $sample_name" >&2; exit 1; }
mv -- "$partial" "$spans_path"
mv -- "$partial.tbi" "$spans_path.tbi"
{ printf '%s\n' "$cram_path"; printf 'alignment_filter\tprimary (-F 0x900)\n'; } > "$spans_path.source.tsv"
echo "Read spans ready: $spans_path"
