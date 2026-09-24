#!/bin/bash
#SBATCH --job-name=coaccess_spans
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --array=1-4
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/coaccess_spans_%A_%a.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/coaccess_spans_%A_%a.err

# Genome-wide aligned read spans per timepoint, from the FIRE CRAM.


# Usage:
#   sbatch 02_read_spans.sh          # all four timepoints
#   bash   02_read_spans.sh LPS_0    # one timepoint

set -uo pipefail

SAMPLES=(LPS_0 LPS_5 LPS_10 LPS_15)
if [ -n "${1:-}" ]; then
    sample_name=$1
else
    sample_name=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}
fi

FIRE_ROOT="/project/spott/lizarraga/pacbio_analysis/macrophage_project/merged_hifi_bams/FIRE"
OUT_ROOT="/project/spott/cshan/fiber-seq/macrophage_project/co-accessibility"
REF="/project/spott/reference/human/GRCh38/hg38.fa"

SAMTOOLS=/project/spott/cshan/envs/dimelo/bin/samtools
BEDTOOLS=/project/spott/cshan/envs/bedtools/bin/bedtools
BGZIP=/project/spott/cshan/envs/dimelo/bin/bgzip
TABIX=/project/spott/cshan/envs/dimelo/bin/tabix

CHROMS=$(echo chr{1..22} chrX chrY)
ncore=8

die() { echo "ERROR: $*" >&2; exit 1; }

cram="${FIRE_ROOT}/${sample_name}/${sample_name}-fire-v0.1-filtered.cram"
[ -s "$cram" ] || die "FIRE CRAM not found: $cram"
[ -s "$REF" ]  || die "reference not found: $REF"

out_dir="${OUT_ROOT}/${sample_name}"
mkdir -p "$out_dir" || die "cannot create $out_dir"
spans="${out_dir}/${sample_name}.read_spans.bed.gz"

echo "Sample     : ${sample_name}"
echo "CRAM       : ${cram}"
echo "Output     : ${spans}"
echo "Started    : $(date)"

if [ -s "$spans" ] && [ -s "${spans}.tbi" ]; then
    echo "Reusing existing read spans ($(zcat "$spans" | wc -l) rows)"
else
    # one chromosome at a time
    for chrom in $CHROMS; do
    # extract primary alignments, convert to BED, and merge all chromosomes together
    # remove secondary and supplementary alignments (-F 0x900) to avoid duplicate read names
        # 0 x 100 = secondary alignment
        # 0 x 800 = supplementary alignment
    # output BAM
        "$SAMTOOLS" view -T "$REF" -F 0x900 -@ ${ncore} -b "$cram" "$chrom" \
            # convert BAM to BED
            | "$BEDTOOLS" bamtobed -i stdin
    done | "$BGZIP" -@ ${ncore} > "$spans" || die "read span extraction failed"
    "$TABIX" -f -p bed "$spans" || die "tabix failed: $spans"
fi

# count the number of primary alignments and distinct read names
n=$(zcat "$spans" | wc -l)
n_uniq=$(zcat "$spans" | cut -f4 | sort -u -T "${SLURM_TMPDIR:-/tmp}" | wc -l)
echo "  ${n} primary alignments, ${n_uniq} distinct read names"
[ "$n" -eq "$n_uniq" ] || echo "  WARNING: $((n - n_uniq)) duplicate read names remain after -F 0x900"

echo "Done -> ${out_dir}"
echo "Finished   : $(date)"
