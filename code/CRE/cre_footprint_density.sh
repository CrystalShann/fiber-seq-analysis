#!/bin/bash
# FiberHMM v2 footprint density at GM12878 active cCREs and Low-DNase cCREs (chr1-22, X, Y), all 31 LCL
# samples. Region files come from cre_fire_freq_dnase.sh regions.
#   bash cre_footprint_density.sh inputs       footprints/{cres.tsv.gz, samples.tsv}
#   bash cre_footprint_density.sh footprints   submits cre_summary/LCL_tf_cooccupancy.sh (count_cre_footprints.R,
#                                              one sample per task, then combine) -> footprints/cre_footprint_summary.tsv.gz

set -euo pipefail

OUT_ROOT=${CRE_ACCESS_ROOT:-/project/spott/cshan/fiber-seq/LCL_project/CRE_access}
SAMPLE_SHEET=/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_samples_31samples.csv
COACCESS_SAMPLES=/project/spott/cshan/fiber-seq/LCL_project/co-accessibility/samples.tsv
CRE_SUMMARY=/project/spott/cshan/fiber-seq/code/co-accessibility/LCL/cre_summary
SETS=(active_cre low_dnase_cre)
CHROM_RE='^chr([1-9]|1[0-9]|2[0-2]|X|Y)$'
REGION_DIR=$OUT_ROOT/regions
FP_DIR=$OUT_ROOT/footprints
export LC_ALL=C

# cCREs sorted by chrom (count_cre_footprints.R expects each chromosome's rows together) and the
# 31 samples' aligned blocks, FIRE elements and FiberHMM v2 calls from the co-access sample table
make_inputs() {
    local set
    mkdir -p "$FP_DIR"
    for set in "${SETS[@]}"; do
        awk -F'\t' -v OFS='\t' -v re="$CHROM_RE" -v set="$set" '$1 ~ re {print $4, $1, $2, $3, set, $10}' "$REGION_DIR/$set.bed"
    done | sort -t$'\t' -k2,2 -k3,3n | { printf 'CRE_ID\tchrom\tstart\tend\tset\tclass\n'; cat; } | gzip -c > "$FP_DIR/cres.tsv.gz"
    awk -F'\t' -v OFS='\t' 'NR == 1 {for (i = 1; i <= NF; i++) c[$i] = i}
        {print $c["sample_name"], $c["blocks_path"], $c["elements_path"], $c["tf_path"]}' "$COACCESS_SAMPLES" > "$FP_DIR/samples.tsv"
    diff <(tail -n +2 "$FP_DIR/samples.tsv" | cut -f1 | sort) \
         <(awk -F, 'NR == 1 {for (i = 1; i <= NF; i++) if ($i == "sample_name") k = i; next} {print $k}' "$SAMPLE_SHEET" | sort) \
        || { echo "samples.tsv does not list the samples of $SAMPLE_SHEET" >&2; exit 1; }
    gzip -cd "$FP_DIR/cres.tsv.gz" | awk -F'\t' 'NR > 1 {n[$5]++} END {for (s in n) print s, n[s]}'
    echo "Samples: $(($(wc -l < "$FP_DIR/samples.tsv") - 1))"
}

submit_footprints() {
    local n a
    n=$(($(wc -l < "$FP_DIR/samples.tsv") - 1))
    export LCL_COACCESS_ROOT=$FP_DIR LCL_TF_OUT=$FP_DIR
    a=$(sbatch --parsable --partition=caslake --qos=caslake --mem=48G --time=6:00:00 --array=1-"$n" \
        "$CRE_SUMMARY/LCL_tf_cooccupancy.sh" footprints)
    sbatch --partition=caslake --qos=caslake --mem=48G --dependency=afterok:"$a" \
        "$CRE_SUMMARY/LCL_tf_cooccupancy.sh" footprints_combine
    echo "footprints array $a"
}

case ${1:-} in
    inputs)     make_inputs ;;
    footprints) submit_footprints ;;
    *)          grep '^#   ' "$0" >&2; exit 1 ;;
esac
