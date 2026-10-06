#!/bin/bash
#SBATCH --job-name=find_msp_accessible_enhancers
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --chdir=/project/spott/cshan/fiber-seq
#SBATCH --output=/project/spott/cshan/fiber-seq/code/enhancer/logs/enhancer_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/code/enhancer/logs/enhancer_%j.err

# Submit: sbatch code/enhancer/shell_scripts/find_accessible_enhancers.sh
# Accessible: ONE MSP block overlaps strictly >50% of the enhancer length.
set -euo pipefail
export LC_ALL=C
[[ $# -eq 0 ]] || { echo "Configure paths below; no command-line options are used." >&2; exit 1; }

project=/project/spott/cshan/fiber-seq
enhancers=${ENHANCERS:-$project/macrophage_project/annotations/enhancer_public_data/fantom5_enhancer_hg38/F5.hg38.enhancers.bed.gz}
msp_root=${MSP_ROOT:-$project/macrophage_project/FiberHMM/extract/firehmm_msp}
out=${OUTPUT_DIR:-$project/macrophage_project/enhancer/TF_co-occ/tables}
bedtools=/project/spott/cshan/envs/bedtools/bin/bedtools
tabix=/project/spott/cshan/envs/dimelo/bin/tabix
samples=(LPS_0 LPS_5 LPS_10 LPS_15)

mkdir -p "$out"
tmp=$(mktemp -d "$out/.enhancer_accessibility.XXXXXX")
trap 'rm -rf -- "$tmp"' EXIT

# Keep enhancer outer coordinates on standard chromosomes, in input order.
if [[ "$enhancers" == *.gz ]]; then gzip -cd -- "$enhancers"; else cat -- "$enhancers"; fi |
awk 'BEGIN { FS=OFS="\t" }
     $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ {
         if (NF < 4 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $3 <= $2 || $4 == "") exit 1
         if (seen[$1,$2,$3]++) next
         if (ids[$4]++) exit 1
         print $1,$2,$3,$4
     }' > "$tmp/enhancers.bed"
[[ -s "$tmp/enhancers.bed" ]] || { echo "No enhancers found." >&2; exit 1; }
cut -f1 "$tmp/enhancers.bed" | sort -u > "$tmp/chromosomes"
while read -r chrom; do
    awk -F '\t' -v c="$chrom" '$1 == c' "$tmp/enhancers.bed" |
        sort -k2,2n -k3,3n > "$tmp/$chrom.bed"
    "$bedtools" merge -i "$tmp/$chrom.bed" > "$tmp/$chrom.query.bed"
done < "$tmp/chromosomes"

# Query nearby fibers, expand each individual MSP block, then test its overlap.
process_sample() {
    trap - EXIT
    local sample=$1 label=${1//_/} chrom input
    : > "$tmp/$sample.ids"
    while read -r chrom; do
        input="$msp_root/$label/${label}_hmm_extracted_msp_$chrom.bed.gz"
        "$tabix" --cache 64 -R "$tmp/$chrom.query.bed" "$input" |
        awk 'BEGIN { FS=OFS="\t" }
             {
                 sub(/,$/, "", $11); sub(/,$/, "", $12)
                 n=split($11, sizes, ","); m=split($12, offsets, ",")
                 if (NF < 12 || n != $10 || m != $10) exit 1
                 for (i=1; i<=n; i++) {
                     start=$2+offsets[i]; end=start+sizes[i]
                     if (sizes[i] !~ /^[0-9]+$/ || sizes[i] <= 0 ||
                         offsets[i] !~ /^[0-9]+$/ || end > $3) exit 1
                     print $1,start,end
                 }
             }' |
        "$bedtools" intersect -a stdin -b "$tmp/$chrom.bed" -F 0.5 -wo |
        awk 'BEGIN { FS="\t" } $8*2 > $6-$5 && !seen[$7]++ { print $7 }' >> "$tmp/$sample.ids"
        # B is the enhancer; the final comparison excludes exactly 50% overlap.
    done < "$tmp/chromosomes"
    echo "Finished: $sample"
}

# One worker per timepoint. Wait for all four before publishing any results.
pids=()
for sample in "${samples[@]}"; do
    process_sample "$sample" &
    pids+=("$!")
done
failed=0
for pid in "${pids[@]}"; do wait "$pid" || failed=1; done
[[ $failed -eq 0 ]] || { echo "MSP processing failed; previous tables retained." >&2; exit 1; }

# Write per-timepoint flags, their overall OR, and accessible enhancer counts.
awk -v samples="${samples[*]}" -v tmp="$tmp" '
BEGIN {
    FS=OFS="\t"; n=split(samples,s," ")
    printf "enhancer_id\tchr\tstart0\tend0"
    for (j=1; j<=n; j++) {
        printf "\t%s",s[j]; path=tmp "/" s[j] ".ids"
        while ((getline id < path)>0) hit[j,id]=1
        close(path)
    }
    printf "\taccessible\n"
}
{
    printf "%s\t%s\t%s\t%s",$4,$1,$2,$3; any=0
    for (j=1; j<=n; j++) {
        yes=((j SUBSEP $4) in hit); any+=yes; counts[j]+=yes
        printf "\t%s",yes ? "TRUE":"FALSE"
    }
    printf "\t%s\n",any ? "TRUE":"FALSE"
}
END {
    path=tmp "/01_accessible_counts.tsv"
    print "sample","n_accessible" > path
    for (j=1; j<=n; j++) print s[j],counts[j]+0 > path
}' "$tmp/enhancers.bed" > "$tmp/01_enhancer_accessibility.tsv"
mv -- "$tmp/01_enhancer_accessibility.tsv" "$tmp/01_accessible_counts.tsv" "$out/"
echo "Wrote accessibility tables: $out"
