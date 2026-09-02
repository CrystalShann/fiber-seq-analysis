#!/bin/bash

# One row per TRANSCRIPT; the interval is 20 bp centered on the
# TSS (matching the example's "TSS coordinates are 20 bp centered on predicted
# TSS", which is why it pairs with bedmap --range 240 for "within 250 bp").
#
# BED6 layout:
#   chrom  tss-10  tss+10  gene_id;transcript_id;gene_name;transcript_type;tags  .  strand
#
# The GTF's transcript tags (Ensembl_canonical, MANE_Select, basic, ...) are
# comma-joined inside column 4, so the example's
#   grep Ensembl_canonical all_tss.bed > canonical_TSS.bed
# works unchanged. Both the all-TSS bed and that canonical subset are written.
#
# Usage:  bash make_gencode_v46_all_tss.sh


# GTF columns:
# 1 chromosome
# 2 source
# 3 feature
# 4 start
# 5 end
# 6 score
# 7 strand
# 8 frame
# 9 attributes

# column 9 contain attributes
# gid    = gene ID - gene_id "ENSG000001";
# tid    = transcript ID - transcript_id "ENST000001";
# gname  = gene name - gene_name "ABC";
# ttype  = transcript type - transcript_type "protein_coding";
# tags   = transcript tags - tag "basic"; tag "Ensembl_canonical";

set -euo pipefail

OUT_DIR="/project/spott/cshan/annotations"
GTF="${OUT_DIR}/gencode.v46.annotation.gtf.gz"
URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz"

all_tss="${OUT_DIR}/gencode.v46.annotation_all_tss.bed"
canonical="${OUT_DIR}/gencodev46_Ensembl_canonical_TSS.bed"

# donwload the gtf file if it does not exist
if [ ! -s "$GTF" ]; then
    echo "Downloading GENCODE v46 GTF..."
    curl -sSf -o "$GTF" "$URL"
fi

echo "Building ${all_tss}..."


##############################
# Create bed6 file from GTF file
##############################
# 1. chr
# 2. start
# 3. end
# 4. annotation information: gene_id;transcript_id;gene_name;transcript_type;tags
# 5. score .
# 6. strand

# output: /project/spott/cshan/annotation/gencode.v46.annotation_all_tss.bed

# decompress the gzipped GTF
# input fields are separated by tabs
# when awk prints multiple fields, separate them using tabs
#
# if column 3 is not a transcript, skip this row and proceed to the next row
#
##############################
# calculate strand aware TSS
##############################
# converts GTF 1 based to BED 0 based
# for + strand: TSS is at GTF start ----> TSS = start -1
# for - strand, TSS is at the GTF end --> TSS = end - 1
# create 20 bp TSS window: TSS window start = TSS start - 10; TSS window end = TSS start + 10
#
##############################
# initialize annotation variables
##############################
# split column9 of GTF to create an array
# loop through every attribute
# collect all tags: a transcript can have multiple tags, and concataenates all tag
# if no tag is present, use NA
#   tag "basic";
#   tag "Ensembl_canonical";
#   tag "MANE_Select";
# remove attribute name and quotation marks
# sort by the bed6 file
#   -k1,1 - sort by chr
#   -k2,2n - numerically by BED start
#   -k3,3n - numerically by BED end

zcat "$GTF" \
    | awk -F'\t' -v OFS='\t' '
        $3 != "transcript" { next }
        {
            # TSS in 0-based coords: GTF start ($4) is 1-based; for minus-strand
            # transcripts the TSS is the GTF end ($5).
            tss = ($7 == "+") ? $4 - 1 : $5 - 1
            s = tss - 10; if (s < 0) s = 0
            e = tss + 10

            gid = tid = gname = ttype = ""; tags = ""
            n = split($9, a, ";")
            for (i = 1; i <= n; i++) {
                f = a[i]
                gsub(/^ +| +$/, "", f)
                if (f ~ /^gene_id /)              { gid   = f }
                else if (f ~ /^transcript_id /)   { tid   = f }
                else if (f ~ /^gene_name /)       { gname = f }
                else if (f ~ /^transcript_type /) { ttype = f }
                else if (f ~ /^tag /)             { tags = (tags == "" ? "" : tags ",") f }
                gsub(/^[a-z_]+ "|"$/, "", gid); gsub(/^[a-z_]+ "|"$/, "", tid)
                gsub(/^[a-z_]+ "|"$/, "", gname); gsub(/^[a-z_]+ "|"$/, "", ttype)
            }
            gsub(/tag "|"/, "", tags)
            if (tags == "") tags = "NA"
            print $1, s, e, gid ";" tid ";" gname ";" ttype ";" tags, ".", $7
        }' \
    | LC_ALL=C sort -k1,1 -k2,2n -k3,3n > "$all_tss"

##############################
# extract canonical transcripts from bed6
##############################

# output: /project/spott/cshan/annotations/gencodev46_Ensembl_canonical_TSS.bed

echo "Filtering to Ensembl_canonical -> ${canonical}..."
grep Ensembl_canonical "$all_tss" > "$canonical"

## Sanity checks
n_all=$(wc -l < "$all_tss")
n_can=$(wc -l < "$canonical")
n_bad=$(awk -F'\t' '$3 - $2 != 20' "$all_tss" | wc -l)
echo "transcript TSSs : ${n_all}"
echo "canonical TSSs  : ${n_can}"
[ "$n_bad" -eq 0 ] || { echo "ERROR: ${n_bad} rows are not 20 bp wide" >&2; exit 1; }
# every canonical row should be a distinct gene
n_can_genes=$(cut -d';' -f1 "$canonical" | cut -f4 | sort -u | wc -l)
echo "canonical genes : ${n_can_genes}"
