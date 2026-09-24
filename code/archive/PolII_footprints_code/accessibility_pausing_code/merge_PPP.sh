#!/bin/bash

#SBATCH --job-name=merge_PPP
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --time=12:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/merge_PPP_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/merge_PPP_%j.err


INPUT_DIR="/project/spott/cshan/fiber-seq/results/PolII/footprint_summary_beds/PPP"
OUTPUT_FILE="${INPUT_DIR}/merged_PPP_footprint_summary.tsv.gz"

echo "Starting merge at $(date)"

cd "$INPUT_DIR"

# create header
echo -e "chrom\tstart\tend\tread_id\tclass\tmidpoint" | gzip > "$OUTPUT_FILE"

# append all files (skip headers since BED doesn't have one)
for f in *.bed.gz; do
    echo "Processing $f"
    zcat "$f" >> temp_merged.tsv
done

# compress final
gzip -c temp_merged.tsv >> "$OUTPUT_FILE"

# cleanup
rm temp_merged.tsv

echo "Finished at $(date)"
