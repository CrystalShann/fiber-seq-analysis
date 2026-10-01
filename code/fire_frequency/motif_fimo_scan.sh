#!/bin/bash
#SBATCH --job-name=fimo_fire_scan
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=30:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/code/fire_frequency/fimo_%A_%a.out
#SBATCH --error=/project/spott/cshan/fiber-seq/code/fire_frequency/fimo_%A_%a.err

set -euo pipefail

FIRE_fasta_dir="/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency/FIRE_region_fasta"
output_root="/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency/tf_motif"

motif="/project/spott/cshan/annotations/JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.txt"
fimo="/project/spott/cshan/tools/meme-5.5.9-install/bin/fimo"

mkdir -p "${output_root}"

for group in differential background; do

    # Input TSV containing chrom, start, end, sequence
    input_tsv="${FIRE_fasta_dir}/${group}_FIRE_regions_with_sequences.tsv"

    # FASTA created from TSV
    fasta="${FIRE_fasta_dir}/${group}_FIRE_regions.fa"

    # FIMO output
    fimo_dir="${output_root}/fimo_${group}"
    fimo_tsv="${fimo_dir}/fimo.tsv"

    echo "========================================"
    echo "Processing ${group} FIRE regions"
    echo "========================================"

    # Check input TSV
    if [[ ! -s "${input_tsv}" ]]; then
        echo "ERROR: input TSV missing or empty:"
        echo "${input_tsv}"
        exit 1
    fi

    # Rebuild FASTA and rescan every group on each submission.
    # --oc below overwrites the existing FIMO output for this group.
    echo "Rebuilding FASTA from ${input_tsv}"
    awk 'BEGIN {FS=OFS="\t"}
         NR > 1 {
             print ">" $1 ":" $2 "-" $3
             print $4
         }' "${input_tsv}" > "${fasta}"

    # Run FIMO
    echo "Running FIMO..."

    "${fimo}" \
        --thresh 1e-4 \
        --max-stored-scores 10000000 \
        --oc "${fimo_dir}" \
        "${motif}" \
        "${fasta}"

    # Verify output
    if [[ ! -s "${fimo_tsv}" ]]; then
        echo "ERROR: FIMO did not create:"
        echo "${fimo_tsv}"
        exit 1
    fi

    echo "FIMO completed successfully:"
    echo "${fimo_tsv}"

done
