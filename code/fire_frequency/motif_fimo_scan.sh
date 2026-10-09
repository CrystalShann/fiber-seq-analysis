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

# Default: skip existing outputs. Override with: sbatch script.sh --force
FORCE_REWRITE=false

for arg in "$@"; do
    case "${arg}" in
        --force)
            FORCE_REWRITE=true
            ;;
        *)
            echo "Unknown argument: ${arg}" >&2
            echo "Usage: sbatch $0 [--force]" >&2
            exit 1
            ;;
    esac
done

FIRE_fasta_dir="/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency/FIRE_region_fasta"
output_root="/project/spott/cshan/fiber-seq/macrophage_project/fire_frequency/tf_motif"

motif="/project/spott/cshan/annotations/JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.txt"
fimo="/project/spott/cshan/tools/meme-5.5.9-install/bin/fimo"

mkdir -p "${output_root}"

for group in differential background; do
    input_tsv="${FIRE_fasta_dir}/${group}_FIRE_regions_with_sequences.tsv"
    fasta="${FIRE_fasta_dir}/${group}_FIRE_regions.fa"
    fimo_dir="${output_root}/fimo_${group}"
    fimo_tsv="${fimo_dir}/fimo.tsv"

    echo "Processing ${group} FIRE regions"

    # Check existence before rebuilding FASTA or running FIMO.
    if [[ -e "${fimo_tsv}" && "${FORCE_REWRITE}" == false ]]; then
        echo "SKIPPING: ${fimo_tsv} already exists."
        echo "Use --force to overwrite."
        continue
    fi

    if [[ ! -s "${input_tsv}" ]]; then
        echo "ERROR: input TSV missing or empty: ${input_tsv}" >&2
        exit 1
    fi

    if [[ "${FORCE_REWRITE}" == true ]]; then
        echo "Force rewrite enabled for ${group}."
    fi

    echo "Rebuilding FASTA from ${input_tsv}"
    awk 'BEGIN {FS=OFS="\t"}
         NR > 1 {
             print ">" $1 ":" $2 "-" $3
             print $4
         }' "${input_tsv}" > "${fasta}"

    echo "Running FIMO..."
    "${fimo}" \
        --thresh 1e-4 \
        --max-stored-scores 10000000 \
        --oc "${fimo_dir}" \
        "${motif}" \
        "${fasta}"

    if [[ ! -s "${fimo_tsv}" ]]; then
        echo "ERROR: FIMO output missing or empty: ${fimo_tsv}" >&2
        exit 1
    fi

    echo "FIMO completed: ${fimo_tsv}"
done
