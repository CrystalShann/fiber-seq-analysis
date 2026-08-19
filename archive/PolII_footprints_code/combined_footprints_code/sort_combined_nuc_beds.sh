#!/bin/bash
#SBATCH --job-name=sort_nucbed
#SBATCH --time=6:00:00
#SBATCH --mem=40G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/sort_nucbed%A.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/sort_nucbed%A.err


module load python/miniforge-25.3.0
module load htslib
source activate /project/spott/cshan/envs/bedtools

BED_DIR="/project/spott/cshan/fiber-seq/results/PolII/FIRE_nuc_by_chr_combined_sample/combined_nuc_beds"
OUT_DIR="${BED_DIR}/sorted_indexed"

mkdir -p "${OUT_DIR}"
shopt -s nullglob

for infile in "${BED_DIR}"/*.bed "${BED_DIR}"/*.bed.gz; do
    if [[ "${infile}" == *_sorted.bed || "${infile}" == *_sorted.bed.gz ]]; then
        continue
    fi

    base="$(basename "${infile}")"
    base="${base%.bed.gz}"
    base="${base%.bed}"
    base="${base%.gz}"

    sorted_bed="${OUT_DIR}/${base}_sorted.bed"
    sorted_gz="${sorted_bed}.gz"

    echo "Sorting ${infile}"

    if [[ "${infile}" == *.gz ]]; then
        bgzip -dc "${infile}" | sort -k1,1 -k2,2n -k3,3n > "${sorted_bed}"
    else
        sort -k1,1 -k2,2n -k3,3n "${infile}" > "${sorted_bed}"
    fi

    bgzip -f "${sorted_bed}"
    tabix -f -p bed "${sorted_gz}"

    echo "Created ${sorted_gz}"
    echo "Created ${sorted_gz}.tbi"
done
