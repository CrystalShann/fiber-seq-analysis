#!/bin/bash
#SBATCH --job-name=sort_parsed_nuc
#SBATCH --time=20:00:00
#SBATCH --mem=100G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/sort_parsed_nuc%A.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/sort_parsed_nuc%A.err

module load python/miniforge-25.3.0
module load htslib
source activate /project/spott/cshan/envs/bedtools

BED_DIR="/project/spott/cshan/fiber-seq/results/PolII/FIRE_nuc_by_chr_combined_sample/parsed_bed_files"
OUT_DIR="${BED_DIR}/sorted_indexed"

mkdir -p "${OUT_DIR}"
shopt -s nullglob

for infile in "${BED_DIR}"/*_parsed.bed; do
    base="$(basename "${infile}" .bed)"
    sorted_gz="${OUT_DIR}/${base}_sorted.bed.gz"

    echo "Sorting ${infile}"
    sort -k1,1 -k2,2n -k3,3n "${infile}" | bgzip -f > "${sorted_gz}"
    tabix -f -p bed "${sorted_gz}"

    echo "Created ${sorted_gz}"
    echo "Created ${sorted_gz}.tbi"
done
