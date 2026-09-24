#!/bin/bash
#SBATCH --job-name=fire_metagene
#SBATCH --time=24:00:00
#SBATCH --mem=180G
#SBATCH --cpus-per-task=16
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --array=0-33
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/fire_metagene_%A_%a.out

module load python/miniforge-25.3.0
# source activate /project/spott/cshan/envs/pb-CpG-tools
source activate /project/spott/cshan/envs/deeptools

FIRE_dir="/project/spott/1_Shared_projects/LCL_Fiber_seq/FIRE/results/"
ref_genome="/project/spott/reference/human/GRCh38/hg38.fa"
gtf="/project/spott/cshan/annotations/gencode.v49.primary_assembly.annotation.gtf"
# cpg_output_base="/project/spott/cshan/fiber-seq/results/fire_CpG"
plot_output_dir="/project/spott/cshan/fiber-seq/results/plots/metagene_heatmaps"
threads=16

# Create output directories
mkdir -p ${cpg_output_base}
mkdir -p ${plot_output_dir}
mkdir -p /project/spott/cshan/fiber-seq/results/logs

# Get list of all sample directories
mapfile -t sample_dirs < <(find ${FIRE_dir} -mindepth 1 -maxdepth 1 -type d | sort)

# Get the sample directory for this array task
sample_dir="${sample_dirs[$SLURM_ARRAY_TASK_ID]}"
sample=$(basename ${sample_dir})

echo "========================================="
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Processing sample: ${sample}"
echo "Sample directory: ${sample_dir}"
echo "========================================="

# Find the CRAM file
# cram_file=$(find ${sample_dir} -maxdepth 1 -name "*-fire-v0.1-filtered.cram" -o -name "*-fire-v0.1.1-filtered.cram" | head -n 1)

# if [ -z "${cram_file}" ]; then
#   echo "ERROR: No CRAM file found in ${sample_dir}"
#    exit 1
# fi

# echo "Found CRAM: ${cram_file}"

# Set up CpG output directory for this sample
# cpg_output_dir="${cpg_output_base}/${sample}"
# mkdir -p ${cpg_output_dir}

# Run aligned_bam_to_cpg_scores if not already done
#if [ ! -f "${cpg_output_dir}/${sample}_CPG.combined.bw" ]; then
#    echo "Running aligned_bam_to_cpg_scores..."
    
#    aligned_bam_to_cpg_scores \
#         --bam ${cram_file} \
#         --output-prefix ${cpg_output_dir}/${sample}_CPG \
#         --ref ${ref_genome} \
#         --threads ${threads}
    
#     echo "CpG scores complete for ${sample}"
# else
#     echo "CpG scores already exist, skipping..."
# fi

# Collect all BigWig files (FIRE trackHub + CpG outputs)
bw_dir="${sample_dir}/trackHub-v0.1/bw"

# Create list of all BigWig files to process
bw_files=()

# Add FIRE BigWigs
if [ -d "${bw_dir}" ]; then
    for bw in ${bw_dir}/*.bw; do
        [ -f "${bw}" ] && bw_files+=("${bw}")
    done
fi

# # Add CpG BigWigs
# for bw in ${cpg_output_dir}/*.bw; do
#     [ -f "${bw}" ] && bw_files+=("${bw}")
# done

# Check if we have any BigWig files
if [ ${#bw_files[@]} -eq 0 ]; then
    echo "ERROR: No BigWig files found for ${sample}"
    exit 1
fi

echo "Found ${#bw_files[@]} BigWig files to process"


# Process each BigWig file
for bw in "${bw_files[@]}"; do
    label=$(basename ${bw} .bw)
    echo "  Processing ${label}..."
    
    # Compute matrix around genes 
    computeMatrix scale-regions \
        -S ${bw} \
        -R ${gtf} \
        --beforeRegionStartLength 3000 \
        --regionBodyLength 5000 \
        --afterRegionStartLength 3000 \
        --skipZeros \
        -p ${threads} \
        -o ${plot_output_dir}/${sample}_${label}.matrix.gz
    
    # Plot heatmap
    plotHeatmap -m ${plot_output_dir}/${sample}_${label}.matrix.gz \
        --plotFileFormat pdf \
        --samplesLabel ${label} \
        --heatmapWidth 8 \
        -out ${plot_output_dir}/${sample}_${label}_heatmap.pdf
    
    echo "  Created heatmap: ${plot_output_dir}/${sample}_${label}_heatmap.pdf"
    
    # # Plot profile
    # plotProfile -m ${plot_output_dir}/${sample}_${label}.matrix.gz \
    #     --plotFileFormat pdf \
    #     --samplesLabel ${label} \
    #     --perGroup \
    #     -out ${plot_output_dir}/${sample}_${label}_profile.pdf
    
    # echo "  Created profile: ${plot_output_dir}/${sample}_${label}_profile.pdf"
done

echo "========================================="
echo "Completed ${sample}"
# echo "CpG BigWigs saved to: ${cpg_output_dir}"
echo "Plots saved to: ${plot_output_dir}"
echo "========================================="