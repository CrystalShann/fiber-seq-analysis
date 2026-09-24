#!/bin/bash
#SBATCH --job-name=fire_metagene
#SBATCH --time=24:00:00
#SBATCH --mem=180G
#SBATCH --cpus-per-task=16
#SBATCH --account=pi-spott
#SBATCH --partition=spott
#SBATCH --array=0-33
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/fire_metagene_%A_%a.out

module load python/miniforge-25.3.0
source activate /project/spott/cshan/envs/deeptools

FIRE_dir="/project/spott/1_Shared_projects/LCL_Fiber_seq/FIRE/results/"
gtf="/project/spott/cshan/annotations/gencode.v49.primary_assembly.annotation.gtf"
plot_output_dir="/project/spott/cshan/fiber-seq/results/plots/metagene_v1"
dnase_bw="/project/spott/cshan/annotations/DNase.ENCFF743ULW.bigWig"
threads=16


mkdir -p ${plot_output_dir}

# Define the sample_dirs array - populate with all sample directories
mapfile -t sample_dirs < <(find ${FIRE_dir} -mindepth 1 -maxdepth 1 -type d | sort)

# Debug: Print total number of samples found
echo "Total samples found: ${#sample_dirs[@]}"
echo "Current array task ID: ${SLURM_ARRAY_TASK_ID}"

# Get sample info
sample_dir="${sample_dirs[$SLURM_ARRAY_TASK_ID]}"
sample=$(basename "${sample_dir}")
bw_dir="${sample_dir}/trackHub-v0.1/bw"

echo "========================================="
echo "Processing sample: ${sample}"
echo "Sample directory: ${sample_dir}"
echo "BigWig directory: ${bw_dir}"
echo "========================================="

# Verify BigWig directory exists
if [ ! -d "${bw_dir}" ]; then
    echo "ERROR: BigWig directory does not exist: ${bw_dir}"
    exit 1
fi

# Arrays to store BigWig files and their labels
declare -a all_bw_files
declare -a all_labels

# Add tracks in order: FIRE peaks, DNase-seq, Linker DNA, Nucleosomes
echo "Checking for BigWig files..."

if [ -f "${bw_dir}/all.fire.coverage.bw" ]; then
    all_bw_files+=("${bw_dir}/all.fire.coverage.bw")
    all_labels+=("FIRE_peaks")
    echo "Found FIRE peaks"
else
    echo "Missing: ${bw_dir}/all.fire.coverage.bw"
fi

if [ -f "${dnase_bw}" ]; then
    all_bw_files+=("${dnase_bw}")
    all_labels+=("DNase-seq")
    echo "Found DNase-seq"
else
    echo "Missing: ${dnase_bw}"
fi

if [ -f "${bw_dir}/all.linker.coverage.bw" ]; then
    all_bw_files+=("${bw_dir}/all.linker.coverage.bw")
    all_labels+=("Linker_DNA")
    echo "Found Linker DNA"
else
    echo "Missing: ${bw_dir}/all.linker.coverage.bw"
fi

if [ -f "${bw_dir}/all.nucleosome.coverage.bw" ]; then
    all_bw_files+=("${bw_dir}/all.nucleosome.coverage.bw")
    all_labels+=("Nucleosomes")
    echo "Found Nucleosomes"
else
    echo "Missing: ${bw_dir}/all.nucleosome.coverage.bw"
fi

# Check if we have all required files
if [ ${#all_bw_files[@]} -lt 4 ]; then
    echo ""
    echo "ERROR: Missing required BigWig files for ${sample}"
    echo "Found ${#all_bw_files[@]} files, need 4 (FIRE, DNase, Linker, Nucleosome)"
    echo ""
    echo "Available files in ${bw_dir}:"
    ls -lh "${bw_dir}"/*.bw 2>/dev/null || echo "No .bw files found"
    exit 1
fi

echo ""
echo "Found all 4 required BigWig files"
echo ""

# Create metagene plot with TSS and TTS (scale-regions mode)
echo "Creating TSS to TTS metagene plot..."
computeMatrix scale-regions \
    -S "${all_bw_files[@]}" \
    -R ${gtf} \
    --beforeRegionStartLength 2000 \
    --regionBodyLength 5000 \
    --afterRegionStartLength 2000 \
    --skipZeros \
    -p ${threads} \
    --samplesLabel "${all_labels[@]}" \
    -o ${plot_output_dir}/${sample}_TSS_TTS.matrix.gz

plotProfile -m ${plot_output_dir}/${sample}_TSS_TTS.matrix.gz \
    --plotFileFormat pdf \
    --perGroup \
    --samplesLabel "${all_labels[@]}" \
    --plotHeight 7 \
    --plotWidth 10 \
    --startLabel "TSS" \
    --endLabel "TTS" \
    --plotTitle "${sample} - Gene metagene" \
    -out ${plot_output_dir}/${sample}_TSS_TTS_profile.pdf

echo "Created: ${plot_output_dir}/${sample}_TSS_TTS_profile.pdf"
echo "========================================="
echo "Completed ${sample}"
echo "Metagene plot saved to: ${plot_output_dir}"
echo "========================================="