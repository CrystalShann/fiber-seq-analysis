#!/bin/bash

#SBATCH --job-name=m6a_quartiles
#SBATCH --account=pi-spott
#SBATCH --partition=bigmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=300G
#SBATCH --time=24:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/m6a_quartiles_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/m6a_quartiles_%j.err

module load R

export SAMPLE=AL10_bc2178_19130
export CHROMS=AUTO
export WINDOW_BP=1000
export BIN_SIZE=10
export MOD_PROB_THRESHOLD=0.9
export MODBAM_TSS_CHUNK_SIZE=200
export TSS_ACCESS_WINDOW=100
export MAX_MODIFIED_HIST_BIN=10
export TSS_SOURCE_FILTER=ALL

export POLII_ROOT=/project/spott/cshan/fiber-seq/results/PolII
export PAUSING_PATH=${POLII_ROOT}/annotations/pausing_index_principal_with_CAGE_TSS_all_genes.tsv
export OUTPUT_ROOT=${POLII_ROOT}/m6a_pausing_quartiles
export BAM=/project/spott/1_Shared_projects/LCL_Fiber_seq/preprocess_final_merged_samples/${SAMPLE}.5mC.6mA.aligned.phased.bam

echo "Sample   : ${SAMPLE}"
echo "Window   : ${WINDOW_BP} bp"
echo "Bin      : ${BIN_SIZE} bp"
echo "AccessWin: ${TSS_ACCESS_WINDOW} bp"
echo "Output   : ${OUTPUT_ROOT}"
echo "Started  : $(date)"

Rscript /project/spott/cshan/fiber-seq/code/PolII_footprints_code/m6a_pausing_quartiles/run_m6a_pausing_quartiles.R

echo "Finished : $(date)"