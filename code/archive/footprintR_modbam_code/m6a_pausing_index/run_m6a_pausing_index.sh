#!/bin/bash
#SBATCH --job-name=pausing_m6a_allchr
#SBATCH --time=24:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=1
#SBATCH --account=pi-spott
#SBATCH --partition=bigmem
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/pausing_m6a_allchr_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/pausing_m6a_allchr_%j.err



module load R

export SAMPLE=AL10_bc2178_19130
export CHROMS=AUTO
export WINDOW_BP=1000
export BIN_SIZE=10
export MOD_PROB_THRESHOLD=0.9
export MODBAM_TSS_CHUNK_SIZE=200
export LOW_Q=0.25
export HIGH_Q=0.75
export TSS_SOURCE_FILTER=ALL
export POLII_ROOT=/project/spott/cshan/fiber-seq/results/PolII
export PAUSING_PATH=/project/spott/cshan/fiber-seq/results/PolII/annotations/pausing_index_principal_with_CAGE_TSS_all_genes.tsv
export OUTPUT_ROOT=/project/spott/cshan/fiber-seq/results/PolII/m6a_pausing_index
export BAM=/project/spott/1_Shared_projects/LCL_Fiber_seq/preprocess_final_merged_samples/AL10_bc2178_19130.5mC.6mA.aligned.phased.bam

export OUTPUT_DIR=${OUTPUT_ROOT}/${SAMPLE}_1kb_bin10_modthresh0.9_all

Rscript /project/spott/cshan/fiber-seq/code/footprintR_modbam_code/m6a_pausing_index/run_m6a_pausing_index.R

Rscript /project/spott/cshan/fiber-seq/code/footprintR_modbam_code/m6a_pausing_index/plot_m6a_pausing_index_outputs.R