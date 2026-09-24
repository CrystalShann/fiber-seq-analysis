#!/bin/bash
#SBATCH --job-name=combine_fps
#SBATCH --time=24:00:00
#SBATCH --mem=120G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/combine_fps%A.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/combine_fps%A.err

module load python/miniforge-25.3.0
source activate /project/spott/cshan/envs/bedtools

python /project/spott/cshan/fiber-seq/code/PolII_footprints_code/combined_diff_size_footprints.py