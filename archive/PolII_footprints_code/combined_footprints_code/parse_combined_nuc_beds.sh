#!/bin/bash
#SBATCH --job-name=parse_combined_nuc_beds
#SBATCH --time=6:00:00
#SBATCH --mem=120G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/parse_combined_nuc_beds%A.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/parse_combined_nuc_beds%A.err

module load python
python3 /project/spott/cshan/fiber-seq/code/PolII_footprints_code/parse_combined_nuc_bed.py
