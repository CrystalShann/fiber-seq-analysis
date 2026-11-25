#!/bin/bash
# Consensus and non-consensus clustered m6A heatmaps of the parameters_test sweep (all k x resolution x
# seed) for every PolII example gene (see parameters_test_plots.R).
#
# sbatch /project/spott/cshan/fiber-seq/code/clustering_methods/Leiden_Manhattan/parameters_test/parameters_test_plots.sh
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --cpus-per-task=48
#SBATCH --mem=160G
#SBATCH --time=24:00:00
#SBATCH --job-name=leiden_parameters_both
#SBATCH --output=/project/spott/cshan/fiber-seq/LCL_project/Leiden_manhattan/parameters_test/logs/plots_both_%j.log
set -euo pipefail
export PARAM_WORKERS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH="/software/gcc-12.2.0-el8-x86_64/lib64:/software/openblas-0.3.29-el8-x86_64/lib:/software/glpk-5.0-el8-x86_64/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export R_LIBS_USER="/project/spott/cshan/Rlibs/x86_64-pc-linux-gnu-library/4.4"
/software/R-4.4.1-el8-x86_64/bin/Rscript --vanilla \
  /project/spott/cshan/fiber-seq/code/clustering_methods/Leiden_Manhattan/parameters_test/parameters_test_plots.R
