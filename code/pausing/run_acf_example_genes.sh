#!/bin/bash
#SBATCH --job-name=acf_example_genes
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/acf_example_genes_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/acf_example_genes_%j.err
# Run acf_example_genes.R (sbatch or bash). rtracklayer/GenomicRanges need a newer libcurl
# ahead of the system one, and reticulate needs a Python with numpy/scipy/pandas.
set -euo pipefail

export LD_LIBRARY_PATH="/project/spott/cshan/envs/dimelo/lib:/software/gcc-12.2.0-el8-x86_64/lib64:/software/openblas-0.3.29-el8-x86_64/lib:/software/glpk-5.0-el8-x86_64/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export R_LIBS_USER="/project/spott/cshan/Rlibs/x86_64-pc-linux-gnu-library/4.4"
export RETICULATE_PYTHON="/project/spott/cshan/envs/Jupyter-notebook/bin/python"
# keep the sourced Python single-threaded, as the parent pipeline does
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

/software/R-4.4.1-el8-x86_64/bin/Rscript --vanilla /project/spott/cshan/fiber-seq/code/pausing/acf_example_genes.R
