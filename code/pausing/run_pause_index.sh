#!/bin/bash
#SBATCH --job-name=pause_index
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/pause_index_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/pause_index_%j.err
# Run pause_index.R (sbatch or bash). rtracklayer will not load unless a newer libcurl
# precedes the system one on LD_LIBRARY_PATH ("undefined symbol: curl_easy_option_next").
set -euo pipefail

export LD_LIBRARY_PATH="/project/spott/cshan/envs/dimelo/lib:/software/gcc-12.2.0-el8-x86_64/lib64:/software/openblas-0.3.29-el8-x86_64/lib:/software/glpk-5.0-el8-x86_64/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export R_LIBS_USER="/project/spott/cshan/Rlibs/x86_64-pc-linux-gnu-library/4.4"

/software/R-4.4.1-el8-x86_64/bin/Rscript --vanilla /project/spott/cshan/fiber-seq/code/pausing/pause_index.R "$@"
