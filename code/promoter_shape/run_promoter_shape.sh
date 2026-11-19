#!/bin/bash
# Run the promoter-shape tables (01_promoter_shape.R); figures come from 02_promoter_shape_figures.Rmd
# (run its chunks in RStudio, or render it). A login node is fine
# (~1 min, < 4 GB); avoid the amd partition (project R crashes there with SIGILL).
set -euo pipefail

export LD_LIBRARY_PATH="/software/gcc-12.2.0-el8-x86_64/lib64:/software/openblas-0.3.29-el8-x86_64/lib:/software/glpk-5.0-el8-x86_64/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export R_LIBS_USER="/project/spott/cshan/Rlibs/x86_64-pc-linux-gnu-library/4.4"
RSCRIPT=/software/R-4.4.1-el8-x86_64/bin/Rscript
CODE=/project/spott/cshan/fiber-seq/code/promoter_shape

if [[ $# -gt 0 ]]; then
  "$RSCRIPT" --vanilla "$@"      # e.g. bash run_promoter_shape.sh code/promoter_shape/test_promoter_shape.R
else
  "$RSCRIPT" --vanilla "$CODE/01_promoter_shape.R"
fi
