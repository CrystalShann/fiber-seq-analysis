#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export LD_LIBRARY_PATH="/software/gcc-12.2.0-el8-x86_64/lib64:/software/openblas-0.3.29-el8-x86_64/lib:/software/glpk-5.0-el8-x86_64/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export R_LIBS_USER="/project/spott/cshan/Rlibs/x86_64-pc-linux-gnu-library/4.4"
exec /software/R-4.4.1-el8-x86_64/bin/Rscript --vanilla "$@"
