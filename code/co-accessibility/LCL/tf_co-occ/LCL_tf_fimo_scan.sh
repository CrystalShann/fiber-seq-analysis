#!/bin/bash
#SBATCH --job-name=LCL_tf_fimo
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/LCL_tf_fimo_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/LCL_tf_fimo_%j.err
# LCL adaptation of code/fire_frequency/motif_fimo_scan.sh; original is untouched.
set -euo pipefail
SCRIPT_DIR=${LCL_TF_CODE:-/project/spott/cshan/fiber-seq/code/co-accessibility/LCL/tf_co-occ}
PYTHON=${LCL_TF_PYTHON:-/project/spott/cshan/envs/Jupyter-notebook/bin/python3}
ROOT=${LCL_COACCESS_ROOT:-/project/spott/cshan/fiber-seq/LCL_project/co-accessibility}
OUT=${LCL_TF_OUT:-$ROOT/tf_co-occ}
FIMO=${LCL_FIMO:-/project/spott/cshan/tools/meme-5.5.9/src/fimo}
MOTIF=${LCL_TF_MOTIFS:-/project/spott/cshan/annotations/JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.txt}
REFERENCE=${LCL_TF_REFERENCE:-/project/spott/reference/human/GRCh38.p14/hg38.fa}
exec "$PYTHON" "$SCRIPT_DIR/LCL_tf_fimo.py" --out "$OUT" --fimo "$FIMO" --motif "$MOTIF" --reference "$REFERENCE"
