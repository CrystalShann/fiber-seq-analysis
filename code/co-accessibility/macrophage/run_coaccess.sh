#!/bin/bash
#SBATCH --job-name=coaccess
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/coaccess_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/coaccess_%j.err

# cCRE co-accessibility across the LPS timecourse.
# Requires 01_make_fire_universe.sh and 02_read_spans.sh to have finished.
#
# Usage:
#   sbatch run_coaccess.sh
#   sbatch run_coaccess.sh --max-dist 10000        # Kevin's single-sample range
#   sbatch run_coaccess.sh --read-rule contain     # require fibers to span the cCRE
#   sbatch run_coaccess.sh --fire-overlap-fraction 0.5  # fraction of cCRE covered by one FIRE
#   sbatch run_coaccess.sh --kevin-compat          # drop zero-accessibility pairs

set -euo pipefail

PYTHON=/project/spott/cshan/envs/Jupyter-notebook/bin/python3
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# SLURM runs a copied shell script from its spool directory, so use the source
# directory when the adjacent engine is unavailable. This can be overridden.
if [ ! -f "${SCRIPT_DIR}/03_coaccess_cres.py" ]; then
    SCRIPT_DIR="${COACCESS_MACROPHAGE_SCRIPT_DIR:-/project/spott/cshan/fiber-seq/code/co-accessibility/macrophage}"
fi

"$PYTHON" "${SCRIPT_DIR}/03_coaccess_cres.py" "$@"
