#!/bin/bash
#SBATCH --job-name=LCL_coaccess
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/project/spott/cshan/fiber-seq/results/logs/LCL_coaccess_%j.out
#SBATCH --error=/project/spott/cshan/fiber-seq/results/logs/LCL_coaccess_%j.err

# Submit after the complete 31-task read-span array succeeds. No --chrom means
# all chr1-22/X/Y and one genome-wide BH family of unique pooled pairs.
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
# SLURM copies the batch file to its spool; resolve the installed sibling then.
if [[ ! -f "$SCRIPT_DIR/LCL_coaccess_cres.py" ]]; then
    SCRIPT_DIR=/project/spott/cshan/fiber-seq/code/co-accessibility/LCL
fi
PYTHON=/project/spott/cshan/envs/Jupyter-notebook/bin/python3
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
exec "$PYTHON" "$SCRIPT_DIR/LCL_coaccess_cres.py" "$@"
