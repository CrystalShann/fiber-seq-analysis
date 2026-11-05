#!/bin/bash
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --array=1-31
#SBATCH --output=compartments-%x-%A_%a.log
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
exec bash "$SCRIPT_DIR/run_R.sh" "$SCRIPT_DIR/02_nucleosome_fraction.R" "$@"
