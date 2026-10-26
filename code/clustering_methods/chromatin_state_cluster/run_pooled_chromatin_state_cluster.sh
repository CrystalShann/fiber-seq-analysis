#!/usr/bin/env bash
# Run after all chromosome preparation tasks succeed. No chromosome array:
# every broad state is fitted using all prepared chromosomes together.
#SBATCH --job-name=chromatin_pooled
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=chromatin_pooled.%j.log
set -euo pipefail
runner=/project/spott/cshan/fiber-seq/code/clustering_methods/chromatin_state_cluster/run_chromatin_state_cluster.sh
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo 'Submit pooled clustering without a job array.' >&2
  exit 2
fi
for argument in "$@"; do
  case "$argument" in
    --chrom|--state|--stage)
      echo "Pooled runner controls $argument; use --chroms to select chromosomes." >&2
      exit 2 ;;
  esac
done
bash "$runner" cluster "$@"
bash "$runner" report "$@"
