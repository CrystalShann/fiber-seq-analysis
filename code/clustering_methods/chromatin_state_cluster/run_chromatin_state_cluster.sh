#!/usr/bin/env bash
# Default: one preparation task per chromosome (1-22, X, Y).
# This script never submits additional jobs.
#SBATCH --job-name=chromatin_states
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --array=1-24
#SBATCH --output=chromatin_states.%A_%a.log
set -euo pipefail

stage=prepare
if (( $# )) && [[ "$1" != --* ]]; then
  stage="$1"
  shift
fi
export R_PROFILE_USER=/dev/null
export LD_LIBRARY_PATH="/project/spott/cshan/envs/dimelo/lib:/software/gcc-12.2.0-el8-x86_64/lib64:/software/openblas-0.3.29-el8-x86_64/lib:/software/glpk-5.0-el8-x86_64/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export R_LIBS_USER="/project/spott/cshan/Rlibs/x86_64-pc-linux-gnu-library/4.4"
rscript=/software/R-4.4.1-el8-x86_64/bin/Rscript
driver=/project/spott/cshan/fiber-seq/code/clustering_methods/chromatin_state_cluster/chromatin_state_cluster.R

if [[ "$stage" == all ]]; then
  if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo 'Use all as a single job, without --array; chromosomes are pooled for clustering.' >&2
    exit 2
  fi
  chromosome_list=chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY
  genes_path=/project/spott/cshan/fiber-seq/LCL_project/pausing/dominant_tss/pause_index_dominant.tsv
  feature_dir=/project/spott/cshan/fiber-seq/LCL_project/pausing/dominant_tss
  shared_args=()
  while (( $# )); do
    if (( $# < 2 )); then
      echo "Expected --option value: $1" >&2
      exit 2
    fi
    case "$1" in
      --chrom|--state|--stage)
        echo "all mode controls $1; use --chroms to select the pooled chromosome set." >&2
        exit 2 ;;
      --chroms) chromosome_list="$2" ;;
      --genes) genes_path="$2"; shared_args+=("$1" "$2") ;;
      --feature-dir) feature_dir="$2"; shared_args+=("$1" "$2") ;;
      *) shared_args+=("$1" "$2") ;;
    esac
    shift 2
  done
  mkdir -p -- "$feature_dir"
  if [[ ! -r "$genes_path" ]]; then
    echo "Missing gene-anchor table: $genes_path. Generate the upstream table or supply --genes with the table used to generate the promoter features." >&2
    exit 2
  fi
  IFS=, read -r -a chromosomes <<< "$chromosome_list"
  if (( ${#chromosomes[@]} == 0 )) || [[ "$chromosome_list" == *, ]]; then
    echo 'Empty chromosome selection.' >&2
    exit 2
  fi
  # Check all feature inputs before starting any preparation.
  for chromosome in "${chromosomes[@]}"; do
    if [[ ! "$chromosome" =~ ^chr([1-9]|1[0-9]|2[0-2]|X|Y)$ ]]; then
      echo "Unsupported chromosome: $chromosome" >&2
      exit 2
    fi
    if [[ ! -r "$feature_dir/fiber_features.$chromosome.tsv.gz" ]]; then
      echo "Missing promoter features: $feature_dir/fiber_features.$chromosome.tsv.gz" >&2
      exit 2
    fi
  done
  for chromosome in "${chromosomes[@]}"; do
    "$rscript" "$driver" "${shared_args[@]}" --stage prepare --chrom "$chromosome"
  done
  "$rscript" "$driver" "${shared_args[@]}" --stage cluster --chroms "$chromosome_list" --state all
  "$rscript" "$driver" "${shared_args[@]}" --stage report --chroms "$chromosome_list"
  exit 0
fi

args=(--stage "$stage")
case "$stage" in
  prepare)
    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
      chromosomes=({1..22} X Y)
      task="$SLURM_ARRAY_TASK_ID"
      [[ "$task" =~ ^[0-9]+$ ]] && (( task >= 1 && task <= 24 )) || exit 2
      args+=(--chrom "chr${chromosomes[task-1]}")
    fi
    ;;
  cluster)
    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
      states=(closed open_fire open_linker other)
      task="$SLURM_ARRAY_TASK_ID"
      [[ "$task" =~ ^[0-9]+$ ]] && (( task >= 1 && task <= 4 )) || exit 2
      args+=(--state "${states[task-1]}")
    fi
    ;;
  report)
    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
      echo 'Run report as a single job, without an array.' >&2
      exit 2
    fi
    ;;
  *) echo "Unknown stage: $stage" >&2; exit 2 ;;
esac

exec "$rscript" "$driver" \
  "${args[@]}" "$@"
