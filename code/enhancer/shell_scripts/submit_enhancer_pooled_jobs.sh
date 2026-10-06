#!/bin/bash
# Prepare once; both methods then read exactly the same saved capped cohort.
set -euo pipefail
PROJECT=${ENHANCER_PROJECT_ROOT:-/project/spott/cshan/fiber-seq}
table_dir="$PROJECT/macrophage_project/enhancer/TF_co-occ/tables"
shared_file=${ENHANCER_SHARED_FIBERS:-$table_dir/enhancer_shared_fibers_pooled_capped.rds}
record="$table_dir/enhancer_pooled_job_ids.tsv"
mkdir -p "$table_dir"
if [[ -f "$record" ]]; then
    previous_ids=$(awk -F '\t' 'NR > 1 { if (n++) printf ","; printf "%s", $2 }' "$record")
    if [[ -n "$previous_ids" ]]; then
        active=$(squeue -h -j "$previous_ids" -o '%i')
        if [[ -n "$active" ]]; then
            echo "Previous pooled jobs are still active: $active. Do not submit duplicates." >&2
            exit 1
        fi
    fi
fi
submit_job() {
    local stage=$1
    shift
    local raw job_id
    raw=$(sbatch --parsable "$@")
    job_id=${raw%%;*}
    [[ $job_id =~ ^[0-9]+$ ]] || { echo "Unexpected sbatch response: $raw" >&2; exit 1; }
    printf '%s\t%s\n' "$stage" "$job_id" >> "$record"
    printf '%s\n' "$job_id"
}
printf 'stage\tjob_id\n' > "$record"
export ENHANCER_SHARED_FIBERS="$shared_file"
prepare=$(submit_job prepare --export=ALL "$PROJECT/code/enhancer/shell_scripts/run_enhancer_shared_sampling.sh")
manhattan=$(submit_job manhattan --dependency="afterok:$prepare" --export=ALL "$PROJECT/code/enhancer/shell_scripts/run_enhancer_manhattan.sh")
acf=$(submit_job acf --dependency="afterok:$prepare" --export=ALL "$PROJECT/code/enhancer/shell_scripts/run_enhancer_acf.sh")
printf 'Shared sample: %s\nManhattan (pooled): %s\nACF (pooled): %s\n' "$prepare" "$manhattan" "$acf"
printf 'Job record: %s\nShared input: %s\n' "$record" "$shared_file"
