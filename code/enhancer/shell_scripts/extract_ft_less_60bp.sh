#!/bin/bash
#SBATCH --job-name=extract_tf_footprints_lt60bp
#SBATCH --account=pi-spott
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=1-4
#SBATCH --output=/project/spott/cshan/fiber-seq/code/enhancer/logs/out_%A_%a.out
#SBATCH --error=/project/spott/cshan/fiber-seq/code/enhancer/logs/err_%A_%a.err

set -euo pipefail
shopt -s failglob
samples=(LPS0 LPS5 LPS10 LPS15)
[[ ${SLURM_ARRAY_TASK_ID:-} =~ ^[1-4]$ ]] || { echo "Submit with sbatch" >&2; exit 1; }
sample=${samples[SLURM_ARRAY_TASK_ID-1]}
root=${TF_EXTRACT_ROOT:-/project/spott/cshan/fiber-seq/macrophage_project/FiberHMM/extract}
bin=/project/spott/cshan/envs/dimelo/bin
outdir="$root/firehmm_tf/ft_by_size/size0-60/$sample"
mkdir -p "$outdir"
workdir=$(mktemp -d "$outdir/.tf_lt60.XXXXXX")
trap 'rm -rf -- "$workdir"' EXIT

for input in "$root/firehmm_tf/$sample/${sample}_hmm_extracted_tf_"*.bed.gz; do
    name=${input##*/}
    output="$outdir/${name/_hmm_extracted_tf_/_tf_size0-60_}"
    [[ -s "$output" && -s "$output.tbi" ]] && continue

    gzip -cd "$input" | awk '
      BEGIN { FS = OFS = "\t" }
      {
        if (NF < 12) exit 1
        sub(/,$/, "", $11); sub(/,$/, "", $12)
        n = split($11, sizes, ","); m = split($12, offsets, ",")
        if ($10 > 0 && (n != $10 || m != $10)) exit 1
        for (i = 1; i <= $10; i++)
          if (sizes[i] > 0 && sizes[i] < 60) {
            start = $2 + offsets[i]
            print $1, start, start + sizes[i], $4, sizes[i], $6
          }
      }' | LC_ALL=C sort -k1,1 -k2,2n -k3,3n -S 2G -T "$workdir" \
         | "$bin/bgzip" -@ "${SLURM_CPUS_PER_TASK:-4}" -c > "$workdir/calls.bed.gz"
    "$bin/tabix" -p bed "$workdir/calls.bed.gz"
    mv "$workdir/calls.bed.gz" "$output"
    mv "$workdir/calls.bed.gz.tbi" "$output.tbi"
    echo "Wrote: $output"
done
