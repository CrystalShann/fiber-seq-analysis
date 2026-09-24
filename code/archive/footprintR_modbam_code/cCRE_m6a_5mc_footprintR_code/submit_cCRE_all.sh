#!/bin/bash
# Submit all three cCRE metaprofile jobs in dependency order:
#   part1 (setup + nucleosomes) → part2 (m6A + 5mC BAM queries) → part3 (plots)
#
# Usage:
#   bash submit_cCRE_all.sh                          # defaults
#   SAMPLE=MY_SAMPLE CHROM=chr2 bash submit_cCRE_all.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Parameters (inherit from env or use defaults) ─────────────────────────────
export SAMPLE="${SAMPLE:-AL10_bc2178_19130}"
export CHROM="${CHROM:-chr1}"
export WINDOW_BP="${WINDOW_BP:-300}"
export BIN_SIZE="${BIN_SIZE:-1}"
export MOD_PROB_THRESHOLD="${MOD_PROB_THRESHOLD:-0.9}"
export CCRE_CHUNK_SIZE="${CCRE_CHUNK_SIZE:-500}"

echo "Submitting cCRE metaprofile jobs"
echo "  SAMPLE            = ${SAMPLE}"
echo "  CHROM             = ${CHROM}"
echo "  WINDOW_BP         = ${WINDOW_BP}"
echo "  BIN_SIZE          = ${BIN_SIZE}"
echo "  MOD_PROB_THRESHOLD= ${MOD_PROB_THRESHOLD}"
echo "  CCRE_CHUNK_SIZE   = ${CCRE_CHUNK_SIZE}"
echo ""

# ── Make sure log directory exists ───────────────────────────────────────────
mkdir -p /project/spott/cshan/fiber-seq/results/logs

# ── Part 1: cCRE setup + nucleosome alignment ─────────────────────────────────
JID1=$(sbatch --parsable "${SCRIPT_DIR}/submit_cCRE_part1.sh")
echo "Part 1 submitted: job ${JID1}"

# ── Part 2: m6A + 5mC BAM queries (depends on part 1) ────────────────────────
JID2=$(sbatch --parsable --dependency=afterok:${JID1} "${SCRIPT_DIR}/submit_cCRE_part2.sh")
echo "Part 2 submitted: job ${JID2} (after ${JID1})"

# ── Part 3: plot (depends on part 2) ─────────────────────────────────────────
JID3=$(sbatch --parsable --dependency=afterok:${JID2} "${SCRIPT_DIR}/submit_cCRE_part3.sh")
echo "Part 3 submitted: job ${JID3} (after ${JID2})"

echo ""
echo "Monitor with:"
echo "  squeue -u \$USER"
echo "  tail -f /project/spott/cshan/fiber-seq/results/logs/cCRE_part2_${JID2}.out"
