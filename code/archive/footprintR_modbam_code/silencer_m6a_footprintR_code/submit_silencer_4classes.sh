#!/bin/bash
set -euo pipefail

SAMPLE="${SAMPLE:-AL10_bc2178_19130}"
CHROM="${CHROM:-chr1}"
WINDOW_BP="${WINDOW_BP:-300}"
BIN_SIZE="${BIN_SIZE:-1}"
MOD_PROB_THRESHOLD="${MOD_PROB_THRESHOLD:-0.9}"
SILENCER_CHUNK_SIZE="${SILENCER_CHUNK_SIZE:-500}"
OUTPUT_SCOPE="${OUTPUT_SCOPE:-sample_chr}"

SCRIPT_DIR="/project/spott/cshan/fiber-seq/code/footprintR_modbam_code/silencer_m6a_footprintR_code"
ALL_SCRIPT="${SCRIPT_DIR}/submit_silencer_all.sh"

if [[ ! -f "${ALL_SCRIPT}" ]]; then
  echo "Missing ${ALL_SCRIPT}" >&2
  exit 1
fi

CLASSES=(
  "REST-Enhancers"
  "REST-Silencers"
  "STARR-Silencers.Robust"
  "STARR-Silencers.Stringent"
)

for class in "${CLASSES[@]}"; do
  echo "Submitting class=${class} sample=${SAMPLE} chrom=${CHROM}"
  SAMPLE="${SAMPLE}" \
  CHROM="${CHROM}" \
  WINDOW_BP="${WINDOW_BP}" \
  BIN_SIZE="${BIN_SIZE}" \
  MOD_PROB_THRESHOLD="${MOD_PROB_THRESHOLD}" \
  SILENCER_CHUNK_SIZE="${SILENCER_CHUNK_SIZE}" \
  SILENCER_CLASS="${class}" \
  OUTPUT_SCOPE="${OUTPUT_SCOPE}" \
  bash "${ALL_SCRIPT}"
done
