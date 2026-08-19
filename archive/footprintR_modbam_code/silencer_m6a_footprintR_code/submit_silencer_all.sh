#!/bin/bash
set -euo pipefail

SAMPLE="${SAMPLE:-AL10_bc2178_19130}"
CHROM="${CHROM:-chr1}"
SILENCER_CLASS="${SILENCER_CLASS:-REST-Silencers}"
WINDOW_BP="${WINDOW_BP:-300}"
BIN_SIZE="${BIN_SIZE:-1}"
MOD_PROB_THRESHOLD="${MOD_PROB_THRESHOLD:-0.9}"
SILENCER_CHUNK_SIZE="${SILENCER_CHUNK_SIZE:-500}"
OUTPUT_SCOPE="${OUTPUT_SCOPE:-sample_chr}"

SCRIPT_DIR="/project/spott/cshan/fiber-seq/code/footprintR_modbam_code/silencer_m6a_footprintR_code"
PART1_SCRIPT="${SCRIPT_DIR}/submit_silencer_part1.sh"
PART2_SCRIPT="${SCRIPT_DIR}/submit_silencer_part2.sh"
PART3_SCRIPT="${SCRIPT_DIR}/submit_silencer_part3.sh"

if [[ ! -f "${PART1_SCRIPT}" || ! -f "${PART2_SCRIPT}" || ! -f "${PART3_SCRIPT}" ]]; then
  echo "Missing one or more part submit scripts in ${SCRIPT_DIR}" >&2
  exit 1
fi

submit_chain_for_chr() {
  local chr="$1"
  local scope="$2"

  local jid1
  local jid2
  local jid3

  jid1=$(sbatch \
    --export=ALL,SAMPLE="${SAMPLE}",CHROM="${chr}",WINDOW_BP="${WINDOW_BP}",BIN_SIZE="${BIN_SIZE}",SILENCER_CLASS="${SILENCER_CLASS}",OUTPUT_SCOPE="${scope}" \
    "${PART1_SCRIPT}" | awk '{print $4}')

  jid2=$(sbatch \
    --dependency=afterok:"${jid1}" \
    --export=ALL,SAMPLE="${SAMPLE}",CHROM="${chr}",WINDOW_BP="${WINDOW_BP}",BIN_SIZE="${BIN_SIZE}",MOD_PROB_THRESHOLD="${MOD_PROB_THRESHOLD}",SILENCER_CHUNK_SIZE="${SILENCER_CHUNK_SIZE}",SILENCER_CLASS="${SILENCER_CLASS}",OUTPUT_SCOPE="${scope}" \
    "${PART2_SCRIPT}" | awk '{print $4}')

  jid3=$(sbatch \
    --dependency=afterok:"${jid2}" \
    --export=ALL,SAMPLE="${SAMPLE}",CHROM="${chr}",WINDOW_BP="${WINDOW_BP}",SILENCER_CLASS="${SILENCER_CLASS}",OUTPUT_SCOPE="${scope}" \
    "${PART3_SCRIPT}" | awk '{print $4}')

  echo "class=${SILENCER_CLASS} sample=${SAMPLE} chr=${chr}: part1=${jid1} part2=${jid2} part3=${jid3}"
}

if [[ "${CHROM}" == "ALL" ]]; then
  ALL_CHRS="${ALL_CHRS:-chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY}"
  IFS=',' read -r -a CHR_LIST <<< "${ALL_CHRS}"
  for chr in "${CHR_LIST[@]}"; do
    submit_chain_for_chr "${chr}" "all_chr"
  done
else
  submit_chain_for_chr "${CHROM}" "${OUTPUT_SCOPE}"
fi
