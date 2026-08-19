# ---------------------------------------------------------------------------
# m6A promoter accessibility tables
# ---------------------------------------------------------------------------
source("/project/spott/cshan/fiber-seq/code/footprintR_modbam_code/run_modbam_promoter_summary_modality_common.R")

run_modality_chunked(
  modbase_code = "a",
  summary_prefix = "accessibility",
  mod_label = "m6A",
  meta_label = "Accessibility (m6A)"
)

message("m6A outputs written to: ", output_dir)
