# ---------------------------------------------------------------------------
# 5mC promoter methylation tables
# ---------------------------------------------------------------------------
source("/project/spott/cshan/fiber-seq/code/footprintR_modbam_code/run_modbam_promoter_summary_modality_common.R")

run_modality_chunked(
  modbase_code = "m",
  summary_prefix = "methylation",
  mod_label = "5mC",
  meta_label = "Methylation (5mC)"
)

message("5mC outputs written to: ", output_dir)
