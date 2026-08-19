# ---------------------------------------------------------------------------
# Join promoter structure + m6A accessibility + 5mC methylation
# ---------------------------------------------------------------------------
source("/project/spott/cshan/fiber-seq/code/footprintR_modbam_code/modbam_footprintR_functions.R")

default_sample <- Sys.getenv("SAMPLE", unset = "AL10_bc2178_19130")
default_chrom <- Sys.getenv("CHROM", unset = "chr1")

results_root_dir <- "/project/spott/cshan/fiber-seq/results/promoter_summary/10kb_TSS"
threshold_dir <- file.path(results_root_dir, "threshold_0.9")
output_dir <- file.path(threshold_dir, default_sample)

tss_table_path <- file.path(output_dir, sprintf("%s_%s_tss_table.tsv.gz", default_sample, default_chrom))
nuc_summary_path <- file.path(output_dir, sprintf("%s_%s_promoter_structure.tsv.gz", default_sample, default_chrom))
m6a_summary_path <- file.path(output_dir, sprintf("%s_%s_accessibility.tsv.gz", default_sample, default_chrom))
m5c_summary_path <- file.path(output_dir, sprintf("%s_%s_methylation.tsv.gz", default_sample, default_chrom))
final_path <- file.path(output_dir, sprintf("%s_%s_promoter_structure_accessibility_methylation.tsv.gz", default_sample, default_chrom))

assert_file_exists(tss_table_path, "part1 TSS table")
assert_file_exists(nuc_summary_path, "part1 promoter structure table")
assert_file_exists(m6a_summary_path, "m6A promoter accessibility table")
assert_file_exists(m5c_summary_path, "5mC promoter methylation table")

tss_table <- data.table::fread(tss_table_path)
nuc_summary <- data.table::fread(nuc_summary_path)
m6a_summary <- data.table::fread(m6a_summary_path)
m5c_summary <- data.table::fread(m5c_summary_path)

promoter_multiomic_summary <- merge(tss_table, nuc_summary, by = "tss_uid", all.x = TRUE)
promoter_multiomic_summary <- merge(promoter_multiomic_summary, m6a_summary, by = "tss_uid", all.x = TRUE)
promoter_multiomic_summary <- merge(promoter_multiomic_summary, m5c_summary, by = "tss_uid", all.x = TRUE)

zero_fill_cols <- intersect(
  c(
    "structure_nucleosome_count", "structure_nucleosome_density_per_kb",
    "structure_core_nucleosome_count",
    "accessibility_modified_calls", "accessibility_unmodified_calls",
    "accessibility_total_calls", "accessibility_positions_with_coverage",
    "methylation_modified_calls", "methylation_unmodified_calls",
    "methylation_total_calls", "methylation_positions_with_coverage"
  ),
  names(promoter_multiomic_summary)
)

for (col in zero_fill_cols) {
  promoter_multiomic_summary[is.na(get(col)), (col) := 0]
}

na_fill_cols <- intersect(
  c(
    "structure_mean_nucleosome_size", "structure_median_nucleosome_size",
    "structure_mean_midpoint_offset", "structure_mean_abs_midpoint_offset",
    "accessibility_fraction_modified", "methylation_fraction_modified"
  ),
  names(promoter_multiomic_summary)
)

for (col in na_fill_cols) {
  promoter_multiomic_summary[is.infinite(get(col)) | is.nan(get(col)), (col) := NA_real_]
}

data.table::setcolorder(
  promoter_multiomic_summary,
  c(
    "chromosome", "tss_coordinate", "strand", "gene_name_or_tss_id",
    "tss_id", "tss_uid", "promoter_start", "promoter_end",
    "structure_nucleosome_count", "structure_mean_nucleosome_size",
    "structure_median_nucleosome_size", "structure_nucleosome_density_per_kb",
    "structure_core_nucleosome_count", "structure_mean_midpoint_offset",
    "structure_mean_abs_midpoint_offset",
    "accessibility_modified_calls", "accessibility_unmodified_calls",
    "accessibility_total_calls", "accessibility_positions_with_coverage",
    "accessibility_fraction_modified",
    "methylation_modified_calls", "methylation_unmodified_calls",
    "methylation_total_calls", "methylation_positions_with_coverage",
    "methylation_fraction_modified"
  )
)

write_gz_tsv(promoter_multiomic_summary, final_path)
message("Joined table written to: ", final_path)
