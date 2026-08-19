# this script generates multi-omic metaprofiles of promoter structure, accessibility, and methylation for a single sample and chromosome, using the modbam_footprintR_functions.R for plotting.

source("/project/spott/cshan/fiber-seq/code/footprintR_modbam_code/modbam_footprintR_functions.R")

output_dir <- "/project/spott/cshan/fiber-seq/results/promoter_summary/10kb_TSS/threshold_0.9/AL10_bc2178_19130"
default_sample <- "AL10_bc2178_19130"
default_chrom <- "chr1"

nuc_meta <- data.table::fread(
  file.path(output_dir, sprintf("%s_%s_nucleosome_meta.tsv.gz", default_sample, default_chrom))
)
m6a_meta <- data.table::fread(
  file.path(output_dir, sprintf("%s_%s_m6A_meta.tsv.gz", default_sample, default_chrom))
)


m5c_meta <- data.table::fread(
  file.path(output_dir, sprintf("%s_%s_5mC_meta.tsv.gz", default_sample, default_chrom))
)


structure_plot <- plot_structure_meta(
  meta_dt         = nuc_meta,
  sample_name     = default_sample,
  chrom           = default_chrom,
  window_bp        = 10000L,       # x-axis range (+/- bp)
  x_break_interval = 1000,         # tick spacing (bp)
  loess_span       = 0.06,        # smoothing bandwidth
  ylim             = NULL,        # e.g. c(0, 0.5) to fix y-axis
  line_color       = "#4D4D4D",
  smooth_color     = "#1f77b4",
  linewidth        = 0.6,
  smooth_linewidth = 0.9,
  base_size        = 12
)


accessibility_plot <- plot_modification_meta(
  meta_dt          = m6a_meta,
  title            = sprintf("Promoter accessibility: %s %s", default_sample, default_chrom),
  color            = "#7B3294",
  window_bp        = 10000L,       # x-axis range (+/- bp)
  x_break_interval = 1000,         # tick spacing (bp)
  loess_span       = 0.06,        # smoothing bandwidth
  ylim             = NULL,        # e.g. c(0, 0.3) to fix y-axis
  linewidth        = 0.6,
  smooth_linewidth = 0.9,
  y_label          = "Modified fraction",
  base_size        = 12
)

methylation_plot <- plot_modification_meta(
  meta_dt          = m5c_meta,
  title            = sprintf("Promoter methylation: %s %s", default_sample, default_chrom),
  color            = "#8C510A",
  window_bp        = 10000L,       # x-axis range (+/- bp)
  x_break_interval = 1000,         # tick spacing (bp)
  loess_span       = 0.06,        # smoothing bandwidth
  ylim             = NULL,        # e.g. c(0, 0.3) to fix y-axis
  linewidth        = 0.6,
  smooth_linewidth = 0.9,
  y_label          = "Modified fraction",
  base_size        = 12
)




# Overlay the plots using patchwork
# structure_plot is included as-is, which allows it to autoscale
combined_plot <- structure_plot / accessibility_plot / methylation_plot + 
  patchwork::plot_layout(heights = c(1, 1, 1))

# Save the plot
ggplot2::ggsave(
  filename = file.path(
    output_dir,
    sprintf("%s_%s_promoter_multiomic_metaprofiles_10kb_threshold_0.9.pdf", default_sample, default_chrom)
  ),
  plot = combined_plot,
  width = 10,
  height = 16
)

# Display the plot
combined_plot