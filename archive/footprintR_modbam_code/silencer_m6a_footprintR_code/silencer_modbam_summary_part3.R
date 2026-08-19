# ---------------------------------------------------------------------------
# Part 3: Plot silencer m6A metaprofiles (optionally with nucleosome structure)
#
# Parameters (override via env vars):
#   SAMPLE           [AL10_bc2178_19130]
#   CHROM            [chr1]
#   WINDOW_BP        [300]
#   SILENCER_CLASS   one of REST-Enhancers, REST-Silencers,
#                    STARR-Silencers.Robust, STARR-Silencers.Stringent,
#                    or ALL [ALL]
# ---------------------------------------------------------------------------

source("/project/spott/cshan/fiber-seq/code/footprintR_modbam_code/modbam_footprintR_functions.R")

default_sample <- Sys.getenv("SAMPLE", unset = "AL10_bc2178_19130")
default_chrom  <- Sys.getenv("CHROM", unset = "chr1")
window_bp      <- as.integer(Sys.getenv("WINDOW_BP", unset = "300"))
silencer_class <- Sys.getenv("SILENCER_CLASS", unset = "ALL")
output_scope   <- Sys.getenv("OUTPUT_SCOPE", unset = "sample_chr")

if (silencer_class == "ALL") {
  class_tag <- "all4classes"
} else {
  class_tag <- gsub("[^A-Za-z0-9._-]", "_", silencer_class)
}

results_root_dir <- file.path(
  "/project/spott/cshan/fiber-seq/results/cCRE_summary/silencer_summary/m6a_modbam",
  output_scope
)
output_dir       <- file.path(results_root_dir, "threshold_0.9", default_sample, class_tag)

nuc_meta_path <- file.path(output_dir, sprintf("%s_%s_%s_nucleosome_meta.tsv.gz", default_sample, default_chrom, class_tag))
m6a_meta_path <- file.path(output_dir, sprintf("%s_%s_%s_m6A_meta.tsv.gz", default_sample, default_chrom, class_tag))

assert_file_exists(m6a_meta_path, "m6A metaprofile (run part2 first)")

class_order <- c("REST-Enhancers", "REST-Silencers", "STARR-Silencers.Robust", "STARR-Silencers.Stringent")

m6a_meta <- data.table::fread(m6a_meta_path)
m6a_meta[, silencer_class := factor(silencer_class, levels = class_order)]

has_nuc <- file.exists(nuc_meta_path)
if (has_nuc) {
  nuc_meta <- data.table::fread(nuc_meta_path)
  nuc_meta[, silencer_class := factor(silencer_class, levels = class_order)]
}

x_breaks <- seq(-window_bp, window_bp, by = 50)
x_labels <- function(x) ifelse(x == 0, "0", sprintf("%+d", x))

base_theme <- ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(
    strip.text       = ggplot2::element_text(size = 11, face = "bold"),
    axis.text.x      = ggplot2::element_text(size = 9),
    panel.grid.minor = ggplot2::element_blank()
  )

vline <- ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.4)

plot_m6a <- function(meta_dt, sample_name, chrom, loess_span = 0.15) {
  ggplot2::ggplot(meta_dt, ggplot2::aes(x = bin, y = fraction_modified)) +
    ggplot2::geom_line(color = "#7B3294", linewidth = 0.4, na.rm = TRUE) +
    ggplot2::geom_smooth(method = "loess", span = loess_span, se = FALSE,
                         color = "#1f77b4", linewidth = 0.8, na.rm = TRUE) +
    ggplot2::scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    vline +
    ggplot2::facet_wrap(~ silencer_class, nrow = 1, scales = "free_y") +
    ggplot2::labs(
      title = sprintf("Silencer accessibility (m6A): %s %s", sample_name, chrom),
      x = "Distance from silencer center (bp)",
      y = "m6A modified fraction"
    ) +
    base_theme
}

plot_structure <- function(meta_dt, sample_name, chrom, loess_span = 0.15) {
  ggplot2::ggplot(meta_dt, ggplot2::aes(x = bin, y = nucleosome_midpoints_per_silencer)) +
    ggplot2::geom_line(color = "#4D4D4D", linewidth = 0.4) +
    ggplot2::geom_smooth(method = "loess", span = loess_span, se = FALSE,
                         color = "#1f77b4", linewidth = 0.8, na.rm = TRUE) +
    ggplot2::scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    vline +
    ggplot2::facet_wrap(~ silencer_class, nrow = 1, scales = "free_y") +
    ggplot2::labs(
      title = sprintf("Silencer structure (nucleosome midpoints): %s %s", sample_name, chrom),
      x = "Distance from silencer center (bp)",
      y = "Nucleosome midpoints per silencer"
    ) +
    base_theme
}

m6a_plot <- plot_m6a(m6a_meta, default_sample, default_chrom)

if (has_nuc) {
  structure_plot <- plot_structure(nuc_meta, default_sample, default_chrom)
  final_plot <- structure_plot / m6a_plot + patchwork::plot_layout(heights = c(1, 1.2))
  plot_height <- 14
} else {
  final_plot <- m6a_plot
  plot_height <- 8
}

out_pdf <- file.path(
  output_dir,
  sprintf("%s_%s_%s_silencer_m6A_metaprofiles_300bp_threshold_0.9.pdf", default_sample, default_chrom, class_tag)
)

ggplot2::ggsave(out_pdf, plot = final_plot, width = 18, height = plot_height)
message("PDF saved to: ", out_pdf)