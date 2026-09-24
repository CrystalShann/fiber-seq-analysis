# ---------------------------------------------------------------------------
# Part 3: Plot multi-omic metaprofiles around cCRE centers
#
# Reads the metaprofile TSVs written by Part 1 (nucleosome_meta) and
# Part 2 (m6A_meta, 5mC_meta), then generates one PDF with three faceted
# panels — structure, accessibility, methylation — each faceted by
# cCRE class.
#
# Analogous to plot_TSS_m6a_5mc.r but for all cCRE classes at once.
# ---------------------------------------------------------------------------

source("/project/spott/cshan/fiber-seq/code/footprintR_modbam_code/modbam_footprintR_functions.R")

default_sample <- Sys.getenv("SAMPLE",    unset = "AL10_bc2178_19130")
default_chrom  <- Sys.getenv("CHROM",     unset = "chr1")
window_bp      <- as.integer(Sys.getenv("WINDOW_BP", unset = "300"))

# ── Paths ─────────────────────────────────────────────────────────────────────
results_root_dir <- "/project/spott/cshan/fiber-seq/results/cCRE_summary/300bp_center"
output_dir       <- file.path(results_root_dir, "threshold_0.9", default_sample)

nuc_meta_path <- file.path(output_dir, sprintf("%s_%s_nucleosome_meta.tsv.gz", default_sample, default_chrom))
m6a_meta_path <- file.path(output_dir, sprintf("%s_%s_m6A_meta.tsv.gz",        default_sample, default_chrom))
m5c_meta_path <- file.path(output_dir, sprintf("%s_%s_5mC_meta.tsv.gz",        default_sample, default_chrom))

assert_file_exists(m6a_meta_path, "m6A metaprofile (run part2 first)")
assert_file_exists(m5c_meta_path, "5mC metaprofile (run part2 first)")

# ── Class order (base classes + synthetic Enhancers) ──────────────────────────
ccre_class_order <- c("PLS","pELS","dELS","CA-TF","CA-CTCF","CA-H3K4me3","CA","TF","Enhancers")

# ── Load metaprofiles ─────────────────────────────────────────────────────────
m6a_meta <- data.table::fread(m6a_meta_path)
m5c_meta <- data.table::fread(m5c_meta_path)

m6a_meta[, ccre_class := factor(ccre_class, levels = ccre_class_order)]
m5c_meta[, ccre_class := factor(ccre_class, levels = ccre_class_order)]

has_nuc <- file.exists(nuc_meta_path)
if (has_nuc) {
  nuc_meta <- data.table::fread(nuc_meta_path)
  nuc_meta[, ccre_class := factor(ccre_class, levels = ccre_class_order)]
}

# ── cCRE count labels for facet strips ────────────────────────────────────────
# Infer n per class from m6A meta (total_calls column present after Part 2)
n_per_class <- m6a_meta[bin == 0L | bin == min(bin),
                        .(n = data.table::uniqueN(ccre_class)), by = ccre_class]
# Simpler: use the label column if available, otherwise just use ccre_class
facet_labeller <- ggplot2::as_labeller(setNames(
  as.character(ccre_class_order),
  ccre_class_order
))

# ── Plot helpers ───────────────────────────────────────────────────────────────
x_breaks <- seq(-window_bp, window_bp, by = 50)
x_labels <- function(x) ifelse(x == 0, "0", sprintf("%+d", x))

base_theme <- ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(
    strip.text       = ggplot2::element_text(size = 11, face = "bold"),
    axis.text.x      = ggplot2::element_text(size = 9),
    panel.grid.minor = ggplot2::element_blank()
  )

vline <- ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                              color = "red", linewidth = 0.4)

# Structure: nucleosome midpoints per cCRE (subset of classes)
plot_structure_ccre <- function(meta_dt, sample_name, chrom,
                                window_bp = 300L, loess_span = 0.15) {
  ggplot2::ggplot(meta_dt,
      ggplot2::aes(x = bin, y = nucleosome_midpoints_per_ccre)) +
    ggplot2::geom_line(color = "#4D4D4D", linewidth = 0.4) +
    ggplot2::geom_smooth(method = "loess", span = loess_span, se = FALSE,
                         color = "#1f77b4", linewidth = 0.8, na.rm = TRUE) +
    ggplot2::scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    vline +
    ggplot2::facet_wrap(~ ccre_class, nrow = 1, scales = "free_y") +
    ggplot2::labs(
      title = sprintf("cCRE structure (nucleosome midpoints): %s %s", sample_name, chrom),
      x     = "Distance from cCRE center (bp)",
      y     = "Nucleosome midpoints per cCRE"
    ) +
    base_theme
}

# Combined m6A (top) / 5mC (bottom) for a subset of cCRE classes
mod_colors <- c("m6A" = "#7B3294", "5mC" = "#8C510A")

plot_combined_mod <- function(dt_subset, title_suffix, loess_span = 0.15) {
  ggplot2::ggplot(dt_subset,
      ggplot2::aes(x = bin, y = fraction_modified, color = mod_type)) +
    ggplot2::geom_line(linewidth = 0.4, na.rm = TRUE) +
    ggplot2::geom_smooth(method = "loess", span = loess_span, se = FALSE,
                         linewidth = 0.8, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = mod_colors) +
    ggplot2::scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    vline +
    ggplot2::facet_grid(mod_type ~ ccre_class, scales = "free_y") +
    ggplot2::labs(
      title = sprintf("cCRE modifications (m6A / 5mC): %s %s %s",
                      default_sample, default_chrom, title_suffix),
      x     = "Distance from cCRE center (bp)",
      y     = "Modified fraction"
    ) +
    base_theme +
    ggplot2::theme(legend.position = "none")
}

# ── Combine m6A and 5mC into one long table ───────────────────────────────────
m6a_meta[, mod_type := "m6A"]
m5c_meta[, mod_type := "5mC"]

combined_mod <- rbind(
  m6a_meta[, .(bin, ccre_class, fraction_modified, mod_type)],
  m5c_meta[, .(bin, ccre_class, fraction_modified, mod_type)]
)
combined_mod[, mod_type := factor(mod_type, levels = c("m6A", "5mC"))]

# ── Split cCRE classes across 2 PDFs ─────────────────────────────────────────
classes_pdf1 <- ccre_class_order[1:5]   # PLS, pELS, dELS, CA-TF, CA-CTCF
classes_pdf2 <- ccre_class_order[6:9]   # CA-H3K4me3, CA, TF, Enhancers

make_pdf_plot <- function(classes, suffix) {
  dt_sub <- combined_mod[ccre_class %in% classes]
  dt_sub[, ccre_class := factor(ccre_class, levels = classes)]
  mod_plot <- plot_combined_mod(dt_sub, title_suffix = suffix)

  if (has_nuc) {
    nuc_sub <- nuc_meta[ccre_class %in% classes]
    nuc_sub[, ccre_class := factor(ccre_class, levels = classes)]
    str_plot <- plot_structure_ccre(nuc_sub, default_sample, default_chrom,
                                    window_bp, loess_span = 0.15)
    str_plot / mod_plot + patchwork::plot_layout(heights = c(1, 2))
  } else {
    mod_plot
  }
}

plot1 <- make_pdf_plot(classes_pdf1, "(classes 1–5)")
plot2 <- make_pdf_plot(classes_pdf2, "(classes 6–9)")

# ── Save two PDFs ─────────────────────────────────────────────────────────────
plot_width  <- 28
plot_height <- if (has_nuc) 18 else 12

out_pdf1 <- file.path(output_dir,
  sprintf("%s_%s_cCRE_multiomic_metaprofiles_300bp_threshold_0.9_part1.pdf",
          default_sample, default_chrom))
out_pdf2 <- file.path(output_dir,
  sprintf("%s_%s_cCRE_multiomic_metaprofiles_300bp_threshold_0.9_part2.pdf",
          default_sample, default_chrom))

ggplot2::ggsave(out_pdf1, plot = plot1, width = plot_width, height = plot_height)
ggplot2::ggsave(out_pdf2, plot = plot2, width = plot_width, height = plot_height)

message("PDF 1 saved to: ", out_pdf1)
message("PDF 2 saved to: ", out_pdf2)
list(plot1 = plot1, plot2 = plot2)
