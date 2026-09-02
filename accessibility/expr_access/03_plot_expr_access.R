#' Metaplots of SAM-seq accessibility (m6A) around protein-coding TSS by
#' expression bin, from the profile table written by 02_tss_m6a_profiles.py.
#'
#' Two views of the same table:
#'   * by_expression: one panel per LPS timepoint, lines colored by expr bin;
#'   * by_timepoint:  one panel per expr bin, lines colored by timepoint.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

TAB_DIR  <- "/project/spott/cshan/fiber-seq/macrophage_project/expr_access/tables"
PLOT_DIR <- "/project/spott/cshan/fiber-seq/macrophage_project/expr_access/plots"
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

prof <- fread(file.path(TAB_DIR, "tss_m6a_profile_by_expr_bin.tsv.gz"))

BIN_LEVELS <- c("not_expressed", "Q1_low", "Q2", "Q3", "Q4_high")
BIN_COLORS <- c("grey55", "#bdd7e7", "#6baed6", "#3182bd", "#08519c")
n_genes <- unique(prof[, .(expr_bin, n_genes)])
BIN_LABELS <- sapply(BIN_LEVELS, function(b) sprintf(
  "%s (n=%s)",
  c(not_expressed = "not expressed (TPM<1)", Q1_low = "Q1 (low)",
    Q2 = "Q2", Q3 = "Q3", Q4_high = "Q4 (high)")[b],
  format(n_genes[expr_bin == b, n_genes], big.mark = ",")))

prof[, expr_bin := factor(expr_bin, levels = BIN_LEVELS, labels = BIN_LABELS)]
prof[, timepoint := factor(sub("LPS_", "", timepoint),
                           levels = c("0", "5", "10", "15"),
                           labels = paste("LPS", c(0, 5, 10, 15), "min"))]
prof[, mid := (bin_start + bin_end) / 2]

base <- list(
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40",
             linewidth = 0.3),
  geom_line(linewidth = 0.5),
  labs(x = "Distance from TSS (bp, sense direction ->)",
       y = "m6A per covered A/T site"),
  theme_bw(base_size = 10),
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")
)

p1 <- ggplot(prof, aes(mid, m6a, color = expr_bin)) +
  base +
  facet_wrap(~timepoint, nrow = 1) +
  scale_color_manual(values = setNames(BIN_COLORS, BIN_LABELS),
                     name = "Expression (mean TPM, 0-15 min)") +
  guides(color = guide_legend(nrow = 2)) +
  ggtitle("accessibility around protein-coding TSS by expression level")

p2 <- ggplot(prof, aes(mid, m6a, color = timepoint)) +
  base +
  facet_wrap(~expr_bin, nrow = 1) +
  scale_color_brewer(palette = "RdYlBu", direction = -1, name = NULL) +
  ggtitle("accessibility around protein-coding TSS by LPS timepoint")

for (nm in c("tss_m6a_metaprofile_by_expression", "tss_m6a_metaprofile_by_timepoint")) {
  p <- if (nm == "tss_m6a_metaprofile_by_expression") p1 else p2
  ggsave(file.path(PLOT_DIR, paste0(nm, ".pdf")), p, width = 12, height = 4.2)
  ggsave(file.path(PLOT_DIR, paste0(nm, ".png")), p, width = 12, height = 4.2,
         dpi = 200)
}
message("Wrote plots to ", PLOT_DIR)
