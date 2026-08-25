# Global codependency distance-decay plots across the LPS timecourse.
#   - per-pair scores re-binned at plot time by floor(log2(dist)), <512 bp merged
#     into the first bin, >= 16 kb pooled (sparse fiber-length tail)
#   - Figure 1: mean binned score, one line per timepoint - the timepoints are
#     shown side by side but not tested against a baseline
#   - per-pair significance comes from fire_codependency.py's Fisher exact test
#     (fisher_p column); the fraction of pairs with p < 0.05 per timepoint and
#     distance bin is exported as a CSV
#   - every figure is saved together with the table of numbers behind it
#
# Usage:  module load R/4.4.1 && Rscript plot_codependency.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
})

codep_dir <- "/project/spott/cshan/fiber-seq/macrophage_project/FIRE_MSP/codependency"
fig_dir   <- file.path(codep_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE)

timepoints <- c("LPS_0", "LPS_5", "LPS_10", "LPS_15")
tp_colors  <- c(LPS_0 = "#2166AC", LPS_5 = "#92C5DE", LPS_10 = "#F4A582",
                LPS_15 = "#B2182B")

# save a figure together with the table behind it

save_with_table <- function(file, plot, data, ...) {
  ggsave(file.path(fig_dir, file), plot, bg = "transparent", ...)
  write_tsv(data, file.path(fig_dir, paste0(tools::file_path_sans_ext(file), ".tbl.gz")))
  message("wrote ", file)
}

## Per-pair scores, all timepoints, log2 distance bins -------------------------
all_pairs <- bind_rows(lapply(timepoints, function(s) {
  read_csv(file.path(codep_dir, paste0("codep_", s, ".csv.gz")), show_col_types = FALSE) |>
    mutate(timepoint = s)
})) |>
  mutate(bin_log2 = pmin(pmax(floor(log2(dist)), 8), 14),   # < 512 bp pooled; > 16 kb pooled
         timepoint = factor(timepoint, levels = timepoints))

bin_breaks <- 8:14
bin_labels <- c("<512bp", ">512bp", ">1kb", ">2kb", ">4kb", ">8kb", "≥16kb")

message("pairs per timepoint:")
print(count(all_pairs, timepoint))

binned <- all_pairs |>
  group_by(timepoint, bin_log2) |>
  summarize(mean_score = mean(score), sem = sd(score) / sqrt(n()), n = n(), .groups = "drop")

## Figure 1: global distance decay, one line per timepoint ---------------------
fig1 <- ggplot(binned, aes(bin_log2, mean_score, color = timepoint)) +
  geom_ribbon(aes(ymin = mean_score - sem, ymax = mean_score + sem, fill = timepoint),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
  scale_color_manual(values = tp_colors) + scale_fill_manual(values = tp_colors, guide = "none") +
  scale_x_continuous(breaks = bin_breaks, labels = bin_labels) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  xlab("Distance bins (bp between peaks)") +
  ylab("Mean binned codependency score") +
  theme_bw(11) + theme(legend.position = "top", legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
save_with_table("codep_distance_decay_global.pdf", fig1, binned, height = 4.5, width = 6)

## Per-pair Fisher significance, summarized per timepoint and distance bin -----
if (!"fisher_p" %in% colnames(all_pairs)) {
  message("no fisher_p column in codep_<s>.csv.gz - re-run fire_codependency.py first")
} else {
  sig <- all_pairs |>
    group_by(timepoint, bin_log2) |>
    summarize(n_pairs = n(),
              n_sig = sum(fisher_p < 0.05),
              frac_sig = mean(fisher_p < 0.05), .groups = "drop")
  write_csv(sig, file.path(fig_dir, "codep_fisher_sig_by_bin.csv"))
  message("wrote codep_fisher_sig_by_bin.csv")
  message("pairs with fisher_p < 0.05 per timepoint:")
  print(all_pairs |> group_by(timepoint) |>
          summarize(n_pairs = n(), n_sig = sum(fisher_p < 0.05),
                    frac_sig = mean(fisher_p < 0.05)))
}
