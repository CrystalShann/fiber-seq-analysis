# Codependency distance-decay plots across the LPS timecourse, following the
# logic of the Stergachis-lab scDAF-seq codependency figure:
#   - per-pair scores re-binned at plot time by floor(log2(dist)), <512 bp merged
#     into the first bin
#   - line plots of mean binned score, one line per condition; here the
#     conditions are the four timepoints (plus the shuffled-independence null),
#     shown separately for the global pair set and for LPS-response loci
#   - difference curves vs baseline (timepoint - LPS_0), the analog of their
#     Same - Opposite / Same - DiffCells panels
#   - per-bin significance exported as a CSV table, not drawn on the figure;
#     here PAIRED t-tests (same pairs tracked at all timepoints) instead of
#     their unpaired Same-vs-Opposite tests
#   - every figure is saved together with the table of numbers behind it
#     (their my_ggsave convention)
#
# LPS-response loci (user-defined):
#   tier 1, immediate-early (~15 min):  EGR1 EGR2 EGR3 FOS FOSB JUN NFKBIA
#       NFKBIZ TNFAIP3 DUSP1 DUSP2 IER3 ATF3 TNF PTGS2 MAFF
#   tier 2, secondary (~1-3 h):         IL1B IL1A CCL3 CCL4 CXCL8 IFNB1 CXCL10
#       SOCS1 SOCS3 IFIT2 IFIT3 IFIT5
# A pair belongs to a tier if EITHER of its peaks lies within 250 bp of a
# canonical TSS of a tier gene (from the FIRE_peaks_canonicalTSS_250bp file).
#
# Usage:  module load R/4.4.1 && Rscript plot_codependency.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
})

codep_dir <- "/project/spott/cshan/fiber-seq/macrophage_project/FIRE_MSP/codependency"
tss_file  <- "/project/spott/cshan/fiber-seq/macrophage_project/FIRE_MSP/LPS_0/LPS_0.FIRE_peaks_canonicalTSS_250bp.bed.gz"
fig_dir   <- file.path(codep_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE)

timepoints <- c("LPS_0", "LPS_5", "LPS_10", "LPS_15")
tp_colors  <- c(LPS_0 = "#2166AC", LPS_5 = "#92C5DE", LPS_10 = "#F4A582",
                LPS_15 = "#B2182B", Shuffled = "grey60")

tier1 <- c("EGR1","EGR2","EGR3","FOS","FOSB","JUN","NFKBIA","NFKBIZ","TNFAIP3",
           "DUSP1","DUSP2","IER3","ATF3","TNF","PTGS2","MAFF")
tier2 <- c("IL1B","IL1A","CCL3","CCL4","CXCL8","IFNB1","CXCL10","SOCS1","SOCS3",
           "IFIT2","IFIT3","IFIT5")

# save a figure together with the table behind it

save_with_table <- function(file, plot, data, ...) {
  ggsave(file.path(fig_dir, file), plot, bg = "transparent", ...)
  write_tsv(data, file.path(fig_dir, paste0(tools::file_path_sans_ext(file), ".tbl.gz")))
  message("wrote ", file)
}

## Peak -> gene map (canonical TSSs within 250 bp) -----------------------------
tss <- read_tsv(tss_file, col_names = c("chrom","start","end","act","info",
                                        "tsspos","strand","distance"),
                show_col_types = FALSE) |>
  mutate(pid  = paste0(chrom, ":", start, "-", end),
         gene = sapply(strsplit(info, ";"), `[`, 3))
peak_genes <- tss |> group_by(pid) |> summarize(genes = list(unique(gene)), .groups = "drop")
tier_of_peak <- sapply(peak_genes$genes, function(g)
  if (any(g %in% tier1)) "tier1" else if (any(g %in% tier2)) "tier2" else NA)
peak_tier <- setNames(tier_of_peak, peak_genes$pid)

## Per-pair scores, all timepoints, log2 distance bins -------------------------
all_pairs <- bind_rows(lapply(timepoints, function(s) {
  read_csv(file.path(codep_dir, paste0("codep_", s, ".csv.gz")), show_col_types = FALSE) |>
    mutate(timepoint = s)
})) |>
  mutate(bin_log2 = pmin(pmax(floor(log2(dist)), 8), 14),   # < 512 bp pooled; > 16 kb pooled (sparse fiber-length tail)
         tier_1 = peak_tier[motif_1],
         tier_2p = peak_tier[motif_2],
         group = case_when(
           (!is.na(tier_1) & tier_1 == "tier1") | (!is.na(tier_2p) & tier_2p == "tier1") ~ "LPS tier 1 (immediate-early)",
           (!is.na(tier_1) & tier_1 == "tier2") | (!is.na(tier_2p) & tier_2p == "tier2") ~ "LPS tier 2 (secondary)",
           TRUE ~ "Global"),
         timepoint = factor(timepoint, levels = timepoints))

bin_breaks <- 8:14
bin_labels <- c("<512bp", ">512bp", ">1kb", ">2kb", ">4kb", ">8kb", "\u226516kb")

message("pairs per group / timepoint:")
print(count(all_pairs, group, timepoint) |> pivot_wider(names_from = timepoint, values_from = n))

## Shuffled null, aggregated into the same log2 bins ---------------------------
## Only 100 bp-binned means exist for the shuffle; a count-weighted aggregation
## of those means into log2 bins is exact for the mean.
shuffle <- bind_rows(lapply(timepoints, function(s) {
  read_csv(file.path(codep_dir, paste0("bin_codep_", s, "_shuffle_avg.csv")),
           show_col_types = FALSE) |> mutate(timepoint = s)
})) |>
  mutate(bin_log2 = pmin(pmax(floor(log2(dist + 50)), 8), 14)) |>
  group_by(bin_log2) |>
  summarize(mean_score = sum(score * count) / sum(count), .groups = "drop") |>
  mutate(timepoint = "Shuffled", group = "Shuffled null")

## Figure 1: distance decay, global vs LPS-response, one line per timepoint ----
binned <- all_pairs |>
  mutate(group = ifelse(group == "Global", "Global", "LPS-response loci")) |>
  group_by(group, timepoint, bin_log2) |>
  summarize(mean_score = mean(score), sem = sd(score) / sqrt(n()), n = n(), .groups = "drop")

null_lines <- bind_rows(shuffle |> mutate(group = "Global"),
                        shuffle |> mutate(group = "LPS-response loci")) |>
  mutate(sem = NA, n = NA)
fig1_df <- bind_rows(binned, null_lines) |>
  mutate(timepoint = factor(timepoint, levels = c(timepoints, "Shuffled")),
         group = factor(group, levels = c("Global", "LPS-response loci")))

fig1 <- ggplot(fig1_df, aes(bin_log2, mean_score, color = timepoint)) +
  geom_ribbon(aes(ymin = mean_score - sem, ymax = mean_score + sem, fill = timepoint),
              alpha = 0.15, color = NA, na.rm = TRUE) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
  facet_wrap(~group) +
  scale_color_manual(values = tp_colors) + scale_fill_manual(values = tp_colors, guide = "none") +
  scale_x_continuous(breaks = bin_breaks, labels = bin_labels) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  xlab("Distance bins (bp between peaks)") +
  ylab("Mean binned codependency score") +
  theme_bw(11) + theme(legend.position = "top", legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
save_with_table("codep_distance_decay_global_vs_LPS.pdf", fig1, fig1_df, height = 4.5, width = 9)

## Figure 2: difference curves vs LPS_0 (analog of their Same - Opposite) ------
diff_df <- binned |>
  select(group, timepoint, bin_log2, mean_score) |>
  pivot_wider(names_from = timepoint, values_from = mean_score) |>
  pivot_longer(all_of(setdiff(timepoints, "LPS_0")), names_to = "timepoint",
               values_to = "mean_score") |>
  mutate(diff = mean_score - LPS_0,
         timepoint = factor(timepoint, levels = timepoints))

fig2 <- ggplot(diff_df, aes(bin_log2, diff, color = timepoint)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
  facet_wrap(~group) +
  scale_color_manual(values = tp_colors) +
  scale_x_continuous(breaks = bin_breaks, labels = bin_labels) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  xlab("Distance bins (bp between peaks)") +
  ylab("Mean binned codependency score\n(timepoint - LPS_0)") +
  theme_bw(11) + theme(legend.position = "top", legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
save_with_table("codep_difference_vs_LPS0.pdf", fig2, diff_df, height = 4.5, width = 9)

## Figure 3: the two LPS tiers separately --------------------------------------
tier_binned <- all_pairs |>
  filter(group != "Global") |>
  group_by(group, timepoint, bin_log2) |>
  summarize(mean_score = mean(score), sem = sd(score) / sqrt(n()), n = n(), .groups = "drop")

fig3 <- ggplot(tier_binned, aes(bin_log2, mean_score, color = timepoint)) +
  geom_ribbon(aes(ymin = mean_score - sem, ymax = mean_score + sem, fill = timepoint),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.6) +
  facet_wrap(~group) +
  scale_color_manual(values = tp_colors) + scale_fill_manual(values = tp_colors, guide = "none") +
  scale_x_continuous(breaks = bin_breaks, labels = bin_labels) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  xlab("Distance bins (bp between peaks)") +
  ylab("Mean binned codependency score") +
  theme_bw(11) + theme(legend.position = "top", legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
save_with_table("codep_distance_decay_LPS_tiers.pdf", fig3, tier_binned, height = 4.5, width = 9)

## Per-bin paired t-tests vs LPS_0 (their per-bin t-test table, paired) --------
paired <- read_csv(file.path(codep_dir, "codep_paired_across_time.csv.gz"),
                   show_col_types = FALSE) |>
  mutate(bin_log2 = pmin(pmax(floor(log2(dist)), 8), 14),
         tier_1 = peak_tier[motif_1],
         tier_2p = peak_tier[motif_2],
         group = ifelse((!is.na(tier_1)) | (!is.na(tier_2p)), "LPS-response loci", "Global"))

tests <- list()
for (grp in c("Global", "LPS-response loci")) {
  gdf <- filter(paired, group == grp)
  for (s in setdiff(timepoints, "LPS_0")) {
    for (b in sort(unique(gdf$bin_log2))) {
      bdf <- filter(gdf, bin_log2 == b)
      if (nrow(bdf) < 5) next
      t <- t.test(bdf[[paste0("score_", s)]], bdf$score_LPS_0, paired = TRUE)
      tests[[length(tests) + 1]] <- data.frame(
        group = grp, timepoint = s, bin_log2 = b, n_pairs = nrow(bdf),
        mean_diff = unname(t$estimate), t_statistic = unname(t$statistic),
        degrees_of_freedom = unname(t$parameter), p_value = t$p.value,
        conf_low = t$conf.int[1], conf_high = t$conf.int[2])
    }
  }
}
tests_df <- bind_rows(tests) |> group_by(group) |>
  mutate(p_adj_BH = p.adjust(p_value, method = "BH")) |> ungroup()
write_csv(tests_df, file.path(fig_dir, "codep_per_bin_paired_ttests_vs_LPS0.csv"))
message("wrote codep_per_bin_paired_ttests_vs_LPS0.csv")
message("significant bins (BH < 0.05):")
print(filter(tests_df, p_adj_BH < 0.05) |> arrange(group, timepoint, bin_log2), n = 50)
