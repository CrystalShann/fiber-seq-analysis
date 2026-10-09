    source("/project/spott/cshan/fiber-seq/code/parsing_functions/plotting_functions.r", local = TRUE)
    for (i in seq_len(nrow(examples))) {
      ex <- examples[i, , drop = FALSE]
      result <- reference_result(ex)
      show_reads(result)
      a <- result$assignments
      cl_cols <- cluster_palette(levels(a$cluster))
      res_ex <- readRDS(file.path(example_dir(ex), "retained_molecule_tracks.rds"))
      stopifnot(setequal(as.character(res_ex$rids_df$RID), a$RID))
      tracks <- molecule_tracks(res_ex)
      rows_C <- a[order(a$cluster, a$time, a$start, a$RID), ]
      res_cl <- add_read_groups(res_ex, labels = setNames(as.character(a$cluster), a$RID),
                                group_col = "cluster", group_levels = levels(a$cluster), verbose = FALSE)
      rows_C2 <- stack_rows(rows_C, "cluster")
      n_C2 <- table(rows_C2$cluster)
      lab_C2 <- setNames(sprintf("%s\nn = %d", names(n_C2), as.integer(n_C2)), names(n_C2))
      pC2_mol <- plot_molecules(rows_C2, tracks, ex, "cluster", lab_C2,
                                bar_cols = c("cluster", "time"),
                                bar_palette = c(cl_cols, time_cols), legend_extra = time_levels)
      pC2_accessibility <- local({
        .profile_args <- list(pileup = res_cl$pileup,
      ex = ex,
      group_col = "cluster",
      group_levels = levels(a$cluster),
      group_cols = cl_cols,
      n_reads = n_C2)
        pileup <- .profile_args$pileup
        ex <- .profile_args$ex
        group_col <- .profile_args$group_col
        group_levels <- .profile_args$group_levels
        group_cols <- .profile_args$group_cols
        n_reads <- .profile_args$n_reads
            stopifnot(all(group_levels %in% names(n_reads)))
            pileup$accessibility_group <- factor(as.character(pileup[[group_col]]),
                                                  levels = group_levels)
            stopifnot(!anyNA(pileup$accessibility_group))
            strip <- setNames(sprintf("%s\nn = %d", group_levels,
                                      as.integer(n_reads[group_levels])), group_levels)
            plot_group_profile(pileup, group_col = "accessibility_group", style = "area",
              y_col = "smooth_frac", mapping = aes(pos, smooth_frac),
              highlight = list(xmin = ex$start, xmax = ex$end, fill = "grey88", position = "under"),
              layers = list(list(geom = "area", mapping = aes(fill = accessibility_group), alpha = 0.9),
                list(geom = "line", linewidth = 0.25, colour = "grey15"))) +
              facet_grid(accessibility_group ~ ., labeller = as_labeller(strip)) +
              scale_fill_manual(values = group_cols, guide = "none") +
              scale_x_continuous(breaks = window_breaks(ex), labels = scales::comma) +
              scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
              labs(x = NULL, y = "m6A fraction (+/- 5 bp mean)") +
              theme_fiberseq("classic", base_size = 9, overrides = example_theme_settings) +
              theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                    axis.line.x = element_blank())
          
      })
      pC2 <- plot_footprint_panels(pC2_accessibility, pC2_mol, ex, n_bars = 2L,
        n_groups = nlevels(a$cluster), n_molecules = nrow(rows_C2),
        subtitle = sprintf(paste0("Accessibility by Leiden cluster; %s motif positions; tested FIRE; ",
                                  "all %d full-span molecules by cluster\n",
                                  "Left bars = cluster, time; molecules ordered by time then read start; ",
                                  "dashed lines = tested region"), ex$selection_TF, nrow(rows_C2)))
      h_C2 <- pC2$height
      path_C2 <- example_figure_path(ex, "panel3C", "cluster_footprints")
      ggsave(path_C2, pC2$plot, width = 8.5, height = h_C2, limitsize = FALSE)
      print(head(res_cl$pileup))
      print(head(rows_C2[, c("RID", "cluster", "time", "start", "y")]))
      print(n_C2)

      h_C1 <- min(14, 3 + 0.02 * nrow(rows_C))
      h_C <- max(h_C1, h_C2)
      g_C1 <- grid::grid.grabExpr(plot_read_heatmap(result, layout = "fire", rows = rows_C, window = ex,
      cluster_colors = cl_cols, timepoint_cols = time_cols, ticks = window_breaks(ex),
      main = example_title(ex), draw = TRUE, use_raster = FALSE),
                                  width = 9, height = h_C)
      pC <- cowplot::plot_grid(
        cowplot::plot_grid(NULL, NULL, nrow = 1, rel_widths = c(9, 8.5), labels = c("C1", "C2")),
        cowplot::plot_grid(g_C1, patchwork::patchworkGrob(pC2$plot), nrow = 1, rel_widths = c(9, 8.5)),
        ncol = 1, rel_heights = c(0.3, h_C))
      path_C <- example_figure_path(ex, "panel3BC", "combined")
      ggsave(path_C, pC, width = 17.5, height = h_C + 0.3, limitsize = FALSE)
      cat("Saved", path_C2, "and", path_C, "\n")
    }
    