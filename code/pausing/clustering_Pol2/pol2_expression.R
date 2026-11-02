#!/usr/bin/env Rscript
# Pausing index and RNA-seq expression of the TSSs in each cluster, after the
# report stage. Expression = ENCSR820PHH RSEM gene TPM, mean of the two
# replicates, joined by Ensembl gene ID (versions ignored, PAR_Y copies dropped).
# Each TSS counts once per cluster, as for the cluster PI. Reads tables/ and
# writes tables/cluster_expression.tsv and plots/pi_expression_by_cluster.pdf.
script <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1L])
here <- dirname(normalizePath(script))
source(file.path(here, "pol2_clustering_functions.R"))
source(file.path(here, "pol2_clustering_plots.R"))
cfg <- yaml::read_yaml(file.path(here, "config.yaml"))
pol2_load_helpers(cfg$project)
tables <- file.path(cfg$output, "tables")

expr <- rbindlist(lapply(cfg$expression, function(f) fread(f, select = c("gene_id", "TPM"))))
expr <- expr[grepl("^ENSG", gene_id) & !grepl("_PAR_Y$", gene_id),
             .(tpm = mean(TPM), n_rep = .N), by = .(gene_key = sub("\\..*$", "", gene_id))]
stopifnot(all(expr$n_rep == length(cfg$expression)))

summary <- fread(file.path(tables, "cluster_summary.tsv"))
m <- fread(file.path(tables, "molecule_assignments.tsv.gz"))
tss <- unique(m[, .(cluster, tss_id, gene_id, gene_name, pi)])
tss[, tpm := expr$tpm[match(sub("\\..*$", "", gene_id), expr$gene_key)]]
tss[, cluster := factor(cluster, levels = summary$cluster)]
message(sum(is.na(tss$tpm)), " of ", nrow(tss), " cluster/TSS pairs have no TPM")

spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(list(rho = NA_real_, p = NA_real_, n = sum(ok)))
  r <- suppressWarnings(stats::cor.test(x[ok], y[ok], method = "spearman"))
  list(rho = unname(r$estimate), p = r$p.value, n = sum(ok))
}
d <- tss[is.finite(tpm)]
out <- d[, c(list(n_tss = .N, median_tpm = median(tpm), frac_tpm_ge1 = mean(tpm >= 1),
                  median_pi = median(pi, na.rm = TRUE)),
             setNames(spearman(pi, tpm), c("rho_pi_tpm", "p_rho", "n_pi_tpm"))), by = cluster]
# TPM of each cluster's TSSs against all other clusters' TSSs.
out[, p_tpm_vs_rest := vapply(as.character(cluster), function(cl)
  stats::wilcox.test(d[cluster == cl, tpm], d[cluster != cl, tpm])$p.value, numeric(1))]
out[, `:=`(fdr_rho = adjust_bh(p_rho), fdr_tpm_vs_rest = adjust_bh(p_tpm_vs_rest))]
setorder(out, cluster)
fwrite(out, file.path(tables, "cluster_expression.tsv"), sep = "\t")
kw <- stats::kruskal.test(tpm ~ cluster, data = d)
all_rho <- spearman(d$pi, d$tpm)
message(sprintf("TPM across clusters: Kruskal-Wallis p = %.3g; PI vs TPM over all clustered TSSs: rho = %.3f (n = %d, p = %.3g)",
                kw$p.value, all_rho$rho, all_rho$n, all_rho$p))
print(out)

colors <- cluster_id_palette(levels(d$cluster))
boxes <- ggplot(d, aes(cluster, log10(tpm + 1), fill = cluster)) +
  geom_boxplot(outlier.size = 0.3) + scale_fill_manual(values = colors, guide = "none") +
  labs(x = "Cluster (median PI, lowest first)", y = "log10(TPM + 1)", title = "Expression of each cluster's TSSs",
       subtitle = sprintf("ENCSR820PHH, mean of 2 replicates; Kruskal-Wallis p = %.2g", kw$p.value)) +
  theme_fiberseq("classic", base_size = 10)
labels <- out[, .(cluster, label = sprintf("rho = %.2f\nn = %d", rho_pi_tpm, n_pi_tpm))]
scatter <- ggplot(d[is.finite(pi) & pi > 0], aes(log10(pi), log10(tpm + 1), colour = cluster)) +
  geom_point(size = 0.5, alpha = 0.5) + facet_wrap(~cluster, nrow = 1) +
  geom_text(data = labels, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = -0.1, vjust = 1.2, size = 2.8) +
  scale_colour_manual(values = colors, guide = "none") +
  labs(x = "log10(pausing index)", y = "log10(TPM + 1)", title = "Pausing index vs expression within each cluster (Spearman)") +
  theme_fiberseq("bw", base_size = 10)
save_figure(cowplot::plot_grid(boxes, scatter, ncol = 1, rel_heights = c(1, 0.9)),
            file.path(cfg$output, "plots", "pi_expression_by_cluster.pdf"), width = 14, height = 8)
