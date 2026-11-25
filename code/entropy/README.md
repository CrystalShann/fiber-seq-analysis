# Promoter entropy and diversity (LCL Fiber-seq)

Three measures of promoter chromatin heterogeneity for the 150 PolII example genes
(housekeeping, cell-type-specific, not expressed), all computed from the saved PolII
full-span cohort (fibers spanning TSS ±1 kb).

## Methods: reference vs. this adaptation

**Nucleosome entropy** — Li et al. 2026, *NAR* 54:gkaf1385 (`references/`).
- Reference: yeast MNase-seq. Nucleosome scores over the 500 bp downstream of the TSS
  are normalized to proportions, $p_i = nuc(i) / \sum nuc$, and their Shannon entropy
  $H = -\sum p_i \log_2 p_i$ is used as a proxy for nucleosome disorganization.
- Here: $nuc(i)$ is single-molecule occupancy, $o_x = n_x / c_x$ (fibers with a
  130-160 bp fibertools nucleosome over base $x$ / fibers covering $x$), at TSS +1..+500 bp,
  1 bp resolution, no binning or smoothing. Computed per gene and per original
  Leiden cluster; requires all 500 bases covered.

**Gini-Simpson diversity** — de Langen et al. 2023, *Cell Genomics* 3:100411 (`references/`).
- Reference: per Pol II consensus peak, $G = 1 - \sum p_i^2$ over the proportions of
  peaks from each tissue/cell annotation; near 1 = shared across conditions, 0 = specific.
- Here: $p_k$ is the fraction of a promoter's fibers in each of its own Leiden m6A
  clusters, so $G$ is the probability that two random fibers are in different
  clusters. Also reported: $G / (1 - 1/K)$ (evenness) and $1/\sum p_k^2$; $G$ is
  correlated with promoter chromatin features (+1 nucleosome, m6A, short footprints).

**Molecule-cluster Shannon entropy** — the same cluster proportions $p_k$ as above,
summarized with $H$, $H / \log_2 K$ and $2^H$ (effective number of clusters).

Cluster Manhattan diversity (`diversity/`) uses abundance-weighted mean absolute
profile differences across consensus Leiden clusters (10 neighbors, resolution 1).
Distances are divided by all 2,001 genomic bases. The pipeline exports this Q as
`manhattan_between` and the main diversity score as
`cluster_quadratic_diversity = 1 / (1 - manhattan_between)`.
The notebook reports gene groups, motifs, chromatin-feature associations and rankings.

## Run order

From `/project/spott/cshan/fiber-seq`, after the PolII read-feature and
recurrent-pattern steps (`code/Pol2/PolII`) have written their tables:

```bash
# 1. nucleosome entropy + short-footprint fractions -> LCL_project/entropy/nucleosome_entropy/tables
sbatch code/entropy/nucleosome/run_entropy.sh --by-cluster
# 2. molecule-cluster entropy -> LCL_project/entropy/molecule_cluster_entropy
bash code/Pol2/PolII/run_gene_examples.sh code/entropy/molecule_cluster/molecule_cluster_entropy.R
# 3. Gini-Simpson + chromatin features (needs 2) -> LCL_project/entropy/gini_simpson_molecule_cluster
bash code/Pol2/PolII/run_gene_examples.sh code/entropy/gini_simpson_molecule_cluster/simpson_functions.R
# 4. optional: cluster Manhattan diversity -> LCL_project/entropy/diversity
#    tables + feature PDFs: run the chunks of diversity_features.rmd
#    (compute-gene-diversity recomputes gene_diversity.tsv when missing or RECOMPUTE_DIVERSITY = TRUE)
# 5. notebook: tables + figures (needs 1-3): run the chunks of entropy_features.rmd
```

Details: [nucleosome/](nucleosome/README.md),
[molecule_cluster/](molecule_cluster/README.md).
