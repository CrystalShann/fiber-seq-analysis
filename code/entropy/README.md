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

Rao diversity (`diversity/`) is a separate analysis using fiber-to-fiber m6A distances.

## Run order

From `/project/spott/cshan/fiber-seq`, after the PolII read-feature and
recurrent-pattern steps (`code/pausing/PolII`) have written their tables:

```bash
# 1. nucleosome entropy + short-footprint fractions -> LCL_project/entropy/tables
sbatch code/entropy/nucleosome/run_entropy.sh --by-cluster
# 2. molecule-cluster entropy -> LCL_project/entropy/molecule_cluster_entropy
bash code/pausing/PolII/run_gene_examples.sh code/entropy/molecule_cluster/molecule_cluster_entropy.R
# 3. Gini-Simpson + chromatin features (needs 2) -> LCL_project/entropy/gini_simpson_molecule_cluster
bash code/pausing/PolII/run_gene_examples.sh code/entropy/gini_simpson_molecule_cluster/simpson_functions.R
# 4. optional: Rao diversity -> LCL_project/entropy/diversity
bash code/pausing/PolII/run_gene_examples.sh code/entropy/diversity/run_diversity.R
# 5. notebook: tables + figures (needs 1-3)
bash code/entropy/nucleosome/run_entropy_features.sh
```

Step 1 also runs step 5 at the end; on a fresh run that pass fails until steps 2-3
exist, so run step 5 again last. Details: [nucleosome/](nucleosome/README.md),
[molecule_cluster/](molecule_cluster/README.md).
