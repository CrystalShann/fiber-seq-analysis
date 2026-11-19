# Molecule-cluster entropy


Run `bash code/pausing/PolII/run_gene_examples.sh code/entropy/molecule_cluster/molecule_cluster_entropy.R`.
Outputs go exclusively to `LCL_project/entropy/molecule_cluster_entropy/{tables,figures}`;
the calculation writes tables only. The shared `../entropy_features.rmd` reads these
tables to generate PDFs. Run `bash code/entropy/nucleosome/run_entropy_features.sh`
after calculation. This calculation does not read or rewrite nucleosome entropy outputs.

The script uses original `recurrent_patterns/tables/<group>/<gene>_read_clusters.tsv`
assignments for the example gene list. These represent the original full-span promoter
clustering cohort (TSS +/-1000 bp). No new clustering or resampling is performed.
The downstream recurrent-family `eligible` flag is not a read eligibility filter:
small original clusters are retained. The saved read-feature cohort is used only
to identify missing assignments; assigned reads outside that cohort are also retained.
Reads are deduplicated within genes using sample-qualified IDs. Identical repeated
labels count once, conflicting valid labels are excluded as ambiguous, and missing
or invalid labels are excluded from assigned-read denominators. A valid label with
an additional invalid duplicate is retained. Unresolvable read identities cannot
enter unique-read denominators. Missing assignment files are not reconstructed from families.

Raw entropy is `-sum(p*log2(p))`, normalized entropy is `H/log2(K)` for K>=2,
and effective states are `2^H`. For K=1, H=0 and normalized entropy is NA; for
no assigned reads all entropy measures are NA. Proportions are validated to sum
to one.

Tables include gene estimates, cluster counts/proportions, group summaries and
group tests. Wilcoxon pairwise p-values are BH-adjusted across the three group
contrasts separately for each outcome; Kruskal-Wallis p-values are unadjusted. Tests use
one observation per gene and do not adjust for depth or clustering resolution.
PDFs show raw/normalized entropy by group and entropy versus cluster count.
No session information or run metadata are generated.

## Findings from the current LCL cohort

149 of 150 genes have valid assignments: 75,276 assigned gene/fiber observations
across 795 original clusters. SLC44A4 has no valid assignments. CLEC12B is included
because this analysis does not require nucleosome calls. One gene has only one
cluster, so normalized entropy is available for 148 genes.

| Group | Genes | Median H (bits) | Median normalized H | Median effective states | Median K | Median assigned fibers |
|---|---:|---:|---:|---:|---:|---:|
| Housekeeping | 50 | 2.4491 | 0.9727 | 5.4607 | 6 | 544 |
| Cell-type-specific | 50 | 2.2519 | 0.9583 | 4.7632 | 5 | 516.5 |
| Nonexpressed | 49 | 2.2573 | 0.9722 | 4.7809 | 5 | 558 |

Raw entropy differs across groups (Kruskal-Wallis p=0.00599). Housekeeping genes
have higher raw entropy than cell-type-specific genes (pairwise BH-adjusted
p=0.0226) and nonexpressed genes (p=0.00976). Cell-type-specific and nonexpressed
genes are not distinguishable in this comparison (p=0.395).

Normalized entropy shows no clear group difference (Kruskal-Wallis p=0.235;
all pairwise BH-adjusted p>=0.295). High median normalized entropies indicate
relatively even read proportions across each gene's observed clusters.

Raw entropy strongly tracks observed cluster count (Spearman rho=0.951), but
has little association with assigned-read depth (rho=0.061, p=0.463).
Normalized entropy has little association with cluster count (rho=-0.047,
p=0.573). The higher housekeeping raw entropy is therefore consistent with
more observed states rather than a strong difference in cluster evenness.
This is descriptive: normalized entropy does not remove clustering-resolution
or read-depth effects, and these cluster states are not direct counts of
biochemically distinct promoter states.
