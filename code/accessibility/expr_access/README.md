# TSS m6A metaprofiles by expression level

Metaplots of chromatin accessibility (m6A) around protein-coding TSS,
stratified by RNA expression, per LPS timepoint (0/5/10/15 min). All TSS are
aligned and strand-oriented (upstream negative, sense direction to the right).

Pipeline (driver: `run_expr_access.sh`, sbatch-able or direct; extra args go to
the python step, e.g. `--chrom chr21` for a test):

1. `01_expression_bins.R` — unfiltered STAR ReadsPerGene counts (reverse
   -stranded col 4) for the 0-15 min samples -> TPM (union-exon lengths,
   gencode v46) -> mean over the 12 samples -> bins: mean TPM < 1 =
   not_expressed, else expressed quartiles Q1..Q4. TSS = Ensembl_canonical
   transcript TSS (v46), protein-coding only, chr1-22/X/Y.
   Writes `tables/tss_expression_bins.tsv`.
2. `02_tss_m6a_profiles.py` — per timepoint, tabix-fetches the ft extract
   m6a_by_chr reads over TSS +/- 1 kb, 10 bp bins on the strand-aware axis;
   m6a = n_mod / n_cov with n_cov

Writes
   `tables/tss_m6a_profile_by_expr_bin.tsv.gz` and
   `tables/tss_m6a_gene_totals.tsv.gz`.
3. `03_plot_expr_access.R` — metaplots colored by expression bin (faceted by
   timepoint) and by timepoint (faceted by bin) into `plots/`.

Outputs under `/project/spott/cshan/fiber-seq/macrophage_project/expr_access/`.
