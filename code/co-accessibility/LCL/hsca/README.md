1. Open the repository and create the log directory.

   ```bash
   cd /project/spott/cshan/fiber-seq
   mkdir -p results/logs
   ```

2. Complete the upstream LCL co-accessibility workflow before HSCA. Check that its required outputs exist.

   ```bash
   for file in \
     LCL_project/co-accessibility/CRE_co-accessibility_FIRE/universe/cre_universe.bed \
     LCL_project/co-accessibility/CRE_co-accessibility_FIRE/coaccess/LCL_coaccess_pairs.tsv.gz \
     LCL_project/co-accessibility/CRE_co-accessibility_FIRE/coaccess/LCL_pair_sample_counts.tsv.gz
   do
     test -s "$file" || { echo "Missing upstream output: $file"; break; }
   done
   ```

3. Submit the full HSCA run. The wrapper runs `LCL_hsca_pool.py` and `LCL_hsca_run.py`. It writes the four REF/ALT test tables plus the 27-sample pooled/noXCI comparisons and cohort/sex manifests under `LCL_project/co-accessibility/hsca/outputs/`.

   ```bash
   hsca_job=$(sbatch --parsable code/co-accessibility/LCL/hsca/LCL_run_hsca.sh)
   echo "$hsca_job"
   ```

4. Wait until the job reports `COMPLETED` with exit code `0:0`. Check the error log if it fails.

   ```bash
   sacct -j "$hsca_job" --format=JobID,State,ExitCode,Elapsed
   tail -n 30 "results/logs/LCL_hsca_${hsca_job}.err"
   ```

5. In an R session on a compute node, render the report after the HSCA job completes. This includes the molecule, clustering, ACF, and domain/CTCF figures. The notebook sources `hap_sele_access_plots.r` automatically. All tests and molecule figures use `LCL_sample_metatable_merged_phased_27samples.csv`; the 31-sample upstream data supply coordinates and per-sample counts only. Re-run the wrapper and render all figures after this cohort change; existing PDFs are not updated by editing the scripts.

   ```r
   rmarkdown::render(
     "/project/spott/cshan/fiber-seq/code/co-accessibility/LCL/hsca/LCL_hsca.Rmd",
     output_dir = "/project/spott/cshan/fiber-seq/docs",
     params = list(make_molecule_plots = TRUE)
   )
   ```

6. Open `/project/spott/cshan/fiber-seq/docs/LCL_hsca.html`.

REF/ALT tests use the same coverage threshold throughout: at least 10 REF reads and
10 ALT reads per donor and cCRE/pair. Pooled tests sum only donor/unit rows meeting
that threshold. H1/H2 tests and FIRE haplotype-peak intersections are no longer run.
Phasing remains available for assigning REF/ALT labels in molecule figures.

Donor genotypes in `pair_genotype_calls()` now come from the donor's `cell_line`
column in `phased_vcf_file`, at the exact chromosome, position, REF and ALT.
Only PASS biallelic SNVs qualify. Diploid 0/1 and 1/0 are heterozygous whether
phased or unphased; 0/0 and 1/1 are homozygous. Missing or ineligible records
are uncalled, never inferred from read counts. REF/ALT read counts remain
diagnostics and read bases still label individual molecules. The statistical
tests retain their separate minimum of 10 reads per allele.

For this genotype-source change, rerun the report
`haplotype_selective_chromatin_accessibility_LCL.Rmd` with
`make_molecule_plots = TRUE`. This refreshes the domain/CTCF figures and genotype
summaries for CLASP1, SLFN5, LZIC and LINC00462, plus the GTEx lead-eQTL m6A
table. In an existing R session, source `hap_sele_access_plots.r` again before
rerunning those chunks. The step-5 statistical tables already use VCF genotypes:
neither `LCL_run_hsca.sh` nor upstream co-accessibility needs rerunning solely
for this change. The separate phased-molecule/clustering panels retain their
VCF phasing requirement for assigning alleles; that is a molecule-assignment
filter, not a read-count genotype call.
