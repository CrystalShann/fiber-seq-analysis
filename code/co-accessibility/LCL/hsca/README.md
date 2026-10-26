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

3. Submit the full HSCA run. The wrapper runs `LCL_hsca_intersect.py` followed by `LCL_hsca_run.py` and writes the six tables under `LCL_project/co-accessibility/hsca/outputs/`.

   ```bash
   hsca_job=$(sbatch --parsable code/co-accessibility/LCL/hsca/LCL_run_hsca.sh)
   echo "$hsca_job"
   ```

4. Wait until the job reports `COMPLETED` with exit code `0:0`. Check the error log if it fails.

   ```bash
   sacct -j "$hsca_job" --format=JobID,State,ExitCode,Elapsed
   tail -n 30 "results/logs/LCL_hsca_${hsca_job}.err"
   ```

5. In an R session on a compute node, render the report after the HSCA job completes. This includes the molecule, clustering, ACF, and domain/CTCF figures. The notebook sources `fire_plots.R` and `hap_sele_access_plots.r` automatically.

   ```r
   rmarkdown::render(
     "/project/spott/cshan/fiber-seq/code/co-accessibility/LCL/hsca/LCL_hsca.Rmd",
     output_dir = "/project/spott/cshan/fiber-seq/docs",
     params = list(make_molecule_plots = TRUE)
   )
   ```

6. Open `/project/spott/cshan/fiber-seq/docs/LCL_hsca.html`.
