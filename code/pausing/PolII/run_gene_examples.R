# Plots every curated housekeeping / cell type specific gene (curated_example_genes() in
# pol2_pausing_functions.r) into LCL_project/pausing/PolII/plots/example_reads/<group>/<gene>/.
# sbatch --account=pi-spott --partition=caslake --cpus-per-task=16 --mem=64G --time=03:00:00 \
#   --output=/project/spott/cshan/fiber-seq/LCL_project/pausing/PolII/run.log \
#   --export=ALL,PAUSING_WORKERS=8 run_gene_examples.sh
dyn.load("/project/spott/cshan/envs/dimelo/lib/libcurl.so.4", local = FALSE)
source("/project/spott/cshan/fiber-seq/code/pausing/PolII/pol2_pausing_functions.r")
data.table::setDTthreads(1L)
run_curated_examples(workers = as.integer(Sys.getenv("PAUSING_WORKERS", "1")))
