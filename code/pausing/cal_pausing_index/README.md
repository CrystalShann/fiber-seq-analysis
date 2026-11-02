# Pol II pausing index

Calculated by `pause_index.R` for GM12878 (GRCh38)

TSSR stands for TSS region: the window of 350 bp around a gene’s transcription start site, from 50 bp upstream to 299 bp downstream

## Calculation steps

1. Load RefSeq transcripts on chr1–22 and chrX; drop genes on multiple chromosomes or strands

2. Collapse identical read 5′ ends and strands within each replicate; pool ChIP and control replicates separately

3. Count read 5′ ends in each strand-aware TSS region (TSS −50 to +299 bp)
  - rpm/bp = reads per million, per base pair = rpm/bp = reads / len / (n / 1e6)
  
4. Calculate `density = reads / region_length / (library_reads / 1e6)` and subtract control density from ChIP density

5. Select the strongest TSS signal per gene, then the longest isoform sharing that TSS.

6. Calculate the same corrected density for the strand-aware gene body (TSS +300 bp to TES +3000 bp)
  - promoter region (TSSR): -50 - TSS + 299
  - gene body: TSS + 300 to TES + 3000
  - TSSR signal is the background-subtracted Pol II density at the promoter
  - Gene-body signal is the background-subtracted Pol II density along the transcribed region
  - PI = TSSR signal / gene-body signal. A high value means Pol II is concentrated at the promoter relative to the body
  
7. Calculate `PI = TSS-region signal / gene-body signal` if TSS signal ≥0.001 rpm/bp and body signal >0; otherwise assign `NA`. Mark PI ≥2 as paused (missing PI is not paused)
  - PI is calculated only when both the promoter and gene body signal is positive 
  - a gene is active if its TSSR signal is at least MIN_DENS = 0.001

8. Average two RNA-seq FPKM replicates, matching RefSeq symbols through GENCODE gene IDs

9. Flag the analysis set: autosomal protein-coding isoforms, FPKM >1, length ≥1 kb, no other RefSeq gene within 3 kb of the gene span, and valid PI. Divide its PI values into low/mid/high tertiles

Genes failing analysis-set filters remain in the output. No H3K4me3 fallback is used for weak TSS signal.

## Input files

Paths below are relative to `/project/spott/cshan/annotations/`.

RefSeq transcripts and GENCODE gene-ID mapping:

```text
hg38.refGene.gtf.gz
gencode.v46.annotation.gtf.gz
```

Pol II ChIP (ENCSR000EAD, seven replicates):

```text
LCL_public_data/Pol2_chip/ENCFF441GAS.bed.gz
LCL_public_data/Pol2_chip/ENCFF306JCA.bed.gz
LCL_public_data/Pol2_chip/ENCFF248WXH.bed.gz
LCL_public_data/Pol2_chip/ENCFF590GEB.bed.gz
LCL_public_data/Pol2_chip/ENCFF944QSC.bed.gz
LCL_public_data/Pol2_chip/ENCFF811KVZ.bed.gz
LCL_public_data/Pol2_chip/ENCFF603TJG.bed.gz
```

Input control (ENCSR000EAF, two replicates):

```text
LCL_public_data/Pol2_chip/ENCFF863XTX.bed.gz
LCL_public_data/Pol2_chip/ENCFF278NPZ.bed.gz
```

RNA-seq expression (ENCSR820PHH, two replicates):

```text
LCL_public_data/gene_quan_ENCSR820PHH/ENCFF873VWU.gene_quan.tsv
LCL_public_data/gene_quan_ENCSR820PHH/ENCFF345SHY.gene_quan.tsv
```

## Output files

- `/project/spott/cshan/fiber-seq/LCL_project/pausing/pause_index.tsv`: one row per retained RefSeq gene, including selected isoform, coordinates (1-based inclusive), read counts, corrected densities, PI, paused status, FPKM, filter flags, and PI tertile.
- Slurm logs: `/project/spott/cshan/fiber-seq/results/logs/pause_index_<jobid>.out` and `.err`.

## Run

From `/project/spott/cshan/fiber-seq`:

```bash
sbatch code/pausing/cal_pausing_index/run_pause_index.sh
```
