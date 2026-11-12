# promoter_annotation

## 04a_puffin_motif_scores.py

For each curated Pol II example gene, this script runs Puffin (Dudnyk et al., *Science* 2024,
"Sequence basis of transcription initiation in the human genome";
https://github.com/jzhoulab/puffin) on the DNA around the TSS your fibers support. It turns
Puffin's per-base breakdown of that sequence into one score per motif at the TSS, and labels
each gene by its strongest motif.

Run:

```
/project/spott/cshan/envs/puffin/bin/python \
    code/pausing/PolII/promoter_annotation/04a_puffin_motif_scores.py
```

Input: `LCL_project/pausing/PolII/read_features/tables/gene_features.tsv` (genes with a fiber
`tss`, 1-based) and `annotations/hg38.fa`. Puffin source and weights: `envs/puffin/puffin_src`.

### Puffin in brief

Puffin predicts transcription initiation signal at every base. It is built from three parts
whose effects add on a log scale, computed on both strands:

| Part | Layer | Size | Meaning |
|---|---|---|---|
| Motifs | `conv` | 10 filters × 2 orientations = 20, 51 bp each | TATA, YY1, SP, ETS, NFY, CREB, NRF1, ZNF143, U1 snRNP, Long Inr |
| Initiators | `conv_inr` | 15 bp | short Inr-like sequence at the TSS |
| Trinucleotides | `conv_sim` | 3 bp | broad composition effects, mostly CpG-related |

How a part makes its prediction:
- A filter's match strength along the sequence is its **activation**.
- A learned **effect curve**, up to 601 bp wide, spreads that activation out to the bases around
  it. Most motifs boost initiation at a set distance from themselves and suppress it close by.
  For example, TATA works about 30 bp upstream and YY1 just downstream.
- The curves from all parts are added, and a softplus turns the sum into the prediction
  (roughly log(signal + 1)).

Puffin needs 325 bp of extra sequence on each side of the region it predicts: the 601-bp curve
plus the 51-bp filter.

### What the script does

1. **Sequence.** Takes hg38 [TSS − 825, TSS + 825) in the gene's orientation. Minus-strand genes
   are reverse-complemented and shifted by 1 bp, as in Puffin's own command-line tool. Puffin
   returns 1000 bp (−500..+499) with the TSS at index 500.
2. **Puffin.** Calls `Puffin.interpret(seq, "FANTOM_CAGE")`. This gives, for every base, the
   prediction, each motif's activation and effect, and the summed motif, initiator,
   trinucleotide and total effects.
3. **Score per motif.** Averages each motif's effect over TSS −10..+9. Each base is weighted by
   Puffin's prediction there (log scale): Σ w·effect / Σ w. This reproduces
   `hg38_top16000_tss_groups.tsv` from the Puffin manuscript (`Motif_contribution.ipynb`) to
   within its rounding (max abs difference 0.0005). The four summed effects are averaged the same
   way.
4. **TATA box.** Finds the strongest TATA+ match between −60 and −10.
5. **Class.** Picks the highest of the 14 class motifs (U1 snRNP left out), drops the strand, and
   renames SP to SP1.
6. **Profiles.** Saves the four summed effects at every position from −500 to +499.

### Outputs

Under `LCL_project/pausing/PolII/promoter_annotation/tables/`.

#### `04a_puffin_motif_scores.tsv`: one row per gene (150)

| Column | Meaning |
|---|---|
| `gene_name` | gene symbol |
| `chrom` | chromosome |
| `puffin_tss0` | TSS used, 0-based (fiber `tss` − 1) |
| `strand` | gene strand; the sequence is oriented so the gene's transcription runs left to right |
| `TATA +`, `YY1 +` | TSS score (step 3) for that motif in the sense orientation |
| `SP +`, `SP -`, `ETS +`, `ETS -`, `NFY +`, `NFY -`, `CREB +`, `CREB -`, `NRF1 +`, `NRF1 -`, `ZNF143 +`, `ZNF143 -` | TSS score for each orientation of these motifs. `+` = motif on the gene's strand, `-` = on the opposite strand |
| `U1 snRNP +` | TSS score for U1 snRNP (5′ splice-site-like). It is reported but never used for the class |
| `total_effect_tss` | weighted TSS average of Puffin's "Sum of total effect" (motifs + initiators + trinucleotides) |
| `motif_effect_tss` | same for "Sum of motif effect" (all 20 motif filters except Long Inr) |
| `initiator_effect_tss` | same for "Sum of initiator effect" (short Inr + Long Inr+, centred on the 1 kb window) |
| `trinucleotide_effect_tss` | same for "Sum of trinucleotide effect" (centred on the 1 kb window) |
| `tata_center` | position relative to the TSS (bp, −60..−10) of the strongest TATA+ activation, at the motif centre |
| `tata_activation` | TATA+ activation at `tata_center` (min-subtracted over the 1 kb window) |
| `dominant_motif_excl_U1` | the highest of the 14 class motif columns (strand kept) |
| `dominant_contribution_excl_U1` | that motif's score |
| `group_excl_U1` | promoter class: `dominant_motif_excl_U1` without strand, SP → SP1 (TATA, YY1, SP1, ETS, NFY, CREB, NRF1, ZNF143) |

All effect and score columns are additive contributions to Puffin's log-scale FANTOM CAGE
prediction. Positive means the motif pushes initiation at the TSS up; negative means it pushes it
down.

#### `04a_puffin_effect_profiles.tsv.gz`: one row per gene × position (150 × 1000)

| Column | Meaning |
|---|---|
| `gene_name` | gene symbol |
| `group` | gene group from `gene_features.tsv` (`housekeeping`, `cell_type_specific`, `not_expressed`) |
| `position` | bp relative to the TSS, −500..+499, in the gene's orientation (negative = upstream) |
| `total`, `motif`, `initiator`, `trinucleotide` | Puffin's per-base summed effects (not weighted), on the same log scale |

