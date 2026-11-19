#!/usr/bin/env python
# Position of each gene's dominant Puffin motif around dominant_TSS, for the motif annotation of
# the representative CAGE peak figures (02_promoter_shape_figures.Rmd, section 9).
#
# Run: /project/spott/cshan/envs/puffin/bin/python code/promoter_shape/03_puffin_dominant_motif_site.py
#
# 04a_puffin_motif_scores.tsv names the dominant motif (largest contribution at the TSS, U1 snRNP
# excluded) but only stores a position for TATA+ and YY1+. Here Puffin is rerun on the same sequence
# (04a's tss_sequence, TSS = dominant_TSS - 1) and the dominant motif's site = TSS-relative position
# of its peak motif activation: TATA+ in -60..-10 and YY1+ in 0..+50 (04a's ranges, so TATA / YY1
# match tata_center / yy1_center), every other motif in MOTIF_SEARCH.
#
# Output: LCL_project/promoter_shape/puffin_dominant_motif_sites.tsv, one row per gene

import importlib.util
import os

import numpy as np
import pandas as pd
import pysam

PROJECT = "/project/spott/cshan/fiber-seq"
SHAPE = f"{PROJECT}/LCL_project/promoter_shape/promoter_shape_annotations.tsv"
SCORES = f"{PROJECT}/LCL_project/pausing/PolII/promoter_annotation/tables/04a_puffin_motif_scores.tsv"
OUT = f"{PROJECT}/LCL_project/promoter_shape/puffin_dominant_motif_sites.tsv"
MOTIF_SEARCH = (-100, 50)  # TSS-relative range for motifs other than TATA+ / YY1+

spec = importlib.util.spec_from_file_location(
    "puffin04a", f"{PROJECT}/code/pausing/PolII/promoter_annotation/04a_puffin_motif_scores.py")
p04a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p04a)  # loads Puffin and chdirs to its source folder
SEARCH = {"TATA +": p04a.TATA_SEARCH, "YY1 +": p04a.YY1_SEARCH}


def main():
    shape = pd.read_csv(SHAPE, sep="\t")
    genes = shape[shape.contains_dominant_TSS].merge(
        pd.read_csv(SCORES, sep="\t")[["gene_name", "dominant_motif_excl_U1", "dominant_contribution_excl_U1"]],
        on="gene_name")
    fa = pysam.FastaFile(p04a.GENOME_FA)
    net = p04a.Puffin()
    rows = []
    for g in genes.itertuples():
        motif = g.dominant_motif_excl_U1
        lo, hi = SEARCH.get(motif, MOTIF_SEARCH)
        df = net.interpret(p04a.tss_sequence(fa, g.chr, g.dominant_TSS - 1, g.strand), "FANTOM_CAGE")
        act = df.loc[motif[:-2] + motif[-1] + " motif activation"].astype(float).values[500 + lo:500 + hi + 1]
        rows.append(dict(gene_name=g.gene_name, dominant_motif=motif,
                         dominant_contribution=g.dominant_contribution_excl_U1,
                         motif_center=int(lo + np.argmax(act)), motif_activation=float(act.max())))
    pd.DataFrame(rows).to_csv(OUT, sep="\t", index=False, float_format="%.4f")
    print(f"{len(rows)} genes -> {OUT}")


if __name__ == "__main__":
    main()
