#!/usr/bin/env python
# Puffin motif contributions at the fiber TSS of every curated PolII example gene.



# Run: /project/spott/cshan/envs/puffin/bin/python \
#        code/pausing/PolII/promoter_annotation/04a_puffin_motif_scores.py
#
# The summed motif / initiator / trinucleotide effects get the same TSS weighting, and their
# per-position values over TSS -500..+499 are kept for the effect profiles. tata_center is the
# TSS-relative position of the peak TATA+ motif activation within TATA_SEARCH.
#
# Outputs under LCL_project/pausing/PolII/promoter_annotation/tables/
#   04a_puffin_motif_scores.tsv        one row per gene
#   04a_puffin_effect_profiles.tsv.gz  one row per gene x position (-500..499)



# This script runs Puffin for all curated examples genes, and turns per base predictions
# into one score per motif at the TSS, and lables each gene by its strongest motif


import os
import sys

import numpy as np
import pandas as pd
import pysam
import torch

PROJECT = "/project/spott/cshan/fiber-seq"
GENES = f"{PROJECT}/LCL_project/pausing/PolII/read_features/tables/gene_features.tsv"
OUT = f"{PROJECT}/LCL_project/pausing/PolII/promoter_annotation/tables/04a_puffin_motif_scores.tsv"
PROFILE_OUT = f"{PROJECT}/LCL_project/pausing/PolII/promoter_annotation/tables/04a_puffin_effect_profiles.tsv.gz"
GENOME_FA = "/project/spott/cshan/annotations/hg38.fa"
PUFFIN_SRC = "/project/spott/cshan/envs/puffin/puffin_src"

MOTIFS = ["TATA +", "YY1 +", "SP +", "SP -", "ETS +", "ETS -", "NFY +", "NFY -", "CREB +", "CREB -",
          "NRF1 +", "NRF1 -", "ZNF143 +", "ZNF143 -", "U1 snRNP +"]
CLASS_MOTIFS = MOTIFS[:-1]
PAD = 825      # 500 bp output each side + Puffin's 325 bp padding
WINDOW = slice(490, 510)  # TSS -10..+9 of the 1000 bp output, TSS at index 500
# Puffin rows of the summed effects (initiator and trinucleotide are centred over the 1 kb output)
TATA_SEARCH = (-60, -10)  # TATA box = peak TATA+ motif activation (motif centre) in this TSS-relative range
SUMS = {"total": "Sum of total effect", "motif": "Sum of motif effect",
        "initiator": "Sum of initiator effect", "trinucleotide": "Sum of trinucleotide effect"}

sys.path.insert(0, PUFFIN_SRC)
os.chdir(PUFFIN_SRC)  # Puffin loads ./resources/puffin.pth
from puffin import Puffin  


def tss_sequence(fa, chrom, tss0, strand):
    # same window as the manuscript: [TSS - 825, TSS + 825) shifted by 1 on the minus strand
    off = 1 if strand == "-" else 0
    s = fa.fetch(chrom, tss0 - PAD + off, tss0 + PAD + off).upper()
    if strand == "-":
        s = s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))
    return s


def tss_effects(net, seq):
    # returns the TSS-weighted scores and the per-position summed effects
    df = net.interpret(seq, "FANTOM_CAGE")
    w = df.loc["Prediction"].astype(float).values[WINDOW]
    row = lambda name: df.loc[name].astype(float).values
    scores = {m: float((w * row(m[:-2] + m[-1] + " motif effect")[WINDOW]).sum() / w.sum()) for m in MOTIFS}
    scores.update({f"{k}_effect_tss": float((w * row(v)[WINDOW]).sum() / w.sum()) for k, v in SUMS.items()})
    act = row("TATA+ motif activation")[500 + TATA_SEARCH[0]:500 + TATA_SEARCH[1] + 1]
    scores["tata_center"] = int(TATA_SEARCH[0] + np.argmax(act))
    scores["tata_activation"] = float(act.max())
    profile = pd.DataFrame({"position": np.arange(-500, 500), **{k: row(v) for k, v in SUMS.items()}})
    return scores, profile


def main():
    torch.set_num_threads(int(os.environ.get("SLURM_CPUS_PER_TASK", 4)))
    genes = pd.read_csv(GENES, sep="\t")
    genes = genes[genes.tss.notna()]
    fa = pysam.FastaFile(GENOME_FA)
    net = Puffin()
    rows, profiles = [], []
    for g in genes.itertuples():
        tss0 = int(g.tss) - 1  # fiber tss is 1-based, Puffin's is 0-based
        scores, profile = tss_effects(net, tss_sequence(fa, g.chrom, tss0, g.strand))
        rows.append(dict(gene_name=g.gene_name, chrom=g.chrom, puffin_tss0=tss0, strand=g.strand, **scores))
        profiles.append(profile.assign(gene_name=g.gene_name, group=g.group))
    res = pd.DataFrame(rows)
    dom = res[CLASS_MOTIFS].idxmax(axis=1)
    res["dominant_motif_excl_U1"] = dom
    res["dominant_contribution_excl_U1"] = res[CLASS_MOTIFS].max(axis=1)
    res["group_excl_U1"] = dom.str.split(" ").str[0].replace({"SP": "SP1"})
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    res.to_csv(OUT, sep="\t", index=False, float_format="%.4f")
    pd.concat(profiles)[["gene_name", "group", "position", *SUMS]].to_csv(
        PROFILE_OUT, sep="\t", index=False, float_format="%.4f")
    print(f"{len(res)} genes -> {OUT}, {PROFILE_OUT}")


if __name__ == "__main__":
    main()
