"""Functions for haplotype- and allele-selective accessibility at LCL cCREs and cCRE pairs."""

import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
from scipy.stats import fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "CRE_co-accessibility_FIRE"))
from LCL_coaccess_cres import BEDTOOLS, CELLS, CHROMS, ROOT, incidence_for_chrom, sample_manifest
from count_and_test_pairs import bh, count_pairs

BCFTOOLS = "/project/spott/cshan/envs/bedtools/bin/bcftools"
PHASED_METATABLE = "/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_phased_27samples.csv"
AUTOSOMES = [c for c in CHROMS if c not in ("chrX", "chrY")]
FDR = 0.05
MIN_READS = 10  # per haplotype or allele, as in the paper's haplotype-selective peak test
CRAM_FIELDS = [b"required_fields=0x23F"]  # QNAME FLAG RNAME POS MAPQ CIGAR SEQ; skips QUAL and aux tags

#################################################
# FIRE haplotype-selective peaks -> cCRE status
#################################################

# <sample>-fire-v0.1-hap-differences.bed.gz: one row per peak that passed FIRE's per-haplotype
# coverage filter, with FIRE's two-sided Fisher test of hap1 vs hap2 accessibility
HAP_DIFF_COLS = ["#chrom", "peak_start", "peak_end", "coverage_H1", "coverage_H2",
                 "hap1_frac_acc", "hap2_frac_acc", "p_value", "p_adjust", "diff"]
REPRESENTATIVE = ["p_value", "p_adjust", "diff", "hap1_frac_acc", "hap2_frac_acc", "coverage_H1", "coverage_H2"]


def read_hap_differences(path):
    """FIRE peaks tested for haplotype selectivity; empty if the file is missing or empty (unphased)."""
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        return pd.DataFrame(columns=["chrom"] + HAP_DIFF_COLS[1:])
    peaks = pd.read_csv(path, sep="\t", usecols=HAP_DIFF_COLS)
    return peaks.rename(columns={"#chrom": "chrom"})


def bedtools_overlaps(a, b):
    """bedtools intersect -wa -wb of two BED-like frames (first 3 columns = chrom, start, end)."""
    with tempfile.TemporaryDirectory(prefix="hsca-") as tmp:
        a_bed, b_bed = Path(tmp) / "a.bed", Path(tmp) / "b.bed"
        a.to_csv(a_bed, sep="\t", header=False, index=False)
        b.to_csv(b_bed, sep="\t", header=False, index=False)
        out = subprocess.run([BEDTOOLS, "intersect", "-a", str(a_bed), "-b", str(b_bed), "-wa", "-wb"],
                             check=True, capture_output=True, text=True).stdout
    cols = [f"a_{c}" for c in a.columns] + [f"b_{c}" for c in b.columns]
    if not out:
        return pd.DataFrame(columns=cols)
    rows = [line.split("\t") for line in out.rstrip("\n").split("\n")]
    return pd.DataFrame(rows, columns=cols)


def cre_status(cres, peaks):
    """One row per cCRE overlapping >=1 haplotype-tested peak.

    status: conflict = significant peaks favour both haplotypes; hap1/hap2 = all significant peaks
    favour one haplotype; ns = none significant. The representative peak is the one with the
    smallest p_adjust (ties: smallest peak_start); its values are copied even for ns cCREs.
    """
    if peaks.empty:
        return pd.DataFrame()
    peaks = peaks.reset_index(drop=True)
    hits = bedtools_overlaps(cres[["chrom", "start", "end", "CRE_ID"]],
                             peaks[["chrom", "peak_start", "peak_end"]].assign(peak=peaks.index))
    if hits.empty:
        return pd.DataFrame()
    hits = pd.DataFrame({"CRE_ID": hits.a_CRE_ID, "peak": hits.b_peak.astype(int)})
    hits = hits.join(peaks, on="peak")
    significant = hits.p_adjust <= FDR
    hits["hap1_hit"] = significant & (hits["diff"] > 0)
    hits["hap2_hit"] = significant & (hits["diff"] < 0)
    grouped = hits.groupby("CRE_ID", sort=False)
    out = pd.DataFrame({"n_hap_peaks": grouped.size(),
                        "n_hsca_peaks": grouped.p_adjust.apply(lambda p: int((p <= FDR).sum())),
                        "hap1": grouped.hap1_hit.any(), "hap2": grouped.hap2_hit.any()})
    out["status"] = np.select([out.hap1 & out.hap2, out.hap1, out.hap2], ["conflict", "hap1", "hap2"], "ns")
    representative = hits.sort_values(["CRE_ID", "p_adjust", "peak_start"]).drop_duplicates("CRE_ID").set_index("CRE_ID")
    out = out.join(representative[REPRESENTATIVE]).reset_index()
    return out[["CRE_ID", "status", "n_hap_peaks", "n_hsca_peaks"] + REPRESENTATIVE]

#################################################
# heterozygous SNVs inside cCREs and the variant used per cCRE / pair
#################################################

def het_snvs(vcf, cell_line, cres):
    """PASS heterozygous biallelic SNVs of one sample that fall inside a cCRE.

    The VCFs hold all 120 YRI samples, and bcftools evaluates -i before -s, so the sample is
    subset first and then filtered (the paper's filter: -f PASS -i 'GT="HET"').
    """
    with tempfile.TemporaryDirectory(prefix="hsca-") as tmp:
        regions = Path(tmp) / "cres.bed"  # .bed suffix: bcftools reads it as 0-based
        cres[["chrom", "start", "end"]].to_csv(regions, sep="\t", header=False, index=False)
        subset = subprocess.Popen([BCFTOOLS, "view", "-R", str(regions), "-s", cell_line, "-f", "PASS", "-Ou", vcf],
                                  stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        het = subprocess.Popen([BCFTOOLS, "view", "-i", 'GT="HET"', "-v", "snps", "-m2", "-M2", "-Ou"],
                               stdin=subset.stdout, stdout=subprocess.PIPE)
        subset.stdout.close()
        query = subprocess.run([BCFTOOLS, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n"],
                               stdin=het.stdout, capture_output=True, text=True, check=True)
        het.stdout.close()
        if het.wait() or subset.wait():
            raise RuntimeError(f"bcftools failed on {vcf}")
    rows = [line.split("\t") for line in query.stdout.splitlines()]
    snvs = pd.DataFrame(rows, columns=["chrom", "pos", "ref", "alt"]).drop_duplicates()
    if snvs.empty:
        return pd.DataFrame(columns=["CRE_ID", "variant_id", "chrom", "pos", "ref", "alt"])
    snvs["pos"] = snvs.pos.astype(int)
    snvs["variant_id"] = snvs.chrom + ":" + snvs.pos.astype(str) + ":" + snvs.ref + ":" + snvs.alt
    hits = bedtools_overlaps(snvs.assign(start=snvs.pos - 1)[["chrom", "start", "pos", "variant_id"]],
                             cres[["chrom", "start", "end", "CRE_ID"]])
    hits = pd.DataFrame({"variant_id": hits.a_variant_id, "CRE_ID": hits.b_CRE_ID})
    return hits.merge(snvs, on="variant_id")[["CRE_ID", "variant_id", "chrom", "pos", "ref", "alt"]]


def choose_cre_variants(het, cres):
    """Per cCRE, the het SNV shared by the most samples; ties go to the SNV nearest the cCRE centre."""
    counts = het.groupby(["CRE_ID", "variant_id", "chrom", "pos", "ref", "alt"], sort=False).sample_name.agg(
        n_het_samples="nunique", het_samples=lambda s: ";".join(sorted(set(s)))).reset_index()
    centre = cres.set_index("CRE_ID")[["start", "end"]].mean(axis=1)
    counts["dist_to_centre"] = (counts.pos - 0.5 - counts.CRE_ID.map(centre)).abs()
    counts = counts.sort_values(["CRE_ID", "n_het_samples", "dist_to_centre", "variant_id"],
                                ascending=[True, False, True, True])
    return counts.drop_duplicates("CRE_ID").reset_index(drop=True)


def choose_pair_variants(pairs, cre_variants):
    """Per pair, the same rule over the het SNVs of CRE1 and CRE2: the better of the two cCRE
    choices (most het samples, then nearest its cCRE centre, then CRE1)."""
    choice = cre_variants.set_index("CRE_ID")
    sides = []
    for side in ("CRE1", "CRE2"):
        s = pairs[["CRE_pair", "CRE1", "CRE2", "chrom"]].join(choice.drop(columns="chrom"), on=side, how="inner")
        sides.append(s.assign(variant_cre=side))
    both = pd.concat(sides).sort_values(["CRE_pair", "n_het_samples", "dist_to_centre", "variant_cre"],
                                        ascending=[True, False, True, True])
    return both.drop_duplicates("CRE_pair").reset_index(drop=True)

#################################################
# per-read haplotype tags and alleles
#################################################

def load_haplotags(phasing_dir):
    """HiPhase read-level-phasing.tsv -> {chrom: {read_name: (phase_block_id, haplotag)}}.

    Keyed by (chrom, read_name) because chimeric reads appear once per chromosome; a read with
    conflicting tags on one chromosome counts as unphased.
    """
    tags = pd.read_csv(Path(phasing_dir) / "read-level-phasing.tsv", sep="\t",
                       usecols=["chrom", "phase_block_id", "read_name", "haplotag"])
    tags = tags[tags.haplotag.isin([1, 2])].drop_duplicates()
    tags = tags[~tags.duplicated(["chrom", "read_name"], keep=False)]
    return {chrom: dict(zip(g.read_name, zip(g.phase_block_id, g.haplotag))) for chrom, g in tags.groupby("chrom")}


def read_alleles(cram, chrom, variants):
    """{variant_id: {read_name: "ref" | "alt"}} from the base each primary read carries.

    Reads with a deletion or reference skip at the SNV, or a base matching neither allele, are dropped.
    """
    calls = {}
    with pysam.AlignmentFile(cram, "rc", format_options=CRAM_FIELDS) as aln:
        for v in variants.itertuples(index=False):
            reads = {}
            for column in aln.pileup(chrom, v.pos - 1, v.pos, truncate=True, stepper="nofilter",
                                     min_base_quality=0, max_depth=1_000_000, multiple_iterators=False):
                for read in column.pileups:
                    a = read.alignment
                    if a.is_secondary or a.is_supplementary or read.is_del or read.is_refskip:
                        continue
                    base = a.query_sequence[read.query_position]
                    if base == v.ref:
                        reads[a.query_name] = "ref"
                    elif base == v.alt:
                        reads[a.query_name] = "alt"
            calls[v.variant_id] = reads
    return calls

#################################################
# counts
#################################################

def block_local_hap_counts(pairs, coverage, accessible, tags):
    """Both-open counts per haplotype on the shared reads of each pair's majority phase block.

    H1/H2 are only consistent within one HiPhase block, so phased reads from other blocks are
    dropped (n_other_block_reads). hap*_acc = both cCREs open, hap*_nacc = not both open.
    """
    rows, empty = [], set()
    for pair in pairs.itertuples(index=False):
        shared = coverage.get(pair.CRE1, empty) & coverage.get(pair.CRE2, empty)
        phased = [(r, tags[r]) for r in shared if r in tags]
        if not phased:
            continue
        blocks = Counter(ps for _, (ps, _) in phased)
        block = min(blocks, key=lambda b: (-blocks[b], b))
        both = accessible.get(pair.CRE1, empty) & accessible.get(pair.CRE2, empty)
        n, acc = [0, 0], [0, 0]
        for r, (ps, hp) in phased:
            if ps == block:
                n[hp - 1] += 1
                acc[hp - 1] += r in both
        rows.append((pair.CRE_pair, block, len(phased) - n[0] - n[1], acc[0], n[0] - acc[0], acc[1], n[1] - acc[1]))
    return pd.DataFrame(rows, columns=["CRE_pair", "phase_block", "n_other_block_reads",
                                       "hap1_acc", "hap1_nacc", "hap2_acc", "hap2_nacc"])


def allele_counts(units, calls):
    """ref_open, ref_closed, alt_open, alt_closed per unit.

    units: (unit_id, variant_id, reads, open_reads). For a cCRE, reads = reads covering it and
    open = FIRE-accessible; for a pair, reads = reads spanning both and open = both accessible.
    Reads cover every base of the cCRE that carries the SNV, so each one is genotyped there.
    """
    rows = []
    for unit, variant, reads, open_reads in units:
        alleles = calls[variant]
        counts = [0, 0, 0, 0]
        for r in reads:
            allele = alleles.get(r)
            if allele is not None:
                counts[2 * (allele == "alt") + (r not in open_reads)] += 1
        rows.append((unit, variant, *counts))
    return pd.DataFrame(rows, columns=["unit", "variant_id", "ref_open", "ref_closed", "alt_open", "alt_closed"])

#################################################
# phase blocks
#################################################

def load_blocks(phasing_dir):
    """HiPhase blocks.tsv: 1-based first/last phased variant of each block."""
    blocks = pd.read_csv(Path(phasing_dir) / "blocks.tsv", sep="\t", usecols=["phase_block_id", "chrom", "start", "end"])
    return {chrom: g.sort_values("start") for chrom, g in blocks.drop_duplicates().groupby("chrom")}


def same_block(chrom, span_start, span_end, blocks):
    """True where the 0-based span [span_start, span_end) lies inside one HiPhase block."""
    g = blocks.get(chrom)
    if g is None:
        return np.zeros(len(span_start), dtype=bool)
    starts, ends = g.start.to_numpy(), g.end.to_numpy()
    i = np.searchsorted(starts, np.asarray(span_start) + 1, side="right") - 1
    return (i >= 0) & (ends[np.clip(i, 0, None)] >= np.asarray(span_end))

#################################################
# tests
#################################################

def fisher_bh(table, cols, group=None, testable=None):
    """Two-sided Fisher on [[c0, c1], [c2, c3]] plus BH within each group (as FIRE hap-diffs.R).

    Identical tables are tested once. Rows outside `testable` get NaN and are left out of BH.
    """
    counts = table[cols].to_numpy(dtype=np.int64)
    unique, inverse = np.unique(counts, axis=0, return_inverse=True)
    p = np.array([fisher_exact([[a, b], [c, d]], alternative="two-sided").pvalue for a, b, c, d in unique])
    table = table.assign(p_value=p[inverse.ravel()] if len(unique) else np.empty(0))
    if testable is not None:
        table.loc[~np.asarray(testable), "p_value"] = np.nan

    def adjust(s):
        out = pd.Series(np.nan, index=s.index)
        ok = s.notna()
        out[ok] = bh(s[ok])
        return out
    table["p_adjust"] = table.p_value.groupby(table[group] if group else np.zeros(len(table))).transform(adjust)
    return table


def call_status(diff, p_adjust, positive, negative):
    """positive/negative label by the sign of diff when p_adjust <= FDR, otherwise ns."""
    return np.where(p_adjust <= FDR, np.where(diff > 0, positive, negative), "ns")

#################################################
# samples
#################################################

def phased_samples(sample_metatable, root):
    """The 31-sample manifest (paths, donor, phasing_dir, phased_vcf_file) restricted to the phased sheet."""
    manifest = sample_manifest(sample_metatable, root)
    phased = pd.read_csv(PHASED_METATABLE, dtype=str, keep_default_na=False)
    samples = manifest[manifest.sample_name.isin(phased.sample_name)].reset_index(drop=True)
    check = samples.merge(phased, on="sample_name", suffixes=("", "_phased"))
    if len(samples) != len(phased) or any(check[c].ne(check[c + "_phased"]).any()
                                          for c in ("cell_line", "fire_dir", "phasing_dir", "phased_vcf_file")):
        raise ValueError("The phased sheet must match the 31-sample manifest")
    return samples
