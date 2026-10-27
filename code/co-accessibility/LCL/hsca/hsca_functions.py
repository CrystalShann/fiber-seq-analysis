"""Functions for allele-selective accessibility at LCL cCREs and cCRE pairs."""

import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
from scipy.stats import fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "CRE_co-accessibility_FIRE"))
from LCL_coaccess_cres import BEDTOOLS, CELLS, CHROMS, ROOT, incidence_for_chrom
from count_and_test_pairs import bh, count_pairs

BCFTOOLS = "/project/spott/cshan/envs/bedtools/bin/bcftools"
PHASED_METATABLE = "/project/spott/1_Shared_projects/LCL_Fiber_seq/Data/LCL_sample_metatable_merged_phased_27samples.csv"
AUTOSOMES = [c for c in CHROMS if c not in ("chrX", "chrY")]
FDR = 0.05
MIN_READS = 10  # REF and ALT reads per donor/unit, before testing or pooling
CRAM_FIELDS = [b"required_fields=0x23F"]  # QNAME FLAG RNAME POS MAPQ CIGAR SEQ; skips QUAL and aux tags

#################################################
# find overlaps in 2 bed files
#################################################
# Find every overlapping pair of intervals in two DataFrames using BED coordinates.
# The first three columns must be chromosome, 0-based start, and exclusive end.
# Return both rows with a_/b_ column prefixes; values parsed from bedtools are strings.
# Temporary BED files are deleted automatically, and a failed bedtools command raises an error.
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

#################################################
# heterozygous SNVs inside cCREs and the variant used per cCRE / pair
#################################################

# find the heterozygous SNVs of one sample that sit inside cCRE -> read genotype from VCF with bcftools, 
# then uses bedtools to link each SNV to its CRE

# Extract this donor’s PASS, heterozygous, biallelic SNVs from the multi-sample VCF.
# Subset the donor before testing heterozygosity so another donor’s genotype cannot qualify a site.
# Intersect the SNVs with cCREs and return one row per variant/cCRE overlap.

def het_snvs(vcf, cell_line, cres):
    """PASS heterozygous biallelic SNVs of one sample that fall inside a cCRE.

    The VCFs hold all 120 YRI samples, and bcftools evaluates -i before -s, so the sample is
    subset first and then filtered (the paper's filter: -f PASS -i 'GT="HET"').
    """
    with tempfile.TemporaryDirectory(prefix="hsca-") as tmp:
        regions = Path(tmp) / "cres.bed"  # .bed suffix: bcftools reads it as 0-based
        cres[["chrom", "start", "end"]].to_csv(regions, sep="\t", header=False, index=False)
        # subset = keep regions inside the CRE, and the sample's genotype column, keep only sites has pass filter
        subset = subprocess.Popen([BCFTOOLS, "view", "-R", str(regions), "-s", cell_line, "-f", "PASS", "-Ou", vcf],
                                  stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        # het (from subset piple) = keep sites where the remaining sample is heterozygous, keep SNPs and ref vs alt allele
        het = subprocess.Popen([BCFTOOLS, "view", "-i", 'GT="HET"', "-v", "snps", "-m2", "-M2", "-Ou"],
                               stdin=subset.stdout, stdout=subprocess.PIPE)
        subset.stdout.close()
        # print chrom, pos, ref, alt, one line per site
        query = subprocess.run([BCFTOOLS, "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\n"],
                               stdin=het.stdout, capture_output=True, text=True, check=True)
        het.stdout.close()
        if het.wait() or subset.wait():
            raise RuntimeError(f"bcftools failed on {vcf}")
    rows = [line.split("\t") for line in query.stdout.splitlines()]
    # Converts pos to an integer and builds variant_id as chrom:pos:ref:alt
    snvs = pd.DataFrame(rows, columns=["chrom", "pos", "ref", "alt"]).drop_duplicates()
    if snvs.empty:
        return pd.DataFrame(columns=["CRE_ID", "variant_id", "chrom", "pos", "ref", "alt"])
      
    snvs["pos"] = snvs.pos.astype(int)
    snvs["variant_id"] = snvs.chrom + ":" + snvs.pos.astype(str) + ":" + snvs.ref + ":" + snvs.alt
    hits = bedtools_overlaps(snvs.assign(start=snvs.pos - 1)[["chrom", "start", "pos", "variant_id"]],
                             cres[["chrom", "start", "end", "CRE_ID"]])
    hits = pd.DataFrame({"variant_id": hits.a_variant_id, "CRE_ID": hits.b_CRE_ID})
    return hits.merge(snvs, on="variant_id")[["CRE_ID", "variant_id", "chrom", "pos", "ref", "alt"]]

#################################################
# picks one anchor SNV per cCRE, the het SNV shared by the most samples
#################################################

# Choose one common anchor SNP per cCRE for REF/ALT comparisons across samples.
# Rank by number of heterozygous samples, then distance to the cCRE centre, then variant ID.
# Retain the supporting sample count and semicolon-separated sample names for later filtering.
# cCREs with no candidate heterozygous SNP do not appear in the result.
def choose_cre_variants(het, cres):
    """Per cCRE, the het SNV shared by the most samples; ties go to the SNV nearest the cCRE centre."""
    counts = het.groupby(["CRE_ID", "variant_id", "chrom", "pos", "ref", "alt"], sort=False).sample_name.agg(
        n_het_samples="nunique", het_samples=lambda s: ";".join(sorted(set(s)))).reset_index()
    centre = cres.set_index("CRE_ID")[["start", "end"]].mean(axis=1)
    counts["dist_to_centre"] = (counts.pos - 0.5 - counts.CRE_ID.map(centre)).abs()
    counts = counts.sort_values(["CRE_ID", "n_het_samples", "dist_to_centre", "variant_id"],
                                ascending=[True, False, True, True])
    return counts.drop_duplicates("CRE_ID").reset_index(drop=True)

#################################################
# picks one anchor SNV for each cCRE pair, taken from whichever of the two cCREs has the better anchor
#################################################

# Choose one anchor SNP for each cCRE pair from the selected SNPs of its two members.
# Prefer the SNP heterozygous in more samples, then closer to its own cCRE centre;
# a remaining tie selects CRE1. variant_cre records which member supplied the SNP.
# A pair can qualify with a SNP in just one member; pairs with neither are omitted.

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
# records which allele (REF or ALT) each read carries at each anchor SNP
#################################################

# Read the actual base at each selected SNP from primary CRAM alignments.
# Return a nested variant-ID/read-name lookup of ref or alt; no HiPhase tag is required.
# Convert 1-based SNP positions to pileup coordinates, exclude deletions/reference skips
# Labels it "ref" if it equals the VCF REF base, "alt" if it equals the ALT base, and drops the read if it matches neither

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

######################################
# counts, for each unit, how many REF and ALT reads are open vs closed
######################################


# CRE level:	One cCRE +	 CRE_ID	+ reads covering that cCRE	reads on which it is FIRE-accessible
# Pair level:	One cCRE pair	+ CRE_pair +	reads spanning both cCREs	reads on which both cCREs are accessible

# Build the REF/ALT by open/closed 2x2 counts for each supplied cCRE or cCRE pair


# Output:
# unit	CRE_ID or CRE_pair.
# variant_id	The anchor SNV used to split the reads.
# ref_open, ref_closed	REF-carrying reads that are open / not open
# alt_open, alt_closed	ALT-carrying reads that are open / not open

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
# Fisher test + BH
#################################################

# 1. REF vs ALT within each donor, for a cCRE or a pair: allele_counts() supplies the cells.
# 2. REF vs ALT pooled across donors: pooled_alleles() sums the eligible donor cells.
# Both analyses require >= MIN_READS REF and ALT reads per donor/unit before testing or pooling.

def fisher_bh(table, cols, group=None, testable=None):
    """Two-sided Fisher on [[c0, c1], [c2, c3]] plus BH within each group
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

#################################################
# lable fisher test result
#################################################

# Positive means group 1 is more open
# not significant results get ns

def call_status(diff, p_adjust, positive, negative):
    """positive/negative label by the sign of diff when p_adjust <= FDR, otherwise ns."""
    return np.where(p_adjust <= FDR, np.where(diff > 0, positive, negative), "ns")

#################################################
# build the manifest directly from the 27 phased samples
#################################################

# 1. Read the phased sheet, requiring exactly 27 distinct, nonempty sample identifiers
# 2. Reject another cohort or mismatched donor/phasing paths instead of silently dropping samples
# 3. Derive the FIRE and upstream per-sample paths used by the counting 

def phased_samples(sample_metatable=PHASED_METATABLE, root=ROOT):
    """Validated 27-sample HSCA manifest; independent of the 31-sample manifest loader."""
    samples = pd.read_csv(sample_metatable, dtype=str, keep_default_na=False)
    required = ["sample_name", "cell_line", "fire_dir", "phasing_dir", "phased_vcf_file"]
    if not set(required).issubset(samples.columns):
        raise ValueError(f"Phased sample CSV requires {required}")
    if len(samples) != 27 or samples.sample_name.nunique() != 27 or samples[required].eq("").any().any():
        raise ValueError("HSCA requires exactly 27 unique phased samples with complete paths")
    canonical = pd.read_csv(PHASED_METATABLE, dtype=str, keep_default_na=False)
    if set(samples.sample_name) != set(canonical.sample_name):
        raise ValueError("HSCA sample names must match the canonical 27 phased samples")
    expected = canonical.set_index("sample_name").loc[samples.sample_name, required[1:]]
    if not np.array_equal(samples[required[1:]].to_numpy(), expected.to_numpy()):
        raise ValueError("HSCA donor and phasing paths must match the phased sheet")
    root = Path(root)
    samples["donor"] = samples.cell_line
    samples["cram_path"] = [str(Path(r.fire_dir) / f"{r.sample_name}-fire-v0.1-filtered.cram") for r in samples.itertuples()]
    samples["peaks_path"] = [str(Path(r.fire_dir) / f"{r.sample_name}-fire-v0.1-peaks.bed.gz") for r in samples.itertuples()]
    samples["elements_path"] = [str(Path(r.fire_dir) / "additional-outputs-v0.1/fire-peaks" / f"{r.sample_name}-v0.1-fire-elements.bed.gz") for r in samples.itertuples()]
    samples["spans_path"] = [str(root / s / f"{s}.read_spans.bed.gz") for s in samples.sample_name]
    samples["blocks_path"] = [str(root / s / f"{s}.aligned_blocks.bed.gz") for s in samples.sample_name]
    return samples
