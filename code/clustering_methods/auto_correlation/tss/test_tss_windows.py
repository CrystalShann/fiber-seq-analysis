#!/usr/bin/env python3
"""Hand-made tests of the strand-oriented window logic in 01_sample_tss_molecules.py."""
import importlib.util
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
spec = importlib.util.spec_from_file_location("tss_sampler", HERE / "01_sample_tss_molecules.py")
sampler = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sampler)


def test_minus_strand_genomic_window():
    ws, we = sampler.genomic_window(["-", "+"], [10000, 10000], -1000, 1000)
    assert ws.tolist() == [9001, 9000] and we.tolist() == [11001, 11000]
    ws, we = sampler.genomic_window(["-", "+"], [10000, 10000], -1000, -100)
    assert ws.tolist() == [10101, 9000] and we.tolist() == [11001, 9900]
    ws, we = sampler.genomic_window(["-", "+"], [10000, 10000], 100, 1000)
    assert ws.tolist() == [9001, 10100] and we.tolist() == [9901, 11000]


def test_minus_strand_read_edges():
    """A minus-strand gene at tss 10000, window [-1000, +1000) -> genomic [9001, 11001)."""
    ws, we = sampler.genomic_window(["-"], [10000], -1000, 1000)
    width = 2000
    accepted = [(9001, 11001), (9000, 11001), (9001, 11002), (8000, 12000)]
    rejected = [(9002, 11001), (9001, 11000), (9002, 11000), (11001, 13000), (7000, 9001)]
    for st, en in accepted:
        lo, hi = sampler.eligible_gene_range(ws, width, st, en)
        assert hi - lo == 1, f"read [{st}, {en}) should be accepted"
    for st, en in rejected:
        lo, hi = sampler.eligible_gene_range(ws, width, st, en)
        assert hi <= lo, f"read [{st}, {en}) should be rejected"
    # One-sided upstream window on the minus strand: [-1000, -100) -> genomic [10101, 11001).
    ws, we = sampler.genomic_window(["-"], [10000], -1000, -100)
    for st, en, ok in [(10101, 11001, True), (10102, 11001, False), (10101, 11000, False), (9000, 12000, True)]:
        lo, hi = sampler.eligible_gene_range(ws, 900, st, en)
        assert (hi - lo == 1) == ok, f"upstream read [{st}, {en})"


def test_oriented_columns_and_reversal():
    tss, s, e = 10000, -1000, 1000
    for strand in ("+", "-"):
        ws, we = sampler.genomic_window([strand], [tss], s, e)
        column = sampler.oriented_column(strand, tss, s, tss)
        assert column == -s == 1000
        assert sampler.column_genomic_position(strand, ws[0], we[0], column) == tss
        # Build a genomic-order row with one call 3 bp downstream of the TSS and reverse for '-'.
        genomic = tss + 3 if strand == "+" else tss - 3
        row = np.zeros(e - s, dtype=np.uint8)
        row[genomic - ws[0]] = 1
        if strand == "-":
            row = row[::-1]
        assert row[1003] == 1 and row.sum() == 1
        assert sampler.oriented_column(strand, tss, s, genomic) == 1003


if __name__ == "__main__":
    for name, test in sorted(globals().items()):
        if name.startswith("test_") and callable(test):
            test()
            print(f"PASS {name}", flush=True)
    print("All tss window tests passed")
