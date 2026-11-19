"""Subset CpG inputs to each gene's chromosome without changing sequence or MM/ML tags."""
import argparse
import csv
from pathlib import Path
import pysam


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--manifest', required=True)
    p.add_argument('--reference', required=True)
    p.add_argument('--output', required=True)
    args = p.parse_args()
    rows = list(csv.DictReader(open(args.manifest), delimiter='\t'))
    genes = {r['cluster']: r['chrom'] for r in rows}
    root = Path(args.output)
    refs = root / 'references'
    refs.mkdir(exist_ok=True)
    with pysam.FastaFile(args.reference) as source:
        for chrom in sorted(set(genes.values())):
            target = refs / (chrom + '.fa')
            if Path(str(target) + '.fai').exists():
                continue
            with open(target, 'w') as out:
                out.write('>' + chrom + '\n')
                for lo in range(0, source.get_reference_length(chrom), 1000000):
                    sequence = source.fetch(chrom, lo, min(lo + 1000000, source.get_reference_length(chrom)))
                    for j in range(0, len(sequence), 100):
                        out.write(sequence[j:j + 100] + '\n')
            pysam.faidx(str(target))
    for gene, chrom in genes.items():
        target = root / (gene + '.chrom.bam')
        if Path(str(target) + '.bai').exists():
            continue
        with pysam.AlignmentFile(str(root / (gene + '.bam')), 'rb') as source:
            header = source.header.to_dict()
            header['SQ'] = [x for x in header['SQ'] if x['SN'] == chrom]
            assert len(header['SQ']) == 1
            with pysam.AlignmentFile(str(target), 'wb', header=header) as out:
                for read in source:
                    assert read.reference_name == chrom
                    mate = read.next_reference_name
                    read.reference_id = 0
                    read.next_reference_id = 0 if mate == chrom else -1
                    out.write(read)
        pysam.index(str(target))
    print('Prepared chromosome references and', len(genes), 'single-chromosome BAMs', flush=True)


if __name__ == '__main__':
    main()
