"""Promoter-only adapter for 03_pol2_cluster_bams.py: fetch each window once.
Preserves its selection, tie ordering, BAM-header handling, and MM/ML tags.
The original helper and all shared functions remain unchanged.
"""
import argparse
import csv
from collections import defaultdict
from pathlib import Path
import pysam


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    with open(args.manifest) as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if len({r['RID'] for r in rows}) != len(rows):
        raise ValueError('Clustered physical molecules must be unique')
    by_bam = defaultdict(list)
    for row in rows:
        by_bam[row['merged_bam_file']].append(row)
    out = Path(args.output)
    with pysam.AlignmentFile(rows[0]['merged_bam_file'], 'rb') as first:
        header = first.header.to_dict()
    source_sq = header['SQ']
    # pb-CpG-tools scans every header contig, including unused decoys absent
    # from hg38.fa. Keep the first sample's header with only requested contigs.
    chroms = {r['chrom'] for r in rows}
    header['SQ'] = [sq for sq in source_sq if sq['SN'] in chroms]
    target_ids = {sq['SN']: i for i, sq in enumerate(header['SQ'])}
    if set(target_ids) != chroms:
        raise ValueError('Requested chromosome absent from BAM header')
    writers = {cl: pysam.AlignmentFile(str(out / (cl + '.unsorted.bam')), 'wb', header=header)
               for cl in sorted({r['cluster'] for r in rows})}
    counts = defaultdict(int)
    try:
        for path, molecules in by_bam.items():
            with pysam.AlignmentFile(path, 'rb') as bam:
                if bam.header.to_dict()['SQ'] != source_sq:
                    raise ValueError('BAM reference dictionaries differ: ' + path)
                windows = defaultdict(list)
                for row in molecules:
                    windows[(row['chrom'], int(row['window_start']), int(row['window_end']))].append(row)
                for (chrom, lo, hi), members in sorted(windows.items(), key=lambda x: (target_ids[x[0][0]], x[0][1])):
                    wanted = {row['original_RID'] for row in members}
                    best = {}
                    for read in bam.fetch(chrom, lo, hi):
                        if read.is_unmapped or read.query_name not in wanted:
                            continue
                        rank = (read.reference_length or 0, not read.is_secondary, not read.is_supplementary)
                        old = best.get(read.query_name)
                        if old is None or rank > old[0]:
                            best[read.query_name] = (rank, read)
                    for row in members:
                        if row['original_RID'] not in best:
                            raise ValueError('Missing alignment: ' + row['RID'])
                        # Copy because the same read may be selected for more than one window.
                        read = pysam.AlignedSegment.fromstring(best[row['original_RID']][1].to_string(), bam.header)
                        ref_name = read.reference_name
                        mate_name = read.next_reference_name
                        read.reference_id = target_ids[ref_name]
                        read.next_reference_id = target_ids.get(mate_name, -1)
                        writers[row['cluster']].write(read)
                        counts[row['cluster']] += 1
            print('Extracted', len(molecules), 'molecules from', Path(path).name, flush=True)
    finally:
        for writer in writers.values():
            writer.close()
    for cl in writers:
        raw = out / (cl + '.unsorted.bam')
        target = out / (cl + '.bam')
        pysam.sort('-o', str(target), str(raw))
        pysam.index(str(target))
        with pysam.AlignmentFile(str(target), 'rb') as bam:
            assert bam.count(until_eof=True) == sum(r['cluster'] == cl for r in rows)
        raw.unlink()
        print(cl, counts[cl], flush=True)


if __name__ == '__main__':
    main()
