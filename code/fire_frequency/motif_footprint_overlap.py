"""Footprint-supported motif presence for FIRE comparisons."""
from pathlib import Path
from collections import defaultdict
import hashlib
import json
import numpy as np
import pysam


def build_footprint_index(project_root, chromosomes, cache_dir, timepoints=(0, 15)):
    """Pool selected timepoints, score>=50, 10..80 bp footprints within same-read FIREs.

    Coordinates are BED (zero-based, half-open). Keep individual intervals;
    merging adjacent footprints could falsely satisfy the motif overlap rule.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    fire_root = Path('/project/spott/lizarraga/pacbio_analysis/macrophage_project/merged_hifi_bams/FIRE')
    tf_root = Path(project_root) / 'macrophage_project/FiberHMM/extract/firehmm_tf'
    result = {}
    for chrom in sorted(chromosomes):
        sources = []
        for time in timepoints:
            sample = f'LPS_{time}'
            sources.append((
                fire_root / sample / f'additional-outputs-v0.1/fire-peaks/{sample}-v0.1-fire-elements.bed.gz',
                tf_root / f'LPS{time}' / f'LPS{time}_hmm_extracted_tf_{chrom}.bed.gz'))
        signature = hashlib.sha256(json.dumps({
            'version': 1, 'score': 50, 'size': [10, 80], 'same_read_fire': True,
            'sources': [(str(p), p.stat().st_size, p.stat().st_mtime_ns) for pair in sources for p in pair]
        }, sort_keys=True).encode()).hexdigest()
        cache = cache_dir / f'{chrom}.npz'
        metadata = cache_dir / f'{chrom}.json'
        if cache.exists() and metadata.exists() and json.loads(metadata.read_text()).get('signature') == signature:
            with np.load(cache) as data:
                result[chrom] = (data['start'], data['end'])
            print(f'{chrom}: reused {len(result[chrom][0]):,} eligible footprint intervals', flush=True)
            continue
        intervals = set()
        for fire_path, tf_path in sources:
            fires = defaultdict(list)
            with pysam.TabixFile(str(fire_path)) as stream:
                for line in stream.fetch(chrom):
                    fields = line.split('\t')
                    if fields[8] not in ('169,169,169', '147,112,219'):
                        fires[fields[3]].append((int(fields[1]), int(fields[2])))
            with pysam.TabixFile(str(tf_path)) as stream:
                for line in stream.fetch(chrom):
                    fields = line.split('\t')
                    segments = fires.get(fields[3])
                    if not segments:
                        continue
                    sizes = fields[10].rstrip(',').split(',')
                    offsets = fields[11].rstrip(',').split(',')
                    scores = fields[12].rstrip(',').split(',')
                    assert len(sizes) == len(offsets) == len(scores) == int(fields[9])
                    base = int(fields[1])
                    for size, offset, score in zip(sizes, offsets, scores):
                        size = int(size)
                        if not 10 <= size <= 80 or float(score) < 50:
                            continue
                        start = base + int(offset)
                        end = start + size
                        if any(a <= start and end <= b for a, b in segments):
                            intervals.add((start, end))
            del fires
        array = np.array(sorted(intervals), dtype=np.int64).reshape(-1, 2)
        result[chrom] = (array[:, 0], array[:, 1])
        np.savez_compressed(cache, start=array[:, 0], end=array[:, 1])
        metadata.write_text(json.dumps({'signature': signature, 'n_intervals': len(array)}))
        print(f'{chrom}: indexed {len(array):,} eligible footprint intervals', flush=True)
    return result


def supported_motifs(left, right, footprints, minimum_fraction=0.5):
    """True iff ONE footprint covers >=minimum_fraction of each motif.

    Equivalent to overlap_bp >= 0.5 * motif_length; no strand constraint.
    Prefix maxima implement interval existence queries without merging footprints.
    """
    left, right = np.asarray(left), np.asarray(right)
    assert np.all(right > left)
    starts, ends = footprints
    required = np.ceil((right - left) * minimum_fraction).astype(np.int64)
    keep = np.zeros(len(left), dtype=bool)
    for length in np.unique(required):
        fp_keep = ends - starts >= length
        eligible_starts = starts[fp_keep]
        if not len(eligible_starts):
            continue
        max_ends = np.maximum.accumulate(ends[fp_keep])
        rows = np.flatnonzero(required == length)
        idx = np.searchsorted(eligible_starts, right[rows] - length, side='right') - 1
        valid = idx >= 0
        keep[rows[valid]] = max_ends[idx[valid]] >= left[rows[valid]] + length
    return keep
