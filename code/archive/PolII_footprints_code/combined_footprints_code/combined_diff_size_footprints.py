import gzip
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

MERGED_DIR = Path('/project/spott/1_Shared_projects/LCL_Fiber_seq/FiberHMM/merged')
OUTPUT_ROOT = Path('/project/spott/cshan/fiber-seq/results/PolII/FIRE_combined_footprints')
COMBINED_PEAKS_DIR = OUTPUT_ROOT / 'joint_trained_peaks'
COMBINED_TRACKS_DIR = OUTPUT_ROOT / 'joint_trained_tracks'
CHROM_SIZES = Path('/project/spott/cshan/annotations/hg38.chrom.sizes')
BEDGRAPH_TO_BIGWIG = '/project/spott/cshan/tools/bedGraphToBigWig'

SIZE_BINS = {
    '80-100': (80, 100),
    '100-120': (100, 120),
    '100-140': (100, 140),
    '120-140': (120, 140),
}


def run(cmd, stdout=None):
    print('+', ' '.join(map(str, cmd)))
    subprocess.run(cmd, check=True, stdout=stdout)


def chrom_key(chrom):
    label = chrom.replace('chr', '')
    if label.isdigit():
        return (0, int(label))
    rank = {'X': 23, 'Y': 24, 'M': 25, 'MT': 25}
    return (1, rank.get(label, 999), label)


def parse_block_list(value):
    value = value.strip()
    if not value or value == '.':
        return []

    try:
        return [int(x) for x in value.rstrip(',').split(',') if x]
    except ValueError as exc:
        raise ValueError(f'Invalid block list {value!r}') from exc


def expand_fp_bed_line(fields):
    if len(fields) < 10:
        raise ValueError(f'Expected at least 10 columns, got {len(fields)}: {fields}')

    chrom = fields[0]
    chrom_start = int(fields[1])
    fiber_name = fields[3]
    block_count = int(fields[7])
    block_starts = parse_block_list(fields[8])
    block_sizes = parse_block_list(fields[9])

    if block_count == 0:
        if block_starts or block_sizes:
            raise ValueError(f'Expected no blocks for {fiber_name}, got {fields}')
        return

    if block_count != len(block_starts) or block_count != len(block_sizes):
        raise ValueError(f'Block count mismatch for {fiber_name}: {fields}')

    for rel_start, block_size in zip(block_starts, block_sizes):
        fp_start = chrom_start + rel_start
        fp_end = fp_start + block_size
        yield chrom, fp_start, fp_end, fiber_name, block_size


def sort_to_file(in_path, out_path):
    with open(out_path, 'w') as out_handle:
        run(['sort', '-k1,1', '-k2,2n', '-k3,3n', str(in_path)], stdout=out_handle)


def bgzip_file(in_path, out_path_gz):
    with open(out_path_gz, 'wb') as out_handle:
        run(['bgzip', '-c', str(in_path)], stdout=out_handle)


def make_bedgraph(sorted_bed_path, out_bedgraph_path):
    with open(out_bedgraph_path, 'w') as out_handle:
        run(
            ['bedtools', 'genomecov', '-bg', '-i', str(sorted_bed_path), '-g', str(CHROM_SIZES)],
            stdout=out_handle,
        )


def ensure_tools_exist():
    for tool_name in ['bgzip', 'tabix', 'bedtools']:
        if shutil.which(tool_name) is None:
            raise RuntimeError(f'Missing required tool on PATH: {tool_name}')
    if not CHROM_SIZES.exists():
        raise FileNotFoundError(f'Chrom sizes file not found: {CHROM_SIZES}')


ensure_tools_exist()
COMBINED_PEAKS_DIR.mkdir(parents=True, exist_ok=True)
COMBINED_TRACKS_DIR.mkdir(parents=True, exist_ok=True)

fp_files = sorted(MERGED_DIR.glob('AL*_bc*_*/joint_trained_results/*_chr*_fp.bed.gz'))
if not fp_files:
    raise FileNotFoundError(f'No sample footprint files found under {MERGED_DIR}')

files_by_chrom = defaultdict(list)
for fp_path in fp_files:
    match = re.search(r'_(chr[^_]+)_fp\.bed\.gz$', fp_path.name)
    if match:
        files_by_chrom[match.group(1)].append(fp_path)

chroms = sorted(files_by_chrom, key=chrom_key)
print(f'Found {len(fp_files)} footprint files across {len(chroms)} chromosomes')
print(f'Chromosomes: {chroms}')

with tempfile.TemporaryDirectory(prefix='combine_joint_trained_fp_bins_') as tmp_dir_name:
    tmp_dir = Path(tmp_dir_name)
    all_chrom_unsorted = {label: tmp_dir / f'combined_all_chrs_{label}bp_fps.unsorted.bed' for label in SIZE_BINS}
    all_chrom_handles = {label: open(path, 'w') for label, path in all_chrom_unsorted.items()}
    all_chrom_counts = defaultdict(int)

    try:
        for chrom in chroms:
            chrom_dir = COMBINED_TRACKS_DIR / chrom
            chrom_dir.mkdir(parents=True, exist_ok=True)

            chrom_unsorted = {
                label: tmp_dir / f'{chrom}_{label}bp_fps.unsorted.bed'
                for label in SIZE_BINS
            }
            chrom_handles = {label: open(path, 'w') for label, path in chrom_unsorted.items()}
            chrom_counts = defaultdict(int)

            try:
                for fp_path in files_by_chrom[chrom]:
                    print(f'Processing {fp_path}')
                    with gzip.open(fp_path, 'rt') as handle:
                        for line in handle:
                            line = line.strip()
                            if not line or line.startswith('#'):
                                continue
                            fields = line.split('\t')
                            for fp_chrom, fp_start, fp_end, fiber_name, fp_len in expand_fp_bed_line(fields):
                                for label, (min_len, max_len) in SIZE_BINS.items():
                                    if min_len <= fp_len <= max_len:
                                        record = '\t'.join([fp_chrom, str(fp_start), str(fp_end), fiber_name]) + '\n'
                                        chrom_handles[label].write(record)
                                        all_chrom_handles[label].write(record)
                                        chrom_counts[label] += 1
                                        all_chrom_counts[label] += 1
            finally:
                for handle in chrom_handles.values():
                    handle.close()

            for label in SIZE_BINS:
                count = chrom_counts[label]
                if count == 0:
                    print(f'Skipping {chrom} {label}: no intervals found')
                    continue

                sorted_bed = tmp_dir / f'{chrom}_{label}bp_fps.sorted.bed'
                final_bed_gz = chrom_dir / f'combined_{chrom}_{label}bp_fps.bed.gz'
                final_bedgraph = tmp_dir / f'combined_{chrom}_{label}bp_fp_cov.bedgraph'
                final_bedgraph_gz = chrom_dir / f'combined_{chrom}_{label}bp_fp_cov.bedgraph.gz'
                final_bigwig = chrom_dir / f'combined_{chrom}_{label}bp_fp_cov.bw'

                sort_to_file(chrom_unsorted[label], sorted_bed)
                bgzip_file(sorted_bed, final_bed_gz)
                run(['tabix', '-f', '-p', 'bed', str(final_bed_gz)])

                make_bedgraph(sorted_bed, final_bedgraph)
                if Path(BEDGRAPH_TO_BIGWIG).exists():
                    run([BEDGRAPH_TO_BIGWIG, str(final_bedgraph), str(CHROM_SIZES), str(final_bigwig)])
                else:
                    print(f'Skipping bigWig creation because {BEDGRAPH_TO_BIGWIG} was not found')
                bgzip_file(final_bedgraph, final_bedgraph_gz)

                print(f'Wrote {chrom} {label}: {count} intervals')

        for handle in all_chrom_handles.values():
            handle.close()

        for label in SIZE_BINS:
            count = all_chrom_counts[label]
            if count == 0:
                print(f'Skipping all-chromosome output for {label}: no intervals found')
                continue

            sorted_bed = tmp_dir / f'combined_all_chrs_{label}bp_fps.sorted.bed'
            final_bed_gz = COMBINED_PEAKS_DIR / f'combined_all_chrs_{label}bp_fps.bed.gz'

            sort_to_file(all_chrom_unsorted[label], sorted_bed)
            bgzip_file(sorted_bed, final_bed_gz)
            run(['tabix', '-f', '-p', 'bed', str(final_bed_gz)])
            print(f'Wrote all-chromosome {label}: {count} intervals')
    finally:
        for handle in all_chrom_handles.values():
            if not handle.closed:
                handle.close()

print('Done. New combined outputs are under:')
print(COMBINED_PEAKS_DIR)
print(COMBINED_TRACKS_DIR)
