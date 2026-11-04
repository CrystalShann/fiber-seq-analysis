#!/usr/bin/env python3
"""Scan explicit promoter windows with the eight DNA-TF family panels (no U1).

Outputs retain every FIMO hit; the report selects the highest score per family.
Coordinates: BED genomic intervals; transcription-oriented relative base positions.
"""
import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import pysam

PROJECT = Path('/project/spott/cshan/fiber-seq')
PANELS = {
    'SP': {'SP1', 'SP2', 'SP3', 'SP4', 'SP5', 'SP8', 'SP9'},
    'ETS': {'ETS1', 'ETS2', 'ELF1', 'ELF2', 'ELF3', 'ELF4', 'ELF5', 'ETV1',
            'ETV2', 'ETV3', 'ETV4', 'ETV5', 'ETV6', 'ETV7', 'ERG', 'FLI1',
            'GABPA', 'SPI1', 'SPIB', 'SPIC', 'SPDEF', 'EHF', 'FEV', 'ERF'},
    'ZNF143': {'ZNF143'}, 'NRF1': {'NRF1'},
    'CREB': {'CREB1', 'CREB3', 'CREB3L1', 'CREB3L2', 'CREB3L3', 'CREB3L4',
             'CREB5', 'ATF1', 'ATF2', 'ATF3', 'ATF4', 'ATF6', 'ATF7'},
    'TATA': {'TBP'}, 'YY1': {'YY1'}, 'NFY': {'NFYA', 'NFYB', 'NFYC'},
}

def family_sites(hits):
    """Non-transitive, best-p-value representative intervals; never union intervals.

    Collapse exact duplicate intervals regardless of strand. Different matrices
    can share a site only with >=80% reciprocal overlap AND centers <=2 bp apart.
    Shifted hits of the same matrix remain distinct. Every member must match the
    representative directly, so neighboring sites cannot merge by chaining.
    """
    rows, membership = [], {}
    for (_, category), group in hits.groupby(['tss_id', 'category'], sort=True):
        representatives = []
        for row in group.sort_values(['p_value', 'matrix_id', 'start', 'end', 'motif_strand']).to_dict('records'):
            match = None
            for rep in representatives:
                same = row['start'] == rep['start'] and row['end'] == rep['end']
                overlap = max(0, min(row['end'], rep['end']) - max(row['start'], rep['start']))
                near = abs(row['start'] + row['end'] - rep['start'] - rep['end']) <= 4
                similar = (row['matrix_id'] != rep['matrix_id'] and near
                           and overlap >= .8 * max(row['end']-row['start'], rep['end']-rep['start']))
                if same or similar:
                    match = rep
                    break
            if match is None:
                match = dict(row)
                match['site_id'] = f'site{len(rows) + 1:08d}'
                match['members'] = []
                rows.append(match)
                representatives.append(match)
            match['members'].append(row)
            membership[row['matrix_hit_id']] = match['site_id']
    result = []
    for rep in rows:
        members = rep.pop('members')
        rep['representative_matrix_id'] = rep.pop('matrix_id')
        rep['matrix_ids'] = ';'.join(sorted({r['matrix_id'] for r in members}))
        rep['matrix_hit_ids'] = ';'.join(r['matrix_hit_id'] for r in members)
        rep['n_matrix_hits'] = len(members)
        rep['orientations_relative_to_transcription'] = ';'.join(sorted({r['orientation_relative_to_transcription'] for r in members}))
        result.append(rep)
    return pd.DataFrame(result), membership


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    default = PROJECT / 'LCL_project/pausing/clustering_Pol2/tables/promoter_motif_cluster_examples/examples'
    parser.add_argument('--promoters', type=Path, default=default/'selected_promoter_windows.tsv')
    parser.add_argument('--output', type=Path, default=default)
    parser.add_argument('--reference', type=Path, default=Path('/project/spott/cshan/annotations/hg38.fa'))
    parser.add_argument('--database', type=Path, default=Path('/project/spott/cshan/annotations/JASPAR2026_CORE_vertebrates_non-redundant_pfms_meme.txt'))
    parser.add_argument('--fimo', type=Path, default=Path('/project/spott/cshan/tools/meme-5.5.9-install/bin/fimo'))
    args = parser.parse_args()
    for path in (args.promoters, args.reference, Path(str(args.reference)+'.fai'), args.database, args.fimo):
        if not path.is_file():
            raise FileNotFoundError(path)
    args.output.mkdir(parents=True, exist_ok=True)
    def write(df, name):
        df.to_csv(args.output/name, sep='\t', index=False, na_rep='NA')
    promoters = pd.read_csv(args.promoters, sep='\t')
    required = ['tss_id', 'chrom', 'window_start', 'window_end', 'strand']
    if not set(required) <= set(promoters) or promoters.tss_id.duplicated().any():
        raise ValueError('Promoters require unique tss_id and explicit chrom/window_start/window_end/strand')
    if not ((promoters.window_end-promoters.window_start).eq(2000).all() and promoters.strand.isin(['+','-']).all()):
        raise ValueError('Expected strand-oriented 2000-bp promoter windows')
    chunks = re.split(r'(?=^MOTIF )', args.database.read_text(), flags=re.MULTILINE)
    manifest, chosen = [], []
    for block in chunks[1:]:
        fields = block.splitlines()[0].split()
        matrix_id, name = fields[1:3]
        for category, names in PANELS.items():
            if name.upper() not in names or (category == 'NRF1' and not matrix_id.startswith('MA0506.')):
                continue
            width = int(re.search(r'letter-probability matrix:.*?\bw\s*=\s*(\d+)', block).group(1))
            manifest.append(dict(category=category, matrix_id=matrix_id, name=name, width=width,
                                 status='available', source_url=f'https://jaspar.elixir.no/matrix/{matrix_id}/'))
            chosen.append(block)
    manifest = pd.DataFrame(manifest)
    assert set(manifest.category) == set(PANELS)
    write(manifest, 'motif_manifest.tsv')
    meme = args.output/'selected_motifs.meme'
    meme.write_text(chunks[0]+''.join(chosen))
    background = args.output/'background.txt'
    background.write_text('A 0.25\nC 0.25\nG 0.25\nT 0.25\n')
    keys = ['chrom','window_start','window_end','strand']
    windows = promoters[keys].drop_duplicates().sort_values(keys).reset_index(drop=True)
    windows['sequence_id'] = [f'prom{i+1:06d}' for i in range(len(windows))]
    complement = str.maketrans('ACGTN','TGCAN')
    seqs = []
    with pysam.FastaFile(str(args.reference)) as ref:
        for row in windows.itertuples():
            seq = ref.fetch(row.chrom, int(row.window_start), int(row.window_end)).upper()
            seqs.append(seq if row.strand == '+' else seq.translate(complement)[::-1])
    windows['usable_sequence'] = [len(s)==2000 and set(s)<=set('ACGT') for s in seqs]
    windows['sequence_status'] = np.where(windows.usable_sequence,'usable','unresolved_length_or_ambiguous_bases')
    lookup = promoters.merge(windows,on=keys,validate='many_to_one')
    write(lookup,'promoter_sequence_lookup.tsv')
    fasta = args.output/'assigned_promoters.fa'
    with fasta.open('w') as handle:
        for row, seq in zip(windows.itertuples(), seqs):
            if row.usable_sequence:
                handle.write(f'>{row.sequence_id}\n{seq}\n')
    if not windows.usable_sequence.any():
        raise ValueError('No usable promoter sequences to scan')
    cmd = [str(args.fimo),'--thresh','1e-4','--max-stored-scores','10000000','--no-qvalue',
           '--motif-pseudo','0.1','--bfile',str(background),'--oc',str(args.output/'fimo'),str(meme),str(fasta)]
    with (args.output/'fimo.log').open('w') as log:
        subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
    raw = pd.read_csv(args.output/'fimo/fimo.tsv',sep='\t',comment='#').rename(columns={
        'motif_id':'matrix_id','sequence_name':'sequence_id','strand':'motif_strand',
        'start':'fimo_start','stop':'fimo_stop','p-value':'p_value'})
    raw = raw.merge(manifest[['matrix_id','category']],on='matrix_id',validate='many_to_one')
    raw = raw.merge(lookup,on='sequence_id',validate='many_to_many')
    plus = raw.strand.eq('+')
    raw['start'] = np.where(plus,raw.window_start+raw.fimo_start-1,raw.window_end-raw.fimo_stop)
    raw['end'] = np.where(plus,raw.window_start+raw.fimo_stop,raw.window_end-raw.fimo_start+1)
    raw['relative_start'] = raw.fimo_start-1001
    raw['relative_end'] = raw.fimo_stop-1001
    raw['relative_center'] = (raw.relative_start+raw.relative_end)/2
    raw['genomic_strand'] = np.where(plus,raw.motif_strand,raw.motif_strand.map({'+':'-','-':'+'}))
    raw['orientation_relative_to_transcription'] = raw.motif_strand
    raw['matrix_hit_id'] = [f'hit{i+1:09d}' for i in range(len(raw))]
    if len(raw):
        sites, membership = family_sites(raw)
        raw['site_id'] = raw.matrix_hit_id.map(membership)
    else:
        raw['site_id'] = pd.Series(dtype=str)
        sites = raw.rename(columns={'matrix_id':'representative_matrix_id'}).copy()
        for column in ['matrix_ids','matrix_hit_ids','n_matrix_hits','orientations_relative_to_transcription']:
            sites[column] = pd.Series(dtype=str)
    write(raw,'matrix_motif_hits.tsv.gz')
    write(sites,'motif_sites.tsv.gz')
    provenance = dict(created_utc=datetime.now(timezone.utc).isoformat(),
        promoter_input=str(args.promoters), promoter_sha256=hashlib.sha256(args.promoters.read_bytes()).hexdigest(),
        reference=str(args.reference),database=str(args.database),
        database_sha256=hashlib.sha256(args.database.read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        fimo_command=cmd,p_value_threshold=1e-4,motif_pseudocount=.1,
        uniform_background=True,both_strands=True,q_values=False,
        n_promoters=len(promoters),n_usable_windows=int(windows.usable_sequence.sum()),
        n_matrices=len(manifest),n_matrix_hits=len(raw),
        coordinates='Genomic BED [start,end); relative first/last bases inclusive -1000..999; minus promoters reverse-complemented.',
        display='Highest FIMO score per family per promoter; ties: p-value, matrix ID, position, strand. All raw hits retained.',
        U1='Omitted; no sequence model supplied.')
    (args.output/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(f'Scanned {len(promoters)} promoters; {len(manifest)} matrices; {len(raw)} hits',flush=True)


if __name__ == '__main__':
    main()
