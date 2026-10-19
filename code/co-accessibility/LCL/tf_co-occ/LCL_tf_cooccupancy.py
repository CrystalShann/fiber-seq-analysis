#!/usr/bin/env python3
"""LCL motif-supported same-fiber occupancy. BED coordinates are zero-based/half-open."""
import argparse
from collections import defaultdict
from pathlib import Path
import sys
import tempfile

import numpy as np
import pandas as pd
import pysam
# Shared co-accessibility helpers remain in the parent LCL directory.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from LCL_coaccess_cres import incidence_for_chrom

ROOT = Path('/project/spott/cshan/fiber-seq/LCL_project/co-accessibility')
TFROOT = Path('/project/spott/1_Shared_projects/LCL_Fiber_seq/FiberHMM/FiberHMM_v2/results/results_v1_trained_model')
REF = '/project/spott/reference/human/GRCh38.p14/hg38.fa'
CELLS = ['co_closed', 'CRE1_access', 'CRE2_access', 'co_access']


def read(path):
    return pd.read_csv(path, sep='\t')


def write(df, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.partial')
    try:
        df.to_csv(tmp, sep='\t', index=False, compression='gzip' if path.suffix == '.gz' else None)
        tmp.replace(path)
    finally:
        tmp.unlink(missing_ok=True)


def label_pairs(pooled, adjusted, minimum=20):
    raw = pooled[CELLS].sum(axis=1)
    p = pooled.copy()
    p['raw_spanning'] = raw
    fields = ['CRE_pair', 'OR_adj', 'fdr_adj', 'co_access', 'expected_co_access_adj', 'fit_status']
    a = adjusted[fields].rename(columns={'co_access': 'co_access_adj'})
    p = p.merge(a, on='CRE_pair', how='left', validate='one_to_one')
    valid = p.fit_status.isin(['ok', 'ok_with_fit_warnings', 'separation_positive', 'separation_negative'])
    valid &= np.isfinite(p.fdr_adj) & np.isfinite(p.OR_adj) & np.isfinite(p.expected_co_access_adj)
    positive = (p.fdr_adj < .05) & (p.OR_adj > 1) & (p.co_access_adj > p.expected_co_access_adj)
    p['set_adjusted'] = np.where(valid, np.where(positive, 'co-access', 'background'), 'excluded')
    pv = np.isfinite(p.fdr) & p.OR.notna() & np.isfinite(p.expected_co_access)
    pp = (p.fdr < .05) & (p.OR > 1) & (p.co_access > p.expected_co_access)
    p['set_pooled'] = np.where(pv, np.where(pp, 'co-access', 'background'), 'excluded')
    p['exclusion'] = np.where(raw < minimum, 'low_raw_coverage', np.where(valid, '', 'invalid_adjusted_fit'))
    p.loc[p.raw_spanning < minimum, ['set_adjusted', 'set_pooled']] = 'excluded'
    return p


def cre_sets(pairs, column):
    valid = pairs[pairs[column] != 'excluded']
    pos = valid[valid[column] == 'co-access']
    yes = set(pos.CRE1) | set(pos.CRE2)
    all_ids = set(valid.CRE1) | set(valid.CRE2)
    return {c: ('co-access' if c in yes else 'background') for c in all_ids}


def prepare(args):
    out = args.out
    paths = [args.root/'coaccess/LCL_coaccess_pairs.tsv.gz',
        args.root/'coaccess/donor_adjusted/all_chr/LCL_coaccess_pairs_donor_adjusted_all_chr.tsv.gz',
        args.root/'universe/cre_universe.bed', args.root/'universe/sample_manifest.tsv',
        args.root/'coaccess/sample_manifest.tsv']
    products = [out/x for x in ['pairs.tsv.gz', 'cres.tsv.gz', 'samples.tsv']]
    out.mkdir(parents=True, exist_ok=True)
    pooled, adj = read(paths[0]), read(paths[1])
    pairs = label_pairs(pooled, adj, args.min_raw)
    write(pairs, products[0])
    retention = []
    for n in [10, 20, 50]:
        labeled = label_pairs(pooled, adj, n)
        for group, count in labeled.groupby('set_adjusted').size().items():
            retention.append((n, group, count))
    write(pd.DataFrame(retention, columns=['min_raw', 'set', 'n_pairs']), out/'coverage_retention.tsv')
    write(pairs.raw_spanning.describe(percentiles=[.01,.05,.25,.5,.75,.95,.99]).rename_axis('statistic').reset_index(), out/'coverage_distribution.tsv')
    cre = pd.read_csv(paths[2], sep='\t', header=None, names=['chrom','start','end','CRE_ID','CRE_class'])
    for analysis in ['adjusted', 'pooled']:
        cre['set_'+analysis] = cre.CRE_ID.map(cre_sets(pairs, 'set_'+analysis)).fillna('excluded')
    cre = cre[(cre.set_adjusted != 'excluded') | (cre.set_pooled != 'excluded')].copy()
    cre = cre.sort_values(['chrom','start','end','CRE_ID']).reset_index(drop=True)
    cre['sequence_name'] = ['region_%07d' % i for i in range(len(cre))]
    write(cre, products[1])
    manifest = read(paths[4])
    manifest['tf_path'] = [str(args.tf_root/s/(s+'.recalled_tf.bed.gz')) for s in manifest.sample_name]
    write(manifest, products[2])


def collapse_hits(hits):
    kept = []
    for _, group in hits.groupby(['CRE_ID', 'motif_id'], sort=True):
        accepted = []
        for r in group.sort_values(['score','start','end','strand'], ascending=[False,True,True,True]).itertuples():
            if not any(r.start < b and a < r.end for a,b in accepted):
                kept.append(r.Index); accepted.append((r.start,r.end))
    return hits.loc[kept].sort_values(['chrom','start','end','motif_id','strand']).reset_index(drop=True)


def annotate_motifs(out, fimo_tsv):
    """Publish annotated motifs; the scanner owns temporary hits."""
    hits = pd.read_csv(fimo_tsv, sep='\t', comment='#')
    hits = hits.dropna(subset=['motif_id', 'sequence_name'])
    cre = read(out/'cres.tsv.gz')
    hits = hits.rename(columns={'start':'relative_start','stop':'relative_stop','motif_alt_id':'TF_label'})
    hits = hits.merge(cre[['sequence_name','CRE_ID','chrom','start','end']], on='sequence_name', validate='many_to_one')
    hits['end'] = hits.start + hits.relative_stop.astype(int)
    hits['start'] = hits.start + hits.relative_start.astype(int)-1
    hits = collapse_hits(hits)
    hits['instance_id'] = ['motif_%09d' % i for i in range(len(hits))]
    write(hits, out/'motifs.tsv.gz')


def footprint_blocks(fields):
    """Parse BED15 without confusing parent-record score and block scores."""
    n = int(fields[9]); base = int(fields[1])
    lengths = [int(x) for x in fields[10].rstrip(',').split(',') if x]
    offsets = [int(x) for x in fields[11].rstrip(',').split(',') if x]
    scores = [float(x) for x in fields[12].rstrip(',').split(',') if x]
    if not len(lengths) == len(offsets) == len(scores) == n: raise ValueError('Malformed BED15')
    return [(base+o,base+o+l) for l,o,s in zip(lengths,offsets,scores) if 10 <= l <= 80 and s >= 50]


def supported(start, end, footprints):
    return any(min(end,b)-max(start,a) >= .5*(end-start) for a,b in footprints)


def count(args):
    out = args.out
    samples = read(out/'samples.tsv')
    if not 1 <= args.sample_index <= len(samples): raise ValueError('sample-index must be 1..31')
    sample = samples.iloc[args.sample_index-1]
    dest = out/'counts'/sample.sample_name
    combined = [out/(x+'_sample_counts.tsv.gz') for x in ['cre','motif','pair']]
    if all(p.is_file() for p in combined): return
    dest.mkdir(parents=True, exist_ok=True)
    outputs = [dest/(x+'.tsv.gz') for x in ['cre','motif','pair']]
    if all(p.is_file() for p in outputs): return
    cre = read(out/'cres.tsv.gz'); motifs_df = read(out/'motifs.tsv.gz'); pairs = read(out/'pairs.tsv.gz')
    pairs = pairs[(pairs.set_adjusted != 'excluded') | (pairs.set_pooled != 'excluded')]
    cre_rows=[]; motif_rows=[]; pair_rows=[]
    with pysam.TabixFile(sample.tf_path) as tf, pysam.TabixFile(sample.elements_path) as fires:
        for chrom,cg in cre.groupby('chrom',sort=True):
            with tempfile.TemporaryDirectory(prefix='lcl_occupancy_') as tmp:
                bed = Path(tmp)/'cres.bed'
                cg[['chrom','start','end','CRE_ID','CRE_class']].to_csv(bed,sep='\t',index=False,header=False)
                cov,acc,_ = incidence_for_chrom(chrom,bed,sample.blocks_path,sample.elements_path,'contain',.5,Path(tmp))
            wanted = set().union(*(cov[c] for c in cg.CRE_ID))
            fire_by_read=defaultdict(list)
            if chrom in fires.contigs:
                for line in fires.fetch(chrom):
                    f=line.split('\t')
                    if f[3] in wanted: fire_by_read[f[3]].append((int(f[1]),int(f[2])))
            observed=set(); fp=defaultdict(list)
            if chrom in tf.contigs:
                for line in tf.fetch(chrom):
                    f=line.split('\t'); rid=f[3]
                    if rid not in wanted: continue
                    observed.add(rid)
                    fp[rid].extend((a,b) for a,b in footprint_blocks(f)
                        if any(x <= a and b <= y for x,y in fire_by_read[rid]))
            bound={}; informed={}
            mg = {c:g for c,g in motifs_df[motifs_df.chrom == chrom].groupby('CRE_ID')}
            for c in cg.CRE_ID:
                informed[c] = acc[c] & observed
                bound[c] = set()
                g = mg.get(c)
                if g is not None:
                    motif_bound=defaultdict(set)
                    for m in g.itertuples():
                        b = {rid for rid in informed[c] if supported(m.start,m.end,fp[rid])}
                        bound[c] |= b
                        motif_bound[m.motif_id] |= b
                    for mid,b in motif_bound.items():
                        motif_rows.append((c,mid,len(acc[c]),len(informed[c]),len(b),len(acc[c]-observed)))
                cre_rows.append((c,len(cov[c]),len(acc[c]),len(informed[c]),len(bound[c]),len(acc[c]-observed)))
            for p in pairs[pairs.chrom == chrom].itertuples():
                shared=cov[p.CRE1]&cov[p.CRE2]; a=acc[p.CRE1]&shared; b=acc[p.CRE2]&shared
                joint=a&b; info=joint&observed
                b1=bound[p.CRE1]&info; b2=bound[p.CRE2]&info
                pair_rows.append((p.CRE_pair,len(joint),len(info),len(b1&b2),len(b1-b2),len(b2-b1),len(info-b1-b2),len(joint-info)))
            print(sample.sample_name,chrom,len(cg),'CREs',flush=True)
    frames = [pd.DataFrame(cre_rows,columns=['CRE_ID','n_covered','n_accessible','n_informative','n_bound','n_missing_tf']),
        pd.DataFrame(motif_rows,columns=['CRE_ID','motif_id','n_accessible','n_informative','n_bound','n_missing_tf']),
        pd.DataFrame(pair_rows,columns=['CRE_pair','n_joint_accessible','n_informative','both','first_only','second_only','neither','n_missing_tf'])]
    for df,path in zip(frames,outputs):
        df['sample_name']=sample.sample_name; df['donor']=sample.donor; write(df,path)


def combine(args):
    out=args.out; samples=read(out/'samples.tsv'); frames={}
    reuse = all((out/(kind+'_sample_counts.tsv.gz')).is_file() for kind in ['cre','motif','pair'])
    for kind in ['cre','motif','pair']:
        if reuse:
            frames[kind] = read(out/(kind+'_sample_counts.tsv.gz'))
            continue
        parts=[]
        for s in samples.sample_name:
            d=out/'counts'/s
            parts.append(read(d/(kind+'.tsv.gz')))
        frames[kind]=pd.concat(parts,ignore_index=True)
        write(frames[kind],out/(kind+'_sample_counts.tsv.gz'))
    cre=read(out/'cres.tsv.gz'); motifs_df=read(out/'motifs.tsv.gz'); pairs=read(out/'pairs.tsv.gz')
    summary=frames['cre'].groupby('CRE_ID')[['n_covered','n_accessible','n_informative','n_missing_tf']].sum().reset_index()
    cre=cre.merge(summary,on='CRE_ID',validate='one_to_one')
    cre['prop_act']=cre.n_accessible/cre.n_covered.replace(0,np.nan)
    cre['length']=cre.end-cre.start
    copies=motifs_df.groupby('CRE_ID').size(); cre['motif_count']=cre.CRE_ID.map(copies).fillna(0).astype(int)
    cre['prop_act_bin']=pd.qcut(cre.prop_act,10,duplicates='drop',labels=False)
    write(cre,out/'cre_covariates.tsv.gz')
    pool=frames['pair'].groupby('CRE_pair')[['n_joint_accessible','n_informative','n_missing_tf']].sum().reset_index()
    pairs=pairs.merge(pool,on='CRE_pair',validate='one_to_one')
    for k in [1,2]:
        c=cre[['CRE_ID','CRE_class','prop_act','length','motif_count']].rename(columns={
            'CRE_ID':'CRE'+str(k),'CRE_class':'class'+str(k),'prop_act':'access'+str(k),
            'length':'length'+str(k),'motif_count':'motifs'+str(k)})
        pairs=pairs.merge(c,on='CRE'+str(k),validate='many_to_one')
    pairs['log_distance']=np.log1p(pairs.dist); pairs['log_joint']=np.log1p(pairs.n_joint_accessible)
    write(pairs,out/'pair_covariates.tsv.gz')
    for analysis in ['adjusted','pooled','autosomes']:
        col='set_pooled' if analysis=='pooled' else 'set_adjusted'
        p=pairs[(pairs[col]!='excluded')&(pairs.n_joint_accessible>=10)&(pairs.n_informative>=10)].copy()
        if analysis=='autosomes': p=p[~p.chrom.isin(['chrX','chrY'])]
        p['set']=p[col]
        dest=out/'analysis'/analysis
        write(p,dest/'pair_eligibility.tsv.gz')
    # Remove staging tables after the combined analysis tables are written.
    for sample in samples.sample_name:
        dest = out/'counts'/sample
        for name in ['cre.tsv.gz','motif.tsv.gz','pair.tsv.gz']:
            (dest/name).unlink(missing_ok=True)
        if dest.exists() and not any(dest.iterdir()): dest.rmdir()
    parent = out/'counts'
    if parent.exists() and not any(parent.iterdir()): parent.rmdir()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('stage',choices=['prepare','count','combine'])
    p.add_argument('--root',type=Path,default=ROOT); p.add_argument('--out',type=Path)
    p.add_argument('--tf-root',type=Path,default=TFROOT)
    p.add_argument('--min-raw',type=int,default=20)
    p.add_argument('--sample-index',type=int,default=1)
    args=p.parse_args(); args.out=args.out or args.root/'tf_co-occ'
    if args.min_raw < 1: p.error('Coverage must be positive')
    if args.stage=='combine':
        for i in range(1,32):
            args.sample_index=i; count(args)
    globals()[args.stage](args)

if __name__=='__main__': main()
