#!/usr/bin/env python3
"""Run FIMO once on unique CREs with opaque sequence identifiers."""
import argparse
from pathlib import Path
import subprocess
import tempfile
import pysam
import sys
from LCL_tf_cooccupancy import REF, read, annotate_motifs


def run_scan(out, fimo, motif, reference):
    if (out/'motifs.tsv.gz').is_file(): return
    # FASTA and every raw FIMO side product live only in scratch and are removed
    # on both successful completion and exceptions.
    with tempfile.TemporaryDirectory(prefix='lcl_fimo_') as directory:
        temp=Path(directory)
        cre=read(out/'cres.tsv.gz')
        with pysam.FastaFile(str(reference)) as fa, (temp/'regions.fa').open('w') as stream:
            for r in cre.itertuples():
                sequence=fa.fetch(r.chrom,r.start,r.end).upper()
                if len(sequence)!=r.end-r.start: raise ValueError('Incomplete reference sequence')
                stream.write(f'>{r.sequence_name}\n{sequence}\n')
        result=subprocess.run([str(fimo),'--thresh','1e-4','--max-stored-scores','10000000',
            '--oc',str(temp/'scan'),str(motif),str(temp/'regions.fa')],capture_output=True,text=True)
        if result.stdout: print(result.stdout, end='')
        if result.stderr: print(result.stderr, end='', file=sys.stderr)
        result.check_returncode()
        if any(word in result.stderr.lower() for word in ['discard','maximum number','max-stored-scores']):
            raise RuntimeError('Possible FIMO truncation; see scheduler stderr')
        annotate_motifs(out,temp/'scan/fimo.tsv')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--fimo',type=Path,required=True)
    p.add_argument('--motif',type=Path,required=True)
    p.add_argument('--reference',type=Path,default=Path(REF));a=p.parse_args()
    run_scan(a.out,a.fimo,a.motif,a.reference)

if __name__=='__main__':main()
