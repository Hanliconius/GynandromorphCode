#!/usr/bin/env python3
"""Optional reference homology search; submit on Pegasus, never on the login node."""
import argparse
import csv
import json
from pathlib import Path
import shlex
import subprocess
import sys
from run_breakpoint import stamp


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-dir',required=True,type=Path)
    p.add_argument('--submit',action='store_true')
    p.add_argument('--partition',default='cpu')
    p.add_argument('--time',default='02:00:00')
    p.add_argument('--blastn',default='blastn')
    a=p.parse_args();out=a.run_dir.resolve()
    cfg=json.loads((out/'run.json').read_text())
    if not (out/'complete.json').is_file():raise RuntimeError('Run the context pipeline first')
    if a.submit:
        cmd=[sys.executable,str(Path(__file__).resolve()),'--run-dir',str(out),'--blastn',a.blastn]
        job=subprocess.check_output(['sbatch','--parsable','--job-name=NDJ_homology',
            '--partition='+a.partition,'--time='+a.time,'--cpus-per-task=1','--mem=12G',
            '--output='+str(out/'logs/homology.%j.out'),'--error='+str(out/'logs/homology.%j.err'),
            '--wrap='+shlex.join(cmd)],text=True).strip()
        print('JOB\t'+job);return
    dest=out/'homology';dest.mkdir(exist_ok=True)
    inputs=[stamp(out/'homology_queries.fa',True),stamp(cfg['reference']),stamp(Path(__file__),True)]
    if (dest/'complete.json').exists():
        old=json.loads((dest/'complete.json').read_text())
        if old['inputs']!=inputs:raise RuntimeError('Homology inputs/code changed; use a fresh run')
        for saved in old['outputs']:
            if stamp(saved['path'],True)!=saved:raise RuntimeError('Homology output changed')
        print('REUSED\t'+str(dest));return
    raw=dest/'all_hits.tsv'
    fields='qseqid sseqid pident length qstart qend sstart send evalue bitscore'
    command=[a.blastn,'-task','blastn','-query',str(out/'homology_queries.fa'),'-subject',cfg['reference'],
        '-dust','yes','-soft_masking','true','-evalue','1e-10','-max_target_seqs','10000',
        '-outfmt','6 '+fields,'-out',str(raw)]
    version=subprocess.check_output([a.blastn,'-version'],text=True)
    subprocess.run(command,check=True)
    rows=[]
    with raw.open() as f:
        for line in f:
            q,s,ident,length,qa,qb,sa,sb,e,bits=line.rstrip().split('\t')
            if float(ident)<90 or int(length)<200:continue
            chrom,span=q.rsplit(':',1);offset=int(span.split('-')[0])
            qs=offset+min(int(qa),int(qb))-1;qe=offset+max(int(qa),int(qb))
            ss=min(int(sa),int(sb))-1;se=max(int(sa),int(sb))
            # Remove the diagonal self match, retain off-diagonal and inverted hits.
            strand='+' if (int(qb)-int(qa))*(int(sb)-int(sa))>=0 else '-'
            if s==chrom and strand=='+' and abs(ss-qs)<=10 and abs(se-qe)<=10:continue
            rows.append([chrom,qs,qe,s,ss,se,strand,ident,length,e,bits])
    with (dest/'nonself_homology.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n')
        w.writerow(['query_chrom','query_start0','query_end0','partner_chrom','partner_start0','partner_end0',
                    'relative_orientation','identity_pct','alignment_bp','evalue','bitscore']);w.writerows(rows)
    (dest/'complete.json').write_text(json.dumps(dict(command=command,version=version,hits=len(rows),inputs=inputs,
        outputs=[stamp(raw,True),stamp(dest/'nonself_homology.tsv',True)],
        interpretation='Candidate reference homologies only; overlapping 5-kb queries duplicate hits. No NAHR call.',
        limitations='>=200 aligned bp, >=90% identity; masked seeding can miss low-complexity homology; no control enrichment test.'),indent=2)+'\n')
    print('COMPLETE\t'+str(dest))


if __name__=='__main__':main()
