#!/usr/bin/env python3
"""Convert existing same-reference annotations to union BED tracks; never download."""
import argparse
import csv
from collections import defaultdict
import gzip
import json
from pathlib import Path
from sequence_features import Track


def op(path):
    return gzip.open(path,'rt') if str(path).endswith('.gz') else open(path)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference',required=True,type=Path)
    p.add_argument('--repeatmasker-out',type=Path)
    p.add_argument('--gff',type=Path)
    p.add_argument('--bed',action='append',default=[],help='NAME=BED; e.g. low_mappability=...')
    p.add_argument('--output',required=True,type=Path)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    if (a.output/'annotations.json').exists():raise RuntimeError('Use a new output directory')
    tracks=defaultdict(lambda:defaultdict(list));sources=[];details=[]
    if a.repeatmasker_out:
        sources.append(str(a.repeatmasker_out.resolve()))
        with op(a.repeatmasker_out) as f:
            for line in f:
                v=line.split()
                if not v or not v[0].isdigit():continue
                if len(v)<11:raise ValueError('Malformed RepeatMasker .out row')
                chrom,start,end,cls=v[4],int(v[5])-1,int(v[6]),v[10].split('/')[0]
                if start<0 or end<=start:raise ValueError('Malformed RepeatMasker coordinates')
                tracks['repeatmasker_all'][chrom].append((start,end))
                name='RM_'+''.join(c if c.isalnum() else '_' for c in cls)
                tracks[name][chrom].append((start,end))
                details.append([chrom,start,end,v[9],v[10],v[8],v[1]])
    if a.gff:
        sources.append(str(a.gff.resolve()))
        with op(a.gff) as f:
            for line in f:
                if line.startswith('##FASTA'):break
                if line.startswith('#') or not line.strip():continue
                v=line.rstrip().split('\t')
                if len(v)!=9:raise ValueError('Malformed GFF3')
                if v[2] in {'gene','exon','CDS'}:
                    start,end=int(v[3])-1,int(v[4])
                    if start<0 or end<=start:raise ValueError('Malformed GFF coordinates')
                    tracks['GFF_'+v[2]][v[0]].append((start,end))
                    if v[2]=='gene' and end-start>=100000:
                        tracks['GFF_gene_ge100kb'][v[0]].append((start,end))
    for arg in a.bed:
        name,path=arg.split('=',1)
        if not name.replace('_','').isalnum() or name in tracks:raise ValueError('Unsafe/duplicate name')
        sources.append(str(Path(path).resolve()))
        with op(path) as f:
            for line in f:
                if line.startswith(('#','track','browser')) or not line.strip():continue
                v=line.split();tracks[name][v[0]].append((int(v[1]),int(v[2])))
    if not tracks:raise ValueError('No annotations supplied')
    if not any('NC_140403.1' in c for c in tracks.values()):
        raise ValueError('No chr15 accession in annotations; wrong assembly/contig names?')
    with (a.output/'repeatmasker_details.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n')
        w.writerow(['chrom','start0','end0','family','class','strand','divergence_pct']);w.writerows(details)
    items=[]
    for name,chroms in sorted(tracks.items()):
        path=a.output/(name+'.bed')
        with path.open('w') as f:
            for chrom,intervals in sorted(chroms.items()):
                t=Track(intervals)
                for start,end in zip(t.starts,t.ends):f.write(f'{chrom}\t{start}\t{end}\n')
        items.append(dict(name=name,path=str(path.resolve()),allow_empty_chr15=True))
    manifest=dict(assembly='GCF_056956725.1',reference=str(a.reference.resolve(strict=True)),
        sources=sources,coordinate_system='0-based half-open BED',tracks=items,
        caveat='Caller must supply annotations from this exact assembly; names alone do not validate liftover')
    (a.output/'annotations.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(a.output/'annotations.json')


if __name__=='__main__':main()
