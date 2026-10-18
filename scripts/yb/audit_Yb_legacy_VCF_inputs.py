#!/usr/bin/env python3
"""Find legacy Yb candidates by VCF contents, without assuming filenames or coordinate equivalence."""
import csv,gzip,json,re
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
def main():
    meta=list(csv.DictReader((ROOT/'data/derived/region_panel_metadata.tsv').open(),delimiter='\t'))
    cydno={r.get('Accession','').strip() for r in meta if re.search(r'cydno|pachinus',r.get('Taxon Name',''),re.I)}
    files=sorted({p for folder in [ROOT/'external/krista_genotype_plots_original',ROOT/'analysis/original_502'] for p in folder.glob('**/*') if p.is_file() and (p.name.endswith('.vcf') or p.name.endswith('.vcf.gz'))})
    out=ROOT/'analysis/yb_legacy_input_audit_v1'
    out.mkdir(exist_ok=False)
    results=[]
    for path in files:
        opener=gzip.open if path.name.endswith('.gz') else open
        samples=[];contigs=set();count=0;spans={};formats=set()
        scan=path.stat().st_size<250_000_000
        with opener(path,'rt') as h:
            for line in h:
                if line.startswith('##FORMAT='):formats.add(line.split('ID=',1)[1].split(',',1)[0])
                elif line.startswith('#CHROM'):
                    samples=line.rstrip().split('\t')[9:]
                    if not scan:break
                elif not line.startswith('#'):
                    f=line.split('\t',9);chrom=f[0];pos=int(f[1]);contigs.add(chrom);count+=1
                    if chrom not in spans:spans[chrom]=[pos,pos]
                    else:spans[chrom]=[min(spans[chrom][0],pos),max(spans[chrom][1],pos)]
        entry=dict(path=str(path),samples=len(samples),sample_IDs=samples,cydno_pachinus_IDs=sorted(set(samples)&cydno),formats=sorted(formats),record_scan_complete=scan,records=count if scan else None,contig_spans=spans)
        results.append(entry)
        if set(samples)&cydno or re.search(r'cortex|ivory|yb',str(path),re.I):
            print('CANDIDATE\t'+json.dumps({k:v for k,v in entry.items() if k not in ('sample_IDs','cydno_pachinus_IDs')})+'\tCYDNO_N='+str(len(entry['cydno_pachinus_IDs'])))
    (out/'vcf_inventory.json').write_text(json.dumps(results,indent=2)+'\n')
    print('REPORT\t'+str(out/'vcf_inventory.json'))
    print('NOTE\tCoordinate/reference compatibility is not yet established; no samples merged or calls changed.')
if __name__=='__main__':main()
