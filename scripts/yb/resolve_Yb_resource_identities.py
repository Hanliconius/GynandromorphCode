#!/usr/bin/env python3
"""Resolve audit BioSample names to metadata accessions; never infer phenotype."""
import argparse,csv,io,json,re,urllib.parse,urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
def read(p):
    with p.open() as h:return list(csv.DictReader(h,delimiter='\t'))
def write(p,rows,fields):
    with p.open('w',newline='') as h:
        w=csv.DictWriter(h,delimiter='\t',fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)
def lookup(s):
    params={'accession':s,'result':'read_run','fields':'run_accession,sample_accession,secondary_sample_accession,sample_alias,scientific_name','format':'tsv'}
    url='https://www.ebi.ac.uk/ena/portal/api/filereport?'+urllib.parse.urlencode(params)
    try:
        with urllib.request.urlopen(url,timeout=90) as h:raw=h.read().decode()
        rows=list(csv.DictReader(io.StringIO(raw),delimiter='\t'))
        if not rows or 'run_accession' not in rows[0]:raise ValueError('No ENA run records')
        return s,rows,raw,''
    except Exception as e:return s,[], '',str(e)
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=ROOT);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    out=a.output;out.mkdir(parents=True,exist_ok=False);(out/'ena').mkdir()
    audit=a.root/'analysis/yb_existing_resources_v1';bams=read(audit/'bam_candidates.tsv');meta=read(a.root/'data/derived/region_panel_metadata.tsv')
    ids=sorted({r['sample'] for r in bams if r['sample'].startswith('SAMEA')})
    index={}
    for row in meta:
        for col in ('Accession','SequenceID','EarthCapeID'):
            for value in re.findall(r'(?:ERR|SRR|DRR|SAMEA|SAMN)\d+',row.get(col,'')):index.setdefault(value,[]).append(row)
    resolved=[];errors=[]
    with ThreadPoolExecutor(max_workers=4) as pool:
        for sample,records,raw,error in pool.map(lookup,ids):
            if raw:(out/'ena'/f'{sample}.tsv').write_text(raw)
            if error:errors.append({'sample':sample,'error':error})
            hits={}
            for rec in records:
                # Run-to-BioSample relationship must explicitly match the requested sample.
                related={rec.get('sample_accession',''),rec.get('secondary_sample_accession','')}
                if sample not in related:continue
                run=rec.get('run_accession','')
                for m in index.get(run,[]):
                    acc=m.get('Accession','')
                    key=(acc,m.get('Taxon Name',''),m.get('Yb-dorsal',''),m.get('Yb-ventral',''))
                    hits[key]=m
            # Multiple runs are never presumed to be distinct butterflies.
            phenos={(k[1],k[2],k[3]) for k in hits}
            status='metadata_matched' if hits and len(phenos)==1 else 'metadata_conflict' if hits else 'unresolved'
            m=next(iter(hits.values())) if status=='metadata_matched' else {}
            paths=sorted({r['resolved_path'] for r in bams if r['sample']==sample and r['reference']=='HmRos' and r['resolved_path'].endswith('.dedup.bam') and '/tmp/' not in r['resolved_path']})
            resolved.append(dict(bam_sample=sample,run_accessions=';'.join(sorted({r.get('run_accession','') for r in records if sample in {r.get('sample_accession',''),r.get('secondary_sample_accession','')}})),metadata_accessions=';'.join(sorted({k[0] for k in hits})),taxon=' '.join(m.get('Taxon Name','').split()),Yb_dorsal=m.get('Yb-dorsal',''),Yb_ventral=m.get('Yb-ventral',''),status=status,final_Rosina_BAMs=';'.join(paths)))
    write(out/'identity_resolution.tsv',resolved,['bam_sample','run_accessions','metadata_accessions','taxon','Yb_dorsal','Yb_ventral','status','final_Rosina_BAMs'])
    write(out/'lookup_errors.tsv',errors,['sample','error'])
    # Preserve local provenance evidence, including NDJ, without dumping it to stdout.
    evidence=[]
    for folder in (a.root/'analysis/genomewide_mixing_corrected_v1/full_sib_validation',a.root/'analysis/NDJ_control_comparison_v1'):
        if not folder.exists():continue
        for path in folder.rglob('*'):
            if not path.is_file() or path.suffix.lower() not in {'.tsv','.csv','.json','.txt'}:continue
            if any(x in path.parts for x in ('tmp','alignment','raw','logs')) or path.stat().st_size>5_000_000:continue
            for n,line in enumerate(path.read_text(errors='replace').splitlines(),1):
                if any(s in line for s in ids) or 'NDJ_' in line:evidence.append(dict(path=str(path),line=n,text=line[:4000]))
    write(out/'local_identity_evidence.tsv',evidence,['path','line','text'])
    summary={'status':'COMPLETE_WITH_UNRESOLVED' if errors or any(r['status']!='metadata_matched' for r in resolved) else 'COMPLETE','BioSamples':len(ids),'metadata_matched':sum(r['status']=='metadata_matched' for r in resolved),'lookup_errors':len(errors),'phenotypes_inferred':False,'alignments_or_calls_changed':False}
    (out/'complete.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    for r in resolved:print('\t'.join(str(r[k]) or '?' for k in ('bam_sample','metadata_accessions','taxon','Yb_dorsal','Yb_ventral','status')))
if __name__=='__main__':main()
