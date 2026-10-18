#!/usr/bin/env python3
"""Audit original sequencing type, ENA strategy, and local reuse before expansion."""
import argparse,csv,io,json,re,time,urllib.parse,urllib.request
from concurrent.futures import ThreadPoolExecutor
from collections import Counter
from pathlib import Path
import pandas as pd
FIELDS='run_accession,sample_accession,library_layout,library_strategy,library_source,library_selection,experiment_title,study_title,library_construction_protocol,experiment_target,target_gene,base_count,read_count,fastq_bytes,bam_ftp,bam_bytes'
def fetch(ids):
 query=' OR '.join('run_accession="'+s+'"' for s in ids)
 url='https://www.ebi.ac.uk/ena/portal/api/search?'+urllib.parse.urlencode(dict(result='read_run',query='('+query+')',fields=FIELDS,format='tsv',limit=0))
 error=''
 for attempt in range(3):
  try:
   with urllib.request.urlopen(url,timeout=120) as h:raw=h.read().decode()
   rows=list(csv.DictReader(io.StringIO(raw),delimiter='\t'))
   if not rows or 'run_accession' not in rows[0]:raise RuntimeError('No valid ENA response')
   return rows,raw,''
  except Exception as e:
   error=str(e)
   if attempt<2:time.sleep(2)
 return [],'',error

def classify(original,row):
 strategy=row.get('library_strategy','').upper();selection=row.get('library_selection','').upper()
 protocol=row.get('library_construction_protocol','')
 evidence=' '.join([original,strategy,selection,protocol,row.get('experiment_target','')]).lower()
 if re.search(r'\brad\b|rad.?seq|restriction|ddrad|gbs|genotyping.by.sequencing',evidence):return 'reduced_representation'
 if re.search(r'capture|targeted|enrich|hybrid.?select|amplicon|amplicon|pcr.?amp',evidence) or strategy in {'AMPLICON','WXS','TARGETED-CAPTURE'}:return 'selective_capture_or_amplicon'
 if re.search(r'\bwgs\b|whole.?genome',original,re.I) or strategy=='WGS':return 'reported_WGS'
 return 'unresolved_sequencing_type'

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();out=a.output;out.mkdir(parents=True,exist_ok=False)
 cfg=json.loads((a.run/'run.json').read_text());tasks=cfg['new_tasks'];ids=[t['sample'] for t in tasks];root=Path(cfg['root']);meta=pd.read_csv(root/'data/derived/region_panel_metadata.tsv',sep='\t',dtype=str).fillna('');meta.columns=[c.strip() for c in meta.columns]
 chunks=[ids[i:i+20] for i in range(0,len(ids),20)];lookup={};errors=[]
 with ThreadPoolExecutor(max_workers=4) as pool:
  for i,(rows,raw,error) in enumerate(pool.map(fetch,chunks)):
   if raw:(out/f'ENA_batch_{i:03d}.tsv').write_text(raw)
   if error:errors.append(dict(samples=';'.join(chunks[i]),error=error))
   for row in rows:
    s=row['run_accession']
    if s in lookup:raise RuntimeError('Duplicate ENA run '+s)
    lookup[s]=row
 result=[]
 for t in tasks:
  s=t['sample'];mm=meta[meta['Accession'].str.contains(r'(?<!\d)'+re.escape(s)+r'(?!\d)',regex=True)]
  cols=[c for c in meta.columns if re.search(r'data.?type|sequenc|library',c,re.I)]
  original='; '.join(sorted({str(v).strip() for c in cols for v in mm[c] if str(v).strip()}))
  row=lookup.get(s,{});category=classify(original,row)
  # Actual reuse selection from the submitted manifest; no assumption that a
  # filename candidate has passed its ENA checksum yet.
  fq=t['fastq_candidates'];paired=bool(fq.get('1',fq.get(1,[]))) and bool(fq.get('2',fq.get(2,[])))
  r=dict(sample=s,taxon=t['taxon'],original_data_type=original,sequencing_category=category,ENA_resolved=bool(row),reuse_Rosina_BAM=bool(t['reuse_bam']),paired_local_FASTQ_candidates=paired,existing_download_complete=(a.run/'raw'/s/'complete.json').is_file(),original_WGS_ENA_non_WGS=bool(re.search(r'wgs|whole.?genome',original,re.I) and row.get('library_strategy','') not in ('','WGS')))
  r.update(row);result.append(r)
 pd.DataFrame(result).to_csv(out/'sequencing_and_reuse.tsv',sep='\t',index=False)
 pd.DataFrame(errors,columns=['samples','error']).to_csv(out/'lookup_errors.tsv',sep='\t',index=False)
 summary=dict(status='COMPLETE' if not errors and len(lookup)==len(ids) else 'COMPLETE_WITH_UNRESOLVED',targets=len(ids),ENA_resolved=len(lookup),original_data_type=dict(Counter(r['original_data_type'] or 'UNKNOWN' for r in result)),sequencing_categories=dict(Counter(r['sequencing_category'] for r in result)),ENA_strategies=dict(Counter(r.get('library_strategy','UNKNOWN') for r in result)),ENA_selections=dict(Counter(r.get('library_selection','UNKNOWN') for r in result)),reused_Rosina_BAMs=sum(r['reuse_Rosina_BAM'] for r in result),paired_local_FASTQ_candidates=sum(r['paired_local_FASTQ_candidates'] and not r['reuse_Rosina_BAM'] for r in result),completed_downloads=sum(r['existing_download_complete'] for r in result),warning='Capture/RAD labels do not establish ventral coverage; inspect targets or observed reads. No tasks excluded or jobs changed.')
 (out/'complete.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
 print('\nGROUPS BY SEQUENCING CATEGORY AND TAXON')
 table=pd.DataFrame(result).groupby(['sequencing_category','taxon']).size()
 print(table.to_string())
if __name__=='__main__':main()
