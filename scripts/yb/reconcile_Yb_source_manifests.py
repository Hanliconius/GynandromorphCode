#!/usr/bin/env python3
"""Read-only manifest provenance reconciliation; no downloads or new analysis jobs.
Follow explicit config paths, including symlink targets outside the project root.
Write a resource ledger and one review decision per pending ventral individual.
"""
import argparse,csv,gzip,json,os,re
from collections import defaultdict,deque
from datetime import datetime,timezone
from pathlib import Path
import pandas as pd
import pysam
CHROMS={'NC_140403.1':('HmRos',(1633041,1635646),(1658079,1659072)),'Hmel215003o':('Hmel2.5',(1403500,1405500),(1429500,1430500)),'Hmr15':('legacy_alias_unverified',(1633041,1635646),(1658079,1659072))}
SKIP={'.git','logs','raw','fastq','bam','bams','alignment','alignments','tmp','windows','code','toolenv','__pycache__'}
def json_refs(obj,field='',context=None):
 context=dict(context or {})
 if isinstance(obj,dict):
  for key in ('sample','sample_id','accession','library_id','library','role'):
   if isinstance(obj.get(key),str):context[key]=obj[key]
  for key,value in obj.items():yield from json_refs(value,field+'.'+str(key),context)
 elif isinstance(obj,list):
  for i,value in enumerate(obj):yield from json_refs(value,field+f'[{i}]',context)
 elif isinstance(obj,str) and obj.startswith('/') and '\n' not in obj:
  yield obj,field,context

def resource_type(path,field):
 name=path.name.lower();f=field.lower()
 if name.endswith(('.bam','.cram')):return 'alignment'
 if name.endswith(('.vcf','.vcf.gz','.bcf')):return 'variant_calls'
 if name.endswith(('.fastq','.fastq.gz','.fq','.fq.gz')):return 'FASTQ'
 if name.endswith(('.bw','.bigwig','.bedgraph','.bedgraph.gz')):return 'depth_profile'
 if name.endswith(('.tsv','.tsv.gz','.csv','.csv.gz','.npz','.parquet')):
  if re.search(r'count|allele|genotype|calls',name+' '+f):return 'allele_table_candidate'
  if re.search(r'coverage|depth|profile|bigwig',name+' '+f):return 'depth_profile'
  if re.search(r'manifest|sample|inventory|metadata|input|reference|source',name):return 'manifest'
 if name.endswith('.json'):return 'configuration'
 return 'other'

def read_table_paths(path):
 opener=gzip.open if path.suffix=='.gz' else open
 with opener(path,'rt',errors='replace') as h:
  rows=csv.DictReader(h,delimiter=',' if '.csv' in path.name else '\t')
  for i,row in enumerate(rows):
   # Only small manifest tables are passed here, not per-site data.
   context={k:str(row[k]) for k in ('sample','sample_id','accession','library_id','role') if row.get(k)}
   for k,v in row.items():
    if isinstance(v,str) and v.startswith('/') and '\n' not in v:yield v,f'row[{i}].{k}',context

def inspect(path,typ):
 result=dict(reference='',contig='',samples='',d_reads='',v_reads='',indexed='',evidence='exists_nonempty')
 if typ=='alignment':
  with pysam.AlignmentFile(str(path),'rc' if path.suffix=='.cram' else 'rb') as bam:
   sm=sorted({r.get('SM') for r in bam.header.to_dict().get('RG',[]) if r.get('SM')});result['samples']=';'.join(sm);result['indexed']=bam.has_index()
   dic=dict(zip(bam.references,bam.lengths));found=[]
   for chrom,(ref,d,v) in CHROMS.items():
    if chrom not in dic:continue
    nr=[]
    for lo,hi in (d,v):nr.append(bam.count(chrom,lo-1,hi,read_callback='all') if bam.has_index() and dic[chrom]>=hi else '')
    found.append(dict(reference=ref,contig=chrom,d_reads=nr[0],v_reads=nr[1]))
   if len(found)==1:result.update(found[0]);result['evidence']='BAM RG/header and indexed interval read-presence check; quality/reference compatibility still requires validation'
   elif not found:result['evidence']='BAM exists; contig naming/reference unresolved'
 if typ=='variant_calls':
  with pysam.VariantFile(str(path)) as v:
   result['samples']=';'.join(v.header.samples);result['contig']=';'.join(v.header.contigs)
   recognized=set(v.header.contigs)&set(CHROMS)
   if len(recognized)==1:
    chrom=next(iter(recognized));ref,d,vspan=CHROMS[chrom];result['reference']=ref
    indexed=any(Path(str(path)+x).is_file() for x in ('.tbi','.csi'));result['indexed']=indexed
    if indexed:
     result['d_records']=sum(1 for _ in v.fetch(chrom,d[0]-1,d[1]));result['v_records']=sum(1 for _ in v.fetch(chrom,vspan[0]-1,vspan[1]))
    elif path.stat().st_size<100_000_000:
     nd=nv=0
     for rec in v:
      if rec.contig!=chrom:continue
      nd+=d[0]<=rec.pos<=d[1];nv+=vspan[0]<=rec.pos<=vspan[1]
     result['d_records']=nd;result['v_records']=nv
  result['evidence']='VCF header and interval record availability where indexed or small; per-sample allele/quality compatibility requires validation'
 if typ in ('allele_table_candidate','depth_profile') and path.name.endswith(('.tsv','.tsv.gz','.csv','.csv.gz')):
  opener=gzip.open if path.suffix=='.gz' else open
  with opener(path,'rt',errors='replace') as h:result['table_header']=h.readline().strip()[:3000]
  result['evidence']='Existing table header; coordinates/filter compatibility not yet established'
 return result

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();cfg=json.loads((a.run/'run.json').read_text());root=Path(cfg['root']);out=a.output;out.mkdir(parents=True,exist_ok=False)
 queue=deque();queued=set();errors=[];edges=[];resources={};sources=0
 def enqueue(path):
  path=Path(path)
  if not path.is_file():return
  real=path.resolve()
  if real not in queued:queue.append(path);queued.add(real)
 # Start with all explicit run/config/manifest files, excluding newly generated
 # expansion inputs as evidence of what existed before this submission.
 for folder in (root/'analysis',root/'data',root/'external'):
  if not folder.exists():continue
  for parent,dirs,names in os.walk(folder):
   dirs[:]=[d for d in dirs if d not in SKIP and not d.startswith('figure_code_') and d not in {'yb_ventral_population_expansion_v1','yb_matched_cohort_v1'}]
   for name in names:
    path=Path(parent)/name
    if name in {'run.json','config.json','audit.json'} or (re.search(r'manifest|inventory|sources|inputs',name,re.I) and name.endswith(('.tsv','.csv','.json'))):
     if path.stat().st_size<10_000_000:enqueue(path)
 # Exact historic configuration documented by the repository; its nested
 # sources may point outside every directory scanned above.
 enqueue(root/'analysis/yb_all_samples_coverage/20261014T194126Z/run.json')
 while queue:
  source=queue.popleft();sources+=1
  try:
   if source.suffix=='.json':refs=list(json_refs(json.loads(source.read_text())))
   else:refs=list(read_table_paths(source))
  except Exception as e:errors.append(dict(path=str(source),error=str(e)));continue
  for value,field,context in refs:
   path=Path(value);typ=resource_type(path,field)
   # Preserve missing explicit paths and directory provenance rather than
   # translating either into a download decision.
   exists=path.is_file();resolved=str(path.resolve());key=(resolved,typ)
   edges.append(dict(source=str(source),field=field,path=value,resolved_path=resolved,context_sample=context.get('sample') or context.get('sample_id') or context.get('accession') or context.get('library_id',''),role=context.get('role',''),resource_type=typ,exists=exists))
   if path.is_dir():
    if typ=='other' and any(k in field.lower() for k in ('source','run','input','coverage','profile','alignment','bam','raw','fastq')):
     for name in ('run.json','config.json','audit.json'):enqueue(path/name)
    continue
   if typ=='configuration' and exists and path.stat().st_size<10_000_000:enqueue(path)
   elif typ=='manifest' and exists and path.stat().st_size<10_000_000:enqueue(path)
   if typ in {'configuration','manifest','other'}:continue
   if key not in resources:
    row=dict(path=value,resolved_path=resolved,resource_type=typ,exists=exists,bytes=path.stat().st_size if exists else 0,source_config=str(source),source_config_field=field,context_samples=set(),role=context.get('role',''),reference='',contig='',samples='',d_reads='',v_reads='',indexed='',evidence='missing_explicit_manifest_path' if not exists else 'exists_nonempty')
    if exists and row['bytes']:
     try:row.update(inspect(path,typ))
     except Exception as e:row['evidence']='inspection_failed';errors.append(dict(path=value,error=str(e)))
    resources[key]=row
   context_sample=context.get('sample') or context.get('sample_id') or context.get('accession') or context.get('library_id')
   if context_sample:resources[key]['context_samples'].add(context_sample)
 # Include the first audit's native artifacts, retaining its limitations.
 old=root/'analysis/yb_existing_resources_v1/bam_candidates.tsv'
 if old.is_file():
  for b in pd.read_csv(old,sep='\t').fillna('').to_dict('records'):
   key=(b['resolved_path'],'alignment')
   if key in resources:continue
   path=Path(b['resolved_path']);r=dict(path=b['path'],resolved_path=b['resolved_path'],resource_type='alignment',exists=path.is_file(),bytes=path.stat().st_size if path.is_file() else 0,source_config=str(old),source_config_field='resolved_path',context_samples={b['sample']})
   try:r.update(inspect(path,'alignment')) if path.is_file() else r.update(evidence='missing_explicit_audit_path')
   except Exception as e:r['evidence']='inspection_failed';errors.append(dict(path=str(path),error=str(e)))
   resources[key]=r
 now=datetime.now(timezone.utc).isoformat();rows=[];by_sample=defaultdict(list)
 for i,r in enumerate(resources.values()):
  r['resource_id']=f'YBR{i:06d}';r['checked_at_utc']=now;r['context_samples']=';'.join(sorted(r['context_samples']));rows.append(r)
  for sample in set(r.get('samples','').split(';'))|set(r['context_samples'].split(';')):
   if sample:by_sample[sample].append(r)
 pd.DataFrame(rows).to_csv(out/'resource_manifest.tsv',sep='\t',index=False);pd.DataFrame(edges).to_csv(out/'provenance_edges.tsv',sep='\t',index=False);pd.DataFrame(errors,columns=['path','error']).to_csv(out/'read_errors.tsv',sep='\t',index=False)
 decisions=[]
 for t in cfg['new_tasks']:
  ss=[r for r in by_sample[t['sample']] if r.get('exists') and r.get('bytes',0)>0]
  bam=[r for r in ss if r['resource_type']=='alignment' and r.get('reference')=='HmRos' and isinstance(r.get('v_reads'),int) and r['v_reads']>0]
  variants=[r for r in ss if r['resource_type']=='variant_calls' and r.get('v_records',0)>0]
  counts=[r for r in ss if r['resource_type']=='allele_table_candidate'];depth=[r for r in ss if r['resource_type']=='depth_profile'];fq=[r for r in ss if r['resource_type']=='FASTQ']
  status='validate_existing_Rosina_alignment' if bam else 'inspect_existing_variant_calls' if variants else 'inspect_existing_allele_tables' if counts else 'resolve_coverage_source' if depth else 'validate_existing_read_pairs' if fq else 'unresolved_not_evidence_of_absence'
  decisions.append(dict(sample=t['sample'],taxon=t['taxon'],interval='v',action=status,Rosina_BAMs_with_ventral_reads=len(bam),variant_calls_with_ventral_records=len(variants),allele_table_candidates=len(counts),depth_profiles=len(depth),FASTQ_resources=len(fq),resource_ids=';'.join(r['resource_id'] for r in ss),previous_manifest_BAM=t['reuse_bam']))
 pd.DataFrame(decisions).to_csv(out/'reuse_decisions.tsv',sep='\t',index=False)
 summary=dict(status='RECONCILIATION_REQUIRES_REVIEW',source_configs_read=sources,resources=len(rows),explicit_missing_paths=sum(not r.get('exists') for r in rows),read_errors=len(errors),targets=len(decisions),decisions={str(k):int(v) for k,v in pd.Series([r['action'] for r in decisions]).value_counts().items()},jobs_changed=False,downloads=0,alignments=0,calling=0,warning='Unresolved entries do not authorize downloads; candidate resources still require reference/allele/filter compatibility checks')
 (out/'complete.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
