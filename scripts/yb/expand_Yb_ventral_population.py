#!/usr/bin/env python3
"""Add ventral data to the full existing dorsal cohort; never shrink it.
Reuse Rosina BAMs/FASTQs, else download paired reads. Align against the full Rosina
assembly, duplicate-mark before retaining chromosome 15, call frozen Yb markers.
"""
import argparse,csv,hashlib,importlib.util,json,os,re,shlex,shutil,subprocess,sys,tempfile,time,urllib.parse,urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import pandas as pd
import pysam
BASE=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
PARENT=BASE/'analysis/yb_expanded_existing_controls_v2'
CHROM='NC_140403.1'
def sha(p,algorithm='sha256'):
 h=hashlib.new(algorithm)
 with Path(p).open('rb') as f:
  for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
 return h.hexdigest()
def save(p,x):
 tmp=Path(str(p)+'.tmp');tmp.write_text(json.dumps(x,indent=2)+'\n');tmp.replace(p)
def run(cmd,**kwargs):return subprocess.run([str(x) for x in cmd],check=True,**kwargs)
def tool_environment(cfg):
 for name,path in cfg.get('tools',{}).items():
  if not Path(path).is_file():raise RuntimeError('Required executable disappeared: '+path)
 dirs=list(dict.fromkeys(str(Path(p).parent) for p in cfg.get('tools',{}).values()))
 os.environ['PATH']=os.pathsep.join(dirs+[os.environ.get('PATH','')])

def engine(cfg):
 p=Path(cfg['output'])/'engine.py'
 if sha(p)!=cfg['engine_sha']:raise RuntimeError('Frozen calling/plotting engine changed')
 spec=importlib.util.spec_from_file_location('yb_ventral_engine',p);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module

def ena(sample):
 params=dict(accession=sample,result='read_run',fields='run_accession,sample_accession,library_layout,fastq_ftp,fastq_md5,fastq_bytes',format='tsv')
 url='https://www.ebi.ac.uk/ena/portal/api/filereport?'+urllib.parse.urlencode(params)
 last=None
 for attempt in range(4):
  try:
   with urllib.request.urlopen(url,timeout=90) as h:raw=h.read().decode()
   rows=list(csv.DictReader(raw.splitlines(),delimiter='\t'));rows=[r for r in rows if r.get('run_accession')==sample and r.get('library_layout')=='PAIRED']
   if len(rows)!=1:raise RuntimeError('No unique paired ENA run '+sample)
   row=rows[0];fields=[row.get(k,'').split(';') for k in ('fastq_ftp','fastq_md5','fastq_bytes')]
   if len({len(x) for x in fields})!=1:raise RuntimeError('Inconsistent ENA mate metadata')
   mates={}
   for url,md5,n in zip(*fields):
    match=re.search(r'_([12])\.(?:fastq|fq)\.gz$',url)
    if not match:continue
    i=int(match[1]);url='https://'+url.removeprefix('ftp://').removeprefix('https://')
    if i in mates or not url.startswith('https://ftp.sra.ebi.ac.uk/') or not re.fullmatch('[a-fA-F0-9]{32}',md5) or not n.isdigit():raise RuntimeError('Invalid/ambiguous ENA pair')
    mates[i]=dict(url=url,md5=md5.lower(),bytes=int(n))
   if set(mates)!={1,2}:raise RuntimeError('Expected exactly two ENA mate files')
   return dict(run=sample,BioSample=row.get('sample_accession',''),mates=mates,raw_response=raw)
  except Exception as e:
   last=e
   if attempt<3:time.sleep(2*(attempt+1))
 raise RuntimeError(str(last))

def existing_fastqs(root,targets):
 found={s:{1:[],2:[]} for s in targets}
 for folder in (root/'analysis',root/'data',root/'external'):
  if not folder.exists():continue
  for parent,dirs,names in os.walk(folder):
   dirs[:]=[d for d in dirs if d not in {'.git','logs','tmp','alignment','bam','bams','windows','code','toolenv','node_modules'} and not d.startswith('figure_code_')]
   for name in names:
    match=re.fullmatch(r'((?:ERR|SRR|DRR)\d+)[_.]([12])\.(?:fastq|fq)\.gz',name)
    if not match or match[1] not in found:continue
    path=Path(parent)/name
    if path.is_file() and path.stat().st_size:found[match[1]][int(match[2])].append(str(path.resolve()))
 return found

def compatible_bam(path,sample,reference):
 try:
  with pysam.FastaFile(str(reference)) as fa:pairs=dict(zip(fa.references,fa.lengths))
  with pysam.AlignmentFile(str(path),'rb') as bam:
   dic=dict(zip(bam.references,bam.lengths));sm={r.get('SM') for r in bam.header.to_dict().get('RG',[]) if r.get('SM')}
   return bam.has_index() and sm=={sample} and CHROM in dic and all(pairs.get(c)==n for c,n in dic.items())
 except (OSError,ValueError):return False

def prepare(a):
 parent=a.parent;out=a.output
 status=json.loads((parent/'complete.json').read_text());old=json.loads((parent/'run.json').read_text())
 if status.get('status')!='COMPLETE' or not status.get('focal_calls_unchanged'):raise RuntimeError('Parent expansion is not complete')
 members=pd.read_csv(parent/'control_membership.tsv',sep='\t');active=members[members['included'].isin([True,'True',1]) & members.reference.eq('HmRos')]
 dorsal=active[active.interval.eq('d')].drop_duplicates('sample');ventral=set(active.loc[active.interval.eq('v'),'sample'])
 missing=dorsal[~dorsal['sample'].isin(ventral)].copy();known=missing[pd.to_numeric(missing['v'],errors='coerce').isin([0,1])].copy()
 unknown=missing[~missing['sample'].isin(known['sample'])]
 if known.empty:raise RuntimeError('No phenotype-annotated dorsal individuals missing ventral calls')
 out.mkdir(parents=True,exist_ok=False)
 for name in ('logs','raw','bam','calls','ref','qc'): (out/name).mkdir()
 shutil.copytree(parent/'frozen',out/'frozen');shutil.copy2(parent/'legacy_reference_check.json',out/'legacy_reference_check.json')
 renderer=parent/'expand_Yb_existing_controls.py'
 if sha(renderer)!=old['script_sha']:raise RuntimeError('Parent engine hash changed')
 shutil.copy2(renderer,out/'engine.py');shutil.copy2(__file__,out/'expand_Yb_ventral_population.py')
 reference=Path(old['references']['HmRos']);shutil.copy2(str(reference)+'.fai',out/'ref/genome.fa.fai');(out/'ref/genome.fa').symlink_to(reference.resolve())
 for suffix in ('.amb','.ann','.bwt','.pac','.sa','.alt'):
  p=Path(str(reference)+suffix)
  if p.is_file():(out/'ref'/('genome.fa'+suffix)).symlink_to(p.resolve())
 with pysam.FastaFile(str(reference)) as fa:
  if len(fa.references)<20:raise RuntimeError('Alignment reference is not the full Rosina assembly')
  length=fa.get_reference_length(CHROM)
 (out/'ref/chromosome15.bed').write_text(f'{CHROM}\t0\t{length}\n')
 for f in parent.glob('HmRos.*'):
  if f.is_file():shutil.copy2(f,out/f.name)
 for task in old['tasks']:
  for suffix in ('.vcf.gz','.vcf.gz.csi','.vcf.gz.tbi','.complete.json'):
   p=parent/'calls'/f'{task["task"]}{suffix}'
   if p.is_file():(out/'calls'/p.name).symlink_to(p.resolve())
 root=Path(old['root'])
 # Also admit ventral-annotated legacy individuals whose dorsal phenotype is unknown.
 helper=engine(dict(output=str(out),engine_sha=sha(out/'engine.py')));getmeta=helper.metadata(root,Path(old['source']))
 with pysam.VariantFile(old['legacy']) as vf:legacy_ids=set(vf.header.samples)
 extra=[]
 for sample in sorted(legacy_ids-set(known['sample'])-ventral):
  row=getmeta(sample)
  if not row['exclusion'] and row['v'] is not None:extra.append(row)
 if extra:known=pd.concat([known,pd.DataFrame(extra)],ignore_index=True)
 fastqs=existing_fastqs(root,set(known['sample']))
 audit=pd.read_csv(root/'analysis/yb_existing_resources_v1/bam_candidates.tsv',sep='\t').fillna('');tasks=[];base_index=len(old['tasks'])
 for i,row in enumerate(known.sort_values('sample').to_dict('records')):
  sample=row['sample'];candidates=audit[(audit['sample']==sample)&audit.reference.eq('HmRos')]
  paths=sorted(set(candidates.resolved_path),key=lambda x:('/tmp/' in x,not x.endswith('.dedup.bam'),x));bam=''
  for p in paths:
   if '/tmp/' not in p and compatible_bam(p,sample,reference):bam=p;break
  tasks.append(dict(index=i,call_task=base_index+i,sample=sample,taxon=row['taxon'],d=int(row['d']) if pd.notna(row['d']) else None,v=int(row['v']),reuse_bam=bam,fastq_candidates=fastqs[sample]))
 tool_dirs=[reference.parent.parent/'toolenv/bin']+list((root/'analysis').glob('*/*/toolenv/bin'))
 tools={}
 for name in ('bwa','samtools','bcftools','curl'):
  found=shutil.which(name)
  if not found:found=next((str(p/name) for p in tool_dirs if (p/name).is_file()),None)
  if not found:raise RuntimeError('No existing executable found for '+name)
  tools[name]=found
 cfg=dict(old);cfg.update(output=str(out),parent=str(parent),engine_sha=sha(out/'engine.py'),script_sha=sha(out/'expand_Yb_ventral_population.py'),new_tasks=tasks,original_tasks=old['tasks'],reference=str(out/'ref/genome.fa'),download_concurrency=a.download_concurrency,align_concurrency=a.align_concurrency,tools=tools)
 cfg['tasks']=list(old['tasks'])+[dict(task=t['call_task'],sample=t['sample'],reference='HmRos',bam=t['reuse_bam'] or str(out/'bam'/f'{t["sample"]}.chr15.dedup.bam')) for t in tasks]
 save(out/'run.json',cfg);dorsal.to_csv(out/'retained_dorsal_individuals.tsv',sep='\t',index=False);unknown.to_csv(out/'ventral_unknown_phenotype.tsv',sep='\t',index=False)
 pd.DataFrame([{k:v for k,v in t.items() if k!='fastq_candidates'} for t in tasks]).to_csv(out/'ventral_expansion_tasks.tsv',sep='\t',index=False)
 scope=dict(legacy_ventral_only_phenotype_candidates=len(extra),dorsal_individuals_preserved=len(dorsal),existing_ventral_individuals=len(ventral),ventral_targets=len(tasks),unknown_ventral_phenotype=len(unknown),reused_Rosina_BAMs=sum(bool(t['reuse_bam']) for t in tasks),paired_local_FASTQ_candidates=sum(bool(t['fastq_candidates'][1]) and bool(t['fastq_candidates'][2]) for t in tasks if not t['reuse_bam']),alignment='Full Rosina genome; duplicate-marked chromosome 15 retained',downloads='Only missing or checksum-invalid mates',PK_and_laboratory_crosses='Excluded through parent cohort')
 save(out/'scope.json',scope);print(json.dumps(scope,indent=2))
 if a.submit:submit(cfg)

def submit(cfg):
 out=Path(cfg['output']);n=len(cfg['new_tasks']);script=out/'expand_Yb_ventral_population.py'
 def sb(stage,partition,cpus,mem,limit,cap=None,deps=None):
  cmd=['sbatch','--parsable','--partition='+partition,'--cpus-per-task='+str(cpus),'--mem='+mem,'--time='+limit,'--job-name=Ybv_'+stage,'--output='+str(out/'logs'/f'{stage}.%A_%a.out'),'--error='+str(out/'logs'/f'{stage}.%A_%a.err')]
  if cap:cmd+=['--array='+f'0-{n-1}%{cap}']
  if deps:cmd+=['--dependency='+deps,'--kill-on-invalid-dep=yes']
  wrapped=shlex.join(['micromamba','run','-n','gynandrocalling','python',str(script),'--output',str(out),'--stage',stage]);cmd+=['--wrap='+wrapped]
  return subprocess.check_output(cmd,text=True).strip().split(';')[0]
 setup=sb('setup','cpu',4,'12G','02:00:00')
 download=sb('download','cpu',2,'4G','1-00:00:00',cfg['download_concurrency'])
 align=sb('align','cpu',16,'32G','1-00:00:00',cfg['align_concurrency'],f'afterok:{setup},aftercorr:{download}')
 call=sb('call','nano',1,'8G','00:30:00',48,f'aftercorr:{align}')
 report=sb('report','nano',1,'8G','00:30:00',deps=f'afterany:{call}')
 jobs=dict(setup=setup,download=download,align=align,call=call,report=report);save(out/'jobs.json',jobs)
 for k,v in jobs.items():print(k.upper()+'\t'+v)
 print('PDF\t'+str(out/'Yb_expanded_existing_controls_Rosina.pdf'))

def setup(cfg):
 out=Path(cfg['output']);ref=out/'ref/genome.fa'
 for tool in ('bwa','samtools','bcftools','curl'):
  if not shutil.which(tool):raise RuntimeError('Missing executable in gynandrocalling: '+tool)
 if not all(Path(str(ref)+suffix).is_file() for suffix in ('.amb','.ann','.bwt','.pac','.sa')):run(['bwa','index',ref])
 save(out/'setup.complete.json',dict(full_genome_reference=str(ref.resolve()),indexed=True,versions={t:shutil.which(t) for t in ('bwa','samtools','bcftools','curl')}))

def download(cfg,t):
 if t['reuse_bam']:print('BAM_REUSED\t'+t['sample']);return
 out=Path(cfg['output']);d=out/'raw'/t['sample'];d.mkdir(exist_ok=True)
 metadata=ena(t['sample']);save(d/'ena.json',metadata);verified={}
 def verified_file(p,vals):
  if not p.is_file() or p.stat().st_size!=vals['bytes']:return False
  stat=p.stat();key=(str(p.resolve()),stat.st_size,stat.st_mtime_ns,vals['md5'])
  if key not in verified:verified[key]=sha(p,'md5')==vals['md5']
  return verified[key]
 def mate(i):
  vals=metadata['mates'][i];target=d/f'{t["sample"]}_{i}.fastq.gz'
  def valid(p):return verified_file(p,vals)
  if not valid(target):
   for path in t['fastq_candidates'].get(str(i),t['fastq_candidates'].get(i,[])):
    source=Path(path)
    if valid(source):
     target.unlink(missing_ok=True);target.symlink_to(source.resolve());break
  if not valid(target):
   tmp=Path(str(target)+'.partial')
   cmd=['curl','--fail','--location','--retry','12','--retry-all-errors','--retry-delay','5','--connect-timeout','30','--speed-time','180','--speed-limit','1024','--continue-at','-','--silent','--show-error','--output',str(tmp),vals['url']]
   run(cmd)
   if not valid(tmp):
    tmp.unlink(missing_ok=True)
    raise RuntimeError(f'FASTQ checksum/size mismatch: {t["sample"]} mate {i}; invalid completed download discarded')
   target.unlink(missing_ok=True);tmp.replace(target)
  return dict(mate=i,bytes=vals['bytes'],md5=vals['md5'],path=str(target.resolve()))
 with ThreadPoolExecutor(max_workers=2) as pool:mates=list(pool.map(mate,(1,2)))
 save(d/'complete.json',dict(sample=t['sample'],BioSample=metadata['BioSample'],mates=mates));print('READS_READY\t'+t['sample'])

def align(cfg,t):
 out=Path(cfg['output']);sample=t['sample']
 if t['reuse_bam']:
  if not compatible_bam(t['reuse_bam'],sample,cfg['reference']):raise RuntimeError('Reused BAM failed reference/identity check')
  print('BAM_REUSED\t'+sample);return
 dst=out/'bam'/f'{sample}.chr15.dedup.bam';done=out/'bam'/f'{sample}.complete.json'
 if done.is_file() and compatible_bam(dst,sample,cfg['reference']):print('EXISTING_ALIGNMENT\t'+sample);return
 d=out/'raw'/sample;readmeta=json.loads((d/'complete.json').read_text())
 for r in readmeta['mates']:
  p=d/f'{sample}_{r["mate"]}.fastq.gz'
  if not p.is_file() or p.stat().st_size!=r['bytes']:raise RuntimeError('Validated FASTQ missing/changed')
 rg='@RG\\tID:'+sample+'\\tSM:'+sample+'\\tPL:ILLUMINA'
 with tempfile.TemporaryDirectory(prefix=sample+'.',dir=out/'bam') as temp:
  tmp=Path(temp);dsttmp=tmp/'chr15.bam'
  commands=[['bwa','mem','-M','-t','8','-R',rg,cfg['reference'],str(d/f'{sample}_1.fastq.gz'),str(d/f'{sample}_2.fastq.gz')],['samtools','sort','-n','-@','2','-m','1G','-T',str(tmp/'name'),'-O','bam','-'],['samtools','fixmate','-m','-u','-','-'],['samtools','sort','-@','2','-m','1G','-T',str(tmp/'coord'),'-O','bam','-'],['samtools','markdup','-s','-u','-','-'],['samtools','view','-b','-L',str(out/'ref/chromosome15.bed'),'-o',str(dsttmp),'-']]
  processes=[]
  try:
   for command in commands:
    prev=processes[-1] if processes else None;p=subprocess.Popen(command,stdin=prev.stdout if prev else None,stdout=subprocess.PIPE if len(processes)<len(commands)-1 else None);processes.append(p)
    if prev:prev.stdout.close()
   codes=[p.wait() for p in reversed(processes)][::-1]
   if any(codes):raise RuntimeError('Alignment pipeline failed: '+str(codes))
  finally:
   for p in processes:
    if p.poll() is None:p.kill();p.wait()
  run(['samtools','quickcheck',dsttmp]);run(['samtools','index',dsttmp]);dsttmp.replace(dst);Path(str(dsttmp)+'.bai').replace(str(dst)+'.bai')
 if not compatible_bam(dst,sample,cfg['reference']):raise RuntimeError('Produced BAM failed reference/sample validation')
 save(done,dict(sample=sample,full_genome_mapping=True,duplicates_marked_before_chromosome_filter=True,retained_chromosome=CHROM,commands=commands));print('BAM_READY\t'+sample)

def call(cfg,t):
 out=Path(cfg['output']);i=t['call_task'];marker=out/'calls'/f'{i}.complete.json';vcf=out/'calls'/f'{i}.vcf.gz'
 if marker.is_file() and vcf.is_file() and json.loads(marker.read_text()).get('vcf_sha256')==sha(vcf):print('CALLS_REUSED\t'+t['sample']);return
 engine(cfg).call(cfg,i)

def report(cfg):
 out=Path(cfg['output']);parent=Path(cfg['parent']);selected=pd.read_csv(parent/'selected_sources.tsv',sep='\t');original_members=pd.read_csv(parent/'control_membership.tsv',sep='\t')
 failures=[];ready=[]
 for t in cfg['new_tasks']:
  marker=out/'calls'/f'{t["call_task"]}.complete.json';vcf=out/'calls'/f'{t["call_task"]}.vcf.gz'
  if not marker.is_file() or not vcf.is_file():failures.append(dict(sample=t['sample'],reason='call_stage_not_completed',index=t['index']));continue
  if json.loads(marker.read_text()).get('vcf_sha256')!=sha(vcf):failures.append(dict(sample=t['sample'],reason='call_hash_mismatch',index=t['index']));continue
  ready.append(t)
  # Replace an uncallable previous ventral source, if present; preserve all dorsal rows.
  selected=selected[~(selected['sample'].eq(t['sample'])&selected.reference.eq('HmRos')&selected.interval.eq('v'))]
  selected=pd.concat([selected,pd.DataFrame([dict(sample=t['sample'],taxon=t['taxon'],d=t['d'],v=t['v'],reference='HmRos',interval='v',source='forced_BAM_call',path=str(vcf))])],ignore_index=True)
 selected.to_csv(out/'selected_sources.tsv',sep='\t',index=False)
 selectedtasks={t['call_task'] for t in ready};local=dict(cfg);local['tasks']=cfg['original_tasks']+[t for t in cfg['tasks'] if t['task'] in selectedtasks]
 module=engine(cfg);module.report(local)
 result=json.loads((out/'complete.json').read_text());(out/'complete.json').unlink()
 new=pd.read_csv(out/'control_membership.tsv',sep='\t')
 def dorsal_members(frame):
  return set(frame.loc[frame.reference.eq('HmRos')&frame.interval.eq('d')&frame['included'].isin([True,'True',1]),'sample'])
 if dorsal_members(new)!=dorsal_members(original_members):raise RuntimeError('Dorsal cohort changed; expansion rejected')
 # Hard guarantee: dorsal plotted frequencies and focal rows are unchanged.
 import numpy as np
 original=pd.read_csv(parent/'HmRos_Ybd_plotted_values.tsv',sep='\t',index_col=0);current=pd.read_csv(out/'HmRos_Ybd_plotted_values.tsv',sep='\t',index_col=0)
 np.testing.assert_equal(original.to_numpy(),current.to_numpy());assert list(original.index)==list(current.index)
 # Independent identity/batch check: legacy dorsal genotypes versus fresh calls, without changing dorsal values.
 table=pd.read_csv(out/'frozen/HmRos_Ybd_unchanged_focal_calls.tsv.gz',sep='\t').sort_values('pos1');baseline_values={};concordance=[]
 old_dorsal=selected[selected.reference.eq('HmRos')&selected.interval.eq('d')&selected['sample'].isin([t['sample'] for t in ready])]
 for (path,kind),group in old_dorsal.groupby(['path','source']):
  vv,_=module.extract(path,group['sample'].tolist(),table,legacy=kind=='legacy_VCF',allow_pl=True);baseline_values.update(vv)
 for t in ready:
  if t['sample'] not in baseline_values:continue
  vv,_=module.extract(out/'calls'/f'{t["call_task"]}.vcf.gz',[t['sample']],table,allow_pl=True)
  old=baseline_values[t['sample']];fresh=vv[t['sample']];both=np.isfinite(old)&np.isfinite(fresh);n=int(both.sum());discord=int((old[both]!=fresh[both]).sum())
  concordance.append(dict(sample=t['sample'],shared_callable_dorsal_markers=n,discordant_genotypes=discord,concordance=(n-discord)/n if n else None,review_flag=n>=10 and discord/n>.1))
 pd.DataFrame(concordance).to_csv(out/'dorsal_legacy_recall_concordance.tsv',sep='\t',index=False,na_rep='NA')
 flags=sum(r['review_flag'] for r in concordance)
 pd.DataFrame(failures,columns=['sample','reason','index']).to_csv(out/'failed_expansion_tasks.tsv',sep='\t',index=False)
 result.update(status='COMPLETE_WITH_FAILED_TASKS' if failures else 'COMPLETE_REQUIRES_CONCORDANCE_REVIEW' if flags else 'COMPLETE',concordance_review_flags=flags,parent=str(parent),dorsal_cohort_preserved=True,dorsal_frequencies_identical=True,ventral_tasks_planned=len(cfg['new_tasks']),ventral_tasks_completed=len(ready),failed_tasks=failures,scope=json.loads((out/'scope.json').read_text()))
 save(out/'complete.json',result);print('EXPANSION_STATUS\t'+result['status']);print('FAILED_TASKS\t'+str(len(failures)))

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--parent',type=Path,default=PARENT);p.add_argument('--output',type=Path,required=True);p.add_argument('--submit',action='store_true');p.add_argument('--retry',action='store_true');p.add_argument('--download-concurrency',type=int,default=48);p.add_argument('--align-concurrency',type=int,default=24);p.add_argument('--stage',choices=['prepare','setup','download','align','call','report'],default='prepare');p.add_argument('--task',type=int);a=p.parse_args()
 if a.stage=='prepare' and not a.retry:prepare(a);return
 cfg=json.loads((a.output/'run.json').read_text());tool_environment(cfg)
 if sha(__file__)!=cfg['script_sha']:raise RuntimeError('Frozen workflow changed')
 if a.retry:submit(cfg);return
 if a.stage in ('setup','report'):globals()[a.stage](cfg)
 else:
  i=a.task if a.task is not None else int(os.environ['SLURM_ARRAY_TASK_ID']);globals()[a.stage](cfg,cfg['new_tasks'][i])
if __name__=='__main__':main()
