#!/usr/bin/env python3
"""Recover downstream Yb stages from validated reads already in hand.
No downloads. Isolated output, individual phenotype ledger, explicit completion gate.
"""
import argparse,csv,hashlib,importlib.util,json,os,shlex,shutil,subprocess,sys
from pathlib import Path
from collections import defaultdict
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
def read(p):
 with Path(p).open() as h:return list(csv.DictReader(h,delimiter='\t'))
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def digest(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def tab(p,rows):
 if not rows:p.write_text('status\nNO_ROWS\n');return
 cols=list(dict.fromkeys(k for r in rows for k in r))
 with p.open('w',newline='') as h:
  w=csv.DictWriter(h,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(rows)
def module(cfg):
 p=Path(cfg['output'])/'code/expansion.py'
 if digest(p)!=cfg['expansion_code_sha']:raise RuntimeError('Expansion code changed')
 spec=importlib.util.spec_from_file_location('existing_expansion',p);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);m.tool_environment(cfg);return m
def ready_reads(source,t):
 p=source/'raw'/t['sample']/'complete.json'
 if not p.is_file():return False,'no_download_completion_marker'
 try:
  x=json.loads(p.read_text())
  if x['sample']!=t['sample'] or len(x['mates'])!=2:return False,'invalid_completion_identity'
  for r in x['mates']:
   f=source/'raw'/t['sample']/f'{t["sample"]}_{r["mate"]}.fastq.gz'
   if not f.is_file() or f.stat().st_size!=r['bytes'] or not r.get('md5'):return False,'completed_read_missing_or_size_changed'
   if f.stat().st_mtime>p.stat().st_mtime+1:return False,'reads_modified_after_checksum_validation'
  return True,'existing_download_checksum_completion_and_unchanged_sizes'
 except Exception as e:return False,str(e)
def prepare(a):
 source=a.source;out=a.output;cfg=json.loads((source/'run.json').read_text());focus=sorted((source/'capture_focus').glob('*/run.json'))
 scope=json.loads(focus[-1].read_text()) if focus else cfg
 out.mkdir(parents=True,exist_ok=False)
 for n in ('code','logs','tasks','bam','calls','raw'):(out/n).mkdir()
 shutil.copy2(__file__,out/'code/recover_Ybv_from_completed_reads.py')
 shutil.copy2(source/'expand_Yb_ventral_population.py',out/'code/expansion.py')
 shutil.copy2(source/'engine.py',out/'engine.py');shutil.copytree(source/'frozen',out/'frozen')
 shutil.copy2(source/'legacy_reference_check.json',out/'legacy_reference_check.json')
 (out/'ref').symlink_to((source/'ref').resolve(),target_is_directory=True)
 for p in source.glob('HmRos.*'):
  if p.is_file():shutil.copy2(p,out/p.name)
 for p in source.glob('Hmel2.5.*'):
  if p.is_file():shutil.copy2(p,out/p.name)
 for p in (source/'calls').iterdir():
  if p.is_file() and not '.partial' in p.name:(out/'calls'/p.name).symlink_to(p.resolve())
 phen={r['sample']:r for r in read(a.ledger/'resolved_phenotypes.tsv')}
 newtasks=[];deferred=[]
 for old in scope['new_tasks']:
  t=dict(old);sample=t['sample'];oldbam=source/'bam'/f'{sample}.chr15.dedup.bam'
  if oldbam.is_file() and (source/'bam'/f'{sample}.complete.json').is_file():t['reuse_bam']=str(oldbam)
  complete,why=ready_reads(source,t)
  if not t['reuse_bam'] and not complete:deferred.append(dict(sample=sample,reason=why));continue
  if (source/'raw'/sample).is_dir():(out/'raw'/sample).symlink_to((source/'raw'/sample).resolve(),target_is_directory=True)
  ph=phen.get(sample,{})
  for k in ('d','v'):
   if ph.get(k) in {'0','1'}:t[k]=int(ph[k])
  newtasks.append(t)
 cfg.update(output=str(out),recovery_source=str(source),new_tasks=newtasks,expansion_code_sha=digest(out/'code/expansion.py'),ledger=str(a.ledger),phenotype_ledger_sha=digest(a.ledger/'resolved_phenotypes.tsv'))
 # Preserve numeric call-task indices used by the frozen caller.
 taskmap={t['call_task']:t for t in newtasks}
 for t in cfg['tasks']:
  if t['task'] in taskmap:
   x=taskmap[t['task']];t['bam']=x['reuse_bam'] or str(out/'bam'/f'{x["sample"]}.chr15.dedup.bam')
 save(out/'run.json',cfg);tab(out/'deferred_without_completed_inputs.tsv',deferred)
 shutil.copy2(a.ledger/'resolved_phenotypes.tsv',out/'resolved_phenotypes.tsv')
 if a.submit:
  python=sys.executable;script=out/'code/recover_Ybv_from_completed_reads.py'
  def sb(stage,extra):
   wrap=shlex.join([python,str(script),'--stage',stage,'--output',str(out)])
   cmd=['sbatch','--parsable','--partition=cpu','--cpus-per-task='+('16' if stage=='worker' else '1'),'--mem='+('32G' if stage=='worker' else '8G'),'--time='+('1-00:00:00' if stage=='worker' else '01:00:00'),'--job-name=Ybv_recover_'+stage,'--output='+str(out/'logs'/f'{stage}.%A_%a.out'),'--error='+str(out/'logs'/f'{stage}.%A_%a.err'),*extra,'--wrap='+wrap]
   return subprocess.check_output(cmd,text=True).strip().split(';')[0]
  jobs=sb('worker',['--array=0-'+str(len(newtasks)-1)+'%48']);reportjob=sb('report',['--dependency=afterany:'+jobs]);save(out/'jobs.json',dict(align_and_call=jobs,report=reportjob))
  print('ALIGN_AND_CALL\t'+jobs+'\nREPORT\t'+reportjob)
 print('EXISTING_INPUT_TASKS\t'+str(len(newtasks))+'\nDEFERRED_INCOMPLETE_INPUTS\t'+str(len(deferred))+'\nOUTPUT\t'+str(out))
def worker(cfg):
 out=Path(cfg['output']);i=int(os.environ['SLURM_ARRAY_TASK_ID']);t=cfg['new_tasks'][i];status=dict(sample=t['sample'],index=i,original_index=t['index'],status='FAILED')
 try:
  m=module(cfg)
  if not t['reuse_bam']:
   okay,why=ready_reads(Path(cfg['recovery_source']),t)
   if not okay:raise RuntimeError(why)
  m.align(cfg,t);m.call(cfg,t)
  status['status']='COMPLETE'
 except Exception as e:
  status['error']=str(e);raise
 finally:save(out/'tasks'/f'{i:04d}.json',status)
def report(cfg):
 import pandas as pd,numpy as np
 out=Path(cfg['output']);statuses=[];ready=[]
 for i,t in enumerate(cfg['new_tasks']):
  p=out/'tasks'/f'{i:04d}.json';s=json.loads(p.read_text()) if p.is_file() else dict(sample=t['sample'],status='MISSING_TASK_OUTPUT')
  statuses.append(s)
  if s['status']=='COMPLETE':ready.append(t)
 tab(out/'stage_completion.tsv',statuses)
 if len(ready)<len(cfg['new_tasks'])*.95:
  save(out/'complete.json',dict(status='INCOMPLETE_NO_FINAL_FIGURE',planned=len(cfg['new_tasks']),completed=len(ready)))
  raise RuntimeError('Fewer than 95% of retained tasks completed; no final PDF generated; inspect stage_completion.tsv and worker logs')
 m=module(cfg);eng=m.engine(cfg);parent=Path(cfg['parent'])
 selected=pd.read_csv(parent/'selected_sources.tsv',sep='\t')
 for t in ready:
  selected=selected[~(selected['sample'].eq(t['sample'])&selected.reference.eq('HmRos')&selected.interval.eq('v'))]
  selected=pd.concat([selected,pd.DataFrame([dict(sample=t['sample'],taxon=t['taxon'],d=t['d'],v=t['v'],reference='HmRos',interval='v',source='forced_BAM_call',path=str(out/'calls'/f'{t["call_task"]}.vcf.gz'))])],ignore_index=True)
 ph={r['sample']:r for r in read(out/'resolved_phenotypes.tsv')}
 excluded=[];conflicts=[]
 for idx,r in selected.iterrows():
  p=ph.get(r['sample'],{})
  if p.get('exclusion'):excluded.append(idx);continue
  for k in ('d','v'):
   if p.get(k) in {'0','1'}:selected.at[idx,k]=int(p[k])
  if p.get('phenotype_conflict')=='True':conflicts.append(dict(sample=r['sample'],d=p.get('d'),v=p.get('v'),note='Documented literature override retained alongside sheet conflict'))
 selected=selected.drop(index=excluded)
 # One row per biological individual; combine same-BioSample calls only when
 # callable dosages agree. Conflicts become unavailable, never extra individuals.
 groups=defaultdict(list)
 for idx,r in selected.iterrows():groups[(ph.get(r['sample'],{}).get('individual_id',r['sample']),r.reference,r.interval)].append((idx,r))
 original_extract=eng.extract;aliases={};drop=[];merges=[]
 for key,group in groups.items():
  if len(group)<2:continue
  group=sorted(group,key=lambda x:x[1]['sample']);representative=group[0][1]['sample']
  if len({tuple(str(r[k]) for k in ('d','v')) for _,r in group})>1:raise RuntimeError('Conflicting phenotypes within biological individual '+str(key))
  aliases[(representative,key[1],key[2])]=[dict(sample=r['sample'],path=r.path,source=r.source) for _,r in group]
  drop.extend(idx for idx,_ in group[1:]);merges.append(dict(individual=key[0],reference=key[1],interval=key[2],representative=representative,runs=';'.join(r['sample'] for _,r in group),rule='union concordant callable genotypes; discordance unavailable'))
 selected=selected.drop(index=drop);selected.to_csv(out/'selected_sources.tsv',sep='\t',index=False)
 def extract_with_individuals(path,samples,table,legacy=False,allow_pl=False):
  result,checks=original_extract(path,samples,table,legacy,allow_pl)
  ref='HmRos' if str(table.contig.iloc[0])=='NC_140403.1' else 'Hmel2.5'
  interval='v' if int(table.pos1.min())>(1650000 if ref=='HmRos' else 1420000) else 'd'
  for s in samples:
   aa=aliases.get((s,ref,interval),[])
   if not aa:continue
   vals=[]
   for a in aa:
    x,_=original_extract(a['path'],[a['sample']],table,legacy=a['source']=='legacy_VCF',allow_pl=True);vals.append(x[a['sample']])
   mat=np.asarray(vals);n=np.isfinite(mat).sum(axis=0)
   minimum=np.min(np.where(np.isfinite(mat),mat,np.inf),axis=0);maximum=np.max(np.where(np.isfinite(mat),mat,-np.inf),axis=0)
   result[s]=np.where((n>0)&(minimum==maximum),minimum,np.nan)
  return result,checks
 eng.extract=extract_with_individuals
 wanted={t['call_task'] for t in ready}|{t['task'] for t in cfg['original_tasks']};local=dict(cfg);local['tasks']=[t for t in cfg['tasks'] if t['task'] in wanted]
 eng.report(local)
 result=json.loads((out/'complete.json').read_text())
 tab(out/'individual_merges.tsv',merges);tab(out/'phenotype_conflicts.tsv',conflicts)
 result.update(status='COMPLETE' if len(ready)==len(cfg['new_tasks']) else 'COMPLETE_WITH_TASK_FAILURES_REQUIRES_REVIEW',retained_tasks=len(cfg['new_tasks']),completed_tasks=len(ready),recovery_source=cfg['recovery_source'],phenotype_ledger=str(out/'resolved_phenotypes.tsv'),phenotype_ledger_sha256=cfg['phenotype_ledger_sha'],biological_individual_merges=merges,no_downloads=True,marker_selection_unchanged=True,focal_calls_unchanged=True)
 save(out/'complete.json',result)
 print('STATUS\t'+result['status'])
 for page in result['pages']:print('PANEL\t'+json.dumps(page))
 print('PDF\t'+str(out/'Yb_expanded_existing_controls_Rosina.pdf'))
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,default=ROOT/'analysis/yb_ventral_population_expansion_v1');p.add_argument('--ledger',type=Path,default=ROOT/'analysis/yb_individual_ledger_v1');p.add_argument('--output',type=Path,required=True);p.add_argument('--submit',action='store_true');p.add_argument('--stage',choices=['prepare','worker','report'],default='prepare');a=p.parse_args()
 if a.stage=='prepare':prepare(a)
 else:globals()[a.stage](json.loads((a.output/'run.json').read_text()))
if __name__=='__main__':
 from collections import defaultdict
 main()
