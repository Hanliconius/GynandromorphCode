#!/usr/bin/env python3
"""Arrayed, read-only validation of storage candidates. Never deletes project data."""
import argparse,csv,gzip,hashlib,json,os,re,shlex,shutil,sqlite3,subprocess,sys,tempfile,time
from collections import defaultdict
from pathlib import Path
import audit_project_storage as audit
def run(cmd):
 p=subprocess.run(list(map(str,cmd)),capture_output=True,text=True)
 if p.returncode:raise RuntimeError('Command failed: '+shlex.join(list(map(str,cmd)))+'; '+p.stderr[-1500:])
 return p.stdout
def rows(p):
 with p.open() as h:return list(csv.DictReader(h,delimiter='\t'))
def dump(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def signature(p):
 s=p.stat();return (s.st_dev,s.st_ino,s.st_size,s.st_mtime_ns)
def digest(p):
 before=signature(p);h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 if signature(p)!=before:raise RuntimeError('File changed during full checksum: '+str(p))
 return h.hexdigest()
def fastq(p,sample):
 before=signature(p);n=0;matching=0
 with gzip.open(p,'rb') as f:
  while True:
   head=f.readline()
   if not head:break
   seq=f.readline().rstrip(b'\r\n');plus=f.readline();quality=f.readline().rstrip(b'\r\n')
   if not head.startswith(b'@') or not plus.startswith(b'+') or len(seq)!=len(quality) or not seq:raise RuntimeError('Invalid FASTQ record: '+str(p))
   n+=1;matching+=sample.encode() in head.split()[0]
 if signature(p)!=before:raise RuntimeError('FASTQ changed during validation')
 return dict(path=str(p),records=n,headers_containing_accession=matching,allocated=p.stat().st_blocks*512)
def prepare(a):
 source=a.source_audit.resolve();info=json.loads((source/'complete.json').read_text());root=Path(info['root'])
 if a.output.exists():raise RuntimeError('Choose a new output')
 out=a.output.resolve();out.mkdir(parents=True);(out/'tasks').mkdir();(out/'logs').mkdir();(out/'code').mkdir()
 for name in ('verify_storage_candidates.py','audit_project_storage.py'):shutil.copy2(Path(__file__).with_name(name),out/'code'/name)
 db=sqlite3.connect(source/'inventory.sqlite');db.row_factory=sqlite3.Row
 files=[dict(r) for r in db.execute("SELECT * FROM files WHERE kind='regular'")];db.close()
 by_path={r['path']:r for r in files};protected={r['path'] for r in files if r['protected']}
 tasks=[];scope=root/'analysis/genomewide_mixing_corrected_v1/full_sib_validation'
 groups=defaultdict(list)
 for r in files:
  p=Path(r['path'])
  if audit.inside(p,scope/'tmp/alignment') and re.search(r'\.(name|fixmate|position)\.bam$',p.name):
   sample=re.sub(r'\.(name|fixmate|position)\.bam$','',p.name);groups[sample].append(r['path'])
 for sample,paths in groups.items():
  final=scope/'alignment'/f'{sample}.dedup.bam'
  tasks.append(dict(kind='bam_intermediates',sample=sample,paths=paths,retained=str(final)))
 duplicate_groups=defaultdict(list)
 for r in rows(source/'duplicate_candidates.tsv'):duplicate_groups[r['group']].append(r['path'])
 for group,paths in duplicate_groups.items():
  if len(paths)>1:tasks.append(dict(kind='duplicate_files',sample=group,paths=paths))
 for r in files:
  p=Path(r['path'])
  if p.suffix=='.sra' and audit.inside(p,root/'analysis/genomewide_mixing_input/sra'):
   sample=p.stem;fq=root/'analysis/genomewide_mixing_input/fastq'/sample
   tasks.append(dict(kind='sra_with_retained_fastq',sample=sample,paths=[str(p)],retained=[str(fq/f'{sample}_{m}.fastq.gz') for m in (1,2)]))
 # Store exact size/inode/mtime expectations; workers refuse changed candidate files.
 for task in tasks:
  task['snapshot']={p:{k:by_path[p][k] for k in ('bytes','mtime','dev','inode','nlink','allocated','protected')} for p in task['paths']}
 cfg=dict(root=str(root),source_audit=str(source),output=str(out),tasks=tasks,protected_runs=info['explicit_protected_runs'],repo=str(a.push_repo) if a.push_repo else '',no_deletions=True)
 dump(out/'run.json',cfg)
 # Preserve run->individual relationships without silently changing running analyses.
 seq=rows(source.parent.parent/'sequencing_audits'/source.name/'sequencing_and_reuse.tsv') if (source.parent.parent/'sequencing_audits'/source.name/'sequencing_and_reuse.tsv').exists() else []
 identities=defaultdict(list)
 for r in seq:identities[r['sample_accession']].append(r['sample'])
 audit.table(out/'shared_BioSamples.tsv',[dict(BioSample=k,runs=';'.join(v),action='one_individual; combine compatible libraries before calling; never count runs as independent individuals') for k,v in identities.items() if k and len(v)>1],['BioSample','runs','action'])
 if not tasks:raise RuntimeError('No candidate tasks')
 if a.submit:
  script=out/'code/verify_storage_candidates.py'
  def sb(stage,extra):
   cmd=['sbatch','--parsable','--job-name=GynVerify_'+stage,'--partition=cpu','--cpus-per-task=2' if stage=='worker' else '--cpus-per-task=1','--mem=8G','--time=06:00:00' if stage=='worker' else '--time=00:30:00','--output='+str(out/'logs'/f'{stage}.%A_%a.out'),'--error='+str(out/'logs'/f'{stage}.%A_%a.err'),*extra,'--wrap=micromamba run -n gynandrocalling python '+shlex.quote(str(script))+' --stage '+stage+' --output '+shlex.quote(str(out))]
   return run(cmd).strip().split(';')[0]
  job=sb('worker',['--array=0-'+str(len(tasks)-1)+'%8'])
  report=sb('report',['--dependency=afterany:'+job])
  dump(out/'jobs.json',dict(validation=job,report=report));print('VALIDATION\t'+job+'\nREPORT\t'+report)
 print('OUTPUT\t'+str(out));print('TASKS\t'+str(len(tasks)))
def worker(cfg):
 out=Path(cfg['output']);idx=int(os.environ['SLURM_ARRAY_TASK_ID']);task=cfg['tasks'][idx];root=Path(cfg['root']);errors=[]
 result=dict(task=idx,kind=task['kind'],sample=task['sample'],status='REQUIRES_REVIEW',candidate_bytes=0,checks={},paths=task['paths'])
 try:
  _,live,unresolved,error=audit.live_jobs(root,os.environ.get('SLURM_ARRAY_JOB_ID',''))
  live.update(Path(p) for p in cfg['protected_runs'])
  protected=audit.protected_inputs(root,live,errors)
  if unresolved or error or errors:raise RuntimeError('Live input protection unresolved: '+str((unresolved,error,errors)))
  for name in task['paths']:
   p=Path(name);s=p.stat();old=task['snapshot'][name]
   if p.is_symlink() or s.st_nlink>1:raise RuntimeError('Symlink or hardlink requires separate review')
   if old['protected'] or audit.match_protected(p,protected):raise RuntimeError('Protected input: '+name)
   if (s.st_dev,s.st_ino,s.st_size)!=(old['dev'],old['inode'],old['bytes']) or abs(s.st_mtime-old['mtime'])>0.001:raise RuntimeError('Candidate changed since inventory')
   if time.time()-s.st_mtime<86400:raise RuntimeError('Candidate modified within one day')
   result['candidate_bytes']+=s.st_blocks*512
  if task['kind']=='duplicate_files':
   sums={p:digest(Path(p)) for p in task['paths']};result['checks']['full_SHA256']=sums
   if len(set(sums.values()))!=1:raise RuntimeError('Edge-matched files differ internally')
   # Retain one canonical file; this is an upper bound before reference repair.
   result['candidate_bytes']-=max(Path(p).stat().st_blocks*512 for p in task['paths'])
   result['status']='EXACT_DUPLICATES_REFERENCE_REPAIR_REQUIRED'
  elif task['kind']=='bam_intermediates':
   p=Path(task['retained']);before=signature(p);run(['samtools','quickcheck','-v',p])
   header=run(['samtools','view','-H',p]);sms=sorted(set(re.findall(r'(?:^|\t)SM:([^\t\n]+)',header)))
   if sms!=[task['sample']]:raise RuntimeError('Unexpected retained BAM sample identity: '+str(sms))
   total=int(run(['samtools','view','-@','1','-c',p]).strip())
   if total<=0:raise RuntimeError('Empty retained BAM')
   if signature(p)!=before:raise RuntimeError('Retained BAM changed during validation')
   index=next((x for x in [Path(str(p)+'.bai'),p.with_suffix('.bai'),Path(str(p)+'.csi')] if x.is_file()),None)
   if index is None:raise RuntimeError('Retained BAM index absent')
   result['checks']=dict(retained=str(p),SM=sms,decoded_records=total,index=str(index),index_stats=run(['samtools','idxstats',p]))
   result['status']='DERIVATIVE_VALIDATED_PROVENANCE_REVIEW_REQUIRED'
  else:
   fqs=[fastq(Path(p),task['sample']) for p in task['retained']]
   if fqs[0]['records']==0 or fqs[0]['records']!=fqs[1]['records']:raise RuntimeError('Empty or unbalanced retained FASTQ pair')
   if any(f['headers_containing_accession']!=f['records'] for f in fqs):raise RuntimeError('FASTQ header accession not established')
   result['checks']=dict(retained_fastqs=fqs,note='gzip CRC and FASTQ structure validated; record counts agree; SRA-to-FASTQ completeness still requires archive spot-count comparison')
   result['status']='RETAINED_READS_VALIDATED_SRA_COMPLETENESS_REVIEW_REQUIRED'
 except Exception as e:result['status']='BLOCKED';result['error']=str(e);result['candidate_bytes']=0
 dump(out/'tasks'/f'{idx:04d}.json',result)
 print(str(idx)+'\t'+result['status'])
def report(cfg):
 out=Path(cfg['output']);results=[]
 for i,task in enumerate(cfg['tasks']):
  p=out/'tasks'/f'{i:04d}.json'
  results.append(json.loads(p.read_text()) if p.is_file() else dict(task=i,kind=task['kind'],sample=task['sample'],status='MISSING_TASK_OUTPUT',candidate_bytes=0))
 audit.table(out/'verification_summary.tsv',results,['task','kind','sample','status','candidate_bytes','error'])
 counts=defaultdict(int);sizes=defaultdict(int)
 for r in results:counts[r['status']]+=1;sizes[r['status']]+=r['candidate_bytes']
 dump(out/'complete.json',dict(status='VALIDATION_COMPLETE_REQUIRES_REVIEW' if all(r['status']!='MISSING_TASK_OUTPUT' for r in results) else 'INCOMPLETE',tasks=len(results),statuses=dict(counts),candidate_bytes_by_status=dict(sizes),candidate_bytes_are_upper_bounds_not_approved_savings=True,no_deletions=True))
 if cfg['repo']:
  repo=Path(cfg['repo']);branch='storage-validation-'+out.name;env=dict(os.environ,GIT_SSH_COMMAND='ssh -o BatchMode=yes')
  def git(*args,cwd=repo):subprocess.run(['git','-C',str(cwd),*map(str,args)],check=True,env=env,stdout=subprocess.DEVNULL)
  with tempfile.TemporaryDirectory(prefix='storage-validation-publish.') as tmp:
   work=Path(tmp)/'checkout';git('fetch','--quiet','origin','main');git('worktree','add','--quiet','--detach',work,'origin/main')
   try:
    dest=work/'results/storage_validation'/out.name;dest.mkdir(parents=True)
    for f in ('complete.json','verification_summary.tsv','shared_BioSamples.tsv'):shutil.copy2(out/f,dest/f)
    shutil.copytree(out/'tasks',dest/'tasks')
    git('add',str(dest.relative_to(work)),cwd=work);git('commit','--quiet','-m','Record read-only storage candidate validation',cwd=work);git('push','--quiet','git@github.com:Hanliconius/GynandromorphCode.git','HEAD:refs/heads/'+branch,cwd=work)
   finally:git('worktree','remove','--force',work)
  print('GIT_BRANCH\t'+branch)
 print(json.dumps(dict(counts)));print('OUTPUT\t'+str(out))
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stage',choices=['prepare','worker','report'],default='prepare');p.add_argument('--source-audit',type=Path);p.add_argument('--output',type=Path,required=True);p.add_argument('--push-repo',type=Path);p.add_argument('--submit',action='store_true');a=p.parse_args()
 if a.stage=='prepare':prepare(a)
 else:
  cfg=json.loads((a.output/'run.json').read_text());globals()[a.stage](cfg)
if __name__=='__main__':main()
