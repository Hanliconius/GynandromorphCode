#!/usr/bin/env python3
"""Read-only Pegasus storage inventory and conservative cleanup review.
Only creates its new report directory. Never deletes, moves, or edits project data.
"""
import argparse,csv,gzip,hashlib,json,os,re,shlex,sqlite3,stat,subprocess,tempfile,time
from collections import defaultdict,deque
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
REPORTS=['complete.json','directory_sizes.tsv','largest_files.tsv','cleanup_candidates.tsv','duplicate_candidates.tsv','protected_paths.tsv','active_jobs.tsv','dependency_summary.tsv','read_errors.tsv','symlink_summary.json','REPORT.md']
def write_json(p,obj):p.write_text(json.dumps(obj,indent=2)+'\n')
def table(path,rows,columns):
 with path.open('w',newline='') as h:
  w=csv.DictWriter(h,fieldnames=columns,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
def inside(p,root):
 try:Path(p).relative_to(root);return True
 except ValueError:return False

def absolute_refs(value):
 if isinstance(value,dict):
  for k,v in value.items():
   if k.lower() in {'root','project_root','output_root'}:continue
   yield from absolute_refs(v)
 elif isinstance(value,list):
  for v in value:yield from absolute_refs(v)
 elif isinstance(value,str) and value.startswith('/') and '\n' not in value and not value.startswith('/full/path/to/'):
  # Existing configurations sometimes store comma-separated mate paths.
  for v in value.split(','):
   if v.startswith('/'):yield Path(v)

def infer_run(path,root):
 p=Path(path)
 if not inside(p,root):return None
 p=p if p.is_dir() else p.parent
 for parent in [p,*p.parents]:
  if parent==root:break
  if any((parent/name).is_file() for name in ('run.json','config.json')):return parent
 return None

def live_jobs(root,exclude):
 jobs=[];runs=set();unresolved=[]
 try:
  proc=subprocess.run(['squeue','--me','--noheader','--format=%i|%T|%j'],capture_output=True,text=True,check=True)
  ids={line.split('|')[0].split('_')[0].strip() for line in proc.stdout.splitlines() if line.strip()}
  ids.discard(exclude)
  for job in sorted(ids):
   result=subprocess.run(['scontrol','show','job',job,'-o'],capture_output=True,text=True)
   if result.returncode:unresolved.append(job);continue
   info=result.stdout
   matches=re.findall(re.escape(str(root))+r'/[^\s\"\'<>;]+',info)
   found=set()
   for m in matches:
    run=infer_run(m,root)
    if run:found.add(str(run));runs.add(run)
   name=re.search(r'\bJobName=(\S+)',info);state=re.search(r'\bJobState=(\S+)',info)
   relevant=bool(matches)
   if relevant and not found:unresolved.append(job)
   jobs.append(dict(job_id=job,name=name[1] if name else '',state=state[1] if state else '',project_related=relevant,protected_runs=';'.join(sorted(found))))
  return jobs,runs,unresolved,''
 except Exception as e:return jobs,runs,unresolved,str(e)

def protected_inputs(root,initial,errors):
 protected={};queue=deque();seen=set()
 def add(path,reason):
  p=Path(path)
  if p==root or not inside(p,root):return
  protected.setdefault(str(p),set()).add(reason)
  try:
   target=p.resolve()
   if target!=p and inside(target,root):protected.setdefault(str(target),set()).add('target_of_protected_symlink')
  except OSError:pass
 for run in initial:
  add(run,'active_or_explicitly_protected_run')
  for n in ('run.json','config.json','audit.json'):
   if (run/n).is_file():queue.append(run/n)
 while queue:
  p=queue.popleft()
  try:
   real=p.resolve()
   if real in seen or p.stat().st_size>10_000_000:continue
   seen.add(real)
   if p.suffix=='.json':obj=json.loads(p.read_text())
   elif p.suffix in {'.tsv','.csv'}:
    with p.open() as h:obj=list(csv.DictReader(h,delimiter=',' if p.suffix=='.csv' else '\t'))
   else:continue
   for ref in absolute_refs(obj):
    if not inside(ref,root) or ref==root:continue
    add(ref,'input_referenced_by_protected_configuration')
    # Follow source run configs and explicit nested config files recursively.
    if ref.is_dir():
     for n in ('run.json','config.json','audit.json'):
      if (ref/n).is_file():queue.append(ref/n)
    elif ref.suffix in {'.json','.tsv','.csv'} and ref.is_file():queue.append(ref)
  except Exception as e:errors.append(dict(path=str(p),error='protected provenance: '+str(e)))
 return protected

def match_protected(path,protected):
 p=Path(path)
 for parent in [p,*p.parents]:
  if str(parent) in protected:return True
 return False

def file_category(name):
 n=name.lower()
 if n.endswith(('.fastq','.fq','.fastq.gz','.fq.gz','.sra')):return 'raw_reads'
 if n.endswith(('.bam','.cram','.sam')):return 'alignment'
 if n.endswith(('.vcf','.vcf.gz','.bcf','.g.vcf.gz')):return 'variant_calls'
 if n.endswith(('.fa','.fasta','.fa.gz','.fasta.gz','.fna','.fna.gz','.bwt','.pac','.sa')):return 'reference_or_index'
 if n.endswith(('.pdf','.png','.svg')):return 'figure'
 return 'other'

def edge_signature(path,size):
 h=hashlib.sha256()
 with path.open('rb') as f:
  h.update(f.read(65536));f.seek(max(0,size-65536));h.update(f.read(65536))
 return h.hexdigest()

def audit(a):
 root=a.root.resolve();out=a.output.resolve()
 if not root.is_dir():raise RuntimeError('Project root not found')
 if out.exists():raise RuntimeError('Use a new output directory')
 out.mkdir(parents=True);started=time.time();errors=[]
 jobs,live,unresolved,queue_error=live_jobs(root,os.environ.get('SLURM_ARRAY_JOB_ID',os.environ.get('SLURM_JOB_ID','')))
 if queue_error:errors.append(dict(path='Slurm',error=queue_error))
 explicit={Path(p).resolve() for p in a.protect_run}
 protected=protected_inputs(root,live|explicit,errors)
 table(out/'active_jobs.tsv',jobs,['job_id','name','state','project_related','protected_runs'])
 table(out/'protected_paths.tsv',[dict(path=p,reasons=';'.join(sorted(v))) for p,v in sorted(protected.items())],['path','reasons'])
 db=sqlite3.connect(out/'inventory.sqlite');db.execute('CREATE TABLE files(path TEXT PRIMARY KEY, bytes INTEGER, allocated INTEGER, mtime REAL, dev INTEGER, inode INTEGER, nlink INTEGER, category TEXT, protected INTEGER, kind TEXT, target TEXT)')
 dirs=defaultdict(lambda:dict(apparent_bytes=0,allocated_unique_bytes=0,regular_files=0,symlinks=0));inode_seen=set();configs=[];symlink_count=broken=outside_targets=0;file_count=0
 for parent,subdirs,names in os.walk(root,followlinks=False,onerror=lambda e:errors.append(dict(path=str(e.filename),error=str(e)))):
  parent=Path(parent);subdirs.sort();names.sort()
  # Never scan report outputs while they are being written.
  subdirs[:]=[d for d in subdirs if parent/d!=out and not inside(parent/d,out)]
  entries=[parent/n for n in names]+[parent/d for d in subdirs if (parent/d).is_symlink()]
  subdirs[:]=[d for d in subdirs if not (parent/d).is_symlink()]
  for p in entries:
   try:
    s=p.lstat();islink=stat.S_ISLNK(s.st_mode);regular=stat.S_ISREG(s.st_mode)
    if not islink and not regular:continue
    allocated=s.st_blocks*512;identity=(s.st_dev,s.st_ino);unique=identity not in inode_seen;inode_seen.add(identity)
    target=''
    if islink:
     symlink_count+=1;target=str(p.resolve(strict=False));broken+=not p.exists();outside_targets+=not inside(target,root)
    prot=match_protected(p,protected) or (islink and match_protected(target,protected))
    db.execute('INSERT INTO files VALUES(?,?,?,?,?,?,?,?,?,?,?)',(str(p),s.st_size,allocated if unique else 0,s.st_mtime,s.st_dev,s.st_ino,s.st_nlink,file_category(p.name),int(prot),'symlink' if islink else 'regular',target))
    file_count+=1
    for ancestor in [parent,*parent.parents]:
     if not inside(ancestor,root):break
     r=dirs[str(ancestor)];r['apparent_bytes']+=s.st_size;r['allocated_unique_bytes']+=allocated if unique else 0;r['regular_files']+=regular;r['symlinks']+=islink
    if regular and s.st_size<10_000_000 and (p.name in {'run.json','config.json','audit.json'} or ('manifest' in p.name.lower() and p.suffix in {'.json','.tsv','.csv'})):configs.append(p)
   except (OSError,ValueError) as e:errors.append(dict(path=str(p),error=str(e)))
  db.commit()
 db.execute('CREATE INDEX file_sizes ON files(bytes)');db.commit()
 # Dependencies include all historical configurations, not only live ones.
 dependencies=defaultdict(set);seen=set();queue=deque(configs)
 while queue:
  p=queue.popleft()
  try:
   real=p.resolve()
   if real in seen or not p.is_file() or p.stat().st_size>10_000_000:continue
   seen.add(real)
   if p.suffix=='.json':refs=absolute_refs(json.loads(p.read_text()))
   else:
    with p.open() as h:rows=list(csv.DictReader(h,delimiter=',' if p.suffix=='.csv' else '\t'))
    refs=(ref for row in rows for ref in absolute_refs(row))
   for ref in refs:
    if ref==root:continue
    dependencies[str(ref)].add(str(p))
    try:dependencies[str(ref.resolve())].add(str(p))
    except OSError:pass
    if ref.is_dir():
     for name in ('run.json','config.json','audit.json'):
      if (ref/name).is_file():queue.append(ref/name)
    elif ref.suffix in {'.json','.tsv','.csv'} and ref.is_file():queue.append(ref)
  except Exception as e:errors.append(dict(path=str(p),error='dependency scan: '+str(e)))
 def references(path):
  return sorted(set().union(*(dependencies.get(str(p),set()) for p in [Path(path),*Path(path).parents] if p!=root)))
 db.row_factory=sqlite3.Row
 with gzip.open(out/'file_inventory.tsv.gz','wt',newline='') as h:
  cols=['path','bytes','allocated','mtime','dev','inode','nlink','category','protected','kind','target'];w=csv.DictWriter(h,delimiter='\t',fieldnames=cols);w.writeheader()
  for r in db.execute('SELECT * FROM files ORDER BY path'):w.writerow(dict(r))
 directory_rows=[]
 for path,r in dirs.items():
  depth=len(Path(path).relative_to(root).parts)
  directory_rows.append(dict(path=path,depth=depth,protected=match_protected(path,protected),**r))
 directory_rows.sort(key=lambda r:r['allocated_unique_bytes'],reverse=True)
 table(out/'directory_sizes.tsv',[r for r in directory_rows if r['depth']<=4],['path','depth','allocated_unique_bytes','apparent_bytes','regular_files','symlinks','protected'])
 largest=[dict(r) for r in db.execute("SELECT * FROM files WHERE kind='regular' ORDER BY bytes DESC LIMIT 200")]
 for r in largest:r['referenced_by_config_count']=len(references(r['path']))
 table(out/'largest_files.tsv',largest,['path','bytes','allocated','category','protected','nlink','referenced_by_config_count'])
 candidates=[]
 for rr in db.execute("SELECT * FROM files WHERE kind='regular' ORDER BY bytes DESC"):
  r=dict(rr);p=Path(r['path']);reason='';age=(started-r['mtime'])/86400;parts=set(p.relative_to(root).parts)
  if re.search(r'\.(name|fixmate|position|sorted|raw)\.bam$',p.name) or parts & {'tmp','temp'}:reason='temporary_or_processing_intermediate_requires_derivative_check'
  elif p.name.endswith(('.partial','.part','.tmp')):reason='partial_file_requires_failed_run_and_integrity_check'
  elif 'toolenv' in parts:reason='analysis_specific_environment_requires_shared_replacement_review'
  if not reason:continue
  refs=references(r['path']);status='protected_live_input' if r['protected'] else 'hardlink_not_independently_reclaimable' if r['nlink']>1 else 'referenced_by_configuration' if refs else 'recent_or_changing_review' if age<1 else 'candidate_for_review'
  if unresolved or queue_error:
   if status=='candidate_for_review':status='review_blocked_by_unresolved_live_job_provenance'
  candidates.append(dict(path=r['path'],bytes=r['bytes'],allocated_bytes=r['allocated'],age_days=round(age,2),reason=reason,status=status,config_references=len(refs),sample_source_configs=';'.join(refs[:3])))
 table(out/'cleanup_candidates.tsv',candidates[:2000],['path','bytes','allocated_bytes','age_days','reason','status','config_references','sample_source_configs'])
 # Size plus first/last 64 KiB identifies candidates, never proven duplicates.
 duplicates=[];budget=64*1024*1024;hashed=0;by_size=defaultdict(list)
 for rr in db.execute("SELECT * FROM files WHERE kind='regular' AND bytes>=33554432 AND nlink=1 AND protected=0 ORDER BY bytes DESC"):
  r=dict(rr)
  if started-r['mtime']>=86400:by_size[r['bytes']].append(r)
 for size,group in by_size.items():
  if len(group)<2:continue
  buckets=defaultdict(list)
  for r in group:
   if hashed+131072>budget:break
   try:
    p=Path(r['path']);before=p.stat();sig=edge_signature(p,size);after=p.stat();hashed+=131072
    if before.st_size!=after.st_size or before.st_mtime_ns!=after.st_mtime_ns:continue
    buckets[sig].append(r)
   except OSError as e:errors.append(dict(path=r['path'],error=str(e)))
  for sig,matches in buckets.items():
   if len(matches)<2:continue
   gid=f'edge_{len(duplicates):05d}'
   for r in matches:duplicates.append(dict(group=gid,path=r['path'],bytes=size,signature=sig,evidence='same_size_and_first_last_64KiB_only_NOT_proven_duplicate',config_references=len(references(r['path']))))
  if hashed>=budget:break
 table(out/'duplicate_candidates.tsv',duplicates[:2000],['group','path','bytes','signature','evidence','config_references'])
 dep_rows=[dict(path=p,config_references=len(src),example_configs=';'.join(sorted(src)[:3])) for p,src in dependencies.items()];dep_rows.sort(key=lambda r:r['config_references'],reverse=True)
 table(out/'dependency_summary.tsv',dep_rows[:2000],['path','config_references','example_configs'])
 table(out/'read_errors.tsv',errors,['path','error']);write_json(out/'symlink_summary.json',dict(symlinks=symlink_count,broken_links=broken,targets_outside_project=outside_targets,symlink_targets_not_traversed=True))
 total=dirs.get(str(root),{});reviewable=[r for r in candidates if r['status']=='candidate_for_review']
 result=dict(status='READ_ONLY_AUDIT_COMPLETE_REQUIRES_REVIEW',root=str(root),output=str(out),started_utc=datetime.fromtimestamp(started,timezone.utc).isoformat(),finished_utc=datetime.now(timezone.utc).isoformat(),files_and_symlinks=file_count,directories=len(dirs),allocated_unique_bytes=total.get('allocated_unique_bytes',0),apparent_bytes=total.get('apparent_bytes',0),active_project_runs=[str(p) for p in sorted(live)],explicit_protected_runs=[str(p) for p in sorted(explicit)],protected_paths=len(protected),unresolved_active_job_ids=unresolved,slurm_snapshot_error=queue_error,read_errors=len(errors),review_candidate_files=len(reviewable),review_candidate_allocated_bytes_upper_bound=sum(r['allocated_bytes'] for r in reviewable),candidate_bytes_are_not_approved_deletion_savings=True,duplicate_evidence='edge signatures only; full-file confirmation required',duplicate_signature_bytes_read_budget=budget,duplicate_signature_bytes_read=hashed,mutations='New report files only; no project files deleted, moved, or edited',live_file_sizes='Point-in-time and provisional; no content hashing for protected or recently modified files',script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 write_json(out/'complete.json',result)
 report=f'''# Gynandromorph storage audit\n\nAllocated regular-file and symlink blocks counted once per inode: {result['allocated_unique_bytes']/2**40:.3f} TiB.\nApparent file sizes: {result['apparent_bytes']/2**40:.3f} TiB. Directory metadata and filesystem overhead are not included.\n\nThis is a read-only snapshot. Active run trees and their explicitly traced input dependencies are protected.\nLive file sizes are provisional. Historical configuration references are retained as dependency evidence.\nAn unreferenced intermediate is a review candidate, not permission to remove it.\n\nUse directory_sizes.tsv and largest_files.tsv to identify storage concentration.\ncleanup_candidates.tsv contains intermediates, partials and per-run environments, with protection/reference status.\nduplicate_candidates.tsv uses size and edge signatures only; it does not prove whole-file identity or reclaimable bytes.\nSymlinks and hardlinks are not counted as additional copies of their targets. Removing one hardlink may reclaim no space.\n\nThe full compressed inventory and SQLite index remain on Pegasus. Only compact reports are pushed to Git.\nNo deletion commands or cleanup operations are generated.\n'''
 (out/'REPORT.md').write_text(report);db.close()
 print('AUDIT\t'+str(out));print('ALLOCATED_TIB\t'+f"{result['allocated_unique_bytes']/2**40:.3f}");print('READ_ERRORS\t'+str(len(errors)))

def publish(a):
 out=a.output.resolve();summary=json.loads((out/'complete.json').read_text())
 if not summary['status'].startswith('READ_ONLY_AUDIT_COMPLETE'):raise RuntimeError('Audit not complete')
 repo=a.push_repo.resolve();branch='storage-audit-'+out.name
 env=dict(os.environ,GIT_SSH_COMMAND='ssh -o BatchMode=yes')
 def git(*args,cwd=repo):return subprocess.run(['git','-C',str(cwd),*map(str,args)],check=True,env=env,stdout=subprocess.DEVNULL)
 # Only this helper's disposable checkout is removed; project data are untouched.
 with tempfile.TemporaryDirectory(prefix='gynandro-storage-publish.') as tmp:
  work=Path(tmp)/'checkout';git('fetch','--quiet','origin','main');git('worktree','add','--quiet','--detach',work,'origin/main')
  try:
   dest=work/'results/storage_audits'/out.name;dest.mkdir(parents=True)
   import shutil
   for f in REPORTS:shutil.copy2(out/f,dest/f)
   if a.sequencing_report:
    seqdest=dest/'sequencing';seqdest.mkdir()
    for f in ('complete.json','sequencing_and_reuse.tsv','lookup_errors.tsv'):
     shutil.copy2(a.sequencing_report/f,seqdest/f)
   git('add',str(dest.relative_to(work)),cwd=work);git('commit','--quiet','-m','Record read-only Pegasus storage audit '+out.name,cwd=work)
   git('push','--quiet','git@github.com:Hanliconius/GynandromorphCode.git','HEAD:refs/heads/'+branch,cwd=work)
  finally:git('worktree','remove','--force',work)
 write_json(out/'git_publish.json',dict(repository='Hanliconius/GynandromorphCode',branch=branch,path='results/storage_audits/'+out.name,files=REPORTS))
 print('GIT_BRANCH\t'+branch)

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=ROOT);p.add_argument('--output',type=Path,required=True);p.add_argument('--protect-run',type=Path,action='append',default=[]);p.add_argument('--push-repo',type=Path);p.add_argument('--sequencing-report',type=Path);p.add_argument('--publish-only',action='store_true');a=p.parse_args()
 if not a.publish_only:audit(a)
 if a.push_repo:publish(a)
if __name__=='__main__':main()
