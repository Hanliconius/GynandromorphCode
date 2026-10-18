#!/usr/bin/env python3
"""Narrow a submitted expansion to capture plus completed/in-flight WGS; preserve data."""
import argparse,csv,hashlib,importlib.util,json,shlex,subprocess,sys,time,os
from pathlib import Path
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def command(args,check=True):return subprocess.run(list(map(str,args)),check=check,capture_output=True,text=True)
def expression(indices):
 vals=sorted(set(indices));groups=[]
 for n in vals:
  if groups and n==groups[-1][1]+1:groups[-1][1]=n
  else:groups.append([n,n])
 return ','.join(str(a) if a==b else f'{a}-{b}' for a,b in groups)
def decisions(tasks,seq,states,root,jobs):
 result=[]
 for t in tasks:
  sample=t['sample'];r=seq[sample];i=t['index'];state=states.get(f"{jobs['download']}_{i}",'')
  capture=r['library_strategy'].upper()=='TARGETED-CAPTURE' and r['library_selection'].lower()=='hybrid selection'
  complete=(root/'raw'/sample/'complete.json').is_file()
  inflight=state in {'RUNNING','COMPLETING','R','CG'}
  fq=t.get('fastq_candidates',{})
  local=bool(t.get('reuse_bam')) and Path(t['reuse_bam']).is_file()
  for mate in (1,2):
   candidates=fq.get(str(mate),fq.get(mate,[]))
   if not any(Path(x).is_file() and Path(x).stat().st_size>0 for x in candidates):break
  else:local=True
  keep=capture or complete or local
  result.append(dict(index=i,sample=sample,keep=keep,reason='confirmed_targeted_capture' if capture else 'WGS_download_completed' if complete else 'WGS_complete_local_inputs' if local else 'excluded_WGS_download_incomplete',download_state=state))
 return result
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,required=True);p.add_argument('--sequencing-report',type=Path);p.add_argument('--apply',action='store_true');p.add_argument('--report-config',type=Path);a=p.parse_args();root=a.run.resolve()
 if a.report_config:
  cfg=json.loads(a.report_config.read_text());spec=importlib.util.spec_from_file_location('frozen_expansion',root/'expand_Yb_ventral_population.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
  m.tool_environment(cfg);m.report(cfg)
  complete=json.loads((root/'complete.json').read_text());scope=json.loads((a.report_config.parent/'scope.json').read_text());complete['capture_focus']=scope;save(root/'complete.json',complete)
  print('CAPTURE_FOCUS_REPORT\t'+str(root/'Yb_expanded_existing_controls_Rosina.pdf'));return
 cfg=json.loads((root/'run.json').read_text());jobs=json.loads((root/'jobs.json').read_text())
 if not a.sequencing_report:raise RuntimeError('--sequencing-report required')
 with a.sequencing_report.open() as h:records=list(csv.DictReader(h,delimiter='\t'))
 seq={r['sample']:r for r in records}
 if len(seq)!=len(records) or set(seq)!=set(t['sample'] for t in cfg['new_tasks']):raise RuntimeError('Sequencing manifest does not exactly match submitted tasks')
 # Hold before snapshotting; running transfers continue.
 if a.apply:
  for stage in ('download','align','call','report'):
   present=command(['squeue','--noheader','--jobs='+str(jobs[stage]),'--format=%i'],False).stdout.strip()
   if not present:continue
   r=command(['scontrol','hold',jobs[stage]],False)
   if r.returncode and not all('already finished' in line.lower() or 'no error' in line.lower() for line in r.stderr.splitlines() if line.strip()):raise RuntimeError('Could not hold '+stage+': '+r.stderr)
 snapshot=command(['squeue','--noheader','--array','--jobs='+','.join(str(jobs[k]) for k in ('download','align','call','report')),'--format=%i|%T']).stdout
 states=dict(line.strip().split('|',1) for line in snapshot.splitlines() if line.strip())
 plan=decisions(cfg['new_tasks'],seq,states,root,jobs)
 if sum(r['reason']=='confirmed_targeted_capture' for r in plan)!=338:raise RuntimeError('Expected 338 confirmed capture tasks; scope needs review')
 stage=root/'capture_focus'/time.strftime('%Y%m%dT%H%M%SZ',time.gmtime());stage.mkdir(parents=True,exist_ok=False)
 (stage/'queue_snapshot.tsv').write_text(snapshot)
 with (stage/'sample_decisions.tsv').open('w',newline='') as h:
  w=csv.DictWriter(h,fieldnames=list(plan[0]),delimiter='\t');w.writeheader();w.writerows(plan)
 keep={r['index'] for r in plan if r['keep']};exclude={r['index'] for r in plan if not r['keep']}
 narrowed=dict(cfg);narrowed['new_tasks']=[t for t in cfg['new_tasks'] if t['index'] in keep]
 callids={t['call_task'] for t in narrowed['new_tasks']};origids={t['task'] for t in cfg['original_tasks']}
 narrowed['tasks']=[t for t in cfg['tasks'] if t['task'] in callids|origids]
 save(stage/'run.json',narrowed)
 scope=dict(policy='capture_or_completed_WGS_or_complete_local_inputs; unfinished WGS cancelled',capture_tasks=338,retained_WGS=sum(r['keep'] and r['reason']!='confirmed_targeted_capture' for r in plan),deferred_WGS=len(exclude),deferred_samples=[r['sample'] for r in plan if not r['keep']],parent_run=str(root),original_run_manifest_unchanged=True,existing_data_deleted=False,selection_independent_of_genotype=True,sequencing_manifest=str(a.sequencing_report.resolve()),sequencing_manifest_sha256=hashlib.sha256(a.sequencing_report.read_bytes()).hexdigest())
 save(stage/'scope.json',scope)
 print('KEEP\t'+str(len(keep))+'\nDEFER\t'+str(len(exclude))+'\nPLAN\t'+str(stage))
 if not a.apply:return
 # Refuse to alter a report already running.
 if any(states.get(str(jobs['report'])+suffix) in {'RUNNING','COMPLETING'} for suffix in ('','_0')):raise RuntimeError('Original report running; arrays remain held')
 # Cancel unfinished WGS work and all superseded focus reports before cleanup.
 for previous in (root/'capture_focus').glob('*/jobs.json'):
  oldjob=json.loads(previous.read_text()).get('report')
  if oldjob:
   state=command(['squeue','--noheader','--jobs='+str(oldjob),'--format=%T'],False).stdout.strip()
   if state in {'RUNNING','COMPLETING'}:raise RuntimeError('Previous focus report running; arrays remain held')
   if state:command(['scancel',str(oldjob)])
 for k in ('download','align','call'):
  pending=[i for i in exclude if states.get(f'{jobs[k]}_{i}')]
  if pending:command(['scancel',str(jobs[k])+'_['+expression(pending)+']'])
 reportstate=states.get(str(jobs['report']),states.get(str(jobs['report'])+'_0',''))
 if reportstate in {'PENDING','PD'}:command(['scancel',jobs['report']])
 elif reportstate:raise RuntimeError('Unexpected original report state '+reportstate)
 # Wait for excluded writers to stop before removing only their raw download files.
 excluded_jobs={f'{jobs[k]}_{i}' for k in ('download','align','call') for i in exclude}
 for attempt in range(30):
  current=command(['squeue','--noheader','--array','--jobs='+','.join(str(jobs[k]) for k in ('download','align','call')),'--format=%i'],False).stdout.split()
  if not excluded_jobs.intersection(current):break
  time.sleep(1)
 else:raise RuntimeError('Excluded tasks still present; cleanup skipped and arrays remain held')
 cleaned=[]
 for t in cfg['new_tasks']:
  if t['index'] not in exclude:continue
  d=root/'raw'/t['sample']
  if (d/'complete.json').exists():raise RuntimeError('WGS finished during scope change; retain data and review '+t['sample'])
  for mate in (1,2):
   for suffix in ('','.part','.partial'):
    path=d/f"{t['sample']}_{mate}.fastq.gz{suffix}"
    if path.is_file() and not path.is_symlink():
     cleaned.append(dict(path=str(path),bytes=path.stat().st_size));path.unlink()
 save(stage/'removed_incomplete_WGS_downloads.json',cleaned)
 # A wrapper applies the filtered manifest without editing frozen code or original run.json.
 script=stage/'focus_Ybv_capture.py';script.write_bytes(Path(__file__).read_bytes())
 wrap=shlex.join(['micromamba','run','-n','gynandrocalling','python',str(script),'--run',str(root),'--report-config',str(stage/'run.json')])
 cmd=['sbatch','--parsable','--partition=nano','--cpus-per-task=1','--mem=8G','--time=00:30:00','--job-name=Ybv_capture_report','--dependency=afterany:'+str(jobs['call']),'--output='+str(stage/'report.%j.out'),'--error='+str(stage/'report.%j.err'),'--wrap='+wrap]
 reportjob=command(cmd).stdout.strip().split(';')[0]
 save(stage/'jobs.json',dict(original=jobs,report=reportjob,operation='capture_focus'))
 # Release only retained pending array elements. Deferred work stays cancelled/held.
 for k in ('download','align','call'):
  pending=[i for i in keep if states.get(f'{jobs[k]}_{i}') in {'PENDING','PD'}]
  if pending:command(['scontrol','release',str(jobs[k])+'_['+expression(pending)+']'])
 # Retain the scheduler's existing throttle; no parent-array update needed.
 print('REPORT\t'+reportjob+'\nPDF\t'+str(root/'Yb_expanded_existing_controls_Rosina.pdf'))
if __name__=='__main__':main()
