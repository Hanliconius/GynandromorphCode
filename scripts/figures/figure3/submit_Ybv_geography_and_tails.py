#!/usr/bin/env python3
"""Submit corrected geography/BA, then the phenotype-first high/low supplement."""
import json,shlex,shutil,subprocess,sys
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis')
SRC=ROOT/'Ybv_sorted_reviewed_BA_geography_v2_20261102'
OUT=ROOT/'Ybv_phenotype_first_high_low_20261102'
def main():
 repo=Path(__file__).resolve().parents[3]
 build=repo/'scripts/figures/figure3/rebuild_Ybv_sorted_reviewed_BA.py'
 tails=repo/'scripts/figures/figure3/add_Ybv_low_BA_genoplot.py'
 meta=repo/'data/figure3_review/Ybv_reviewed_metadata_20261102_geography_v2.tsv'
 for p in [build,tails,meta]:
  if not p.is_file():raise FileNotFoundError(p)
 if OUT.exists():raise RuntimeError('Final output already exists; inspect its jobs/logs before resubmitting: '+str(OUT))
 command=[sys.executable,str(build),'--submit','--metadata',str(meta),'--output',str(SRC)]
 if not SRC.exists():subprocess.run(command,check=True)
 elif not (SRC/'complete.json').is_file():
  if not (SRC/'jobs.json').is_file():raise RuntimeError('Source directory exists without a recorded job; inspect '+str(SRC))
  job=str(json.loads((SRC/'jobs.json').read_text())['plot'])
  q=subprocess.run(['squeue','-h','-j',job,'-o','%T'],capture_output=True,text=True)
  if not q.stdout.strip():
   states=subprocess.check_output(['sacct','-X','-n','-j',job,'--format=State','--parsable2'],text=True).strip()
   if 'COMPLETED' in states:raise RuntimeError('Source job reports completed but complete.json is absent; inspect source log')
   subprocess.run(command+['--retry'],check=True)
 render=[sys.executable,str(tails),'--submit','--order','phenotype-geography','--output',str(OUT)]
 if (SRC/'complete.json').is_file():
  subprocess.run(render,check=True);return
 job=str(json.loads((SRC/'jobs.json').read_text())['plot'])
 pending=SRC/'tails_submission_job.json'
 if pending.is_file():
  old=str(json.loads(pending.read_text())['job'])
  active=subprocess.run(['squeue','-h','-j',old],capture_output=True,text=True)
  if active.stdout.strip():print('FOLLOWUP_ALREADY_QUEUED\t'+old);return
 # Freeze the plotting submitter so later repo edits cannot change queued work.
 frozen=SRC/'high_low_submitter.py';shutil.copy2(tails,frozen);render[1]=str(frozen)
 follow=subprocess.check_output(['sbatch','--parsable','--partition=nano','--time=00:03:00','--mem=2G','--cpus-per-task=1',
  '--job-name=Ybv_next','--dependency=afterok:'+job,
  '--output='+str(SRC/'logs/tails_submit.%j.out'),'--error='+str(SRC/'logs/tails_submit.%j.err'),
  '--wrap='+shlex.join(render)],text=True).strip().split(';')[0]
 pending.write_text(json.dumps(dict(job=follow,source_job=job,output=str(OUT)),indent=2)+'\n')
 print('GEOGRAPHY_BA\t'+job+'\nAUTOMATIC_PLOT_SUBMISSION\t'+follow+'\nFINAL_OUTPUT\t'+str(OUT))
if __name__=='__main__':main()
