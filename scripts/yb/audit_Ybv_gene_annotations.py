#!/usr/bin/env python3
"""Find existing cortex/ivory annotations; export evidence without changing figures."""
import argparse,csv,gzip,json,os,re,shlex,shutil,subprocess,sys,tempfile
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
def run(a):
 a.output.mkdir(parents=True,exist_ok=False)
 rows=[];files=[];errors=[]
 # Existing annotation resources only. No sequence downloads or changes.
 for directory,dirs,names in os.walk(ROOT):
  dirs[:]=[d for d in dirs if d not in {'.git','__pycache__','node_modules'}]
  for name in names:
   if not name.lower().endswith(('.gff','.gff3','.gtf','.gff.gz','.gff3.gz','.gtf.gz')):continue
   path=Path(directory)/name
   if not path.is_file():continue
   files.append(str(path))
   try:
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt',errors='replace') as h:
     for number,line in enumerate(h,1):
      if line.startswith('#'):continue
      fields=line.rstrip().split('\t')
      if len(fields)!=9 or not re.search(r'\b(?:ivory|cortex)\b',fields[8],re.I):continue
      rows.append(dict(path=str(path),line=number,seqid=fields[0],feature=fields[2],start=fields[3],end=fields[4],strand=fields[6],attributes=fields[8]))
   except (OSError,EOFError) as e:errors.append(dict(path=str(path),error=str(e)))
 def table(name,data,cols):
  with (a.output/name).open('w') as h:
   w=csv.DictWriter(h,delimiter='\t',fieldnames=cols);w.writeheader();w.writerows(data)
 table('gene_annotation_candidates.tsv',rows,['path','line','seqid','feature','start','end','strand','attributes'])
 table('read_errors.tsv',errors,['path','error'])
 (a.output/'annotation_files_examined.txt').write_text('\n'.join(files)+'\n')
 (a.output/'complete.json').write_text(json.dumps(dict(status='ANNOTATION_EVIDENCE_COLLECTED',files=len(files),matches=len(rows),errors=len(errors),coordinates_not_yet_approved_for_display=True,reference='Require native NC_140403.1 or verified Hmel2.5 transfer; no constant coordinate offset'),indent=2)+'\n')
 if a.push_repo:
  repo=a.push_repo;env=dict(os.environ,GIT_SSH_COMMAND='ssh -o BatchMode=yes')
  def git(*args,cwd=repo):subprocess.run(['git','-C',str(cwd),*map(str,args)],check=True,env=env,stdout=subprocess.DEVNULL)
  with tempfile.TemporaryDirectory(prefix='yb-gene-publish.') as tmp:
   work=Path(tmp)/'repo';git('fetch','--quiet','origin','main');git('worktree','add','--quiet','--detach',work,'origin/main')
   try:
    dest=work/'results/yb_gene_annotation_audit'/a.output.name;dest.mkdir(parents=True)
    for p in a.output.iterdir():shutil.copy2(p,dest/p.name)
    git('add',str(dest.relative_to(work)),cwd=work);git('commit','--quiet','-m','Record existing cortex ivory annotation evidence',cwd=work);git('push','--quiet','git@github.com:Hanliconius/GynandromorphCode.git','HEAD:refs/heads/yb-gene-annotation-audit-v1',cwd=work)
   finally:git('worktree','remove','--force',work)
  print('BRANCH\tyb-gene-annotation-audit-v1')
 print('ANNOTATION_FILES\t'+str(len(files)));print('MATCHES\t'+str(len(rows)));print('OUTPUT\t'+str(a.output))
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'analysis/Ybv_gene_annotation_audit_v1');p.add_argument('--push-repo',type=Path);p.add_argument('--submit',action='store_true');a=p.parse_args()
 if not a.submit:return run(a)
 stage=ROOT/'analysis/Ybv_gene_annotation_code';stage.mkdir(parents=True,exist_ok=True);f=stage/'audit_Ybv_gene_annotations.py';shutil.copy2(__file__,f)
 cmd=[sys.executable,str(f),'--output',str(a.output)]
 if a.push_repo:cmd+=['--push-repo',str(a.push_repo)]
 print('JOB\t'+subprocess.check_output(['sbatch','--parsable','--partition=nano','--cpus-per-task=1','--mem=4G','--time=00:30:00','--job-name=Ybv_genes','--output='+str(stage/'audit.%j.out'),'--error='+str(stage/'audit.%j.err'),'--wrap='+shlex.join(cmd)],text=True).strip())
if __name__=='__main__':main()
