#!/usr/bin/env python3
"""Inventory publication code, validate fresh clones, prepare isolated recovery branches.

Default: collect and test locally on Pegasus. --push uploads recovery branches only.
No existing checkout is changed. No sequencing or full analysis pipeline is launched.
"""
import argparse
import ast
import csv
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

PROJECTS = {
 'GynandromorphCode': ('main','gynandrocalling',['GynandromorphCode_repo','GynandroCalling']),
 'ColiasGenomics_2026': ('main','colias_w_popgen',['ColiasGenomics_2026','colias_W_popgen']),
 'SynTangle': ('pegasus-genespace-pilot','syntangle_test',['SynTangle_pipeline_test','SynTangle_pegasus_test']),
 'MerianCode': ('main',None,['MerianCode']),
 'ivory-mir193_paper': ('main',None,[]),
 'genome-observatory-live': ('main',None,['genome-observatory-live']),
}
EXTENSIONS={'.py','.r','.sh','.bash','.slurm','.sbatch','.pl','.awk','.nf','.smk','.jl','.rs','.js','.ts'}
SKIP={'.git','.ssh','.cache','node_modules','site-packages','__pycache__','.venv','venv','envs','micromamba','miniconda3','software'}
SECRET=re.compile(rb'-----BEGIN (?:OPENSSH |RSA |EC )?PRIVATE KEY-----|(?:gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,})')

def command(args,cwd=None,timeout=1800,env=None):
    return subprocess.run([str(a) for a in args],cwd=cwd,env=env,capture_output=True,text=True,timeout=timeout)

def git(repo,*args):
    c=command(['git','-C',repo,*args])
    if c.returncode: raise RuntimeError(c.stderr.strip())
    return c.stdout

def digest(raw):
    return hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()

def candidate(path):
    return (path.suffix.lower() in EXTENSIONS
            or path.name in {'Snakefile','Makefile','Dockerfile','pyproject.toml','requirements.txt','environment.yml','environment.yaml','renv.lock'}
            or (path.suffix.lower() in {'.json','.toml','.yaml','.yml','.tsv'} and any(x in {'config','configs','tests','fixtures','.github'} for x in path.parts)))

def scan(root,issues):
    def error(e): issues.append({'path':str(root),'error':str(e)})
    for folder,dirs,files in os.walk(root,followlinks=False,onerror=error):
        dirs[:]=[d for d in dirs if d not in SKIP and (not d.startswith('.') or d=='.github') and not (Path(folder)/d).is_symlink()]
        for name in files:
            p=Path(folder)/name
            if p.is_symlink() or not candidate(p): continue
            try:
                if p.stat().st_size>2_000_000:
                    issues.append({'path':str(p),'error':'code candidate exceeds 2 MB; not copied'});continue
                raw=p.read_bytes()
                raw.decode('utf-8')
                if SECRET.search(raw):
                    issues.append({'path':str(p),'error':'credential-like content; not copied'});continue
                yield p,raw
            except (OSError,UnicodeError) as e: issues.append({'path':str(p),'error':str(e)})

def remote_project(root):
    c=command(['git','-C',root,'remote','get-url','origin'],timeout=15)
    if c.returncode: return None
    match=re.search(r'(?:github\.com[:/])(?:Hanliconius|DNAcrobatics)/([^/\s]+?)(?:\.git)?$',c.stdout.strip())
    return match.group(1) if match else None

def run_validation(repo,project,environment,out):
    records=[]
    def run(label,args,timeout=900):
        try:
            c=command(args,cwd=repo,timeout=timeout,env=dict(os.environ,PYTHONPATH=str(repo/'src')+os.pathsep+str(repo),MPLBACKEND='Agg'))
            (out/(label+'.log')).write_text(c.stdout+'\n'+c.stderr)
            records.append({'check':label,'status':'PASS' if c.returncode==0 else 'FAIL','exit_code':c.returncode,'command':args})
        except (OSError,subprocess.TimeoutExpired) as e:
            (out/(label+'.log')).write_text(str(e));records.append({'check':label,'status':'BLOCKED','detail':str(e)})
    prefix=['micromamba','run','-n',environment] if environment else []
    if environment and not shutil.which('micromamba'):
        records.append({'check':'environment','status':'BLOCKED','detail':'micromamba unavailable'});prefix=None
    syntax=[]
    for p in repo.rglob('*'):
        if '.git' in p.parts or not p.is_file(): continue
        try:
            if p.suffix=='.py': ast.parse(p.read_bytes(),filename=str(p))
            elif p.suffix in {'.sh','.slurm','.sbatch'}:
                c=command(['bash','-n',p],timeout=30)
                if c.returncode: syntax.append([str(p.relative_to(repo)),c.stderr])
        except Exception as e: syntax.append([str(p.relative_to(repo)),str(e)])
    records.append({'check':'python_shell_syntax','status':'FAIL' if syntax else 'PASS','errors':syntax,'python':sys.version})
    if prefix is not None:
        run('python_environment',prefix+['python','-c','import sys; print(sys.version)'])
        if (repo/'tests').is_dir(): run('unit_tests',prefix+['python','-m','unittest','discover','-s','tests','-v'])
        else: records.append({'check':'unit_tests','status':'NOT_AVAILABLE'})
        if project=='SynTangle':
            run('real_data_fixture_tests',prefix+['python','-m','unittest','discover','-s','validation/real_data','-p','test_*.py','-v'])
            run('cli_help',prefix+['python','-m','syntangle','--help'])
        if project=='MerianCode':
            run('cli_help',['python3','-m','meriancode','--help'])
            run('rank_example',['python3','-m','meriancode','rank','--genes','examples/genes.tsv','--output',str(out/'rank_example.tsv')])
        if project=='GynandromorphCode':
            run('analysis_dependencies',prefix+['python','-c','import numpy,pandas,scipy,matplotlib,pysam,pyBigWig; print("imports OK")'])
        if project=='ColiasGenomics_2026':
            run('analysis_dependencies',prefix+['python','-c','import numpy,pandas,matplotlib; print("imports OK")'])
    rfiles=[str(p) for p in repo.rglob('*') if p.is_file() and p.suffix.lower()=='.r' and '.git' not in p.parts]
    if rfiles:
        if shutil.which('Rscript'):
            run('r_parse',['Rscript','-e','for (p in commandArgs(TRUE)) parse(file=p)',*rfiles])
        else: records.append({'check':'r_parse','status':'BLOCKED','detail':'Rscript unavailable in current environment'})
    if project=='ivory-mir193_paper' and shutil.which('Rscript'):
        run('r_dependencies',['Rscript','-e','p<-c("Seurat","sctransform","SeuratDisk"); stopifnot(all(vapply(p,requireNamespace,logical(1),quietly=TRUE))); sessionInfo()'])
    records.append({'check':'full_data_reproduction','status':'NOT_RUN','detail':'Requires frozen input manifests and accepted output comparisons; fixture tests do not establish scientific reproduction.'})
    (out/'validation.json').write_text(json.dumps(records,indent=2))
    return records

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base',type=Path,default=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe'))
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--extra-root',action='append',default=[],metavar='PROJECT=/absolute/path')
    parser.add_argument('--project',action='append',choices=list(PROJECTS))
    parser.add_argument('--push',action='store_true',help='Upload only newly collected code to a new recovery branch per project')
    parser.add_argument('--skip-tests',action='store_true')
    args=parser.parse_args()
    out=args.out.expanduser().resolve();out.mkdir(parents=True,exist_ok=False)
    key=Path.home()/'.ssh/id_ed25519_github'
    if not key.is_file(): parser.error('Existing Pegasus GitHub key not found: '+str(key))
    os.environ['GIT_SSH_COMMAND']='ssh -i '+shlex.quote(str(key))+' -o IdentitiesOnly=yes -o BatchMode=yes -o StrictHostKeyChecking=yes'
    stamp=dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    roots={p:[args.base/x for x in cfg[2]] for p,cfg in PROJECTS.items()}
    # Discover only immediate project checkouts by their verified remote identity.
    for p in args.base.iterdir():
        if p.is_dir() and not p.is_symlink() and (p/'.git').exists():
            identity=remote_project(p)
            if identity in roots and p not in roots[identity]: roots[identity].append(p)
    for item in args.extra_root:
        name,path=item.split('=',1)
        if name not in roots or not Path(path).is_absolute(): parser.error('Invalid --extra-root')
        roots[name].append(Path(path))
    summary=[]
    for project in args.project or list(PROJECTS):
        project_out=out/project;project_out.mkdir();repo=project_out/'github';reports=project_out/'checks';reports.mkdir()
        branch,environment,_=PROJECTS[project]
        row={'project':project,'branch':branch,'issues':[]}
        try:
            url=f'git@github.com:Hanliconius/{project}.git'
            c=command(['git','clone','--quiet','--branch',branch,url,repo])
            if c.returncode: raise RuntimeError(c.stderr)
            row['commit']=git(repo,'rev-parse','HEAD').strip()
            (project_out/'remote_refs.tsv').write_text(git(repo,'for-each-ref','--format=%(refname)\t%(objectname)','refs/remotes/origin'))
            # All reachable remote history, including intake branches, counts as preserved.
            objects=git(repo,'rev-list','--objects','--remotes=origin')
            present={line.split(' ',1)[0] for line in objects.splitlines()}
            index=[]; recovered=[]
            for n,root in enumerate(roots[project]):
                if not root.exists():
                    row['issues'].append({'path':str(root),'error':'source root absent'});continue
                identity=remote_project(root) if (root/'.git').exists() else None
                if (root/'.git').exists() and identity!=project:
                    row['issues'].append({'path':str(root),'error':'wrong project remote; skipped'});continue
                if (root/'.git').exists():
                    (project_out/f'checkout_{n}_status.txt').write_text(git(root,'status','--short'))
                    (project_out/f'checkout_{n}_head.txt').write_text(git(root,'rev-parse','HEAD'))
                for path,raw in scan(root,row['issues']):
                    sha=digest(raw);known=sha in present
                    dest=Path('recovery')/stamp/f'source_{n}'/path.relative_to(root)
                    record={'source':str(path),'sha256':hashlib.sha256(raw).hexdigest(),'git_blob':sha,'status':'ALREADY_ON_GITHUB' if known else 'RECOVERED','destination':'' if known else str(dest)}
                    if not known:
                        target=repo/dest;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
                        target.chmod(0o755 if path.stat().st_mode&0o111 else 0o644)
                        recovered.append(record)
                    index.append(record)
            with (project_out/'inventory.tsv').open('w') as handle:
                writer=csv.DictWriter(handle,fieldnames=['source','sha256','git_blob','status','destination'],delimiter='\t');writer.writeheader();writer.writerows(index)
            row.update(scanned=len(index),recovered=len(recovered),source_roots=[str(p) for p in roots[project]])
            if not index: row['issues'].append({'error':'No cluster code inventoried; source location unresolved or absent'})
            # Validate the clean GitHub revision, excluding archival recovery candidates.
            if not args.skip_tests:
                validation_repo=project_out/'validation_checkout'
                c=command(['git','clone','--quiet','--no-hardlinks',repo,validation_repo])
                if c.returncode: raise RuntimeError(c.stderr)
                row['validation']=run_validation(validation_repo,project,environment,reports)
            if recovered:
                recovery=repo/'recovery'/stamp
                (recovery/'MANIFEST.json').write_text(json.dumps(recovered,indent=2))
                (recovery/'README.md').write_text('Recovered cluster code, not a validated or canonical pipeline. Preserve source paths and checksums. Colias Alba code remains quarantined and scientifically invalidated. Validate each selected workflow before promotion.\n')
                git(repo,'checkout','-b','pegasus-publication-recovery-'+stamp)
                git(repo,'add','--','recovery/'+stamp)
                git(repo,'-c','user.name=Joseph Hanly','-c','user.email=43384926+Hanliconius@users.noreply.github.com','commit','-m','Preserve missing Pegasus code with source checksums')
                row['recovery_commit']=git(repo,'rev-parse','HEAD').strip()
                row['recovery_branch']='pegasus-publication-recovery-'+stamp
                if args.push:
                    git(repo,'push','origin','HEAD:refs/heads/'+row['recovery_branch']);row['pushed']=True
                else: row['pushed']=False
            if project=='ivory-mir193_paper':
                c=command(['git','ls-remote','git@github.com:DNAcrobatics/ivory-mir193_paper.git','refs/heads/*','refs/tags/*'])
                (project_out/'publication_refs.txt').write_text(c.stdout+c.stderr)
                row['publication_copy_access']='PASS' if c.returncode==0 else 'BLOCKED'
                if c.returncode==0:
                    src={line.split()[1]:line.split()[0] for line in command(['git','ls-remote',url,'refs/heads/*','refs/tags/*']).stdout.splitlines()}
                    dst={line.split()[1]:line.split()[0] for line in c.stdout.splitlines()}
                    row['publication_ref_differences']=[{'ref':k,'source':src.get(k),'destination':dst.get(k)} for k in sorted(src.keys()|dst.keys()) if src.get(k)!=dst.get(k)]
        except Exception as e: row['error']=str(e)
        summary.append(row)
        (out/'summary.json').write_text(json.dumps(summary,indent=2))
        print(project+': '+('BLOCKED' if 'error' in row else f"{row.get('scanned',0)} candidates; {row.get('recovered',0)} recovered; pushed={row.get('pushed',False)}"),flush=True)
    print('Report: '+str(out/'summary.json'))
    return 1 if any('error' in r for r in summary) else 0

if __name__=='__main__': sys.exit(main())
