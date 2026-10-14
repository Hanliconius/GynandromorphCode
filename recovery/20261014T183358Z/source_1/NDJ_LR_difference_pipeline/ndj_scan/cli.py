"""Preflight, SLURM orchestration and local execution of the NDJ scan."""
import argparse
import importlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
import pysam
from .common import *

DEFAULT_PROJECT='/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling'
DEFAULT_SAMTOOLS='/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin/samtools'


def indexed(path, suffixes):
    for suffix in suffixes:
        p=Path(str(path)+suffix)
        if p.is_file(): return p
    # BAM index can alternatively replace .bam with .bai.
    if str(path).endswith('.bam') and Path(path).with_suffix('.bai').exists():
        return Path(path).with_suffix('.bai')
    raise ValueError(f'Missing index for {path}')


def configure(args):
    project=Path(args.project).resolve()
    root=project/'analysis/NDJ_control_comparison_v1'
    vcf=Path(args.vcf or root/'variant_calls/genomewide/NDJ_controls.joint_at_control_sites.bcf').resolve()
    reference=Path(args.reference or (project/'analysis/genomewide_mixing_input/metadata/reference_fasta.path').read_text().strip()).resolve()
    bam_l=Path(args.bam_left or root/'alignment/NDJ_left/NDJ_left.dedup.bam').resolve()
    bam_r=Path(args.bam_right or root/'alignment/NDJ_right/NDJ_right.dedup.bam').resolve()
    if bam_l==bam_r or args.left==args.right: raise ValueError('Difference scan requires distinct libraries and sample names')
    inputs=[vcf,indexed(vcf,['.csi','.tbi']),reference,Path(str(reference)+'.fai'),bam_l,indexed(bam_l,['.bai','.csi']),bam_r,indexed(bam_r,['.bai','.csi'])]
    with pysam.FastaFile(str(reference)) as fa:
        refdict=dict(zip(fa.references,fa.lengths))
    if args.contig:
        chromosomes=[dict(contig=c,label=c,background=c not in args.exclude_background) for c in args.contig]
    else:
        cohort=root/'metadata/cohort.json'
        inputs.append(cohort)
        original=json.loads(cohort.read_text())
        chromosomes=[dict(contig=x['contig'],label=str(x['chromosome']),background=int(x['chromosome'])!=21) for x in original['chromosomes']]
    if len({c['contig'] for c in chromosomes})!=len(chromosomes): raise ValueError('Duplicate contigs')
    if sum(c['background'] for c in chromosomes)<2: raise ValueError('Need at least two background chromosomes for held-out calibration')
    for c in chromosomes:
        c['length']=refdict[c['contig']]
        if c['contig'] in args.exclude_background: c['background']=False
    for bam,sample in [(bam_l,args.left),(bam_r,args.right)]:
        with pysam.AlignmentFile(str(bam),'rb',require_index=True) as b:
            if dict(zip(b.references,b.lengths))!=refdict: raise ValueError(f'Wrong BAM reference dictionary: {bam}')
            names={row.get('SM') for row in b.header.to_dict().get('RG',[])}
            if names!={sample}: raise ValueError(f'BAM SM mismatch: expected {sample}, saw {names}')
            next(b.fetch(chromosomes[0]['contig']),None)
    with pysam.VariantFile(str(vcf)) as b:
        if args.left not in b.header.samples or args.right not in b.header.samples: raise ValueError('NDJ sample name missing from BCF')
        ad=b.header.formats.get('AD')
        if ad is None or ad.number!='R' or ad.type!='Integer': raise ValueError('FORMAT/AD must be Number=R,Type=Integer')
        for c in chromosomes:
            if c['contig'] not in b.header.contigs or b.header.contigs[c['contig']].length!=c['length']:
                raise ValueError(f'VCF reference dictionary mismatch for {c["contig"]}')
            next(b.fetch(c['contig']),None)
    samtools=shutil.which(args.samtools) or args.samtools
    version=subprocess.run([samtools,'--version'],capture_output=True,text=True,check=True).stdout.splitlines()[0]
    subprocess.run([samtools,'quickcheck',str(bam_l),str(bam_r)],check=True)
    if args.repeats:
        inputs.extend([Path(args.repeats).resolve(),indexed(args.repeats,['.tbi','.csi'])])
        with pysam.TabixFile(args.repeats) as b:
            if not set(b.contigs)&set(refdict): raise ValueError('Repeat BED uses a different contig naming scheme')
    software = {}
    for package in ['numpy', 'pandas', 'scipy', 'matplotlib', 'pysam']:
        importlib.import_module(package)
        software[package] = importlib.metadata.version(package)
    cfg=dict(version=VERSION,project=str(project),vcf=str(vcf),reference=str(reference),bam_left=str(bam_l),bam_right=str(bam_r),
        left=args.left,right=args.right,chromosomes=chromosomes,samtools=str(Path(samtools).resolve()),
        repeats=str(Path(args.repeats).resolve()) if args.repeats else None,
        min_depth=args.min_depth,max_depth=args.max_depth,filter=args.filter,window=args.window,
        min_window_sites=args.min_window_sites,mapq=args.mapq,baseq=args.baseq,
        candidate_delta=args.candidate_delta,candidates=args.candidates,depth_shift=args.depth_shift,
        seed=20260919,background_sites=args.background_sites,software=software,samtools_version=version,
        inputs=[fingerprint(p) for p in inputs],created=datetime.now(timezone.utc).isoformat(),
        ascertainment='Existing control-discovered SNP panel; not de novo NDJ-private discovery',
        coordinate_convention='SNP pos 1-based; window start/end and BED 0-based half-open')
    print(f'PREFLIGHT_PASS\tchromosomes={len(chromosomes)}\tsamples={args.left},{args.right}\nINPUT\t{vcf}',flush=True)
    return cfg


def dispatch(stage, run, index=None):
    if stage=='scan':
        from .scan import scan
        scan(run,index)
    elif stage=='background':
        from .scan import calibrate
        calibrate(run)
    elif stage=='score':
        from .scan import score
        score(run,index)
    elif stage=='shortlist':
        from .scan import shortlist
        shortlist(run)
    elif stage=='qc':
        from .qc import qc
        qc(run,index)
    elif stage=='report':
        from .report import report
        report(run)


def submit_chain(run, cfg, args):
    log=run/'submissions.jsonl'
    if log.exists():
        ids=[json.loads(line)['job_id'] for line in log.read_text().splitlines()]
        check=subprocess.run(['squeue','-h','-j',','.join(ids),'-o','%i %T'],capture_output=True,text=True)
        # A purged set of IDs may yield nonzero; only tolerate explicit invalid-ID errors.
        if check.returncode and 'Invalid job id' not in check.stderr:
            raise RuntimeError('Could not check old jobs: '+check.stderr)
        if check.stdout.strip(): raise RuntimeError('Existing jobs remain active; no duplicate submission:\n'+check.stdout)
    dependency=None
    script=run/'code'/'run_ndj_scan.py'
    stages=[('scan',True,'cpu','24:00:00','16G'),('background',False,'nano','00:30:00','8G'),
            ('score',True,'cpu','12:00:00','16G'),('shortlist',False,'nano','00:30:00','4G'),
            ('qc',True,'cpu','04:00:00','4G'),('report',False,'cpu','02:00:00','8G')]
    for name,array_job,partition,time,mem in stages:
        wrapper=run/f'{name}.slurm'
        command=[sys.executable,str(script),'worker','--run',str(run),'--stage',name]
        body='#!/bin/bash\nset -euo pipefail\nexport OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1\n'
        body+='export MPLCONFIGDIR='+shlex.quote(str(run/'mplcache'))+'\n'
        body+=shlex.join(command)
        if array_job: body+=' --index "${SLURM_ARRAY_TASK_ID}"'
        wrapper.write_text(body+'\n')
        subprocess.run(['bash','-n',str(wrapper)],check=True)
        log_token = '%A_%a' if array_job else '%j'
        cmd=['sbatch','--parsable','--job-name=NDJ_'+name,'--partition='+partition,'--time='+time,'--mem='+mem,
            '--cpus-per-task=1','--output='+str(run/'logs'/f'{name}.{log_token}.out'),
            '--error='+str(run/'logs'/f'{name}.{log_token}.err')]
        if array_job: cmd+=['--array=0-'+str(len(cfg['chromosomes'])-1)+'%'+str(args.parallel)]
        if dependency: cmd+=['--dependency=afterok:'+dependency]
        cmd.append(str(wrapper))
        result=subprocess.run(cmd,capture_output=True,text=True)
        if result.returncode: raise RuntimeError(f'Submission failed at {name}: {result.stderr}\nSubmitted jobs were preserved; see {log}')
        job=result.stdout.strip().split(';')[0]
        if not re.fullmatch(r'\d+',job): raise RuntimeError('Unexpected sbatch response '+result.stdout)
        with log.open('a') as f:
            f.write(json.dumps(dict(stage=name,job_id=job,dependency=dependency,command=cmd))+'\n')
        dependency=job
        print(f'SUBMITTED\t{name}\t{job}',flush=True)
    print(f'RUN_DIRECTORY\t{run}\nEXPECTED_PDF\t{run}/figures/NDJ_LR_difference_narrative.pdf',flush=True)


def main():
    p=argparse.ArgumentParser(description='NDJ Left/Right exploratory AF and regional difference scan')
    p.add_argument('action',choices=['preflight','submit','local','worker'])
    p.add_argument('--project',default=DEFAULT_PROJECT)
    p.add_argument('--run',type=Path,help='New run directory, or existing run for --resume/worker')
    p.add_argument('--resume',action='store_true',help='Reuse immutable config and completed tasks; refuse if submitted jobs still active')
    p.add_argument('--stage',choices=['scan','background','score','shortlist','qc','report'])
    p.add_argument('--index',type=int)
    p.add_argument('--vcf');p.add_argument('--reference');p.add_argument('--bam-left');p.add_argument('--bam-right')
    p.add_argument('--left',default='NDJ_left');p.add_argument('--right',default='NDJ_right')
    p.add_argument('--samtools',default=DEFAULT_SAMTOOLS)
    p.add_argument('--contig',action='append');p.add_argument('--exclude-background',action='append',default=[])
    p.add_argument('--repeats',help='Optional bgzip/tabix BED, matching reference; absence remains NOT_ASSESSED')
    p.add_argument('--min-depth',type=int,default=10);p.add_argument('--max-depth',type=int)
    p.add_argument('--filter',choices=['all','pass'],default='all')
    p.add_argument('--window',type=int,default=100000);p.add_argument('--min-window-sites',type=int,default=30)
    p.add_argument('--mapq',type=int,default=20);p.add_argument('--baseq',type=int,default=20)
    p.add_argument('--candidate-delta',type=float,default=.15);p.add_argument('--candidates',type=int,default=200)
    p.add_argument('--depth-shift',type=float,default=.5);p.add_argument('--background-sites',type=int,default=100000)
    p.add_argument('--parallel',type=int,default=21)
    args=p.parse_args()
    if args.action=='worker':
        if not args.run or not args.stage: p.error('worker requires --run and --stage')
        if args.stage in ['scan','score','qc'] and args.index is None: p.error('array worker requires --index')
        dispatch(args.stage,args.run,args.index);return
    if min(args.min_depth,args.window,args.min_window_sites,args.parallel,args.background_sites)>0 and args.candidates>=2:
        pass
    else: p.error('Depth/window/parallel/background must be positive; candidates >=2')
    if not 0<=args.candidate_delta<=1 or args.depth_shift<=0: p.error('Invalid effect-size threshold')
    if args.resume:
        if not args.run: p.error('--resume requires --run')
        run,cfg=load_config(args.run)
        if (run/'COMPLETE.json').exists():
            print('ALREADY_COMPLETE\t'+str(run/'figures/NDJ_LR_difference_narrative.pdf'));return
    else:
        cfg=configure(args)
        if args.action=='preflight':return
        run=(args.run or Path(args.project)/'analysis/NDJ_LR_difference_scan_v1'/datetime.now(timezone.utc).strftime('run_%Y%m%dT%H%M%SZ')).resolve()
        run.mkdir(parents=True,exist_ok=False)
        (run/'logs').mkdir();(run/'chromosomes').mkdir()
        write_json(run/'config.json',cfg)
        source=Path(__file__).resolve().parent.parent
        shutil.copytree(source/'ndj_scan',run/'code'/'ndj_scan',ignore=shutil.ignore_patterns('__pycache__'))
        shutil.copy2(source/'run_ndj_scan.py',run/'code'/'run_ndj_scan.py')
        write_json(run/'code.sha256.json',{str(f.relative_to(run/'code')):hashlib.sha256(f.read_bytes()).hexdigest() for f in (run/'code').rglob('*.py')})
    if args.action=='submit': submit_chain(run,cfg,args)
    elif args.action=='local':
        for stage in ['scan','background','score','shortlist','qc','report']:
            for index in (range(len(cfg['chromosomes'])) if stage in ['scan','score','qc'] else [None]):
                dispatch(stage,run,index)
    print('RUN_DIRECTORY\t'+str(run),flush=True)


if __name__=='__main__':
    main()
