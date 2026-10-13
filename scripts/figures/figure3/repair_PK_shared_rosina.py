#!/usr/bin/env python3
"""Repair the already-submitted PK run, keeping downloads and using one Rosina alignment."""
import argparse
import gzip
import json
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pysam

HERE=Path(__file__).resolve().parent
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')

def submit_job(base, worker, jobs, stage, indices, cpu, mem, limit, dependency=None, partition='cpu', label=None):
    if not indices:
        return None
    args=['sbatch','--parsable',f'--partition={partition}',
          '--array='+','.join(map(str,indices)),f'--cpus-per-task={cpu}',
          f'--mem={mem}',f'--time={limit}',
          f'--job-name=PKfix_{label or stage}',
          f'--output={base}/logs/{label or stage}.fixed.%A_%a.out',
          f'--error={base}/logs/{label or stage}.fixed.%A_%a.err']
    if dependency:
        args.append('--dependency='+dependency)
    jid=subprocess.check_output(args+[str(worker),stage],text=True).strip().split(';')[0]
    jobs['fixed_'+(label or stage)]=jid
    (base/'jobs.json').write_text(json.dumps(jobs,indent=2)+'\n')
    print('SUBMITTED',label or stage,jid,flush=True)
    return jid

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--run',type=Path,required=True)
    args=parser.parse_args()
    base=args.run.resolve()
    if base.parent != ROOT/'analysis/pk_dark_pale_optix_yb':
        raise RuntimeError('Refusing to modify a run outside the PK optix/Yb project')
    cfgfile=base/'run.json'
    jobsfile=base/'jobs.json'
    cfg=json.loads(cfgfile.read_text())
    jobs=json.loads(jobsfile.read_text())
    if cfg.get('shared_rosina_optix') or any(k.startswith('fixed_') for k in jobs):
        raise RuntimeError('Already repaired: refusing duplicate Slurm jobs')
    if not jobs.get('download'):
        raise RuntimeError('Original download job not found; refusing unsafe dependency changes')
    if len(cfg['samples'])!=21:
        raise RuntimeError('Expected 21 PK samples')
    ros=next(r for r in cfg['references'] if r['name']=='HmRos')
    opt=next(r for r in cfg['references'] if r['name']=='Optix')
    if Path(ros['fasta']).resolve()!=Path(opt['fasta']).resolve():
        raise RuntimeError('Optix is not the same FASTA as Rosina; cannot share alignment')
    if opt.get('source_contig')!='NC_140406.1' or opt['contig']!='Hmr18' or ros['contig']!='NC_140403.1':
        raise RuntimeError('Unexpected chromosome alias or reference')
    if not (base/'optix_sites.bed').is_file() or not (base/'optix_alleles.tsv.gz').is_file():
        raise RuntimeError('Missing frozen Optix allele inputs')
    if not (base/'worker.sh').is_file() or not (base/'code/scripts/figures/figure3').is_dir():
        raise RuntimeError('Missing frozen worker or report code')

    # Alias only the region and allele coordinates, never the read alignment.
    chrom=opt['source_contig']
    bed=pd.read_csv(base/'optix_sites.bed',sep='\t',header=None)
    if set(bed[0].astype(str))!={'Hmr18'}:
        raise RuntimeError('Unexpected Optix BED chromosome')
    bed[0]=chrom
    bed.to_csv(base/'optix_sites_rosina.bed',sep='\t',index=False,header=False)
    with gzip.open(base/'optix_alleles.tsv.gz','rt') as inp, (base/'optix_alleles_rosina.tsv').open('w') as out:
        for line in inp:
            if not line.startswith('Hmr18\t'):
                raise RuntimeError('Unexpected Optix allele chromosome')
            out.write(chrom+line[len('Hmr18'):])
    target=base/'optix_alleles_rosina.tsv'
    pysam.tabix_compress(str(target),str(target)+'.gz',force=True)
    pysam.tabix_index(str(target)+'.gz',seq_col=0,start_col=1,end_col=1,zerobased=False,force=True)
    target.unlink()
    (base/'rosina_to_optix_chrom.tsv').write_text(chrom+'\tHmr18\n')

    # Retain all original task IDs (important for the already-submitted array).
    # New alignments are solely Hmel2.5 and HmRos; Optix profiles
    # consume an extracted BAM from the corresponding HmRos task.
    align_ids=[]
    reused_ids=[]
    optix_ids=[]
    for i,t in enumerate(cfg['tasks']):
        ref=t['reference']['name']
        if ref=='Optix':
            optix_ids.append(i)
            t['reuse']={'bam':str(base/'HmRos'/t['sample']/'optix.bam'),'dependency':None}
        elif t.get('reuse'):
            reused_ids.append(i)
        else:
            align_ids.append(i)
    if not align_ids or not optix_ids:
        raise RuntimeError('Unexpected task inventory')
    # An old Rosina alignment may have produced just the Yb chromosome.
    # Such a task must be realigned if Optix calls are needed for that sample.
    optix_samples={cfg['tasks'][i]['sample'] for i in optix_ids}
    for i in reused_ids[:]:
        t=cfg['tasks'][i]
        if t['reference']['name']=='HmRos' and t['sample'] in optix_samples:
            t['reuse']=None
            reused_ids.remove(i)
            align_ids.append(i)
    align_ids.sort()
    cfg['shared_rosina_optix']={'contig':chrom,'fasta':ros['fasta']}

    # Preflight before cancelling downstream work.
    engine=HERE/'extend_PK_optix_Yb.py'
    compile(engine.read_text(),str(engine),'exec')
    if 'shared_rosina_optix' not in engine.read_text():
        raise RuntimeError('Repository engine does not contain the shared-alignment fix')
    if (base/'run.before_shared_rosina.json').exists():
        raise RuntimeError('Repair backup already exists')

    shutil.copy2(cfgfile,base/'run.before_shared_rosina.json')
    shutil.copy2(jobsfile,base/'jobs.before_shared_rosina.json')
    shutil.copy2(base/'extend_PK_optix_Yb.py',base/'worker.before_shared_rosina.py')
    cfgfile.write_text(json.dumps(cfg,indent=2)+'\n')
    shutil.copy2(engine,base/'extend_PK_optix_Yb.py')

    # Deliberately DO NOT cancel the download array, or external downloads.
    for key in ('report','new_profiles','reused_profiles','align','prepare'):
        if jobs.get(key):
            subprocess.run(['scancel',str(jobs[key])],check=True)
            print('CANCEL_OBSOLETE',key,jobs[key],flush=True)

    prereqs=['afterok:'+str(jobs['download'])]
    for jid in set(cfg.get('existing_download_dependencies',{}).values()):
        prereqs.append('afterok:'+str(jid))
    worker=base/'worker.sh'
    alignment=submit_job(base,worker,jobs,'align',align_ids,16,'64G','2-00:00:00',','.join(prereqs),label='align')
    new_profiles=submit_job(base,worker,jobs,'profile',align_ids,2,'8G','04:00:00','aftercorr:'+alignment,label='yb_profiles')
    deps=[new_profiles]
    if reused_ids:
        dependencies=sorted({cfg['tasks'][i]['reuse']['dependency'] for i in reused_ids if cfg['tasks'][i]['reuse'].get('dependency')})
        dep='afterok:'+':'.join(dependencies) if dependencies else None
        deps.append(submit_job(base,worker,jobs,'profile',reused_ids,2,'8G','04:00:00',dep,label='reused_yb_profiles'))
    deps.append(submit_job(base,worker,jobs,'profile',optix_ids,2,'8G','04:00:00','afterok:'+alignment,label='optix_profiles'))
    submit_job(base,worker,jobs,'report',[0],1,'8G','00:30:00','afterok:'+':'.join(deps),partition='nano',label='report')
    print('DOWNLOAD_PRESERVED',jobs['download'])
    print('ALIGNMENTS',len(align_ids),'Optix_extra_alignments',0,'Optix_profiles',len(optix_ids))
    print('RUN_DIR',base)

if __name__=='__main__':
    main()
