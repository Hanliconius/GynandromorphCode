#!/usr/bin/env python3
"""Retry transient PK ENA failures, validate all FASTQs, rewire pending alignments."""
import argparse, csv, json, subprocess, sys
from pathlib import Path

EXPECTED=(3,11,12,13,14)
def call(args):
    return subprocess.check_output([str(x) for x in args],text=True).strip()
def main():
    p=argparse.ArgumentParser()
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--verify',action='store_true')
    a=p.parse_args()
    root=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/pk_dark_pale_optix_yb')
    base=a.run.resolve()
    if base.parent!=root:raise RuntimeError('Wrong project directory')
    cfg=json.loads((base/'run.json').read_text())
    jobs=json.loads((base/'jobs.json').read_text())
    if len(cfg['samples'])!=21:raise RuntimeError('Unexpected sample count')
    if a.verify:
        missing=[]
        for s in cfg['samples']:
            acc=s['sample']
            folder=Path(cfg['fastq_root'])/acc
            marker=folder/'download.complete.json'
            if not marker.is_file():missing.append(acc+': no completion marker');continue
            data=json.loads(marker.read_text())
            if data['ENA']['sample_accession']!=s['biosample']:
                missing.append(acc+': mismatched BioSample')
                continue
            for mate in (1,2):
                fq=folder/f'{acc}_{mate}.fastq.gz'
                if not fq.is_file() or fq.stat().st_size==0:
                    missing.append(acc+': missing mate '+str(mate))
        if missing:
            raise RuntimeError('FASTQ validation failed: '+', '.join(missing))
        print('VALIDATED 21 paired FASTQ samples',flush=True)
        return
    if 'retry_download' in jobs:
        raise RuntimeError('Recovery already submitted')
    original=jobs.get('download')
    alignment=jobs.get('fixed_align')
    if original!='138462689' or alignment!='138462830':
        raise RuntimeError('Unexpected original job IDs; refusing to rewire')
    # Assert only these failed elements, and all other tasks still active or completed.
    states=call(['sacct','-X','-n','-j',original,'--format=JobID,State','--parsable2'])
    failed=set()
    for line in states.splitlines():
        cols=line.split('|')
        if len(cols)<2:continue
        jid,state=cols[:2]
        if jid.startswith(original+'_') and state.startswith(('FAILED','CANCELLED','TIMEOUT','OUT_OF_MEMORY')):
            suffix=jid.split('_',1)[1]
            if suffix.isdigit():failed.add(int(suffix))
    if failed!=set(EXPECTED):
        raise RuntimeError('Failed array elements changed: '+repr(sorted(failed)))
    state=call(['squeue','-h','-j',alignment,'-o','%T'])
    if not state or any(x not in ('PENDING','PD') for x in state.splitlines()):
        raise RuntimeError('Alignment job is no longer wholly pending')
    worker=base/'worker.sh'
    if not worker.is_file():raise RuntimeError('Frozen worker missing')
    # Same downloader, same run configuration and FASTQ storage. Only failed samples.
    retry=call(['sbatch','--parsable','--partition=cpu','--array='+','.join(map(str,EXPECTED)),
                '--cpus-per-task=1','--mem=4G','--time=12:00:00','--job-name=PKretry_download',
                f'--output={base}/logs/retry_download.%A_%a.out',
                f'--error={base}/logs/retry_download.%A_%a.err',str(worker),'download']).split(';')[0]
    jobs['retry_download']=retry
    (base/'jobs.json').write_text(json.dumps(jobs,indent=2)+'\n')
    # A gate executes only when both original array is finished (irrespective
    # of its known five failures) and all retries succeed.
    validator=Path(__file__).resolve()
    dep=f'afterany:{original},afterok:{retry}'
    gate=call(['sbatch','--parsable','--partition=nano','--cpus-per-task=1','--mem=2G',
               '--time=00:15:00','--job-name=PKverify_fastq','--dependency='+dep,
               f'--output={base}/logs/verify_fastq.%j.out',
               f'--error={base}/logs/verify_fastq.%j.err',
               '--wrap='+f'{sys.executable} {validator} --run {base} --verify']).split(';')[0]
    jobs['verify_download']=gate
    (base/'jobs.json').write_text(json.dumps(jobs,indent=2)+'\n')
    # Preserve existing 42-element alignment array and all its downstream jobs.
    call(['scontrol','update',f'JobId={alignment}',f'Dependency=afterok:{gate}'])
    print('RETRY_DOWNLOAD',retry)
    print('VALIDATION_GATE',gate)
    print('REWIRED_ALIGNMENT',alignment)
    print('PRESERVED_DOWNSTREAM',jobs.get('fixed_yb_profiles'),jobs.get('fixed_optix_profiles'),jobs.get('fixed_report'))
if __name__=='__main__':main()
