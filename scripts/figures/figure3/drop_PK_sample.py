#!/usr/bin/env python3
"""Exclude an unavailable PK and resubmit only pending corrected locus work."""
import argparse,json,shutil,subprocess
from pathlib import Path
from repair_PK_shared_rosina import submit_job
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
def validate(cfg,excluded):
    errors=[]
    for s in cfg['samples']:
        if s['sample'] in excluded:continue
        d=Path(cfg['fastq_root'])/s['sample'];marker=d/'download.complete.json'
        if not marker.is_file():errors.append(s['sample']+': no completion marker');continue
        record=json.loads(marker.read_text())
        if record.get('ENA',{}).get('sample_accession')!=s['biosample']:errors.append(s['sample']+': BioSample mismatch')
        if any(not (d/f"{s['sample']}_{mate}.fastq.gz").is_file() or (d/f"{s['sample']}_{mate}.fastq.gz").stat().st_size==0 for mate in (1,2)):errors.append(s['sample']+': missing FASTQ mate')
    if errors:raise RuntimeError('Retained samples not ready: '+', '.join(errors))
def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--exclude',required=True);a=p.parse_args();base=a.run.resolve()
    if base.parent!=ROOT/'analysis/pk_dark_pale_optix_yb':raise RuntimeError('Wrong project run')
    cfg=json.loads((base/'run.json').read_text());jobs=json.loads((base/'jobs.json').read_text())
    if not cfg.get('shared_rosina_optix'):raise RuntimeError('Requires corrected shared Rosina workflow')
    if cfg.get('excluded_samples'):raise RuntimeError('Exclusion already applied; inspect jobs.json before rerunning')
    if a.exclude!='SRR11243345' or len(cfg['samples'])!=21:raise RuntimeError('This recovery excludes only P4009 from the 21-PK run')
    excluded={a.exclude};validate(cfg,excluded)
    align=[i for i,t in enumerate(cfg['tasks']) if t['sample'] not in excluded and t['reference']['name']!='Optix' and not t.get('reuse')]
    yb=[i for i,t in enumerate(cfg['tasks']) if t['sample'] not in excluded and t['reference']['name']!='Optix']
    optix=[i for i,t in enumerate(cfg['tasks']) if t['sample'] not in excluded and t['reference']['name']=='Optix']
    if len(align)!=40 or len(yb)!=40 or len(optix)!=15:raise RuntimeError('Unexpected retained task inventory')
    keys=['fixed_align','fixed_yb_profiles','fixed_optix_profiles','fixed_report']
    for key in keys:
        jid=jobs[key];states=subprocess.check_output(['squeue','-h','-j',str(jid),'-o','%T'],text=True).splitlines()
        if not states or any(x.strip()!='PENDING' for x in states):raise RuntimeError(key+' is not wholly pending; refusing to cancel completed/running work')
    engine=base/'extend_PK_optix_Yb.py';source=engine.read_text()
    old="base=Path(cfg['run_dir']);samples=cfg['samples'];values="
    if source.count(old)!=1 or "'samples':21," not in source:raise RuntimeError('Frozen report shape changed')
    source=source.replace(old,"base=Path(cfg['run_dir']);samples=[s for s in cfg['samples'] if s['sample'] not in cfg.get('excluded_samples',[])];values=")
    source=source.replace("'samples':21,","'samples':len(samples),'excluded_samples':cfg.get('excluded_samples',[]),")
    compile(source,str(engine),'exec')
    for name in ['run.json','jobs.json','extend_PK_optix_Yb.py']:shutil.copy2(base/name,base/(name+'.before_P4009_exclusion'))
    cfg['excluded_samples']=[a.exclude];cfg['exclusion_reason']='User requested exclusion: P4009 FASTQ download incomplete; P2010 download succeeded'
    (base/'run.json').write_text(json.dumps(cfg,indent=2)+'\n');engine.write_text(source)
    (base/'excluded_samples.tsv').write_text('sample\tlabel\treason\nSRR11243345\tP4009\tFASTQ download incomplete; excluded by user\n')
    for key in keys:subprocess.run(['scancel',str(jobs[key])],check=True)
    worker=base/'worker.sh'
    alignment=submit_job(base,worker,jobs,'align',align,16,'64G','2-00:00:00',label='align')
    yp=submit_job(base,worker,jobs,'profile',yb,2,'8G','04:00:00','aftercorr:'+alignment,label='yb_profiles')
    op=submit_job(base,worker,jobs,'profile',optix,2,'8G','04:00:00','afterok:'+alignment,label='optix_profiles')
    submit_job(base,worker,jobs,'report',[0],1,'8G','00:30:00','afterok:'+yp+':'+op,partition='nano',label='report')
    print('SCOPE\t20 PK; 40 whole-genome alignments; 15 new optix profiles; P2010 retained')
    print('PDF\t'+str(base/'PK_dark_pale_optix_Yb_posthoc.pdf'))
if __name__=='__main__':main()
