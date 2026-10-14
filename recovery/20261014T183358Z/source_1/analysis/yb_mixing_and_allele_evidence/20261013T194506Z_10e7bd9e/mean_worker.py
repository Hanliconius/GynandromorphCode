#!/usr/bin/env python3
"""Submit matched coverage and race overlays on both existing Yb references.

Uses samtools depth rather than Mosdepth: 50-bp median depth / full-contig
mean depth, MAPQ >=30, paired overlaps counted once, default flag exclusions.
This is an explicit implementation of the normalization strategy, not an
exact reproduction of Mosdepth. Existing unfiltered bigWigs are also plotted.
"""
import argparse,json,os,shlex,shutil,subprocess,sys,uuid
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import pandas as pd
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
WINDOWS={'Hmel2.5':('Hmel215003o',1395000,1440000),
         'HmRos':('NC_140403.1',1620000,1670000)}
LEGACY={'Hmel2.5':[(1403500,1405500),(1429500,1430500)],
        'HmRos':[(1633041,1635646),(1658079,1659072)]}
def require(path):
    p=Path(path)
    if not p.is_absolute() or not p.is_file() or not p.stat().st_size:
        raise RuntimeError(f'Missing nonempty absolute input: {p}')
    return p
def submit(control):
    cfg=json.loads(require(control/'run.json').read_text())
    focal=Path(cfg['focal_run'])
    fc=json.loads(require(focal/'run.json').read_text())
    entries=[]
    for source,run in ((cfg,control),(fc,focal)):
        for task in source['tasks']:
            name=task['reference']['name']; sample=task['sample']
            directory=run/name/'coverage'
            bam=require(directory/(sample+'.dedup.bam'))
            require(str(bam)+'.bai')
            require(directory/(sample+'.complete.json'))
            bw=require(directory/(sample+'.dedup.50bp.bw'))
            entries.append({**task,'bam':str(bam),'bigwig':str(bw),
                'role':'control' if run==control else 'focal'})
    if len(entries)!=28:raise RuntimeError('Expected 20 control and 8 focal alignments')
    import matplotlib,pyBigWig
    samtools=str(require(cfg['tools']['samtools']))
    subprocess.run([samtools,'--version'],stdout=subprocess.DEVNULL,check=True)
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8]
    base=ROOT/'analysis/yb_control_overlays'/stamp
    (base/'logs').mkdir(parents=True)
    frozen=base/Path(__file__).name;shutil.copy2(__file__,frozen)
    config={'run_dir':str(base),'tasks':entries,'samtools':samtools,
            'control_run':str(control),'focal_run':str(focal),
            'method':'samtools depth -aa -s -q 0 -Q 30 -G 1796; 50-bp median / full-contig mean including zeros',
            'legacy_intervals':'previous Yb(d)/Yb(v) annotations; not verified Luca functional CRE boundaries'}
    configpath=base/'run.json';configpath.write_text(json.dumps(config,indent=2)+'\n')
    worker=base/'worker.sh'
    worker.write_text('#!/usr/bin/env bash\nset -euo pipefail\nexec '+
        shlex.join([sys.executable,str(frozen),'--config',str(configpath)])+
        ' --stage "$1" --task "$SLURM_ARRAY_TASK_ID"\n')
    jobs=[]
    try:
        for stage,array,cpus,mem,limit in [('coverage','0-27%28',2,'8G','04:00:00'),
                                        ('bundle','0-0%1',1,'8G','00:30:00')]:
            command=['sbatch','--parsable','--partition=cpu',f'--array={array}',
                f'--cpus-per-task={cpus}',f'--mem={mem}',f'--time={limit}',
                f'--job-name=Yb_{stage}',f'--chdir={base}',
                f'--output={base}/logs/{stage}.%A_%a.out',f'--error={base}/logs/{stage}.%A_%a.err']
            if jobs:command.append('--dependency=afterok:'+jobs[-1][1])
            job=subprocess.check_output(command+[str(worker),stage],text=True).strip().split(';')[0]
            jobs.append((stage,job))
            (base/'jobs.json').write_text(json.dumps(dict(jobs),indent=2)+'\n')
            print(f'{stage.upper()}_ARRAY\t{job}',flush=True)
    except Exception:
        for _,job in jobs:subprocess.run(['scancel',job],check=False)
        raise
    print(f'RUN_DIR\t{base}')
    print('SAMPLES\t10 controls plus 4 focal libraries per reference; no new alignments or downloads')
def coverage(config,task):
    e=config['tasks'][task];ref=e['reference']['name']
    contig,start,end=WINDOWS[ref]
    bam=require(e['bam']);tool=config['samtools']
    header=subprocess.check_output([tool,'view','-H',str(bam)],text=True)
    lengths={}
    for line in header.splitlines():
        if line.startswith('@SQ\t'):
            fields=dict(x.split(':',1) for x in line.split('\t')[1:])
            lengths[fields['SN']]=int(fields['LN'])
    if contig not in lengths or lengths[contig]<end:raise RuntimeError('Wrong/missing Yb contig in BAM')
    size=lengths[contig];depth=np.zeros(size,dtype=np.int32)
    command=[tool,'depth','-aa','-s','-q','0','-Q','30','-G','1796','-r',contig,str(bam)]
    proc=subprocess.Popen(command,stdout=subprocess.PIPE,text=True)
    count=0
    try:
        for line in proc.stdout:
            chrom,pos,d=line.rstrip().split('\t')
            pos=int(pos)
            if chrom!=contig or pos!=count+1:raise RuntimeError('Noncontinuous depth stream')
            depth[count]=int(d);count+=1
        if proc.wait()!=0:raise RuntimeError('samtools depth failed')
    finally:
        if proc.poll() is None:proc.kill();proc.wait()
    if count!=size:raise RuntimeError(f'Depth omitted positions: {count}/{size}')
    mean=float(depth.mean())
    if not np.isfinite(mean) or mean<=0:raise RuntimeError('Zero/nonfinite scaffold mean')
    starts=np.arange(((start-1)//50)*50,end,50)
    medians=np.array([np.mean(depth[s:min(s+50,size)]) for s in starts])
    import pyBigWig
    with pyBigWig.open(e['bigwig']) as bw:
        if bw.chroms().get(contig)!=size:raise RuntimeError('BigWig/BAM contig length mismatch')
        raw=[]
        for s in starts:
            vals=np.asarray(bw.values(contig,int(s),int(min(s+50,size)),numpy=True))
            raw.append(float(np.nan_to_num(vals,nan=0).mean()))
    directory=Path(config['run_dir'])/ref/'profiles';directory.mkdir(parents=True,exist_ok=True)
    out=directory/(e['sample']+'.tsv.gz');tmp=directory/(e['sample']+'.tmp.tsv.gz')
    pd.DataFrame({'position':starts+25.5,'mean_depth':medians,
        'normalized_depth':medians/mean,'raw_bigwig':raw}).to_csv(tmp,sep='\t',index=False)
    tmp.replace(out)
    (directory/(e['sample']+'.audit.json')).write_text(json.dumps({
        **e,'contig':contig,'contig_length':size,'mean_scaffold_depth':mean,
        'command':command,'window_size':50,'window_step':50,
        'mean_includes_zero_positions':True,'engine':'samtools depth; not Mosdepth',
        'partial_terminal_bin':'median of available contig bases'},indent=2)+'\n')
    print(f'PROFILE_COMPLETE\t{ref}\t{e["sample"]}\tmean_depth={mean:.4f}')
def plot(config,task,bundle=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'pdf.fonttype':42,'font.size':8})
    ref=list(WINDOWS)[task];contig,start,end=WINDOWS[ref]
    base=Path(config['run_dir'])/ref
    entries=[e for e in config['tasks'] if e['reference']['name']==ref]
    profiles={e['sample']:pd.read_csv(require(base/'profiles'/(e['sample']+'.tsv.gz')),sep='\t') for e in entries}
    groups={}
    for e in entries:
        if e['role']=='control':groups.setdefault((e['species'],e['subspecies'],e['yellow_bar']),[]).append(e['sample'])
    focal=[e for e in entries if e['role']=='focal']
    x=next(iter(profiles.values())).position.to_numpy()/1e6
    for sample,t in profiles.items():
        if not np.array_equal(t.position.to_numpy()/1e6,x):raise RuntimeError('Mismatched profile coordinates')
    outdir=base/'plots';outdir.mkdir(exist_ok=True)
    for column,suffix,ylabel in [('normalized_depth','normalized','50-bp median / scaffold mean'),
                                ('raw_bigwig','raw_log1p','ln(1 + raw bigWig coverage)')]:
        fig,axes=plt.subplots(len(groups),3,figsize=(16,10),sharex='col',sharey=True,
            gridspec_kw={'width_ratios':[2.4,1,1]},squeeze=False)
        fig.subplots_adjust(left=.08,right=.98,top=.88,bottom=.13,hspace=.32,wspace=.12)
        for row_axes,((species,race,phenotype),samples) in zip(axes,groups.items()):
            values=np.vstack([profiles[s][column].to_numpy() for s in samples])
            # Normalize each individual first, then average, giving equal weights.
            avg=values.mean(axis=0)
            transform=np.log1p if column=='raw_bigwig' else lambda a:a
            row_axes[0].set_title(f'H. {species} {race} | yellow bar {phenotype} | n={len(samples)}',loc='left',fontsize=9)
            for ax in row_axes:
                ax.plot(x,transform(avg),color='black',lw=1.1,label=f'{race}: n={len(samples)}')
                for e in focal:
                    side=e['side'];leg=e['leg']
                    color='#1764AB' if side.lower()=='left' else '#D04B35'
                    ax.plot(x,transform(profiles[e['sample']][column].to_numpy()),color=color,
                        ls='--' if str(leg)=='T3' else '-',lw=.8,alpha=.8,label=f'{side} {leg}')
                for lo,hi in LEGACY[ref]:ax.axvspan(lo/1e6,hi/1e6,color='#999999',alpha=.12,lw=0)
                ax.spines[['top','right']].set_visible(False)
                if column=='normalized_depth':
                    ax.set_ylim(0,2)
                    ax.set_yticks([0,.5,1,1.5,2])
                    ax.axhline(1,color='#999999',ls=':',lw=.6)
        handles,labels=axes[0,0].get_legend_handles_labels()
        labels[0]='Race mean (n=2)'
        fig.legend(handles,labels,ncol=5,fontsize=8,loc='upper center',bbox_to_anchor=(.55,.955))
        ranges=[(start,end)]+[(lo-1000,hi+1000) for lo,hi in LEGACY[ref]]
        for col,(lo,hi) in enumerate(ranges):
            axes[-1,col].set_xlim(lo/1e6,hi/1e6)
            axes[-1,col].set_xlabel(f'{contig} position (Mb)')
            if col:
                name='Yb(d)' if col==1 else 'Yb(v)'
                axes[0,col].set_title(f'{name} interval zoom\n1 kb flanks',fontsize=10)
            else:
                axes[0,col].set_title(axes[0,col].get_title()+'\nFull locus',loc='left',fontsize=9)
            axes[-1,col].ticklabel_format(axis='x',style='plain',useOffset=False)
        fig.supylabel(ylabel,x=.02)
        fig.suptitle(f'Yb coverage on {ref}: focal Left/Right overlaid with paper controls',fontsize=11)
        note=('MAPQ >=30; duplicate/secondary/QC-failed/unmapped reads excluded; paired overlaps counted once.\n'
              'Each individual normalized by its full-contig mean, then averaged within race (n=2). Samtools implementation; 50-bp nonoverlapping bins. Display clipped to 0–2; values outside this range are hidden.'
              if column=='normalized_depth' else
              'Existing unfiltered 50-bp bamCoverage tracks; raw control values averaged within race, then transformed with ln(1 + coverage).')
        fig.text(.13,.035,note+'\nIndividual control curves omitted for clarity. Grey shading: legacy Yb(d)/Yb(v) intervals; Luca functional CRE boundaries remain unverified. Focal legs are libraries from one individual.',fontsize=7)
        out=outdir/f'Yb_{ref}_controls_LR_{suffix}.pdf'
        fig.savefig(out);fig.savefig(out.with_suffix('.png'),dpi=200)
        if bundle is not None:bundle.savefig(fig)
        plt.close(fig)
        print(f'PDF\t{out}')
    (outdir/'overlay.audit.json').write_text(json.dumps({'reference':ref,'method':config['method'],
        'groups':[{'species':s,'race':r,'yellow_bar':p,'n':len(ids),'accessions':ids} for (s,r,p),ids in groups.items()],
        'focal_individuals':1,'focal_libraries':4,'luca_CRE_verified':False},indent=2)+'\n')
def bundle_plots(config):
    from matplotlib.backends.backend_pdf import PdfPages
    out=Path(config['run_dir'])/'Yb_both_references_controls_LR_normalized_0to2_with_interval_zooms_and_raw_log1p.pdf'
    tmp=out.with_name(out.stem+'.tmp.pdf')
    with PdfPages(tmp) as bundle:
        for task in range(2):plot(config,task,bundle=bundle)
        bundle.infodict().update(Title='Yb control overlays: Hmel2.5 and HmRos',
            Subject='Race means n=2 with focal Left/Right legs; normalized 0–2 and raw log coverage')
    tmp.replace(out)
    print(f'BUNDLE_PDF\\t{out}')
def replot(base):
    require(base/'run.json')
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    frozen=base/f'replot_{stamp}.py'
    shutil.copy2(__file__,frozen)
    command=shlex.join([sys.executable,str(frozen),'--config',str(base/'run.json'),'--stage','bundle','--task','0'])
    job=subprocess.check_output(['sbatch','--parsable','--partition=nano','--time=00:30:00','--mem=8G',
        '--job-name=YbBundle',f'--chdir={base}',f'--output={base}/logs/bundle.%j.out',
        f'--error={base}/logs/bundle.%j.err','--wrap',command],text=True).strip().split(';')[0]
    print(f'REPLOT_JOB\\t{job}')
    print(f'BUNDLE_PDF\\t{base}/Yb_both_references_controls_LR_normalized_0to2_with_interval_zooms_and_raw_log1p.pdf')
def main():
    p=argparse.ArgumentParser();p.add_argument('--control-run',type=Path);p.add_argument('--replot-run',type=Path)
    p.add_argument('--config',type=Path);p.add_argument('--stage',choices=['coverage','plot','bundle']);p.add_argument('--task',type=int)
    a=p.parse_args()
    if a.replot_run:replot(a.replot_run)
    elif a.control_run:submit(a.control_run)
    elif a.config and a.stage and a.task is not None:
        config=json.loads(a.config.read_text())
        if a.stage=='bundle':bundle_plots(config)
        else:(coverage if a.stage=='coverage' else plot)(config,a.task)
    else:p.error('Use --control-run or --config/--stage/--task')
if __name__=='__main__':
    try:main()
    except Exception as e:print(f'ERROR: {e}',file=sys.stderr);sys.exit(2)
