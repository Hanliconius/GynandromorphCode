#!/usr/bin/env python3
"""Combine compatible linear mean profiles; add 21 PK and four race controls."""
import argparse,json,sys,shlex,shutil,subprocess,importlib.util,hashlib
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import pandas as pd
REFS=['Hmel2.5','HmRos']
INTERVALS={'Hmel2.5':[(1403500,1405500),(1429500,1430500)],'HmRos':[(1633041,1635646),(1658079,1659072)]}
PRIMERS={'Hmel2.5':(1404112,1405075),'HmRos':(1633657,1635225)}
EXTRA=['ERR2298205','ERR2298206','ERR2298204','ERR2298214']
def require(p):
    p=Path(p)
    if not p.is_absolute() or not p.is_file():raise RuntimeError('Missing absolute file: '+str(p))
    return p
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def submit(a):
    old=json.loads(require(a.mean_run/'run.json').read_text());pk=json.loads(require(a.pk_run/'run.json').read_text());extra=json.loads(require(a.extra_run/'run.json').read_text())
    base=a.mean_run.parent.parent/'yb_all_samples_coverage'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ');(base/'logs').mkdir(parents=True)
    samples=[];tasks=[]
    for e in old['tasks']:
        ref=e['reference']['name'];s=e['sample'];profile=require(a.mean_run/ref/'profiles'/(s+'.tsv.gz'))
        samples.append(dict(reference=ref,sample=s,group='focal' if e['role']=='focal' else e['subspecies'],label=s,profile=str(profile),side=e.get('side'),leg=e.get('leg')))
    if len(pk['samples'])!=21:raise RuntimeError('Expected 21 manifest PKs')
    for ref in REFS:
        for e in pk['samples']:
            d=a.pk_run/ref/e['sample'];require(d/'profile.complete.json');profile=require(d/'coverage.tsv.gz');norm=require(d/'normalization.json')
            n=json.loads(norm.read_text())
            if n.get('QC')!='samtools depth -aa -s -q0 -Q30 -G1796; mean including zero-depth bases':raise RuntimeError('Unexpected PK coverage method')
            samples.append(dict(reference=ref,sample=e['sample'],group='PK '+e['phenotype'],label=e['label']+' '+e['sample']+' '+e['phenotype'],profile=str(profile),normalization=str(norm)))
    for e in extra['tasks']:
        if e['sample'] not in EXTRA:continue
        ref=e['reference']['name'];s=e['sample'];d=a.extra_run/ref/'coverage';require(d/(s+'.complete.json'))
        t={**e,'bam':str(require(d/(s+'.dedup.bam'))),'bigwig':str(require(d/(s+'.dedup.50bp.bw')))};tasks.append(t)
        samples.append(dict(reference=ref,sample=s,group=e['subspecies'],label=s,profile=str(base/ref/'profiles'/(s+'.tsv.gz'))))
    if len(tasks)!=8 or len(samples)!=78 or len({(e['reference'],e['sample']) for e in samples})!=78:raise RuntimeError('Expected 39 unique samples per reference')
    cfg={'run_dir':str(base),'tasks':tasks,'samtools':old['samtools'],'samples':samples,'sources':{'mean_run':str(a.mean_run),'pk_run':str(a.pk_run),'extra_run':str(a.extra_run)},'fit_summary':str(require(a.fit_summary))}
    config=base/'run.json';config.write_text(json.dumps(cfg,indent=2)+'\n');code=base/Path(__file__).name;shutil.copy2(__file__,code);shutil.copy2(require(a.mean_run/'mean_worker.py'),base/'mean_worker.py')
    common=['sbatch','--parsable','--partition=nano','--time=00:30:00','--mem=8G','--cpus-per-task=2',f'--chdir={base}'];cmd=[sys.executable,str(code),'--config',str(config)]
    job=subprocess.check_output(common+['--job-name=YbExtraMean','--array=0-7%8',f'--output={base}/logs/mean.%A_%a.out',f'--error={base}/logs/mean.%A_%a.err','--wrap='+shlex.join(cmd+['--stage','mean'])+' --task "$SLURM_ARRAY_TASK_ID"'],text=True).strip().split(';')[0]
    report=subprocess.check_output(common+['--job-name=YbAllCoverage','--dependency=afterok:'+job,f'--output={base}/logs/report.%j.out',f'--error={base}/logs/report.%j.err','--wrap='+shlex.join(cmd+['--stage','report'])],text=True).strip()
    (base/'jobs.json').write_text(json.dumps({'mean':job,'report':report},indent=2)+'\n');print('MEAN_ARRAY',job,'REPORT_JOB',report,'RUN_DIR',base,sep='\n')
def report(cfg):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    base=Path(cfg['run_dir']);groups=['rosina','melpomene','PK dark','PK pale','amaryllis','thelxinoe','florencia'];colors={'rosina':'#662A83','melpomene':'#555555','PK dark':'#663399','PK pale':'#CC70A7','amaryllis':'#5B8C45','thelxinoe':'#C79028','florencia':'#777777'}
    fit=pd.read_csv(cfg['fit_summary'],sep='\t').set_index('statistic').value;cl=float(fit['estimated_c_left_from_right']);cr=float(fit['estimated_c_right_from_left']);det=1-cl-cr
    if det<=0:raise RuntimeError('Invalid mixing inversion')
    audit={'sources':cfg['sources'],'samples':cfg['samples'],'profile_hashes':{e['profile']:sha(e['profile']) for e in cfg['samples']},'method':'50bp arithmetic mean / own full-contig mean including zeros; MAPQ30 BQ0 samtools -s -G1796','baseline':'per-bin mean of four melpomene on Hmel2.5; four rosina on HmRos; mask baseline below 0.25','mixing':{'left_from_right':cl,'right_from_left':cr,'source':cfg['fit_summary'],'fit_summary_sha256':sha(cfg['fit_summary']),'interpretation':'optix-fit transfer, equal underlying contig dosage and comparable mapping response assumed; no fitted coverage model'},'PK_dark_pale':'manifest phenotype; not inferred Yb genotype','intervals':'legacy annotations; separate Luca primer anchor','grid':'exact shared bin centres; no interpolation','panels':[]};rows=[]
    output=base/'Yb_all_21PK_4rosina_4melpomene_other_races_LR_coverage.pdf'
    with PdfPages(output) as pdf:
        for ref in REFS:
            entries=[e for e in cfg['samples'] if e['reference']==ref];tables={e['sample']:pd.read_csv(e['profile'],sep='\t').set_index('position') for e in entries}
            positions=sorted(set.intersection(*(set(t.index) for t in tables.values())));x=np.array(positions)
            if not positions:raise RuntimeError('No shared profile grid')
            for t in tables.values():
                if t.index.duplicated().any() or not {'mean_depth','normalized_depth'}.issubset(t.columns):raise RuntimeError('Need unique arithmetic mean profiles')
            y={s:t.reindex(x).normalized_depth.to_numpy() for s,t in tables.items()};matched='melpomene' if ref=='Hmel2.5' else 'rosina';ids=[e['sample'] for e in entries if e['group']==matched]
            if len(ids)!=4:raise RuntimeError('Expected four baseline controls')
            baseline=np.mean([y[s] for s in ids],axis=0);den=np.where(np.isfinite(baseline)&(baseline>=.25),baseline,np.nan)
            focal={e['sample']:e for e in entries if e['group']=='focal'};corrected={}
            for leg in ('T2','T3'):
                l=next(s for s,e in focal.items() if e['side']=='Left' and e['leg']==leg);r=next(s for s,e in focal.items() if e['side']=='Right' and e['leg']==leg)
                corrected[l]=((1-cr)*y[l]-cl*y[r])/det;corrected[r]=((1-cl)*y[r]-cr*y[l])/det
            for scale in ('contig normalized','reference-control relative'):
                divisor=np.ones(len(x)) if scale=='contig normalized' else den
                for adjusted in (False,True):
                    fig,axes=plt.subplots(len(groups),3,figsize=(16,15),sharey=True);windows=[(x.min(),x.max())]+[(lo-1000,hi+1000) for lo,hi in INTERVALS[ref]]
                    for i,group in enumerate(groups):
                        es=[e for e in entries if e['group']==group];arr=np.array([y[e['sample']]/divisor for e in es]);mean=np.mean(arr,axis=0)
                        for j,(lo,hi) in enumerate(windows):
                            ax=axes[i,j];mask=(x>=lo)&(x<=hi)
                            for curve in arr:ax.plot(x[mask]/1e6,curve[mask],color=colors[group],alpha=.23,lw=.5)
                            ax.plot(x[mask]/1e6,mean[mask],color=colors[group],lw=1.5,label=group+' mean')
                            for s,e in focal.items():ax.plot(x[mask]/1e6,((corrected[s] if adjusted else y[s])/divisor)[mask],color='#1764AB' if e['side']=='Left' else '#D9472D',ls='-' if e['leg']=='T2' else '--',lw=.9,label=e['side']+' '+e['leg'])
                            for aa,bb in INTERVALS[ref]:ax.axvspan(aa/1e6,bb/1e6,color='grey',alpha=.15)
                            aa,bb=PRIMERS[ref];ax.axvspan(aa/1e6,bb/1e6,color='#662A83',alpha=.08);ax.axhline(1,color='grey',lw=.5);ax.set_xlim(lo/1e6,hi/1e6);ax.set_ylim(0,2)
                            if i==0:ax.set_title(['Full plotted range','Yb(d): 1kb flanks','Yb(v): 1kb flanks'][j]);
                            if j==0:ax.set_ylabel(group+f' (n={len(es)})')
                            if i==len(groups)-1:ax.set_xlabel('Position (Mb)')
                    axes[0,2].legend(fontsize=7,loc='upper right');fig.suptitle(ref+' | '+scale+' | '+('LR adjusted using saved optix mixing fit' if adjusted else 'Observed LR')+' | 39 libraries/individuals')
                    fig.text(.08,.018,'Thin curves: individual controls/PKs; thick curve: group mean. Blue/red: focal Left/Right; solid T2, dashed T3. Display clipped to 0–2.\nGrey: legacy Yb intervals; purple: Luca primer span. Relative baseline: '+matched+' n=4; low baseline bins masked.\nAdjusted coverage is a linear inversion after own-contig scaling, assuming equal dosage/comparable mapping; negatives retained in TSV.',fontsize=8);fig.tight_layout(rect=[0,.07,1,.96]);pdf.savefig(fig);plt.close(fig)
                    audit['panels'].append({'reference':ref,'scale':scale,'LR_adjusted':adjusted,'baseline_masked_bins':int(np.isnan(den).sum())})
                for e in entries:
                    vals=y[e['sample']]/divisor
                    for k,(pos,value) in enumerate(zip(x,vals)):rows.append(dict(reference=ref,sample=e['sample'],group=e['group'],label=e['label'],position=pos,scale=scale,observed=value,adjusted_LR=corrected[e['sample']][k]/divisor[k] if e['sample'] in corrected else np.nan))
    pd.DataFrame(rows).to_csv(base/'coverage_plot_values.tsv.gz',sep='\t',index=False);(base/'audit.json').write_text(json.dumps(audit,indent=2)+'\n');print('BUNDLE_PDF',output,sep='\t')
def main():
    p=argparse.ArgumentParser();p.add_argument('--mean-run',type=Path);p.add_argument('--pk-run',type=Path);p.add_argument('--extra-run',type=Path);p.add_argument('--fit-summary',type=Path);p.add_argument('--config',type=Path);p.add_argument('--stage',choices=['mean','report']);p.add_argument('--task',type=int);a=p.parse_args()
    if not a.config:submit(a);return
    cfg=json.loads(a.config.read_text())
    if a.stage=='mean':
        spec=importlib.util.spec_from_file_location('mean_worker',a.config.parent/'mean_worker.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);m.coverage(cfg,a.task)
    else:report(cfg)
if __name__=='__main__':main()
