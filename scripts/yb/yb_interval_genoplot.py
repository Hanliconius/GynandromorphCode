#!/usr/bin/env python3
"""Control-selected SNP frequency tracks at legacy Yb(d)/Yb(v) intervals."""
import argparse,json,sys,subprocess
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
import pandas as pd
from assign_Yb_local_SNP_fragments import observations
REFS=['Hmel2.5','HmRos']
REGIONS={'Hmel2.5':('Hmel215003o',[(1403500,1405500),(1429500,1430500)]),'HmRos':('NC_140403.1',[(1633041,1635646),(1658079,1659072)])}
PRIMERS={'Hmel2.5':(1404112,1405075),'HmRos':(1633657,1635225)}
FOCAL=['TN2L_S1','TN2R_S2','TN3L_S3','TN3R_S4']
def analyse(base,task):
    cfg=json.loads((base/'run.json').read_text());ref=REFS[task//2];interval=['d','v'][task%2]
    chrom,ranges=REGIONS[ref];lo,hi=ranges[task%2];start,end=lo-1000,hi+1000
    entries=[e for e in cfg['tasks'] if e['reference']['name']==ref and (e['role']=='focal' or (e.get('species')=='melpomene' and e.get('subspecies') in ('rosina','melpomene')))]
    counts={};controls=defaultdict(list)
    for e in entries:
        path=Path(e['bam'])
        if not path.is_absolute() or not path.is_file():raise RuntimeError('Missing absolute BAM: '+str(path))
        c=defaultdict(Counter)
        for calls in observations(path,chrom,start,end).values():
            for pos,b in calls.items():c[pos][b]+=1
        counts[e['sample']]=c
        if e['role']=='control':controls[e['subspecies']].append(e['sample'])
    if any(len(controls[r])<2 for r in ('rosina','melpomene')):raise RuntimeError('Need two controls per race')
    rows=[]
    for pos in sorted(set().union(*(set(counts[s]) for ids in controls.values() for s in ids))):
        alleles={};valid=True
        for race in ('rosina','melpomene'):
            cons=[]
            for s in controls[race]:
                c=counts[s][pos];n=sum(c.values())
                if n<8:valid=False;break
                b,k=c.most_common(1)[0]
                if k/n<.9:valid=False;break
                cons.append(b)
            if not valid or len(set(cons))!=1:valid=False;break
            alleles[race]=cons[0]
        if not valid or alleles['rosina']==alleles['melpomene']:continue
        for e in entries:
            c=counts[e['sample']][pos];n=sum(c.values());r=c[alleles['rosina']];m=c[alleles['melpomene']]
            rows.append(dict(reference=ref,interval=interval,contig=chrom,pos1=pos+1,inside_interval=lo<=pos<hi,sample=e['sample'],role=e['role'],race=e.get('subspecies','focal'),rosina_base=alleles['rosina'],melpomene_base=alleles['melpomene'],depth=n,rosina_count=r,melpomene_count=m,other_count=n-r-m,rosina_frequency=r/(r+m) if r+m>=8 else np.nan))
    cols=['reference','interval','contig','pos1','inside_interval','sample','role','race','rosina_base','melpomene_base','depth','rosina_count','melpomene_count','other_count','rosina_frequency']
    pd.DataFrame(rows,columns=cols).to_csv(base/f'{ref}_Yb{interval}_SNP_counts.tsv.gz',sep='\t',index=False)
    print('INTERVAL_COMPLETE',ref,interval,len(rows)//len(entries),sep='\t')
def plot(base):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.patches import ConnectionPatch
    output=base/'Yb_d_v_control_selected_SNP_frequency_genoplots_both_references.pdf'
    corrections=[]
    with PdfPages(output) as pdf:
        for ref in REFS:
            for idx,interval in enumerate(('d','v')):
                df=pd.read_csv(base/f'{ref}_Yb{interval}_SNP_counts.tsv.gz',sep='\t');chrom,ranges=REGIONS[ref];lo,hi=ranges[idx]
                fig,axes=plt.subplots(3,1,figsize=(12,7),gridspec_kw={'height_ratios':[3,2,1]})
                fig.suptitle(f'{ref}: Yb({interval}) | control-selected SNP frequencies')
                if df.empty:
                    axes[0].text(.5,.5,'No SNPs passed the control filters',ha='center',transform=axes[0].transAxes)
                    axes[1].axis('off')
                else:
                    positions=np.sort(df.pos1.unique());p=df.pivot(index='sample',columns='pos1',values='rosina_frequency').reindex(columns=positions)
                    ctrl=[df.loc[(df.role=='control')&(df.race==r)].groupby('pos1').rosina_frequency.mean().reindex(positions).to_numpy() for r in ('rosina','melpomene')]
                    f=p.reindex(FOCAL).to_numpy();matrix=np.vstack(ctrl+[row for row in f]);cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#dddddd')
                    im=axes[0].imshow(matrix,aspect='auto',vmin=0,vmax=1,cmap=cmap,extent=(-.5,len(positions)-.5,5.5,-.5))
                    axes[0].set_yticks(range(6),['rosina controls','melpomene controls','Left T2','Right T2','Left T3','Right T3'])
                    corrected=f.copy()
                    for l,r in ((0,1),(2,3)):
                        corrected[l]=(.8*f[l]-.2*f[r])/.6;corrected[r]=(.8*f[r]-.2*f[l])/.6
                    axes[1].imshow(corrected,aspect='auto',vmin=0,vmax=1,cmap=cmap,extent=(-.5,len(positions)-.5,3.5,-.5))
                    axes[1].set_yticks(range(4),['Left T2','Right T2','Left T3','Right T3']);axes[1].set_title('Assumed symmetric 20% mixing correction; equal local dosage',fontsize=10)
                    for sample,values in zip(FOCAL,corrected):
                        for pos,value in zip(positions,values):corrections.append(dict(reference=ref,interval=interval,sample=sample,pos1=int(pos),corrected_frequency=value))
                    for ax in axes[:2]:ax.set_xticks([]);ax.set_xlim(-.5,len(positions)-.5)
                    axes[0].set_title('Observed; each column is one SNP, ordered along the locus',fontsize=10)
                    fig.colorbar(im,ax=axes[:2],label='Rosina-associated allele frequency',fraction=.025,pad=.02)
                    for j,pos in enumerate(positions):
                        fig.add_artist(ConnectionPatch(xyA=(j,3.5),coordsA=axes[1].transData,xyB=(pos,0),coordsB=axes[2].transData,color='grey',alpha=.18,lw=.4))
                ax=axes[2];ax.axhline(0,color='black',lw=1);ax.axvspan(lo,hi,color='grey',alpha=.25,label='Legacy interval')
                if interval=='d':
                    a,b=PRIMERS[ref];ax.plot([a,b],[.3,.3],color='purple',lw=4);ax.text((a+b)/2,.4,'Luca primer span',ha='center',fontsize=8)
                ax.set_xlim(lo-1000,hi+1000);ax.set_ylim(-.2,.8);ax.set_yticks([]);ax.ticklabel_format(axis='x',style='plain',useOffset=False);ax.set_xlabel(f'{chrom} position (bp)');ax.legend(loc='upper right',fontsize=8)
                fig.text(.08,.025,'MAPQ/baseQ >=30; overlapping mates counted once. Control depth >=8, consensus >=90%, agreement within race.\nGrey cells: fewer than 8 race-allele observations. Other alleles retained in counts TSV. Controls also selected the markers.\nIntervals are legacy annotations; primer span is an anchor, not a proven CRE boundary. Yb(v) markers are race-associated, not validated ventral markers.\nCorrected values outside 0–1 remain in TSV and are colour-saturated; focal libraries are from one individual.',fontsize=7)
                fig.subplots_adjust(left=.16,right=.86,bottom=.23,top=.88,hspace=.45);pdf.savefig(fig);plt.close(fig)
    pd.DataFrame(corrections).to_csv(base/'assumed20pct_corrected_SNP_frequencies.tsv.gz',sep='\t',index=False)
    print('BUNDLE_PDF',output,sep='\t')
def main():
    p=argparse.ArgumentParser();p.add_argument('--source-run',type=Path);p.add_argument('--config',type=Path);p.add_argument('--stage',choices=['counts','plot']);p.add_argument('--task',type=int);a=p.parse_args()
    if a.source_run:
        from datetime import datetime,timezone
        import shutil
        base=a.source_run/'interval_SNP_genoplot'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ');(base/'logs').mkdir(parents=True)
        cfg=json.loads((a.source_run/'run.json').read_text());(base/'run.json').write_text(json.dumps(cfg,indent=2))
        for name in (Path(__file__).name,'assign_Yb_local_SNP_fragments.py'):shutil.copy2(Path(__file__).parent/name,base/name)
        cmd=[sys.executable,str(base/Path(__file__).name),'--config',str(base/'run.json')]
        common=['sbatch','--parsable','--partition=nano','--time=00:30:00','--mem=8G','--cpus-per-task=1',f'--chdir={base}']
        import shlex
        job=subprocess.check_output(common+['--job-name=YbSNPinterval','--array=0-3%4',f'--output={base}/logs/counts.%A_%a.out',f'--error={base}/logs/counts.%A_%a.err','--wrap='+shlex.join(cmd+['--stage','counts'])+' --task "$SLURM_ARRAY_TASK_ID"'],text=True).strip().split(';')[0]
        plotjob=subprocess.check_output(common+['--job-name=YbSNPgenoplot','--dependency=afterok:'+job,f'--output={base}/logs/plot.%j.out',f'--error={base}/logs/plot.%j.err','--wrap='+shlex.join(cmd+['--stage','plot'])],text=True).strip()
        print('COUNTS_ARRAY',job,'PLOT_JOB',plotjob,'RUN_DIR',base,sep='\n')
    elif a.stage=='counts':analyse(a.config.parent,a.task)
    else:plot(a.config.parent)
if __name__=='__main__':main()
