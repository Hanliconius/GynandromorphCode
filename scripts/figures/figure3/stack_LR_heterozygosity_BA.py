#!/usr/bin/env python3
"""Stack existing LR heterozygosity with independent comparative phenotype BA.

No calling or mixture fitting. Dennis scores are reused. Yb scores are new
exploratory scores on all cached SNPs and the frozen retained v7 cohort.
"""
import argparse,getpass,hashlib,importlib.util,json,shlex,shutil,subprocess,sys
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
RUN=ROOT/'analysis/rosina_genomewide_joint_recall/20261014T223341Z_f956c8bc'
HO=ROOT/'analysis/LR_window_heterozygosity_v1'
DEN=RUN/'Optix_Arnaud_review_nonPK_v2/merged'
YB=ROOT/'analysis/ybv_REF_ALT_rebuild_v1'
COHORT=ROOT/'analysis/ybv_short_ALT_geography_v7'


# Author exclusion, 2026-10-22: all calculations, irrespective of BA flags.
EXCLUDED_TAXON_PATTERN=r'numata|pardalinus'
EXCLUDED_SPECIMEN_IDS=set(['ERR1303880','ERR1303882','ERR1303884','ERR1303885','ERR260303','ERR260304'])
def excluded_specimen(record):
    import re
    text=' | '.join(str(v) for v in record.values()) if isinstance(record,dict) else str(record)
    return bool(re.search(EXCLUDED_TAXON_PATTERN,text,re.I)) or any(
        re.search(r'(?<![A-Za-z0-9])'+re.escape(s)+r'(?![A-Za-z0-9])',text) for s in EXCLUDED_SPECIMEN_IDS)
def require_allowed_specimens(records,context):
    bad=[r for r in records if excluded_specimen(r)]
    if bad:
        names=[str(r.get('sample',r.get('specimen_id',r))) if isinstance(r,dict) else str(r) for r in bad]
        raise RuntimeError('EXCLUDED_TAXA_IN_CALCULATION: '+context+': '+', '.join(names)+
            '; numata/pardalinus are excluded from all calculations. Rebuild affected cohorts, BA and frequencies before reusing results.')


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):Path(p).write_text(json.dumps(v,indent=2)+'\n')

def score(d,y,rule,min_rate=.8):
    import numpy as np
    d=np.asarray(d,float);y=np.asarray(y,int);called=np.isfinite(d)
    rates=[called[y==c].mean() if (y==c).any() else 0 for c in (0,1)]
    if min(rates)<min_rate:return dict(balanced_accuracy=np.nan,status='insufficient_callability',n_absent=int(sum(called&(y==0))),n_present=int(sum(called&(y==1))))
    dd=d[called];yy=y[called]
    if rule=='carrier':preds=[dd>0,dd<1]
    elif rule=='homozygote':preds=[dd==1,dd==0]
    else:raise ValueError(rule)
    scores=[.5*(pred[yy==1].mean()+(~pred[yy==0]).mean()) for pred in preds]
    tied=abs(scores[0]-scores[1])<1e-12
    return dict(balanced_accuracy=np.nan if tied else float(max(scores)),
                status='orientation_tie' if tied else 'scored',pattern_is_alt=None if tied else bool(scores[0]>scores[1]),
                ALT_orientation_BA=float(scores[0]),REF_orientation_BA=float(scores[1]),
                call_rate_absent=float(rates[0]),call_rate_present=float(rates[1]),
                n_absent=int(sum(called&(y==0))),n_present=int(sum(called&(y==1))))

def prepare(a):
    import pandas as pd
    files=[HO/'complete.json',HO/'worker.py',HO/'run.json',
           HO/'Ivory_Yb_exploratory_fixed_parameter_pairs.tsv.gz',
           RUN/'optix_consistent_v2/reciprocal_mixture/optix_reciprocal_L_R_genotype_pairs.tsv.gz',
           DEN/'new_comparative_SNP_scores.tsv.gz',DEN/'Dennis_BA90_intervals_sites.tsv',
           DEN/'Dennis_frequency_group_membership.json',YB/'SNPs.tsv',YB/'individual_ALT_dosages.tsv.gz',
           COHORT/'audit.json',COHORT/'membership_geography.tsv']
    for f in files:
        if not f.is_file():raise RuntimeError('Missing documented input: '+str(f))
    if json.loads((HO/'complete.json').read_text())['status']!='COMPLETE':raise RuntimeError('Ho incomplete')
    if json.loads((COHORT/'audit.json').read_text())['status']!='COMPLETE':raise RuntimeError('Yb cohort incomplete')
    out=a.output
    if a.retry:
        if not out.is_dir() or not (out/'run.json').is_file() or not (out/'jobs.json').is_file():
            raise RuntimeError('--retry requires an existing prepared run')
        if (out/'complete.json').exists():
            raise RuntimeError('Completed run will not be overwritten')
        previous=json.loads((out/'run.json').read_text())
        if 'Dennis_BA' not in previous or 'Yb_BA' not in previous:
            raise RuntimeError('Output directory is not owned by this workflow')
        jobs=json.loads((out/'jobs.json').read_text())
        active=set(subprocess.check_output(['squeue','--noheader','--user',getpass.getuser(),'--format=%A'],text=True).split())
        recorded={str(j).split(';')[0] for j in jobs.values()}
        if active & recorded:
            raise RuntimeError('Run has an active Slurm job; wait before retrying')
        (out/'logs').mkdir(exist_ok=True)
    else:
        out.mkdir(parents=True,exist_ok=False);(out/'logs').mkdir()
    shutil.copy2(__file__,out/'worker.py')
    save(out/'run.json',dict(input_hashes={str(f):sha(f) for f in files},
         script_sha256=sha(__file__),retry_failed_run=a.retry,SNP_calling_changed=False,mixture_fitting=False,window_bp=10000,step_bp=2500,
         Dennis_BA='Existing dominant carrier rule, original non-PK cohort',
         Yb_BA='New all-cached-SNP scores, frozen v7 membership; carrier and homozygote rules separately; >=80% called in each class',
         unscored_flanks='No expanded-cohort BA outside cached Yb SNP window; shaded, never set to zero'))
    job=subprocess.check_output(['sbatch','--parsable','--partition=nano','--time=00:30:00','--mem=8G',
         '--job-name=LR_Ho_BA','--output='+str(out/'logs/plot.%j.out'),
         '--error='+str(out/'logs/plot.%j.err'),'--wrap='+shlex.join([sys.executable,str(out/'worker.py'),'--output',str(out),'--run'])],text=True).strip().split(';')[0]
    save(out/'jobs.json',dict(plot=job));print('JOB\t'+job+'\nOUTPUT\t'+str(out))

def load_scores(out):
    import pandas as pd,numpy as np
    d=pd.read_csv(DEN/'new_comparative_SNP_scores.tsv.gz',sep='\t')
    d=d.loc[d.module.eq('Dennis')].copy()
    if d.empty or d.duplicated(['chrom','pos','ref','alt']).any():raise RuntimeError('Empty/duplicate Dennis scores')
    frozen=pd.read_csv(DEN/'Dennis_BA90_intervals_sites.tsv',sep='\t')
    check=frozen.merge(d,on=['chrom','pos','ref','alt'],suffixes=('_frozen','_all'),validate='one_to_one')
    if len(check)!=len(frozen) or not np.allclose(check.balanced_accuracy_frozen,check.balanced_accuracy_all,atol=1e-9):
        raise RuntimeError('Dennis broad scores disagree with frozen selected scores')
    groups=json.loads((DEN/'Dennis_frequency_group_membership.json').read_text())
    for group in groups:require_allowed_specimens(group['samples'],'frozen Dennis BA scores')
    dn={int(g['phenotype']):g['n'] for g in groups}
    if dn!={1:142,0:281}:raise RuntimeError('Dennis cohort changed')
    members=pd.read_csv(COHORT/'membership_geography.tsv',sep='\t').fillna('')
    require_allowed_specimens(members.to_dict('records'),'Yb(v) cached BA')
    members['v']=pd.to_numeric(members.v,errors='coerce')
    members=members.loc[members.v.isin([0,1])].copy()
    if members['sample'].duplicated().any():raise RuntimeError('Duplicate Yb individual')
    # PK and focal individuals must never enter phenotype statistics.
    words=members.astype(str).agg(' | '.join,axis=1).str.lower()
    bad=words.str.contains(r'piano|bwpk|pk_dark|pk_pale|gynandromorph|tn_left|tn_right|tn2l|tn2r|tn3l|tn3r',regex=True)
    if bad.any():raise RuntimeError('PK/focal candidate in frozen Yb cohort; resolve rather than silently recode')
    matrix=pd.read_csv(YB/'individual_ALT_dosages.tsv.gz',sep='\t',index_col=0)
    if not set(members['sample'])<=set(matrix.index) or matrix.index.duplicated().any():
        raise RuntimeError('Yb cached individual identities do not match retained cohort')
    matrix=matrix.loc[members['sample']]
    sites=pd.read_csv(YB/'SNPs.tsv',sep='\t');sites['pos1']=sites.pos1.astype(int)
    if sites.pos1.duplicated().any() or set(map(int,matrix.columns))!=set(sites.pos1):
        raise RuntimeError('Yb cached SNP columns mismatch')
    y=members.v.to_numpy(int);rows=[]
    for r in sites.itertuples():
        for rule in ['carrier','homozygote']:
            v=matrix[str(r.pos1)].to_numpy(float)
            if np.any(np.isfinite(v)&~np.isin(v,[0,.5,1])):raise RuntimeError('Non-diploid cached dosage')
            rows.append(dict(chrom=r.contig,pos=int(r.pos1),ref=r.ref,alt=r.alt,rule=rule,**score(v,y,rule)))
    ys=pd.DataFrame(rows);ys.to_csv(out/'Ybv_all_cached_SNP_BA.tsv.gz',sep='\t',index=False)
    members.to_csv(out/'Ybv_BA_membership.tsv',sep='\t',index=False)
    d.to_csv(out/'Dennis_existing_all_SNP_BA.tsv.gz',sep='\t',index=False)
    return d,ys,dict(Dennis=dn,Ybv={1:int(sum(y==1)),0:int(sum(y==0))})

def plot_page(pdf,h,s,locus,lo,hi,window,cohort,rule=None):
    import numpy as np,matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    fig,axes=plt.subplots(5,1,figsize=(10.5,8.4),sharex=True,
                         gridspec_kw={'height_ratios':[.6,2.2,1.8,1.15,1.]})
    fig.subplots_adjust(left=.11,right=.97,top=.89,bottom=.15,hspace=.14)
    intervals=[('D1',1215210,1216124),('D2',1218557,1221494)] if locus=='Dennis' else [('Yb(d)',1633041,1635646),('Yb(v)',1658079,1659072)]
    for label,start,end in intervals:
        if end<lo or start>hi:continue
        axes[0].plot([start/1e6,end/1e6],[.3,.3],color='#555555',lw=4)
        axes[0].text((start+end)/2e6,.65,label,ha='center',fontsize=8)
        for ax in axes[1:]:ax.axvspan(start/1e6,end/1e6,color='#555555',alpha=.07)
    axes[0].set_ylim(0,1);axes[0].axis('off')
    for cutoff,style in [(.95,'-'),(.99,'--')]:
        w=h.loc[h.window_kind.eq('bp')&h.window_size.eq(window)&np.isclose(h.posterior_cutoff,cutoff)]
        for side,col in [('left','#1764AB'),('right','#D04B35')]:
            axes[1].plot(w.mid/1e6,w[side+'_Ho'],c=col,ls=style,lw=1.3 if cutoff==.95 else .8,
                         label=side.title()+(f' (P >= {cutoff})'))
        axes[3].plot(w.mid/1e6,w.right_minus_left,color='#5B3A83',ls=style,lw=1.2)
        if cutoff==.95:
            axes[4].plot(w.mid/1e6,w.n_common,color='#333333',lw=1,label='Accepted both sides')
            axes[4].plot(w.mid/1e6,w.n_discovered,color='#999999',lw=.7,label='Discovery SNPs')
    x=s.pos/1e6;b=s.balanced_accuracy
    axes[2].scatter(x,b,s=7,color='#333333',alpha=.6,linewidths=0,rasterized=False)
    high=b>=.9;axes[2].scatter(x[high],b[high],s=9,color='#B15422',linewidths=0)
    axes[2].axhline(.5,color='#aaaaaa',lw=.6,ls=':');axes[2].axhline(.9,color='#B15422',lw=.6,ls='--')
    axes[2].set_ylabel('Phenotype BA');axes[2].set_ylim(0,1.02)
    if locus=='Yb(v)':
        a,z=int(s.pos.min()),int(s.pos.max())
        if lo<a:axes[2].axvspan(lo/1e6,a/1e6,facecolor='#eeeeee',hatch='///',edgecolor='#cccccc',lw=0)
        if hi>z:axes[2].axvspan(z/1e6,hi/1e6,facecolor='#eeeeee',hatch='///',edgecolor='#cccccc',lw=0)
        axes[2].text(.01,.08,'Hatched: no expanded-cohort BA data',transform=axes[2].transAxes,fontsize=7)
    axes[1].set_ylim(0,.6);axes[1].set_ylabel('Heterozygous fraction\namong callable SNPs')
    axes[1].legend(ncol=2,fontsize=7,loc='upper right')
    axes[3].axhline(0,color='#aaaaaa',lw=.6);axes[3].set_ylim(-.5,.5);axes[3].set_ylabel('Right - Left')
    axes[4].set_ylabel('SNP count');axes[4].legend(fontsize=7,loc='upper right')
    axes[4].set_xlabel('Rosina chromosome '+('18' if locus=='Dennis' else '15')+' position (Mb)')
    axes[4].set_xlim(lo/1e6,hi/1e6);axes[4].xaxis.set_major_locator(MaxNLocator(6))
    for ax in axes[1:]:ax.spines[['top','right']].set_visible(False)
    desc='existing carrier BA' if locus=='Dennis' else rule+' BA (new exploratory scores)'
    fig.suptitle(f'{locus} | heterozygosity and phenotype association\n{window:,} bp Ho windows / 2,500 bp slide; {desc}',fontsize=12)
    footer=f'BA cohort: phenotype+ n={cohort[1]}, phenotype- n={cohort[0]}; focal individual excluded. All eligible scored SNPs shown, no BA cutoff.\n'
    footer+='BA requires >=80% called in each phenotype class; tied allele orientations unscored. BA and Ho have different ascertainment.\n'
    footer+='Ho: same SNP denominator for both sides, depth >=10 each, >=20 SNPs/window; existing calls and mixture parameters unchanged.'
    fig.text(.11,.035,footer,fontsize=7)
    pdf.savefig(fig);plt.close(fig)

def ten_kb_statistics(out):
    """Re-window existing genotype pairs; never count reads or infer calls again."""
    import pandas as pd
    worker=HO/'worker.py'
    spec=importlib.util.spec_from_file_location('frozen_ho_windows',worker)
    old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    tables={}
    for locus,path in [
        ('Optix',RUN/'optix_consistent_v2/reciprocal_mixture/optix_reciprocal_L_R_genotype_pairs.tsv.gz'),
        ('Ivory_Yb',HO/'Ivory_Yb_exploratory_fixed_parameter_pairs.tsv.gz')]:
        t=old.normalize(pd.read_csv(path,sep='\t'),locus)
        _,lo,hi=old.REGIONS[locus]
        parts=[]
        for cutoff in [.95,.99]:
            w=old.windows(t,lo,hi,'bp',10000,cutoff)
            if w.empty or not (w.end-w.start+1).eq(10000).all() or not w.start.diff().dropna().eq(2500).all():
                raise RuntimeError('Unexpected 10kb/2.5kb window geometry')
            parts.append(w.assign(locus=locus))
        tables[locus]=pd.concat(parts,ignore_index=True)
        tables[locus].to_csv(out/(locus+'_10kb_step2500_window_statistics.tsv.gz'),sep='\t',index=False)
    return tables['Optix'],tables['Ivory_Yb']

def run(out):
    import pandas as pd,matplotlib
    matplotlib.use('Agg');matplotlib.rcParams.update({'pdf.fonttype':42,'font.size':8})
    from matplotlib.backends.backend_pdf import PdfPages
    d,y,cohorts=load_scores(out)
    dh,yh=ten_kb_statistics(out)
    paths=[]
    for locus,h,s in [('Dennis',dh,d),('Yb(v)',yh,y)]:
        dest=out/('Dennis_Ho_BA_stacked.pdf' if locus=='Dennis' else 'Ybv_Ho_BA_stacked.pdf')
        with PdfPages(dest) as pdf:
            rules=[None] if locus=='Dennis' else ['carrier','homozygote']
            for rule in rules:
                ss=s if rule is None else s.loc[s.rule.eq(rule)]
                bounds=[(1060000,1240000),(1208000,1225000)] if locus=='Dennis' else [(1630000,1692000),(1650000,1662000)]
                sizes=[10000]
                for lo,hi in bounds:
                    for size in sizes:plot_page(pdf,h,ss,locus,lo,hi,size,cohorts['Dennis' if locus=='Dennis' else 'Ybv'],rule)
        paths.append(str(dest));print('PDF\t'+str(dest),flush=True)
    save(out/'complete.json',dict(status='COMPLETE',PDFs=paths,cohorts=cohorts,
         Dennis_BA_reused=True,Ybv_BA_new_exploratory=True,calling_changed=False,mixture_refitted=False,
         window_bp=10000,step_bp=2500,existing_Yb_genotype_pairs_reused=True,
         Ybv_BA_span=[int(y.pos.min()),int(y.pos.max())],warning='Hatched Yb flanks have no BA data for this cohort'))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ROOT/'analysis/LR_Ho_BA_stacked_10kb_v2')
    p.add_argument('--submit',action='store_true');p.add_argument('--run',action='store_true')
    p.add_argument('--retry',action='store_true',help='Overwrite this failed run in place; refuses complete or active runs')
    a=p.parse_args()
    if a.retry and not a.submit:p.error('--retry requires --submit')
    if a.submit:prepare(a)
    elif a.run:run(a.output)
    else:p.error('Use --submit')

if __name__=='__main__':main()
