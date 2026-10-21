#!/usr/bin/env python3
"""Unselected SNP heterozygosity: reuse Optix calls; extend Yb with frozen parameters.

All computation runs on Pegasus. No downloads, alignment, optimization, or edits
to previous outputs. This measures heterozygosity among SNPs, not per-base diversity.
"""
import argparse, hashlib, importlib.util, json, shlex, shutil, subprocess, sys
from pathlib import Path

ROOT = Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
PARENT = ROOT/'analysis/rosina_genomewide_joint_recall/20261014T223341Z_f956c8bc'
YB_FIT = ROOT/'analysis/original_502/reciprocal_L_R_mixture/optix_reciprocal_L_R_mixture.summary.tsv'
MODEL = '134_fit_reciprocal_L_R_mixture_genotypes.py'
REGIONS = {'Optix': ('NC_140406.1',1060000,1240000),
           'Ivory_Yb': ('NC_140403.1',1630000,1692000)}
INTERVALS = {'Optix': [('Band',1112651,1114314),('Ray1',1156746,1160220),
                      ('Ray2',1174547,1187833),('Dennis1',1215210,1216124),('Dennis2',1218557,1221494)],
             'Ivory_Yb': [('Yb(d)',1633041,1635646),('Yb(v)',1658079,1659072)]}

def save(path, value):
    Path(path).write_text(json.dumps(value,indent=2)+'\n')

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def readcfg(out):
    return json.loads((out/'run.json').read_text())

def command(args):
    subprocess.run(list(map(str,args)),check=True)

def snp(r):
    return len(r.ref)==1 and r.ref in 'ACGT' and len(r.alts or ())==1 and len(r.alts[0])==1 and r.alts[0] in 'ACGT'

def prepare(a):
    import pandas as pd, pysam
    cfg=json.loads((PARENT/'optix_consistent_v2/run.json').read_text())
    completion=json.loads((PARENT/'optix_consistent_v2/corrected_calls.complete.json').read_text())
    if completion['status']!='AUDITED_CALLS_AND_MIXTURE_COMPLETE':
        raise RuntimeError('Optix correction incomplete')
    names=cfg['focal']; samples=[r for r in cfg['samples'] if r['sample'] in names]
    if len(samples)!=4 or sorted(r['side'] for r in samples)!=['Left','Left','Right','Right']:
        raise RuntimeError('Expected two libraries per side')
    model=a.model_script or Path(__file__).resolve().parents[2]/'reciprocal_mixing'/MODEL
    pairs=PARENT/'optix_consistent_v2/reciprocal_mixture/optix_reciprocal_L_R_genotype_pairs.tsv.gz'
    yb=PARENT/'integration_Yb/Yb_chromosome.vcf.gz'
    needed=[pairs,yb,YB_FIT,model,Path(cfg['fasta']),*[Path(r['bam']) for r in samples]]
    for p in needed:
        if not p.is_file(): raise RuntimeError('Missing documented input: '+str(p))
    # Verify reference identity, library identity and indexed BAM availability.
    with pysam.FastaFile(cfg['fasta']) as fa:
        dictionary=dict(zip(fa.references,fa.lengths))
        op=pd.read_csv(pairs,sep='\t')
        if set(op.chrom)!={'NC_140406.1'}:
            raise RuntimeError('Optix corrected calls use an unexpected contig')
        for r in op.itertuples():
            if fa.fetch(r.chrom,int(r.pos)-1,int(r.pos)).upper()!=r.ref:
                raise RuntimeError('Optix REF/FASTA disagreement at '+str(r.pos))
    for r in samples:
        with pysam.AlignmentFile(r['bam'],'rb') as b:
            sm={x.get('SM') for x in b.header.to_dict().get('RG',[]) if x.get('SM')}
            if not b.has_index() or sm!={r['sample']} or dict(zip(b.references,b.lengths))!=dictionary:
                raise RuntimeError('Focal BAM reference/index/sample validation failed: '+r['sample'])
    t=pd.read_csv(YB_FIT,sep='\t').set_index('statistic')['value']
    keys=['estimated_c_left_from_right','estimated_c_right_from_left','estimated_error_rate',
          'ALT_mapping_bias_log_odds','beta_binomial_concentration']
    params={k:float(t[k]) for k in keys}
    if not (0<=params[keys[0]]<.5 and 0<=params[keys[1]]<.5 and
            0<params[keys[2]]<.5 and params[keys[4]]>0):
        raise RuntimeError('Invalid frozen mixture parameters')
    chrom,lo,hi=REGIONS['Ivory_Yb']; sites=[]
    with pysam.VariantFile(str(yb)) as vf,pysam.FastaFile(cfg['fasta']) as fa:
        for r in vf.fetch(chrom,lo-1,hi):
            if not snp(r): continue
            if fa.fetch(chrom,r.pos-1,r.pos).upper()!=r.ref:
                raise RuntimeError('Yb REF/FASTA disagreement at '+str(r.pos))
            sites.append(dict(chrom=chrom,pos=r.pos,ref=r.ref,alt=r.alts[0]))
    sites=pd.DataFrame(sites)
    if sites.empty or sites.pos.duplicated().any():
        raise RuntimeError('Empty or duplicate-position broad biallelic SNP set')
    out=a.output; out.mkdir(parents=True,exist_ok=False); (out/'logs').mkdir()
    shutil.copy2(__file__,out/'worker.py'); shutil.copy2(model,out/'model.py')
    sites.to_csv(out/'Ivory_Yb_unselected_sites.tsv',sep='\t',index=False)
    config=dict(output=str(out),samples=samples,fasta=cfg['fasta'],bcftools=cfg['bcftools'],
                optix_pairs=str(pairs),yb_vcf=str(yb),parameters=params,fit_summary=str(YB_FIT),
                policy='Yb exploratory inference with archived Yb parameters fixed; no refit; Optix existing inference unchanged',
                phenotype_or_BA_selection=False,downloads=False,alignments=False,
                hashes={str(p):sha(p) for p in [model,YB_FIT,pairs]},script_sha256=sha(__file__))
    save(out/'run.json',config)
    jobs={}
    def launch(stage,array=None,dependency=None):
        args=[sys.executable,str(out/'worker.py'),'--output',str(out),'--stage',stage]
        wrap=shlex.join(args)
        if array: wrap+=' --task "$SLURM_ARRAY_TASK_ID"'
        cmd=['sbatch','--parsable','--partition='+('cpu' if array else 'nano'),
             '--time='+('02:00:00' if array else '00:30:00'),'--cpus-per-task=1','--mem=8G',
             '--job-name=LR_Ho_'+stage,'--output='+str(out/'logs'/(stage+'.%A_%a.out')),
             '--error='+str(out/'logs'/(stage+'.%A_%a.err'))]
        if array: cmd+=['--array='+array]
        if dependency: cmd+=['--dependency=afterok:'+dependency]
        jid=subprocess.check_output(cmd+['--wrap='+wrap],text=True).strip().split(';')[0]
        jobs[stage]=jid;save(out/'jobs.json',jobs);print(stage.upper()+'\t'+jid,flush=True)
        return jid
    counts=launch('counts','0-3');launch('report',dependency=counts)
    print('OUTPUT\t'+str(out))

def counts(out,task):
    import pandas as pd,pysam
    cfg=readcfg(out);row=cfg['samples'][task];chrom,lo,hi=REGIONS['Ivory_Yb']
    target=out/(row['sample']+'.pileup.bcf')
    args=[cfg['bcftools'],'mpileup','-f',cfg['fasta'],'-r',f'{chrom}:{lo}-{hi}',
          '-q','30','-Q','20','-d','250','-a','FORMAT/AD,FORMAT/DP','-Ob','-o',target,row['bam']]
    command(args);command([cfg['bcftools'],'index','-f',target])
    sites=pd.read_csv(out/'Ivory_Yb_unselected_sites.tsv',sep='\t'); wanted=sites.set_index('pos')
    found={}
    with pysam.VariantFile(str(target)) as vf:
        if list(vf.header.samples)!=[row['sample']]:raise RuntimeError('Pileup sample identity mismatch')
        for r in vf:
            if r.pos not in wanted.index or r.info.get('INDEL',False):continue
            s=wanted.loc[r.pos]
            if r.ref!=s.ref or r.pos in found:raise RuntimeError('Pileup allele/duplicate error at '+str(r.pos))
            c=r.samples[row['sample']];ad=c.get('AD') or ()
            d={b:(n or 0) for b,n in zip(r.alleles,ad)}
            nr,na=d.get(s.ref,0),d.get(s.alt,0)
            found[r.pos]=(nr,na,sum(d.values())-nr-na)
    rows=[]
    for r in sites.itertuples():
        nr,na,other=found.get(r.pos,(0,0,0))
        rows.append(dict(pos=r.pos,ref_depth=nr,alt_depth=na,other_depth=other))
    pd.DataFrame(rows).to_csv(out/(row['sample']+'.counts.tsv.gz'),sep='\t',index=False)
    save(out/(row['sample']+'.complete.json'),dict(status='COMPLETE',command=list(map(str,args)),sites=len(rows)))

def infer_yb(out,cfg):
    import numpy as np,pandas as pd
    from scipy.special import logsumexp
    spec=importlib.util.spec_from_file_location('archived_model',out/'model.py')
    model=importlib.util.module_from_spec(spec);spec.loader.exec_module(model)
    t=pd.read_csv(out/'Ivory_Yb_unselected_sites.tsv',sep='\t')
    for side in ('Left','Right'):
        totals=np.zeros((len(t),3))
        for r in cfg['samples']:
            if r['side']!=side:continue
            json.loads((out/(r['sample']+'.complete.json')).read_text())
            c=pd.read_csv(out/(r['sample']+'.counts.tsv.gz'),sep='\t')
            if not np.array_equal(c.pos,t.pos):raise RuntimeError('Counts/site ordering mismatch')
            totals+=c[['ref_depth','alt_depth','other_depth']].to_numpy()
        t[side+'_total']=totals[:,0]+totals[:,1];t[side+'_alt']=totals[:,1]
        t[side+'_other']=totals[:,2]
    p=cfg['parameters']
    ll=model.state_log_likelihoods(t.Left_alt.to_numpy(),t.Left_total.to_numpy(),
           t.Right_alt.to_numpy(),t.Right_total.to_numpy(),p['estimated_c_left_from_right'],
           p['estimated_c_right_from_left'],p['estimated_error_rate'],
           p['ALT_mapping_bias_log_odds'],p['beta_binomial_concentration'])
    posterior=np.exp(ll-logsumexp(ll,axis=1)[:,None]);best=posterior.argmax(axis=1)
    t['pair_posterior']=posterior[np.arange(len(t)),best]
    t['mixture_left_GT']=model.GTL[best];t['mixture_right_GT']=model.GTR[best]
    t.to_csv(out/'Ivory_Yb_exploratory_fixed_parameter_pairs.tsv.gz',sep='\t',index=False)
    return t

def normalize(t,locus):
    import numpy as np,pandas as pd
    t=t.copy();chrom,lo,hi=REGIONS[locus]
    if 'chrom' not in t:raise RuntimeError('Missing chromosome')
    t=t.loc[t.chrom.eq(chrom)&t.pos.between(lo,hi)].copy()
    if t.empty or t.pos.duplicated().any():raise RuntimeError('Empty/duplicate focal SNPs: '+locus)
    t=t.loc[t.ref.isin(list('ACGT'))&t.alt.isin(list('ACGT'))].sort_values('pos').copy()
    valid={'0/0':0,'0/1':1,'1/0':1,'1/1':0}
    for side in ('left','right'):
        t[side+'_het']=t['mixture_'+side+'_GT'].astype(str).str.replace('|','/',regex=False).map(valid)
        col=side+'_DP' if side+'_DP' in t else side.title()+'_total'
        t[side+'_depth']=pd.to_numeric(t[col],errors='coerce')
    t['other_ok']=True
    if 'Left_other' in t:
        t['other_ok']=(t.Left_other/(t.Left_total+t.Left_other).clip(lower=1)<=.05)&(t.Right_other/(t.Right_total+t.Right_other).clip(lower=1)<=.05)
    return t

def windows(t,lo,hi,kind,size,cutoff,min_sites=20):
    import numpy as np,pandas as pd
    common=(t.pair_posterior>=cutoff)&(t.left_depth>=10)&(t.right_depth>=10)&t.other_ok&t.left_het.notna()&t.right_het.notna()
    good=t.loc[common]; pos=t.pos.to_numpy(); gp=good.pos.to_numpy(); rows=[]
    if kind=='bp':
        starts=range(lo,hi-size+2,max(1,size//4)); bounds=[(s,s+size-1) for s in starts]
    else:
        # SNP-count windows defined by eligible discovery sites before confidence filtering.
        bounds=[(int(pos[i]),int(pos[i+size-1])) for i in range(0,len(pos)-size+1,max(1,size//4))]
    for start,end in bounds:
        a=np.searchsorted(gp,start,side='left');b=np.searchsorted(gp,end,side='right')
        z=good.iloc[a:b];n=len(z);discovered=int(((pos>=start)&(pos<=end)).sum())
        l=float(z.left_het.mean()) if n>=min_sites else np.nan
        r=float(z.right_het.mean()) if n>=min_sites else np.nan
        rows.append(dict(start=start,end=end,mid=(start+end)/2,n_common=n,n_discovered=discovered,
                         fraction_accepted=n/discovered if discovered else np.nan,left_Ho=l,right_Ho=r,
                         right_minus_left=r-l,window_kind=kind,window_size=size,posterior_cutoff=cutoff))
    return pd.DataFrame(rows)

def draw(t,locus,out,pdf):
    import numpy as np,pandas as pd,matplotlib.pyplot as plt
    chrom,lo,hi=REGIONS[locus]; summaries=[]; widths=[1000,5000,10000] if locus=='Optix' else [1000,2000,5000]
    for kind,size in [('bp',w) for w in widths]+[('snps',50),('snps',100)]:
        fig,axes=plt.subplots(4,1,figsize=(11,7.8),sharex=True,gridspec_kw={'height_ratios':[.65,3,1.3,1.2]})
        fig.subplots_adjust(left=.10,right=.97,top=.88,bottom=.15,hspace=.15)
        ax=axes[0]
        for label,a,b in INTERVALS[locus]:
            ax.plot([a/1e6,b/1e6],[.35,.35],color='#555555',lw=4);ax.text((a+b)/2e6,.62,label,ha='center',fontsize=8)
            for other in axes[1:]:other.axvspan(a/1e6,b/1e6,color='#666666',alpha=.08)
        ax.set_ylim(0,1);ax.axis('off')
        for cutoff,style in [(.95,'-'),(.99,'--')]:
            w=windows(t,lo,hi,kind,size,cutoff);summaries.append(w.assign(locus=locus))
            if w.empty:continue
            for side,color in [('left','#1764AB'),('right','#D04B35')]:
                axes[1].plot(w.mid/1e6,w[side+'_Ho'],color=color,ls=style,lw=1.5 if cutoff==.95 else .9,
                             label=f'{side.title()} | posterior >= {cutoff}')
            axes[2].plot(w.mid/1e6,w.right_minus_left,color='#5B3A83',ls=style,lw=1.3)
            if cutoff==.95:
                axes[3].plot(w.mid/1e6,w.n_common,color='#333333',lw=1,label='Accepted on both sides')
                axes[3].plot(w.mid/1e6,w.n_discovered,color='#999999',lw=.8,label='Discovery SNPs')
        axes[1].set_ylim(0,1);axes[1].set_ylabel('Heterozygous fraction\namong callable SNPs')
        axes[1].legend(ncol=2,fontsize=8,loc='upper right')
        axes[2].axhline(0,color='#aaaaaa',lw=.7);axes[2].set_ylim(-1,1);axes[2].set_ylabel('Right - Left')
        axes[3].set_ylabel('SNP count');axes[3].legend(fontsize=7,loc='upper right')
        axes[3].set_xlabel(chrom+' position (Mb)');axes[3].set_xlim(lo/1e6,hi/1e6)
        for a in axes[1:]:a.spines[['top','right']].set_visible(False)
        title=f'{size:,} bp windows; {size//4:,} bp step' if kind=='bp' else f'{size} discovery SNPs per window; {size//4} SNP step'
        fig.suptitle(locus.replace('_','/')+' | Left/Right heterozygosity\n'+title,fontsize=13)
        method='Existing corrected Optix genotypes' if locus=='Optix' else 'Exploratory Yb genotypes: archived Yb parameters fixed; no refit'
        fig.text(.10,.055,method+'\nSame callable sites per side; depth >=10 each; at least 20 accepted SNPs/window. Gaps are masked windows.\nAll broad-window biallelic SNPs; no phenotype/BA filtering. Overlapping windows are not independent.',fontsize=8)
        pdf.savefig(fig);fig.savefig(out/(locus+'_'+kind+'_'+str(size)+'.png'),dpi=140);plt.close(fig)
    pd.concat(summaries,ignore_index=True).to_csv(out/(locus+'_window_statistics.tsv.gz'),sep='\t',index=False)

def report(out):
    import pandas as pd
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    cfg=readcfg(out);plt.rcParams.update({'pdf.fonttype':42,'font.size':9})
    optix=normalize(pd.read_csv(cfg['optix_pairs'],sep='\t'),'Optix')
    yb=normalize(infer_yb(out,cfg),'Ivory_Yb')
    pdfs=[]
    for locus,t in [('Optix',optix),('Ivory_Yb',yb)]:
        dest=out/(locus+'_LR_heterozygosity_versions.pdf')
        with PdfPages(dest) as pdf:draw(t,locus,out,pdf)
        pdfs.append(str(dest));print('PDF\t'+str(dest),flush=True)
    save(out/'complete.json',dict(status='COMPLETE',pdfs=pdfs,sites={'Optix':len(optix),'Ivory_Yb':len(yb)},
         scale='Heterozygous fraction among shared callable discovery SNPs; not per callable base',
         primary_posterior=.95,sensitivity_posterior=.99,min_AD_depth_each_side=10,min_SNPs_per_window=20,
         Yb_new_exploratory_calls=True,mixture_parameters_refitted=False,previous_outputs_changed=False,
         ascertainment='Existing broad comparative biallelic SNP discovery; no BA or phenotype selection',
         caution='Optix and Yb use their respective historical parameter fits; no direct cross-locus absolute diversity comparison'))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ROOT/'analysis/LR_window_heterozygosity_v1')
    p.add_argument('--model-script',type=Path);p.add_argument('--submit',action='store_true')
    p.add_argument('--stage',choices=['counts','report']);p.add_argument('--task',type=int)
    a=p.parse_args()
    if a.submit:prepare(a)
    elif a.stage=='counts':counts(a.output,a.task)
    elif a.stage=='report':report(a.output)
    else:p.error('Use --submit')

if __name__=='__main__':main()
