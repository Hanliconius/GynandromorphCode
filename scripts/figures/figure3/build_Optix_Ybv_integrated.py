#!/usr/bin/env python3
"""Full-span Ho, gene annotation and frozen publication SNP panels; no calling."""
import argparse,getpass,gzip,hashlib,importlib.util,json,re,shlex,shutil,subprocess,sys
from pathlib import Path
from urllib.parse import unquote
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
RUN=ROOT/'analysis/rosina_genomewide_joint_recall/20261014T223341Z_f956c8bc'
DEN=RUN/'Optix_Arnaud_review_nonPK_v2/merged'
HO=ROOT/'analysis/LR_window_heterozygosity_v1'
DF=ROOT/'analysis/Dennis_figure_set_v1'
YF=ROOT/'analysis/Ybv_figure_set_v1'
V7=ROOT/'analysis/ybv_short_ALT_geography_v7'
OLD=ROOT/'analysis/yb_completed_reads_recovery_v1/frozen/HmRos_Ybv_unchanged_focal_calls.tsv.gz'
GFF=ROOT/'analysis/genomewide_mixing_input/reference/GCF_056956725.1/GCF_056956725.1_ASM5695672v1_genomic.gff.gz'
SPANS={'Dennis':('NC_140406.1',1060000,1240000),'Yb(v)':('NC_140403.1',1630000,1692000)}
INTERVALS={'Dennis':[('B',1111500,1117000),('R1',1156000,1161000),('R2',1174000,1189000),('D1',1215000,1217000),('D2',1218000,1223000)],
           'Yb(v)':[('Yb(d)',1633041,1635646),('Yb(v)',1658079,1659072)]}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):Path(p).write_text(json.dumps(v,indent=2)+'\n')
def boolean(x):
    if str(x).lower() in ('true','1','1.0'):return True
    if str(x).lower() in ('false','0','0.0'):return False
    raise ValueError('Invalid Boolean '+str(x))

def prepare(a):
    features=Path(__file__).resolve().parents[2]/'yb/Ybv_gene_features.json'
    ba_helper=Path(__file__).with_name('stack_LR_heterozygosity_BA.py')
    inputs=[ba_helper,DEN/'new_comparative_SNP_scores.tsv.gz',DEN/'Dennis_BA90_intervals_sites.tsv',
        DEN/'Dennis_frequency_group_membership.json',ROOT/'analysis/ybv_REF_ALT_rebuild_v1/SNPs.tsv',
        ROOT/'analysis/ybv_REF_ALT_rebuild_v1/individual_ALT_dosages.tsv.gz',HO/'complete.json',HO/'worker.py',HO/'Ivory_Yb_exploratory_fixed_parameter_pairs.tsv.gz',
        RUN/'optix_consistent_v2/reciprocal_mixture/optix_reciprocal_L_R_genotype_pairs.tsv.gz',
        DF/'all_individual_associated_dosages.tsv.gz',DF/'individual_membership.tsv',
        DF/'frozen_markers_and_BA.tsv',DEN/'Dennis_BA90_intervals_frequencies.tsv',
        YF/'individual_ALT_dosages.tsv.gz',V7/'membership_geography.tsv',
        V7/'region_race_ALT_frequencies.tsv',OLD,GFF,features]
    for p in inputs:
        if not p.is_file():raise RuntimeError('Missing documented input: '+str(p))
    if json.loads((HO/'complete.json').read_text())['status']!='COMPLETE':raise RuntimeError('Ho not complete')
    out=a.output
    if out.exists():
        if not a.overwrite:raise RuntimeError('Use --overwrite to update this figure in place')
        prior=json.loads((out/'run.json').read_text())
        if prior.get('workflow')!='Optix_Ybv_integrated':raise RuntimeError('Different workflow owns output')
        jobs=json.loads((out/'jobs.json').read_text()) if (out/'jobs.json').exists() else {}
        active=set(subprocess.check_output(['squeue','-h','-u',getpass.getuser(),'-o','%A'],text=True).split())
        if active & set(map(str,jobs.values())):raise RuntimeError('Figure job is still active')
    out.mkdir(parents=True,exist_ok=True);(out/'logs').mkdir(exist_ok=True)
    (out/'complete.json').unlink(missing_ok=True)
    shutil.copy2(__file__,out/'worker.py');shutil.copy2(features,out/'Ybv_gene_features.json');shutil.copy2(ba_helper,out/'ba_worker.py')
    save(out/'run.json',dict(workflow='Optix_Ybv_integrated',input_hashes={str(p):sha(p) for p in inputs if p!=GFF},
        native_GFF=dict(path=str(GFF),bytes=GFF.stat().st_size),spans=SPANS,window_bp=10000,slide_bp=2500,
        genotypes_and_cohorts_changed=False,mixture_refitted=False,code_sha256=sha(__file__)))
    job=subprocess.check_output(['sbatch','--parsable','--partition=nano','--time=00:30:00','--mem=8G',
        '--job-name=OxYb_panel','--output='+str(out/'logs/plot.%j.out'),'--error='+str(out/'logs/plot.%j.err'),
        '--wrap='+shlex.join([sys.executable,str(out/'worker.py'),'--output',str(out),'--run'])],text=True).strip().split(';')[0]
    save(out/'jobs.json',dict(plot=job));print('JOB\t'+job+'\nOUTPUT\t'+str(out))

def publication_data():
    import numpy as np,pandas as pd
    ds=pd.read_csv(DF/'frozen_markers_and_BA.tsv',sep='\t')
    if len(ds)!=24:raise RuntimeError('Dennis marker count changed')
    keys=[':'.join(map(str,[r.chrom,r.pos,r.ref,r.alt])) for r in ds.itertuples()]
    dm=pd.read_csv(DF/'all_individual_associated_dosages.tsv.gz',sep='\t',index_col=0)
    if set(dm.columns)!=set(keys):raise RuntimeError('Dennis cache marker mismatch')
    dm=dm[keys];meta=pd.read_csv(DF/'individual_membership.tsv',sep='\t').fillna('')
    meta['D']=pd.to_numeric(meta.D,errors='coerce');meta['included_in_BA']=meta.included_in_BA.map(boolean)
    frozen=pd.read_csv(DEN/'Dennis_BA90_intervals_frequencies.tsv',sep='\t',index_col=0)
    if frozen.shape!=(4,24) or list(frozen.columns)!=keys:raise RuntimeError('Frozen Dennis order differs')
    def dpool(ids):
        x=dm.loc[ids].to_numpy(float);n=np.isfinite(x).sum(axis=0)
        return np.divide(np.nansum(x,axis=0),n,out=np.full(24,np.nan),where=(n>0)&(n>=.8*len(ids)))
    dr={}
    for state,n in [(1,142),(0,281)]:
        mm=meta.loc[meta.included_in_BA&meta.D.eq(state)]
        if len(mm)!=n:raise RuntimeError('Dennis cohort changed')
        if not np.allclose(dpool(mm['sample']),frozen.iloc[0 if state else 3],equal_nan=True,atol=1e-8):raise RuntimeError('Dennis consensus differs')
        dr[state]=[]
        ranks=mm.groupby('race').size().sort_values(ascending=False,kind='stable')
        for race in ranks.index[:4]:
            ids=mm.loc[mm.race.eq(race),'sample'];dr[state].append((f'{race} | D{state} | n={len(ids)}',dpool(ids),'frequency'))
    rows=dr[1]+[('Dennis present consensus (n=142)',frozen.iloc[0].to_numpy(float),'consensus'),
        ('Left (corrected)',frozen.iloc[1].to_numpy(float),'genotype'),('Right (corrected)',frozen.iloc[2].to_numpy(float),'genotype'),
        ('Dennis absent consensus (n=281)',frozen.iloc[3].to_numpy(float),'consensus')]+dr[0]
    associated=[r.alt if boolean(r.pattern_is_alt) else r.ref for r in ds.itertuples()]
    other=[r.ref if boolean(r.pattern_is_alt) else r.alt for r in ds.itertuples()]
    den=dict(sites=ds.pos.to_numpy(),rows=rows,a=associated,b=other,polarity='Purple: Dennis-associated allele; orange: other allele')
    old=pd.read_csv(OLD,sep='\t').sort_values('pos1');pos=old.pos1.astype(int).tolist()
    ym=pd.read_csv(YF/'individual_ALT_dosages.tsv.gz',sep='\t',index_col=0)
    mm=pd.read_csv(V7/'membership_geography.tsv',sep='\t').fillna('');mm['v']=pd.to_numeric(mm.v,errors='coerce')
    if len(old)!=27 or list(map(int,ym.columns))!=pos or set(ym.index)!=set(mm['sample']):raise RuntimeError('Yb cache identity changed')
    ym=ym.loc[mm['sample']];freq=pd.read_csv(V7/'region_race_ALT_frequencies.tsv',sep='\t')
    def ypool(ids):return ym.loc[ids].mean(axis=0).to_numpy(float)
    yr={};bread={}
    for state in [1,0]:
        group=mm.loc[mm.v.eq(state)];v=ypool(group['sample'])
        label='All retained populations | v'+('+' if state else '-')
        saved=freq.loc[freq.row.eq(label)].set_index('pos1').reindex(pos).ALT_frequency.to_numpy(float)
        if not np.allclose(v,saved,equal_nan=True,atol=1e-8):raise RuntimeError('Yb consensus differs')
        bread[state]=(f'Yb(v){"+" if state else "-"} consensus (n={len(group)})',v,'consensus');yr[state]=[]
        ranks=group.groupby(['geography','race'],sort=True).size().sort_values(ascending=False,kind='stable')
        for geo,race in ranks.index[:4]:
            ids=group.loc[group.geography.eq(geo)&group.race.eq(race),'sample']
            yr[state].append((f'{race} | v{state} | n={len(ids)}',ypool(ids),'frequency'))
    focal=[]
    for side in ['left','right']:
        v=old[side+'_rosina_dosage'].to_numpy(float,copy=True);flip=old.rosina_base.eq(old.ref).to_numpy();v[flip]=1-v[flip]
        focal.append((side.title()+' (corrected)',v,'genotype'))
    yb=dict(sites=np.asarray(pos),rows=yr[1]+[bread[1]]+focal+[bread[0]]+yr[0],a=old.alt.tolist(),b=old.ref.tolist(),palette_type='ALT',
            polarity='Purple: ALT; orange: REF (colours do not encode phenotype)')
    return {'Dennis':den,'Yb(v)':yb}

def genes():
    wanted={v[0]:v[1:] for v in SPANS.values()};features={};parents={};exons=[]
    with gzip.open(GFF,'rt') as handle:
        for line in handle:
            if line.startswith('#'):continue
            f=line.rstrip().split('\t')
            if len(f)!=9 or f[0] not in wanted:continue
            lo,hi=wanted[f[0]];start,end=int(f[3]),int(f[4])
            if end<lo or start>hi:continue
            attrs=dict(z.split('=',1) for z in f[8].split(';') if '=' in z)
            attrs={k:unquote(v) for k,v in attrs.items()};ident=attrs.get('ID','')
            if f[2]=='gene':features[ident]=dict(chrom=f[0],start=start,end=end,strand=f[6],label=attrs.get('gene',attrs.get('Name',ident)),exons=[])
            if ident and 'Parent' in attrs:parents[ident]=attrs['Parent'].split(',')[0]
            if f[2]=='exon':exons.append((attrs.get('Parent','').split(',')[0],start,end))
    for parent,start,end in exons:
        seen=set()
        while parent not in features and parent in parents and parent not in seen:seen.add(parent);parent=parents[parent]
        if parent in features:features[parent]['exons'].append([start,end])
    return list(features.values())

def draw_letters(ax,data):
    import numpy as np
    from matplotlib.patches import Rectangle
    from matplotlib.colors import LinearSegmentedColormap,to_rgba
    cmap=LinearSegmentedColormap.from_list('alleles',['#D55E00','#F4ECDD','#54278F'])
    rows=data['rows'];heights=np.array([2.7 if r[2]=='genotype' else 1. for r in rows]);centres=np.cumsum(heights)-heights/2-.5
    ax.set_xlim(-.5,len(data['sites'])-.5);ax.set_ylim(heights.sum()-.5,-.5)
    ax.set_yticks(centres,[r[0] for r in rows],fontsize=7);ax.set_xticks([]);ax.tick_params(length=0,pad=6)
    for i,(label,vals,kind) in enumerate(rows):
        cy,rh=centres[i],heights[i];bold=kind=='genotype'
        for j,x in enumerate(vals):
            if not np.isfinite(x):ax.text(j,cy,'-',color='#999999',ha='center',va='center');continue
            a,b=data['a'][j],data['b'][j]
            if bold and x not in (0,.5,1):raise RuntimeError('Non-discrete focal genotype')
            if bold and x==.5:
                yb=data.get('palette_type')=='ALT'
                ax.add_patch(Rectangle((j-.47,cy-.43*rh),.94,.43*rh,facecolor='#FAE4D4' if yb else '#E4DDED',edgecolor='none'))
                ax.add_patch(Rectangle((j-.47,cy),.94,.43*rh,facecolor='#E4DDED' if yb else '#FAE4D4',edgecolor='none'))
                pairs=[(b,'#D55E00'),(a,'#54278F')] if yb else [(a,'#54278F'),(b,'#D55E00')]
                for yy,(base,col) in zip([cy-.22*rh,cy+.22*rh],pairs):ax.text(j,yy,base,ha='center',va='center',color=col,weight='bold',fontsize=13)
            else:
                base=a if x>.5 else b if x<.5 else a+'/'+b
                background=to_rgba(cmap(x),.18) if data.get('palette_type')=='ALT' else cmap(x)
                if kind!='consensus':ax.add_patch(Rectangle((j-.47,cy-.43*rh),.94,.86*rh,facecolor=('#E4DDED' if x>.5 else '#FAE4D4') if bold else background,edgecolor='white',lw=.3))
                letter_colour=('#54278F' if x>.5 else '#D55E00') if kind!='frequency' or data.get('palette_type')=='ALT' else ('white' if x>.8 else '#333333')
                ax.text(j,cy,base,ha='center',va='center',fontsize=14 if bold else 9 if kind=='consensus' else 8,weight='bold',color=letter_colour)
            if bold:ax.add_patch(Rectangle((j-.47,cy-.43*rh),.94,.86*rh,facecolor='none',edgecolor='#777777',lw=.65))
        if bold:ax.get_yticklabels()[i].set_weight('bold');ax.get_yticklabels()[i].set_fontsize(9)
    for sp in ax.spines.values():sp.set_visible(False)

def run(out):
    import numpy as np,pandas as pd,matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle,ConnectionPatch
    spec=importlib.util.spec_from_file_location('old_ho',HO/'worker.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
    plt.rcParams.update({'pdf.fonttype':42,'font.size':8,'font.family':'DejaVu Sans'})
    data=publication_data();native=genes();yf=json.loads((out/'Ybv_gene_features.json').read_text())
    spec=importlib.util.spec_from_file_location('frozen_ba',out/'ba_worker.py');ba_helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(ba_helper)
    den_ba,yb_ba,ba_cohorts=ba_helper.load_scores(out)
    # One tighter common Ho scale, determined before drawing either locus.
    window_tables={};max_ho=0.
    for locus in data:
        path=RUN/'optix_consistent_v2/reciprocal_mixture/optix_reciprocal_L_R_genotype_pairs.tsv.gz' if locus=='Dennis' else HO/'Ivory_Yb_exploratory_fixed_parameter_pairs.tsv.gz'
        table=old.normalize(pd.read_csv(path,sep='\t'),'Optix' if locus=='Dennis' else 'Ivory_Yb')
        _,lo,hi=SPANS[locus]
        for cutoff in [.95,.99]:
            w=old.windows(table,lo,hi,'bp',10000,cutoff);window_tables[(locus,cutoff)]=w
            finite=w[['left_Ho','right_Ho']].to_numpy().ravel();finite=finite[np.isfinite(finite)]
            if len(finite):max_ho=max(max_ho,float(finite.max()))
    ho_top=min(1.,max(.3,float(np.ceil(max_ho*1.1*10)/10)))
    fig=plt.figure(figsize=(14,17.5));audit=[]
    for index,(locus,d) in enumerate(data.items()):
        chrom,lo,hi=SPANS[locus];base=.515 if index==0 else .035
        left,width=.35,.62
        ba=fig.add_axes([left,base+.37,width,.075]);ho=fig.add_axes([left,base+.265,width,.09]);ann=fig.add_axes([left,base+.185,width,.065]);heat=fig.add_axes([left,base+.010,width,.16])
        fig.text(.025,base+.46,('A  Optix / Dennis' if locus=='Dennis' else 'B  Ivory / Yb(v)'),fontsize=13,weight='bold')
        scores=den_ba if locus=='Dennis' else yb_ba
        if locus=='Dennis':
            ba.scatter(scores.pos,scores.balanced_accuracy,s=6,c='#333333',alpha=.65,linewidths=0)
        else:
            for rule,marker,col in [('carrier','o','#333333'),('homozygote','x','#8A8A8A')]:
                take=scores.loc[scores.rule.eq(rule)]
                ba.scatter(take.pos,take.balanced_accuracy,s=8,c=col,marker=marker,alpha=.65,linewidths=.5,label=rule+' rule')
            lower,upper=int(scores.pos.min()),int(scores.pos.max())
            for x1,x2 in [(lo,lower),(upper,hi)]:
                if x2>x1:ba.axvspan(x1,x2,facecolor='#eeeeee',hatch='///',edgecolor='#cccccc',lw=0)
            ba.legend(fontsize=6,ncol=2,loc='upper right')
            ba.text(.01,.08,'Hatched: no expanded-cohort BA data',transform=ba.transAxes,fontsize=6)
        ba.axhline(.5,c='#aaaaaa',lw=.5,ls=':');ba.axhline(.9,c='#B15422',lw=.6,ls='--')
        ba.set_xlim(lo,hi);ba.set_ylim(0,1.02);ba.set_ylabel('Phenotype BA')
        ba.tick_params(axis='x',labelbottom=False);ba.spines[['top','right']].set_visible(False)
        for cutoff,style in [(.95,'-'),(.99,'--')]:
            w=window_tables[(locus,cutoff)]
            for side,color in [('left','#1764AB'),('right','#D04B35')]:ho.plot(w.mid,w[side+'_Ho'],color=color,ls=style,lw=1.25 if cutoff==.95 else .8,label=side.title()+f' | P >= {cutoff}')
        ho.set_xlim(lo,hi);ho.set_ylim(0,ho_top);ho.set_ylabel('Heterozygous fraction\namong callable SNPs');ho.tick_params(axis='x',labelbottom=False)
        ho.legend(ncol=2,fontsize=7,loc='upper right');ho.spines[['top','right']].set_visible(False)
        ann.set_xlim(lo,hi);ann.set_ylim(-.2,4.5);ann.axhline(0,color='#999999',lw=.5)
        for label,a,b in INTERVALS[locus]:
            ann.add_patch(Rectangle((a,.1),b-a,.24,facecolor='#555555' if label.startswith(('D','Yb')) else '#BBBBBB'))
            ann.text((a+b)/2,.5 if label!='D2' else .9,label,ha='center',fontsize=7)
        gs=[g for g in native if g['chrom']==chrom]
        if locus=='Yb(v)':
            c=yf['cortex'];gs=[g for g in gs if not (g['start']==c['start'] and g['end']==c['end'])]
            gs+=[dict(chrom=chrom,label='cortex',**{k:c[k] for k in ['start','end','strand','exons']})]
            iv=yf['ivory'];gs+=[dict(chrom=chrom,start=min(x[0] for x in iv['exons']),end=max(x[1] for x in iv['exons']),strand='+',label='ivory (continues)',exons=iv['exons'])]
        lanes=[lo-1]*3
        for g in sorted(gs,key=lambda z:z['start']):
            lane=next((i for i,end in enumerate(lanes) if g['start']>end),int(np.argmin(lanes)))
            lanes[lane]=g['end']+(hi-lo)*.04;y=1.5+lane*.95;a=max(lo,g['start']);b=min(hi,g['end'])
            if b<a:continue
            col='#267B82' if g['label'].startswith('ivory') else '#444444'
            ann.plot([a,b],[y,y],c=col,lw=.65)
            for a1,b1 in g['exons']:
                if b1<lo or a1>hi:continue
                ann.add_patch(Rectangle((max(lo,a1),y-.11),min(hi,b1)-max(lo,a1)+1,.22,facecolor=col,lw=0))
            x1,x2=(a,b) if g['strand']=='+' else (b,a)
            ann.annotate('',xy=(x2,y),xytext=(x1,y),arrowprops=dict(arrowstyle='->',color=col,lw=.6))
            ann.text((a+b)/2,y+.2,g['label'],fontsize=6.5,ha='center',color=col)
        ticks=np.linspace(lo,hi,5);ann.set_xticks(ticks,[f'{v/1e6:.3f}' for v in ticks]);ann.set_yticks([])
        ann.set_xlabel(('Chromosome 18' if locus=='Dennis' else 'Chromosome 15')+' (Rosina, Mb)',fontsize=7,labelpad=1)
        for sp in ann.spines.values():sp.set_visible(False)
        draw_letters(heat,d)
        for j,pos in enumerate(d['sites']):fig.add_artist(ConnectionPatch(xyA=(int(pos),0),coordsA=ann.transData,xyB=(j,-.5),coordsB=heat.transData,color='#BFC3C7',lw=.35,alpha=.5,clip_on=False))
        fig.text(.35,base-.004,d['polarity'],fontsize=7)
        audit.append(dict(locus=locus,chrom=chrom,span=[lo,hi],markers=len(d['sites']),genes=gs,rows=[r[0] for r in d['rows']]))
    fig.text(.025,.014,'Ho: 10-kb windows / 2.5-kb slide; same callable SNPs on both sides. Gene tracks and connector origins retain full genomic spans.\nSNP panels: frozen publication cohorts, allele orientations and corrected genotypes. Split cells are unphased heterozygotes; consensus is a per-site majority base.',fontsize=7)
    dest=out/'Optix_Ybv_heterozygosity_annotation_genotypes.pdf';fig.savefig(dest);plt.close(fig)
    save(out/'complete.json',dict(status='COMPLETE',PDF=str(dest),panels=audit,window_bp=10000,slide_bp=2500,
         spans_changed=False,genotypes_or_cohorts_changed=False,ivory_annotation_partial=True,annotation_source=str(GFF),
         BA_cohorts=ba_cohorts,Yb_BA_rules=['carrier','homozygote'],BA_above_Ho=True,
         Ho_ylim=[0,ho_top],heterozygote_style='top-bottom rectangles; unphased'))
    print('PDF\t'+str(dest));print('STATUS\tCOMPLETE')

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=ROOT/'analysis/Optix_Ybv_integrated_figure')
    p.add_argument('--submit',action='store_true');p.add_argument('--overwrite',action='store_true');p.add_argument('--run',action='store_true');a=p.parse_args()
    if a.submit:prepare(a)
    elif a.run:run(a.output)
    else:p.error('Use --submit')
if __name__=='__main__':main()
