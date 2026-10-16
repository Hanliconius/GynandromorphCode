#!/usr/bin/env python3
"""Score new Optix calls independently of PK/focal samples and draw genoplots.
All genotypes and scores come from the new callset. Legacy tables supply
phenotype/cohort annotations only; their SNP scores and GTs are never read.
"""
import argparse
import csv
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

MODULES=('Band','Dennis','Ray')
KEY=['chrom','pos','ref','alt']
TRAITS={'Band':'B','Dennis':'D','Ray':'R'}
INTERVALS={'Band':[(1111500,1117000)],'Ray':[(1156000,1161000),(1174000,1189000)],
           'Dennis':[(1215000,1217000),(1218000,1223000)]}
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')


def binary(value):
    try:
        n=float(value)
        return int(n) if n in (0,1) else None
    except (ValueError,TypeError):return None


def aliases(row):
    return {str(row.get(k,'')).strip() for k in ('sample','specimen','label','ena_run','ena_sample')
            if str(row.get(k,'')).strip()}


def table_ids(row):
    ids=set()
    for k in ('SequenceID','EarthCapeID','Accession','sample'):
        value=str(row.get(k,'' )).strip()
        if value:
            ids.add(value)
            ids.update(re.findall(r'(?:SRR|ERR|DRR|SAMN|SAMEA|ERS)\d+',value))
    return ids


def unique_value(values,label):
    good={v for v in values if v is not None}
    if len(good)>1:raise RuntimeError('Conflicting phenotype annotations: '+label)
    return next(iter(good)) if good else None


def phenotype_manifest(cfg,root,pk_path,overrides=None):
    meta_path=root/'data/derived/region_panel_metadata.tsv'
    meta=pd.read_csv(meta_path,sep='\t',dtype=str).fillna('')
    meta.columns=[c.strip() for c in meta.columns]
    records=meta.to_dict('records')
    lookup={}
    for r in records:
        for token in table_ids(r):lookup.setdefault(token,[]).append(r)
    oldpath=root/'analysis/figures/key_regions_population_context_genoplots/exact_BA_comparative_sample_manifest.tsv'
    old=pd.read_csv(oldpath,sep='\t',dtype=str).fillna('')
    legacy={}
    for r in old.to_dict('records'):
        module=r.get('BA_pattern')
        if module in MODULES:legacy.setdefault((module,r['sample']),[]).append(binary(r.get('phenotype')))
    pk=pd.read_csv(pk_path,sep='\t',dtype=str).fillna('')
    pkmap={r['sample']:r for r in pk.to_dict('records')}
    bandpath=root/'analysis/original_502/checkpoints/optix_band_original_region_reciprocal_genoplot.long.tsv.gz'
    bands={}
    if bandpath.is_file():
        bt=pd.read_csv(bandpath,sep='\t',usecols=['sample','plot_group']).drop_duplicates()
        for r in bt.to_dict('records'):
            b={'Band absent (B=0)':0,'Band present (B=1)':1}.get(r['plot_group'])
            if b is not None:bands.setdefault(r['sample'],[]).append(b)
    overrides_map={}
    if overrides:
        ot=pd.read_csv(overrides,sep='\t',dtype=str).fillna('')
        if ot['sample'].duplicated().any():raise RuntimeError('Duplicated phenotype override sample')
        overrides_map={r['sample']:r for r in ot.to_dict('records')}
        unknown=set(overrides_map)-set(cfg['comparative'])
        if unknown:raise RuntimeError('Overrides outside retained cohort: '+','.join(sorted(unknown)))
    out=[]
    for r in cfg['samples']:
        if r['sample'] not in cfg['comparative']:continue
        ids=aliases(r);hits=[]
        for token in ids:hits.extend(lookup.get(token,[]))
        ispk=r['sample'] in pkmap or str(r['group']).startswith('PK_')
        if ispk and r['sample'] not in pkmap:raise RuntimeError('Unknown PK specimen in cohort')
        entry=dict(sample=r['sample'],label=r['label'],group=r['group'],is_PK=ispk,
                   phenotype_source='canonical metadata / legacy phenotype manifest',
                   taxon='; '.join(sorted({h.get('Taxon Name','') for h in hits if h.get('Taxon Name')})))
        for module,col in TRAITS.items():
            value=unique_value([binary(h.get(col)) for h in hits],r['sample']+' '+col)
            historical=[]
            for token in ids:historical.extend(legacy.get((module,token),[]))
            if module=='Band':
                # Preserve the original western comparison eligibility; do not
                # silently broaden B=0 to every band-absent taxon.
                fallback=[]
                for token in ids:fallback.extend(bands.get(token,[]))
                eligible=any((module,token) in legacy for token in ids) and bool(fallback)
                if value is None:value=unique_value(fallback,r['sample']+' Band')
            else:
                eligible=True
                oldvalue=unique_value(historical,r['sample']+' '+module)
                value=unique_value([value,oldvalue],r['sample']+' '+module)
            if ispk:value=binary(pkmap[r['sample']][col])
            if re.search(r'cydno|pachinus',entry['taxon']+' '+str(r['group']),re.I):eligible=False
            ov=overrides_map.get(r['sample'],{})
            if ov.get(col,'')!='':
                value=binary(ov[col])
                if value is None:raise RuntimeError('Nonbinary phenotype override')
                entry['phenotype_source']='explicit phenotype override'
            if module=='Band' and ov.get('Band_eligible','')!='':
                e=binary(ov['Band_eligible'])
                if e is None:raise RuntimeError('Nonbinary Band eligibility override')
                eligible=bool(e)
            entry[col]=value
            entry[module+'_selection_eligible']=not ispk and eligible and value is not None
        if ispk:
            entry['label']=pkmap[r['sample']]['label']
            entry['PK_colour']=pkmap[r['sample']]['phenotype']
        out.append(entry)
    return pd.DataFrame(out),[str(meta_path),str(oldpath),str(bandpath),str(pk_path)]


def valid_dosage(call):
    gt=call.get('GT') or ()
    if len(gt)!=2 or any(a not in (0,1) for a in gt):return np.nan
    if (call.get('DP') or 0)<2 or call.get('GQ') is None or call.get('GQ')<20:return np.nan
    return sum(gt)/2


def score(dosages,phenotypes,min_rate=.8):
    """Dominant balanced accuracy, considering ALT and REF carrier orientation."""
    d=np.asarray(dosages,dtype=float);y=np.asarray(phenotypes,dtype=int)
    called=np.isfinite(d)
    rates=[float(np.mean(called[y==c])) for c in (0,1)]
    if min(rates)<min_rate:return None
    eligible=called
    alt_carrier=d[eligible]>0;ref_carrier=d[eligible]<1;yy=y[eligible]
    def ba(pred):return .5*(float(np.mean(pred[yy==1]))+float(np.mean(~pred[yy==0])))
    a,b=ba(alt_carrier),ba(ref_carrier)
    # Equal orientation support does not identify an associated allele.
    if abs(a-b)<1e-12:return None
    return dict(balanced_accuracy=max(a,b),pattern_is_alt=a>b,
                call_rate_absent=rates[0],call_rate_present=rates[1])


def build(cfg,meta,post,out):
    import pysam
    names=cfg['comparative'];index={n:i for i,n in enumerate(names)}
    cohorts={}
    for module in MODULES:
        t=meta.loc[meta[module+'_selection_eligible']]
        if t['is_PK'].any():raise RuntimeError('PK sample entered selection cohort')
        counts=t[TRAITS[module]].value_counts()
        if any(counts.get(c,0)<2 for c in (0,1)):
            raise RuntimeError(f'{module}: need >=2 explicitly annotated non-PK specimens per class; see phenotype_manifest.tsv')
        cohorts[module]=([index[n] for n in t['sample']],t[TRAITS[module]].astype(int).to_numpy())
    path=Path(cfg['output'])/'integration_Optix/Optix_window.vcf.gz'
    keys=[];dosages=[];scored=[]
    with pysam.VariantFile(str(path)) as vf:
        if set(vf.header.samples)!=set(names):raise RuntimeError('VCF cohort differs from run manifest')
        for rec in vf:
            if not (len(rec.ref)==1 and rec.alts and len(rec.alts)==1 and len(rec.alts[0])==1
                    and rec.ref.upper() in 'ACGT' and rec.alts[0].upper() in 'ACGT'):continue
            key=(rec.contig,rec.pos,rec.ref,rec.alts[0])
            d=np.array([valid_dosage(rec.samples[n]) for n in names])
            keys.append(key);dosages.append(d)
            for module,(indices,phen) in cohorts.items():
                result=score(d[indices],phen)
                if result:scored.append(dict(zip(KEY,key),module=module,site_index=len(keys)-1,**result))
    if len(keys)!=len(set(keys)):raise RuntimeError('Duplicated physical SNP keys')
    scores=pd.DataFrame(scored)
    if scores.empty:raise RuntimeError('No scoreable SNPs')
    scores.to_csv(out/'new_comparative_SNP_scores.tsv.gz',sep='\t',index=False)
    recalled=pd.read_csv(post/'reciprocal_mixture/optix_reciprocal_L_R_genotype_pairs.tsv.gz',sep='\t')
    if recalled.duplicated(KEY).any():raise RuntimeError('Duplicate reciprocal mixture SNP')
    recall={tuple(r[k] for k in KEY):r for r in recalled.to_dict('records')}
    return scores,np.array(dosages),recall,cohorts


def display_sites(scores,recall,cutoff,interval_only,module,max_sites):
    t=scores.loc[(scores.module==module)&(scores.balanced_accuracy>=cutoff)].copy()
    n_scored=len(t)
    if interval_only:
        inside=np.zeros(len(t),dtype=bool)
        for lo,hi in INTERVALS[module]:inside|=t.pos.between(lo,hi).to_numpy()
        t=t.loc[inside]
    n_interval=len(t)
    # Explicitly a display callability filter, never a desired-state filter.
    def callable_pair(row):
        q=recall.get(tuple(row[k] for k in KEY))
        return q is not None and q['pair_posterior']>=.95 and all(q.get('mixture_'+side+'_GT') in ('0/0','0/1','1/1') for side in ('left','right'))
    t=t.loc[t.apply(callable_pair,axis=1)] if len(t) else t
    n_callable=len(t)
    t=t.sort_values(KEY)
    if len(t)>max_sites:
        t=t.iloc[np.unique(np.rint(np.linspace(0,len(t)-1,max_sites)).astype(int))]
    return t,dict(scored_above_cutoff=n_scored,in_scope=n_interval,focal_callable=n_callable,displayed=len(t))


def individual_label(row):
    """Format pandas itertuples rows, preserving unknown phenotypes as ?."""
    vals=' '.join(f'{c}{"?" if pd.isna(getattr(row,c)) else int(getattr(row,c))}' for c in ('B','D','R'))
    group=('PK '+str(row.PK_colour)) if row.is_PK else str(row.group)
    return f'{row.label} {group} | {vals} | {row.sample}'


def draw(t,dosages,recall,meta,names,cohort,module,cutoff,interval_only,out,pdf,counts):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.patches import Rectangle,ConnectionPatch
    pk=meta.loc[meta.is_PK].sort_values(['PK_colour','label'])
    controls=meta.loc[~meta.is_PK].sort_values(['group','label'])
    individuals=pd.concat([controls,pk],ignore_index=True)
    rows=['present','TN_Left','TN_Right']+individuals['sample'].tolist()+['absent']
    labels=[f'{module}+ controls','TNLeft corrected','TNRight corrected']+[individual_label(r) for r in individuals.itertuples()]+[f'{module}− controls']
    nrows=len(rows);fig=plt.figure(figsize=(13,max(8,3+nrows*.29)))
    fig.suptitle(f'optix — {module}',x=.36,ha='left',fontsize=15)
    fig.text(.36,.93,f'BA ≥ {cutoff:.2f} | '+('labelled intervals' if interval_only else 'full window')+
             f' | {counts["focal_callable"]} callable SNPs; {len(t)} displayed',fontsize=9)
    heat=fig.add_axes([.36,.13,.59,.59]);locus=fig.add_axes([.36,.82,.59,.065])
    if len(t):
        cols=[];indices,y=cohort
        lookup={n:i for i,n in enumerate(names)}
        for r in t.to_dict('records'):
            d=dosages[int(r['site_index'])].copy()
            if not r['pattern_is_alt']:d=1-d
            focal=recall[tuple(r[k] for k in KEY)]
            f=[]
            for side in ('left','right'):
                gt=focal['mixture_'+side+'_GT'];dose=sum(int(x) for x in gt.split('/'))/2
                f.append(dose if r['pattern_is_alt'] else 1-dose)
            co=d[indices]
            cols.append([float(np.nanmean(co[y==1])),*f,*[d[lookup[n]] for n in individuals['sample']],float(np.nanmean(co[y==0]))])
        matrix=np.array(cols).T
        cmap=LinearSegmentedColormap.from_list('associated_dosage',['#D7EDF5','#ECCB63','#662A83']);cmap.set_bad('#C8C8C8')
        im=heat.pcolormesh(np.arange(len(t)+1)-.5,np.arange(nrows+1)-.5,matrix,
                          cmap=cmap,vmin=0,vmax=1,shading='flat',rasterized=False)
        heat.set_xlim(-.5,len(t)-.5)
        bar=fig.colorbar(im,cax=fig.add_axes([.965,.24,.009,.35]));bar.set_ticks([0,.5,1]);bar.ax.tick_params(labelsize=8)
        bar.set_label('Trait-associated allele dosage / control frequency',fontsize=8)
        for i,r in enumerate(t.itertuples()):
            fig.add_artist(ConnectionPatch(xyA=(r.pos/1e6,-.12),coordsA=locus.transData,
                xyB=(i,-.5),coordsB=heat.transData,color='#929aa3',alpha=.16,lw=.35,clip_on=False))
        heat.axhline(.5,color='white',lw=1);heat.axhline(2.5,color='white',lw=2)
        if len(controls):heat.axhline(2.5+len(controls),color='white',lw=2)
        heat.axhline(nrows-1.5,color='white',lw=1)
        heat.set_yticks(range(nrows),labels,fontsize=7);heat.set_ylim(nrows-.5,-.5)
        # Save exact plotted values, including every retained individual and no-call.
        pd.DataFrame(matrix,index=rows,columns=t.pos.astype(str)).to_csv(
            out/f'{module}_BA{int(cutoff*100)}_{"intervals" if interval_only else "full"}_matrix.tsv',sep='\t',na_rep='NA')
    else:
        heat.text(.5,.5,'No SNPs pass these filters',transform=heat.transAxes,ha='center');heat.set_yticks([])
    heat.set_xticks([]);heat.tick_params(axis='y',length=0)
    locus.set_xlim(1.06,1.24);locus.set_ylim(-.2,1.2);locus.axhline(0,color='#666',lw=.8)
    for mod in MODULES:
        for i,(lo,hi) in enumerate(INTERVALS[mod]):
            color='#555' if mod==module else '#c8c8c8'
            locus.add_patch(Rectangle((lo/1e6,.1),(hi-lo)/1e6,.25,color=color,lw=0))
            locus.text((lo+hi)/2e6,.63 if not (mod=='Dennis' and i==1) else 1.05,
                       {'Band':'B','Ray':'R','Dennis':'D'}[mod]+('' if mod=='Band' else str(i+1)),ha='center',fontsize=8,color=color)
    locus.set_yticks([]);locus.set_xlabel('Rosina chromosome 18 (NC_140406.1), Mb; analysis windows',fontsize=8)
    for ax in (heat,locus):
        for spine in ax.spines.values():spine.set_visible(False)
    n0=int(np.sum(cohort[1]==0));n1=int(np.sum(cohort[1]==1))
    fig.text(.36,.065,f'Selection: non-PK controls only (absent n={n0}; present n={n1}); DP≥2, GQ≥20; ≥80% callable per class.\n'
        'Focal display: reciprocal-mixture posterior ≥0.95; both sides callable. Conflicting states retained.\n'
        'Columns are consecutive SNPs; if capped, evenly spaced by genomic rank. Grey individual cells = no-call; no phasing inferred.',fontsize=8)
    pdf.savefig(fig);plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--output-dir',type=Path)
    p.add_argument('--root',type=Path,default=ROOT)
    p.add_argument('--pk-catalog',type=Path,required=True)
    p.add_argument('--phenotype-overrides',type=Path)
    p.add_argument('--max-sites',type=int,default=80)
    a=p.parse_args()
    if a.max_sites<2:p.error('--max-sites must be at least 2')
    import matplotlib
    matplotlib.use('Agg');matplotlib.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'font.size':9})
    from matplotlib.backends.backend_pdf import PdfPages
    cfg=json.loads((a.run/'run.json').read_text());post=a.run/'optix_consistent_v2'
    completion=json.loads((post/'corrected_calls.complete.json').read_text())
    if completion['status']!='AUDITED_CALLS_AND_MIXTURE_COMPLETE':raise RuntimeError('Corrected calls incomplete')
    out=a.output_dir or a.run/('Optix_genoplots_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ'))
    out.mkdir(parents=True,exist_ok=False)
    meta,sources=phenotype_manifest(cfg,a.root,a.pk_catalog,a.phenotype_overrides)
    meta.to_csv(out/'phenotype_manifest.tsv',sep='\t',index=False)
    scores,dosages,recall,cohorts=build(cfg,meta,post,out)
    pages=[];pdfpath=out/'Optix_uniform_Band_Dennis_Ray_genoplots.pdf'
    with PdfPages(pdfpath) as pdf:
        for cutoff in (.80,.90):
            for scope in (False,True):
                for module in MODULES:
                    t,counts=display_sites(scores,recall,cutoff,scope,module,a.max_sites)
                    t.to_csv(out/f'{module}_BA{int(cutoff*100)}_{"intervals" if scope else "full"}_sites.tsv',sep='\t',index=False)
                    draw(t,dosages,recall,meta,cfg['comparative'],cohorts[module],module,cutoff,scope,out,pdf,counts)
                    pages.append(dict(module=module,BA_min=cutoff,labelled_intervals=scope,**counts))
    audit=dict(input=str(a.run),comparative_genotypes='New uniformly recalled VCF only',
        phenotype_sources=sources,phenotype_overrides=str(a.phenotype_overrides) if a.phenotype_overrides else None,
        scoring='Dominant balanced accuracy: ALT-carrier versus REF-carrier; higher orientation retained; ties omitted',
        PK_and_focal_used_for_scoring=False,legacy_SNP_scores_used=False,
        focal_callability_display_filter=True,focal_expected_state_filter=False,
        max_sites_per_panel=a.max_sites,phase_inferred=False,
        unknown_phenotypes=meta.loc[meta[['B','D','R']].isna().any(axis=1),'sample'].tolist(),
        cohorts={m:meta.loc[meta[m+'_selection_eligible'],['sample',TRAITS[m]]].to_dict('records') for m in MODULES},
        pages=pages,pdf=str(pdfpath),pdf_sha256=hashlib.sha256(pdfpath.read_bytes()).hexdigest())
    (out/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print('PDF\t'+str(pdfpath));print('AUDIT\t'+str(out/'audit.json'))
    print('UNANNOTATED_SAMPLES\t'+','.join(audit['unknown_phenotypes']))
    for row in pages:print('SITES\t'+json.dumps(row))

if __name__=='__main__':main()
