#!/usr/bin/env python3
"""Uniform Rosina Yb SNPs: independent d/v scoring and corrected L/R display.

2026-10-16 adaptation, not a replacement for archived historical results.
No comparative SNP recall. Focal counts use the same bcftools mpileup settings
as the corrected Optix workflow. The completed NEW Optix fit is transferred
without refitting parameters on selected Yb sites.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logsumexp

KEY=['chrom','pos','ref','alt']
CHROM='NC_140403.1'
WINDOW=(1600000,1700000)
INTERVALS={'Ybd':(1633041,1635646),'Ybv':(1658079,1659072)}
PRIMER=(1633657,1635225)
MODEL='134_fit_reciprocal_L_R_mixture_genotypes.py'
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
PARAMS=['estimated_c_left_from_right','estimated_c_right_from_left',
        'estimated_error_rate','ALT_mapping_bias_log_odds','beta_binomial_concentration']


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def dump(path,obj):
    def scalar(x):
        if isinstance(x,np.generic):return x.item()
        if isinstance(x,Path):return str(x)
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(obj,indent=2,default=scalar)+'\n')


def load_module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def binary(x):
    try:return int(float(x)) if float(x) in (0,1) else None
    except (ValueError,TypeError):return None


def phenotype_manifest(cfg,metadata,out):
    """Explicit traits only; missing annotations never become absence."""
    import re
    raw=pd.read_csv(metadata,sep='\t',dtype=str).fillna('')
    raw.columns=[c.strip() for c in raw.columns]
    lookup={}
    for row in raw.to_dict('records'):
        ids=set()
        for k in ('SequenceID','EarthCapeID','Accession','sample'):
            text=row.get(k,'').strip()
            if text:ids.add(text);ids.update(re.findall(r'(?:SRR|ERR|DRR|ERS|SAMN|SAMEA)\d+',text))
        for token in ids:lookup.setdefault(token,[]).append(row)
    rows=[]
    for r in cfg['samples']:
        if r['sample'] not in cfg['comparative']:continue
        tokens={str(r.get(k,'')) for k in ('sample','specimen','label','ena_run','ena_sample')}
        hits=[]
        for t in tokens:hits.extend(lookup.get(t,[]))
        group=str(r['group']);is_pk=group.startswith('PK_')
        entry=dict(sample=r['sample'],label=r['label'],group=group,is_PK=is_pk)
        for mod,col in [('Ybd','Yb-dorsal'),('Ybv','Yb-ventral')]:
            vals={binary(h.get(col)) for h in hits}-{None}
            # User corrections take precedence, with their source recorded.
            if group.lower()=='cythera':
                value=0;source='Joe Hanly 2026-10-16: no Yb on cythera (dorsal/ventral)'
            elif group.lower()=='vulcanus':
                value=int(mod=='Ybv');source='Joe Hanly 2026-10-16: dorsal absent, ventral present'
            else:
                if len(vals)>1:raise RuntimeError('Conflicting '+col+' annotations: '+r['sample'])
                value=next(iter(vals)) if vals else None
                source='region_panel_metadata.tsv:'+col if vals else 'unresolved; excluded from trait scoring'
            entry[mod]=value;entry[mod+'_source']=source
        rows.append(entry)
    result=pd.DataFrame(rows)
    result.to_csv(out/'phenotype_manifest.tsv',sep='\t',index=False,na_rep='NA')
    return result


def valid_dosage(call):
    gt=call.get('GT') or ()
    if len(gt)!=2 or any(x not in (0,1) for x in gt):return np.nan
    if (call.get('DP') or 0)<2 or call.get('GQ') is None or call['GQ']<20:return np.nan
    return sum(gt)/2


def score(d,y):
    called=np.isfinite(d)
    if any(np.mean(called[y==c])<.8 for c in (0,1)):return None
    yy=y[called];x=d[called]
    def ba(p):return .5*(np.mean(p[yy==1])+np.mean(~p[yy==0]))
    a,b=ba(x>0),ba(x<1)
    if abs(a-b)<1e-12:return None
    return dict(balanced_accuracy=float(max(a,b)),pattern_is_alt=bool(a>b))


def posterior_calls(model,counts,parameters):
    """Nine ordered genotype pairs, equal priors, depth on each pooled side."""
    lr=counts['L_ref'].to_numpy(float);la=counts['L_alt'].to_numpy(float)
    rr=counts['R_ref'].to_numpy(float);ra=counts['R_alt'].to_numpy(float)
    ll=model.state_log_likelihoods(la,lr+la,ra,rr+ra,*parameters)
    post=np.exp(ll-logsumexp(ll,axis=1,keepdims=True))
    idx=np.argmax(post,axis=1);best=post[np.arange(len(idx)),idx]
    good=(best>=.95)&(lr+la>=2)&(rr+ra>=2)
    result=counts.copy();result['pair_posterior']=best
    result['left_dosage']=np.where(good,model.GL[idx],np.nan)
    result['right_dosage']=np.where(good,model.GR[idx],np.nan)
    result['left_GT']=np.where(good,model.GTL[idx],'./.')
    result['right_GT']=np.where(good,model.GTR[idx],'./.')
    return result


def pool_counts(out,cfg):
    tables=[]
    for name in cfg['focal']:
        marker=json.loads((out/(name+'.complete.json')).read_text())
        t=pd.read_csv(marker['counts'],sep='\t')
        if t.duplicated(KEY).any():raise RuntimeError('Duplicate focal SNP keys')
        row=next(r for r in cfg['samples'] if r['sample']==name)
        t['side']=row['side'];tables.append(t)
    first=tables[0][KEY].reset_index(drop=True)
    if any(not t[KEY].reset_index(drop=True).equals(first) for t in tables):
        raise RuntimeError('Focal library site keys differ')
    result=first.copy()
    for side,short in [('Left','L'),('Right','R')]:
        ts=[t for t in tables if t['side'].iloc[0]==side]
        if len(ts)!=2:raise RuntimeError('Need T2+T3 per side')
        for key,suffix in [('ref_depth','ref'),('alt_depth','alt'),('other_depth','other')]:
            result[short+'_'+suffix]=sum(t[key].to_numpy() for t in ts)
    return result


def group_mean(values):
    values=np.asarray(values,float);n=int(np.isfinite(values).sum())
    return (float(np.nanmean(values)) if len(values) and n/len(values)>=.8 else np.nan),n


def group_rows(meta,mod,names):
    index={n:i for i,n in enumerate(names)}
    groups=[]
    for state in (1,0):
        rows=meta.loc[~meta.is_PK & meta[mod].eq(state)]
        groups.append(dict(label=f'{mod} {"present" if state else "absent"} controls | n={len(rows)}',
            samples=rows['sample'].tolist(),indices=[index[n] for n in rows['sample']],kind='control'))
    for key,t in meta.loc[meta.is_PK].groupby(['group','Ybd','Ybv'],dropna=False,sort=True):
        def label(v):return '?' if pd.isna(v) else str(int(v))
        groups.append(dict(label=f'{key[0].replace("PK_","PK ")} | d{label(key[1])}/v{label(key[2])} | n={len(t)}',
            samples=t['sample'].tolist(),indices=[index[n] for n in t['sample']],kind='PK'))
    return groups


def plot_panel(sites,dosages,calls,meta,names,mod,cutoff,scope,out,pdf):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.patches import Rectangle,ConnectionPatch
    groups=group_rows(meta,mod,names)
    # One positive track, corrected L/R, separate PK classes, one negative track.
    ordered=[groups[0]]+[g for g in groups if g['kind']=='PK']+[groups[1]]
    labels=[ordered[0]['label'],'TNLeft corrected','TNRight corrected']+[g['label'] for g in ordered[1:]]
    matrix=np.full((len(labels),len(sites)),np.nan);callability=[]
    for j,r in enumerate(sites.to_dict('records')):
        d=dosages[int(r['site_index'])].copy()
        if not r['pattern_is_alt']:d=1-d
        focal=calls.loc[tuple(r[k] for k in KEY)]
        matrix[1:3,j]=[focal.left_dosage,focal.right_dosage]
        if not r['pattern_is_alt']:matrix[1:3,j]=1-matrix[1:3,j]
        for i,g in enumerate(ordered):
            row=0 if i==0 else i+2
            freq,n=group_mean(d[g['indices']]);matrix[row,j]=freq
            callability.append(dict(group=g['label'],**{k:r[k] for k in KEY},
                n_called=n,n_total=len(g['indices']),frequency=None if np.isnan(freq) else freq))
    stem=f'{mod}_BA{int(cutoff*100)}_{scope}'
    sites.to_csv(out/(stem+'_sites.tsv'),sep='\t',index=False)
    columns=[':'.join(str(r[k]) for k in KEY) for r in sites.to_dict('records')]
    pd.DataFrame(matrix,index=labels,columns=columns).to_csv(out/(stem+'_matrix.tsv'),sep='\t',na_rep='NA')
    pd.DataFrame(callability).to_csv(out/(stem+'_group_callability.tsv.gz'),sep='\t',index=False)
    fig=plt.figure(figsize=(12,max(6,3+len(labels)*.42)))
    fig.suptitle(f'{mod} — uniform Rosina SNP frequencies',x=.32,ha='left',fontsize=14)
    fig.text(.32,.92,f'BA ≥ {cutoff:.2f} | {scope} | {len(sites)} SNPs displayed',fontsize=9)
    locus=fig.add_axes([.32,.79,.61,.08]);heat=fig.add_axes([.32,.23,.61,.45])
    if len(sites):
        cmap=LinearSegmentedColormap.from_list('associated_dosage',['#D7EDF5','#ECCB63','#662A83'])
        cmap.set_bad('#C8C8C8')
        im=heat.pcolormesh(np.arange(len(sites)+1)-.5,np.arange(len(labels)+1)-.5,matrix,
            cmap=cmap,vmin=0,vmax=1,shading='flat')
        heat.set_xlim(-.5,len(sites)-.5)
        heat.set_ylim(len(labels)-.5,-.5);heat.set_yticks(range(len(labels)),labels,fontsize=8)
        heat.axhline(.5,color='white',lw=2);heat.axhline(2.5,color='white',lw=2)
        heat.axhline(len(labels)-1.5,color='white',lw=2)
        cb=fig.colorbar(im,cax=fig.add_axes([.95,.28,.012,.30]));cb.set_ticks([0,.5,1])
        cb.set_label('Associated allele frequency / corrected focal dosage',fontsize=8)
        for j,r in enumerate(sites.itertuples()):
            fig.add_artist(ConnectionPatch(xyA=(r.pos/1e6,-.1),coordsA=locus.transData,
                xyB=(j,-.5),coordsB=heat.transData,color='#929aa3',lw=.35,alpha=.20,clip_on=False))
    else:heat.text(.5,.5,'No SNPs pass these filters',transform=heat.transAxes,ha='center');heat.set_yticks([])
    heat.set_xticks([]);heat.tick_params(axis='y',length=0)
    limits=WINDOW if scope=='full_window' else INTERVALS[mod]
    pad=0 if scope=='full_window' else 300
    locus.set_xlim((limits[0]-pad)/1e6,(limits[1]+pad)/1e6)
    locus.set_ylim(-.2,1.2);locus.axhline(0,color='#666',lw=.8)
    for name,(lo,hi) in INTERVALS.items():
        if hi<limits[0]-pad or lo>limits[1]+pad:continue
        locus.add_patch(Rectangle((lo/1e6,.1),(hi-lo)/1e6,.25,color='#888' if name==mod else '#ccc',lw=0))
        locus.text((lo+hi)/2e6,.7,name+' analysis interval',ha='center',fontsize=8)
    if mod=='Ybd':
        lo,hi=PRIMER;locus.plot([lo/1e6,hi/1e6],[.50,.50],color='#222',lw=2)
        locus.text((lo+hi)/2e6,1.02,'Luca primer span (anchors)',ha='center',fontsize=7)
    locus.set_yticks([]);locus.set_xlabel('Rosina chromosome 15 (NC_140403.1), Mb',fontsize=8)
    for ax in (locus,heat):
        for spine in ax.spines.values():spine.set_visible(False)
    fig.text(.32,.07,'Selection: non-PK comparative genotypes only; DP≥2, GQ≥20; ≥80% callable per phenotype.\n'
        'Controls/PK: equal-individual allele frequencies; ≥80% callable per group. Grey = missing.\n'
        'Left/Right: T2+T3 pooled by side; archived nine-state model; new Optix fit transferred; posterior≥0.95.\n'
        'No desired L/R genotype filter; columns are consecutive SNPs. Labelled intervals are analysis windows, not deletion boundaries.',fontsize=7)
    pdf.savefig(fig);plt.close(fig)
    return dict(module=mod,BA_min=cutoff,scope=scope,displayed=len(sites))


def analyse(cfg,out):
    import pysam
    import matplotlib
    matplotlib.use('Agg');matplotlib.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42})
    from matplotlib.backends.backend_pdf import PdfPages
    for name,expected in cfg['frozen_sha256'].items():
        if digest(out/name)!=expected:raise RuntimeError('Frozen file changed: '+name)
    if digest(cfg['yb_vcf'])!=cfg['yb_vcf_sha256']:raise RuntimeError('Source Yb VCF changed')
    meta=pd.read_csv(out/'phenotype_manifest.tsv',sep='\t')
    names=cfg['comparative'];index={n:i for i,n in enumerate(names)}
    cohorts={}
    for mod in INTERVALS:
        t=meta.loc[~meta.is_PK & meta[mod].notna()]
        if any((t[mod]==c).sum()<2 for c in (0,1)):
            raise RuntimeError(mod+': need ≥2 non-PK known specimens in each class; see phenotype_manifest.tsv')
        cohorts[mod]=([index[n] for n in t['sample']],t[mod].astype(int).to_numpy())
    model=load_module(out/MODEL,'archived_model')
    fit=pd.read_csv(out/'mixture.summary.tsv',sep='\t').set_index('statistic')['value']
    parameters=[float(fit[k]) for k in PARAMS]
    focal=posterior_calls(model,pool_counts(out,cfg),parameters)
    focal.to_csv(out/'Yb_uniform_corrected_LR_all_sites.tsv.gz',sep='\t',index=False,na_rep='NA')
    calls=focal.set_index(KEY,verify_integrity=True)
    dosages=[];scored=[];keys=[];qc={n:dict(covered_called_GQ_missing=0,callable_sites=0) for n in names}
    with pysam.VariantFile(cfg['yb_vcf']) as vf:
        if set(vf.header.samples)!=set(names):raise RuntimeError('Comparative VCF cohort differs')
        for rec in vf:
            if not (rec.contig==CHROM and len(rec.ref)==1 and rec.alts and len(rec.alts)==1
                    and len(rec.alts[0])==1 and rec.ref.upper() in 'ACGT' and rec.alts[0].upper() in 'ACGT'):continue
            key=(rec.contig,rec.pos,rec.ref,rec.alts[0]);keys.append(key)
            d=np.array([valid_dosage(rec.samples[n]) for n in names]);dosages.append(d)
            for i,n in enumerate(names):
                c=rec.samples[n];gt=c.get('GT') or ()
                qc[n]['covered_called_GQ_missing']+=int(len(gt)==2 and all(v is not None for v in gt)
                    and (c.get('DP') or 0)>=2 and c.get('GQ') is None)
                qc[n]['callable_sites']+=int(np.isfinite(d[i]))
            for mod,(idx,phen) in cohorts.items():
                result=score(d[idx],phen)
                if result:scored.append(dict(zip(KEY,key),module=mod,site_index=len(dosages)-1,**result))
    if len(set(keys))!=len(keys):raise RuntimeError('Duplicate comparative SNP keys')
    if set(keys)!=set(calls.index):raise RuntimeError('Focal/comparative site keys differ')
    dump(out/'comparative_qc.json',qc)
    problems=[n for n,q in qc.items() if q['covered_called_GQ_missing'] or not q['callable_sites']]
    if problems:raise RuntimeError('Comparative GQ/callability failure: '+','.join(problems))
    scores=pd.DataFrame(scored,columns=KEY+['module','site_index','balanced_accuracy','pattern_is_alt'])
    scores.to_csv(out/'Yb_all_comparative_SNP_scores.tsv.gz',sep='\t',index=False)
    dosages=np.asarray(dosages);pages=[];pdfpath=out/'Yb_d_v_uniform_Rosina_frequency_corrected_LR_and_PK.pdf'
    with PdfPages(pdfpath) as pdf:
        for cutoff in (.8,.9):
            for scope in ('full_window','interval_zoom'):
                for mod in INTERVALS:
                    t=scores.loc[(scores.module==mod)&(scores.balanced_accuracy>=cutoff)].copy()
                    above=len(t)
                    if scope=='interval_zoom':t=t.loc[t.pos.between(*INTERVALS[mod])]
                    in_scope=len(t)
                    # Display callability only: never choose the desired L/R state.
                    t=t.loc[[np.isfinite(calls.loc[tuple(r[k] for k in KEY)].left_dosage)
                             for r in t.to_dict('records')]]
                    callable_count=len(t);t=t.sort_values(KEY)
                    if len(t)>cfg['max_sites']:
                        t=t.iloc[np.unique(np.rint(np.linspace(0,len(t)-1,cfg['max_sites'])).astype(int))]
                    page=plot_panel(t,dosages,calls,meta,names,mod,cutoff,scope,out,pdf)
                    pages.append(dict(page,above_cutoff=above,in_scope=in_scope,focal_callable=callable_count))
    dump(out/'audit.json',dict(status='COMPLETE',reference='Rosina',source_run=cfg['source_run'],
        comparative_vcf=cfg['yb_vcf'],comparative=len(names),focal_libraries=cfg['focal'],
        fit_source=cfg['fit_source'],fit_policy='Transfer completed new Optix fit unchanged; no Yb refit',
        fitted_parameters=dict(zip(PARAMS,parameters)),frozen_sha256=cfg['frozen_sha256'],
        PK_and_focal_used_for_selection=False,desired_focal_state_filter=False,
        selection='Dominant balanced accuracy, ALT/REF carrier orientation; ties omitted',
        min_DP=2,min_comparative_GQ=20,min_class_and_group_call_rate=.8,posterior_cutoff=.95,
        phenotype_manifest_sha256=digest(out/'phenotype_manifest.tsv'),
        unannotated=meta.loc[meta[['Ybd','Ybv']].isna().any(axis=1),'sample'].tolist(),
        cohorts={m:meta.loc[~meta.is_PK & meta[m].notna(),['sample',m]].to_dict('records') for m in INTERVALS},
        intervals=INTERVALS,primer_span=PRIMER,interval_interpretation='analysis windows; primer span is anchored, not a deletion call',
        pages=pages,pdf=str(pdfpath),pdf_sha256=digest(pdfpath)))
    print('BUNDLE_PDF\t'+str(pdfpath));print('AUDIT\t'+str(out/'audit.json'))


def submit(a):
    import pysam
    repo=Path(__file__).resolve().parents[2]
    source=a.run.resolve(strict=True);cfg=json.loads((source/'run.json').read_text())
    completion=json.loads((source/'optix_consistent_v2/corrected_calls.complete.json').read_text())
    if completion['status']!='AUDITED_CALLS_AND_MIXTURE_COMPLETE':raise RuntimeError('New Optix fit incomplete')
    archived=repo/'scripts/reciprocal_mixing'/MODEL
    if digest(archived)!=completion['model_sha256']:
        raise RuntimeError('Archived model differs from the completed new Optix fit')
    if cfg.get('no_hmel25') is not True:raise RuntimeError('Expected Rosina-only workflow')
    if len(cfg['focal'])!=4 or not cfg['comparative']:raise RuntimeError('Bad cohort')
    from collections import Counter
    if Counter(r['side'] for r in cfg['samples'] if r['sample'] in cfg['focal'])!={'Left':2,'Right':2}:
        raise RuntimeError('Expected two libraries per side')
    yb=source/'integration_Yb/Yb_window.vcf.gz'
    with pysam.VariantFile(str(yb)) as vf:
        if set(vf.header.samples)!=set(cfg['comparative']):raise RuntimeError('Yb VCF cohort mismatch')
        if not any(r.contig==CHROM and WINDOW[0]<=r.pos<=WINDOW[1] for r in vf):
            raise RuntimeError('No variants in Yb window')
    for tool in (cfg['bcftools'],'sbatch'):
        if not shutil.which(tool):raise RuntimeError('Unavailable executable: '+tool)
    for r in cfg['samples']:
        if r['sample'] in cfg['focal']:
            with pysam.AlignmentFile(r['bam'],'rb') as bam:
                bam.check_index()
                if {x.get('SM') for x in bam.header.to_dict().get('RG',[])}!={r['sample']}:
                    raise RuntimeError('Focal BAM SM mismatch')
    fit=source/'optix_consistent_v2/reciprocal_mixture/optix_reciprocal_L_R_mixture.summary.tsv'
    parameters=pd.read_csv(fit,sep='\t').set_index('statistic')['value']
    if not all(np.isfinite(float(parameters[k])) for k in PARAMS):raise RuntimeError('Invalid fit parameters')
    c1,c2,err,bias,kappa=[float(parameters[k]) for k in PARAMS]
    if not (0<=c1<.5 and 0<=c2<.5 and 0<err<.5 and kappa>=2):
        raise RuntimeError('Fitted parameters outside model bounds')
    timestamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    out=a.output_dir or source/('Yb_uniform_frequency_'+timestamp)
    out=out.absolute();out.mkdir(parents=True,exist_ok=False);(out/'logs').mkdir()
    # The established Optix worker is reused via a read-only Yb input alias.
    alias=out/'source_view/integration_Optix';alias.mkdir(parents=True)
    (alias/'Optix_window.vcf.gz').symlink_to(yb)
    for ext in ('.csi','.tbi'):
        p=Path(str(yb)+ext)
        if p.is_file():Path(str(alias/'Optix_window.vcf.gz')+ext).symlink_to(p)
    metadata=a.root/'data/derived/region_panel_metadata.tsv'
    phenotype_manifest(cfg,metadata,out)
    phen=pd.read_csv(out/'phenotype_manifest.tsv',sep='\t')
    for mod in INTERVALS:
        for state in (0,1):
            if ((~phen.is_PK)&phen[mod].eq(state)).sum()<2:
                raise RuntimeError('Insufficient '+mod+' control phenotypes; see '+str(out/'phenotype_manifest.tsv'))
    files={Path(__file__):out/Path(__file__).name,
        repo/'scripts/figures/figure3/extend_optix_consistent_postcall.py':out/'optix_worker.py',
        repo/'scripts/reciprocal_mixing'/MODEL:out/MODEL,
        fit:out/'mixture.summary.tsv',metadata:out/'metadata_source.tsv'}
    for src,dest in files.items():shutil.copy2(src,dest)
    # Also freeze the audited comparative phenotype source and fitted provenance.
    shutil.copy2(source/'optix_consistent_v2/corrected_calls.complete.json',out/'optix_fit.complete.json')
    cfg.update(source_run=str(source),parent=str(out/'source_view'),output=str(out),
        yb_vcf=str(yb),yb_vcf_sha256=digest(yb),limits=list(WINDOW),fit_source=str(fit),max_sites=a.max_sites,
        frozen_sha256={p.name:digest(p) for p in files.values()})
    dump(out/'run.json',cfg)
    def launch(stage,array,dep=None):
        command=[sys.executable,str(out/Path(__file__).name),'--run',str(out),'--stage',stage]
        worker=out/(stage+'.sh')
        worker.write_text('#!/usr/bin/env bash\nset -euo pipefail\nexec '+shlex.join(command)+' --task "${SLURM_ARRAY_TASK_ID:-0}"\n')
        args=['sbatch','--parsable','--partition=nano','--time=00:30:00','--mem=12G','--cpus-per-task=2',
            '--job-name=YbUniform_'+stage,'--output='+str(out/'logs'/(stage+'.%A_%a.out')),
            '--error='+str(out/'logs'/(stage+'.%A_%a.err')),'--array='+array]
        if dep:args+=['--dependency=afterok:'+dep]
        args+=[str(worker)]
        return subprocess.check_output(args,text=True).strip().split(';')[0]
    counts=launch('counts','0-3');dump(out/'jobs.json',{'counts':counts})
    analysis=launch('analyse','0',counts);dump(out/'jobs.json',{'counts':counts,'analyse':analysis})
    print('COUNTS_ARRAY\t'+counts);print('PLOT_JOB\t'+analysis);print('RUN_DIR\t'+str(out))
    print('PDF\t'+str(out/'Yb_d_v_uniform_Rosina_frequency_corrected_LR_and_PK.pdf'))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True)
    p.add_argument('--submit',action='store_true')
    p.add_argument('--root',type=Path,default=ROOT)
    p.add_argument('--output-dir',type=Path)
    p.add_argument('--max-sites',type=int,default=120)
    p.add_argument('--stage',choices=['counts','analyse'])
    p.add_argument('--task',type=int,default=0)
    a=p.parse_args()
    if a.max_sites<2:p.error('--max-sites must be ≥2')
    if a.submit:submit(a)
    elif a.stage:
        cfg=json.loads((a.run/'run.json').read_text())
        if a.stage=='counts':
            worker=load_module(a.run/'optix_worker.py','uniform_counts_worker');worker.counts(cfg,a.task)
        else:analyse(cfg,a.run)
    else:p.error('--submit or --stage required')


if __name__=='__main__':main()
