#!/usr/bin/env python3
"""Optix review package; no calling, phenotype recoding or focal-state SNP selection."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd
import pysam
import plot_optix_uniform_genoplots as base
import plot_optix_merged_race_frequencies as merged
import plot_optix_publication_nonPK as pooled

CHECK_PK={'SRR11243353':'P1007','SRR11243344':'P1008'}

def prepare(cfg,root,pk):
    meta,sources=base.phenotype_manifest(cfg,root,pk)
    known=meta.group.str.lower().isin(['cythera','vulcanus']) & ~meta.is_PK
    meta.loc[known,['B','D','R']]=[1,0,0]
    for mod in base.MODULES:meta.loc[known,mod+'_selection_eligible']=True
    return meta,sources

def pair_doses(r,alt):
    result=[]
    for side in ('left','right'):
        gt=r.get('mixture_'+side+'_GT')
        d=sum(map(int,gt.split('/')))/2 if gt in ('0/0','0/1','1/1') else np.nan
        result.append(d if alt else 1-d)
    return result

def conditional_scores(d,y,stratum):
    rows=[]
    for value in (0,1):
        take=(stratum==value)&np.isfinite(y)
        yy=y[take].astype(int)
        counts=[int(sum(yy==c)) for c in (0,1)]
        result=base.score(d[take],yy) if min(counts)>=2 else None
        rows.append(dict(condition_value=value,n_absent=counts[0],n_present=counts[1],
            status='scored' if result else 'insufficient class size, callability, or tied orientation',**(result or {})))
    return rows

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True);p.add_argument('--pk-catalog',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();out=a.output;out.mkdir(parents=True,exist_ok=False)
    import matplotlib
    matplotlib.use('Agg');matplotlib.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42})
    from matplotlib.backends.backend_pdf import PdfPages
    original=json.loads((a.run/'run.json').read_text());post=a.run/'optix_consistent_v2'
    marker=json.loads((post/'corrected_calls.complete.json').read_text())
    if marker['status']!='AUDITED_CALLS_AND_MIXTURE_COMPLETE':raise RuntimeError('Corrected focal calls incomplete')
    datasets={};audit={}
    for name in ('merged','recalled_only'):
        folder=out/name;folder.mkdir()
        cfg=dict(original)
        path=a.run/('Optix_legacy_merge_v1/Optix_legacy_plus_recalled.vcf.gz' if name=='merged' else 'integration_Optix/Optix_window.vcf.gz')
        with pysam.VariantFile(str(path)) as vf:names=list(vf.header.samples)
        existing={r['sample']:r for r in original['samples']}
        cfg['samples']=[existing.get(n,dict(sample=n,label=n,group='Unresolved taxon')) for n in names]
        cfg['comparative']=names;cfg['merged_vcf']=str(path)
        meta,sources=prepare(cfg,base.ROOT,a.pk_catalog)
        meta.to_csv(folder/'phenotype_manifest.tsv',sep='\t',index=False)
        scores,dosages,recall,cohorts=merged.build(cfg,meta,post,folder)
        datasets[name]=(cfg,meta,scores,dosages,recall,cohorts)
        pages=[]
        for pdfname,mods in [('Dennis_main_candidate_BA90.pdf',['Dennis']),('Band_Ray_supplement_BA90.pdf',['Band','Ray'])]:
            with PdfPages(folder/pdfname) as pdf:
                for scope in (False,True):
                    for mod in mods:
                        t,c=base.display_sites(scores,recall,.9,scope,mod,80)
                        pooled.draw_frequency(t,dosages,recall,meta,names,cohorts[mod],mod,.9,scope,folder,pdf,c)
                        t.to_csv(folder/f'{mod}_BA90_{"intervals" if scope else "full"}_sites.tsv',sep='\t',index=False)
                        pages.append(dict(module=mod,interval_only=scope,**c))
        memberships={m:pooled.frequency_groups(meta,names,m) for m in base.MODULES}
        (folder/'control_membership.json').write_text(json.dumps(memberships,indent=2)+'\n')
        audit[name]=dict(vcf=str(path),samples=len(names),pages=pages,phenotype_sources=sources,
            unknown_phenotypes=meta.loc[meta[['B','D','R']].isna().any(axis=1),'sample'].tolist())
    shared_comparison=datasets['merged'][2].merge(datasets['recalled_only'][2],
        on=base.KEY+['module'],suffixes=('_merged','_recalled'))
    shared_comparison['BA_difference']=shared_comparison.balanced_accuracy_merged-shared_comparison.balanced_accuracy_recalled
    shared_comparison['orientation_agrees']=shared_comparison.pattern_is_alt_merged==shared_comparison.pattern_is_alt_recalled
    shared_comparison.to_csv(out/'same_SNPs_merged_vs_recalled_scores.tsv.gz',sep='\t',index=False)
    cfg,meta,scores,dosages,recall,cohorts=datasets['merged'];names=cfg['comparative']
    # Band audit includes every BA>=.80 SNP BEFORE any focal-state consideration.
    band=scores.loc[(scores.module=='Band')&(scores.balanced_accuracy>=.8)].copy()
    rows=[];keys=set()
    for r in band.to_dict('records'):
        key=tuple(r[k] for k in base.KEY);q=recall.get(key,{})
        l,rr=pair_doses(q,r['pattern_is_alt'])
        rows.append(dict(**r,left_associated_dosage=l,right_associated_dosage=rr,
            pair_posterior=q.get('pair_posterior'),expected_left0_righthet=bool(l==0 and rr==.5),
            left_GT=q.get('mixture_left_GT'),right_GT=q.get('mixture_right_GT')))
        keys.add(key)
    pd.DataFrame(rows).to_csv(out/'Band_all_BA80_sites_focal_audit.tsv',sep='\t',index=False,na_rep='NA')
    # Individual PK diagnostics use independent SNP scores and ordinary display filters.
    pk=meta.loc[meta.is_PK].copy();pk['label']=pk.apply(lambda r:r['label']+(' [VERIFY PHOTO: R0]' if r['sample'] in CHECK_PK else ''),axis=1)
    with PdfPages(out/'PK_individual_diagnostics_Band_Dennis_Ray.pdf') as pdf:
        for mod in base.MODULES:
            for scope in (False,True):
                t,c=base.display_sites(scores,recall,.9,scope,mod,80)
                pooled.draw_individual(t,dosages,recall,pk,names,cohorts[mod],mod,.9,scope,out,pdf,c)
                keys.update(tuple(r[k] for k in base.KEY) for r in t.to_dict('records'))
    # Save per-individual GT quality and AD, never infer biological B/b from REF/ALT.
    calls=[]
    with pysam.VariantFile(cfg['merged_vcf']) as vf:
        for rec in vf:
            key=(rec.contig,rec.pos,rec.ref,rec.alts[0])
            if key not in keys:continue
            for n in names:
                c=rec.samples[n];gt=c.get('GT') or ()
                calls.append(dict(zip(base.KEY,key),sample=n,GT='/'.join('.' if x is None else str(x) for x in gt),
                    DP=c.get('DP'),GQ=c.get('GQ'),AD=','.join(str(x) for x in (c.get('AD') or ())),
                    callable_DP2_GQ20=bool(np.isfinite(base.valid_dosage(c)))))
    pd.DataFrame(calls).to_csv(out/'Band_and_PK_diagnostic_calls.tsv.gz',sep='\t',index=False,na_rep='NA')
    # Test D interval association with B within D classes, and D within B classes.
    nonpk=meta.loc[~meta.is_PK].copy();idx={n:i for i,n in enumerate(names)}
    cross=nonpk.groupby(['B','D','R'],dropna=False).size().rename('n').reset_index()
    cross.to_csv(out/'nonPK_BDR_class_sizes.tsv',sep='\t',index=False)
    conditional=[]
    for r in scores.loc[scores.module=='Dennis'].to_dict('records'):
        if not any(lo<=r['pos']<=hi for lo,hi in base.INTERVALS['Dennis']):continue
        for target,condition,mod in [('B','D','Band'),('D','B','Dennis')]:
            t=nonpk.loc[nonpk[mod+'_selection_eligible']]
            d=dosages[int(r['site_index'])][[idx[n] for n in t['sample']]]
            for result in conditional_scores(d,t[target].to_numpy(float),t[condition].to_numpy(float)):
                conditional.append(dict(zip(base.KEY,[r[k] for k in base.KEY]),target=target,condition=condition,**result))
    pd.DataFrame(conditional).to_csv(out/'Dennis_interval_conditional_B_D_association.tsv.gz',sep='\t',index=False)
    # PK removal cannot affect any control summary in this version.
    pd.DataFrame([dict(sample=n,specimen=s,B=1,D=1,R=0,
        verification_status='Awaiting original photograph or specimen scoring sheet; genotype is not phenotype evidence',
        source='User supplied specimen table') for n,s in CHECK_PK.items()]).to_csv(out/'PK_phenotype_verification.tsv',sep='\t',index=False)
    # Geography must be explicit; do not classify races east/west from their SNPs.
    canonical=pd.read_csv(base.ROOT/'data/derived/region_panel_metadata.tsv',sep='\t',dtype=str).fillna('')
    canonical.columns=[c.strip() for c in canonical.columns]
    lookup={}
    for r in canonical.to_dict('records'):
        for token in base.table_ids(r):lookup.setdefault(token,[]).append(r)
    geography=[]
    for r in nonpk.to_dict('records'):
        hits=lookup.get(r['sample'],[])
        geography.append(dict(sample=r['sample'],taxon=r['taxon'],group=r['group'],B=r['B'],D=r['D'],R=r['R'],
            Country='; '.join(sorted({h.get('Country','') for h in hits})),
            Latitude='; '.join(sorted({h.get('Latitude','') for h in hits})),
            Longitude='; '.join(sorted({h.get('Longitude','') for h in hits})),
            Andes_side='',classification_source=''))
    pd.DataFrame(geography).to_csv(out/'Band_geographic_membership_to_verify.tsv',sep='\t',index=False)
    pd.DataFrame([dict(dataset=n,**r) for n,v in audit.items() for r in v['pages']]).to_csv(out/'merged_vs_recalled_site_counts.tsv',sep='\t',index=False)
    audit.update(phenotypes_changed=False,SNP_selection_uses_focal_state=False,
        PK_scoring_excluded=True,PK_control_summaries_excluded=True,
        styling="compact 7.2-inch vector panels; 6–10 pt sans-serif",geographic_association_status='Pending explicit verified east/west membership',
        PK_photo_verification_status='Pending primary specimen evidence',
        expected_Band_state_count=sum(r['expected_left0_righthet'] for r in rows),
        expected_Band_state_count_posterior95=sum(r['expected_left0_righthet'] and (r['pair_posterior'] or 0)>=.95 for r in rows))
    (out/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print('REVIEW_DIR\t'+str(out));print('BAND_EXPECTED_STATE_SITES\t'+str(audit['expected_Band_state_count_posterior95']))
    print('PENDING\tVerified geographic memberships; P1007/P1008 specimen photographs or original scoring')
if __name__=='__main__':main()
