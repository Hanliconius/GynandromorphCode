#!/usr/bin/env python3
"""Exploratory merged Optix controls: phenotype-pooled frequency tracks including PK."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
import pysam
import plot_optix_uniform_genoplots as base
import plot_optix_all_phenotype_controls as freq

groups=freq.frequency_groups

def build(cfg,meta,post,out):
    import pysam
    names=cfg['comparative'];index={n:i for i,n in enumerate(names)}
    cohorts={}
    for module in base.MODULES:
        t=meta.loc[meta[module+'_selection_eligible']]
        if t['is_PK'].any():raise RuntimeError('PK sample entered selection cohort')
        counts=t[base.TRAITS[module]].value_counts()
        if any(counts.get(c,0)<2 for c in (0,1)):
            raise RuntimeError(f'{module}: need >=2 explicitly annotated non-PK specimens per class; see phenotype_manifest.tsv')
        cohorts[module]=([index[n] for n in t['sample']],t[base.TRAITS[module]].astype(int).to_numpy())
    path=Path(cfg['merged_vcf'])
    keys=[];dosages=[];scored=[]
    with pysam.VariantFile(str(path)) as vf:
        if set(vf.header.samples)!=set(names):raise RuntimeError('VCF cohort differs from run manifest')
        for rec in vf:
            if not (len(rec.ref)==1 and rec.alts and len(rec.alts)==1 and len(rec.alts[0])==1
                    and rec.ref.upper() in 'ACGT' and rec.alts[0].upper() in 'ACGT'):continue
            key=(rec.contig,rec.pos,rec.ref,rec.alts[0])
            d=np.array([base.valid_dosage(rec.samples[n]) for n in names])
            keys.append(key);dosages.append(d)
            for module,(indices,phen) in cohorts.items():
                result=base.score(d[indices],phen)
                if result:scored.append(dict(zip(base.KEY,key),module=module,site_index=len(keys)-1,**result))
    if len(keys)!=len(set(keys)):raise RuntimeError('Duplicated physical SNP keys')
    scores=pd.DataFrame(scored)
    if scores.empty:raise RuntimeError('No scoreable SNPs')
    scores.to_csv(out/'new_comparative_SNP_scores.tsv.gz',sep='\t',index=False)
    recalled=pd.read_csv(post/'reciprocal_mixture/optix_reciprocal_L_R_genotype_pairs.tsv.gz',sep='\t')
    if recalled.duplicated(base.KEY).any():raise RuntimeError('Duplicate reciprocal mixture SNP')
    recall={tuple(r[k] for k in base.KEY):r for r in recalled.to_dict('records')}
    return scores,np.array(dosages),recall,cohorts

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True);p.add_argument('--pk-catalog',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    a=p.parse_args();out=a.output_dir;out.mkdir(parents=True,exist_ok=False)
    merged=a.run/'Optix_legacy_merge_v1/Optix_legacy_plus_recalled.vcf.gz'
    cfg=json.loads((a.run/'run.json').read_text())
    with pysam.VariantFile(str(merged)) as vf:names=list(vf.header.samples)
    byname={r['sample']:r for r in cfg['samples']}
    cfg['samples']=[byname.get(n,dict(sample=n,label=n,group='Unresolved taxon')) for n in names]
    cfg['comparative']=names;cfg['merged_vcf']=str(merged)
    meta,sources=base.phenotype_manifest(cfg,base.ROOT,a.pk_catalog)
    known=meta['group'].str.lower().isin(['cythera','vulcanus']) & ~meta.is_PK
    meta.loc[known,['B','D','R']]=[1,0,0]
    for mod in base.MODULES:meta.loc[known,mod+'_selection_eligible']=True
    # Exact metadata taxon names keep original and recalled individuals together.
    has_taxon=meta.taxon.fillna('').str.strip().ne('') & ~meta.is_PK
    meta.loc[has_taxon,'group']=meta.loc[has_taxon,'taxon'].str.strip()
    unresolved=meta.group.str.startswith('Unresolved')
    meta.loc[unresolved,'group']='Unresolved: '+meta.loc[unresolved,'sample']
    # No inferred phenotypes: unknowns are labelled ?, and excluded from that trait's scoring.
    meta.to_csv(out/'phenotype_manifest.tsv',sep='\t',index=False)
    post=a.run/'optix_consistent_v2'
    completion=json.loads((post/'corrected_calls.complete.json').read_text())
    if completion['status']!='AUDITED_CALLS_AND_MIXTURE_COMPLETE':raise RuntimeError('Focal calls incomplete')
    scores,dosages,recall,cohorts=build(cfg,meta,post,out)
    freq.frequency_groups=groups
    import matplotlib
    matplotlib.use('Agg');matplotlib.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42,'font.size':9})
    from matplotlib.backends.backend_pdf import PdfPages
    pdfpath=out/'Optix_merged_Band_Dennis_Ray_phenotype_frequencies.pdf';pages=[]
    with PdfPages(pdfpath) as pdf:
        for cutoff in (.8,.9):
            for scope in (False,True):
                for module in base.MODULES:
                    t,counts=base.display_sites(scores,recall,cutoff,scope,module,80)
                    t.to_csv(out/f'{module}_BA{int(cutoff*100)}_{"intervals" if scope else "full"}_sites.tsv',sep='\t',index=False)
                    freq.draw_frequency(t,dosages,recall,meta,names,cohorts[module],module,cutoff,scope,out,pdf,counts)
                    pages.append(dict(module=module,cutoff=cutoff,interval_only=scope,**counts))
    membership={module:groups(meta,names,module) for module in base.MODULES}
    (out/'frequency_group_membership.json').write_text(json.dumps(membership,indent=2)+'\n')
    audit=dict(pdf=str(pdfpath),pdf_sha256=hashlib.sha256(pdfpath.read_bytes()).hexdigest(),
        merged_vcf=str(merged),calling_histories_uniform=False,phenotype_sources=sources,
        pages=pages,display='one positive and one negative control frequency per panel including matching PK; corrected focal tracks',
        group_callability_min=.8,DP_min=2,GQ_min=20,PK_and_focal_used_for_scoring=False,
        unknown_phenotypes=meta.loc[meta[['B','D','R']].isna().any(axis=1),'sample'].tolist(),
        unresolved_taxa=meta.loc[meta.group.str.startswith('Unresolved'),'sample'].tolist(),
        cohorts={m:meta.loc[meta[m+'_selection_eligible'],['sample',base.TRAITS[m]]].to_dict('records') for m in base.MODULES})
    (out/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print('PDF\t'+str(pdfpath));print('CONTROL_ROWS_PER_PANEL\t2')
    for m,(_,y) in cohorts.items():print(f'{m}\tabsent={sum(y==0)}\tpresent={sum(y==1)}')
    print('UNKNOWN_PHENOTYPES\t'+str(len(audit['unknown_phenotypes'])))
    print('UNRESOLVED_TAXA\t'+str(len(audit['unresolved_taxa'])))
if __name__=='__main__':main()
