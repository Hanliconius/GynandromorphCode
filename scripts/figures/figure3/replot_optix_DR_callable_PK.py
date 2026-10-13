#!/usr/bin/env python3
"""Pages 2/4 refinement: drop missing focal calls, optionally restore PK individuals."""
import argparse,gzip,json,re
from pathlib import Path
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import Fig3_optix_full_range_frequency as full
import Fig3_optix_BA099_both_callable as base
from replot_optix_pruned import select
from replot_optix_DR_genoplot import draw,KEY

def callable_sites(v):
    focal=v.loc[v.row.isin(['TN_Left','TN_Right'])]
    key=['module']+KEY
    matrix=focal.pivot(index=key,columns='row',values='frequency').reindex(columns=['TN_Left','TN_Right'])
    good=matrix.loc[matrix.notna().all(axis=1)].reset_index()[key]
    return v.merge(good,on=key,validate='many_to_one')

def main():
    p=argparse.ArgumentParser();p.add_argument('--input-dir',type=Path,required=True)
    p.add_argument('--pk-samples',nargs='+')
    p.add_argument('--metadata',type=Path,default=base.PROJECT/'data/derived/region_panel_metadata.tsv')
    p.add_argument('--vcf',type=Path,default=base.PROJECT/'external/krista_genotype_plots_original/Genotype_Plots/Optix_Hmr18-938195-1239421_biallelic.vcf')
    a=p.parse_args()
    opener=gzip.open if a.vcf.suffix=='.gz' else open
    with opener(a.vcf,'rt') as h:
        for line in h:
            if line.startswith('#CHROM'):
                header=line.rstrip().split('\t')[9:];break
        else:raise RuntimeError('VCF header not found')
    if a.pk_samples:
        samples=a.pk_samples
        matched_metadata=[]
    else:
        meta=pd.read_csv(a.metadata,sep='\t',dtype=str).fillna('')
        # Recover the same biological group as analysis 151, not a sample-ID prefix.
        text=meta.astype(str).agg(' | '.join,axis=1)
        hit=text.str.contains(r'Butterfly[ _-]*World[ _-]*piano[ _-]*keys|piano[ _-]*keys|BWPK',case=False,regex=True)
        matched=meta.loc[hit].copy()
        if matched.empty:
            raise RuntimeError('No piano-key group found in canonical metadata. Metadata columns: '+', '.join(meta.columns))
        samples=[];matched_metadata=[]
        valid=set(header)
        for _,row in matched.iterrows():
            candidates=sorted(set(str(x).strip() for x in row.values)&valid)
            if len(candidates)!=1:
                raise RuntimeError('BWPK metadata row must resolve to exactly one VCF sample; matches='+repr(candidates)+'; row='+repr(row.to_dict()))
            samples.append(candidates[0]);matched_metadata.append(row.to_dict())
        samples=sorted(set(samples))
        if len(samples)!=5:
            raise RuntimeError(f'Expected the five analysis-151 BWPK individuals; resolved {len(samples)}: {samples}')
    if len(set(samples))!=len(samples):raise RuntimeError('Duplicate PK IDs')
    # Phenotypes are independently recorded metadata, never inferred from genotype calls.
    phenotype_meta=pd.read_csv(a.metadata,sep='\t',dtype=str).fillna('')
    phenotype_fields=['B','D','R','Yb-dorsal','Yb-ventral','Br','N','Sb']
    labels={};phenotypes={}
    def state(value):
        try:
            number=float(value)
            return str(int(number)) if number in (0,1) else str(value)
        except (ValueError,TypeError):return '?'
    for sample in samples:
        candidates=phenotype_meta.loc[phenotype_meta.apply(lambda col:col.astype(str).str.strip().eq(sample)).any(axis=1)]
        if len(candidates)!=1:raise RuntimeError(f'{sample}: phenotype metadata must match exactly one row')
        row=candidates.iloc[0]
        phenotypes[sample]={field:state(row.get(field,'')) for field in phenotype_fields}
        ph=phenotypes[sample]
        def sign(field):return '+' if ph[field]=='1' else ('-' if ph[field]=='0' else '?')
        labels[sample]=f"{sample} | B{sign('B')} D{sign('D')} R{sign('R')} | Yb d{sign('Yb-dorsal')}/v{sign('Yb-ventral')}"
    samples=sorted(samples,key=lambda sample:(phenotypes[sample]['D']!='1',sample))
    values=pd.read_csv(a.input_dir/'plot_values.tsv.gz',sep='\t')
    values=values.loc[values.module.isin(['Dennis','Ray'])]
    subsets={c:callable_sites(select(values,c,True)) for c in [.8,.9]}
    union=pd.concat(list(subsets.values())).drop_duplicates(['module']+KEY+['row'])
    sites=union[KEY].drop_duplicates()
    calls,seen=full.read_calls(a.vcf,sites,samples)
    out=a.input_dir/'DR_callable_PK';out.mkdir(exist_ok=True)
    pd.DataFrame([{'Accession':sample,**phenotypes[sample]} for sample in samples]).to_csv(out/'PK_recorded_phenotypes.tsv',sep='\t',index=False)
    audit=[]
    with PdfPages(out/'optix_DR_BA80_BA90_callable_with_and_without_PK.pdf') as combined:
        for cutoff,v in subsets.items():
            candidate=select(values,cutoff,True)
            counts={m:len(v.loc[v.module.eq(m),KEY].drop_duplicates()) for m in ['Ray','Dennis']}
            draw(v,cutoff,True,out,combined,focal_callable=True)
            pk=[]
            for _,site in v.loc[v.row.eq('TN_Left')].iterrows():
                is_alt=base.pattern_is_alt(site.selected_orientation)
                key=tuple(site[k] for k in KEY)
                for sample in samples:
                    dose=calls.get((key,sample),np.nan)
                    pk.append({**site.to_dict(),'row':sample,'frequency':dose if is_alt else 1-dose,
                        'source':'individual VCF GT; DP>=2 GQ>=20','pair_posterior':np.nan})
            with_pk=pd.concat([v,pd.DataFrame(pk)],ignore_index=True)
            draw(with_pk,cutoff,True,out,combined,pk_samples=samples,focal_callable=True,row_labels=labels)
            with_pk.to_csv(out/f'BA{int(cutoff*100)}_values_with_PK.tsv.gz',sep='\t',index=False)
            audit.append({'BA_min':cutoff,'before':len(candidate[['module']+KEY].drop_duplicates()),
                'after':len(v[['module']+KEY].drop_duplicates()),'counts':counts})
    (out/'audit.json').write_text(json.dumps({'PK_samples':samples,'PK_phenotypes':phenotypes,'PK_row_labels':labels,'PK_metadata':str(a.metadata),'matched_metadata':matched_metadata,'VCF':str(a.vcf),'comparisons':audit,
        'focal_callability_filter':True,'focal_genotype_state_filter':False,'PK_used_for_site_filter':False,
        'PK_QC':'DP>=2 GQ>=20','PK_missing_retained':True},indent=2)+'\n')
    print('PK_SAMPLES\t'+','.join(samples))
    print('PDF\t'+str(out/'optix_DR_BA80_BA90_callable_with_and_without_PK.pdf'))
if __name__=='__main__':main()
