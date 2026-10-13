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
    p.add_argument('--vcf',type=Path,default=base.PROJECT/'external/krista_genotype_plots_original/Genotype_Plots/Optix_Hmr18-938195-1239421_biallelic.vcf')
    a=p.parse_args()
    opener=gzip.open if a.vcf.suffix=='.gz' else open
    with opener(a.vcf,'rt') as h:
        for line in h:
            if line.startswith('#CHROM'):
                header=line.rstrip().split('\t')[9:];break
        else:raise RuntimeError('VCF header not found')
    samples=a.pk_samples or [s for s in header if re.match(r'^PK',s,re.I)]
    if not samples:raise RuntimeError('No PK-prefixed samples in VCF. Supply exact sample IDs with --pk-samples.')
    if len(set(samples))!=len(samples):raise RuntimeError('Duplicate PK IDs')
    values=pd.read_csv(a.input_dir/'plot_values.tsv.gz',sep='\t')
    values=values.loc[values.module.isin(['Dennis','Ray'])]
    subsets={c:callable_sites(select(values,c,True)) for c in [.8,.9]}
    union=pd.concat(list(subsets.values())).drop_duplicates(['module']+KEY+['row'])
    sites=union[KEY].drop_duplicates()
    calls,seen=full.read_calls(a.vcf,sites,samples)
    out=a.input_dir/'DR_callable_PK';out.mkdir(exist_ok=True)
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
            draw(with_pk,cutoff,True,out,combined,pk_samples=samples,focal_callable=True)
            with_pk.to_csv(out/f'BA{int(cutoff*100)}_values_with_PK.tsv.gz',sep='\t',index=False)
            audit.append({'BA_min':cutoff,'before':len(candidate[['module']+KEY].drop_duplicates()),
                'after':len(v[['module']+KEY].drop_duplicates()),'counts':counts})
    (out/'audit.json').write_text(json.dumps({'PK_samples':samples,'VCF':str(a.vcf),'comparisons':audit,
        'focal_callability_filter':True,'focal_genotype_state_filter':False,'PK_used_for_site_filter':False,
        'PK_QC':'DP>=2 GQ>=20','PK_missing_retained':True},indent=2)+'\n')
    print('PK_SAMPLES\t'+','.join(samples))
    print('PDF\t'+str(out/'optix_DR_BA80_BA90_callable_with_and_without_PK.pdf'))
if __name__=='__main__':main()
