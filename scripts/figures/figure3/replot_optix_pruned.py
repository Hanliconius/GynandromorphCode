#!/usr/bin/env python3
"""Compare independent BA cutoffs and annotated interval restrictions from cached values."""
import argparse,json
from pathlib import Path
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import plot_optix_simple_evidence as plot

INTERVALS={'Band':[(1111500,1117000)],'Dennis':[(1215000,1217000),(1218000,1223000)],'Ray':[(1156000,1161000),(1174000,1189000)]}

def select(values,cutoff,interval_only):
    keep=values.balanced_accuracy.ge(cutoff)
    if interval_only:
        inside=pd.Series(False,index=values.index)
        for module,regions in INTERVALS.items():
            for start,end in regions:
                inside |= values.module.eq(module)&values.pos.between(start,end)
        keep &= inside
    return values.loc[keep].copy()

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--input-dir',type=Path,required=True)
    a=p.parse_args()
    values=pd.read_csv(a.input_dir/'plot_values.tsv.gz',sep='\t')
    out=a.input_dir/'pruned_comparison';out.mkdir(exist_ok=True)
    details=[]
    with PdfPages(out/'optix_BA_thresholds_full_region_vs_labelled_intervals.pdf') as combined:
        for cutoff in [.80,.90,.95]:
            for interval_only in [False,True]:
                scope='labelled_intervals' if interval_only else 'full_region'
                v=select(values,cutoff,interval_only)
                target=out/f'BA{int(cutoff*100):02d}_{scope}';target.mkdir(exist_ok=True)
                counts={m:len(v.loc[v.module.eq(m),plot.base.KEY].drop_duplicates()) for m in plot.full.MODULES}
                if v.empty:
                    details.append({'BA_min':cutoff,'scope':scope,'counts':counts,'empty':True})
                    continue
                v.to_csv(target/'plot_values.tsv.gz',sep='\t',index=False)
                for compact in [False,True]:
                    plot.draw(v,target,compact=compact,combined=combined,
                        subtitle=f'BA >= {cutoff:.2f} | '+('Trait-specific labelled intervals' if interval_only else 'Full region')+f" | B {counts['Band']}, D {counts['Dennis']}, R {counts['Ray']} SNPs")
                details.append({'BA_min':cutoff,'scope':scope,'counts':counts,'empty':False})
    audit={'input':str(a.input_dir/'plot_values.tsv.gz'),'comparative_selection_only':True,
        'focal_state_selection':False,'cluster_summary':False,'intervals':INTERVALS,
        'interval_definition':'Existing analysis windows, not inferred boundaries','pages':details}
    (out/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print('PDF\t'+str(out/'optix_BA_thresholds_full_region_vs_labelled_intervals.pdf'))
    for d in details:print(f"BA>={d['BA_min']:.2f}\t{d['scope']}\t"+json.dumps(d['counts']))
if __name__=='__main__':main()
