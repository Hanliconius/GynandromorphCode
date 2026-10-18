#!/usr/bin/env python3
"""Replot completed Yb expansion using the same eligible individuals at d and v."""
import argparse,hashlib,importlib.util,json,shutil
from pathlib import Path
import pandas as pd
DEFAULT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/yb_expanded_existing_controls_v2')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def matched_ids(members):
    included=members[members['included'].isin([True,'True','true',1])]
    result={}
    for ref,group in included.groupby('reference'):
        good=group[pd.to_numeric(group['d'],errors='coerce').isin([0,1]) & pd.to_numeric(group['v'],errors='coerce').isin([0,1])]
        result[ref]=set(good.loc[good.interval.eq('d'),'sample']) & set(good.loc[good.interval.eq('v'),'sample'])
    return result

def build(parent,out):
    complete=json.loads((parent/'complete.json').read_text());cfg=json.loads((parent/'run.json').read_text())
    if complete.get('status')!='COMPLETE' or not complete.get('focal_calls_unchanged'):raise RuntimeError('Input expansion is not complete with preserved focal calls')
    script=parent/'expand_Yb_existing_controls.py'
    if sha(script)!=cfg['script_sha']:raise RuntimeError('Input renderer hash changed')
    members=pd.read_csv(parent/'control_membership.tsv',sep='\t');selected=pd.read_csv(parent/'selected_sources.tsv',sep='\t');ids=matched_ids(members)
    keep=selected.apply(lambda r:r['sample'] in ids.get(r['reference'],set()),axis=1)
    selected=selected.loc[keep].copy()
    for ref,group in selected.groupby('reference'):
        if not ids[ref]:raise RuntimeError('No common cohort for '+ref)
        for interval in ('d','v'):
            g=group[group.interval.eq(interval)]
            if g['sample'].duplicated().any() or set(g['sample'])!=ids[ref]:raise RuntimeError('Unequal/duplicated matched cohort')
            if set(g[interval])!={0,1}:raise RuntimeError('Empty phenotype class for '+ref+' '+interval)
    if set(ids)!={'HmRos','Hmel2.5'}:raise RuntimeError('Missing reference cohort')
    out.mkdir(parents=True,exist_ok=False);shutil.copytree(parent/'frozen',out/'frozen')
    # Completed calls are read in place. No copies, recalls, or changes to parent data.
    (out/'calls').symlink_to((parent/'calls').resolve(),target_is_directory=True)
    shutil.copy2(parent/'legacy_reference_check.json',out/'legacy_reference_check.json')
    selected.to_csv(out/'selected_sources.tsv',sep='\t',index=False)
    cfg['output']=str(out);cfg['matched_cohort_parent']=str(parent)
    (out/'run.json').write_text(json.dumps(cfg,indent=2)+'\n')
    spec=importlib.util.spec_from_file_location('yb_completed_renderer',script);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    original_draw=module.draw
    def draw(table,matrix,labels,ref,interval,pdf,title):
        return original_draw(table,matrix,labels,ref,interval,pdf,'matched individuals; '+title)
    module.draw=draw;module.report(cfg)
    result=json.loads((out/'complete.json').read_text());(out/'complete.json').unlink()
    after=pd.read_csv(out/'control_membership.tsv',sep='\t');actual=matched_ids(after)
    for ref in ids:
        if actual[ref]!=ids[ref]:raise RuntimeError('Matched membership changed during rendering')
        for interval in ('d','v'):
            subset=after[(after.reference==ref)&(after.interval==interval)]
            if set(subset['sample'])!=ids[ref] or not subset['included'].all():raise RuntimeError('A matched individual was lost')
    result['status']='COMPLETE';result['matched_cohort']=True;result['parent_expanded_run']=str(parent);result['same_individuals_in_both_panels']=True
    result['cohort_members']={ref:sorted(ss) for ref,ss in ids.items()};result['cohort_N']={ref:len(ss) for ref,ss in ids.items()}
    phenocounts=[]
    for ref,group in selected.groupby('reference'):
        individual=group.drop_duplicates('sample')
        for (d,v),g in individual.groupby(['d','v']):phenocounts.append(dict(reference=ref,d=int(d),v=int(v),n=len(g),samples=';'.join(sorted(g['sample']))))
    pd.DataFrame(phenocounts).to_csv(out/'joint_phenotype_membership.tsv',sep='\t',index=False)
    for old,new in [('Yb_expanded_existing_controls_Rosina.pdf','Yb_matched_cohort_Rosina.pdf'),('Yb_existing_controls_Hmel25.pdf','Yb_matched_cohort_Hmel25.pdf')]:
        (out/old).rename(out/new)
    shutil.copy2(__file__,out/'plot_Yb_matched_cohort.py');result['matched_renderer_sha256']=sha(__file__)
    result['Rosina_pdf']=str(out/'Yb_matched_cohort_Rosina.pdf');result['Hmel25_pdf']=str(out/'Yb_matched_cohort_Hmel25.pdf')
    (out/'complete.json').write_text(json.dumps(result,indent=2)+'\n')
    print('MATCHED_N\t'+json.dumps(result['cohort_N']));print('FINAL_PDF\t'+result['Rosina_pdf'])

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--parent',type=Path,default=DEFAULT);p.add_argument('--output',type=Path,required=True);a=p.parse_args();build(a.parent,a.output)
if __name__=='__main__':main()
