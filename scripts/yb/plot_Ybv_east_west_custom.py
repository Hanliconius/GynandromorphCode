#!/usr/bin/env python3
"""Explicitly labelled Ybv cohort sensitivity views; no genotype-based exclusions."""
import argparse,hashlib,json,shlex,shutil,subprocess,sys
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
def choose(r,view):
 import re
 words=set(re.findall('[a-z]+',r['race'].lower()))
 cyd=bool(words&{'cydno','pachinus'})
 mel='melpomene' in words or r['race'].lower() in ('cythera','vulcanus')
 if view=='full':return True
 if view=='western_melpomene':return mel and bool(words&{'rosina','cythera','vulcanus'})
 if view=='cydno_pachinus':return cyd
 if view=='exploratory_exclusions':return not cyd and not bool(words&{'amaryllis','nanna'})
 raise ValueError(view)
def run(a):
 import pandas as pd
 locality=pd.read_csv(ROOT/'analysis/ybv_melpomene_localities_v1/melpomene_individual_localities.tsv',sep='\t').fillna('')
 region={};unknown=set()
 for r in locality.to_dict('records'):
  countries={' '.join(m.get('Country','').split()) for m in json.loads(r['local_metadata']) if m.get('Country','').strip()}
  if len(countries)!=1:unknown.add(r['sample']);continue
  region[r['sample']]='East (French Guiana)' if next(iter(countries))=='French Guiana' else 'West (other known localities)'
 helper=Path(__file__).with_name('plot_Ybv_race_frequencies.py');base=helper.read_text()
 anchor=" for name,column in [('TNLeft','left_rosina_dosage'),('TNRight','right_rosina_dosage')]:"
 assert base.count(anchor)==1
 injection=''' # Full-cohort source assertions above run before any sensitivity restriction.
 original_membership=list(membership)
 membership=[dict(r) for r in membership if r['sample'] not in UNKNOWN]
 for r in membership:
  r['original_race']=r['race'];r['original_v']=r['v']
  if r['sample'] in REGION:r['race']='H. melpomene melpomene | '+REGION[r['sample']]
  if VIEW=='custom_pool':
   import re
   words=set(re.findall('[a-z]+',r['original_race'].lower()))
   western_positive=r['original_v']==1 and ('melpomene' in words or r['original_race'].lower() in ('cythera','vulcanus')) and bool(words&{'rosina','cythera','vulcanus'})
   r['v']=int(western_positive or REGION.get(r['sample'])=='West (other known localities)')
   r['display_group']='Western melpomene v+ + West melpomene (mixed phenotype)' if r['v'] else 'All remaining individuals'
   r['race']=r['display_group'];r['d']=None
 if not membership:raise RuntimeError('Empty requested cohort')
 labels=[];freq=[];callable_fraction=[];styles=[];counts=[];rows=[]
 for v in (1,0):
  ss=[r['sample'] for r in membership if r['v']==v]
  if ss:add((('Western v+ + West melpomene (mixed phenotype)' if v else 'All remaining individuals') if VIEW=='custom_pool' else ('All recorded v+' if v else 'All recorded v-'))+f' (n={len(ss)})',ss,'summary',v=v)
 groups=defaultdict(list)
 for r in membership:
  if VIEW!='custom_pool':groups[(r['race'],r['d'],r['v'])].append(r['sample'])
'''
 base=base.replace(anchor,injection+anchor)
 # Remove historical summary rows, which would obscure these explicit cohorts.
 start=base.index(' if base.is_file():');end=base.index(' for r in membership:r[',start)
 base=base[:start]+base[end:]
 base=base.replace("'Yb(v) | frequencies by race'","'Yb(v) | '+TITLE")
 base=base.replace("Rosina allele frequency / focal dosage","Rosina allele frequency / focal dosage")
 a.output.mkdir(parents=True,exist_ok=False)
 titles={'race':'Races: melpomene East / West','custom_pool':'Custom pool: western v+ + West melpomene vs others'}
 for view,title in titles.items():
  ns={'__name__':'yb_sensitivity','__file__':str(helper),'choose':choose,'VIEW':view,'TITLE':title,'REGION':region,'UNKNOWN':unknown}
  exec(compile(base,str(helper),'exec'),ns)
  folder=a.output/view
  ns['build'](argparse.Namespace(source=a.source,baseline=ROOT/'analysis/yb_matched_cohort_v1',output=folder))
  audit=json.loads((folder/'audit.json').read_text());audit.update(view=view,cohort_title=title,genotype_based_individual_exclusions=False,selection_rule={'race':'Full eligible cohort, unknown melpomene locality omitted; melpomene split FG East versus other known localities West. User-defined geographic bins.','custom_pool':'Custom mixed-phenotype pool: western melpomene v-positive plus West melpomene, versus all remaining eligible individuals including French Guiana. Original phenotypes retained in original_v; v is display group only.'}[view])
  audit['pooled_frequencies_match_completed_figure']=False; audit['original_phenotypes_preserved_in_membership_column']='original_v'; audit['unknown_locality_excluded']=sorted(unknown);audit['full_source_cohort_verified_before_subsetting']=True
  (folder/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
  (folder/'Ybv_frequency_by_race_Rosina.pdf').rename(a.output/f'Ybv_{view}_Rosina.pdf')
  audit['PDF']=str(a.output/f'Ybv_{view}_Rosina.pdf');(folder/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
 shutil.copy2(__file__,a.output/Path(__file__).name)
 print('STATUS\tCOMPLETE')
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,default=ROOT/'analysis/yb_completed_reads_recovery_v1');p.add_argument('--output',type=Path,required=True);p.add_argument('--submit',action='store_true');a=p.parse_args()
 if not a.submit:return run(a)
 helper=Path(__file__).with_name('plot_Ybv_race_frequencies.py');digest=hashlib.sha256(Path(__file__).read_bytes()+helper.read_bytes()).hexdigest()[:12]
 stage=ROOT/'analysis/yb_cohort_sensitivity_code'/digest;stage.mkdir(parents=True,exist_ok=True)
 for f in (Path(__file__),helper):shutil.copy2(f,stage/f.name)
 args=[sys.executable,str(stage/Path(__file__).name),'--source',str(a.source),'--output',str(a.output)]
 job=subprocess.check_output(['sbatch','--parsable','--partition=nano','--cpus-per-task=1','--mem=8G','--time=00:30:00','--job-name=Ybv_cohort','--output='+str(stage/'plot.%j.out'),'--error='+str(stage/'plot.%j.err'),'--wrap='+shlex.join(args)],text=True).strip()
 print('JOB\t'+job+'\nOUTPUT\t'+str(a.output))
if __name__=='__main__':main()
