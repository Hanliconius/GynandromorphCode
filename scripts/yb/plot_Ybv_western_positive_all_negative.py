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
 if view=='western_positive_all_negative':return r['v']==0 or (r['v']==1 and mel and bool(words&{'rosina','cythera','vulcanus'}))
 if view=='full':return True
 if view=='western_melpomene':return mel and bool(words&{'rosina','cythera','vulcanus'})
 if view=='cydno_pachinus':return cyd
 if view=='exploratory_exclusions':return not cyd and not bool(words&{'amaryllis','nanna'})
 raise ValueError(view)
def run(a):
 helper=Path(__file__).with_name('plot_Ybv_race_frequencies.py');base=helper.read_text()
 anchor=" for name,column in [('TNLeft','left_rosina_dosage'),('TNRight','right_rosina_dosage')]:"
 assert base.count(anchor)==1
 injection=''' # Full-cohort source assertions above run before any sensitivity restriction.
 original_membership=list(membership)
 membership=[r for r in membership if choose(r,VIEW)]
 if not membership:raise RuntimeError('Empty requested cohort')
 labels=[];freq=[];callable_fraction=[];styles=[];counts=[];rows=[]
 for v in (1,0):
  ss=[r['sample'] for r in membership if r['v']==v]
  if ss:add(f'{"Western melpomene v+" if v else "All species / regions v-"} (n={len(ss)})',ss,'summary',v=v)
 groups=defaultdict(list)
 for r in membership:groups[(r['race'],r['d'],r['v'])].append(r['sample'])
'''
 base=base.replace(anchor,injection+anchor)
 # Remove historical summary rows, which would obscure these explicit cohorts.
 start=base.index(' if base.is_file():');end=base.index(' for r in membership:r[',start)
 base=base[:start]+base[end:]
 base=base.replace("'Yb(v) | frequencies by race'","'Yb(v) | '+TITLE")
 a.output.mkdir(parents=True,exist_ok=False)
 titles={'western_positive_all_negative':'Western melpomene v+ vs all v-'}
 for view,title in titles.items():
  ns={'__name__':'yb_sensitivity','__file__':str(helper),'choose':choose,'VIEW':view,'TITLE':title}
  exec(compile(base,str(helper),'exec'),ns)
  folder=a.output/view
  ns['build'](argparse.Namespace(source=a.source,baseline=ROOT/'analysis/yb_matched_cohort_v1',output=folder))
  audit=json.loads((folder/'audit.json').read_text());audit.update(view=view,cohort_title=title,genotype_based_individual_exclusions=False,selection_rule={'western_positive_all_negative':'Positive: all eligible rosina/cythera/vulcanus melpomene. Negative: all eligible v-negative individuals across species and regions. PK excluded; no genotype-based exclusions. Geographic/species composition differs between phenotype groups.','full':'All eligible individuals','western_melpomene':'H. melpomene rosina, cythera, vulcanus; all eligible individuals; no western negative melpomene comparison available','cydno_pachinus':'All eligible cydno and pachinus; matching and nonmatching retained','exploratory_exclusions':'Exclude amaryllis, nanna, all cydno and pachinus. Post hoc illustration; not independent association evidence'}[view])
  audit['pooled_frequencies_match_completed_figure']=view=='full';audit['full_source_cohort_verified_before_subsetting']=True
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
