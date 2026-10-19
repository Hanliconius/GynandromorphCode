#!/usr/bin/env python3
"""Two geography-stratified views of the frozen, completed Yb(v) cohort."""
import argparse,hashlib,json,shlex,shutil,subprocess,sys
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
# Conservative taxon-range assignments. Ambiguous/wide-ranging taxa stay unresolved.
WEST={'rosina','cythera','vulcanus','chioneus','zelinde','weymeri','pachinus'}
EAST={'amaryllis','aglaope','malleti','ecuadoriensis','bellula','plesseni','thelxinoe','florencia','lativitta','malletti'}
def geography(taxon):
 import re
 words=set(re.findall(r'[a-z]+',taxon.lower()))
 if words&WEST:return 'Western / trans-Andean','taxon-range assignment'
 if words&EAST:return 'Eastern / cis-Andean','taxon-range assignment'
 if 'nanna' in words or 'besckei' in words:return 'Atlantic Brazil','taxon-range assignment'
 return 'Unresolved geography','no unambiguous curated range assignment'
def run(a):
 helper=Path(__file__).with_name('plot_Ybv_race_frequencies.py')
 code=helper.read_text()
 # Reuse the audited extractor, individual deduplication and pooled-frequency assertions.
 old="membership.append(dict(sample=s,individual_id=ind,race=race,d=d,v=v,aliases=';'.join(aliases),source=r['source'],path=r['path'],callable_markers=int(np.isfinite(vector).sum())))"
 new="membership.append(dict(sample=s,individual_id=ind,race=race,d=d,v=v,geography=geography(race)[0],geography_evidence=geography(race)[1],aliases=';'.join(aliases),source=r['source'],path=r['path'],callable_markers=int(np.isfinite(vector).sum())))"
 if code.count(old)!=1:raise RuntimeError('Helper changed: membership patch requires review')
 code=code.replace(old,new)
 start=code.index(' for (race,d,v),ss in sorted(groups.items()')
 end=code.index(' out.mkdir(',start)
 loop=''' geo_groups=defaultdict(list)
 for r in membership:
  key=(r['geography'],r['race'] if VIEW=='race' else '',r['v'])
  geo_groups[key].append(r['sample'])
 for (region,race,v),ss in sorted(geo_groups.items(),key=lambda x:(x[0][0],x[0][1],-x[0][2])):
  name=region+(' | '+race.replace('Heliconius ','H. ').replace('ssp. ','') if race else '')
  label=f'{name} | v{"+" if v else "-"} (n={len(ss)})'
  add(label,ss,'positive' if v else 'negative',race,None,v)
'''
 code=code[:start]+loop+code[end:]
 code=code.replace("figsize=(9.4,height)","figsize=(12.4,height)").replace('left=3.9;right=.3','left=6.5;right=.3').replace('/9.4','/12.4').replace('(9.4-left-right)','(12.4-left-right)')
 code=code.replace("'Yb(v) | frequencies by race'","'Yb(v) | geography + race' if VIEW=='race' else 'Yb(v) | geography + phenotype'")
 a.output.mkdir(parents=True,exist_ok=False)
 for view in ('race','phenotype'):
  ns={'__name__':'yb_geography_helper','__file__':str(helper),'geography':geography,'VIEW':view}
  exec(compile(code,str(helper),'exec'),ns)
  args=argparse.Namespace(source=a.source,baseline=ROOT/'analysis/yb_matched_cohort_v1',output=a.output/view)
  ns['build'](args)
  folder=args.output
  (folder/'Ybv_frequency_by_race_Rosina.pdf').rename(a.output/f'Ybv_geography_{view}_Rosina.pdf')
  audit=json.loads((folder/'audit.json').read_text());audit['geography_rule']='Conservative taxon-range assignments; ambiguous ranges unresolved; no genotype-based assignments';audit['view']=view; audit['grouping']='geography then race then ventral phenotype' if view=='race' else 'geography then ventral phenotype'
  audit['PDF']=str(a.output/f'Ybv_geography_{view}_Rosina.pdf')
  (folder/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
 shutil.copy2(__file__,a.output/Path(__file__).name)
 (a.output/'geography_assignments.json').write_text(json.dumps({'western_taxa':sorted(WEST),'eastern_taxa':sorted(EAST),'Atlantic_Brazil':['nanna','besckei'],'other_taxa':'unresolved; notably melpomene is not assigned from its name alone'},indent=2)+'\n')
 print('STATUS\tCOMPLETE')
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,default=ROOT/'analysis/yb_completed_reads_recovery_v1');p.add_argument('--output',type=Path,required=True);p.add_argument('--submit',action='store_true');a=p.parse_args()
 if not a.submit:return run(a)
 helper=Path(__file__).with_name('plot_Ybv_race_frequencies.py')
 digest=hashlib.sha256(Path(__file__).read_bytes()+helper.read_bytes()).hexdigest()[:12]
 stage=ROOT/'analysis/yb_geography_code'/digest;stage.mkdir(parents=True,exist_ok=True)
 for f in (Path(__file__),helper):shutil.copy2(f,stage/f.name)
 args=[sys.executable,str(stage/Path(__file__).name),'--source',str(a.source),'--output',str(a.output)]
 job=subprocess.check_output(['sbatch','--parsable','--partition=nano','--cpus-per-task=1','--mem=8G','--time=00:30:00','--job-name=Ybv_geo','--output='+str(stage/'plot.%j.out'),'--error='+str(stage/'plot.%j.err'),'--wrap='+shlex.join(args)],text=True).strip()
 print('JOB\t'+job)
 for view in ('race','phenotype'):print('PDF\t'+str(a.output/f'Ybv_geography_{view}_Rosina.pdf'))
if __name__=='__main__':main()
