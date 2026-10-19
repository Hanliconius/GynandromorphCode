#!/usr/bin/env python3
"""Locality-only diagnostic for plotted melpomene melpomene; unknown locality omitted."""
import argparse,hashlib,json,shlex,shutil,subprocess,sys
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
def run(a):
 import pandas as pd
 local=pd.read_csv(a.localities,sep='\t').fillna('');lookup={};excluded=[]
 for r in local.to_dict('records'):
  records=json.loads(r['local_metadata']);locations=set()
  for m in records:
   country=' '.join(str(m.get('Country','')).split()).replace('PanamÃ¡','Panama').replace('Panamá','Panama')
   lat=str(m.get('Latitude','')).strip();lon=str(m.get('Longitude','')).strip()
   if country and lat and lon:locations.add((country,lat,lon))
  if len(locations)==1:
   country,lat,lon=next(iter(locations));lookup[r['sample']]=f'{country} | {lat}, {lon}'
  else:excluded.append({'sample':r['sample'],'reason':'missing or conflicting coordinate locality','records':records})
 helper=Path(__file__).with_name('plot_Ybv_race_frequencies.py');base=helper.read_text()
 anchor=" for name,column in [('TNLeft','left_rosina_dosage'),('TNRight','right_rosina_dosage')]:"
 assert base.count(anchor)==1
 injection=''' # Verify original full-cohort frequencies before restricting this diagnostic.
 membership=[r for r in membership if r['sample'] in LOCALITIES]
 if len(membership)!=len(LOCALITIES):raise RuntimeError('Locality cohort does not match callable plotted individuals')
 for r in membership:
  r['original_taxon']=r['race'];r['race']=LOCALITIES[r['sample']]
 membership.sort(key=lambda r:(r['race'],r['sample']))
 labels=[];freq=[];callable_fraction=[];styles=[];counts=[];rows=[]
 groups=defaultdict(list)
 for r in membership:groups[(r['race'],r['d'],r['v'])].append(r['sample'])
'''
 base=base.replace(anchor,injection+anchor)
 start=base.index(' if base.is_file():');end=base.index(' for r in membership:r[',start);base=base[:start]+base[end:]
 base=base.replace("'Yb(v) | frequencies by race'","'Yb(v) | melpomene by collection locality'").replace("'Yb(v) | individual cydno genotypes'","'Yb(v) | melpomene individuals by locality'").replace("'Yb(v) | callable individuals by race'","'Yb(v) | callability by locality'")
 anchor2='  draw(np.asarray(callable_fraction),labels,styles,positions,out,pdf,True)'
 assert base.count(anchor2)==1
 base=base.replace(anchor2,anchor2+'''
  individual_labels=[r['race']+' | '+r['sample'] for r in membership]
  draw(np.asarray([values[r['sample']] for r in membership]),individual_labels,['individual']*len(membership),positions,out,pdf,individual=True)
''')
 ns={'__name__':'yb_locality_plot','__file__':str(helper),'LOCALITIES':lookup}
 exec(compile(base,str(helper),'exec'),ns)
 ns['build'](argparse.Namespace(source=a.source,baseline=ROOT/'analysis/yb_matched_cohort_v1',output=a.output))
 out=a.output;pdf=out/'Ybv_melpomene_localities_Rosina.pdf';(out/'Ybv_frequency_by_race_Rosina.pdf').rename(pdf)
 audit=json.loads((out/'audit.json').read_text());audit.update(PDF=str(pdf),full_source_cohort_verified_before_subsetting=True,pooled_frequencies_match_completed_figure=False,selection='Exact melpomene melpomene locality audit membership with unique known country and coordinates; no genotype-based exclusions',excluded_unknown_or_conflicting=excluded,locality_groups=sorted(set(lookup.values())),pages=['locality frequencies with focal rows','locality callability','individual genotypes by locality'])
 (out/'audit.json').write_text(json.dumps(audit,indent=2)+'\n');shutil.copy2(__file__,out/Path(__file__).name)
 print('LOCALITY_INDIVIDUALS\t'+str(len(lookup))+'\nEXCLUDED\t'+str(len(excluded))+'\nPDF\t'+str(pdf))
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,default=ROOT/'analysis/yb_completed_reads_recovery_v1');p.add_argument('--localities',type=Path,default=ROOT/'analysis/ybv_melpomene_localities_v1/melpomene_individual_localities.tsv');p.add_argument('--output',type=Path,default=ROOT/'analysis/ybv_melpomene_locality_plot_v1');p.add_argument('--submit',action='store_true');a=p.parse_args()
 if not a.submit:return run(a)
 helper=Path(__file__).with_name('plot_Ybv_race_frequencies.py');digest=hashlib.sha256(Path(__file__).read_bytes()+helper.read_bytes()).hexdigest()[:12];stage=ROOT/'analysis/yb_locality_plot_code'/digest;stage.mkdir(parents=True,exist_ok=True)
 for f in (Path(__file__),helper):shutil.copy2(f,stage/f.name)
 args=[sys.executable,str(stage/Path(__file__).name),'--source',str(a.source),'--localities',str(a.localities),'--output',str(a.output)]
 job=subprocess.check_output(['sbatch','--parsable','--partition=nano','--cpus-per-task=1','--mem=8G','--time=00:30:00','--job-name=Ybv_localplot','--output='+str(stage/'plot.%j.out'),'--error='+str(stage/'plot.%j.err'),'--wrap='+shlex.join(args)],text=True).strip();print('JOB\t'+job)
if __name__=='__main__':main()
