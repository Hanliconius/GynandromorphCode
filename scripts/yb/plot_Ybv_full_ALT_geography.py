#!/usr/bin/env python3
"""Unselected full-window ALT frequencies by region, phenotype, and race."""
import argparse,hashlib,json,shlex,shutil,subprocess,sys
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
WEST={'rosina','cythera','vulcanus','chioneus','zelinde','weymeri','pachinus'}
EAST={'amaryllis','aglaope','malleti','ecuadoriensis','bellula','plesseni','thelxinoe','florencia','lativitta','nanna','besckei','burchelli','meriana','xenoclea','heurippa','linaresi'}
def region(taxon):
 import re
 words=set(re.findall('[a-z]+',taxon.lower()))
 if words&WEST:return 'West','taxon range'
 if words&EAST:return 'East','taxon range; Atlantic Brazil included in East'
 return 'Unresolved','ambiguous taxon range'
def run(a):
 import numpy as np,pandas as pd
 source=a.source;complete=json.loads((source/'complete.json').read_text())
 if complete['status']!='COMPLETE' or complete['scale']!='0=REF, 1=ALT':raise RuntimeError('Require complete REF/ALT rebuild')
 members=pd.read_csv(source/'membership.tsv',sep='\t').fillna('');dos=pd.read_csv(source/'individual_ALT_dosages.tsv.gz',sep='\t',index_col=0)
 dos=dos.loc[members['sample']];positions=list(map(int,dos.columns));matrix=dos.to_numpy(float)
 if len(members)!=complete['individuals'] or len(positions)!=complete['SNPs']:raise RuntimeError('Source count mismatch')
 local={}
 for r in pd.read_csv(ROOT/'analysis/ybv_melpomene_localities_v1/melpomene_individual_localities.tsv',sep='\t').fillna('').to_dict('records'):
  cc={' '.join(m.get('Country','').split()) for m in json.loads(r['local_metadata']) if m.get('Country','').strip()}
  if len(cc)==1:local[r['sample']]='East' if next(iter(cc))=='French Guiana' else 'West'
 for i,r in members.iterrows():
  geo,evidence=region(r['race'])
  if r['sample'] in local:geo=local[r['sample']];evidence='User-defined melpomene bins: FG East; other known localities West (not a literal Andes-side classification)'
  members.loc[i,'geography']=geo;members.loc[i,'geography_evidence']=evidence
  if r['sample'] in local:members.loc[i,'race']=r['race']+' | '+geo
 out=a.output;out.mkdir(parents=True,exist_ok=False)
 members.to_csv(out/'membership_geography.tsv',sep='\t',index=False)
 members.groupby(['geography','race','d','v'],dropna=False).size().reset_index(name='individuals').to_csv(out/'surface_phenotype_counts.tsv',sep='\t',index=False)
 members[(pd.to_numeric(members.d,errors='coerce')==1)&(pd.to_numeric(members.v,errors='coerce')==0)].to_csv(out/'d_positive_v_negative_candidates.tsv',sep='\t',index=False)
 labels=[];freq=[];fractions=[];styles=[];records=[]
 def add(label,ids,style):
  x=matrix[ids];n=np.isfinite(x).sum(axis=0);f=np.divide(np.nansum(x,axis=0),n,out=np.full(len(positions),np.nan),where=n>0)
  labels.append(label+f' (n={len(ids)})');freq.append(f);fractions.append(n/len(ids) if ids else np.full(len(positions),np.nan));styles.append(style)
  for j,pos in enumerate(positions):records.append(dict(row=label,pos1=pos,n_individuals=len(ids),n_called=int(n[j]),ALT_frequency=f[j]))
 for geo in ('West','East','Unresolved'):
  for v in (1,0):
   ids=members.index[members.geography.eq(geo)&members.v.eq(v)].tolist()
   if ids or geo!='Unresolved':add(f'All {geo.lower()} | v{"+" if v else "-"}',ids,'summary')
 for (geo,race,d,v),g in members.groupby(['geography','race','d','v'],sort=True,dropna=False):
  ds='?' if str(d)=='' else '+' if float(d)==1 else '-'
  add(f'{geo} | {race.replace("Heliconius ","H. ").replace("ssp. ","")} | d{ds}/v{"+" if int(v) else "-"}',g.index.tolist(),'positive' if int(v) else 'negative')
 pd.DataFrame(records).to_csv(out/'region_race_ALT_frequencies.tsv',sep='\t',index=False,na_rep='NA')
 helper=Path(__file__).with_name('plot_Ybv_race_frequencies.py');code=helper.read_text().replace('figsize=(9.4,height)','figsize=(11.4,height)').replace('left=3.9;right=.3','left=5.9;right=.3').replace('/9.4','/11.4').replace('(9.4-left-right)','(11.4-left-right)').replace('Rosina allele frequency / focal dosage','ALT allele frequency').replace('27 unchanged SNP columns','All full-window callable SNP columns').replace('Yb(v) | frequencies by race','Yb(v) | ALT frequencies by geography and race').replace('Yb(v) | callable individuals by race','Yb(v) | callability by geography and race')
 ns={'__name__':'geo_display','__file__':str(helper)};exec(compile(code,str(helper),'exec'),ns)
 import matplotlib;matplotlib.use('Agg');matplotlib.rcParams.update({'pdf.fonttype':42})
 from matplotlib.backends.backend_pdf import PdfPages
 with PdfPages(out/'Ybv_full_ALT_geography_race.pdf') as pdf:
  ns['draw'](np.asarray(freq),labels,styles,positions,out,pdf)
  ns['draw'](np.asarray(fractions),labels,styles,positions,out,pdf,True)
 shutil.copy2(__file__,out/Path(__file__).name)
 audit=dict(status='COMPLETE',source=str(source),individuals=len(members),SNPs=len(positions),allele_scale='0=REF,1=ALT',no_new_calls=True,no_marker_selection=True,no_phenotype_changes=True,unknown_geography_retained=True,geography_bins='West/East are broad population bins; melpomene user-defined FG versus other collection localities; see membership evidence',d_positive_v_negative=int(((pd.to_numeric(members.d,errors='coerce')==1)&(pd.to_numeric(members.v,errors='coerce')==0)).sum()))
 (out/'audit.json').write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps(audit,indent=2));print('PDF\t'+str(out/'Ybv_full_ALT_geography_race.pdf'))
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,default=ROOT/'analysis/ybv_REF_ALT_rebuild_v1');p.add_argument('--output',type=Path,default=ROOT/'analysis/ybv_full_ALT_geography_v1');p.add_argument('--submit',action='store_true');a=p.parse_args()
 if not a.submit:return run(a)
 helper=Path(__file__).with_name('plot_Ybv_race_frequencies.py');tag=hashlib.sha256(Path(__file__).read_bytes()+helper.read_bytes()).hexdigest()[:12];stage=ROOT/'analysis/yb_full_ALT_geography_code'/tag;stage.mkdir(parents=True,exist_ok=True)
 for f in (Path(__file__),helper):shutil.copy2(f,stage/f.name)
 cmd=[sys.executable,str(stage/Path(__file__).name),'--source',str(a.source),'--output',str(a.output)]
 job=subprocess.check_output(['sbatch','--parsable','--partition=nano','--cpus-per-task=1','--mem=8G','--time=00:30:00','--job-name=Ybv_ALTgeo','--output='+str(stage/'plot.%j.out'),'--error='+str(stage/'plot.%j.err'),'--wrap='+shlex.join(cmd)],text=True).strip();print('JOB\t'+job)
if __name__=='__main__':main()
