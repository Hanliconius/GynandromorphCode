#!/usr/bin/env python3
"""Unselected full-window ALT frequencies by region, phenotype, and race."""
import argparse,hashlib,json,shlex,shutil,subprocess,sys
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
WEST={'rosina','cythera','vulcanus','chioneus','zelinde','weymeri','pachinus','cydnides'}
EAST={'amaryllis','aglaope','malleti','ecuadoriensis','bellula','plesseni','thelxinoe','florencia','lativitta','nanna','besckei','burchelli','meriana','xenoclea','heurippa','linaresi','elevatus','vicina','timareta'}
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
 dos=dos.loc[members['sample']]
 full_matrix=dos.to_numpy(float)
 full_n=np.isfinite(full_matrix).sum(axis=1)
 full_callable=full_n/full_matrix.shape[1]
 qc=members[['sample','individual_id','race','d','v']].copy()
 qc['full_window_SNPs_called']=full_n
 qc['full_window_callable_fraction']=full_callable
 qc['full_window_REF_homozygote_fraction_among_called']=np.divide((full_matrix==0).sum(axis=1),full_n,out=np.full(len(members),np.nan),where=full_n>0)
 qc['flag_callable_below_80pct']=full_callable<.8
 qc['REF_homozygotes_are_valid_calls']=True
 old=pd.read_csv(ROOT/'analysis/yb_completed_reads_recovery_v1/frozen/HmRos_Ybv_unchanged_focal_calls.tsv.gz',sep='\t').sort_values('pos1')
 fullsites=pd.read_csv(source/'SNPs.tsv',sep='\t').set_index('pos1')
 for r in old.itertuples():
  if r.pos1 not in fullsites.index or str(fullsites.loc[r.pos1,'ref'])!=r.ref or str(fullsites.loc[r.pos1,'alt'])!=r.alt:raise RuntimeError('Original marker not identically represented in joint calls: '+str(r.pos1))
 dos=dos[[str(int(x)) for x in old.pos1]];positions=list(map(int,dos.columns));matrix=dos.to_numpy(float)
 if len(members)!=complete['individuals'] or len(positions)!=len(old):raise RuntimeError('Source count mismatch')
 local={};unknown_local=set();countries_by_sample={}
 for r in pd.read_csv(ROOT/'analysis/ybv_melpomene_localities_v1/melpomene_individual_localities.tsv',sep='\t').fillna('').to_dict('records'):
  cc={' '.join(m.get('Country','').split()) for m in json.loads(r['local_metadata']) if m.get('Country','').strip()}
  if len(cc)!=1:unknown_local.add(r['sample'])
  else:countries_by_sample[r['sample']]=next(iter(cc))
  if len(cc)==1:local[r['sample']]='East' if next(iter(cc))=='French Guiana' else 'West'
 for i,r in members.iterrows():
  geo,evidence=region(r['race'])
  if r['sample'] in local:geo=local[r['sample']];evidence='User-defined melpomene bins: FG East; other known localities West (not a literal Andes-side classification)'
  members.loc[i,'geography']=geo;members.loc[i,'geography_evidence']=evidence
  if r['sample'] in local:
   country=countries_by_sample[r['sample']]
   members.loc[i,'race']=r['race']+' | '+('Guianan (FG)' if country=='French Guiana' else 'Colombian' if country=='Colombia' else country.replace('PanamÃ¡','Panama'))
  if 'nov. ECU' in r['race']:
   members.loc[i,'geography_evidence']='Provisional East: timareta eastern Ecuador range; exact accessions need locality confirmation'
 out=a.output;out.mkdir(parents=True,exist_ok=False)
 qc.to_csv(out/'full_window_sample_quality.tsv',sep='\t',index=False)
 # Primary retains all genotype states. Separate sensitivity uses only callable
 # fraction; the REF genotype itself never triggers exclusion.
 keep=np.asarray([r['sample'] not in unknown_local for r in members.to_dict('records')]) & (full_callable>=a.min_callable)
 excluded=members.loc[~keep].copy();excluded['reason']=['unknown melpomene locality' if x in unknown_local else 'full-window callability below '+str(a.min_callable) for x in excluded['sample']]
 excluded.to_csv(out/'excluded_samples.tsv',sep='\t',index=False)
 members=members.loc[keep].reset_index(drop=True);matrix=matrix[keep]
 # Export the actual sheet records for provisional Ecuadorian assignments.
 metadata=pd.read_csv(ROOT/'data/derived/region_panel_metadata.tsv',sep='\t',dtype=str).fillna('')
 ecu=members[members.race.str.contains('nov. ECU',regex=False)]
 ids={alias for text in ecu.aliases for alias in text.split(';')}
 metadata[metadata.Accession.isin(ids)].to_csv(out/'nov_ECU_collection_metadata.tsv',sep='\t',index=False)
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
   if ids or geo!='Unresolved':add(f'All {geo.lower()} population group | v{"+" if v else "-"}',ids,'summary')
 for label,col in [('TNLeft','left_rosina_dosage'),('TNRight','right_rosina_dosage')]:
  x=old[col].to_numpy(float);flip=old.rosina_base.eq(old.ref).to_numpy();x[flip]=1-x[flip]
  labels.append(label);freq.append(x);fractions.append(np.where(np.isfinite(x),1.,np.nan));styles.append('focal')
 for (geo,race,d,v),g in members.groupby(['geography','race','d','v'],sort=True,dropna=False):
  ds='?' if str(d)=='' else '+' if float(d)==1 else '-'
  add(f'{geo} | {race.replace("Heliconius ","H. ").replace("ssp. ","")} | d{ds}/v{"+" if int(v) else "-"}',g.index.tolist(),'positive' if int(v) else 'negative')
 pd.DataFrame(records).to_csv(out/'region_race_ALT_frequencies.tsv',sep='\t',index=False,na_rep='NA')
 helper=Path(__file__).with_name('plot_Ybv_race_frequencies.py');code=helper.read_text().replace('figsize=(9.4,height)','figsize=(11.4,height)').replace('left=3.9;right=.3','left=5.9;right=.3').replace('/9.4','/11.4').replace('(9.4-left-right)','(11.4-left-right)').replace('Rosina allele frequency / focal dosage','ALT allele frequency').replace('27 unchanged SNP columns','27 original marker columns; comparative calls from expanded joint rebuild').replace('Yb(v) | frequencies by race','Yb(v) | ALT frequencies by geography and race').replace('Yb(v) | callable individuals by race','Yb(v) | callability by geography and race')
 ns={'__name__':'geo_display','__file__':str(helper)};exec(compile(code,str(helper),'exec'),ns)
 import matplotlib;matplotlib.use('Agg');matplotlib.rcParams.update({'pdf.fonttype':42})
 from matplotlib.backends.backend_pdf import PdfPages
 with PdfPages(out/'Ybv_short_ALT_geography_race.pdf') as pdf:
  ns['draw'](np.asarray(freq),labels,styles,positions,out,pdf)
  ns['draw'](np.asarray(fractions),labels,styles,positions,out,pdf,True)
 shutil.copy2(__file__,out/Path(__file__).name)
 audit=dict(status='COMPLETE',source=str(source),individuals=len(members),SNPs=len(positions),allele_scale='0=REF,1=ALT',no_new_calls=True,marker_set='Original 27 markers; exact REF/ALT checked against new joint calls',focal='Original corrected focal calls converted to ALT dosage',no_phenotype_changes=True,unknown_geography_retained=True,min_full_window_callable_fraction=a.min_callable,REF_zero_not_exclusion_criterion=True,unknown_melpomene_locality_omitted=sorted(unknown_local),provisional_ECU_assignments=True,geography_bins='West/East are broad population bins; melpomene user-defined FG versus other collection localities; see membership evidence',d_positive_v_negative=int(((pd.to_numeric(members.d,errors='coerce')==1)&(pd.to_numeric(members.v,errors='coerce')==0)).sum()))
 (out/'audit.json').write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps(audit,indent=2));print('PDF\t'+str(out/'Ybv_short_ALT_geography_race.pdf'))
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,default=ROOT/'analysis/ybv_REF_ALT_rebuild_v1');p.add_argument('--output',type=Path,default=ROOT/'analysis/ybv_short_ALT_geography_v2');p.add_argument('--min-callable',type=float,default=0.0);p.add_argument('--submit',action='store_true');a=p.parse_args()
 if not a.submit:
  run(a)
  if a.min_callable==0:
   run(argparse.Namespace(source=a.source,output=a.output/'callability80_sensitivity',min_callable=.8))
  return
 helper=Path(__file__).with_name('plot_Ybv_race_frequencies.py');tag=hashlib.sha256(Path(__file__).read_bytes()+helper.read_bytes()).hexdigest()[:12];stage=ROOT/'analysis/yb_full_ALT_geography_code'/tag;stage.mkdir(parents=True,exist_ok=True)
 for f in (Path(__file__),helper):shutil.copy2(f,stage/f.name)
 cmd=[sys.executable,str(stage/Path(__file__).name),'--source',str(a.source),'--output',str(a.output),'--min-callable',str(a.min_callable)]
 job=subprocess.check_output(['sbatch','--parsable','--partition=nano','--cpus-per-task=1','--mem=8G','--time=00:30:00','--job-name=Ybv_ALTgeo','--output='+str(stage/'plot.%j.out'),'--error='+str(stage/'plot.%j.err'),'--wrap='+shlex.join(cmd)],text=True).strip();print('JOB\t'+job)
if __name__=='__main__':main()
