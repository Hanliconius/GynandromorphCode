#!/usr/bin/env python3
"""Inspect cydno/pachinus identities, plot inclusion and individual frozen-marker calls."""
import argparse,csv,hashlib,importlib.util,json,shlex,shutil,subprocess,sys
from collections import defaultdict,Counter
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
SOURCE=ROOT/'analysis/yb_completed_reads_recovery_v1'
def read(p):
 with Path(p).open() as h:return list(csv.DictReader(h,delimiter='\t'))
def tab(p,rows):
 if not rows:p.write_text('status\nNO_ROWS\n');return
 cols=list(dict.fromkeys(k for r in rows for k in r))
 with p.open('w',newline='') as h:
  w=csv.DictWriter(h,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(rows)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p,name):
 s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def build(a):
 import numpy as np,pandas as pd,pysam
 source=a.source;out=a.output;out.mkdir(parents=True,exist_ok=False)
 cfg=json.loads((source/'run.json').read_text());status=json.loads((source/'complete.json').read_text())
 if status['status']!='COMPLETE':raise RuntimeError('Recovery incomplete')
 helper=load(Path(__file__).with_name('plot_Ybv_race_frequencies.py'),'race_helper')
 if sha(source/'engine.py')!=cfg['engine_sha']:raise RuntimeError('Engine changed')
 eng=load(source/'engine.py','engine');focal=source/'frozen/HmRos_Ybv_unchanged_focal_calls.tsv.gz'
 if sha(focal)!=cfg['focal_hashes'][focal.name]:raise RuntimeError('Focal marker/call hash mismatch')
 sites=pd.read_csv(focal,sep='\t').sort_values('pos1');positions=sites.pos1.astype(int).tolist()
 phen={r['sample']:r for r in read(source/'resolved_phenotypes.tsv')}
 candidates={s:r for s,r in phen.items() if 'cydno' in r['taxon'].lower() or 'pachinus' in r['taxon'].lower()}
 selected=[r for r in read(source/'selected_sources.tsv') if r['reference']=='HmRos' and r['interval']=='v'];chosen={r['sample']:r for r in selected}
 members={r['sample']:r for r in read(source/'control_membership.tsv') if r['reference']=='HmRos' and r['interval']=='v'}
 paths={}
 for parent in (Path(cfg['parent']),Path(cfg['recovery_source'])):
  p=parent/'selected_sources.tsv'
  if p.is_file():
   for r in read(p):
    if r['reference']=='HmRos' and r['interval']=='v':paths[r['sample']]=(r['path'],r['source'])
 for t in cfg['new_tasks']:
  p=source/'calls'/f'{t["call_task"]}.vcf.gz'
  if p.is_file() and (source/'calls'/f'{t["call_task"]}.complete.json').is_file():paths[t['sample']]=(str(p),'forced_BAM_call')
 for r in selected:paths[r['sample']]=(r['path'],r['source'])
 merges={}
 for r in read(source/'individual_merges.tsv'):
  if r.get('reference')=='HmRos' and r.get('interval')=='v':merges[r['representative']]=r['runs'].split(';')
 vectors={};raw=[];quality=[]
 def get(s):
  if s in vectors:return vectors[s]
  if s not in paths:vectors[s]=np.full(len(sites),np.nan);return vectors[s]
  path,kind=paths[s];vv,q=eng.extract(path,[s],sites,legacy=kind=='legacy_VCF',allow_pl=True);value=vv[s];quality.append(dict(sample=s,path=path,**q))
  if s in merges:
   aa=[]
   for alias in merges[s]:
    if alias not in paths:raise RuntimeError('Merged alias source not found')
    p,k=paths[alias];v,_=eng.extract(p,[alias],sites,legacy=k=='legacy_VCF',allow_pl=True);aa.append(v[alias])
   value=helper.combine(aa)
  vectors[s]=value;return value
 # Independently confirm that reconstructing the selected cohort gives the same
 # pooled frequencies as the finished panel, before interpreting subgroups.
 old=pd.read_csv(source/'HmRos_Ybv_plotted_values.tsv',sep='\t',index_col=0)
 active=[]
 for r in selected:
  s=r['sample'];v=get(s)
  expected=str(members.get(s,{}).get('included','')).lower() in {'true','1'}
  if bool(np.isfinite(v).any())!=expected:raise RuntimeError('Reconstruction mismatch '+s)
  if expected:active.append(r)
 for state,index in ((1,0),(0,3)):
  ss=[r['sample'] for r in active if int(float(r['v']))==state];mat=np.array([get(s) for s in ss]);n=np.isfinite(mat).sum(0);mean=np.divide(np.nansum(mat,0),n,out=np.full(len(sites),np.nan),where=n>0)
  np.testing.assert_allclose(mean,old.iloc[index].to_numpy(float),equal_nan=True,rtol=0,atol=1e-12)
 inventory=[];seenindividuals=set()
 for s,r in sorted(candidates.items(),key=lambda x:(x[1]['taxon'],x[0])):
  v=get(s);included=str(members.get(s,{}).get('included','')).lower() in {'true','1'}
  reason='included_in_completed_Ybv_plot' if included else 'explicit_PK_or_cross_exclusion' if r.get('exclusion') else 'surface_phenotype_unresolved' if r.get('v') not in {'0','1'} else 'no_usable_calls_at_frozen_markers' if s in paths else 'not_selected_and_no_native_ventral_VCF_in_completed_sources'
  inventory.append(dict(sample=s,individual_id=r['individual_id'],taxon=r['taxon'],d=r['d'],v=r['v'],d_source=r['d_evidence'],v_source=r['v_evidence'],included=included,reason=reason,callable_markers=int(np.isfinite(v).sum()),total_markers=len(sites),path=paths.get(s,('',''))[0]))
  if s not in paths:continue
  path,kind=paths[s];indexed=any(Path(path+suffix).is_file() for suffix in ('.csi','.tbi'));records={}
  with pysam.VariantFile(path) as vf:
   chrom='Hmr15' if kind=='legacy_VCF' else str(sites.contig.iloc[0])
   if chrom in vf.header.contigs:
    for rec in (vf.fetch(chrom,min(positions)-1,max(positions)) if indexed else iter(vf)):
     if rec.contig==chrom and rec.pos in positions:records.setdefault(rec.pos,[]).append(rec.copy())
   for j,marker in enumerate(sites.itertuples()):
    exact=[rr for rr in records.get(int(marker.pos1),[]) if rr.ref==marker.ref and tuple(rr.alts or ())==(marker.alt,)]
    if len(exact)>1:raise RuntimeError('Duplicated exact VCF marker')
    rec=exact[0] if exact else None;c=rec.samples[s] if rec else None
    gt=c.get('GT') if c is not None else None
    raw.append(dict(sample=s,taxon=r['taxon'],d=r['d'],v=r['v'],included=included,chrom=str(marker.contig),pos1=int(marker.pos1),ref=marker.ref,alt=marker.alt,rosina_base=marker.rosina_base,raw_GT='/'.join('.' if x is None else str(x) for x in (gt or ())),DP=c.get('DP') if c is not None else None,GQ=c.get('GQ') if c is not None else None,AD=str(c.get('AD')) if c is not None else '',PL=str(c.get('PL')) if c is not None else '',qualified_Rosina_dosage=v[j],record_state='exact_REF_ALT' if rec else 'no_exact_REF_ALT_record',path=path))
 tab(out/'cydno_individual_inclusion.tsv',inventory);tab(out/'cydno_individual_genotypes.tsv',raw)
 tab(out/'VCF_extraction_quality.tsv',quality)
 labels=['TNLeft','TNRight'];matrix=[sites.left_rosina_dosage.to_numpy(float),sites.right_rosina_dosage.to_numpy(float)];styles=['focal','focal']
 # Show melpomene controls as explicit cross-species anchors.
 rows=[]
 def mean(ids):
  if not ids:return np.full(len(sites),np.nan),np.zeros(len(sites),dtype=int)
  arr=np.array([get(s) for s in ids]);n=np.isfinite(arr).sum(0)
  return np.divide(np.nansum(arr,0),n,out=np.full(len(sites),np.nan),where=n>0),n
 for state in (1,0):
  ss=[r['sample'] for r in active if 'melpomene' in str(r['taxon']).lower() and int(float(r['v']))==state];value,n=mean(ss)
  labels.append(f'H. melpomene controls | v{"+" if state else "-"} (n={len(ss)})');matrix.append(value);styles.append('positive' if state else 'negative')
 groups=defaultdict(list)
 for r in inventory:
  if r['included']:groups[(r['taxon'],r['d'],r['v'])].append(r['sample'])
 for (taxon,d,v),ss in sorted(groups.items()):
  value,n=mean(ss);labels.append(f'{taxon.replace("Heliconius ","H. ").replace("ssp. ","")} | d{d}/v{v} (n={len(ss)})');matrix.append(value);styles.append('positive' if v=='1' else 'negative')
  for j,pos in enumerate(positions):rows.append(dict(taxon=taxon,d=d,v=v,n_individuals=len(ss),pos1=pos,n_called=int(n[j]),Rosina_frequency=value[j]))
 tab(out/'cydno_race_site_frequencies.tsv',rows)
 import matplotlib;matplotlib.use('Agg');matplotlib.rcParams.update({'pdf.fonttype':42,'font.family':'DejaVu Sans'})
 from matplotlib.backends.backend_pdf import PdfPages
 with PdfPages(out/'Ybv_cydno_genotypes_Rosina.pdf') as pdf:
  helper.draw(np.asarray(matrix),labels,styles,positions,out,pdf)
  # Unknown/unselected/no-call specimens remain visible as grey rows.
  ilabels=['TNLeft','TNRight'];imatrix=matrix[:2];istyles=['focal','focal']
  for r in inventory:
   text=f'{r["taxon"].replace("Heliconius ","H. ").replace("ssp. ","")} | {r["sample"]} | d{r["d"] or "?"}/v{r["v"] or "?"}'
   if not r['included']:text+=' [out]'
   ilabels.append(text);imatrix.append(get(r['sample']));istyles.append('positive' if r['v']=='1' else 'negative' if r['v']=='0' else 'unknown')
  helper.draw(np.asarray(imatrix),ilabels,istyles,positions,out,pdf,individual=True)
 audit=dict(status='COMPLETE',source=str(source),source_complete_sha256=sha(source/'complete.json'),focal_sha256=sha(focal),unchanged_focal_calls_and_markers=True,pooled_frequencies_reproduced=True,candidates=len(inventory),currently_included=sum(r['included'] for r in inventory),counts_by_taxon={taxon:dict(candidates=len([r for r in inventory if r['taxon']==taxon]),included=sum(r['included'] for r in inventory if r['taxon']==taxon)) for taxon in sorted({r['taxon'] for r in inventory})},literature=['https://doi.org/10.1111/mec.14277','https://doi.org/10.7554/eLife.68549'],reference_base_not_asserted_causal=True,no_downloads_or_calling=True,individual_page_is_dosage_not_race_frequency=True,PDF=str(out/'Ybv_cydno_genotypes_Rosina.pdf'))
 (out/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
 print('SUMMARY\t'+json.dumps(audit))
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,default=SOURCE);p.add_argument('--output',type=Path,required=True);p.add_argument('--submit',action='store_true');a=p.parse_args()
 if a.submit:
  stage=ROOT/'analysis/yb_cydno_code'/sha(__file__)[:12];stage.mkdir(parents=True,exist_ok=True)
  for name in ('check_Ybv_cydno_genotypes.py','plot_Ybv_race_frequencies.py'):shutil.copy2(Path(__file__).with_name(name),stage/name)
  args=[sys.executable,str(stage/'check_Ybv_cydno_genotypes.py'),'--source',str(a.source),'--output',str(a.output)]
  cmd=['sbatch','--parsable','--partition=nano','--cpus-per-task=1','--mem=8G','--time=00:30:00','--job-name=Ybv_cydno','--output='+str(stage/'check.%j.out'),'--error='+str(stage/'check.%j.err'),'--wrap='+shlex.join(args)]
  print('JOB\t'+subprocess.check_output(cmd,text=True).strip());print('EXPECTED_PDF\t'+str(a.output/'Ybv_cydno_genotypes_Rosina.pdf'))
 else:build(a)
if __name__=='__main__':main()
