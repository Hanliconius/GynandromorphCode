#!/usr/bin/env python3
"""Audit old allele polarity, plot ALT dosage, and jointly discover full-window SNPs from existing BAMs."""
import argparse,csv,hashlib,importlib.util,json,os,shlex,shutil,subprocess,sys
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
CHROM='NC_140403.1';LO=1657079;HI=1660072
SOURCE=ROOT/'analysis/yb_completed_reads_recovery_v1'
def read(p):
 with Path(p).open() as h:return list(csv.DictReader(h,delimiter='\t'))
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def tab(p,rows):
 import pandas as pd
 pd.DataFrame(rows).to_csv(p,sep='\t',index=False)
def sb(out,stage,extra=()):
 cmd=[sys.executable,str(out/'rebuild_Ybv_REF_ALT.py'),'--output',str(out),'--stage',stage]
 return subprocess.check_output(['sbatch','--parsable','--partition=cpu','--cpus-per-task='+('8' if stage=='joint' else '1'),'--mem='+('32G' if stage=='joint' else '8G'),'--time=04:00:00','--job-name=Ybv_'+stage,'--output='+str(out/'logs'/f'{stage}.%j.out'),'--error='+str(out/'logs'/f'{stage}.%j.err'),*extra,'--wrap='+shlex.join(cmd)],text=True).strip().split(';')[0]
def prepare(a):
 out=a.output;out.mkdir(parents=True,exist_ok=False);(out/'logs').mkdir();shutil.copy2(__file__,out/Path(__file__).name)
 helper=Path(__file__).with_name('plot_Ybv_race_frequencies.py');shutil.copy2(helper,out/helper.name)
 save(out/'run.json',dict(source=str(a.source),output=str(out),region=f'{CHROM}:{LO}-{HI}',bounds=[LO,HI]))
 print('AUDIT\t'+sb(out,'audit'))
def audit(cfg):
 import pandas as pd,numpy as np,pysam
 source=Path(cfg['source']);out=Path(cfg['output']);sc=json.loads((source/'run.json').read_text());fa=sc['references']['HmRos']
 sites=pd.read_csv(source/'frozen/HmRos_Ybv_unchanged_focal_calls.tsv.gz',sep='\t').sort_values('pos1');rows=[]
 with pysam.FastaFile(fa) as f:
  for r in sites.itertuples():
   base=f.fetch(CHROM,int(r.pos1)-1,int(r.pos1)).upper()
   rows.append(dict(pos1=r.pos1,FASTA=base,REF=r.ref,ALT=r.alt,old_rosina_control_base=r.rosina_base,old_purple_is='ALT' if r.rosina_base==r.alt else 'REF',REF_matches_FASTA=base==r.ref))
 tab(out/'old_marker_polarity.tsv',rows)
 raw=[]
 for entry in read(source/'selected_sources.tsv'):
  if entry['reference']!='HmRos' or entry['interval']!='v' or 'rosina' not in entry['taxon'].lower():continue
  with pysam.VariantFile(entry['path']) as vf:
   wanted={int(r.pos1):r for r in sites.itertuples()}
   records=vf.fetch(CHROM,int(sites.pos1.min())-1,int(sites.pos1.max())) if vf.index is not None else iter(vf)
   for rec in records:
    if rec.contig!=CHROM or rec.pos not in wanted:continue
    r=wanted[rec.pos];c=rec.samples[entry['sample']]
    raw.append(dict(sample=entry['sample'],pos1=rec.pos,REF=rec.ref,ALT=','.join(rec.alts or ()),GT=str(c.get('GT')),DP=c.get('DP'),GQ=c.get('GQ'),AD=str(c.get('AD')),PL=str(c.get('PL')),old_rosina_control_base=r.rosina_base,source=entry['path']))
 tab(out/'rosina_raw_genotypes_old_markers.tsv',raw)

 if not all(r['REF_matches_FASTA'] for r in rows):raise RuntimeError('FASTA/REF mismatch; stopping before plotting or new calling')
 # Reconstruct and verify the existing full cohort first; then convert every
 # individual and both focal rows to raw ALT dosage without changing calls.
 code=(out/'plot_Ybv_race_frequencies.py').read_text();anchor=" for name,column in [('TNLeft','left_rosina_dosage'),('TNRight','right_rosina_dosage')]:";assert code.count(anchor)==1
 inject=''' flip=(sites.rosina_base.to_numpy()==sites.ref.to_numpy())
 for s in values:values[s][flip]=1-values[s][flip]
 for i in range(len(freq)):freq[i][flip]=1-freq[i][flip]
 for row in rows:
  if int(row['pos1']) in set(sites.loc[flip,'pos1']):row['frequency']=1-row['frequency']
 for column in ('left_rosina_dosage','right_rosina_dosage'):
  sites.loc[flip,column]=1-sites.loc[flip,column]
'''
 code=code.replace(anchor,inject+anchor).replace('Rosina allele frequency / focal dosage','ALT allele frequency / focal ALT dosage').replace('Yb(v) | frequencies by race','Yb(v) | ALT frequency, original marker set')
 ns={'__name__':'old_alt','__file__':str(out/'plot_Ybv_race_frequencies.py')};exec(compile(code,'ALT_display','exec'),ns)
 ns['build'](argparse.Namespace(source=source,baseline=ROOT/'analysis/yb_matched_cohort_v1',output=out/'old_markers_ALT'))
 aa=json.loads((out/'old_markers_ALT/audit.json').read_text());aa.update(allele_scale='0=REF, 1=ALT',original_cohort_verified_before_conversion=True,orientation_audit=str(out/'old_marker_polarity.tsv'));save(out/'old_markers_ALT/audit.json',aa)
 # Exact biological membership and all run aliases are preserved.
 membership=pd.read_csv(out/'old_markers_ALT/race_membership.tsv',sep='\t').fillna('');membership.to_csv(out/'membership.tsv',sep='\t',index=False)
 wanted={s for x in membership.aliases for s in x.split(';')};candidates={s:[] for s in wanted}
 for t in sc.get('tasks',[]):
  if t.get('reference')=='HmRos' and t.get('sample') in wanted:candidates[t['sample']].append(t['bam'])
 inventory=ROOT/'analysis/yb_existing_resources_v1/bam_candidates.tsv'
 for r in read(inventory):
  if r['sample'] in wanted and r['reference']=='HmRos':candidates[r['sample']].append(r['resolved_path'])
 with pysam.FastaFile(fa) as f:dictionary=dict(zip(f.references,f.lengths))
 decisions=[];bams=[]
 for sample in sorted(wanted):
  chosen=None;checks=[]
  for path in dict.fromkeys(candidates[sample]):
   p=Path(path)
   if not p.is_file():continue
   try:
    with pysam.AlignmentFile(str(p),'rb') as b:
     sm={r.get('SM') for r in b.header.to_dict().get('RG',[]) if r.get('SM')};bd=dict(zip(b.references,b.lengths))
     okay=b.has_index() and sm=={sample} and CHROM in bd and all(dictionary.get(k)==v for k,v in bd.items())
     if okay:chosen=str(p.resolve());break
     checks.append(dict(path=path,reason='index, sample, or dictionary incompatible'))
   except Exception as e:checks.append(dict(path=path,reason=str(e)))
  decisions.append(dict(sample=sample,action='recount_existing_alignment' if chosen else 'UNRESOLVED_NO_NEW_DOWNLOAD',bam=chosen or '',evidence=json.dumps(checks)))
  if chosen:bams.append(chosen)
 tab(out/'full_interval_reuse_manifest.tsv',decisions)
 cfg.update(reference=fa,bams=bams,unresolved=[r['sample'] for r in decisions if not r['bam']],old_marker_count=len(sites),individuals=len(membership));save(out/'run.json',cfg)
 if cfg['unresolved']:
  save(out/'audit.complete.json',dict(status='AUDIT_COMPLETE_FULL_INTERVAL_BLOCKED_BY_UNRESOLVED_RESOURCES',unresolved=cfg['unresolved'],no_downloads=True,old_markers_ALT_plot_complete=True));print('UNRESOLVED\t'+str(len(cfg['unresolved']))+'\nMANIFEST\t'+str(out/'full_interval_reuse_manifest.tsv'));return
 (out/'bams.list').write_text('\n'.join(bams)+'\n');save(out/'audit.complete.json',dict(status='COMPLETE',old_marker_count=len(sites),rosina_control_base_is_REF=sum(r['old_purple_is']=='REF' for r in rows),rosina_control_base_is_ALT=sum(r['old_purple_is']=='ALT' for r in rows),individuals=len(membership),BAMs=len(bams),no_downloads=True,no_alignments=True))
 j=sb(out,'joint');r=sb(out,'report',['--dependency=afterok:'+j]);save(out/'jobs.json',dict(joint=j,report=r));print('JOINT\t'+j+'\nREPORT\t'+r)
def joint(cfg):
 out=Path(cfg['output']);tmp=out/'full_window.partial.vcf.gz'
 mp=['bcftools','mpileup','--threads','8','-f',cfg['reference'],'-r',cfg['region'],'-b',str(out/'bams.list'),'-q','30','-Q','20','-d','10000','-a','FORMAT/AD,FORMAT/DP','-Ou']
 call=['bcftools','call','--threads','8','-m','-A','-f','GQ','-Oz','-o',str(tmp)]
 p=subprocess.Popen(mp,stdout=subprocess.PIPE);q=subprocess.run(call,stdin=p.stdout);p.stdout.close();rc=p.wait()
 if rc or q.returncode:raise RuntimeError('Joint calling failed')
 tmp.rename(out/'full_window.vcf.gz');subprocess.run(['bcftools','index','-f',str(out/'full_window.vcf.gz')],check=True);save(out/'joint.complete.json',dict(status='COMPLETE',commands=[mp,call],reference=cfg['reference']))
def report(cfg):
 import numpy as np,pandas as pd,pysam
 out=Path(cfg['output']);source=Path(cfg['source']);spec=importlib.util.spec_from_file_location('quality_engine',source/'engine.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 members=pd.read_csv(out/'membership.tsv',sep='\t').fillna('');sites=[];vectors=[];quality_counts={}
 with pysam.FastaFile(cfg['reference']) as fa,pysam.VariantFile(str(out/'full_window.vcf.gz')) as vf:
  for rec in vf:
   if len(rec.ref)!=1 or not rec.alts or len(rec.alts)!=1 or len(rec.alts[0])!=1:continue
   if fa.fetch(rec.contig,rec.pos-1,rec.pos).upper()!=rec.ref:raise RuntimeError('Joint REF mismatch')
   vector=[]
   for r in members.to_dict('records'):
    vals=[]
    for alias in r['aliases'].split(';'):
     val,kind=m.quality(rec.samples[alias],allow_pl=True);quality_counts[kind]=quality_counts.get(kind,0)+1
     if val is not None:vals.append(val)
    vector.append(vals[0] if vals and len(set(vals))==1 else np.nan)
   if not np.isfinite(vector).any():continue
   sites.append(dict(contig=rec.contig,pos1=rec.pos,ref=rec.ref,alt=rec.alts[0]));vectors.append(vector)
 if not sites:raise RuntimeError('No callable full-window biallelic SNPs')
 matrix=np.asarray(vectors).T;positions=[r['pos1'] for r in sites];pd.DataFrame(matrix,index=members['sample'],columns=positions).to_csv(out/'individual_ALT_dosages.tsv.gz',sep='\t',compression='gzip',na_rep='NA');tab(out/'SNPs.tsv',sites)
 rows=[];freq=[];labels=[];styles=[];callable_fractions=[]
 groups=[('All v+',members.index[members.v.eq(1)].tolist()),('All v-',members.index[members.v.eq(0)].tolist())]
 for (race,v),g in members.groupby(['race','v'],sort=True):groups.append((f'{race} | v{int(v)}',g.index.tolist()))
 for label,ids in groups:
  if not ids:continue
  x=matrix[ids];n=np.isfinite(x).sum(axis=0);f=np.divide(np.nansum(x,axis=0),n,out=np.full(len(sites),np.nan),where=n>0)
  labels.append(label+f' (n={len(ids)})');freq.append(f);styles.append('summary');callable_fractions.append(n/len(ids))
  for i,site in enumerate(sites):rows.append(dict(row=label,**site,n_individuals=len(ids),n_called=int(n[i]),ALT_frequency=f[i]))
 tab(out/'race_ALT_frequencies.tsv',rows)
 # Report every SNP, not only successful phenotype markers. Selection metrics
 # are exploratory in the same cohort; TN is excluded throughout discovery.
 metrics=[]
 for j,site in enumerate(sites):
  z=dict(site)
  for v in (1,0):
   ids=members.index[members.v.eq(v)].tolist();x=matrix[ids,j];z[f'n_v{v}']=int(np.isfinite(x).sum());z[f'callable_fraction_v{v}']=float(np.isfinite(x).mean());z[f'ALT_AF_v{v}']=float(np.nanmean(x)) if np.isfinite(x).any() else None
  z['delta_AF']=z['ALT_AF_v1']-z['ALT_AF_v0'] if z['ALT_AF_v1'] is not None and z['ALT_AF_v0'] is not None else None;metrics.append(z)
 tab(out/'all_SNP_phenotype_metrics.tsv',metrics)
 import matplotlib;matplotlib.use('Agg')
 spec=importlib.util.spec_from_file_location('display',out/'plot_Ybv_race_frequencies.py');display=importlib.util.module_from_spec(spec);spec.loader.exec_module(display)
 # Use the established layout but label its allele scale correctly.
 code=(out/'plot_Ybv_race_frequencies.py').read_text().replace('Rosina allele frequency / focal dosage','ALT allele frequency').replace('27 unchanged SNP columns','All callable biallelic full-window SNP columns').replace('Yb(v) | frequencies by race','Yb(v) | full-window ALT frequencies').replace(' / focal dosage','')
 ns={'__name__':'full_ALT','__file__':str(out/'plot_Ybv_race_frequencies.py')};exec(compile(code,'full_ALT','exec'),ns)
 from matplotlib.backends.backend_pdf import PdfPages
 with PdfPages(out/'Ybv_full_window_ALT_by_race.pdf') as pdf:
  ns['draw'](np.asarray(freq),labels,styles,positions,out,pdf)
  ns['draw'](np.asarray(callable_fractions),labels,styles,positions,out,pdf,True)
 save(out/'complete.json',dict(status='COMPLETE',individuals=len(members),SNPs=len(sites),scale='0=REF, 1=ALT',discovery='All callable biallelic SNPs in the window; no phenotype or focal selection',focal_new_SNP_calls='Not inferred: old focal calls only appear in original-marker ALT diagnostic',quality_evaluations=quality_counts,no_downloads=True,no_alignments=True,source=str(source)))
 print('STATUS\tCOMPLETE\nSNPS\t'+str(len(sites))+'\nPDF\t'+str(out/'Ybv_full_window_ALT_by_race.pdf'))
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,default=SOURCE);p.add_argument('--output',type=Path,default=ROOT/'analysis/ybv_REF_ALT_rebuild_v1');p.add_argument('--stage',choices=['prepare','audit','joint','report'],default='prepare');p.add_argument('--submit',action='store_true');a=p.parse_args()
 if a.stage=='prepare':prepare(a)
 else:globals()[a.stage](json.loads((a.output/'run.json').read_text()))
if __name__=='__main__':main()
