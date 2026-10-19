#!/usr/bin/env python3
"""Trace phenotype evidence, BAM coverage, exact marker calls and plotted membership.
Read-only with respect to all source data. No downloads, alignments or new calls.
"""
import argparse,csv,hashlib,importlib.util,json,os,re,shlex,shutil,subprocess,tempfile
from collections import Counter,defaultdict
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
LIT='https://pmc.ncbi.nlm.nih.gov/articles/PMC8289415/'
def read(p):
 if not p.is_file():return []
 if p.suffix.lower() in {'.xlsx','.xls'}:
  import pandas as pd
  sheets=pd.read_excel(p,sheet_name=None,dtype=str)
  return [dict(frame_row,source_sheet=sheet) for sheet,frame in sheets.items() for frame_row in frame.fillna('').to_dict('records')]
 with p.open() as h:return list(csv.DictReader(h,delimiter='\t' if p.suffix!='.csv' else ','))
def tab(p,rr):
 if not rr:p.write_text('status\nNO_ROWS\n');return
 cols=list(dict.fromkeys(k for r in rr for k in r))
 with p.open('w',newline='') as h:
  w=csv.DictWriter(h,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(rr)
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def binary(v):
 try:return int(float(v)) if float(v) in (0,1) else None
 except (TypeError,ValueError):return None
def ids(text):return set(re.findall(r'(?:ERR|SRR|DRR|SAMEA|SAMN|ERS)\d+',str(text)))
def resolve(evidence,taxon,interval):
 vals={binary(r.get('value')) for r in evidence if r['interval']==interval}-{None}
 # Apply only explicitly documented surface-specific states, never d->v copying.
 text=taxon.lower();override=None;why=''
 if 'cythera' in text or 'vulcanus' in text:override=0 if interval=='d' else 1;why='documented_surface_specific_melpomene'
 elif 'weymeri' in text or 'pachinus' in text:override=1;why='Livraghi2021_both_surface_bar_comparison'
 elif 'chioneus' in text or 'zelinde' in text:override=0;why='Livraghi2021_both_surface_bar_absent_comparison'
 if override is not None:return override,why,bool(vals-{override})
 return (next(iter(vals)),'recorded_surface_phenotype',False) if len(vals)==1 else (None,'conflicting_records' if vals else 'no_recorded_surface_phenotype',len(vals)>1)
def prepare(a):
 root=a.root;run=a.run;out=a.output;out.mkdir(parents=True,exist_ok=False)
 for n in ('code','tasks','logs'):(out/n).mkdir()
 shutil.copy2(__file__,out/'code/reconcile_Yb_individuals.py')
 cfg=json.loads((run/'run.json').read_text());meta=root/'data/derived/region_panel_metadata.tsv'
 # Canonical sheet plus small historically recorded metadata/phenotype tables.
 sourcepaths={meta}
 for directory in (root/'data',root/'metadata',root/'external/krista_genotype_plots_original',root/'analysis/shareable_metadata_for_Arnaud'):
  if not directory.exists():continue
  for parent,dirs,names in os.walk(directory):
   dirs[:]=[d for d in dirs if d not in {'.git','logs','bam','bams','fastq','raw','toolenv'}]
   for name in names:
    p=Path(parent)/name
    if p.suffix.lower() in {'.tsv','.csv','.xlsx','.xls'} and re.search('metadata|phenotype|sample',name,re.I) and p.stat().st_size<10_000_000:sourcepaths.add(p)
 evidence=[];taxa=defaultdict(set);errors=[]
 for p in sorted(sourcepaths):
  try:
   for number,r in enumerate(read(p),2):
    sampleids=set()
    for key,value in r.items():
     if re.search('accession|sequenceid|earthcapeid|sample|run',key,re.I):sampleids.update(ids(value))
    taxon=next((r[k] for k in r if re.fullmatch('taxon.?name|taxon|scientific.?name',k,re.I)), '')
    for s in sampleids:
     if taxon:taxa[s].add(' '.join(taxon.split()))
     for key,value in r.items():
      norm=re.sub('[^a-z]','',key.lower())
      interval='d' if norm in {'ybdorsal','ybd','dorsalyb'} else 'v' if norm in {'ybventral','ybv','ventralyb'} else None
      if interval and binary(value) is not None:evidence.append(dict(sample=s,interval=interval,value=binary(value),source=str(p),line=number,column=key,taxon=taxon))
  except Exception as e:errors.append(dict(path=str(p),error=str(e)))
 for entry in json.loads((Path(cfg['source'])/'run.json').read_text()).get('entries',[]):
  if entry.get('group'):taxa[entry['sample']].add(str(entry['group']))
 members=read(run/'control_membership.tsv');old=read(Path(cfg['parent'])/'control_membership.tsv')
 selected=read(run/'selected_sources.tsv')+read(Path(cfg['parent'])/'selected_sources.tsv')
 sampleids=set(taxa)|{t['sample'] for t in cfg['new_tasks']}|{r['sample'] for r in members+old}
 markers={}
 import pandas as pd
 for ref in ('HmRos','Hmel2.5'):
  for interval in ('d','v'):
   path=run/'frozen'/f'{ref}_Yb{interval}_unchanged_focal_calls.tsv.gz'
   frame=pd.read_csv(path,sep='\t').sort_values('pos1')
   markers[ref+':'+interval]=[dict(contig=str(r.contig),pos=int(r.pos1),ref=str(r.ref),alt=str(r.alt)) for r in frame.itertuples()]
 tasks=[];seen=set()
 def add(kind,path,reference):
  p=Path(path)
  if not p.is_file():return
  key=(kind,str(p.resolve()),reference)
  if key in seen:return
  seen.add(key);tasks.append(dict(task=len(tasks),kind=kind,path=str(p),reference=reference))
 for r in read(root/'analysis/yb_existing_resources_v1/bam_candidates.tsv'):
  if r.get('reference') in {'HmRos','Hmel2.5'} and '/tmp/' not in r.get('resolved_path','') and not re.search(r'\.(name|fixmate|position)\.bam$',r.get('resolved_path','')):add('bam',r['resolved_path'],r['reference'])
 for p in (run/'bam').glob('*.bam'):add('bam',p,'HmRos')
 for t in cfg['tasks']:
  if t.get('bam'):add('bam',t['bam'],t['reference'])
 for r in selected:add('vcf',r['path'],r['reference'])
 for t in cfg['tasks']:add('vcf',run/'calls'/f'{t["task"]}.vcf.gz',t['reference'])
 bio={}
 for p in sorted((root/'analysis/sequencing_audits').glob('*/sequencing_and_reuse.tsv')):
  for r in read(p):
   if r.get('sample_accession'):bio[r['sample']]=r['sample_accession']
 phen=[]
 for s in sorted(sampleids):
  taxon='; '.join(sorted(taxa[s]));ee=[r for r in evidence if r['sample']==s]
  d,ds,dc=resolve(ee,taxon,'d');v,vs,vc=resolve(ee,taxon,'v')
  exclusion='PK_or_laboratory_cross_or_focal' if re.search('piano|bwpk|PK_|gynandromorph|hybrid|backcross|\\bf1\\b|\\bf2\\b',taxon,re.I) else ''
  phen.append(dict(sample=s,individual_id=bio.get(s,s),taxon=taxon,d=d,v=v,d_evidence=ds,v_evidence=vs,phenotype_conflict=dc or vc,phenotype_class=f'd{d}/v{v}',cydno_priority=bool(re.search('cydno|pachinus',taxon,re.I)),exclusion=exclusion))
 tab(out/'phenotype_evidence.tsv',evidence);tab(out/'resolved_phenotypes.tsv',phen);tab(out/'input_read_errors.tsv',errors)
 bundle=dict(root=str(root),run=str(run),output=str(out),markers=markers,tasks=tasks,phenotypes=phen,old_members=old,members=members,repo=str(a.push_repo) if a.push_repo else '',literature=LIT,source_tables=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(sourcepaths)],no_new_downloads_or_alignments_or_calls=True)
 save(out/'run.json',bundle)
 if a.submit:
  script=out/'code/reconcile_Yb_individuals.py'
  def sb(stage,extra):
   wrap=shlex.join(['micromamba','run','-n','gynandrocalling','python',str(script),'--stage',stage,'--output',str(out)])
   cmd=['sbatch','--parsable','--partition=cpu','--cpus-per-task=1','--mem=6G','--time=02:00:00','--job-name=Yb_ledger_'+stage,'--output='+str(out/'logs'/f'{stage}.%A_%a.out'),'--error='+str(out/'logs'/f'{stage}.%A_%a.err'),*extra,'--wrap='+wrap]
   return subprocess.check_output(cmd,text=True).strip().split(';')[0]
  job=sb('worker',['--array=0-'+str(len(tasks)-1)+'%12']);report=sb('report',['--dependency=afterany:'+job])
  save(out/'jobs.json',dict(evidence_array=job,report=report));print('EVIDENCE_ARRAY\t'+job+'\nREPORT\t'+report)
 print('OUTPUT\t'+str(out));print('EXISTING_RESOURCE_TASKS\t'+str(len(tasks)))
def worker(cfg):
 import pysam,math
 out=Path(cfg['output']);task=cfg['tasks'][int(os.environ['SLURM_ARRAY_TASK_ID'])];p=Path(task['path']);records=[];ref=task['reference'];error=''
 try:
  if task['kind']=='bam':
   with pysam.AlignmentFile(p,'rb') as bam:
    samples=sorted({r.get('SM','') for r in bam.header.to_dict().get('RG',[]) if r.get('SM')})
    if len(samples)!=1:raise RuntimeError('Not a single sample BAM: '+repr(samples))
    sample=samples[0]
    if not bam.has_index():raise RuntimeError('BAM index absent')
    for interval in ('d','v'):
     mm=cfg['markers'][ref+':'+interval];chrom=mm[0]['contig'];lo=min(r['pos'] for r in mm)-1;hi=max(r['pos'] for r in mm)
     if chrom not in bam.references:raise RuntimeError('Expected native contig absent: '+chrom)
     count=bam.count(chrom,lo,hi,read_callback='all')
     cov=bam.count_coverage(chrom,lo,hi,quality_threshold=20,read_callback=lambda r:not(r.is_unmapped or r.is_secondary or r.is_qcfail or r.is_duplicate) and r.mapping_quality>=30)
     depths=[sum(c[r['pos']-1-lo] for c in cov) for r in mm]
     records.append(dict(sample=sample,reference=ref,interval=interval,kind='bam',path=str(p),reads_overlapping=count,markers_DP2=sum(d>=2 for d in depths),markers_total=len(mm),mean_marker_DP=sum(depths)/len(depths),note='BAM SM identity and native contig required; MAPQ30 baseQ20; no genotype inference'))
  else:
   with pysam.VariantFile(str(p)) as vf:
    for interval in ('d','v'):
     mm=cfg['markers'][ref+':'+interval];chrom='Hmr15' if p.name.startswith('Ivory_Hmr15') else mm[0]['contig'];markers={r['pos']:r for r in mm}
     if chrom not in vf.header.contigs:continue
     lo=min(markers)-1;hi=max(markers)
     counters={s:Counter() for s in vf.header.samples}
     indexed=any(Path(str(p)+suffix).is_file() for suffix in ('.csi','.tbi'))
     it=vf.fetch(chrom,lo,hi) if indexed else iter(vf)
     for rec in it:
      if rec.contig!=chrom or rec.pos not in markers:continue
      wanted=markers[rec.pos];exact=rec.ref==wanted['ref'] and tuple(rec.alts or ())==(wanted['alt'],)
      for s,c in counters.items():
       c['position_records']+=1
       if not exact:c['allele_mismatch']+=1;continue
       c['exact_records']+=1;call=rec.samples[s];gt=call.get('GT') or ();dp=call.get('DP') or 0
       if dp>=2:c['DP2']+=1
       if len(gt)!=2 or any(x not in (0,1) for x in gt) or dp<2:c['missing_or_low_DP']+=1;continue
       q=call.get('GQ')
       if q is None:
        pl=call.get('PL')
        if pl and len(pl)==3 and all(x is not None for x in pl):
         w=[10**(-(x-min(pl))/10) for x in pl];state=sum(gt)
         if w[state]<max(w)-1e-10:c['GT_PL_disagree']+=1;continue
         err=1-w[state]/sum(w);q=-10*math.log10(max(1e-99,err));c['PL_derived_Q']+=1
        else:c['GQ_and_PL_missing']+=1;continue
       if q>=20:c['callable']+=1
       else:c['low_quality']+=1
     for s,c in counters.items():records.append(dict(sample=s,reference=ref,interval=interval,kind='vcf',path=str(p),markers_total=len(mm),**dict(c)))
 except Exception as e:error=str(e)
 save(out/'tasks'/f'{task["task"]:04d}.json',dict(task=task,records=records,error=error))
def report(cfg):
 out=Path(cfg['output']);allrows=[];errors=[]
 for t in cfg['tasks']:
  p=out/'tasks'/f'{t["task"]:04d}.json'
  if not p.exists():errors.append(dict(path=t['path'],error='worker_output_missing'));continue
  x=json.loads(p.read_text());allrows.extend(x['records'])
  if x['error']:errors.append(dict(path=t['path'],error=x['error']))
 tab(out/'resource_evidence.tsv',allrows);tab(out/'resource_errors.tsv',errors)
 lookup=defaultdict(list)
 for r in allrows:lookup[(r['sample'],r['reference'],r['interval'])].append(r)
 prior=defaultdict(list)
 for r in cfg['members']:prior[(r['sample'],r['reference'],r['interval'])].append(r)
 ledger=[]
 for ph in cfg['phenotypes']:
  for ref in ('HmRos','Hmel2.5'):
   for interval in ('d','v'):
    rr=lookup[(ph['sample'],ref,interval)];called=max([r.get('callable',0) for r in rr if r['kind']=='vcf'] or [0]);covered=max([r.get('markers_DP2',0) for r in rr if r['kind']=='bam'] or [0])
    plotted=any(str(r.get('included','')).lower() in {'true','1'} for r in prior[(ph['sample'],ref,interval)])
    action='excluded_PK_or_cross' if ph['exclusion'] else 'resolve_phenotype_conflict' if ph['phenotype_conflict'] and ph[interval] is None else 'resolve_surface_phenotype' if ph[interval] is None else 'already_plotted' if plotted else 'restore_existing_callable_VCF_to_panel' if called else 'recount_existing_BAM_at_frozen_markers' if covered else 'inspect_available_BAM_no_marker_DP2' if any(r['kind']=='bam' for r in rr) else 'inspect_existing_VCF_filter_and_allele_evidence' if rr else 'trace_historical_resource_provenance_not_proven_absent'
    ledger.append(dict(**ph,reference=ref,interval=interval,existing_callable_markers=called,existing_BAM_markers_DP2=covered,currently_plotted=plotted,action=action,existing_paths=';'.join(sorted({r['path'] for r in rr}))))
 tab(out/'individual_interval_ledger.tsv',ledger);tab(out/'cydno_priority.tsv',[r for r in ledger if r['cydno_priority']]);tab(out/'restore_or_recount.tsv',[r for r in ledger if r['action'] in {'restore_existing_callable_VCF_to_panel','recount_existing_BAM_at_frozen_markers'}])
 # Actual Slurm failures and current run status, not silently interpreted as missing coverage.
 run=Path(cfg['run']);runcfg=json.loads((run/'run.json').read_text());latest=sorted((run/'capture_focus').glob('*/run.json'))
 selectedcfg=json.loads(latest[-1].read_text()) if latest else runcfg
 stages=[]
 for t in selectedcfg['new_tasks']:
  s=t['sample'];i=t['index'];call=run/'calls'/f'{t["call_task"]}.complete.json'
  stages.append(dict(sample=s,index=i,download_complete=(run/'raw'/s/'complete.json').is_file(),BAM_complete=(run/'bam'/f'{s}.complete.json').is_file(),call_complete=call.is_file(),log_evidence=';'.join(str(p) for p in sorted((run/'logs').glob(f'*.*_{i}.err')))))
 tab(out/'expansion_stage_accounting.tsv',stages)
 examples=[]
 for r in [r for r in stages if not r['call_complete']][:12]:
  for stage in ('download','align','call'):
   pp=sorted((run/'logs').glob(f"{stage}.*_{r['index']}.err"))
   if pp:examples.append(dict(sample=r['sample'],index=r['index'],stage=stage,path=str(pp[-1]),last_lines='\n'.join(pp[-1].read_text(errors='replace').splitlines()[-18:])))
 tab(out/'expansion_log_examples.tsv',examples)
 for p in (run/'complete.json',run/'failed_expansion_tasks.tsv',run/'source_harmonization.json'):
  if p.is_file():shutil.copy2(p,out/('source_'+p.name))
 counts=Counter(r['action'] for r in ledger if r['reference']=='HmRos' and r['interval']=='v')
 summary=dict(status='LEDGER_REQUIRES_REVIEW',phenotype_individual_ids=len({r['individual_id'] for r in cfg['phenotypes']}),resource_tasks=len(cfg['tasks']),read_errors=len(errors),HmRos_ventral_actions=dict(counts),surface_classes=dict(Counter(r['phenotype_class'] for r in cfg['phenotypes'] if not r['exclusion'])),stage_completion=dict(download=sum(r['download_complete'] for r in stages),BAM=sum(r['BAM_complete'] for r in stages),call=sum(r['call_complete'] for r in stages),retained_tasks=len(stages)),phenotypes_not_inferred_from_genotypes=True,no_new_downloads_or_alignments_or_calls=True,literature=cfg['literature'])
 save(out/'complete.json',summary)
 if cfg['repo']:
  repo=Path(cfg['repo']);branch='yb-individual-ledger-'+out.name;env=dict(os.environ,GIT_SSH_COMMAND='ssh -o BatchMode=yes')
  def git(*args,cwd=repo):subprocess.run(['git','-C',str(cwd),*map(str,args)],env=env,check=True,stdout=subprocess.DEVNULL)
  with tempfile.TemporaryDirectory(prefix='yb-ledger.') as tmp:
   work=Path(tmp)/'checkout';git('fetch','--quiet','origin','main');git('worktree','add','--quiet','--detach',work,'origin/main')
   try:
    dest=work/'results/yb_individual_ledger'/out.name;dest.mkdir(parents=True)
    for p in out.iterdir():
     if p.is_file() and p.name!='run.json':shutil.copy2(p,dest/p.name)
    git('add',str(dest.relative_to(work)),cwd=work);git('commit','--quiet','-m','Trace Yb surface phenotypes and existing interval evidence',cwd=work);git('push','--quiet','git@github.com:Hanliconius/GynandromorphCode.git','HEAD:refs/heads/'+branch,cwd=work)
   finally:git('worktree','remove','--force',work)
  print('GIT_BRANCH\t'+branch)
 print('SUMMARY\t'+json.dumps(summary))
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=ROOT);p.add_argument('--run',type=Path,default=ROOT/'analysis/yb_ventral_population_expansion_v1');p.add_argument('--output',type=Path,required=True);p.add_argument('--push-repo',type=Path);p.add_argument('--submit',action='store_true');p.add_argument('--stage',choices=['prepare','worker','report'],default='prepare');a=p.parse_args()
 if a.stage=='prepare':prepare(a)
 else:globals()[a.stage](json.loads((a.output/'run.json').read_text()))
if __name__=='__main__':main()
