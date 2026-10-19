#!/usr/bin/env python3
"""Resolve collection localities for the exact plotted melpomene melpomene individuals."""
import argparse,csv,json,shlex,subprocess,sys,time,urllib.request,xml.etree.ElementTree as ET
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
def ena(accession,folder):
 def get(acc):
  p=folder/(acc+'.xml')
  if not p.exists():
   for attempt in range(3):
    try:
     req=urllib.request.Request('https://www.ebi.ac.uk/ena/browser/api/xml/'+acc,headers={'User-Agent':'Ybv-locality-audit/1.0'})
     with urllib.request.urlopen(req,timeout=45) as h:data=h.read()
     ET.fromstring(data);p.write_bytes(data);break
    except Exception:
     if attempt==2:raise
     time.sleep(attempt+1)
  return ET.fromstring(p.read_bytes())
 try:
  tree=get(accession);ref=tree.find('.//SAMPLE_DESCRIPTOR');sample=ref.get('accession') if ref is not None else accession
  if ref is not None:tree=get(sample)
  attrs={}
  for x in tree.findall('.//SAMPLE_ATTRIBUTE'):
   key=x.findtext('TAG','');val=x.findtext('VALUE','');attrs[key]=val
  return dict(accession=accession,ena_sample=sample,ena_title=tree.findtext('.//SAMPLE_NAME/SCIENTIFIC_NAME',''),attributes=attrs,error='')
 except Exception as e:return dict(accession=accession,error=str(e),attributes={})
def run(a):
 import pandas as pd
 out=a.output;out.mkdir(parents=True,exist_ok=False);raw=out/'ena';raw.mkdir()
 members=pd.read_csv(a.membership,sep='\t').fillna('')
 def target(t):return ' '.join(t.lower().replace('spp.','ssp.').split())=='heliconius melpomene ssp. melpomene' or ' '.join(t.lower().split())=='heliconius melpomene melpomene'
 subset=members[members.race.map(target)].copy()
 if subset.empty:raise RuntimeError('No exact melpomene melpomene row found')
 meta=pd.read_csv(ROOT/'data/derived/region_panel_metadata.tsv',sep='\t',dtype=str).fillna('')
 ids=sorted({x for r in subset.to_dict('records') for x in (r.get('aliases','') or r['sample']).split(';')})
 results=list(ThreadPoolExecutor(max_workers=6).map(lambda s:ena(s,raw),ids));lookup={r['accession']:r for r in results}
 rows=[]
 for r in subset.to_dict('records'):
  aliases=(r.get('aliases','') or r['sample']).split(';')
  matches=meta[meta.Accession.isin(aliases)].to_dict('records')
  public=[lookup[x] for x in aliases]
  rows.append(dict(sample=r['sample'],individual_id=r.get('individual_id',''),aliases=';'.join(aliases),taxon=r['race'],v=r['v'],local_metadata=json.dumps(matches,ensure_ascii=False),ENA_locality_attributes=json.dumps(public,ensure_ascii=False),geography_assignment='UNASSIGNED: inspect collection records'))
 pd.DataFrame(rows).to_csv(out/'melpomene_individual_localities.tsv',sep='\t',index=False)
 # A compact, readable report preserves all attributes separately in the TSV/XML.
 lines=[]
 for r in rows:
  lines.append('SAMPLE\t'+r['sample'])
  for m in json.loads(r['local_metadata']):
   lines.append('LOCAL\t'+' | '.join(f'{k}={v}' for k,v in m.items() if v and any(w in k.lower() for w in ('country','lat','long','cape','sequence','taxon','note'))))
  for e in json.loads(r['ENA_locality_attributes']):
   attrs=e['attributes'];loc={k:v for k,v in attrs.items() if any(w in k.lower() for w in ('geograph','country','local','lat','lon','collect','isolation','strain','sample name'))}
   lines.append('ENA\t'+e['accession']+'\t'+json.dumps(loc,ensure_ascii=False)+(('\tERROR='+e['error']) if e['error'] else ''))
 (out/'localities.txt').write_text('\n'.join(lines)+'\n')
 (out/'complete.json').write_text(json.dumps(dict(status='COMPLETE',individuals=len(rows),accessions=len(ids),lookup_errors=sum(bool(r['error']) for r in results),membership=str(a.membership),genotypes_changed=False,phenotypes_changed=False),indent=2)+'\n')
 print('INDIVIDUALS\t'+str(len(rows))+'\nREPORT\t'+str(out/'localities.txt')+'\nLOOKUP_ERRORS\t'+str(sum(bool(r['error']) for r in results)))
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--membership',type=Path,default=ROOT/'analysis/ybv_western_positive_all_negative_v1/western_positive_all_negative/race_membership.tsv');p.add_argument('--output',type=Path,default=ROOT/'analysis/ybv_melpomene_localities_v1');p.add_argument('--submit',action='store_true');a=p.parse_args()
 if not a.submit:return run(a)
 stage=ROOT/'analysis/yb_locality_code';stage.mkdir(parents=True,exist_ok=True)
 import hashlib,shutil
 frozen=stage/('audit_'+hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:12]+'.py');shutil.copy2(__file__,frozen)
 args=[sys.executable,str(frozen),'--membership',str(a.membership),'--output',str(a.output)]
 job=subprocess.check_output(['sbatch','--parsable','--partition=nano','--cpus-per-task=1','--mem=4G','--time=00:20:00','--job-name=Ybv_locality','--output='+str(stage/'audit.%j.out'),'--error='+str(stage/'audit.%j.err'),'--wrap='+shlex.join(args)],text=True).strip();print('JOB\t'+job)
if __name__=='__main__':main()
