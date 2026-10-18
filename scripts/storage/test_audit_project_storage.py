#!/usr/bin/env python3
import importlib.util,json,os,tempfile,time
from pathlib import Path
from types import SimpleNamespace
spec=importlib.util.spec_from_file_location('storage',Path(__file__).with_name('audit_project_storage.py'));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
def read_tsv(p):
 import csv
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
with tempfile.TemporaryDirectory() as temp:
 base=Path(temp);root=base/'project';root.mkdir();live=root/'analysis/live';live.mkdir(parents=True);old=root/'analysis/old';old.mkdir();inputs=root/'data';inputs.mkdir()
 protected_input=inputs/'source.position.bam';protected_input.write_bytes(b'protected')
 (live/'run.json').write_text(json.dumps({'root':str(root),'source':str(old)}));(old/'run.json').write_text(json.dumps({'manifest':str(old/'inputs.tsv')}));(old/'inputs.tsv').write_text('input\n'+str(protected_input)+'\n')
 work=root/'analysis/stale/tmp';work.mkdir(parents=True);intermediate=work/'sample.fixmate.bam';intermediate.write_bytes(b'intermediate')
 hard=inputs/'a.dat';hard.write_bytes(b'hardlink bytes');os.link(hard,inputs/'b.dat')
 outside=base/'outside.dat';outside.write_bytes(b'x'*5000);(inputs/'outside_link').symlink_to(outside);(inputs/'broken_link').symlink_to(base/'missing')
 stamp=time.time()-3*86400
 for p in (protected_input,intermediate):os.utime(p,(stamp,stamp))
 for name in ('candidateA.dat','candidateB.dat'):
  with (inputs/name).open('wb') as h:h.truncate(34*1024*1024)
  os.utime(inputs/name,(stamp,stamp))
 # Different internal contents with identical edges must not be called proven duplicates.
 with (inputs/'candidateB.dat').open('r+b') as h:h.seek(2*1024*1024);h.write(b'changed middle')
 os.utime(inputs/'candidateB.dat',(stamp,stamp))
 m.live_jobs=lambda root,exclude:([],set(),[],'')
 out=root/'analysis/storage_audits/test';m.audit(SimpleNamespace(root=root,output=out,protect_run=[live]))
 summary=json.loads((out/'complete.json').read_text());assert summary['protected_paths']>=3
 files=read_tsv(out/'largest_files.tsv');assert any(r['path']==str(protected_input) and r['protected']=='1' for r in files)
 candidates=read_tsv(out/'cleanup_candidates.tsv');assert next(r for r in candidates if r['path']==str(protected_input))['status']=='protected_live_input';assert next(r for r in candidates if r['path']==str(intermediate))['status']=='candidate_for_review'
 dup=read_tsv(out/'duplicate_candidates.tsv');assert len(dup)==2 and all('NOT_proven_duplicate' in r['evidence'] for r in dup)
 syms=json.loads((out/'symlink_summary.json').read_text());assert syms['broken_links']==1 and syms['targets_outside_project']==2
 assert not any(str(out) in r['path'] for r in files)
 assert intermediate.read_bytes()==b'intermediate' and protected_input.read_bytes()==b'protected'
 assert outside.stat().st_size==5000
 print('PASS: recursive live-input protection, immutable data, symlink handling, report exclusion, conservative duplicate evidence')
