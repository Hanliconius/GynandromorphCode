#!/usr/bin/env python3
import gzip,importlib.util,json,os,tempfile,time
from pathlib import Path
spec=importlib.util.spec_from_file_location('verify',Path(__file__).with_name('verify_storage_candidates.py'));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
with tempfile.TemporaryDirectory() as tmp:
 root=Path(tmp);(root/'tasks').mkdir();a=root/'a';b=root/'b';a.write_bytes(b'abc');b.write_bytes(b'abc')
 stamp=time.time()-2*86400
 for p in (a,b):os.utime(p,(stamp,stamp))
 def task():
  paths=list(map(str,(a,b)));snap={}
  for p in (a,b):
   s=p.stat();snap[str(p)]=dict(bytes=s.st_size,mtime=s.st_mtime,dev=s.st_dev,inode=s.st_ino,nlink=s.st_nlink,allocated=s.st_blocks*512,protected=False)
  return dict(kind='duplicate_files',sample='test',paths=paths,snapshot=snap)
 cfg=dict(root=str(root),output=str(root),protected_runs=[],tasks=[task()])
 m.audit.live_jobs=lambda *args:([],set(),[],'')
 os.environ['SLURM_ARRAY_TASK_ID']='0';m.worker(cfg)
 assert json.loads((root/'tasks/0000.json').read_text())['status']=='EXACT_DUPLICATES_REFERENCE_REPAIR_REQUIRED'
 b.write_bytes(b'xyz');os.utime(b,(stamp,stamp));cfg['tasks']=[task()];m.worker(cfg)
 assert json.loads((root/'tasks/0000.json').read_text())['status']=='BLOCKED'
 cfg['protected_runs']=[str(root/'a')];cfg['tasks']=[task()];m.worker(cfg)
 assert json.loads((root/'tasks/0000.json').read_text())['status']=='BLOCKED'
 fq=root/'reads.gz'
 with gzip.open(fq,'wb') as h:h.write(b'@ERR1.1\nAC\n+\n!!\n')
 assert m.fastq(fq,'ERR1')['records']==1
 with gzip.open(fq,'wb') as h:h.write(b'@ERR1.1\nAC\n+\n!\n')
 try:m.fastq(fq,'ERR1');raise AssertionError('Malformed FASTQ accepted')
 except RuntimeError:pass
 assert a.read_bytes()==b'abc' and b.read_bytes()==b'xyz'
 print('PASS: full checksums reject unequal content, live inputs block validation, FASTQ structural checks, no input changes')
