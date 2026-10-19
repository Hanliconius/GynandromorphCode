import importlib.util,json,tempfile,os,time
from pathlib import Path
s=importlib.util.spec_from_file_location('recover',Path(__file__).with_name('recover_Ybv_from_completed_reads.py'));m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
with tempfile.TemporaryDirectory() as tmp:
 root=Path(tmp);d=root/'raw/ERR1';d.mkdir(parents=True);t=dict(sample='ERR1')
 assert m.ready_reads(root,t)[0] is False
 for n in (1,2):(d/f'ERR1_{n}.fastq.gz').write_bytes(b'valid-size')
 marker=d/'complete.json'
 marker.write_text(json.dumps(dict(sample='ERR1',mates=[dict(mate=n,bytes=10,md5='recorded_download_checksum') for n in (1,2)])))
 assert m.ready_reads(root,t)[0] is True
 (d/'ERR1_2.fastq.gz').write_bytes(b'changed-size');assert m.ready_reads(root,t)[0] is False
 assert m.digest(d/'ERR1_1.fastq.gz')==m.digest(d/'ERR1_1.fastq.gz')
 print('PASS: completed input marker required; both mates and recorded sizes checked; changed reads rejected')
