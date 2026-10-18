from pathlib import Path
import importlib.util,tempfile
spec=importlib.util.spec_from_file_location('focus',Path(__file__).with_name('focus_Ybv_capture.py'));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
assert m.expression([4,0,1,2,7,8])=='0-2,4,7-8'
with tempfile.TemporaryDirectory() as tmp:
 root=Path(tmp);(root/'raw/done').mkdir(parents=True);(root/'raw/done/complete.json').write_text('{}')
 samples=['capture','pending','running','done'];tasks=[dict(sample=s,index=i) for i,s in enumerate(samples)]
 seq={s:dict(library_strategy='WGS',library_selection='RANDOM') for s in samples}
 seq['capture']=dict(library_strategy='Targeted-Capture',library_selection='Hybrid Selection')
 plan=m.decisions(tasks,seq,{'10_2':'RUNNING'},root,dict(download=10))
 assert [r['keep'] for r in plan]==[True,False,True,True]
 assert plan[1]['reason']=='deferred_WGS_no_new_download'
 assert plan[2]['reason']=='WGS_transfer_already_running'
 print('PASS: retain capture and available/in-flight WGS, defer pending WGS, preserve array indices')
