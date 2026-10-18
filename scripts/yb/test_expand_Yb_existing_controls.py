import importlib.util,tempfile,json,gzip
from pathlib import Path
import numpy as np,pandas as pd,pysam
spec=importlib.util.spec_from_file_location('expand',Path(__file__).with_name('expand_Yb_existing_controls.py').resolve());m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
assert m.quality({'GT':(0,0),'DP':10,'GQ':None})[0] is None
assert m.quality({'GT':(0,0),'DP':10,'GQ':None,'PL':(0,40,80)},True)[0]==0
assert m.quality({'GT':(0,1),'DP':10,'GQ':None,'PL':(0,40,80)},True)[0] is None
assert m.quality({'GT':(0,2),'DP':10,'GQ':99})[0] is None
with tempfile.TemporaryDirectory() as td:
 root=Path(td);source=root/'source';source.mkdir();(root/'data/derived').mkdir(parents=True);audit=root/'analysis/yb_existing_resources_v1';audit.mkdir(parents=True)
 meta=[{'Accession':'ERR1','Taxon Name':'test rosina','Yb-dorsal':1,'Yb-ventral':1},{'Accession':'ERR2','Taxon Name':'test melpomene','Yb-dorsal':0,'Yb-ventral':0},{'Accession':'SRR1','Taxon Name':'Butterfly world piano keys','Yb-dorsal':1,'Yb-ventral':1}]
 pd.DataFrame(meta).to_csv(root/'data/derived/region_panel_metadata.tsv',sep='\t',index=False)
 m.save(source/'audit.json',{'status':'COMPLETE','original_SNP_selection_unchanged':True});m.save(source/'run.json',{'entries':[]})
 for ref,chrom in m.CHROMS.items():
  fa=root/'data'/f'{ref}.fa';fa.write_text('>'+chrom+'\n'+'A'*1700000+'\n');pysam.faidx(str(fa))
  for interval in ('d','v'):
   pos=m.BOUNDS[ref][interval][0]+1
   pd.DataFrame([dict(contig=chrom,pos1=pos,ref='A',alt='G',rosina_base='G',left_rosina_dosage=.5,right_rosina_dosage=0)]).to_csv(source/f'{ref}_Yb{interval}_unchanged_focal_calls.tsv.gz',sep='\t',index=False)
 def vcf(path,chrom,positions,names):
  path.parent.mkdir(parents=True,exist_ok=True)
  text='##fileformat=VCFv4.2\n##contig=<ID='+chrom+',length=1700000>\n'
  for k,num,typ in [('GT','1','String'),('DP','1','Integer'),('GQ','1','Integer')]:text+=f'##FORMAT=<ID={k},Number={num},Type={typ},Description="test">\n'
  text+='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(names)+'\n'
  for pos in positions:text+=f'{chrom}\t{pos}\t.\tA\tG\t.\tPASS\t.\tGT:DP:GQ\t'+'\t'.join('1/1:10:40' if n in ('ERR1','SRR1') else '0/0:10:40' for n in names)+'\n'
  path.write_text(text)
 legacy=root/'external/krista_genotype_plots_original/Genotype_Plots/Ivory_Hmr15-1627000-1637000_biallelic.vcf';vcf(legacy,'Hmr15',[1633042],['ERR1','ERR2','SRR1','TN_Left'])
 rows=[]
 for ref,chrom in m.CHROMS.items():
  native=root/f'{ref}.vcf';vcf(native,chrom,[m.BOUNDS[ref][i][0]+1 for i in ('d','v')],['ERR1','ERR2','SRR1'])
  for i in ('d','v'):
   for sample in (['ERR1'] if ref=='HmRos' and i=='d' else ['ERR1','ERR2','SRR1']):rows.append(dict(sample=sample,reference=ref,interval=i,path=str(native),callable_DP2_GQ20=1))
 # ERR2 dorsal is legacy-only; other sources are native.
 pd.DataFrame(rows).to_csv(audit/'sample_callability.tsv',sep='\t',index=False)
 pd.DataFrame(columns=['sample','reference','interval','reads_overlapping','resolved_path']).to_csv(audit/'bam_candidates.tsv',sep='\t',index=False)
 original_run=m.run
 def fake_run(cmd,**kwargs):
  if cmd[0]=='bgzip':
   f=str(cmd[-1]);pysam.tabix_compress(f,f+'.gz',force=True);Path(f).unlink()
  elif cmd[0]=='tabix':pysam.tabix_index(str(cmd[-1]),seq_col=0,start_col=1,end_col=1,force=True)
  else:return original_run(cmd,**kwargs)
 m.run=fake_run
 from types import SimpleNamespace
 out=root/'output';m.prepare(SimpleNamespace(root=root,output=out,source=source,rosina_reference=None,hmel_reference=None,submit=False))
 selected=pd.read_csv(out/'selected_sources.tsv',sep='\t');assert not selected['sample'].eq('SRR1').any();assert len(selected)==8
 assert selected.loc[(selected['sample']=='ERR2')&(selected.reference=='HmRos')&(selected.interval=='d'),'source'].iloc[0]=='legacy_VCF'
 m.report(json.loads((out/'run.json').read_text()))
 result=json.loads((out/'complete.json').read_text());assert result['focal_calls_unchanged']
 members=pd.read_csv(out/'control_membership.tsv',sep='\t');assert len(members)==8 and members.included.all()
 # Check REF mismatch aborts rather than merging by position.
 legacy.write_text(legacy.read_text().replace('\tA\tG\t','\tC\tG\t'))
 try:m.prepare(SimpleNamespace(root=root,output=root/'bad',source=source,rosina_reference=None,hmel_reference=None,submit=False))
 except RuntimeError as e:assert 'compatibility failed' in str(e)
 else:raise AssertionError('Reference mismatch accepted')
 print('PASS: native/legacy merge, PK exclusion, focal preservation, quality rules, all four panels, REF mismatch rejection')
