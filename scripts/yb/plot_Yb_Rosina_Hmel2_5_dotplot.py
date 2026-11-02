#!/usr/bin/env python3
"""Replot the archived Yb assembly alignment over the current figure coordinates.

No new downloads, alignments, genotypes, BA calculations or marker selection.
PAF query: Hmel2.5; target: Rosina. All alignment records are retained.
"""
import argparse,csv,hashlib,json,re,shlex,shutil,subprocess,sys
from pathlib import Path
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis')
OUTPUT=ROOT/'Ybv_figure_assembly_dotplot_20261102'
ARCHIVE=Path('results/yb_exact_coordinate_transfer/job_138421491')
PAF='Hmel2.5_to_HmRos_Yb.paf'
FEATURES='Hmel2.5_to_HmRos_Yb_features.tsv'

def save(p,value):p.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write_tsv(p,rows,fields):
 with p.open('w') as h:
  w=csv.DictWriter(h,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
def region(name,length):
 m=re.fullmatch(r'(.+):(\d+)-(\d+)',name)
 if not m:raise ValueError('Extract header lacks documented region: '+name)
 contig,start,end=m.group(1),int(m.group(2)),int(m.group(3))
 if end-start+1!=length:raise ValueError('Extract header/length mismatch: '+name)
 return contig,start-1,end

def parse_paf(path):
 segments=[];gaps=[];blocks=[];names=set()
 for number,line in enumerate(path.read_text().splitlines(),1):
  if not line.strip():continue
  f=line.split('\t')
  if len(f)<12:raise ValueError(f'Malformed PAF record {number}')
  qn,ql,qs,qe,strand,tn,tl,ts,te,nmatch,alen,mq=f[:12]
  ql,qs,qe,tl,ts,te,nmatch,alen,mq=map(int,[ql,qs,qe,tl,ts,te,nmatch,alen,mq])
  qc,qo,qend=region(qn,ql);tc,to,tend=region(tn,tl);names.add((qc,qo,qend,tc,to,tend))
  if not 0<=qs<qe<=ql or not 0<=ts<te<=tl or strand not in ['+','-']:raise ValueError('Invalid PAF bounds')
  tags={s.split(':',2)[0]:s.split(':',2)[2] for s in f[12:] if len(s.split(':',2))==3}
  cigar=tags.get('cg','');ops=re.findall(r'(\d+)([MIDN=X])',cigar)
  if not ops or ''.join(n+op for n,op in ops)!=cigar:raise ValueError(f'Missing/unsupported CIGAR at record {number}')
  sign=1 if strand=='+' else -1;q=qs if sign==1 else qe;t=ts
  block=dict(block=number,strand=strand,mapq=mq,identity=nmatch/alen,
   Rosina_start0=to+ts,Rosina_end0=to+te,Hmel_start0=qo+qs,Hmel_end0=qo+qe)
  blocks.append(block)
  for text,op in ops:
   n=int(text);qnext=q+sign*n if op in 'MI=X' else q;tnext=t+n if op in 'MDN=X' else t
   if op in 'M=X':segments.append(dict(block=number,strand=strand,mapq=mq,identity=nmatch/alen,
    Rosina_start0=to+t,Rosina_end0=to+tnext,Hmel_start0=qo+q,Hmel_end0=qo+qnext))
   elif n>=20:gaps.append(dict(block=number,strand=strand,operation=op,length_bp=n,
    Rosina_start0=to+t,Rosina_end0=to+tnext,Hmel_start0=qo+q,Hmel_end0=qo+qnext))
   q,t=qnext,tnext
  if q!=(qe if sign==1 else qs) or t!=te:raise ValueError(f'CIGAR endpoints disagree with PAF record {number}')
 if len(names)!=1:raise ValueError('Expected exactly one documented query/target extraction pair')
 return segments,gaps,blocks,next(iter(names))

def clipped(segment,lo,hi):
 x0,x1=segment['Rosina_start0'],segment['Rosina_end0'];y0,y1=segment['Hmel_start0'],segment['Hmel_end0']
 a,b=max(x0,lo),min(x1,hi)
 if a>=b:return None
 return ((a,y0+(y1-y0)*(a-x0)/(x1-x0)),(b,y0+(y1-y0)*(b-x0)/(x1-x0)))

def prepare(a):
 repo=Path(__file__).resolve().parents[2];paf=a.paf or repo/ARCHIVE/PAF;features=repo/ARCHIVE/FEATURES
 for p in [paf,features]:
  if not p.is_file():raise FileNotFoundError('Missing archived alignment input: '+str(p))
 # Structural validation before submission; no data-driven plotting on login node.
 segments,gaps,blocks,names=parse_paf(paf)
 if a.output.exists():raise RuntimeError('Preserving existing output: '+str(a.output))
 a.output.mkdir(parents=True);(a.output/'logs').mkdir();shutil.copy2(__file__,a.output/'worker.py')
 save(a.output/'run.json',dict(paf=str(paf),features=str(features),source_job='138421491',
  hashes={str(p):digest(p) for p in [paf,features]},figure_span_Rosina_1based=[1630000,1692000],
  SNP_strip_context_Rosina_1based=[1656500,1660500],alignment_records=len(blocks),
  aligned_segments=len(segments),coordinate_convention='Exported geometry: 0-based interbase coordinates; axes: Mb',
  new_downloads=False,new_alignments=False,all_alignment_records_retained=True))
 cmd=[sys.executable,str(a.output/'worker.py'),'--run','--output',str(a.output)]
 job=subprocess.check_output(['sbatch','--parsable','--partition=nano','--time=00:10:00','--mem=4G','--cpus-per-task=1',
  '--job-name=Ybv_dotplot','--output='+str(a.output/'logs/plot.%j.out'),'--error='+str(a.output/'logs/plot.%j.err'),
  '--wrap='+shlex.join(cmd)],text=True).strip().split(';')[0]
 save(a.output/'jobs.json',dict(plot=job));print(f'JOB\t{job}\nOUTPUT\t{a.output}')

def run(out):
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from matplotlib.collections import LineCollection
 from matplotlib.lines import Line2D
 from matplotlib.ticker import MultipleLocator,FormatStrFormatter
 cfg=json.loads((out/'run.json').read_text());paf=Path(cfg['paf']);fp=Path(cfg['features'])
 for path,expected in cfg['hashes'].items():
  if digest(Path(path))!=expected:raise RuntimeError('Source changed since submission: '+path)
 segments,gaps,blocks,names=parse_paf(paf);qc,qo,qend,tc,to,tend=names
 fields=['block','strand','mapq','identity','Rosina_start0','Rosina_end0','Hmel_start0','Hmel_end0']
 write_tsv(out/'alignment_blocks.tsv',blocks,fields);write_tsv(out/'aligned_segments.tsv',segments,fields)
 write_tsv(out/'within_alignment_gaps_ge20bp.tsv',gaps,['block','strand','operation','length_bp','Rosina_start0','Rosina_end0','Hmel_start0','Hmel_end0'])
 with fp.open() as h:features=list(csv.DictReader(h,delimiter='\t'))
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
 views=[('01_Yb_current_figure_dotplot.pdf',1630000-1,1692000,'Yb: current figure span',10000),
  ('02_Ybv_SNP_strip_dotplot.pdf',1656500-1,1660500,'Yb(v): SNP-strip context',1000),
  ('03_Yb_broad_alignment_context.pdf',to,tend,'Yb: existing assembly-alignment context',100000)]
 outputs=[]
 for filename,lo,hi,title,step in views:
  lines=[(s,clipped(s,lo,hi)) for s in segments];lines=[(s,p) for s,p in lines if p is not None]
  if not lines:raise RuntimeError('No alignment segments in requested view: '+title)
  ys=[xy[1] for _,pair in lines for xy in pair];pad=max(100,(max(ys)-min(ys))*.04)
  yl=(min(ys)-pad,max(ys)+pad) if filename!=views[2][0] else (qo,qend)
  fig,ax=plt.subplots(figsize=(8.8,5.8));fig.subplots_adjust(left=.14,bottom=.17,right=.97,top=.86)
  for strand,color in [('+','#253746'),('-','#C44E52')]:
   data=[[(x/1e6,y/1e6) for x,y in pair] for s,pair in lines if s['strand']==strand]
   if data:ax.add_collection(LineCollection(data,colors=color,linewidths=.8,alpha=.8))
  for f in features:
   if f['feature'] not in ['Yb_d','Yb_v']:continue
   if f['start_status']!='exact' or f['end_status']!='exact':continue
   x0,x1=float(f['HmRos_start_exact'])-1,float(f['HmRos_end_exact'])
   y0,y1=float(f['Hmel2_5_start'])-1,float(f['Hmel2_5_end'])
   color='#9ECAE1' if f['feature']=='Yb_d' else '#F4CF58'
   if lo<x1 and hi>x0:
    ax.axvspan(x0/1e6,x1/1e6,color=color,alpha=.23,lw=0,zorder=0)
    ax.text((max(lo,x0)+min(hi,x1))/2e6,1.025,f['feature'].replace('_','(')+')',transform=ax.get_xaxis_transform(),ha='center',fontsize=8)
   if yl[0]<y1 and yl[1]>y0:ax.axhspan(y0/1e6,y1/1e6,color=color,alpha=.23,lw=0,zorder=0)
  ax.set(xlim=(lo/1e6,hi/1e6),ylim=(yl[0]/1e6,yl[1]/1e6),xlabel=f'Rosina · {tc} (Mb)',ylabel=f'Hmel2.5 melpomene · {qc} (Mb)',title=title)
  ax.xaxis.set_major_locator(MultipleLocator(step/1e6));ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
  ax.yaxis.set_major_formatter(FormatStrFormatter('%.3f'));ax.tick_params(length=3)
  ax.legend(handles=[Line2D([0],[0],color='#253746',label='Forward alignment'),Line2D([0],[0],color='#C44E52',label='Reverse alignment')],loc='lower right',frameon=False,fontsize=8)
  fig.text(.14,.045,'Existing minimap2 alignment; gaps are not joined. Assembly differences do not establish population variation.',fontsize=7,color='#555555')
  fig.savefig(out/filename);plt.close(fig)
  outputs.append(dict(PDF=str(out/filename),Rosina_interbase_span=[lo,hi],Hmel_interbase_span=list(yl),alignment_segments=len(lines)))
 save(out/'complete.json',dict(status='COMPLETE',outputs=outputs,alignment_records=len(blocks),
  source_job='138421491',new_downloads=False,new_alignments=False,figure_or_calls_changed=False,
  warning='M CIGAR runs include matches and mismatches. Blank areas may reflect alignment fragmentation/repeats/assembly differences, not proven structural variants.'))
 for o in outputs:print('PDF\t'+o['PDF'])
 print('STATUS\tCOMPLETE')

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--submit',action='store_true');p.add_argument('--run',action='store_true');p.add_argument('--output',type=Path,default=OUTPUT);p.add_argument('--paf',type=Path)
 a=p.parse_args()
 if a.submit==a.run:p.error('Choose exactly one of --submit or --run')
 if a.submit:prepare(a)
 else:run(a.output)
if __name__=='__main__':main()
