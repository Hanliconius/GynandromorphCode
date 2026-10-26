#!/usr/bin/env python3
"""Render NDJ Figure 4 from the frozen scan; run data-driven work on Pegasus."""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
SOURCE = ROOT/'analysis/NDJ_LR_difference_scan_v1/run_20260920T012415Z'
OUTPUT = ROOT/'analysis/Figure4_NDJ_panels'
CHROM = 'NC_140403.1'
BLUE, ORANGE = '#4C78A8', '#F58518'


def write_json(path, obj):
    tmp = path.with_suffix('.partial.json')
    tmp.write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n')
    tmp.replace(path)


def stamp(path):
    path = Path(path)
    return dict(path=str(path.resolve()), bytes=path.stat().st_size,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def submit(a):
    out = a.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if (out/'complete.json').exists() and not a.overwrite:
        raise RuntimeError('Completed output exists; use --overwrite to rerender in place')
    # Locate repo dependencies explicitly; never discover result inputs by filename.
    repo = Path(__file__).resolve().parents[3]
    deps = {
        'fig2_common.py': repo/'scripts/figures/figure2/common.py',
        'fig2_AF.py': repo/'scripts/figures/figure2/Fig2H_observed_AF.py',
        'Ybv_gene_features.json': repo/'scripts/yb/Ybv_gene_features.json',
    }
    for p in [a.source/'config.json', a.source/'COMPLETE.json', a.source/'windows.tsv',
              ROOT/'palettes/barfbow/luka_bleute.json', *deps.values()]:
        if not p.is_file():
            raise FileNotFoundError(p)
    logs = out/'logs'; logs.mkdir(exist_ok=True)
    shutil.copy2(__file__, out/'worker.py')
    for name, src in deps.items():
        shutil.copy2(src, out/name)
    run = dict(source=str(a.source.resolve()), output=str(out), ratio=a.ratio,
               zoom_start=a.zoom_start, zoom_end=a.zoom_end if a.zoom_end is not None else 3100000,
               anomaly_start=200000, anomaly_end=2100000, mixing_max=a.mixing_max,
               project=str(ROOT), code_sha256=stamp(out/'worker.py')['sha256'])
    if a.refresh_03_04:
        previous = json.loads((out/'run.json').read_text())
        for k in ['source','ratio','zoom_start','zoom_end','anomaly_start','anomaly_end','mixing_max','project']:
            run[k] = previous[k]
        if a.zoom_end is not None:run['zoom_end']=a.zoom_end
        run['refresh_03_04'] = True
        run['protected_outputs'] = [stamp(p) for p in sorted(out.glob('0[12]_*.pdf'))]
        if not run['protected_outputs']:
            raise RuntimeError('No existing 01/02 outputs found to preserve')
        if (out/'complete.json').exists():
            shutil.copy2(out/'complete.json',out/'previous_complete.json')
            (out/'complete.json').unlink()
    if a.full_chr15:
        cfg=json.loads((Path(run['source'])/'config.json').read_text())
        matches=[c for c in cfg['chromosomes'] if c['contig']==CHROM]
        if len(matches)!=1:raise RuntimeError('Expected exactly one chromosome 15 entry')
        run['zoom_start']=0;run['zoom_end']=int(matches[0]['length'])
        run['full_chr15']=True
    run['stack_C_Cprime']=a.stack
    write_json(out/'run.json', run)
    dependencies = []
    if a.refresh_03_04:
        # The frozen coverage cache is 100 kb; sub-window values cannot be inferred.
        # Reuse a validated cached prefix, counting only any newly requested extent.
        cfg = json.loads((Path(run['source'])/'config.json').read_text())
        required = [Path(cfg[k]) for k in ['bam_left','bam_right','reference']]
        required += [Path(run['source'])/'background.json']
        for path in required:
            if not path.is_file(): raise FileNotFoundError(path)
        for side in ['left','right']:
            bam=Path(cfg['bam_'+side])
            if not any(p.is_file() for p in [Path(str(bam)+'.bai'),bam.with_suffix('.bai'),Path(str(bam)+'.csi')]):
                raise FileNotFoundError('Missing BAM index: '+str(bam))
        write_json(out/'refinement_reuse_plan.json',dict(
            samples=['NDJ_left','NDJ_right'],source_config=stamp(Path(run['source'])/'config.json'),
            action='recount_existing_alignment',reason='100-kb source windows cannot resolve 10-kb depth',
            interval=[CHROM,run['zoom_start'],run['zoom_end']],window_bp=10000,
            BAMs=[str(cfg['bam_left']),str(cfg['bam_right'])],downloads=0,alignments=0,calling=0))
        depth_cmd=[sys.executable,str(out/'worker.py'),'--depth-worker','--output',str(out)]
        depth_job=subprocess.check_output(['sbatch','--parsable','--job-name=NDJ_10kb',
            '--partition=nano','--time=00:30:00','--cpus-per-task=1','--mem=4G','--array=0-1',
            '--output='+str(logs/'depth.%A_%a.out'),'--error='+str(logs/'depth.%A_%a.err'),
            '--wrap='+shlex.join(depth_cmd)],text=True).strip().split(';')[0]
        dependencies=['--dependency=afterok:'+depth_job]
        print('CHR15_DEPTH_ARRAY\t'+depth_job)
    cmd = [sys.executable, str(out/'worker.py'), '--run', '--output', str(out)]
    job = subprocess.check_output([
        'sbatch', '--parsable', '--job-name=NDJ_Fig4', '--partition=nano',
        '--time=00:30:00', '--cpus-per-task=1', '--mem=12G',
        '--output='+str(logs/'plot.%j.out'), '--error='+str(logs/'plot.%j.err'),
        *dependencies, '--wrap='+shlex.join(cmd)], text=True).strip().split(';')[0]
    write_json(out/'jobs.json', dict(render=job, **({'depth_array':depth_job} if a.refresh_03_04 else {})))
    print('JOB\t'+job+'\nOUTPUT\t'+str(out))


def depth_worker(out):
    import numpy as np
    import pandas as pd
    import pysam
    run=json.loads((out/'run.json').read_text())
    cfg=json.loads((Path(run['source'])/'config.json').read_text())
    side=['left','right'][int(os.environ['SLURM_ARRAY_TASK_ID'])]
    bam=Path(cfg['bam_'+side]); ref=Path(cfg['reference'])
    def meta(p):
        p=Path(p);st=p.stat()
        return dict(path=str(p.resolve()),bytes=st.st_size,mtime_ns=st.st_mtime_ns)
    signature=dict(bam=meta(bam),reference=meta(ref),chrom=CHROM,start=run['zoom_start'],
        end=run['zoom_end'],window_bp=10000,baseq=cfg['baseq'],mapq=cfg['mapq'],exclude='0xF04',
        mate_overlap_removed=True,samtools=cfg['samtools'])
    table=out/f'chr15_10kb_{side}.tsv';done=out/f'chr15_10kb_{side}.complete.json'
    prefix=None;prefix_audit=None
    if table.exists() and done.exists():
        old=json.loads(done.read_text())
        if old.get('signature')==signature and old.get('table')==stamp(table):
            print('REUSED\t'+str(table));return
        oldsig=old.get('signature',{})
        if (old.get('table')==stamp(table) and oldsig.get('end',0)<signature['end']
            and oldsig.get('end',0)%10000==0
            and {k:v for k,v in oldsig.items() if k!='end'}=={k:v for k,v in signature.items() if k!='end'}):
            prefix=pd.read_csv(table,sep='\t');prefix_audit=old
            signature_start=int(oldsig['end'])
    done.unlink(missing_ok=True)
    start,end,width=run['zoom_start'],run['zoom_end'],10000
    if prefix is not None:start=signature_start
    with pysam.FastaFile(str(ref)) as fa, pysam.AlignmentFile(str(bam),'rb') as b:
        if b.get_reference_length(CHROM)!=fa.get_reference_length(CHROM):
            raise RuntimeError('BAM/reference length mismatch')
        seq=fa.fetch(CHROM,start,end).upper()
    if len(seq)!=end-start:raise RuntimeError('Reference interval truncated')
    n=(end-start+width-1)//width
    totals=np.zeros(n,dtype=np.int64);covered=np.zeros(n,dtype=np.int64)
    command=[cfg['samtools'],'depth','-aa','-s','-q',str(cfg['baseq']),'-Q',str(cfg['mapq']),
        '-G','0xF04','-r',f'{CHROM}:{start+1}-{end}',str(bam)]
    previous=start
    with (out/f'chr15_10kb_{side}.stderr.log').open('w') as err:
        proc=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=err,text=True)
        try:
            for line in proc.stdout:
                chrom,pos,d=line.rstrip().split('\t');pos=int(pos);d=int(d)
                if chrom!=CHROM or pos!=previous+1:raise RuntimeError('Non-consecutive depth output')
                previous=pos;idx=(pos-start-1)//width;totals[idx]+=d;covered[idx]+=d>0
            if proc.wait()!=0:raise RuntimeError('samtools depth failed')
        except BaseException:
            proc.kill();proc.wait();raise
    if previous!=end:raise RuntimeError('Incomplete depth stream')
    rows=[]
    for i in range(n):
        a=start+i*width;b=min(a+width,end);s=seq[a-start:b-start]
        rows.append(dict(chrom=CHROM,start=a,end=b,depth=totals[i]/(b-a),
            covered=covered[i]/(b-a),acgt_fraction=sum(s.count(x) for x in 'ACGT')/(b-a)))
    data=pd.DataFrame(rows)
    if prefix is not None:data=pd.concat([prefix,data],ignore_index=True)
    if data.duplicated(['chrom','start','end']).any():raise RuntimeError('Overlapping cache extension')
    temp=table.with_suffix('.partial.tsv');data.to_csv(temp,sep='\t',index=False);temp.replace(table)
    write_json(done,dict(signature=signature,command=command,table=stamp(table),
        reused_prefix=prefix_audit,new_interval=[start,end]))
    print('DEPTH_COMPLETE\t'+str(table))


def dense_chr15(out,run,normalization):
    import pandas as pd
    import numpy as np
    bg=json.loads((Path(run['source'])/'background.json').read_text())['depth_normalizers']
    for side in ['left','right']:
        done=json.loads((out/f'chr15_10kb_{side}.complete.json').read_text())
        if done['table']!=stamp(out/f'chr15_10kb_{side}.tsv'):
            raise RuntimeError('Dense coverage table changed after completion: '+side)
    left=pd.read_csv(out/'chr15_10kb_left.tsv',sep='\t')
    right=pd.read_csv(out/'chr15_10kb_right.tsv',sep='\t')
    d=left.merge(right,on=['chrom','start','end'],suffixes=('_left','_right'),validate='one_to_one')
    use=(d.covered_left.ge(.5)&d.covered_right.ge(.5)&d.acgt_fraction_left.ge(.9)&d.acgt_fraction_right.ge(.9))
    d=d.loc[use].copy();d['mid']=(d.start+d.end)/2
    for side in ['left','right']:
        d[side]=d['depth_'+side]/bg[side]/normalization[side]
    d['delta']=d.right-d.left
    fit=fit_segments(d)
    d.attrs['segment_fit']=fit
    for side in ['left','right']:
        d[side+'_segment_median']=np.where(d.mid<fit['split_bp'],fit['levels'][side][0],fit['levels'][side][1])
    write_json(out/'chr15_descriptive_segments.json',fit)
    d.to_csv(out/'chr15_10kb_plotted.tsv',sep='\t',index=False)
    return d


def fit_segments(d):
    """One shared descriptive transition; robust median levels, no copy-number fit."""
    import numpy as np
    if not d.mid.is_monotonic_increasing:raise RuntimeError('Unsorted depth windows')
    values=d[['left','right']].to_numpy(float)
    null=float(np.abs(values-np.median(values,axis=0)).sum())
    choices=[]
    for i in range(20,len(d)-19):
        if d.mid.iloc[i-1]-d.mid.iloc[0]<500000 or d.mid.iloc[-1]-d.mid.iloc[i]<500000:continue
        a,b=values[:i],values[i:]
        ma,mb=np.median(a,axis=0),np.median(b,axis=0)
        loss=float(np.abs(a-ma).sum()+np.abs(b-mb).sum())
        choices.append((loss,i,ma,mb))
    if not choices:raise RuntimeError('Insufficient usable coverage for descriptive segments')
    loss,i,a,b=min(choices,key=lambda item:item[0])
    split=float((d.mid.iloc[i-1]+d.mid.iloc[i])/2)
    return dict(method='Shared one-transition median segmentation; minimum total absolute deviation across both sides',
        interpretation='Descriptive summary only; not a ploidy call, significance test, or base-resolution breakpoint',
        split_bp=split,transition_bracket_bp=[int(d.start.iloc[i-1]),int(d.end.iloc[i])],
        window_bp=10000,minimum_segment_span_bp=500000,minimum_usable_windows_per_segment=20,
        levels={s:[float(a[j]),float(b[j])] for j,s in enumerate(['left','right'])},
        absolute_loss=loss,one_segment_absolute_loss=null,
        start_bp=int(d.start.iloc[0]),end_bp=int(d.end.iloc[-1]),
        left_colour='#5A4FA3',right_colour='#E56B00')


def bools(s):
    if s.dtype == bool:
        return s
    t = s.astype(str).str.lower()
    if not t.isin(['true', 'false', '1', '0']).all():
        raise ValueError('Unrecognized boolean field')
    return t.isin(['true', '1'])


def load_data(run):
    import numpy as np
    import pandas as pd
    src = Path(run['source'])
    cfg = json.loads((src/'config.json').read_text())
    done = json.loads((src/'COMPLETE.json').read_text())
    if done.get('status') != 'exploratory_scan_complete_not_biological_validation':
        raise RuntimeError('Source NDJ scan is not marked complete')
    identities = [str(cfg[k]).lower() for k in ('left', 'right')]
    if identities != ['ndj_left', 'ndj_right']:
        raise RuntimeError('Wrong specimen/sample order: '+repr(identities))
    w = pd.read_csv(src/'windows.tsv', sep='\t')
    required = {'chrom','start','end','usable_depth','contig_background',
                'left_depth_relative','right_depth_relative'}
    if required - set(w.columns):
        raise RuntimeError('Missing coverage fields: '+repr(required-set(w.columns)))
    if w.duplicated(['chrom','start','end']).any():
        raise RuntimeError('Duplicate coverage windows')
    chr15 = [(i,c) for i,c in enumerate(cfg['chromosomes']) if c['contig'] == CHROM]
    if len(chr15) != 1:
        raise RuntimeError('Exactly one Rosina chromosome 15 entry required')
    length = int(chr15[0][1]['length'])
    if not 0 <= run['zoom_start'] < run['anomaly_start'] < run['anomaly_end'] < run['zoom_end'] <= length:
        raise RuntimeError('Zoom must include both sides of the earlier 0.2-2.1 Mb interval')
    use = bools(w.usable_depth)
    for side in ['left','right']:
        use &= np.isfinite(w[side+'_depth_relative']) & w[side+'_depth_relative'].gt(0)
    w = w.loc[use].copy()
    # Re-estimate only a scalar normalization from already cached windows.
    # Exclude chr15 and sex chromosomes from the autosomal background.
    background = bools(w.contig_background) & w.chrom.ne(CHROM)
    if background.sum() < 20:
        raise RuntimeError('Insufficient autosomal background windows')
    normalization = {}
    for side in ['left','right']:
        scale = float(w.loc[background, side+'_depth_relative'].median())
        normalization[side] = scale
        w[side] = w[side+'_depth_relative']/scale
    w['delta'] = w.right-w.left
    w['mid'] = (w.start+w.end)/2
    bounds = dict(delta=list(map(float, np.quantile(w.loc[background,'delta'], [.025,.975]))))
    # Union of the side-specific reference envelopes, not a ploidy confidence band.
    q = [np.quantile(w.loc[background,s], [.025,.975]) for s in ['left','right']]
    bounds['relative'] = [float(min(x[0] for x in q)), float(max(x[1] for x in q))]
    return cfg, w, normalization, bounds


def style():
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':8,
        'axes.labelsize':8, 'xtick.labelsize':7, 'ytick.labelsize':7,
        'axes.linewidth':.65, 'pdf.fonttype':42, 'ps.fonttype':42,
        'savefig.facecolor':'white'})


def clean(ax):
    for s in ['top','right']:
        ax.spines[s].set_visible(False)
    ax.tick_params(length=2.5, width=.65, pad=2)


def save(fig, path, sizes):
    # Fixed canvas, including labels: deliberately do NOT use bbox_inches='tight'.
    import matplotlib.pyplot as plt
    w,h = map(float,fig.get_size_inches())
    fig.savefig(path)
    sizes[path.name] = dict(width_inches=w,height_inches=h,outer_ratio=w/h)
    plt.close(fig)


def expectations(ax, metric, bounds, mix, xmax):
    """Model ranges and empirical background are distinct in label and provenance."""
    if metric == 'delta':
        regions = [(1-2*mix,1,'L 1n / R 3n'),(-1,-1+2*mix,'L 3n / R 1n')]
        base = 0
    else:
        regions = [(.5,.5+mix,'1n'),(1.5-mix,1.5,'3n')]
        base = 1
    for low,high,label in regions:
        ax.axhspan(low,high,facecolor='#D8D8D8',alpha=.55,zorder=1)
        ax.axhline(low,color='#999999',lw=.4,zorder=1)
        ax.axhline(high,color='#999999',lw=.4,zorder=1)
        ax.text(xmax,(low+high)/2,label,fontsize=6.5,ha='right',va='center',
                color='#555555',bbox=dict(facecolor='white',alpha=.85,edgecolor='none',pad=.7),zorder=6)
    lo,hi = bounds[metric]
    ax.axhspan(lo,hi,facecolor='#BFCAD2',alpha=.32,zorder=0)
    ax.axhline(base,color='#444444',lw=.6,zorder=2)


def coverage(out, run, cfg, w, bounds, features, sizes, dense=None, stacked=False):
    import matplotlib.pyplot as plt
    import numpy as np
    offset,ticks,labels = {},[],[]
    cursor = 0
    for c in cfg['chromosomes']:
        offset[c['contig']] = cursor
        ticks.append((cursor+c['length']/2)/1e6)
        labels.append('Z' if str(c['label'])=='21' else str(c['label']))
        cursor += c['length']+300000
    w['genome_mb'] = (w.chrom.map(offset)+w.mid)/1e6
    if stacked:
        canvas=plt.figure(figsize=(10,20/run['ratio']))
        panels=canvas.subfigures(2,1,hspace=0)
    for zoom in [False,True]:
        if run.get('refresh_03_04') and not zoom and not stacked:continue
        for metric in ['delta','relative']:
            if stacked and metric!=('relative' if zoom else 'delta'):continue
            fig = panels[int(zoom)] if stacked else plt.figure(figsize=(10,10/run['ratio']))
            ax = fig.add_axes([.075,.29,.91,.56] if zoom else [.075,.25,.91,.66])
            fig.text(.013,.98,"C'" if zoom else 'C',fontsize=13,fontweight='bold',va='top')
            if zoom:
                data = w.loc[w.chrom.eq(CHROM) & w.end.gt(run['zoom_start']) & w.start.lt(run['zoom_end'])]
                if dense is not None:data=dense
                lo,hi = run['zoom_start']/1e6,run['zoom_end']/1e6
                x = data.mid/1e6
                ax.set_xlabel('Chromosome 15 position (Mb)',labelpad=3)
                # Earlier exploratory boundaries; not newly estimated breakpoints.
                for p in [run['anomaly_start'],run['anomaly_end']]:
                    ax.axvline(p/1e6,color='#777777',lw=.65,ls=':',zorder=2)
                exons = features['ivory']['exons']
                a,b = min(x[0] for x in exons)/1e6,max(x[1] for x in exons)/1e6
                ax.plot([a,b],[1.06,1.06],transform=ax.get_xaxis_transform(),c='#222222',lw=1,clip_on=False)
                ax.plot([(a+b)/2], [1.06],marker='v',ms=4,c='#222222',transform=ax.get_xaxis_transform(),clip_on=False)
                ax.text((a+b)/2,1.11,'ivory*',transform=ax.get_xaxis_transform(),ha='center',fontsize=7)
                suffix = 'chr15'
            else:
                data = w; x = w.genome_mb
                lo,hi = 0,(cursor-300000)/1e6
                for i,c in enumerate(cfg['chromosomes']):
                    if i%2:
                        ax.axvspan(offset[c['contig']]/1e6,(offset[c['contig']]+c['length'])/1e6,
                                   color='#F0F0F0',zorder=-2)
                ax.set_xticks(ticks,labels)
                ax.set_xlabel('Chromosome',labelpad=3)
                suffix = 'genome'
            expectations(ax,metric,bounds,run['mixing_max'],hi-.005*(hi-lo))
            point_slice=slice(None,None,3) if zoom and dense is not None else slice(None)
            if metric == 'delta':
                ax.scatter(x.iloc[point_slice],data.delta.iloc[point_slice],s=7 if zoom else 4,c='black',linewidths=0,zorder=3)
                if not zoom:
                    selected = data.chrom.eq(CHROM) & data.mid.between(run['anomaly_start'],run['anomaly_end'])
                    ax.scatter(x.loc[selected],data.loc[selected,'delta'],s=7,c='#C92A2A',linewidths=0,zorder=4)
                ax.set_ylabel('Relative coverage\nRight - Left')
                vals = data.delta.to_numpy()
                # Retain every point; no quantile clipping or hidden extreme windows.
                ax.set_ylim(min(-1.15,float(vals.min())-.08,bounds['delta'][0]-.05),
                            max(1.15,float(vals.max())+.08,bounds['delta'][1]+.05))
            else:
                point_colours=[('left','#5A4FA3'),('right','#E56B00')] if zoom and dense is not None else [('left',BLUE),('right',ORANGE)]
                for side,col in point_colours:
                    fine=zoom and dense is not None
                    ax.scatter(x.iloc[point_slice],data[side].iloc[point_slice],s=4 if fine else (9 if zoom else 4),c=col,linewidths=0,
                               alpha=.85 if fine else .8,zorder=3,label=side.title())
                    if fine:
                        fit=data.attrs['segment_fit'];level=fit['levels'][side]
                        for j,(a,b) in enumerate([(fit['start_bp'],fit['split_bp']),(fit['split_bp'],fit['end_bp'])]):
                            ax.hlines(level[j],a/1e6,b/1e6,colors='#222222' if side=='left' else '#666666',lw=1.5,zorder=5)
                vals = data[['left','right']].to_numpy()
                ax.set_ylim(max(0,min(.35,float(vals.min())-.08)),max(1.65,float(vals.max())+.08,bounds['relative'][1]+.05))
                ax.set_ylabel('Relative coverage\nDiploid background = 1')
                ax.legend(frameon=False,loc='upper left',ncol=2,fontsize=6.5,handletextpad=.3)
            ax.set_xlim(lo,hi); clean(ax)
            fig.text(.075,.045,f'Blue-grey: central 95% autosomal background; grey: conditional 0-{run["mixing_max"]:.0%} mixing model.',fontsize=6.3)
            if zoom:
                note=('Every third 10-kb point shown; segment medians use all windows (black/grey).' if metric=='relative' else 'Every third 10-kb coverage point shown.') if dense is not None else 'Cached source depth windows.'
                fig.text(.075,.085,note,fontsize=6.1)
                fig.text(.075,.012,'Dotted: earlier exploratory interval. *Mapped ivory exons; annotation is partial.',fontsize=6.1)
            if not stacked:save(fig,out/f'{"03" if zoom else "02"}_{suffix}_{metric}.pdf',sizes)
    if stacked:save(canvas,out/'05_C_Cprime_stacked.pdf',sizes)


def bars(out,cfg,w,sizes):
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.ticker import FuncFormatter
    rows = []
    for c in cfg['chromosomes']:
        d = w.loc[w.chrom.eq(c['contig'])]
        if d.empty:
            continue
        weights = d.end-d.start
        rows.append(dict(label='Z' if str(c['label'])=='21' else str(c['label']),
                         left=float(np.average(d.left,weights=weights)),
                         right=float(np.average(d.right,weights=weights)),usable_windows=len(d)))
    import pandas as pd
    pd.DataFrame(rows).to_csv(out/'chromosome_coverage.tsv',sep='\t',index=False)
    fig = plt.figure(figsize=(4.25,6.15)); ax = fig.add_axes([.10,.13,.85,.80])
    fig.text(.018,.988,'B',fontsize=15,fontweight='bold',va='top')
    y = np.arange(len(rows)); mx = max(1.1,max(max(r['left'],r['right']) for r in rows)*1.1)
    for side,col,sign in [('left',BLUE,-1),('right',ORANGE,1)]:
        ax.barh(y,[sign*r[side] for r in rows],height=.48,color=col,edgecolor='black',linewidth=.35)
    for p in [-1,1]:
        ax.axvline(p,color='#999999',ls='--',lw=.8,zorder=0)
    ax.axvspan(-.105*mx,.105*mx,color='white',zorder=4)
    ax.axvline(0,color='#444444',lw=.8,zorder=5)
    for i,r in enumerate(rows):
        ax.text(0,i,r['label'],ha='center',va='center',fontsize=7.4,zorder=6)
    ax.set(xlim=(-mx,mx),ylim=(len(rows)-.4,-.6),yticks=[],xlabel='Relative coverage')
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,p:f'{abs(x):.1f}'))
    ax.grid(axis='x',color='#E5E5E5',lw=.65);ax.set_axisbelow(True);clean(ax)
    ax.spines['left'].set_visible(False)
    fig.text(.25,.065,'Left',ha='center',color=BLUE)
    fig.text(.78,.065,'Right',ha='center',color=ORANGE)
    fig.text(.5,.95,'Chromosome',ha='center')
    fig.text(.5,.015,'Usable-window means; W not assessed by source scan.',ha='center',fontsize=6.2)
    save(fig,out/'01_B_chromosome_coverage.pdf',sizes)


def af(out,run,cfg,sizes):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap,LogNorm
    from matplotlib.collections import PolyCollection
    from matplotlib.backends.backend_pdf import PdfPages
    # Reuse Figure 2 functions; neutral marginals because no original-GT classification
    # for NDJ is present in this cache. Do not substitute TN genotype histograms.
    common = module(out/'fig2_common.py','common');sys.modules['common']=common
    aesthetic = module(out/'fig2_AF.py','ndj_fig2_af')
    palette = Path(run['project'])/'palettes/barfbow/luka_bleute.json'
    if not palette.is_file():
        raise FileNotFoundError('Exact Figure 2 palette required: '+str(palette))
    cmap = ListedColormap(json.loads(palette.read_text()))
    nx,ny = 45,int(45/np.sqrt(3));sx,sy=(1+2e-9)/nx,1/ny
    centres = np.concatenate([np.array([(i,j) for i in range(nx+1) for j in range(ny+1)]),
                              np.array([(i+.5,j+.5) for i in range(nx) for j in range(ny)])])*[sx,sy]
    centres[:,0]-=1e-9
    def hexids(l,r):
        x,y=(l+1e-9)/sx,r*ny; ix,iy=np.rint(x).astype(int),np.rint(y).astype(int)
        jx,jy=np.floor(x).astype(int),np.floor(y).astype(int)
        first=(x-ix)**2+3*(y-iy)**2<(x-jx-.5)**2+3*(y-jy-.5)**2
        return np.where(first,ix*(ny+1)+iy,(nx+1)*(ny+1)+jx*ny+jy)
    stats={k:dict(hex=np.zeros(len(centres),dtype=np.int64),lh=np.zeros(50,dtype=np.int64),rh=np.zeros(50,dtype=np.int64)) for k in ['genome','chr15_region']}
    inputs=[]
    for i,c in enumerate(cfg['chromosomes']):
        path=Path(run['source'])/'chromosomes'/f'{i:02d}'/'counts.npz'
        inputs.append(stamp(path))
        with np.load(path,allow_pickle=False) as z:
            pos=z['pos'];nl=z['left_ref']+z['left_alt'];nr=z['right_ref']+z['right_alt']
            if np.any(nl<cfg['min_depth']) or np.any(nr<cfg['min_depth']):
                raise RuntimeError('Counts violate recorded depth threshold')
            l=z['left_alt']/nl;r=z['right_alt']/nr
            for name,mask in [('genome',np.ones(len(pos),dtype=bool)),
                              ('chr15_region',(pos>run['anomaly_start'])&(pos<=run['anomaly_end']))]:
                if name=='chr15_region' and c['contig']!=CHROM:continue
                stats[name]['hex']+=np.bincount(hexids(l[mask],r[mask]),minlength=len(centres))
                stats[name]['lh']+=np.histogram(l[mask],bins=50,range=(0,1))[0]
                stats[name]['rh']+=np.histogram(r[mask],bins=50,range=(0,1))[0]
    audit=[]
    with PdfPages(out/'04_AF_Figure2_style.pdf') as pdf:
        for name in stats:
            d=stats[name];n=int(d['hex'].sum())
            if n==0 or d['lh'].sum()!=n or d['rh'].sum()!=n:raise RuntimeError('AF count total mismatch')
            counts=d['hex'];keep=counts>0
            compressed=common.informative_density_vmax(counts,centres[:,0],centres[:,1])
            vmax=float(np.sqrt(compressed*counts.max()))
            fig=plt.figure(figsize=(5.3,5.3));ax=fig.add_axes([.25,.29,.60,.60])
            top=fig.add_axes([.25,.925,.60,.035]);side=fig.add_axes([.09,.29,.045,.60]);cbax=fig.add_axes([.89,.43,.02,.29])
            unit=np.array([[0,-.5],[.5,-.25],[.5,.25],[0,.5],[-.5,.25],[-.5,-.25]])
            # Voronoi metric dx^2+3dy^2: width sx, FULL height 2*sy/3.
            # The earlier factor 4*sy/3 doubled the height and overlapped bins.
            verts=centres[keep,None,:]+unit[None,:,:]*[sx,2*sy/3]
            coll=PolyCollection(verts,array=counts[keep].astype(float),cmap=cmap,
                                norm=LogNorm(1,max(1.0001,vmax),clip=True),edgecolors='none')
            ax.add_collection(coll);ax.set(xlim=(-.025,1.025),ylim=(-.025,1.025),aspect='equal',
                xlabel='Left ALT allele fraction',ylabel='Right ALT allele fraction')
            ticks=[0,.25,.5,.75,1];ax.set_xticks(ticks);ax.set_yticks(ticks);clean(ax)
            ax.xaxis.label.set_size(8);ax.yaxis.label.set_size(8)
            ax.tick_params(labelsize=7)
            ax.plot([0,1],[0,1],color='white',lw=.5,alpha=.45)
            x=(np.arange(50)+.5)/50
            aesthetic.marginal(top,d['lh'],x,(-.025,1.025))
            aesthetic.marginal(side,d['rh'],x,(-.025,1.025),vertical=True)
            # Keep endpoint arrows, move their percentages into the footer.
            for marginal_ax in [top,side]:
                for txt in list(marginal_ax.texts):txt.remove()
            cb=fig.colorbar(coll,cax=cbax,extend='max');cb.set_label('Sites per hexagon',fontsize=7)
            cb.ax.tick_params(labelsize=6);cb.outline.set_linewidth(.5)
            title='Whole genome' if name=='genome' else 'Chromosome 15: 0.2-2.1 Mb'
            fig.text(.25,.985,title,fontsize=9,va='top')
            fig.text(.25,.19,f'{n:,} SNP records; AD depth >= {cfg["min_depth"]} on both sides',fontsize=7)
            fig.text(.25,.15,f'AF=0 / 1 bins: Left {d["lh"][0]/n:.1%} / {d["lh"][-1]/n:.1%}; Right {d["rh"][0]/n:.1%} / {d["rh"][-1]/n:.1%}',fontsize=6.2)
            fig.text(.25,.11,'Observed counts; no mixing correction or fitted ploidy overlay.',fontsize=6.3)
            fig.text(.25,.075,'Control-discovered SNP ascertainment retained.',fontsize=6.3)
            pdf.savefig(fig);plt.close(fig)
            audit.append(dict(scope=name,SNP_records=n,display_vmax=vmax))
            pd.DataFrame(dict(left_centre=centres[:,0],right_centre=centres[:,1],count=counts)).to_csv(out/(name+'_AF_hexagons.tsv'),sep='\t',index=False)
    sizes['04_AF_Figure2_style.pdf']=dict(width_inches=5.3,height_inches=5.3,pages=2)
    return audit,inputs,stamp(palette)


def run_render(out):
    import pandas as pd
    run=json.loads((out/'run.json').read_text())
    (out/'complete.json').unlink(missing_ok=True)
    style();cfg,w,normalization,bounds=load_data(run)
    # Validate all required AF inputs before producing any figures.
    for i,c in enumerate(cfg['chromosomes']):
        path=Path(run['source'])/'chromosomes'/f'{i:02d}'/'counts.npz'
        if not path.is_file() or not path.stat().st_size:
            raise FileNotFoundError(path)
    features=json.loads((out/'Ybv_gene_features.json').read_text())
    if features['assembly']!='GCF_056956725.1' or features['contig']!=CHROM:
        raise RuntimeError('Wrong ivory reference')
    decisions=[dict(biological_sample_id=s,interval='whole scan / chr15',reference=cfg['reference'],action='reuse_calls_and_coverage',selected_resource_id=str(Path(run['source'])/'windows.tsv'),evidence='completed NDJ scan; cached counts and depth',unresolved_questions='W not in source scan') for s in ['NDJ_left','NDJ_right']]
    pd.DataFrame(decisions).to_csv(out/'reuse_decisions.tsv',sep='\t',index=False)
    refresh=run.get('refresh_03_04',False)
    prior=json.loads((out/'previous_complete.json').read_text()) if refresh else {}
    sizes=dict(prior.get('PDFs',{}))
    dense=dense_chr15(out,run,normalization) if refresh else None
    if not refresh:bars(out,cfg,w,sizes)
    coverage(out,run,cfg,w,bounds,features,sizes,dense)
    if run.get('stack_C_Cprime'):
        coverage(out,run,cfg,w,bounds,features,sizes,dense,stacked=True)
    if refresh and run.get('stack_C_Cprime'):
        af_audit,count_inputs,palette=prior['AF'],prior['count_inputs'],prior['palette']
    else:
        af_audit,count_inputs,palette=af(out,run,cfg,sizes)
    if not refresh:w.to_csv(out/'coverage_windows_plotted.tsv.gz',sep='\t',index=False)
    used=[stamp(Path(run['source'])/n) for n in ['config.json','COMPLETE.json','windows.tsv']]
    audit=dict(status='COMPLETE',specimen=['NDJ_left','NDJ_right'],source=run['source'],
        inputs=used,count_inputs=count_inputs,palette=palette,run=run,PDFs=sizes,
        coverage_window_bp=cfg['window'],source_filters={k:cfg.get(k) for k in ['mapq','baseq','min_depth','max_depth','filter','ascertainment']},
        coverage_normalization='Existing side-specific relative depth divided by its median in usable autosomal windows excluding chr15 and sex chromosomes',
        normalization_scalars=normalization,background_envelopes=bounds,
        background_interpretation='Empirical central 95% observed autosomal variation, not confidence interval or calibrated null',
        mixing_model='c is opposite-side DNA/read contribution measured at a diploid background; D_L=(1-c_L)*p_L/2+c_L*p_R/2; D_R=c_R*p_L/2+(1-c_R)*p_R/2',
        diploid_expectation=1,one_n_band=[.5,.5+run['mixing_max']],three_n_band=[1.5-run['mixing_max'],1.5],
        reciprocal_delta_band=[1-2*run['mixing_max'],1],
        model_interpretation='Conditional reciprocal 1n/3n versus 2n background, not independent evidence of ploidy; mixing is not estimated here',
        W='not assessed by source scan; not shown as zero',ivory=features,AF=af_audit,
        AF_aesthetic='Figure2H Luka-Bleute palette, vector hexagons, log density midpoint compression, neutral LOESS marginals; no unsupported control-GT colours or TN predictive rings',
        downloads=0,alignments=0,calling=0,new_coverage=2 if refresh else 0,
        chr15_refinement=dict(window_bp=10000,source='same NDJ BAMs and samtools filters',
            style='Purple Left/orange Right raw points; black/grey descriptive median segments',
            display_stride=3,segment_fit_uses_all_usable_windows=True,
            segment_fit=dense.attrs['segment_fit']) if refresh else None)
    for old in run.get('protected_outputs',[]):
        if stamp(Path(old['path']))!=old:raise RuntimeError('Protected 01/02 PDF changed: '+old['path'])
    audit['protected_outputs_unchanged']=run.get('protected_outputs',[])
    write_json(out/'complete.json',audit)
    print('STATUS\tCOMPLETE')
    for name in sizes:print('PDF\t'+str(out/name))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=SOURCE)
    p.add_argument('--output',type=Path,default=OUTPUT)
    p.add_argument('--submit',action='store_true');p.add_argument('--run',action='store_true')
    p.add_argument('--overwrite',action='store_true')
    p.add_argument('--refresh-03-04',action='store_true',help='Recount chr15 at 10 kb; update only 03/04 and preserve 01/02 hashes')
    p.add_argument('--depth-worker',action='store_true')
    p.add_argument('--stack',action='store_true',help='Stack genome delta above chr15 relative coverage; keep existing AF PDF on targeted refresh')
    p.add_argument('--full-chr15',action='store_true',help='Display all of reference chromosome 15, resolving its exact length from the frozen source config')
    p.add_argument('--ratio',type=float,choices=[4,5],default=4)
    p.add_argument('--zoom-start',type=int,default=0);p.add_argument('--zoom-end',type=int,default=None)
    p.add_argument('--mixing-max',type=float,default=.2)
    a=p.parse_args()
    if not 0<=a.mixing_max<.5:p.error('Mixing must lie in [0,0.5)')
    if a.submit:submit(a)
    elif a.depth_worker:depth_worker(a.output.resolve())
    elif a.run:run_render(a.output.resolve())
    else:p.error('Choose --submit (Pegasus) or internal --run')


if __name__=='__main__':main()
