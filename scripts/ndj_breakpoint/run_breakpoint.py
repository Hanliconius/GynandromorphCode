#!/usr/bin/env python3
"""NDJ chr15 breakpoint context: audit, sequence screens, backgrounds and BAM evidence."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import random
import shutil
import shlex
import subprocess
import sys
from datetime import datetime, timezone
from sequence_features import Track, scan, score

PROJECT = Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
CHROM = 'NC_140403.1'
ASSEMBLY = 'GCF_056956725.1'


def js(path, data):
    tmp = Path(str(path)+'.partial')
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
    tmp.replace(path)


def table(path, rows, fields=None):
    rows = list(rows)
    fields = fields or (list(rows[0]) if rows else [])
    with Path(path).open('w') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)


def stamp(path, checksum=False):
    p = Path(path).resolve(strict=True); st = p.stat()
    d = dict(path=str(p), bytes=st.st_size, mtime_ns=st.st_mtime_ns)
    if checksum:
        h = hashlib.sha256()
        with p.open('rb') as f:
            for block in iter(lambda:f.read(1048576), b''): h.update(block)
        d['sha256'] = h.hexdigest()
    return d


def read_json(path):
    return json.loads(Path(path).read_text())


def index_path(path):
    return next((p for p in [Path(str(path)+'.bai'), Path(path).with_suffix('.bai'),
                             Path(str(path)+'.csi')] if p.is_file()), None)


def resolve(a):
    import pysam
    cfg = read_json(a.source/'config.json')
    if [str(cfg.get(k, '')).lower() for k in ['left','right']] != ['ndj_left','ndj_right']:
        raise ValueError('Source config must identify NDJ_left and NDJ_right')
    complete = read_json(a.source/'COMPLETE.json')
    if complete.get('status') != 'exploratory_scan_complete_not_biological_validation':
        raise ValueError('Source scan completion marker not recognized')
    ref = Path(cfg['reference'])
    if not Path(str(ref)+'.fai').is_file(): raise FileNotFoundError(str(ref)+'.fai')
    match = [c for c in cfg['chromosomes'] if c['contig'] == CHROM]
    if len(match) != 1: raise ValueError('Expected exactly one chr15 entry')
    with pysam.FastaFile(str(ref)) as fa:
        length = fa.get_reference_length(CHROM)
    if length != int(match[0]['length']): raise ValueError('Source/reference length mismatch')
    # Never substitute historical 0.2/2.1-Mb plotting guide lines for a fitted transition.
    if a.interval:
        start, end = a.interval
        origin = 'explicit user-supplied 0-based half-open uncertainty interval'
        segment_inputs = []
    else:
        fitpath = a.figure4/'chr15_descriptive_segments.json'
        fit = read_json(fitpath)
        run = read_json(a.figure4/'run.json')
        if Path(run['source']).resolve() != a.source.resolve():
            raise ValueError('Figure4 and requested source runs differ')
        start, end = map(int, fit['transition_bracket_bp'])
        origin = 'Figure4 descriptive transition bracket; NOT a confidence interval or resolved junction'
        segment_inputs = [stamp(fitpath, True), stamp(a.figure4/'run.json', True)]
    if not 0 <= start < end <= length: raise ValueError('Invalid interval')
    if end-start > 1000000: raise ValueError('Interval >1 Mb: resolve the coverage transition first')
    bams = {'NDJ_left':cfg['bam_left'], 'NDJ_right':cfg['bam_right']}
    for text in a.control_bam:
        name, path = text.split('=',1)
        if name in bams: raise ValueError('Duplicate sample name')
        bams[name] = path
    records = []
    dependencies = [stamp(a.source/n, True) for n in ['config.json','COMPLETE.json']]
    dependencies += [stamp(ref), stamp(str(ref)+'.fai', True)] + segment_inputs
    for name, path in bams.items():
        idx = index_path(path)
        if idx is None: raise FileNotFoundError('Missing BAM index: '+path)
        with pysam.AlignmentFile(path, 'rb') as bam:
            bam.check_index()
            if bam.get_reference_length(CHROM) != length: raise ValueError('BAM reference mismatch: '+name)
            sq = next(x for x in bam.header.to_dict()['SQ'] if x['SN']==CHROM)
            rg = sorted({x.get('SM','UNKNOWN') for x in bam.header.to_dict().get('RG',[])})
        dependencies += [stamp(path), stamp(idx)]
        records.append(dict(biological_sample_id=name, interval=f'{CHROM}:{start}-{end}', reference=str(ref),
            action='recount_existing_alignment' if not a.skip_reads else 'reuse_reference_only',
            selected_resource_id=str(Path(path).resolve()), rejected_resource_ids='',
            evidence='frozen source config + BAM index and chr15 length; RG='+','.join(rg),
            unresolved_questions='junction unresolved; BAM M5 checked in worker if available',
            reference_md5=sq.get('M5',''), read_group_samples=','.join(rg)))
    annotations = None
    if a.annotations:
        annotations = read_json(a.annotations)
        if annotations['assembly'] != ASSEMBLY or Path(annotations['reference']).resolve() != ref.resolve():
            raise ValueError('Annotation manifest assembly/reference mismatch')
        dependencies.append(stamp(a.annotations, True))
        names = set()
        for track in annotations['tracks']:
            if track['name'] in names or not track['name'].replace('_','').isalnum():
                raise ValueError('Duplicate/unsafe annotation name')
            names.add(track['name']); dependencies.append(stamp(track['path'], True))
    repo = Path(__file__).resolve().parents[2]
    registry = repo/'manifests/known_analysis_sources.tsv'
    if registry.is_file(): dependencies.append(stamp(registry, True))
    ledger = a.project/'data/derived/resource_manifest.tsv'
    ledger_status = 'missing; frozen source config and BAMs independently validated'
    if ledger.is_file():
        # Snapshot the ledger, but do not make a mutable global registry a compute-cache dependency.
        with ledger.open() as f:
            matches = [r for r in csv.DictReader(f, delimiter='\t')
                       if r.get('path') in bams.values() or r.get('resolved_path') in
                       [str(Path(p).resolve()) for p in bams.values()]]
        ledger_status = dict(path=str(ledger),matching_rows=matches)
    return dict(project=str(a.project.resolve()), assembly=ASSEMBLY, chrom=CHROM, length=length, reference=str(ref.resolve()),
        start=start, end=end, origin=origin, bams=bams, reuse=records, annotations=annotations,
        ledger=ledger_status, inputs=dependencies, controls=a.controls, seed=a.seed,
        padding=a.padding, min_controls=a.min_controls, gc_caliper=a.gc_caliper,
        valid_caliper=.02, mapq=a.mapq, clip_min=20, skip_reads=a.skip_reads,
        versions=dict(python=sys.version, pysam=pysam.__version__),
        code=[stamp(Path(__file__),True),stamp(Path(__file__).with_name('sequence_features.py'),True)])


def load_tracks(manifest, length):
    tracks = {}
    for item in manifest['tracks'] if manifest else []:
        intervals = []
        with Path(item['path']).open() as f:
            for line in f:
                if not line.strip() or line.startswith(('#','track','browser')): continue
                fields = line.rstrip().split('\t')
                if fields[0] != CHROM: continue
                a,b = map(int,fields[1:3])
                if not 0 <= a < b <= length: raise ValueError('Annotation BED outside reference')
                intervals.append((a,b))
        if not intervals and not item.get('allow_empty_chr15',False):
            raise ValueError('No chr15 annotation records: '+item['name']+'; verify contig names')
        tracks['annotation_'+item['name']] = Track(intervals)
    return tracks


def background(s, tracks, a, b, cfg, rng):
    import numpy as np
    width = b-a; n = len(s)
    bases = np.frombuffer(s.upper().encode(), dtype='S1')
    gc_prefix = np.r_[np.int64(0), np.cumsum((bases==b'G') | (bases==b'C'))]
    valid_prefix = np.r_[np.int64(0), np.cumsum(np.isin(bases, [b'A', b'C', b'G', b'T']))]
    target = score(s, tracks, a, b)
    # Uniform starts, unique but possibly overlapping windows. These are genomic
    # reference ranks, not independent biological replicates.
    excluded_a = max(0, cfg['start']-max(cfg['padding'])-100000)
    excluded_b = min(n, cfg['end']+max(cfg['padding'])+100000)
    seen = set(); unconditioned = []; matched = []
    limit = max(10000, cfg['controls']*100)
    for _ in range(limit):
        x = rng.randrange(n-width+1)
        if x in seen: continue
        seen.add(x)
        if x < excluded_b and x+width > excluded_a: continue
        valid = int(valid_prefix[x+width]-valid_prefix[x])
        if valid/width < .9: continue
        gc = int(gc_prefix[x+width]-gc_prefix[x])/valid
        is_match = (target['GC'] is not None and abs(gc-target['GC']) <= cfg['gc_caliper']
            and abs(valid/width-target['valid_fraction']) <= cfg['valid_caliper'])
        if len(unconditioned)<cfg['controls'] or (is_match and len(matched)<cfg['controls']):
            v = score(s, tracks, x,x+width)
            row = dict(start=x,end=x+width,**v)
            if len(unconditioned)<cfg['controls']: unconditioned.append(row)
            if is_match and len(matched)<cfg['controls']: matched.append(row)
        if len(unconditioned) == cfg['controls'] and len(matched) == cfg['controls']: break
    return target, {'same_chr':unconditioned,'GC_valid_matched':matched}


def compare(target, controls, cfg, label):
    import numpy as np
    rows=[]
    for cohort, values in controls.items():
        for feature, observed in target.items():
            conditioned = cohort=='GC_valid_matched' and feature in {'GC','valid_fraction','ambiguous_fraction'}
            status = ('conditioning_variable' if conditioned else
                      'target_low_callability' if target['valid_fraction']<.9 else
                      'insufficient_controls' if len(values)<cfg['min_controls'] else 'exploratory')
            v = [r[feature] for r in values if r[feature] is not None]
            valid = status=='exploratory' and observed is not None and len(v)>=cfg['min_controls']
            hi = (1+sum(x>=observed for x in v))/(len(v)+1) if valid else None
            lo = (1+sum(x<=observed for x in v))/(len(v)+1) if valid else None
            rows.append(dict(scale=label,background=cohort,feature=feature,observed=observed,
                controls=len(v),median=float(np.median(v)) if v else None,
                percentile=100*(sum(x<observed for x in v)+.5*sum(x==observed for x in v))/len(v)
                    if v and observed is not None else None,
                p_high=hi,p_low=lo,p_two_sided=min(1,2*min(hi,lo)) if valid else None,
                q_BH_exploratory=None,status=status))
    return rows


def bh(rows):
    eligible = sorted([r for r in rows if r['p_two_sided'] is not None],key=lambda r:r['p_two_sided'])
    q=1.; total=len(eligible)
    for i in range(total-1,-1,-1):
        q=min(q,eligible[i]['p_two_sided']*total/(i+1))
        eligible[i]['q_BH_exploratory']=q


def read_evidence(cfg, out, seq):
    import pysam
    a=max(0,cfg['start']-100000); b=min(len(seq),cfg['end']+100000)
    totals=[]; events=[]; bins=[]
    # Query names retained only in local outputs, never committed. Supplementary
    # records are excluded; SA tags on primary records carry split-read evidence.
    for sample,path in cfg['bams'].items():
        counts={}; seen_events=set()
        with pysam.AlignmentFile(path,'rb') as bam:
            sq=next(x for x in bam.header.to_dict()['SQ'] if x['SN']==CHROM)
            if sq.get('M5') and sq['M5'].lower()!=hashlib.md5(seq.upper().encode()).hexdigest():
                raise ValueError('BAM M5 differs from reference: '+sample)
            for r in bam.fetch(CHROM,a,b):
                if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail: continue
                ix=max(0,(max(a,r.reference_start)-a)//1000)
                c=counts.setdefault(ix,dict(primary=0,low_mapq=0,SA=0,clipped=0,improper=0))
                c['primary']+=1
                if r.mapping_quality<cfg['mapq']:
                    c['low_mapq']+=1; continue
                sa=r.get_tag('SA') if r.has_tag('SA') else ''
                c['SA']+=bool(sa)
                c['improper']+=r.is_paired and not r.is_proper_pair
                cigar=r.cigartuples or []
                clipped=False
                for edge,ci,pos in [('left',0,r.reference_start),('right',-1,r.reference_end)]:
                    if cigar and cigar[ci][0]==4 and cigar[ci][1]>=cfg['clip_min'] and a<=pos<=b:
                        key=(r.query_name,r.is_read1,edge,pos)
                        if key in seen_events:continue
                        seen_events.add(key);clipped=True
                        clip=r.query_sequence[:cigar[ci][1]] if edge=='left' else r.query_sequence[-cigar[ci][1]:]
                        events.append(dict(sample=sample,read=r.query_name,read1=int(r.is_read1),
                            boundary_0based=pos,edge=edge,strand='-' if r.is_reverse else '+',
                            clip_bp=len(clip),mapq=r.mapping_quality,SA=sa,
                            mate_chrom=r.next_reference_name if r.next_reference_id>=0 else '*',
                            mate_start_0based=r.next_reference_start,clip_sequence=clip))
                c['clipped']+=clipped
                if sa and not clipped:
                    events.append(dict(sample=sample,read=r.query_name,read1=int(r.is_read1),
                        boundary_0based=r.reference_start,edge='SA_only_anchor_not_breakpoint',
                        strand='-' if r.is_reverse else '+',clip_bp=0,mapq=r.mapping_quality,SA=sa,
                        mate_chrom=r.next_reference_name if r.next_reference_id>=0 else '*',
                        mate_start_0based=r.next_reference_start,clip_sequence=''))
        for ix,c in sorted(counts.items()):
            bins.append(dict(sample=sample,start=a+ix*1000,end=min(b,a+(ix+1)*1000),**c))
        totals.append(dict(sample=sample,primary=sum(c['primary'] for c in counts.values()),
                           note='Counts are alignments, not independent molecules; bins use alignment start'))
    table(out/'read_evidence_1kb.tsv',bins,['sample','start','end','primary','low_mapq','SA','clipped','improper'])
    table(out/'candidate_reads.tsv',events,['sample','read','read1','boundary_0based','edge','strand',
          'clip_bp','mapq','SA','mate_chrom','mate_start_0based','clip_sequence'])
    clusters={}
    for r in events:
        if not r['clip_bp']:continue
        key=(r['sample'],r['boundary_0based'],r['edge'])
        clusters.setdefault(key,set()).add(r['read'])
    table(out/'clip_clusters.tsv',[dict(sample=k[0],boundary_0based=k[1],edge=k[2],fragments=len(v))
        for k,v in sorted(clusters.items(),key=lambda kv:-len(kv[1]))],
        ['sample','boundary_0based','edge','fragments'])
    return totals


def plot(out, rows, local, cfg):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    plt.rcParams.update({'pdf.fonttype':42,'font.size':8})
    with PdfPages(out/'breakpoint_context.pdf') as pdf:
        for bg in ['same_chr','GC_valid_matched']:
            selected=[r for r in rows if r['background']==bg and r['scale']=='pad_0' and r['percentile'] is not None]
            fig,ax=plt.subplots(figsize=(9,max(4,.24*len(selected))))
            ax.barh([r['feature'] for r in selected],[r['percentile'] for r in selected],color='#4C78A8')
            ax.axvline(50,color='black',lw=.6);ax.set_xlim(0,100)
            ax.set_xlabel('Empirical percentile (ties midrank; exploratory reference comparison)')
            ax.set_title(f'NDJ chr15 transition bracket | {bg}\nNo sequence-resolved junction inferred')
            fig.tight_layout();pdf.savefig(fig);plt.close(fig)
        keys=['GC','G4_canonical_fraction','inverted10_fraction','tandem2_6_fraction','ambiguous_fraction']
        fig,axes=plt.subplots(len(keys),1,figsize=(10,7),sharex=True)
        for ax,key in zip(axes,keys):
            ax.plot([(r['start']+r['end'])/2e6 for r in local],[r[key] for r in local],lw=.7)
            ax.axvspan(cfg['start']/1e6,cfg['end']/1e6,color='#F58518',alpha=.2)
            ax.set_ylabel(key.replace('_fraction',''),fontsize=7)
        axes[-1].set_xlabel('Chr15 position (Mb); 1-kb descriptive windows')
        fig.tight_layout();pdf.savefig(fig);plt.close(fig)



def register_outputs(cfg, out, outputs):
    """Append/update only this run's artifacts, preserving existing ledger rows."""
    import fcntl
    now=datetime.now(timezone.utc).isoformat()
    rows=[]
    for item in outputs:
        rows.append(dict(resource_id='NDJ_breakpoint:'+item['path'],biological_sample_id='NDJ_left;NDJ_right',
            run_accessions='',sample_aliases='',library_id='',role='breakpoint_context',
            resource_type='analysis_output',path=item['path'],resolved_path=item['path'],exists='true',
            bytes=item['bytes'],mtime_utc=datetime.fromtimestamp(item['mtime_ns']/1e9,timezone.utc).isoformat(),
            checksum=item['sha256'],checksum_type='sha256',checked_at_utc=now,
            source_config=str(out/'run.json'),source_config_field='inputs',producing_code_commit='see frozen code hashes in run.json',
            completion_marker=str(out/'complete.json'),validation_status='exploratory_context_complete',
            reference_accession=ASSEMBLY,contig=CHROM,start_1based=cfg['start']+1,end_1based=cfg['end'],
            interval_coverage_evidence='reference context and local alignment screen; see per-module status',
            sequencing_type='inherited source BAMs',sequencing_type_source='frozen source config',
            filters='see run.json',statistic='see feature_comparisons.tsv',read_group_samples='see reuse_decisions.tsv',indexed='NA'))
    table(out/'resource_manifest.tsv',rows)
    if not cfg.get('project'):return  # synthetic tests do not touch shared state
    ledger=Path(cfg['project'])/'data/derived/resource_manifest.tsv'
    ledger.parent.mkdir(parents=True,exist_ok=True)
    with Path(str(ledger)+'.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        previous=[];fields=list(rows[0])
        if ledger.exists():
            with ledger.open() as f:
                reader=csv.DictReader(f,delimiter='\t');previous=list(reader)
                fields=list(dict.fromkeys(list(reader.fieldnames or [])+fields))
            shutil.copy2(ledger,out/'resource_ledger_before_update.tsv')
        paths={r['path'] for r in rows}
        previous=[r for r in previous if r.get('path') not in paths]
        temp=Path(str(ledger)+'.'+str(os.getpid())+'.partial')
        table(temp,previous+rows,fields);temp.replace(ledger)


def work(out):
    import pysam
    cfg=read_json(out/'run.json')
    for old in cfg['inputs']:
        now=stamp(old['path'],'sha256' in old)
        if now!=old:raise RuntimeError('Input changed since submission: '+old['path'])
    for old in cfg['code']:
        if stamp(old['path'],True)!=old:raise RuntimeError('Worker code changed')
    with pysam.FastaFile(cfg['reference']) as fa: seq=fa.fetch(CHROM)
    print('STAGE\tsequence screens',flush=True)
    tracks=scan(seq); tracks.update(load_tracks(cfg['annotations'],len(seq)))
    distances=[]
    for name,t in tracks.items():
        nearest=min(((max(cfg['start']-y,x-cfg['end'],0),x,y) for x,y in zip(t.starts,t.ends)),default=None)
        distances.append(dict(feature=name,distance_to_bracket_bp=nearest[0] if nearest else None,
            nearest_start0=nearest[1] if nearest else None,nearest_end0=nearest[2] if nearest else None,
            status='present_on_chr15' if nearest else 'no_hits_on_chr15'))
    table(out/'feature_distances.tsv',distances)
    a=max(0,cfg['start']-100000);b=min(len(seq),cfg['end']+100000)
    with (out/'feature_intervals.bed').open('w') as f:
        for k,t in tracks.items():
            for x,y in t.intervals(a,b):f.write(f'{CHROM}\t{x}\t{y}\t{k}\n')
    rng=random.Random(cfg['seed']);rows=[];balance=[]
    for pad in cfg['padding']:
        x=max(0,cfg['start']-pad);y=min(len(seq),cfg['end']+pad)
        label=f'pad_{pad}'
        target,controls=background(seq,tracks,x,y,cfg,rng)
        rows.extend(compare(target,controls,cfg,label))
        table(out/(label+'_controls.tsv'),[dict(background=k,**r) for k,v in controls.items() for r in v])
        balance.append(dict(scale=label,start=x,end=y,**{k:len(v) for k,v in controls.items()}))
    bh(rows);table(out/'feature_comparisons.tsv',rows);table(out/'control_balance.tsv',balance)
    local=[dict(start=x,end=min(x+1000,b),**score(seq,tracks,x,min(x+1000,b))) for x in range(a,b,1000)]
    table(out/'local_1kb_features.tsv',local)
    reads=read_evidence(cfg,out,seq) if not cfg['skip_reads'] else 'not_assessed'
    # Reference queries for optional long-homology search, with explicit coordinates.
    with (out/'homology_queries.fa').open('w') as f:
        for x in range(a,b,2500):
            y=min(x+5000,b)
            if y-x>=200:f.write(f'>{CHROM}:{x}-{y}\n{seq[x:y]}\n')
    plot(out,rows,local,cfg)
    modules=dict(sequence='complete',annotations='provided' if cfg['annotations'] else 'not_assessed',
        reads=reads,long_homology='not_assessed; optional homology_search.py',
        junction_microhomology='not_assessed; requires assembled/read-supported alternate junction',
        replication_timing='not_assessed; cannot infer from reference sequence',
        chromatin_Rloops='not_assessed; motif and gene context are not measurements')
    js(out/'module_status.json',modules)
    report=f'''# NDJ chr15 breakpoint context\n\nReference: {ASSEMBLY}, {CHROM}.\nInterval (0-based half-open): [{cfg['start']}, {cfg['end']}).\n{cfg['origin']}.\n\nSee feature_comparisons.tsv for observed fractions, background medians, midrank\npercentiles and exploratory two-sided empirical p/BH values across all scales and\nbackgrounds. They are reference-sequence comparisons, not a test of a causal\nmechanism or a population of breakpoints. Overlapping controls are not biological\nreplicates. GC-matched results condition on GC; use same_chr to assess GC itself.\n\nLocal 1-kb peaks are descriptive, not independently tested. A motif anywhere in\nthe bracket cannot be assigned to the break until the junction is resolved.\nLow MAPQ, gaps, duplicated sequence and alignment artifacts must be checked\nbefore interpreting coverage or clipping. Repeated reads can derive from shared\nfragments. Left/right evidence may reflect tissue mixture.\n\nRead module: {'skipped' if cfg['skip_reads'] else 'complete'}.\nAnnotation tracks: {'supplied' if cfg['annotations'] else 'NOT ASSESSED; softmask is not a TE annotation'}.\nLong homology and junction mechanisms remain unassessed in this first pass.\n'''
    (out/'REPORT.md').write_text(report)
    outputs=[stamp(p,True) for p in sorted(out.iterdir()) if p.is_file() and p.name not in
             {'complete.json','run.json','jobs.json'} and p.suffix in {'.tsv','.bed','.pdf','.md','.fa','.json'}]
    register_outputs(cfg,out,outputs)
    outputs.append(stamp(out/'resource_manifest.tsv',True))
    js(out/'complete.json',dict(status='exploratory_context_complete',modules=modules,outputs=outputs))
    print('COMPLETE\t'+str(out),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project',type=Path,default=PROJECT)
    p.add_argument('--source',type=Path)
    p.add_argument('--figure4',type=Path)
    p.add_argument('--output',type=Path)
    p.add_argument('--interval',type=int,nargs=2,metavar=('START0','END0'))
    p.add_argument('--annotations',type=Path)
    p.add_argument('--control-bam',action='append',default=[],metavar='NAME=PATH')
    p.add_argument('--controls',type=int,default=999)
    p.add_argument('--min-controls',type=int,default=100)
    p.add_argument('--gc-caliper',type=float,default=.03)
    p.add_argument('--padding',type=int,nargs='+',default=[0,5000,50000])
    p.add_argument('--seed',type=int,default=150116)
    p.add_argument('--mapq',type=int,default=30)
    p.add_argument('--skip-reads',action='store_true')
    p.add_argument('--partition',default='nano');p.add_argument('--time',default='00:30:00')
    p.add_argument('--submit',action='store_true');p.add_argument('--run',action='store_true')
    a=p.parse_args()
    if a.run:
        if not a.output:p.error('--run needs --output')
        work(a.output.resolve());return
    if not a.submit:p.error('Use --submit on Pegasus')
    if a.controls<a.min_controls or a.min_controls<10 or min(a.padding)<0 or 0 not in a.padding or not 0<a.gc_caliper<=1:
        p.error('Invalid controls, padding or GC caliper')
    a.source=a.source or a.project/'analysis/NDJ_LR_difference_scan_v1/run_20260920T012415Z'
    a.figure4=a.figure4 or a.project/'analysis/Figure4_NDJ_panels'
    cfg=resolve(a)
    digest=hashlib.sha256(json.dumps({k:v for k,v in cfg.items() if k!='ledger'},sort_keys=True).encode()).hexdigest()[:12]
    out=(a.output or a.project/'analysis/NDJ_chr15_breakpoint_context_v1'/digest).resolve()
    out.mkdir(parents=True,exist_ok=True)
    if (out/'complete.json').exists():
        previous=read_json(out/'run.json')
        # Code path differs after freezing, but code hash and input signatures remain.
        previous.pop('worker_code',None)
        if previous.get('request_digest')!=digest:raise RuntimeError('Output belongs to different inputs/options')
        for old in read_json(out/'complete.json')['outputs']:
            if stamp(old['path'],True)!=old:raise RuntimeError('Completed output changed: '+old['path'])
        print('REUSED\t'+str(out));return
    if (out/'run.json').exists():raise RuntimeError('Incomplete/existing run: inspect logs and use a new --output directory')
    cfg['request_digest']=digest
    (out/'code').mkdir(exist_ok=True);(out/'logs').mkdir(exist_ok=True)
    for name in ['run_breakpoint.py','sequence_features.py']:
        shutil.copy2(Path(__file__).with_name(name),out/'code'/name)
    cfg['code']=[stamp(p,True) for p in sorted((out/'code').glob('*.py'))]
    js(out/'run.json',cfg)
    table(out/'reuse_decisions.tsv',cfg['reuse'])
    js(out/'resource_ledger_snapshot.json',cfg['ledger'])
    cmd=[sys.executable,str(out/'code/run_breakpoint.py'),'--run','--output',str(out)]
    job=subprocess.check_output(['sbatch','--parsable','--job-name=NDJ_breakpoint','--partition='+a.partition,
        '--time='+a.time,'--cpus-per-task=1','--mem=12G','--output='+str(out/'logs/%j.out'),
        '--error='+str(out/'logs/%j.err'),'--wrap='+shlex.join(cmd)],text=True).strip().split(';')[0]
    js(out/'jobs.json',dict(job=job,submitted=datetime.now(timezone.utc).isoformat()))
    print('JOB\t'+job+'\nOUTPUT\t'+str(out))


if __name__=='__main__':main()
