#!/usr/bin/env python3
"""Audit local clip candidates and SA partners; no automatic junction/mechanism call."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

def read(path):
    return json.loads(Path(path).read_text())

def stamp(path):
    p=Path(path).resolve(strict=True); st=p.stat()
    return dict(path=str(p),bytes=st.st_size,mtime_ns=st.st_mtime_ns)

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write_table(path, rows, fields):
    with Path(path).open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n')
        w.writeheader();w.writerows(rows)

def geometry(cigar, strand):
    """Query alignment interval in original-read orientation, including hard clips."""
    ops=re.findall(r'(\d+)([MIDNSHP=X])',cigar)
    if not ops or ''.join(n+op for n,op in ops)!=cigar:
        raise ValueError('Malformed CIGAR: '+cigar)
    ops=[(int(n),op) for n,op in ops]
    length=sum(n for n,op in ops if op in 'MIS H=X'.replace(' ',''))
    leading=0
    for n,op in ops:
        if op not in 'SH':break
        leading+=n
    trailing=0
    for n,op in reversed(ops):
        if op not in 'SH':break
        trailing+=n
    start,end=leading,length-trailing
    if strand=='-':start,end=length-end,length-start
    return start,end,length,sum(n for n,op in ops if op in 'MDN=X')

def sa_support(cigar, strand, edge, tag):
    pstart,pend,length,_=geometry(cigar,strand)
    toward_low=(edge=='left') == (strand=='+')
    clip=(0,pstart) if toward_low else (pend,length)
    boundary=pstart if toward_low else pend
    records=[]
    for entry in tag.rstrip(';').split(';'):
        if not entry:continue
        chrom,pos,orient,scigar,mq,nm=entry.split(',')
        a,b,total,span=geometry(scigar,orient)
        if total!=length:continue
        overlap=max(0,min(b,clip[1])-max(a,clip[0]))
        if int(mq)<30 or overlap<20 or overlap<.5*(clip[1]-clip[0]):continue
        # Nearest query endpoint to primary alignment; map it to a reference boundary.
        low_endpoint=abs(a-boundary)<=abs(b-boundary)
        partner=int(pos)-1
        if low_endpoint == (orient=='-'):partner+=span
        gap=(pstart-b) if toward_low else (a-pend)
        records.append(dict(partner_chrom=chrom,partner_boundary0=partner,
            partner_strand=orient,primary_strand=strand,query_gap_bp=gap,
            clipped_query_overlap_bp=overlap,SA_mapq=int(mq),SA_nm=int(nm)))
    return records

def selftest():
    assert geometry('20S80M','+')==(20,100,100,80)
    assert geometry('20S80M','-')==(0,80,100,80)
    assert geometry('10H10S80M','+')==(20,100,100,80)
    x=sa_support('80M20S','+','right','chrB,501,+,80S20M,60,0;')[0]
    assert x['partner_boundary0']==500 and x['query_gap_bp']==0
    x=sa_support('80M20S','+','right','chrB,501,-,20M80S,60,0;')[0]
    assert x['partner_boundary0']==520
    assert not sa_support('80M20S','+','right','chrB,501,+,80M20S,60,0;')
    assert not sa_support('80M20S','+','right','chrB,501,+,80S20M,10,0;')

def candidates(source,cfg):
    with (source/'clip_clusters.tsv').open() as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    counts=defaultdict(dict)
    for r in rows:
        counts[(int(r['boundary_0based']),r['edge'])][r['sample']]=int(r['fragments'])
    inside=[];outside=[]
    for (pos,edge),s in counts.items():
        if max(s.values())<3:continue
        row=dict(boundary0=pos,edge=edge,original_max=max(s.values()),
                 scope='bracket' if cfg['start']<=pos<cfg['end'] else 'flank_comparison')
        (inside if row['scope']=='bracket' else outside).append(row)
    # Keep exact original sites; nearby rows can reuse the same fragments.
    outside=sorted(outside,key=lambda r:(-r['original_max'],r['boundary0']))[:20]
    return sorted(inside+outside,key=lambda r:(r['boundary0'],r['edge']))

def inspect(cfg,sites,out):
    import pysam
    summaries=[];partners=[]
    for sample,path in cfg['bams'].items():
        with pysam.AlignmentFile(path,'rb') as bam:
            if bam.get_reference_length(cfg['chrom'])!=cfg['length']:
                raise ValueError('BAM/reference contig mismatch')
            for site in sites:
                pos=site['boundary0']; edge=site['edge']
                eligible=set();clipped=set();good=set();tagged=set()
                groups=defaultdict(set);starts=defaultdict(set);gaps=defaultdict(list)
                sequences=Counter()
                for r in bam.fetch(cfg['chrom'],max(0,pos-1000),min(cfg['length'],pos+1000)):
                    if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail:
                        continue
                    if r.mapping_quality<30:continue
                    # Unique RG+QNAME fragments intersecting the +/-5 bp anchor neighbourhood.
                    if r.reference_start>pos+5 or r.reference_end<pos-5:continue
                    ident=(r.get_tag('RG') if r.has_tag('RG') else '',r.query_name)
                    eligible.add(ident)
                    cigar=r.cigartuples or []
                    ci=0 if edge=='left' else -1
                    endpoint=r.reference_start if edge=='left' else r.reference_end
                    if not cigar or cigar[ci][0]!=4 or cigar[ci][1]<20 or abs(endpoint-pos)>5:continue
                    clipped.add(ident);n=cigar[ci][1]
                    qualities=r.query_qualities
                    q=qualities[:n] if edge=='left' and qualities is not None else (
                        qualities[-n:] if qualities is not None else [])
                    if not len(q) or sum(q)/len(q)<20:continue
                    good.add(ident)
                    seq=r.query_sequence[:n] if edge=='left' else r.query_sequence[-n:]
                    sequences[hashlib.sha256(seq.encode()).hexdigest()]+=1
                    if not r.has_tag('SA'):continue
                    tagged.add(ident)
                    for hit in sa_support(r.cigarstring,'-' if r.is_reverse else '+',edge,r.get_tag('SA')):
                        # Exact partner endpoint first; nearby endpoints remain explicit.
                        key=(hit['partner_chrom'],hit['partner_boundary0'],
                             hit['partner_strand'],hit['primary_strand'])
                        groups[key].add(ident)
                        starts[key].add((r.reference_start,r.reference_end,r.is_reverse))
                        gaps[key].append(hit['query_gap_bp'])
                for key,fragments in sorted(groups.items(),key=lambda kv:-len(kv[1])):
                    partners.append(dict(sample=sample,**site,partner_chrom=key[0],
                        partner_boundary0=key[1],partner_strand=key[2],primary_strand=key[3],
                        fragments=len(fragments),distinct_anchor_geometries=len(starts[key]),
                        query_gap_min=min(gaps[key]),query_gap_max=max(gaps[key])))
                summaries.append(dict(sample=sample,**site,eligible_fragments=len(eligible),
                    clipped_fragments=len(clipped),quality_clip_fragments=len(good),
                    clip_fraction=len(clipped)/len(eligible) if eligible else '',
                    SA_tagged_clip_fragments=len(tagged),distinct_exact_partner_groups=len(groups),
                    largest_exact_partner_group=max(map(len,groups.values()),default=0),
                    distinct_clip_sequence_hashes=len(sequences)))
    write_table(out/'candidate_summary.tsv',summaries,
        ['sample','boundary0','edge','original_max','scope','eligible_fragments','clipped_fragments',
         'quality_clip_fragments','clip_fraction','SA_tagged_clip_fragments','distinct_exact_partner_groups',
         'largest_exact_partner_group','distinct_clip_sequence_hashes'])
    write_table(out/'SA_partner_groups.tsv',partners,
        ['sample','boundary0','edge','original_max','scope','partner_chrom','partner_boundary0',
         'partner_strand','primary_strand','fragments','distinct_anchor_geometries','query_gap_min','query_gap_max'])
    return summaries,partners

def worker(out):
    plan=read(out/'plan.json');source=Path(plan['source']);cfg=read(source/'run.json')
    for item in plan['inputs']:
        if stamp(item['path'])!=item:raise RuntimeError('Input changed: '+item['path'])
    if sha(Path(__file__))!=plan['code_sha256']:raise RuntimeError('Worker changed')
    selftest()
    summaries,partners=inspect(cfg,plan['sites'],out)
    core=[r for r in partners if r['scope']=='bracket']
    count=sum(r['fragments']>=3 and r['distinct_anchor_geometries']>=2 for r in core)
    text=f"""# NDJ junction follow-up
Screened {len(plan['sites'])} exact clipping sites, including up to 20 strong flank comparison sites.
Bracket partner groups with >=3 unique fragments and >=2 anchor geometries: {count}.
This is a descriptive review threshold, not a significance test or a confirmed junction.
If zero, the appropriate conclusion is: no obvious coherent alternate junction in this SA-tag screen.
A positive count only prioritizes manual review. Inspect orientation, partner positions and query gaps.

Clipping fractions use RG+read-name fragments intersecting the +/-5-bp anchor neighbourhood.
Primary MAPQ >=30; duplicate/QC-fail/secondary/supplementary records excluded.
Quality clips require >=20 bases and mean clip base quality >=20.
SA partners require MAPQ >=30 and alignment of >=20 bases and >=half the clipped query region.
Partner groups use exact endpoints; nearby alternatives and overlapping candidate rows are not independent.
Flank comparison sites were selected for strong clipping, not random null controls.
Negative query gaps are alignment overlaps, NOT demonstrated microhomology.
Missing SA tags can reflect aligner limitations; no SA support does not prove no rearrangement.
Short clips can map ambiguously and small indels can yield locally consistent SA groups.
No de novo assembly, TE annotation or new whole-genome alignment is performed.
Raw read sequences/names are not included in this review package.
"""
    (out/'REPORT.md').write_text(text)
    # Register small review artifacts without replacing unrelated ledger records.
    import fcntl
    ledger=Path(cfg['project'])/'data/derived/resource_manifest.tsv'
    ledger.parent.mkdir(parents=True,exist_ok=True)
    new=[]
    for name in ['REPORT.md','candidate_summary.tsv','SA_partner_groups.tsv']:
        path=out/name
        new.append(dict(resource_id='NDJ_junction:'+str(path),biological_sample_id='NDJ_left;NDJ_right',
            resource_type='analysis_output',path=str(path),resolved_path=str(path),exists='true',
            bytes=path.stat().st_size,checksum=sha(path),checksum_type='sha256',
            source_config=str(out/'plan.json'),completion_marker=str(out/'complete.json'),
            validation_status='review_ready_not_validated_junction',reference_accession=cfg['assembly'],
            contig=cfg['chrom'],start_1based=cfg['start']+1,end_1based=cfg['end'],
            statistic='descriptive clipping and SA partner groups'))
    write_table(out/'resource_manifest.tsv',new,list(new[0]))
    with Path(str(ledger)+'.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        old=[];fields=list(new[0])
        if ledger.exists():
            with ledger.open() as f:
                reader=csv.DictReader(f,delimiter='\t');old=list(reader)
                fields=list(dict.fromkeys(list(reader.fieldnames or [])+fields))
            shutil.copy2(ledger,out/'ledger_before_update.tsv')
        paths={r['path'] for r in new}
        temp=out/'ledger.partial.tsv'
        write_table(temp,[r for r in old if r.get('path') not in paths]+new,fields)
        # Output and ledger are under the same project filesystem.
        temp.replace(ledger)
    outputs={p.name:sha(p) for p in out.iterdir() if p.suffix in {'.tsv','.md'} and p.name!='ledger_before_update.tsv'}
    (out/'complete.json').write_text(json.dumps(dict(status='review_ready_not_validated_junction',
        outputs=outputs,qualifying_descriptive_groups=count),indent=2)+'\n')
    print('COMPLETE\t'+str(out))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling/analysis/NDJ_chr15_breakpoint_context_v1/ec8e51556aa7'))
    p.add_argument('--output',type=Path)
    p.add_argument('--submit',action='store_true')
    p.add_argument('--worker',action='store_true')
    p.add_argument('--self-test',action='store_true')
    a=p.parse_args()
    selftest()
    if a.self_test:print('SELF_TESTS_PASS');return
    if a.worker:
        if a.output is None:p.error('--worker requires --output')
        worker(a.output.resolve());return
    if not a.submit:p.error('Use --submit on Pegasus')
    import pysam
    source=a.source.resolve();cfg=read(source/'run.json')
    if read(source/'complete.json')['status']!='exploratory_context_complete':
        raise ValueError('Source not complete')
    inputs=[stamp(source/n) for n in ['run.json','complete.json','clip_clusters.tsv']]
    # Reuse exactly the already-audited BAMs; do not discover/download replacements.
    for path in cfg['bams'].values():
        old=next(x for x in cfg['inputs'] if x['path']==str(Path(path).resolve()))
        if stamp(path)!={k:old[k] for k in ['path','bytes','mtime_ns']}:
            raise ValueError('BAM changed since first pass')
        inputs.append(stamp(path))
        idx=next((x for x in [Path(path+'.bai'),Path(path).with_suffix('.bai'),Path(path+'.csi')] if x.is_file()),None)
        if idx is None:raise FileNotFoundError('Missing BAM index')
        inputs.append(stamp(idx))
        with pysam.AlignmentFile(path,'rb') as b:
            b.check_index()
            if b.get_reference_length(cfg['chrom'])!=cfg['length']:raise ValueError('Contig mismatch')
    plan=dict(source=str(source),inputs=inputs,code_sha256=sha(__file__),sites=candidates(source,cfg),
              pysam=pysam.__version__,python=sys.version)
    digest=hashlib.sha256(json.dumps(plan,sort_keys=True).encode()).hexdigest()[:12]
    out=(a.output or source/('junction_followup_'+digest)).resolve()
    if (out/'complete.json').exists():
        if read(out/'plan.json')!=plan:raise RuntimeError('Different existing run')
        for name,value in read(out/'complete.json')['outputs'].items():
            if sha(out/name)!=value:raise RuntimeError('Output changed')
        print('REUSED\t'+str(out));return
    if out.exists():raise RuntimeError('Incomplete output exists; inspect logs then use a fresh --output')
    out.mkdir(parents=True);shutil.copy2(__file__,out/'worker.py')
    (out/'plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    write_table(out/'reuse_decisions.tsv',[dict(sample=s,interval=f"{cfg['chrom']}:{cfg['start']}-{cfg['end']}",
        action='recount_existing_alignment',path=p,evidence='unchanged audited source BAM; same reference; indexed')
        for s,p in cfg['bams'].items()],['sample','interval','action','path','evidence'])
    command=[sys.executable,str(out/'worker.py'),'--worker','--output',str(out)]
    job=subprocess.check_output(['sbatch','--parsable','--partition=nano','--time=00:30:00',
        '--cpus-per-task=1','--mem=4G','--job-name=NDJ_junction',
        '--output='+str(out/'job.%j.out'),'--error='+str(out/'job.%j.err'),
        '--wrap='+shlex.join(command)],text=True).strip().split(';')[0]
    (out/'job.json').write_text(json.dumps(dict(job=job))+'\n')
    print('JOB\t'+job+'\nOUTPUT\t'+str(out))

if __name__=='__main__':main()
