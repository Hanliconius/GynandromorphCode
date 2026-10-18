#!/usr/bin/env python3
"""Read-only project-wide Yb BAM/VCF/BCF audit. No downloads, alignments or calling."""
import argparse,csv,json,os,re
from pathlib import Path
import pysam
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
CHROMS={'NC_140403.1':('HmRos',[(1633041,1635646),(1658079,1659072)]),'Hmel215003o':('Hmel2.5',[(1403500,1405500),(1429500,1430500)]),'Hmr15':('legacy_Rosina_alias_UNVERIFIED',[(1633041,1635646),(1658079,1659072)])}
def tsv(path,rows,cols):
    with path.open('w',newline='') as h:
        w=csv.DictWriter(h,delimiter='\t',fieldnames=cols,extrasaction='ignore');w.writeheader();w.writerows(rows)
def callable_gt(c):
    gt=c.get('GT') or ()
    return len(gt)==2 and all(x in (0,1) for x in gt) and (c.get('DP') or 0)>=2 and c.get('GQ') is not None and c['GQ']>=20

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    out=a.output;out.mkdir(parents=True,exist_ok=False)
    meta=list(csv.DictReader((ROOT/'data/derived/region_panel_metadata.tsv').open(),delimiter='\t'))
    lookup={}
    for r in meta:
        for k in ('Accession','SequenceID','EarthCapeID'):
            for s in re.findall(r'(?:ERR|SRR|DRR)\d+',r.get(k,'')):lookup.setdefault(s,[]).append(r)
    bamrows=[];vcfrows=[];samples=[];errors=[];seen=set();files=[]
    for folder in (ROOT/'analysis',ROOT/'external',ROOT/'data'):
        if not folder.exists():continue
        for parent,dirs,names in os.walk(folder,followlinks=False):
            dirs[:]=[d for d in dirs if d not in {'logs','.git','__pycache__','raw','fastq','code','toolenv','windows'}]
            for n in names:
                if n.endswith(('.bam','.cram','.vcf','.vcf.gz','.bcf')):
                    path=Path(parent)/n
                    try:
                        real=path.resolve(strict=True)
                        if real in seen:continue
                        seen.add(real);files.append(path)
                    except OSError as e:errors.append(dict(path=str(path),error=str(e)))
    # BAMs: RG identity, dictionary, usable index, actual reads in both intervals.
    for path in files:
        if path.suffix not in ('.bam','.cram'):continue
        try:
            with pysam.AlignmentFile(str(path),'rc' if path.suffix=='.cram' else 'rb') as bam:
                sm=sorted({r.get('SM','') for r in bam.header.to_dict().get('RG',[]) if r.get('SM')})
                identities=sm if len(sm)==1 else ['UNRESOLVED_RG']
                dic=dict(zip(bam.references,bam.lengths));indexed=bam.has_index()
                for chrom,(reference,intervals) in CHROMS.items():
                    if chrom not in dic:continue
                    for label,(lo,hi) in zip(('d','v'),intervals):
                        covered=dic[chrom]>=hi
                        reads=bam.count(chrom,lo-1,hi,read_callback='all') if indexed and covered else None
                        for sample in identities:
                            bamrows.append(dict(sample=sample,path=str(path),resolved_path=str(path.resolve()),reference=reference,chrom=chrom,interval=label,contig_length=dic[chrom],indexed=indexed,reads_overlapping=reads,RG_samples=','.join(sm),status='candidate_requires_reference_dictionary_verification' if indexed and reads else 'unusable_or_no_reads'))
        except Exception as e:errors.append(dict(path=str(path),error=str(e)))
    # Native indexed callsets and small legacy locus VCFs. Never stream genome VCFs.
    for path in files:
        if not path.name.endswith(('.vcf','.vcf.gz','.bcf')):continue
        try:
            with pysam.VariantFile(str(path)) as vf:
                names=list(vf.header.samples)
                if not names:continue
                chroms=set(vf.header.contigs)&set(CHROMS)
                if not chroms:continue
                indexed=any(Path(str(path)+suffix).is_file() for suffix in ('.csi','.tbi'))
                small=path.stat().st_size<100_000_000
                stats={(c,i):dict(records=0,called={s:0 for s in names}) for c in chroms for i in ('d','v')}
                if indexed:
                    def records():
                        for c in chroms:
                            for lo,hi in CHROMS[c][1]:
                                for r in vf.fetch(c,lo-1,hi):yield r
                    iterator=records()
                elif small:iterator=iter(vf)
                else:
                    vcfrows.append(dict(path=str(path),reference='unknown',chrom=','.join(sorted(chroms)),interval='not_scanned_unindexed_large',records='',samples=len(names),indexed=False));continue
                for rec in iterator:
                    if rec.contig not in chroms or len(rec.ref)!=1 or not rec.alts or len(rec.alts)!=1 or len(rec.alts[0])!=1:continue
                    for label,(lo,hi) in zip(('d','v'),CHROMS[rec.contig][1]):
                        if not lo<=rec.pos<=hi:continue
                        stat=stats[(rec.contig,label)];stat['records']+=1
                        for s in names:stat['called'][s]+=callable_gt(rec.samples[s])
                for (chrom,label),stat in stats.items():
                    vcfrows.append(dict(path=str(path),reference=CHROMS[chrom][0],chrom=chrom,interval=label,records=stat['records'],samples=len(names),indexed=indexed))
                    if not stat['records']:continue
                    for s,n in stat['called'].items():
                        rr=lookup.get(s,[])
                        samples.append(dict(sample=s,reference=CHROMS[chrom][0],interval=label,path=str(path),SNPs=stat['records'],callable_DP2_GQ20=n,callable_fraction=n/stat['records'],taxon='; '.join(sorted({r.get('Taxon Name','') for r in rr})),Yb_dorsal='; '.join(sorted({r.get('Yb-dorsal','') for r in rr})),Yb_ventral='; '.join(sorted({r.get('Yb-ventral','') for r in rr}))))
        except Exception as e:errors.append(dict(path=str(path),error=str(e)))
    tsv(out/'bam_candidates.tsv',bamrows,['sample','path','resolved_path','reference','chrom','interval','contig_length','indexed','reads_overlapping','RG_samples','status'])
    tsv(out/'callset_inventory.tsv',vcfrows,['path','reference','chrom','interval','records','samples','indexed'])
    tsv(out/'sample_callability.tsv',samples,['sample','reference','interval','path','SNPs','callable_DP2_GQ20','callable_fraction','taxon','Yb_dorsal','Yb_ventral'])
    tsv(out/'read_errors.tsv',errors,['path','error'])
    # One row per biological ID and reference/interval; retain source alternatives separately.
    plan=[]
    allids=set(lookup)|{r['sample'] for r in bamrows if r['sample']!='UNRESOLVED_RG'}|{r['sample'] for r in samples}
    for s in sorted(allids):
        for label in ('d','v'):
            calls=[r for r in samples if r['sample']==s and r['interval']==label and r['callable_DP2_GQ20']>0]
            bams=[r for r in bamrows if r['sample']==s and r['interval']==label and r['indexed'] and r['reads_overlapping']]
            rr=lookup.get(s,[]);taxon='; '.join(sorted({r.get('Taxon Name','') for r in rr}))
            plan.append(dict(sample=s,interval=label,taxon=taxon,existing_callable_VCFs=len(calls),existing_BAMs_with_reads=len(bams),
                action='harmonize_existing_calls' if calls else 'reuse_existing_alignment_for_counts' if bams else 'no_usable_data_found_in_searched_roots',
                VCFs='; '.join(r['path'] for r in calls),BAMs='; '.join(r['path'] for r in bams)))
    tsv(out/'sample_reuse_plan.tsv',plan,['sample','interval','taxon','existing_callable_VCFs','existing_BAMs_with_reads','action','VCFs','BAMs'])
    summary=dict(status='COMPLETE',root=str(ROOT),files_examined=len(files),errors=len(errors),reference_aliases_verified=False,downloads=0,alignments=0,calling=0,
        warning='Candidate calls must be harmonized against reference/alleles; preserve restored focal calling independently.',counts={})
    for label in ('d','v'):
        q=[r for r in plan if r['interval']==label]
        summary['counts'][label]={action:sum(r['action']==action for r in q) for action in sorted({r['action'] for r in q})}
    (out/'complete.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('OUTPUT\t'+str(out))
    for label,c in summary['counts'].items():print('Yb('+label+')\t'+json.dumps(c))
    print('ERRORS_RECORDED\t'+str(len(errors)))
if __name__=='__main__':main()
