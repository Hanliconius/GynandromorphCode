#!/usr/bin/env python3
"""Exploratory Optix merge: exact shared SNPs, recalled samples take precedence.
Does not recall genotypes or change the corrected focal mixture results.
"""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import re

ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
PARENT=ROOT/'analysis/rosina_genomewide_joint_recall/20261014T223341Z_f956c8bc'
CHROM='NC_140406.1'

def focal(name):
    return bool(re.match(r'^TN(?:[23][LR](?:_|$)|_(?:Left|Right)(?:_|$))',name,re.I))

def records(path):
    op=gzip.open if str(path).endswith('.gz') else open
    with op(path,'rt') as h:
        for line in h:
            if line.startswith('#CHROM'):
                yield 'samples',line.rstrip().split('\t')[9:]
            elif not line.startswith('#'):
                f=line.rstrip().split('\t')
                if f[0] not in ('Hmr18',CHROM,'18'):continue
                pos=int(f[1])
                if not 1060000<=pos<=1240000:continue
                if len(f[3])!=1 or len(f[4])!=1 or f[3] not in 'ACGT' or f[4] not in 'ACGT':continue
                fmt=f[8].split(':')
                yield (pos,f[3],f[4]),[dict(zip(fmt,x.split(':'))) for x in f[9:]]

def gt(call):
    s=call.get('GT','.')
    a=re.split(r'[/|]',s)
    return tuple(sorted(map(int,a))) if len(a)==2 and all(x in ('0','1') for x in a) else None

def passed(call):
    try:return gt(call) is not None and int(call.get('DP','.'))>=2 and float(call.get('GQ','.'))>=20
    except (ValueError,TypeError):return False

def run(a):
    a.output.mkdir(parents=True,exist_ok=False)
    recalled={};newnames=None
    for key,row in records(a.recalled):
        if key=='samples':newnames=row;continue
        if key in recalled:raise RuntimeError(f'Duplicate recalled SNP: {key}')
        recalled[key]=row
    if not newnames or len(newnames)!=len(set(newnames)):raise RuntimeError('Invalid recalled sample IDs')
    if any(focal(n) for n in newnames):raise RuntimeError('Unexpected focal sample in comparative recalled VCF')
    original={};oldnames=None;seen=set();original_sites=0
    for key,row in records(a.original):
        if key=='samples':oldnames=row;continue
        if key in seen:raise RuntimeError(f'Duplicate original SNP: {key}')
        seen.add(key);original_sites+=1
        if key in recalled:original[key]=row
    if not oldnames or len(oldnames)!=len(set(oldnames)):raise RuntimeError('Invalid original sample IDs')
    overlap=sorted(set(oldnames)&set(newnames))
    extra=[n for n in oldnames if n not in newnames and not focal(n)]
    names=newnames+extra
    oi={n:i for i,n in enumerate(oldnames)};ni={n:i for i,n in enumerate(newnames)}
    stats={n:Counter() for n in names};agreement=Counter();shared=sorted(original)
    if not shared:raise RuntimeError('No exact shared SNPs; chromosome/coordinate/allele compatibility unresolved')
    # Keep raw GT/DP/GQ/AD; absent fields remain absent. No fabricated reference calls.
    out=a.output/'Optix_legacy_plus_recalled.vcf'
    with out.open('w') as h:
        h.write('##fileformat=VCFv4.2\n##source=exploratory_legacy_plus_recalled_Optix\n')
        h.write(f'##contig=<ID={CHROM}>\n')
        for name,num,typ,desc in [('GT','1','String','Genotype'),('DP','1','Integer','Read depth'),('GQ','1','Float','Genotype quality'),('AD','R','Integer','Allele depths')]:
            h.write(f'##FORMAT=<ID={name},Number={num},Type={typ},Description="{desc}">\n')
        h.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(names)+'\n')
        for key in shared:
            old=original[key];new=recalled[key]
            calls=new+[old[oi[n]] for n in extra]
            for n,c in zip(names,calls):
                stats[n]['sites']+=1;stats[n]['GT_present']+=gt(c) is not None
                stats[n]['DP2_GQ20_callable']+=passed(c)
                stats[n]['GQ_missing']+=c.get('GQ','.')=='.'
            for n in overlap:
                x,y=old[oi[n]],new[ni[n]]
                if gt(x) is not None and gt(y) is not None:
                    agreement['both_GT_present']+=1;agreement['discordant_GT']+=gt(x)!=gt(y)
                if passed(x) and passed(y):
                    agreement['both_DP2_GQ20']+=1;agreement['discordant_DP2_GQ20']+=gt(x)!=gt(y)
            fields=[':'.join(c.get(k,'.') for k in ('GT','DP','GQ','AD')) for c in calls]
            h.write(f'{CHROM}\t{key[0]}\t.\t{key[1]}\t{key[2]}\t.\t.\t.\tGT:DP:GQ:AD\t'+'\t'.join(fields)+'\n')
    import pysam
    dest=str(out)+'.gz';pysam.tabix_compress(str(out),dest,force=False);pysam.tabix_index(dest,preset='vcf',force=False);out.unlink()
    with (a.output/'sample_callability.tsv').open('w') as h:
        h.write('sample\tcall_source\tsites\tGT_present\tDP2_GQ20_callable\tGQ_missing\n')
        for n,c in stats.items():h.write(n+'\t'+('recalled' if n in ni else 'legacy')+'\t'+'\t'.join(str(c[k]) for k in ('sites','GT_present','DP2_GQ20_callable','GQ_missing'))+'\n')
    audit=dict(status='EXPLORATORY_MERGE_COMPLETE',original=str(a.original),recalled=str(a.recalled),
        recalled_samples=len(newnames),original_samples=len(oldnames),overlapping_sample_IDs=len(overlap),
        added_legacy_samples=len(extra),merged_samples=len(names),excluded_focal=[n for n in oldnames if focal(n)],
        original_biallelic_window_SNPs=original_sites,recalled_biallelic_window_SNPs=len(recalled),
        exact_shared_SNPs=len(shared),recalled_site_retention=len(shared)/len(recalled),
        genotype_agreement=dict(agreement),calling_histories_uniform=False,
        site_policy='Exact position REF ALT intersection; Hmr18 renamed NC_140406.1; no allele swaps or coordinate shifts',
        sample_policy='Exact-ID deduplication; recalled takes precedence; focal tracks excluded; phenotype eligibility not inferred',
        corrected_focal_calls='Unchanged; not included in this comparative VCF',
        vcf=dest,vcf_sha256=hashlib.sha256(Path(dest).read_bytes()).hexdigest())
    (a.output/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(f"SAMPLES\t{len(names)} ({len(newnames)} recalled + {len(extra)} legacy)")
    print(f"SHARED_SNPS\t{len(shared)} / {len(recalled)} recalled ({len(shared)/len(recalled):.1%})")
    print('VCF\t'+dest);print('AUDIT\t'+str(a.output/'audit.json'))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--original',type=Path,default=ROOT/'external/krista_genotype_plots_original/Genotype_Plots/Optix_Hmr18-938195-1239421_biallelic.vcf')
    p.add_argument('--recalled',type=Path,default=PARENT/'integration_Optix/Optix_window.vcf.gz')
    p.add_argument('--output',type=Path,default=PARENT/'Optix_legacy_merge_v1')
    run(p.parse_args())
