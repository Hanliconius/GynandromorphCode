"""Small synthetic fixture, never presented as real NDJ results."""
from pathlib import Path
import json
import numpy as np
import pysam


def make(root):
    root=Path(root);root.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(71)
    chroms=['chrA','chrB','chrZ'];length=5000
    seq='ACGT'*(length//4)
    fasta=root/'reference.fa'
    fasta.write_text(''.join('>'+c+'\n'+seq+'\n' for c in chroms));pysam.faidx(str(fasta))
    header=pysam.VariantHeader()
    for c in chroms:header.contigs.add(c,length=length)
    header.formats.add('AD','R','Integer','Allelic depths')
    header.add_sample('NDJ_left');header.add_sample('NDJ_right')
    bcf=root/'paired.vcf'
    with pysam.VariantFile(str(bcf),'w',header=header) as out:
        for ci,c in enumerate(chroms):
            for pos in range(5,length-5,5):
                ref=seq[pos-1];alt={'A':'G','G':'A','C':'T','T':'C'}[ref]
                p=rng.choice([.0,.5,1.0],p=[.35,.5,.15]);pl=pr=p
                if ci==0 and 1500<=pos<2500:pl,pr=.8,.2
                n1,n2=60,55
                a,b=rng.binomial(n1,pl),rng.binomial(n2,pr)
                rec=out.new_record(contig=c,start=pos-1,alleles=(ref,alt))
                rec.samples['NDJ_left']['AD']=(n1-a,a)
                rec.samples['NDJ_right']['AD']=(n2-b,b)
                if pos==10:rec.samples['NDJ_left']['AD']=(None,None)
                out.write(rec)
    bgz=root/'paired.vcf.gz'
    pysam.tabix_compress(str(bcf),str(bgz),force=True);pysam.tabix_index(str(bgz),preset='vcf',force=True)
    bams=[]
    for side in ['left','right']:
        path=root/(side+'.bam')
        hdr={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':c,'LN':length} for c in chroms], 'RG':[{'ID':'rg1','SM':'NDJ_'+side}]}
        with pysam.AlignmentFile(str(path),'wb',header=hdr) as out:
            for ci,c in enumerate(chroms):
                for start in range(0,length-150,25):
                    n=3 if side=='right' and ci==0 and 1500<=start<2500 else 8
                    for rep in range(n):
                        read=pysam.AlignedSegment();read.query_name=f'{ci}_{start}_{rep}'
                        read.query_sequence=seq[start:start+150];read.flag=16 if rep%2 else 0
                        read.reference_id=ci;read.reference_start=start;read.mapping_quality=60;read.cigar=[(0,150)]
                        read.query_qualities=pysam.qualitystring_to_array('I'*150);read.set_tag('RG','rg1');out.write(read)
        pysam.index(str(path));bams.append(path)
    samtools=root/'samtools'
    import sys
    samtools.write_text('#!'+sys.executable+'\nimport sys,pysam\nif sys.argv[1]=="--version":\n print("samtools 1.23 (pysam test wrapper)")\nelse:\n sys.stdout.write(getattr(pysam.samtools,sys.argv[1])(*sys.argv[2:]))\n')
    samtools.chmod(0o755)
    return ['--vcf',str(bgz),'--reference',str(fasta),'--bam-left',str(bams[0]),'--bam-right',str(bams[1]),'--samtools',str(samtools),
            '--contig','chrA','--contig','chrB','--contig','chrZ','--exclude-background','chrZ','--window','500',
            '--min-window-sites','30','--background-sites','3000','--candidates','20']

if __name__=='__main__':
    import sys
    print(json.dumps(make(Path(sys.argv[1]).resolve())))
