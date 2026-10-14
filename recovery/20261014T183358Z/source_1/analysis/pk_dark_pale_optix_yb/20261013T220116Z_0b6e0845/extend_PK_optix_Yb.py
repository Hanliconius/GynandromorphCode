#!/usr/bin/env python3
"""Post hoc dark/pale PK extension: immutable site selection, arrays, full-genome mapping."""
import argparse,csv,gzip,importlib.util,json,os,shlex,shutil,subprocess,sys,tempfile,uuid
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import pandas as pd
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
HERE=Path(__file__).resolve().parent
WINDOWS={'Hmel2.5':('Hmel215003o',1395000,1440000),'HmRos':('NC_140403.1',1620000,1670000),'Optix':('Hmr18',1060000,1240000)}
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def run(args,**kw):return subprocess.run([str(x) for x in args],check=True,**kw)
def header(path):
    with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as h:
        for line in h:
            if line.startswith('#CHROM'):return line.rstrip().split('\t')[9:]
    raise RuntimeError('VCF header missing')
def download_new(a):
    """Start raw-read acquisition independently of reference preparation."""
    samples=list(csv.DictReader((HERE/'pk_dark_pale_samples.tsv').open(),delimiter='\t'))
    if len(samples)!=21 or len({s['sample'] for s in samples})!=21:
        raise RuntimeError('Expected exactly 21 distinct dark/pale PKs')
    vcf=ROOT/'external/krista_genotype_plots_original/Genotype_Plots/Optix_Hmr18-938195-1239421_biallelic.vcf'
    existing=set(header(vcf))
    samples=[s for s in samples if s['sample'] not in existing]
    if not samples:
        print('DOWNLOAD\tNo additional PK samples');return
    base=ROOT/'analysis/pk_dark_pale_downloads'/(datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8])
    (base/'logs').mkdir(parents=True)
    shutil.copy2(__file__,base/Path(__file__).name)
    shutil.copy2(HERE.parents[1]/'yb/submit_Yb_extra_controls.py',base/'download_worker.py')
    cfg={'run_dir':str(base),'samples':samples,'fastq_root':str(ROOT/'analysis/yb_extra_control_fastq')}
    (base/'run.json').write_text(json.dumps(cfg,indent=2)+'\n')
    with (base/'samples.tsv').open('w') as h:
        writer=csv.DictWriter(h,fieldnames=list(samples[0]),delimiter='\t');writer.writeheader();writer.writerows(samples)
    worker=base/'download.sh'
    worker.write_text('#!/usr/bin/env bash\nset -euo pipefail\nexec '+shlex.join([sys.executable,str(base/Path(__file__).name),'--config',str(base/'run.json'),'--stage','download'])+' --task "$SLURM_ARRAY_TASK_ID"\n')
    jobid=subprocess.check_output(['sbatch','--parsable','--partition=cpu',f'--array=0-{len(samples)-1}','--cpus-per-task=1','--mem=4G','--time=12:00:00','--job-name=PK_download',f'--output={base}/logs/download.%A_%a.out',f'--error={base}/logs/download.%A_%a.err',str(worker)],text=True).strip().split(';')[0]
    (base/'jobs.json').write_text(json.dumps({'download':jobid},indent=2)+'\n')
    (base.parent/'latest_run.path').write_text(str(base)+'\n')
    print('DOWNLOAD_ARRAY\t'+jobid);print('NEW_PK_SAMPLES\t'+str(len(samples)));print('RUN_DIR\t'+str(base))
def submit(a):
    samples=list(csv.DictReader((HERE/'pk_dark_pale_samples.tsv').open(),delimiter='\t'))
    if len(samples)!=21 or len({s['sample'] for s in samples})!=21:raise RuntimeError('Expected exactly 21 distinct dark/pale PKs')
    focal=a.focal_run or Path((ROOT/'analysis/yb_local_both_references/latest_run.path').read_text().strip())
    fc=json.loads((focal/'run.json').read_text());tools=fc['tools']
    helper=module('align_tools',HERE.parents[1]/'yb/rebuild_Yb_local_both_references.py')
    tools={**tools,'bcftools':helper.tool('bcftools')}
    import pysam
    vcf=ROOT/'external/krista_genotype_plots_original/Genotype_Plots/Optix_Hmr18-938195-1239421_biallelic.vcf'
    existing=set(header(vcf));missing=[s['sample'] for s in samples if s['sample'] not in existing]
    values=pd.read_csv(a.optix_run/'plot_values.tsv.gz',sep='\t')
    keys=['chrom','pos','ref','alt'];sites=values[keys].drop_duplicates().sort_values(keys)
    refs=[dict(name=r['name'],fasta=str(focal/r['name']/'reference/genome.fa'),contig=r['contig']) for r in fc['references']]
    if missing:
        # Hmr18 is the legacy chromosome name; test aliases in prepared rosina first.
        candidates=([a.optix_reference] if a.optix_reference else
                    [Path(r['fasta']) for r in refs if r['name']=='HmRos'])
        if not a.optix_reference:
            for searchroot in [ROOT/'data',ROOT/'reference',ROOT/'references',ROOT/'external',ROOT/'analysis/genomewide_mixing_input/reference']:
                if searchroot.exists():
                    candidates.extend(Path(str(fai)[:-4]) for fai in searchroot.rglob('*.fai'))
        optixref,source_contig,audit=resolve_optix_alias(candidates,sites)
        print('OPTIX_ALIAS\t'+source_contig+' -> Hmr18; exact frozen REF agreement',flush=True)
        refs.append(dict(name='Optix',fasta=str(optixref),contig='Hmr18',source_contig=source_contig,alias_audit=audit))
    base=ROOT/'analysis/pk_dark_pale_optix_yb'/(datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')+'_'+uuid.uuid4().hex[:8]);(base/'logs').mkdir(parents=True)
    shutil.copy2(__file__,base/Path(__file__).name);shutil.copy2(HERE.parents[1]/'yb/submit_Yb_extra_controls.py',base/'download_worker.py')
    shutil.copy2(HERE/'pk_dark_pale_samples.tsv',base/'samples.tsv')
    frozen_reports=base/'code/scripts/figures/figure3';frozen_reports.mkdir(parents=True)
    for name in ['Fig3_optix_BA099_both_callable.py','Fig3_optix_full_range_frequency.py','replot_optix_pruned.py','replot_optix_DR_callable_PK.py','replot_optix_DR_genoplot.py']:
        shutil.copy2(HERE/name,frozen_reports/name)

    values.to_csv(base/'frozen_optix_values.tsv.gz',sep='\t',index=False)
    sites[['chrom','pos']].assign(start=sites.pos-1)[['chrom','start','pos']].drop_duplicates().to_csv(base/'optix_sites.bed',sep='\t',index=False,header=False)
    alleles=[]
    for (chrom,pos,ref),group in sites.groupby(['chrom','pos','ref'],sort=True):
        alleles.append(f"{chrom}\t{pos}\t{ref},"+','.join(sorted(set(group.alt)))+'\n')
    (base/'optix_alleles.tsv').write_text(''.join(alleles))
    pysam.tabix_compress(str(base/'optix_alleles.tsv'),str(base/'optix_alleles.tsv.gz'),force=True)
    pysam.tabix_index(str(base/'optix_alleles.tsv.gz'),seq_col=0,start_col=1,end_col=1,zerobased=False,force=True)
    reuse={};waitjobs=set();fastq=ROOT/'analysis/yb_extra_control_fastq' 
    # Reuse existing five-PK alignments, including queued tasks, instead of duplicating them.
    for old in sorted((ROOT/'analysis/yb_extra_controls').glob('*/run.json')):
        cfg=json.loads(old.read_text());jobs=json.loads((old.parent/'jobs.json').read_text())
        for oldindex,task in enumerate(cfg.get('tasks',[])):
            acc=task['sample'];ref=task['reference']['name'];bam=old.parent/ref/'coverage'/(acc+'.dedup.bam')
            if acc not in {s['sample'] for s in samples}:continue
            if ref not in {r['name'] for r in refs if r['name']!='Optix'}:continue
            marker=bam.parent/(acc+'.complete.json')
            if marker.is_file() and bam.is_file():reuse[(ref,acc)]={'bam':str(bam),'dependency':None}
            elif jobs.get('align') and (ref,acc) not in reuse:
                state=subprocess.check_output(['sacct','-X','-n','-j',jobs['align']+'_'+str(oldindex),'--format=State','--parsable2'],text=True)
                if any(x.strip().split('|')[0] in ['PENDING','RUNNING','CONFIGURING','COMPLETING'] for x in state.splitlines()):reuse[(ref,acc)]={'bam':str(bam),'dependency':jobs['align']+'_'+str(oldindex)}
    tasks=[]
    for ref in refs:
        for s in samples:
            if ref['name']=='Optix' and s['sample'] not in missing:continue
            old=reuse.get((ref['name'],s['sample']))
            if old and old['dependency']:waitjobs.add(old['dependency'])
            tasks.append(dict(reference=ref,sample=s['sample'],reuse=old))
    needed={t['sample'] for t in tasks if not t['reuse']}
    for acc in needed:
        for mate in [1,2]:
            target=fastq/acc/f'{acc}_{mate}.fastq.gz'
            local=ROOT/'analysis/genomewide_mixing_input/fastq'/acc/f'{acc}_{mate}.fastq.gz'
            if not target.exists() and local.is_file() and local.stat().st_size:
                target.parent.mkdir(parents=True,exist_ok=True);target.symlink_to(local)

    download_wait={}
    for source in sorted((ROOT/'analysis/pk_dark_pale_downloads').glob('*/run.json')):
        previous=json.loads(source.read_text());previous_jobs=source.parent/'jobs.json'
        if not previous_jobs.exists() or previous.get('fastq_root')!=str(fastq):continue
        jid=json.loads(previous_jobs.read_text()).get('download')
        if not jid:continue
        states=subprocess.check_output(['sacct','-X','-n','-j',jid,'--format=JobID,State','--parsable2'],text=True)
        active={line.split('|')[0].strip() for line in states.splitlines() if len(line.split('|'))>1 and line.split('|')[1].strip() in {'PENDING','RUNNING','CONFIGURING','COMPLETING'}}
        for index,sample in enumerate(previous['samples']):
            element=jid+'_'+str(index)
            if sample['sample'] in needed and element in active:download_wait[sample['sample']]=element
    cfg={'run_dir':str(base),'repo':str(base/'code'),'samples':samples,'tasks':tasks,'references':refs,'tools':tools,'fastq_root':str(fastq),'optix_vcf':str(vcf),'optix_existing_samples':[s['sample'] for s in samples if s['sample'] in existing],'existing_download_dependencies':download_wait,
        'manifest_source':'user-provided dark/pale PK table; cydno, crosses, backcrosses, ivory-only excluded','site_selection_changed':False,'PK_used_for_ascertainment':False,'coverage_QC':'samtools depth -aa -s -q0 -Q30 -G1796; mean including zero-depth bases'}
    (base/'run.json').write_text(json.dumps(cfg,indent=2)+'\n')
    worker=base/'worker.sh';worker.write_text('#!/usr/bin/env bash\nset -euo pipefail\nexec '+shlex.join([sys.executable,str(base/Path(__file__).name),'--config',str(base/'run.json')])+' --stage "$1" --task "${SLURM_ARRAY_TASK_ID:-0}"\n')
    jobs={}
    def job(stage,ids,cpus,mem,limit,dependency=None,partition='cpu',key=None):
        args=['sbatch','--parsable',f'--partition={partition}',f'--array={ids}',f'--cpus-per-task={cpus}',f'--mem={mem}',f'--time={limit}',f'--job-name=PK_{stage}',f'--output={base}/logs/{stage}.%A_%a.out',f'--error={base}/logs/{stage}.%A_%a.err']
        if dependency:args.append('--dependency='+dependency)
        j=subprocess.check_output(args+[str(worker),stage],text=True).strip().split(';')[0];jobs[key or stage]=j;(base/'jobs.json').write_text(json.dumps(jobs,indent=2)+'\n');print((key or stage).upper()+'\t'+j,flush=True);return j
    new=[i for i,t in enumerate(tasks) if not t['reuse']];old=[i for i,t in enumerate(tasks) if t['reuse']]
    try:
        if missing:
            prep=job('prepare',str(len(refs)-1),8,'32G','06:00:00')
        dl=[i for i,s in enumerate(samples) if s['sample'] in needed and s['sample'] not in download_wait]
        prerequisites=sorted(set(download_wait.values()))
        if dl:prerequisites.append(job('download',','.join(map(str,dl)),1,'4G','12:00:00'))
        if missing:prerequisites.append(prep)
        deps=[]
        if new:
            d='afterok:'+':'.join(prerequisites) if prerequisites else None
            align=job('align',','.join(map(str,new)),16,'64G','2-00:00:00',d)
            deps.append(job('profile',','.join(map(str,new)),2,'8G','04:00:00','aftercorr:'+align,key='new_profiles'))
        if old:deps.append(job('profile',','.join(map(str,old)),2,'8G','04:00:00',('afterok:'+':'.join(sorted(waitjobs))) if waitjobs else None,key='reused_profiles'))
        job('report','0',1,'8G','00:30:00','afterok:'+':'.join(deps),partition='nano')
    except Exception:
        for j in jobs.values():subprocess.run(['scancel',j],check=False)
        raise
    (base.parent/'latest_run.path').write_text(str(base)+'\n')
    print('RUN_DIR\t'+str(base));print(f'SCOPE\t21 PK; {len(new)} new locus alignment tasks; {len(old)} reused; {len(missing)} optix samples require new calls')
def resolve_optix_alias(candidates,sites):
    import pysam
    if set(sites.chrom.astype(str))!={'Hmr18'}:raise RuntimeError('Expected legacy Hmr18 optix sites')
    bases=sites[['pos','ref']].drop_duplicates()
    if bases.pos.duplicated().any():raise RuntimeError('Conflicting frozen REF bases')
    lo=int(bases.pos.min())-1;hi=int(bases.pos.max())
    offsets=bases.pos.to_numpy(dtype=int)-1-lo;expected=bases.ref.str.upper().to_numpy()
    if any(len(x)!=1 for x in expected):raise RuntimeError('Alias check requires SNP REF bases')
    seen=set();audit=[]
    for candidate in candidates:
        candidate=Path(candidate).resolve()
        if candidate in seen or not Path(str(candidate)+'.fai').is_file():continue
        seen.add(candidate)
        with pysam.FastaFile(str(candidate)) as fa:
            if len(fa.references)<10:continue
            matches=[]
            for contig,length in zip(fa.references,fa.lengths):
                if length<max(1250000,hi):continue
                sequence=fa.fetch(contig,lo,hi).upper()
                count=sum(sequence[int(i)]==base for i,base in zip(offsets,expected))
                audit.append(dict(fasta=str(candidate),contig=contig,matched=count,total=len(expected),length=length))
                if count==len(expected):matches.append(contig)
            if len(matches)==1:return candidate,matches[0],audit
            if len(matches)>1:raise RuntimeError('Ambiguous optix chromosome aliases: '+repr(matches))
    best=sorted(audit,key=lambda x:x['matched']/x['total'],reverse=True)[:3]
    raise RuntimeError('No chromosome matches all frozen optix REF bases at unchanged coordinates. Best candidates='+repr(best)+'; chromosome names alone cannot establish coordinate equivalence.')
def prepare(cfg,task):
    import pysam
    ref=cfg['references'][task];d=Path(cfg['run_dir'])/ref['name']/'reference';d.mkdir(parents=True,exist_ok=True);fa=d/'genome.fa'
    renamed=ref.get('source_contig',ref['contig'])!=ref['contig']
    if not fa.exists():
        if renamed:
            temporary=d/'genome.fa.partial'
            with pysam.FastaFile(ref['fasta']) as source,temporary.open('w') as out:
                if ref['contig'] in source.references:raise RuntimeError('Chromosome rename would create duplicate FASTA IDs')
                for name,length in zip(source.references,source.lengths):
                    out.write('>'+ (ref['contig'] if name==ref['source_contig'] else name)+'\n')
                    for start in range(0,length,60000):
                        sequence=source.fetch(name,start,min(length,start+60000))
                        out.write('\n'.join(sequence[i:i+60] for i in range(0,len(sequence),60))+'\n')
            temporary.replace(fa)
        else:fa.symlink_to(ref['fasta'])
    if not renamed:
        for suffix in ['.fai','.amb','.ann','.bwt','.pac','.sa']:
            source=Path(ref['fasta']+suffix)
            if source.is_file() and not Path(str(fa)+suffix).exists():Path(str(fa)+suffix).symlink_to(source)
    if not Path(str(fa)+'.fai').exists():run([cfg['tools']['samtools'],'faidx',fa])
    if not Path(str(fa)+'.bwt').exists():run([cfg['tools']['bwa'],'index',fa])
    (d/'chromosome_alias.json').write_text(json.dumps(ref,indent=2)+'\n')
def align(cfg,task):
    t=cfg['tasks'][task];ref=t['reference'];acc=t['sample'];tools=cfg['tools'];base=Path(cfg['run_dir']);d=base/ref['name']/acc;d.mkdir(parents=True,exist_ok=True)
    if (d/'alignment.complete.json').is_file():
        if ref['name']!='HmRos' or not cfg.get('shared_rosina_optix') or (d/'optix.bam').is_file():return
    fa=base/'Optix/reference/genome.fa' if ref['name']=='Optix' else Path(ref['fasta'])
    inputs=[Path(cfg['fastq_root'])/acc/f'{acc}_{i}.fastq.gz' for i in [1,2]]
    with tempfile.TemporaryDirectory(prefix=acc+'.',dir=d) as temp:
        tmp=Path(temp);raw=tmp/'full.sorted.bam';marked=tmp/'full.marked.bam'
        cmd=shlex.join([tools['bwa'],'mem','-M','-t','12','-R',f'@RG\\tID:{acc}\\tSM:{acc}',str(fa),*map(str,inputs)])+' | '+shlex.join([tools['samtools'],'sort','-@','3','-m','2G','-o',str(raw),'-'])
        run(['bash','-o','pipefail','-c',cmd])
        picard=[tools['picard']] if 'picard' in tools else [tools['java'],'-Xmx48g','-jar',tools['picard_jar']]
        run(picard+['MarkDuplicates',f'I={raw}',f'O={marked}',f'M={d}/duplicate_metrics.txt','VALIDATION_STRINGENCY=LENIENT'])
        run([tools['samtools'],'index','-@','3',marked])
        # Retain whole Yb contig for defensible contig-mean normalization; only optix span for SNP calling.
        region=ref['contig'] if ref['name']!='Optix' else 'Hmr18:1050000-1250000'
        retained=tmp/'locus.bam';run([tools['samtools'],'view','-b','-o',retained,marked,region]);run([tools['samtools'],'index',retained])
        retained.replace(d/'locus.bam');Path(str(retained)+'.bai').replace(d/'locus.bam.bai')
        if ref['name']=='HmRos' and cfg.get('shared_rosina_optix'):
            optix=tmp/'optix.bam'
            chrom=cfg['shared_rosina_optix']['contig']
            run([tools['samtools'],'view','-b','-o',optix,marked,f'{chrom}:1050000-1250000'])
            run([tools['samtools'],'index',optix])
            optix.replace(d/'optix.bam')
            Path(str(optix)+'.bai').replace(d/'optix.bam.bai')
    (d/'alignment.complete.json').write_text(json.dumps({'sample':acc,'full_reference':ref['fasta'],'retained_region':region,'full_genome_alignment':True,'duplicates_marked_before_region_extraction':True},indent=2)+'\n')
def profile(cfg,task):
    t=cfg['tasks'][task];ref=t['reference']['name'];acc=t['sample'];base=Path(cfg['run_dir']);d=base/ref/acc;d.mkdir(parents=True,exist_ok=True)
    if (d/'profile.complete.json').is_file():return
    bam=Path(t['reuse']['bam']) if t['reuse'] else d/'locus.bam';tools=cfg['tools'];contig,lo,hi=WINDOWS[ref]
    if ref=='Optix':
        out=d/'calls.vcf'
        if cfg.get('shared_rosina_optix'):
            fa=Path(cfg['shared_rosina_optix']['fasta'])
            cmd=shlex.join([tools['bcftools'],'mpileup','-Ou','-f',str(fa),'-R',str(base/'optix_sites_rosina.bed'),'-a','FORMAT/DP,FORMAT/AD','-q','30','-Q','20',str(bam)])+' | '+shlex.join([tools['bcftools'],'call','-m','-A','-C','alleles','-T',str(base/'optix_alleles_rosina.tsv.gz'),'-f','GQ','-Ov','-o',str(d/'calls.rosina.vcf')])
            run(['bash','-o','pipefail','-c',cmd])
            run([tools['bcftools'],'annotate','--rename-chrs',base/'rosina_to_optix_chrom.tsv','-Ov','-o',out,d/'calls.rosina.vcf'])
            (d/'calls.rosina.vcf').unlink()
        else:
            fa=base/'Optix/reference/genome.fa'
            command=shlex.join([tools['bcftools'],'mpileup','-Ou','-f',str(fa),'-R',str(base/'optix_sites.bed'),'-a','FORMAT/DP,FORMAT/AD','-q','30','-Q','20',str(bam)])+' | '+shlex.join([tools['bcftools'],'call','-m','-A','-C','alleles','-T',str(base/'optix_alleles.tsv.gz'),'-f','GQ','-Ov','-o',str(out)])
            run(['bash','-o','pipefail','-c',command])
    else:
        hdr=subprocess.check_output([tools['samtools'],'view','-H',str(bam)],text=True)
        length=next(int(dict(z.split(':',1) for z in line.split('\t')[1:])['LN']) for line in hdr.splitlines() if line.startswith('@SQ\t') and 'SN:'+contig in line.split('\t'))
        depth=np.zeros(length,dtype=np.int32);proc=subprocess.Popen([tools['samtools'],'depth','-aa','-s','-q','0','-Q','30','-G','1796','-r',contig,str(bam)],stdout=subprocess.PIPE,text=True)
        count=0
        try:
            for line in proc.stdout:
                chrom,pos,n=line.rstrip().split('\t')
                if chrom!=contig or int(pos)!=count+1:raise RuntimeError('Depth stream not continuous')
                depth[count]=int(n);count+=1
            if proc.wait()!=0 or count!=length:raise RuntimeError('Depth task incomplete')
        finally:
            if proc.poll() is None:proc.kill();proc.wait()
        allele_span=f'{contig}:{lo+1}-{hi}'
        callcmd=shlex.join([tools['bcftools'],'mpileup','-Ou','-f',t['reference']['fasta'],'-r',allele_span,'-a','FORMAT/DP,FORMAT/AD','-q','30','-Q','20',str(bam)])+' | '+shlex.join([tools['bcftools'],'call','-m','-A','-f','GQ','-Oz','-o',str(d/'Yb_calls.vcf.gz')])
        run(['bash','-o','pipefail','-c',callcmd]);run([tools['bcftools'],'index','-t',d/'Yb_calls.vcf.gz'])
        mean=float(depth.mean())
        if mean<=0:raise RuntimeError('Zero Yb contig depth')
        starts=np.arange(lo,hi,50);means=np.array([depth[x:x+50].mean() for x in starts])
        pd.DataFrame({'position':starts+25.5,'mean_depth':means,'normalized_depth':means/mean}).to_csv(d/'coverage.tsv.gz',sep='\t',index=False)
        (d/'normalization.json').write_text(json.dumps({'contig':contig,'length':length,'mean_contig_depth':mean,'source_bam':str(bam),'QC':cfg['coverage_QC']},indent=2)+'\n')
    (d/'profile.complete.json').write_text(json.dumps(t,indent=2)+'\n')
def calls_from_new(path,wanted):
    result={}
    indexed={}
    for key in wanted:indexed.setdefault(key[:3],[]).append(key)
    with path.open() as h:
        for line in h:
            if line.startswith('#'):continue
            f=line.rstrip().split('\t');pos=int(f[1]);key=(f[0],pos,f[3]);fmt=dict(zip(f[8].split(':'),f[9].split(':')))
            gt=fmt.get('GT','.')
            try:dp=float(fmt.get('DP','0'))
            except ValueError:dp=0
            try:gq=float(fmt.get('GQ','0'))
            except ValueError:gq=0
            for chrom,sitepos,ref,alt in indexed.get(key,[]):
                dose=np.nan
                if dp>=2 and gq>=20:
                    alleles=f[4].split(',')
                    bases=[ref]+alleles
                    try:
                        observed=[bases[int(x)] for x in gt.replace('|','/').split('/')]
                        if len(observed)==2 and all(x in (ref,alt) for x in observed):dose=observed.count(alt)/2
                    except (ValueError,IndexError):pass
                result[(chrom,sitepos,ref,alt)]=dose
    return result
def report(cfg):
    base=Path(cfg['run_dir']);samples=[s for s in cfg['samples'] if s['sample'] not in cfg.get('excluded_samples',[])];values=pd.read_csv(base/'frozen_optix_values.tsv.gz',sep='\t');keys=['chrom','pos','ref','alt']
    repo=Path(cfg['repo']);sys.path.insert(0,str(repo/'scripts/figures/figure3'))
    import Fig3_optix_full_range_frequency as full
    from replot_optix_pruned import select
    from replot_optix_DR_callable_PK import callable_sites
    from replot_optix_DR_genoplot import draw
    from matplotlib.backends.backend_pdf import PdfPages
    selected=values[keys].drop_duplicates();wanted=list(selected.itertuples(index=False,name=None))
    existing=cfg['optix_existing_samples'];calls,_=full.read_calls(Path(cfg['optix_vcf']),selected,existing) if existing else ({},set())
    for acc in [s['sample'] for s in samples if s['sample'] not in existing]:
        calls.update({(key,acc):dose for key,dose in calls_from_new(base/'Optix'/acc/'calls.vcf',wanted).items()})
    labels={s['sample']:f"{s['label']} {s['phenotype']} | B{'+' if s['B']=='1' else '-'} D{'+' if s['D']=='1' else '-'} R{'+' if s['R']=='1' else '-'}" for s in samples}
    ids=[s['sample'] for s in samples];plotvalues=[]
    with PdfPages(base/'PK_dark_pale_optix_Yb_posthoc.pdf') as pdf:
        for modules in [['Ray','Dennis'],['Band']]:
            for cutoff in [.8,.9]:
                v=callable_sites(select(values.loc[values.module.isin(modules)],cutoff,True));records=[]
                for _,site in v.loc[v.row.eq('TN_Left')].iterrows():
                    key=tuple(site[k] for k in keys);is_alt=full.base.pattern_is_alt(site.selected_orientation)
                    for acc in ids:
                        dose=calls.get((key,acc),np.nan);records.append({**site.to_dict(),'row':acc,'frequency':dose if is_alt else 1-dose,'source':'post hoc individual genotype; DP2 GQ20; no ascertainment'})
                expanded=pd.concat([v,pd.DataFrame(records)],ignore_index=True);plotvalues.append(expanded.assign(BA_cutoff=cutoff))
                draw(expanded,cutoff,True,base,pdf,pk_samples=ids,focal_callable=True,row_labels=labels,modules=modules)
        import matplotlib.pyplot as plt
        for ref in ['Hmel2.5','HmRos']:
            for phenotype in ['dark','pale']:
                group=[s for s in samples if s['phenotype']==phenotype]
                for start in range(0,len(group),6):
                    subset=group[start:start+6];fig,axes=plt.subplots(len(subset),1,figsize=(11,2+len(subset)*1.1),sharex=True,squeeze=False)
                    fig.subplots_adjust(left=.24,right=.96,top=.9,bottom=.12,hspace=.35)
                    for ax,s in zip(axes[:,0],subset):
                        t=pd.read_csv(base/ref/s['sample']/'coverage.tsv.gz',sep='\t');ax.plot(t.position/1e6,t.normalized_depth,lw=.8,color='#662A83' if phenotype=='dark' else '#CC70A7');ax.axhline(1,color='#aaa',lw=.5);ax.set_ylim(0,2);ax.set_ylabel(s['label']+' '+s['sample'],rotation=0,ha='right',va='center',fontsize=8)
                        ax.spines[['top','right']].set_visible(False)
                    axes[-1,0].set_xlabel(WINDOWS[ref][0]+' position (Mb)');fig.suptitle(f'Yb | {ref} | PK {phenotype} | individual normalized mean coverage')
                    fig.text(.24,.03,'50-bp mean / own full-contig mean; MAPQ>=30; duplicate exclusion; paired overlaps counted once. Range 0-2.',fontsize=8);pdf.savefig(fig);plt.close(fig)
    pd.concat(plotvalues,ignore_index=True).to_csv(base/'optix_posthoc_plot_values.tsv.gz',sep='\t',index=False)
    (base/'complete.json').write_text(json.dumps({'samples':len(samples),'excluded_samples':cfg.get('excluded_samples',[]),'PK_used_for_ascertainment':False,'PDF':str(base/'PK_dark_pale_optix_Yb_posthoc.pdf'),'Yb_phenotypes':'not inferred from ivory or optix states; coverage rows shown individually','optix_existing':existing,'new_optix_QC':'Known-allele-constrained diploid bcftools calls; DP>=2 GQ>=20; MAPQ30 BQ20','focal_mixing':'existing reciprocal-mixture calls; no second correction'},indent=2)+'\n')
    print('PDF\t'+str(base/'PK_dark_pale_optix_Yb_posthoc.pdf'))
def main():
    p=argparse.ArgumentParser();p.add_argument('--submit',action='store_true');p.add_argument('--download-new',action='store_true',help='Submit new PK raw-read download array without requiring reference genomes');p.add_argument('--optix-run',type=Path);p.add_argument('--focal-run',type=Path);p.add_argument('--optix-reference',type=Path)
    p.add_argument('--config',type=Path);p.add_argument('--stage',choices=['prepare','download','align','profile','report']);p.add_argument('--task',type=int,default=0);a=p.parse_args()
    if a.download_new:download_new(a)
    elif a.submit:submit(a)
    else:
        cfg=json.loads(a.config.read_text())
        if a.stage=='download':
            module('download_worker',Path(cfg['run_dir'])/'download_worker.py').download(cfg,a.task)
            sample=cfg['samples'][a.task]
            record=json.loads((Path(cfg['fastq_root'])/sample['sample']/'download.complete.json').read_text())
            if record['ENA']['sample_accession']!=sample['biosample']:raise RuntimeError('ENA BioSample mismatch')
        elif a.stage=='prepare':prepare(cfg,a.task)
        elif a.stage=='align':align(cfg,a.task)
        elif a.stage=='profile':profile(cfg,a.task)
        else:report(cfg)
if __name__=='__main__':main()
