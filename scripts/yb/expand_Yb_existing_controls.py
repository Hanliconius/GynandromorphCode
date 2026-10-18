#!/usr/bin/env python3
"""Expand frozen Yb panels using existing comparative VCFs and arrayed BAM recalls.
Focal calls and marker ascertainment are frozen. Comparative frequency is GT dosage/2.
Legacy/native sources remain traceable and are also plotted separately.
"""
import argparse,csv,hashlib,json,math,os,re,shlex,shutil,subprocess,sys
from pathlib import Path
import numpy as np
import pandas as pd
import pysam
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
CHROMS={'HmRos':'NC_140403.1','Hmel2.5':'Hmel215003o'}
BOUNDS={'HmRos':{'d':(1633041,1635646),'v':(1658079,1659072)},'Hmel2.5':{'d':(1403500,1405500),'v':(1429500,1430500)}}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,obj):p.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def read(p):
    with Path(p).open() as h:return list(csv.DictReader(h,delimiter='\t'))
def tab(p,rows,fields):
    with p.open('w',newline='') as h:
        w=csv.DictWriter(h,delimiter='\t',fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)
def run(args,**kwargs):return subprocess.run([str(x) for x in args],check=True,**kwargs)
def binary(x):
    try:return int(float(x)) if float(x) in (0,1) else None
    except (TypeError,ValueError):return None

def find_source(root):
    valid=[]
    for p in (root/'analysis/yb_previous_calls_new_samples').glob('*/audit.json'):
        a=json.loads(p.read_text());groups={str(r.get('group','')).lower() for r in a.get('new_samples',[])}
        if a.get('status')=='COMPLETE' and a.get('original_SNP_selection_unchanged') is True and {'cythera','vulcanus'}<=groups:valid.append(p.parent)
    if not valid:raise RuntimeError('No completed restored Yb source including cythera/vulcanus')
    return sorted(valid)[-1]

def quality(c,allow_pl=False):
    gt=c.get('GT') or ();dp=c.get('DP') or 0
    if len(gt)!=2 or any(x not in (0,1) for x in gt) or dp<2:return None,'missing_or_low_DP'
    q=c.get('GQ')
    if q is not None:return (sum(gt)/2 if q>=20 else None),'GQ'
    # For biallelic calls with explicit PL evidence. This is a likelihood-derived
    # quality, explicitly distinguished from the caller's reported GQ.
    pl=c.get('PL')
    if not allow_pl or not pl or len(pl)!=3 or any(v is None for v in pl):return None,'GQ_missing'
    x=np.asarray(pl,dtype=float);w=10**(-(x-x.min())/10);probs=w/w.sum();state=sum(gt)
    if probs[state]<probs.max()-1e-10:return None,'GT_PL_disagree'
    error=max(1e-99,float(probs.sum()-probs[state]));q=-10*math.log10(error)
    return (state/2 if q>=20 else None),'PL_derived_Q'

def metadata(root,source):
    rows=read(root/'data/derived/region_panel_metadata.tsv');index={}
    for r in rows:
        ids=set()
        for col in ('Accession','SequenceID','EarthCapeID'):
            text=r.get(col,'');ids.update(re.findall(r'(?:ERR|SRR|DRR|SAMEA|SAMN|ERS)\d+',text))
            if text.strip():ids.add(text.strip())
        for s in ids:index.setdefault(s,[]).append(r)
    cfg=json.loads((source/'run.json').read_text());entries={r['sample']:r for r in cfg.get('entries',[])}
    def get(s):
        rr=index.get(s,[]);group=str(entries.get(s,{}).get('group',''));taxon='; '.join(sorted({' '.join(r.get('Taxon Name','').split()) for r in rr}))
        text=(taxon+' '+group).lower()
        exclusion=''
        if not re.fullmatch(r'(?:ERR|SRR|DRR)\d+',s):exclusion='focal_or_unresolved_identity'
        elif s in {'SRR11243339','SRR11243345','ERR1307105'}:exclusion='previously_excluded_or_cydno_backcross'
        elif re.search(r'piano|bwpk|pk_|gynandromorph|hybrid|backcross|f1|f2',text):exclusion='PK_focal_or_cross'
        states={}
        for k,col in [('d','Yb-dorsal'),('v','Yb-ventral')]:
            vals={binary(r.get(col)) for r in rr}-{None}
            if len(vals)>1:raise RuntimeError('Conflicting phenotype '+s+' '+k)
            states[k]=next(iter(vals)) if vals else None
        if 'cythera' in text or 'vulcanus' in text:states={'d':0,'v':1}
        return dict(sample=s,taxon=taxon or group,d=states['d'],v=states['v'],exclusion=exclusion)
    return get

def fasta_candidates(root):
    for folder in (root/'data',root/'references',root/'external',root/'analysis'):
        if not folder.exists():continue
        for parent,dirs,names in os.walk(folder):
            dirs[:]=[d for d in dirs if d not in {'raw','fastq','logs','tmp','alignment','bams','bam','windows','.git','code','toolenv'} and not d.startswith('figure_code_')]
            for name in names:
                if name.endswith('.fai'):
                    p=Path(parent)/name[:-4]
                    if p.is_file():yield p

def compatible_fasta(path,ref,tables):
    try:
        with pysam.FastaFile(str(path)) as fa:
            if CHROMS[ref] not in fa.references:return False
            for t in tables:
                for r in t.itertuples():
                    if fa.fetch(CHROMS[ref],int(r.pos1)-1,int(r.pos1)).upper()!=r.ref.upper():return False
        return True
    except (OSError,ValueError):return False

def prepare(a):
    root=a.root;out=a.output;out.mkdir(parents=True,exist_ok=False);(out/'logs').mkdir();(out/'calls').mkdir();(out/'frozen').mkdir()
    source=a.source or find_source(root);source_audit=json.loads((source/'audit.json').read_text())
    if source_audit.get('status')!='COMPLETE' or source_audit.get('original_SNP_selection_unchanged') is not True:raise RuntimeError('Source is not a completed restored run')
    getmeta=metadata(root,source);tables={};hashes={}
    for ref in CHROMS:
        for interval in ('d','v'):
            name=f'{ref}_Yb{interval}_unchanged_focal_calls.tsv.gz';p=source/name
            t=pd.read_csv(p,sep='\t').sort_values('pos1')
            if t.empty or t.pos1.duplicated().any():raise RuntimeError('Empty/duplicate frozen marker table '+str(p))
            if set(t.contig)!={CHROMS[ref]}:raise RuntimeError('Unexpected frozen contig')
            if any(len(str(x))!=1 for x in list(t.ref)+list(t.alt)):raise RuntimeError('Frozen non-SNP alleles')
            if not all(x in (r,alt) for x,r,alt in zip(t.rosina_base,t.ref,t.alt)):raise RuntimeError('Invalid Rosina allele polarity')
            shutil.copy2(p,out/'frozen'/name);hashes[name]=sha(p);tables[(ref,interval)]=t
    candidates=list(dict.fromkeys(fasta_candidates(root)));references={}
    for ref,explicit in [('HmRos',a.rosina_reference),('Hmel2.5',a.hmel_reference)]:
        choices=[explicit] if explicit else candidates
        choices=[p for p in choices if compatible_fasta(p,ref,[tables[(ref,i)] for i in ('d','v')])]
        if not choices:raise RuntimeError('No indexed FASTA with exact frozen REF agreement for '+ref+'; supply --rosina-reference / --hmel-reference')
        def fasta_rank(p):
            with pysam.FastaFile(str(p)) as fa:return (-sum(fa.lengths),len(str(p)),str(p))
        references[ref]=str(sorted(choices,key=fasta_rank)[0])
    calls=read(root/'analysis/yb_existing_resources_v1/sample_callability.tsv');bams=read(root/'analysis/yb_existing_resources_v1/bam_candidates.tsv')
    exclusions=[];tasks=[];selected=[]
    allids=sorted({r['sample'] for r in calls+bams})
    for s in allids:
        m=getmeta(s)
        for ref in CHROMS:
            for interval in ('d','v'):
                if m['exclusion'] or m[interval] is None:
                    exclusions.append(dict(sample=s,reference=ref,interval=interval,reason=m['exclusion'] or 'unknown_phenotype',taxon=m['taxon']));continue
                vv=[r for r in calls if r['sample']==s and r['reference']==ref and r['interval']==interval and int(r['callable_DP2_GQ20'])>0]
                # Never silently use legacy aliases as native VCFs.
                vv.sort(key=lambda r:('integration_Yb/' not in r['path'],'rosina_genomewide_joint_recall/' not in r['path'],r['path']))
                if vv:
                    selected.append(dict(**m,reference=ref,interval=interval,source='existing_native_VCF',path=vv[0]['path']));continue
                bb=[r for r in bams if r['sample']==s and r['reference']==ref and float(r['reads_overlapping'] or 0)>0 and '/tmp/' not in r['resolved_path'] and not re.search(r'\.(name|fixmate|position)\.bam$',r['resolved_path'])]
                bb.sort(key=lambda r:(not r['resolved_path'].endswith('.dedup.bam'),r['resolved_path']))
                chosen=None
                with pysam.FastaFile(references[ref]) as fa:dic=dict(zip(fa.references,fa.lengths))
                for b in bb:
                    with pysam.AlignmentFile(b['resolved_path'],'rb') as bam:
                        bd=dict(zip(bam.references,bam.lengths));sm={r.get('SM') for r in bam.header.to_dict().get('RG',[]) if r.get('SM')}
                        if sm!={s} or not bam.has_index() or any(dic.get(c)!=n for c,n in bd.items()):continue
                        chosen=b;break
                if chosen:
                    key=(s,ref);task=next((x for x in tasks if (x['sample'],x['reference'])==key),None)
                    if task is None:
                        task=dict(task=len(tasks),sample=s,reference=ref,bam=chosen['resolved_path']);tasks.append(task)
                    selected.append(dict(**m,reference=ref,interval=interval,source='forced_BAM_call',path=str(out/'calls'/f'{task["task"]}.vcf.gz')))
                else:exclusions.append(dict(sample=s,reference=ref,interval=interval,reason='no_compatible_native_resource',taxon=m['taxon']))
    # Audit the full legacy locus; require exact REF agreement throughout the used dorsal interval.
    legacy=root/'external/krista_genotype_plots_original/Genotype_Plots/Ivory_Hmr15-1627000-1637000_biallelic.vcf'
    checked=0;mismatches=[];target_checked=0;target_mismatches=[]
    dlo,dhi=BOUNDS['HmRos']['d']
    with pysam.FastaFile(references['HmRos']) as fa,pysam.VariantFile(str(legacy)) as vf:
        legacy_names=list(vf.header.samples)
        for r in vf:
            if r.contig!='Hmr15':raise RuntimeError('Unexpected legacy contig')
            checked+=1
            disagrees=fa.fetch(CHROMS['HmRos'],r.pos-1,r.pos-1+len(r.ref)).upper()!=r.ref.upper()
            if disagrees:mismatches.append(r.pos)
            if r.pos<=dhi and r.pos+len(r.ref)-1>=dlo:
                target_checked+=1
                if disagrees:target_mismatches.append(r.pos)
    save(out/'legacy_reference_check.json',dict(path=str(legacy),records_checked=checked,mismatches=len(mismatches),first_mismatches=mismatches[:20],target_interval=[dlo,dhi],target_records_checked=target_checked,target_mismatches=target_mismatches,outside_interval_mismatches_excluded=[x for x in mismatches if x not in target_mismatches],scope='Exact REF agreement required across used dorsal interval; outside-interval discrepancies audited and excluded; no whole-chromosome identity claim'))
    if not target_checked or target_mismatches:raise RuntimeError('Legacy/native coordinate compatibility failed within used dorsal interval; no merge submitted')
    for s in legacy_names:
        m=getmeta(s)
        if m['exclusion'] or m['d'] is None:continue
        # Native source takes precedence at the individual level, never double count.
        if not any(x['sample']==s and x['reference']=='HmRos' and x['interval']=='d' for x in selected):selected.append(dict(**m,reference='HmRos',interval='d',source='legacy_VCF',path=str(legacy)))
    tab(out/'selected_sources.tsv',selected,['sample','taxon','d','v','reference','interval','source','path'])
    tab(out/'excluded_sources.tsv',exclusions,['sample','reference','interval','reason','taxon'])
    for ref in CHROMS:
        t=pd.concat([tables[(ref,i)] for i in ('d','v')]).drop_duplicates('pos1').sort_values('pos1')
        target=out/f'{ref}.alleles.tsv'
        with target.open('w') as h:
            for r in t.itertuples():h.write(f'{r.contig}\t{r.pos1}\t{r.ref},{r.alt}\n')
        run(['bgzip','-f',target]);run(['tabix','-s','1','-b','2','-e','2',str(target)+'.gz'])
        with (out/f'{ref}.regions.bed').open('w') as h:
            for r in t.itertuples():h.write(f'{r.contig}\t{r.pos1-1}\t{r.pos1}\n')
    frozen_script=out/'expand_Yb_existing_controls.py';shutil.copy2(__file__,frozen_script)
    cfg=dict(root=str(root),output=str(out),source=str(source),source_audit_sha=sha(source/'audit.json'),focal_hashes=hashes,metadata_sha256=sha(root/'data/derived/region_panel_metadata.tsv'),references=references,tasks=tasks,legacy=str(legacy),script_sha=sha(frozen_script),selected_rows=len(selected))
    save(out/'run.json',cfg)
    def submit(stage,dependency=None,array=None):
        cmd=['sbatch','--parsable','--partition=nano','--cpus-per-task=1','--mem=8G','--time=00:30:00','--job-name=Yb_'+stage,'--output='+str(out/'logs'/f'{stage}.%A_%a.out'),'--error='+str(out/'logs'/f'{stage}.%A_%a.err')]
        if dependency:cmd+=['--dependency=afterok:'+dependency]
        if array:cmd+=['--array='+array]
        wrapped=shlex.join(['micromamba','run','-n','gynandrocalling','python',str(frozen_script),'--output',str(out),'--stage',stage])
        cmd+=['--wrap='+wrapped];return subprocess.check_output(cmd,text=True).strip().split(';')[0]
    if a.submit:
        c=submit('call',array=f'0-{len(tasks)-1}%12') if tasks else None
        f=submit('report',dependency=c);save(out/'jobs.json',dict(calls=c,report=f))
        print('CALL_ARRAY\t'+str(c));print('REPORT\t'+f)
    print('BAM_TASKS\t'+str(len(tasks)));print('OUTPUT\t'+str(out));print('PDF\t'+str(out/'Yb_expanded_existing_controls_Rosina.pdf'))

def call(cfg,task):
    out=Path(cfg['output']);t=cfg['tasks'][task];ref=t['reference'];dst=out/'calls'/f'{task}.vcf.gz';tmp=out/'calls'/f'{task}.partial.vcf.gz'
    mp=['bcftools','mpileup','-f',cfg['references'][ref],'-R',str(out/f'{ref}.regions.bed'),'-q','30','-Q','20','-d','10000','-a','FORMAT/AD,FORMAT/DP','-Ou',t['bam']]
    ca=['bcftools','call','-m','-A','-C','alleles','-T',str(out/f'{ref}.alleles.tsv.gz'),'-f','GQ','-Oz','-o',str(tmp)]
    p=subprocess.Popen(mp,stdout=subprocess.PIPE)
    try:
        result=subprocess.run(ca,stdin=p.stdout);p.stdout.close();rc=p.wait()
        if rc or result.returncode:raise RuntimeError(f'Calling failed: mpileup={rc}, call={result.returncode}')
    finally:
        if p.poll() is None:p.kill();p.wait()
    with pysam.VariantFile(str(tmp)) as v:
        if list(v.header.samples)!=[t['sample']]:raise RuntimeError('Called sample identity mismatch')
    tmp.rename(dst);run(['bcftools','index','-f',dst]);save(out/'calls'/f'{task}.complete.json',dict(**t,commands=[mp,ca],vcf_sha256=sha(dst)))

def extract(path,samples,table,legacy=False,allow_pl=False):
    # Map only exact REF/ALT pairs. No strand flips or coordinate-only merging.
    markers={int(r.pos1):r for r in table.itertuples()};result={s:np.full(len(table),np.nan) for s in samples};indices={int(p):i for i,p in enumerate(table.pos1)};stats={};seen=set()
    with pysam.VariantFile(str(path)) as vf:
        missing=set(samples)-set(vf.header.samples)
        if missing:raise RuntimeError('VCF sample missing: '+str(missing))
        chrom='Hmr15' if legacy else str(table.contig.iloc[0])
        indexed=any(Path(str(path)+s).is_file() for s in ('.tbi','.csi'))
        records=vf.fetch(chrom,int(table.pos1.min())-1,int(table.pos1.max())) if indexed else iter(vf)
        for rec in records:
            if rec.contig!=chrom or rec.pos not in markers:continue
            r=markers[rec.pos]
            if rec.ref!=r.ref or tuple(rec.alts or ())!=(r.alt,):continue
            if rec.pos in seen:raise RuntimeError('Duplicate matching VCF marker '+str(path)+' '+str(rec.pos))
            seen.add(rec.pos)
            for s in samples:
                value,kind=quality(rec.samples[s],allow_pl);stats[kind]=stats.get(kind,0)+1
                if value is not None:result[s][indices[rec.pos]]=value if r.rosina_base==r.alt else 1-value
    return result,dict(exact_matching_markers=len(seen),quality_evaluations=stats)

def draw(table,matrix,labels,ref,interval,pdf,title):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.patches import ConnectionPatch,Rectangle
    fig=plt.figure(figsize=(7.2,2.65));fig.text(.025,.94,'A' if interval=='d' else 'B',weight='bold',fontsize=10);fig.text(.24,.94,f'Yb({interval})',weight='bold',fontsize=9);fig.text(.94,.94,ref,ha='right',fontsize=6.5)
    ax=fig.add_axes([.24,.28,.69,.36]);locus=fig.add_axes([.24,.76,.69,.07]);n=len(table)
    cmap=LinearSegmentedColormap.from_list('yb',['#D7EDF5','#ECCB63','#662A83']);cmap.set_bad('#C8C8C8')
    im=ax.pcolormesh(np.arange(n+1)-.5,np.arange(5)-.5,np.ma.masked_invalid(matrix),vmin=0,vmax=1,cmap=cmap)
    ax.set_ylim(3.5,-.5);ax.set_xlim(-.5,n-.5);ax.set_yticks(range(4),labels,fontsize=7);ax.set_xticks([]);ax.tick_params(length=0)
    lo,hi=BOUNDS[ref][interval];locus.set_xlim(min(lo-300,table.pos1.min()-20),max(hi+300,table.pos1.max()+20));locus.set_ylim(-.1,1);locus.axhline(0,color='#888',lw=.5);locus.add_patch(Rectangle((lo,.1),hi-lo,.3,color='#777',lw=0));locus.set_yticks([]);locus.ticklabel_format(axis='x',style='plain',useOffset=False);locus.tick_params(labelsize=6,length=2)
    for j,pos in enumerate(table.pos1):fig.add_artist(ConnectionPatch(xyA=(pos,0),coordsA=locus.transData,xyB=(j,-.5),coordsB=ax.transData,color='#9ca2a8',alpha=.22,lw=.25,clip_on=False))
    cb=fig.colorbar(im,cax=fig.add_axes([.68,.12,.24,.025]),orientation='horizontal');cb.set_ticks([0,.5,1]);cb.ax.tick_params(labelsize=6,length=2);cb.outline.set_visible(False);cb.set_label('Rosina allele frequency / focal dosage',fontsize=6,labelpad=1)
    fig.text(.24,.04,f'{n} frozen markers · {title} · grey: unavailable',fontsize=5.8)
    for axis in (ax,locus):
        for spine in axis.spines.values():spine.set_visible(False)
    pdf.savefig(fig);plt.close(fig)

def report(cfg):
    out=Path(cfg['output']);source=Path(cfg['source'])
    if sha(source/'audit.json')!=cfg['source_audit_sha']:raise RuntimeError('Source audit changed')
    if sha(Path(cfg['root'])/'data/derived/region_panel_metadata.tsv')!=cfg['metadata_sha256']:raise RuntimeError('Phenotype metadata changed after preparation')
    for name,digest in cfg['focal_hashes'].items():
        if sha(source/name)!=digest or sha(out/'frozen'/name)!=digest:raise RuntimeError('Frozen focal input changed')
    for t in cfg['tasks']:
        c=json.loads((out/'calls'/f'{t["task"]}.complete.json').read_text())
        if sha(out/'calls'/f'{t["task"]}.vcf.gz')!=c['vcf_sha256']:raise RuntimeError('Task VCF hash mismatch')
    selected=pd.read_csv(out/'selected_sources.tsv',sep='\t');audits=[];members=[];sitecounts=[];sourcechecks=[];pages=[]
    import matplotlib
    matplotlib.use('Agg');matplotlib.rcParams.update({'font.family':'sans-serif','font.size':7,'pdf.fonttype':42})
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(out/'Yb_expanded_existing_controls_Rosina.pdf') as ros,PdfPages(out/'Yb_existing_controls_Hmel25.pdf') as hmel,PdfPages(out/'Yb_call_source_sensitivity.pdf') as sensitivity:
        for ref in CHROMS:
            for interval in ('d','v'):
                table=pd.read_csv(out/'frozen'/f'{ref}_Yb{interval}_unchanged_focal_calls.tsv.gz',sep='\t').sort_values('pos1');ss=selected[(selected.reference==ref)&(selected.interval==interval)];values={}
                for (path,kind),group in ss.groupby(['path','source']):
                    result,check=extract(path,group['sample'].tolist(),table,legacy=kind=='legacy_VCF',allow_pl=True);values.update(result);sourcechecks.append(dict(reference=ref,interval=interval,path=path,source=kind,**check))
                if ss['sample'].duplicated().any():raise RuntimeError('Duplicated biological sample in panel')
                # n counts only individuals with at least one usable frozen-marker call.
                active=ss[ss['sample'].map(lambda s:np.isfinite(values[s]).any())].copy()
                for r in ss.to_dict('records'):members.append(dict(**r,markers_called=int(np.isfinite(values[r['sample']]).sum()),included=bool(np.isfinite(values[r['sample']]).any())))
                def matrix_for(group):
                    rows=[];ns=[];ids_by_state=[]
                    for state in (1,0):
                        ids=group.loc[group[interval]==state,'sample'].tolist();ids_by_state.append(ids);ns.append(len(ids))
                        mat=np.array([values[s] for s in ids]) if ids else np.empty((0,len(table)));count=np.isfinite(mat).sum(axis=0)
                        mean=np.divide(np.nansum(mat,axis=0),count,out=np.full(len(table),np.nan),where=count>0);rows.append(mean)
                    return np.vstack([rows[0],table.left_rosina_dosage,table.right_rosina_dosage,rows[1]]),ns,ids_by_state
                matrix,ns,ids=matrix_for(active)
                np.testing.assert_equal(matrix[1:3],table[['left_rosina_dosage','right_rosina_dosage']].to_numpy().T)
                labels=[f'Yb({interval})+ (n = {ns[0]})','TNLeft','TNRight',f'Yb({interval})− (n = {ns[1]})']
                draw(table,matrix,labels,ref,interval,ros if ref=='HmRos' else hmel,'non-PK genotype frequencies')
                pd.DataFrame(matrix,index=labels,columns=table.pos1).to_csv(out/f'{ref}_Yb{interval}_plotted_values.tsv',sep='\t',na_rep='NA')
                for state,individuals in zip((1,0),ids):
                    for j,pos in enumerate(table.pos1):sitecounts.append(dict(reference=ref,interval=interval,pos1=pos,phenotype=state,n_individuals_called=sum(np.isfinite(values[s][j]) for s in individuals)))
                for kind,group in active.groupby('source'):
                    mm,nn,_=matrix_for(group);ll=[f'Yb({interval})+ (n = {nn[0]})','TNLeft','TNRight',f'Yb({interval})− (n = {nn[1]})'];draw(table,mm,ll,ref,interval,sensitivity,kind)
                # Explicit per-marker native-vs-legacy contrasts, without claiming batch equivalence.
                for state in (1,0):
                    for kind,group in active[active[interval]==state].groupby('source'):
                        for j,pos in enumerate(table.pos1):
                            vv=[values[s][j] for s in group['sample'] if np.isfinite(values[s][j])]
                            audits.append(dict(reference=ref,interval=interval,pos1=pos,phenotype=state,source=kind,n=len(vv),frequency=float(np.mean(vv)) if vv else None))
                pages.append(dict(reference=ref,interval=interval,markers=len(table),positive_n=ns[0],negative_n=ns[1],focal_identical=True))
                print(f'PANEL\t{ref}\tYb({interval})\tpositive={ns[0]}\tnegative={ns[1]}')
    pd.DataFrame(members).to_csv(out/'control_membership.tsv',sep='\t',index=False,na_rep='NA');pd.DataFrame(sitecounts).to_csv(out/'per_site_callable_individuals.tsv',sep='\t',index=False);pd.DataFrame(audits).to_csv(out/'source_frequency_sensitivity.tsv',sep='\t',index=False,na_rep='NA');save(out/'source_harmonization.json',sourcechecks)
    save(out/'complete.json',dict(status='COMPLETE',source=str(source),focal_calls_unchanged=True,marker_selection_unchanged=True,PK_in_controls=False,comparative_statistic='Equal individual genotype dosage/2, DP>=2 and reported GQ>=20; all sources may use explicitly recorded PL-derived likelihood Q>=20 when GQ absent',legacy_and_native_calling_methods_identical=False,legacy_reference_check=json.loads((out/'legacy_reference_check.json').read_text()),source_sensitivity_pdf=str(out/'Yb_call_source_sensitivity.pdf'),pages=pages,focal_hashes=cfg['focal_hashes'],phenotype_overrides={'cythera':{'d':0,'v':1},'vulcanus':{'d':0,'v':1}},references=cfg['references']))
    print('PDF\t'+str(out/'Yb_expanded_existing_controls_Rosina.pdf'))

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=ROOT);p.add_argument('--source',type=Path);p.add_argument('--output',type=Path,required=True);p.add_argument('--rosina-reference',type=Path);p.add_argument('--hmel-reference',type=Path);p.add_argument('--submit',action='store_true');p.add_argument('--stage',choices=['prepare','call','report'],default='prepare');p.add_argument('--task',type=int);a=p.parse_args()
    if a.stage=='prepare':prepare(a)
    else:
        cfg=json.loads((a.output/'run.json').read_text())
        if sha(__file__)!=cfg['script_sha']:raise RuntimeError('Frozen worker changed')
        if a.stage=='call':call(cfg,a.task if a.task is not None else int(os.environ['SLURM_ARRAY_TASK_ID']))
        else:report(cfg)
if __name__=='__main__':main()
