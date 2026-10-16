#!/usr/bin/env python3
"""Add display samples to frozen historical Yb SNPs and focal calls.
No SNP selection, focal recounting, focal inference, or parameter refitting.
"""
import argparse
import hashlib
import importlib.util
import json
import shlex
import shutil
import subprocess
import sys
from collections import Counter,defaultdict
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import logsumexp

ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
COLORS=['#D7EDF5','#ECCB63','#662A83']
INTERVALS={'Hmel2.5':[(1403500,1405500),(1429500,1430500)],'HmRos':[(1633041,1635646),(1658079,1659072)]}
PRIMERS={'Hmel2.5':(1404112,1405075),'HmRos':(1633657,1635225)}
PARAMS=['estimated_c_left_from_right','estimated_c_right_from_left','estimated_error_rate','ALT_mapping_bias_log_odds','beta_binomial_concentration']

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,obj):Path(p).write_text(json.dumps(obj,indent=2)+'\n')
def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def parse_baseline(log):
    hits=[line.split('\t',1)[1] for line in Path(log).read_text().splitlines() if line.startswith('BUNDLE_PDF\t')]
    if len(hits)!=1:raise RuntimeError('Need one completed restoration BUNDLE_PDF in '+str(log))
    pdf=Path(hits[0])
    if not pdf.is_file():raise RuntimeError('Restored PDF missing')
    return pdf.parent

def mean_available(values):
    arr=np.asarray(values,float);n=np.isfinite(arr).sum(axis=0)
    return np.divide(np.nansum(arr,axis=0),n,out=np.full(arr.shape[1],np.nan),where=n>0)

def pk_dosage(model,pars,nr,na,ros_is_ref):
    if nr+na<model.MIN_CALL_DP:return np.nan,np.nan
    means=model.map_fraction(model.DOSAGE,pars[2],pars[3])
    ll=model.beta_binomial_logpmf(na,nr+na,means,pars[4])
    post=np.exp(ll-logsumexp(ll));best=int(np.argmax(post));p=float(post[best])
    value=float(model.DOSAGE[best]) if p>=model.POSTERIOR_CUTOFF else np.nan
    return (1-value if ros_is_ref else value),p

def prepare(a):
    import pysam,matplotlib,scipy
    baseline=parse_baseline(a.restore_log)
    shared=json.loads((a.shared_run/'run.json').read_text())
    if shared.get('no_hmel25') is not True:raise RuntimeError('Expected Rosina BAM manifest')
    entries=[r for r in shared['samples'] if r['sample'] in shared['comparative']]
    if len(entries)!=len(shared['comparative']):raise RuntimeError('Incomplete comparative manifest')
    if len({r['sample'] for r in entries})!=len(entries):raise RuntimeError('Duplicate samples')
    for group in ('cythera','vulcanus'):
        if not any(str(r['group']).lower()==group for r in entries):
            raise RuntimeError('Missing requested new race: '+group)
    # Explicit shared reference dictionary equality; never transfer Rosina coordinates to Hmel.
    focal=ROOT/'analysis/yb_local_both_references/20261012T194933Z_7fa9de96'
    with pysam.FastaFile(str(focal/'HmRos/reference/genome.fa')) as f:
        expected=list(zip(f.references,f.lengths))
    for r in entries:
        with pysam.AlignmentFile(r['bam'],'rb') as bam:
            bam.check_index()
            if list(zip(bam.references,bam.lengths))!=expected:raise RuntimeError('Wrong reference: '+r['sample'])
            if {rg.get('SM') for rg in bam.header.to_dict().get('RG',[])}!={r['sample']}:
                raise RuntimeError('Wrong BAM sample identity: '+r['sample'])
    out=a.output_dir or ROOT/'analysis/yb_previous_calls_new_samples'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    out=out.absolute();out.mkdir(parents=True,exist_ok=False);(out/'logs').mkdir()
    frozen=out/'baseline';frozen.mkdir()
    inputs=[]
    for ref in INTERVALS:
        for interval in ('d','v'):
            for suffix in ('genotype_pairs.tsv.gz','PK_genotypes.tsv.gz','extra_race_frequencies.tsv.gz'):
                p=baseline/f'{ref}_Yb{interval}_{suffix}'
                if suffix.startswith('genotype') and not p.is_file():raise RuntimeError('Missing frozen focal calls: '+str(p))
                if p.is_file():shutil.copy2(p,frozen/p.name);inputs.append(dict(source=str(p),frozen=p.name,sha256=sha(p)))
            p=a.counts_run/f'{ref}_Yb{interval}_SNP_counts.tsv.gz'
            shutil.copy2(p,frozen/p.name);inputs.append(dict(source=str(p),frozen=p.name,sha256=sha(p)))
    for name in ('fit_summary.tsv','original_caller.py','audit.json'):
        p=baseline/name;shutil.copy2(p,frozen/name);inputs.append(dict(source=str(p),frozen=name,sha256=sha(p)))
    # Obtain the original fragment counter from the pinned restoration snapshot.
    restoration=Path(a.restore_log).parent.parent
    helper=restoration/'code/scripts/yb/assign_Yb_local_SNP_fragments.py'
    shutil.copy2(helper,out/'original_fragment_counter.py')
    shutil.copy2(__file__,out/Path(__file__).name)
    cfg=dict(output=str(out),baseline_source=str(baseline),shared_run=str(a.shared_run),
        entries=entries,inputs=inputs,helper_sha256=sha(out/'original_fragment_counter.py'),
        worker_sha256=sha(out/Path(__file__).name),focal_calls_reused=True,
        SNP_selection_changed=False,parameters_refitted=False,
        new_samples_reference='HmRos only; no equivalent new Hmel2.5 BAMs supplied')
    dump(out/'run.json',cfg)
    command=[sys.executable,str(out/Path(__file__).name),'--work',str(out)]
    worker=out/'worker.sh';worker.write_text('#!/usr/bin/env bash\nset -euo pipefail\nexec '+shlex.join(command)+'\n')
    job=subprocess.check_output(['sbatch','--parsable','--partition=nano','--time=00:30:00',
        '--cpus-per-task=2','--mem=12G','--job-name=YbAddSamples',
        '--output='+str(out/'logs/add_samples.%j.out'),'--error='+str(out/'logs/add_samples.%j.err'),str(worker)],text=True).strip().split(';')[0]
    dump(out/'jobs.json',dict(add_samples=job))
    print('JOBID\t'+job);print('PDF\t'+str(out/'Yb_previous_calls_with_vulcanus_cythera_and_all_ready_PK.pdf'))
    print('LOG\t'+str(out/'logs'/f'add_samples.{job}.out'))

def work(out):
    import matplotlib
    matplotlib.use('Agg');matplotlib.rcParams.update({'pdf.fonttype':42,'ps.fonttype':42})
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.patches import ConnectionPatch
    cfg=json.loads((out/'run.json').read_text());base=out/'baseline'
    if sha(Path(__file__))!=cfg['worker_sha256']:raise RuntimeError('Frozen display adapter changed')
    for r in cfg['inputs']:
        if sha(base/r['frozen'])!=r['sha256']:raise RuntimeError('Frozen baseline changed')
    if sha(out/'original_fragment_counter.py')!=cfg['helper_sha256']:raise RuntimeError('Original counter changed')
    model=module(base/'original_caller.py','original_model')
    observations=module(out/'original_fragment_counter.py','original_counter').observations
    fit=pd.read_csv(base/'fit_summary.tsv',sep='\t').set_index('statistic')['value']
    pars=[float(fit[k]) for k in PARAMS]
    cmap=LinearSegmentedColormap.from_list('original',COLORS);cmap.set_bad('#C8C8C8')
    pages=[];membership=[];audit_samples=[]
    pdfpath=out/'Yb_previous_calls_with_vulcanus_cythera_and_all_ready_PK.pdf'
    with PdfPages(pdfpath) as pdf:
        for ref in ('Hmel2.5','HmRos'):
            for interval_index,interval in enumerate(('d','v')):
                frozen_calls=base/f'{ref}_Yb{interval}_genotype_pairs.tsv.gz'
                table=pd.read_csv(frozen_calls,sep='\t').sort_values('pos1')
                if table.empty:raise RuntimeError('No original SNP records')
                if table.pos1.duplicated().any():raise RuntimeError('Duplicate original marker')
                raw=pd.read_csv(base/f'{ref}_Yb{interval}_SNP_counts.tsv.gz',sep='\t')
                positions=table.pos1.tolist();records=table.to_dict('records');chrom=records[0]['contig']
                labels=['rosina controls (original)','melpomene controls (original)','TNLeft corrected (unchanged)','TNRight corrected (unchanged)']
                values=[];original_control_ids=set()
                for race in ('rosina','melpomene'):
                    selected=raw.loc[(raw.role=='control') & raw.race.eq(race)]
                    original_control_ids.update(selected['sample'])
                    values.append(selected.groupby('pos1').rosina_frequency.mean().reindex(positions).to_numpy())
                original_focal=table[['left_rosina_dosage','right_rosina_dosage']].to_numpy().T
                values.extend(original_focal)
                seen=set(original_control_ids)
                # Preserve original PK calls byte-for-byte in value, then add missing PKs.
                pkfile=base/f'{ref}_Yb{interval}_PK_genotypes.tsv.gz'
                pk_values={}
                if pkfile.is_file():
                    pks=pd.read_csv(pkfile,sep='\t')
                    for sample,t in pks.groupby('sample'):
                        pk_values[sample]=t.set_index('pos1').accepted_rosina_dosage.reindex(positions).to_numpy()
                racefile=base/f'{ref}_Yb{interval}_extra_race_frequencies.tsv.gz'
                extra={}
                if racefile.is_file():
                    old=pd.read_csv(racefile,sep='\t')
                    for (race,sample),t in old.groupby(['race','sample']):
                        extra.setdefault(race,[]).append((sample,t.set_index('pos1').frequency.reindex(positions).to_numpy()))
                        seen.add(sample)
                new_rows=[];pk_rows=[]
                entries=cfg['entries'] if ref=='HmRos' else []
                sites={p-1 for p in positions}
                for r in entries:
                    sample=r['sample'];group=str(r['group']);is_pk=group.startswith('PK_')
                    if sample in seen:continue
                    if is_pk and sample in pk_values:continue
                    counts=defaultdict(Counter)
                    for fragment in observations(r['bam'],chrom,min(sites),max(sites)+1,sites).values():
                        for pos,b in fragment.items():counts[pos][b]+=1
                    vv=[]
                    for rec in records:
                        c=counts[rec['pos1']-1];nr=c[rec['ref']];na=c[rec['alt']]
                        ros=c[rec['rosina_base']];total=nr+na
                        if is_pk:
                            val,post=pk_dosage(model,pars,nr,na,rec['rosina_base']==rec['ref'])
                            pk_rows.append(dict(sample=sample,pos1=rec['pos1'],ref_count=nr,alt_count=na,
                                other_count=sum(c.values())-total,posterior=post,rosina_dosage=val))
                        else:val=ros/total if total>=8 else np.nan
                        vv.append(val)
                        new_rows.append(dict(sample=sample,group=group,pos1=rec['pos1'],ref=rec['ref'],alt=rec['alt'],
                            ref_count=nr,alt_count=na,other_count=sum(c.values())-total,display_value=val))
                    if is_pk:pk_values[sample]=np.array(vv)
                    else:extra.setdefault(group,[]).append((sample,np.array(vv)))
                    audit_samples.append(dict(reference=ref,interval=interval,sample=sample,group=group,bam=r['bam'],status='added at original markers'))
                entries_by_id={r['sample']:r for r in cfg['entries']}
                if ref=='HmRos':
                    expected_PK={r['sample'] for r in cfg['entries'] if str(r['group']).startswith('PK_')}
                    if not expected_PK.issubset(pk_values):raise RuntimeError('Retained PKs omitted')
                for sample,v in sorted(pk_values.items()):
                    r=entries_by_id.get(sample,{})
                    label=f'{r.get("label",sample)} | {r.get("group","PK").replace("PK_","PK ")} | d+/v+ | {sample}'
                    labels.append(label);values.append(v)
                    membership.append(dict(reference=ref,interval=interval,label=label,samples=[sample],
                        kind='PK genotype dosage',phenotype_source='Joe Hanly: all retained PKs Yb+, dorsal and ventral',
                        calls='frozen previous PK calls' if sample in (set(pks['sample']) if pkfile.is_file() else set()) else 'same original unmixed posterior95 model at frozen markers'))
                for group,items in sorted(extra.items()):
                    phenotype={'cythera':'d−/v−','vulcanus':'d−/v+'}.get(group.lower(),'')
                    label=f'{group} {phenotype} | n={len(items)}'
                    labels.append(label);values.append(mean_available([v for sample,v in items]))
                    membership.append(dict(reference=ref,interval=interval,label=label,samples=[sample for sample,v in items],
                        kind='mean individual read frequency',minimum_individual_biallelic_depth=8,
                        phenotype_source='Joe Hanly 2026-10-16' if phenotype else 'original race label'))
                mat=np.asarray(values)
                # This is a hard invariant, not an expectation about biological states.
                np.testing.assert_equal(mat[2:4],original_focal)
                pd.DataFrame(mat,index=labels,columns=positions).to_csv(out/f'{ref}_Yb{interval}_plotted_values.tsv',sep='\t',na_rep='NA')
                table.to_csv(out/f'{ref}_Yb{interval}_unchanged_focal_calls.tsv.gz',sep='\t',index=False)
                pd.DataFrame(new_rows).to_csv(out/f'{ref}_Yb{interval}_new_sample_counts.tsv.gz',sep='\t',index=False)
                if pk_rows:pd.DataFrame(pk_rows).to_csv(out/f'{ref}_Yb{interval}_new_PK_calls.tsv.gz',sep='\t',index=False)
                height=len(labels);fig=plt.figure(figsize=(15,max(7,3.8+height*.29)))
                heat=fig.add_axes([.37,.30,.54,.55]);locus=fig.add_axes([.37,.13,.54,.09])
                im=heat.imshow(mat,cmap=cmap,vmin=0,vmax=1,aspect='auto',interpolation='nearest')
                heat.set_yticks(range(height),labels,fontsize=7);heat.set_xticks([])
                heat.axhline(3.5,color='white',lw=2)
                cb=fig.colorbar(im,cax=fig.add_axes([.93,.34,.012,.44]));cb.set_ticks([0,.5,1])
                cb.set_label('Rosina-associated frequency / genotype dosage',fontsize=8)
                lo,hi=INTERVALS[ref][interval_index];locus.axvspan(lo,hi,color='grey',alpha=.25)
                locus.axhline(0,color='grey',lw=.8)
                if interval=='d':
                    aa,bb=PRIMERS[ref];locus.plot([aa,bb],[.3,.3],color=COLORS[2],lw=3)
                    locus.text((aa+bb)/2,.4,'Luca primer span',ha='center',fontsize=8)
                for j,p in enumerate(positions):
                    fig.add_artist(ConnectionPatch(xyA=(j,height-.5),coordsA=heat.transData,
                        xyB=(p,.5),coordsB=locus.transData,color='#89939d',alpha=.18,lw=.35))
                locus.set_xlim(min(lo-1000,min(positions)-100),max(hi+1000,max(positions)+100))
                locus.set_ylim(-.1,.65);locus.set_yticks([]);locus.ticklabel_format(axis='x',style='plain',useOffset=False)
                locus.set_xlabel(chrom+' position (bp)',fontsize=8)
                fig.suptitle(f'{ref} Yb({interval}) — original {len(table)} SNPs and unchanged Left/Right calls',fontsize=13)
                note='New cythera/vulcanus/PK BAMs are Rosina-only; this Hmel2.5 page retains the earlier display cohort.\n' if ref=='Hmel2.5' else ''
                fig.text(.37,.025,note+'Original rosina/melpomene marker selection, focal counts and mixing fit preserved.\n'
                    'Added races: equal-individual mean read frequency, depth≥8; PK: posterior95 genotype dosage.\n'
                    'Grey = missing. New samples do not enter marker selection or focal inference. No phasing implied.',fontsize=7)
                pdf.savefig(fig);plt.close(fig)
                pages.append(dict(reference=ref,interval=interval,original_markers=len(table),rows=height,
                    focal_values_identical=True,new_samples_available=ref=='HmRos',PK_included=len(pk_values)))
    dump(out/'membership.json',membership)
    dump(out/'audit.json',dict(status='COMPLETE',baseline_source=cfg['baseline_source'],
        original_inputs=cfg['inputs'],parameters=dict(zip(PARAMS,pars)),parameters_refitted=False,
        original_SNP_selection_unchanged=True,focal_calls_reused_without_inference=True,
        original_counter_sha256=cfg['helper_sha256'],new_samples=audit_samples,pages=pages,
        new_samples_reference=cfg['new_samples_reference'],pdf=str(pdfpath),pdf_sha256=sha(pdfpath)))
    print('BUNDLE_PDF\t'+str(pdfpath));print('AUDIT\t'+str(out/'audit.json'))
    for page in pages:print('PANEL\t'+json.dumps(page))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--submit',action='store_true');p.add_argument('--restore-log',type=Path)
    p.add_argument('--counts-run',type=Path);p.add_argument('--shared-run',type=Path)
    p.add_argument('--output-dir',type=Path);p.add_argument('--work',type=Path)
    a=p.parse_args()
    if a.submit:
        if not all((a.restore_log,a.counts_run,a.shared_run)):p.error('Need --restore-log, --counts-run and --shared-run')
        for path in (a.restore_log,a.counts_run,a.shared_run):
            if not path.is_absolute():p.error('Use absolute paths')
        prepare(a)
    elif a.work:work(a.work)
    else:p.error('Use --submit or --work')

if __name__=='__main__':main()
