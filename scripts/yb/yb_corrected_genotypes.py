#!/usr/bin/env python3
"""Yb display adapter: archived optix likelihood, fixed saved fit, no refitting."""
import argparse,hashlib,importlib.util,json,shutil
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import logsumexp
REFS=['Hmel2.5','HmRos']
REGIONS={'Hmel2.5':('Hmel215003o',[(1403500,1405500),(1429500,1430500)]),'HmRos':('NC_140403.1',[(1633041,1635646),(1658079,1659072)])}
PRIMERS={'Hmel2.5':(1404112,1405075),'HmRos':(1633657,1635225)}
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

PKS=['SRR11243344','SRR11243353','SRR11243364','SRR11243375','SRR11243376']
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
def available_pk():
    meta_path=ROOT/'data/derived/region_panel_metadata.tsv'
    meta=pd.read_csv(meta_path,sep='\t',dtype=str).fillna('')
    candidates=sorted((ROOT/'analysis/yb_extra_controls').glob('*/run.json'),reverse=True)
    for config in candidates:
        jobs=config.parent/'jobs.json'
        if jobs.is_file() and str(json.loads(jobs.read_text()).get('align'))=='138462316':
            run=config.parent;break
    else:raise RuntimeError('Cannot locate original PK alignment run 138462316')
    entries={};labels={}
    for sample in PKS:
        rows=meta.loc[meta['Accession'].eq(sample)]
        if rows.empty:raise RuntimeError('Missing PK phenotype metadata: '+sample)
        phen=rows[['Taxon Name','Yb-dorsal','Yb-ventral']].drop_duplicates()
        if len(phen)!=1:raise RuntimeError('Ambiguous PK phenotype metadata: '+sample)
        row=rows.iloc[0];text=' '.join(row.astype(str)).lower()
        import re
        shade='dark' if re.search(r'\b(?:pk|piano keys)[ ,_-]*dark\b',text) else ('pale' if re.search(r'\b(?:pk|piano keys)[ ,_-]*pale\b',text) else 'shade unverified')
        labels[sample]=sample+' | '+shade+' | Yb d='+row['Yb-dorsal']+', v='+row['Yb-ventral']
        entries[sample]={ref:run/ref/'coverage'/(sample+'.dedup.bam') for ref in REFS}
    return run,entries,labels,meta_path
def pk_matrix(entries,labels,ref,records,model,pars,out,interval):
    from assign_Yb_local_SNP_fragments import observations
    from collections import Counter,defaultdict
    matrices=[];names=[];audit=[];rows=[]
    for sample,paths in entries.items():
        bam=paths[ref];complete=bam.with_name(sample+'.complete.json')
        if not bam.is_file() or not Path(str(bam)+'.bai').is_file() or not complete.is_file():
            audit.append({'sample':sample,'reference':ref,'status':'alignment not complete; omitted'});continue
        sites={r['pos1']-1 for r in records};contig=records[0]['contig']
        counts=defaultdict(Counter)
        for calls in observations(bam,contig,min(sites),max(sites)+1,sites).values():
            for pos,base in calls.items():counts[pos][base]+=1
        vals=[]
        for r in records:
            cnt=counts[r['pos1']-1];alt=cnt[r['alt']];n=alt+cnt[r['ref']];other=sum(cnt.values())-n;value=np.nan;posterior=np.nan;gt='missing'
            if n>=model.MIN_CALL_DP:
                expected=model.map_fraction(model.DOSAGE,pars[2],pars[3])
                ll=model.beta_binomial_logpmf(alt,n,expected,pars[4]);post=np.exp(ll-logsumexp(ll));best=int(np.argmax(post));posterior=float(post[best]);gt=str(model.GT[best])
                if posterior>=model.POSTERIOR_CUTOFF:
                    value=float(model.DOSAGE[best]);value=1-value if r['rosina_base']==r['ref'] else value
            vals.append(value);rows.append(dict(sample=sample,reference=ref,interval=interval,pos1=r['pos1'],ref=r['ref'],alt=r['alt'],ref_count=cnt[r['ref']],alt_count=alt,other_count=other,posterior=posterior,best_GT=gt,accepted_rosina_dosage=value,phenotype_label=labels[sample]))
        matrices.append(vals);names.append(labels[sample]);audit.append({'sample':sample,'reference':ref,'bam':str(bam),'complete_marker_sha256':digest(complete),'status':'included'})
    pd.DataFrame(rows).to_csv(out/f'{ref}_Yb{interval}_PK_genotypes.tsv.gz',sep='\t',index=False)
    return matrices,names,audit


def race_matrix(source,ref,records,out,interval):
    from assign_Yb_local_SNP_fragments import observations
    from collections import Counter,defaultdict
    cfg=json.loads((source/'run.json').read_text())
    groups=[('melpomene','amaryllis','H. m. amaryllis | Yb+'),('timareta','thelxinoe','H. t. thelxinoe | Yb+'),('timareta','florencia','H. t. florencia | Yb-')]
    values=[];labels=[];audit=[];allrows=[]
    sites={r['pos1']-1 for r in records};chrom=records[0]['contig']
    for species,race,label in groups:
        entries=[e for e in cfg['tasks'] if e['reference']['name']==ref and e.get('species')==species and e.get('subspecies')==race and e['role']=='control']
        if len(entries)!=2:raise RuntimeError('Expected two existing controls for '+race)
        individuals=[]
        for e in entries:
            bam=Path(e['bam'])
            if not bam.is_file() or not Path(str(bam)+'.bai').is_file():raise RuntimeError('Missing race BAM '+str(bam))
            counts=defaultdict(Counter)
            for calls in observations(bam,chrom,min(sites),max(sites)+1,sites).values():
                for pos,b in calls.items():counts[pos][b]+=1
            vals=[]
            for r in records:
                cnt=counts[r['pos1']-1];ros=cnt[r['rosina_base']];otherallele=r['alt'] if r['rosina_base']==r['ref'] else r['ref'];other=cnt[otherallele];n=ros+other
                freq=ros/n if n>=8 else np.nan
                vals.append(freq);allrows.append(dict(reference=ref,interval=interval,species=species,race=race,sample=e['sample'],pos1=r['pos1'],rosina_count=ros,other_selected_allele_count=other,third_allele_count=sum(cnt.values())-n,frequency=freq))
            individuals.append(vals)
            audit.append({'reference':ref,'species':species,'race':race,'sample':e['sample'],'bam':str(bam)})
        arr=np.asarray(individuals,dtype=float);valid=np.isfinite(arr);denom=valid.sum(axis=0);mean=np.divide(np.nansum(arr,axis=0),denom,out=np.full(arr.shape[1],np.nan),where=denom>0)
        values.append(mean);labels.append(label+' (n=2)')
    pd.DataFrame(allrows).to_csv(out/f'{ref}_Yb{interval}_extra_race_frequencies.tsv.gz',sep='\t',index=False)
    return values,labels,audit

def main():
    p=argparse.ArgumentParser();p.add_argument('--counts-run',required=True,type=Path);p.add_argument('--fit-summary',required=True,type=Path);p.add_argument('--focal-run',required=True,type=Path);p.add_argument('--add-race-strips',action='store_true');p.add_argument('--add-ready-pks',action='store_true');a=p.parse_args()
    import pysam
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.patches import ConnectionPatch
    caller=Path(__file__).resolve().parents[1]/'reciprocal_mixing/134_fit_reciprocal_L_R_mixture_genotypes.py'
    for path in (a.counts_run,a.fit_summary,a.focal_run,caller):
        if not path.is_absolute() or not path.exists():raise RuntimeError('Missing absolute input: '+str(path))
    spec=importlib.util.spec_from_file_location('archived_caller',caller);model=importlib.util.module_from_spec(spec);spec.loader.exec_module(model)
    summary=pd.read_csv(a.fit_summary,sep='\t').set_index('statistic').value
    keys=['estimated_c_left_from_right','estimated_c_right_from_left','estimated_error_rate','ALT_mapping_bias_log_odds','beta_binomial_concentration'];pars=[float(summary[k]) for k in keys]
    if not all(np.isfinite(pars)) or not all(0<=x<.5 for x in pars[:2]) or not 0<pars[2]<.5 or pars[4]<2:raise RuntimeError('Invalid fit parameters')
    from datetime import datetime,timezone
    out=a.counts_run/'corrected_genotypes'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ');out.mkdir(parents=True)
    shutil.copy2(caller,out/'original_caller.py');shutil.copy2(a.fit_summary,out/'fit_summary.tsv');shutil.copy2(__file__,out/Path(__file__).name)
    cmap=LinearSegmentedColormap.from_list('optix',[(0,'#D7EDF5'),(.5,'#ECCB63'),(1,'#662A83')]);cmap.set_bad('#C8C8C8');plt.rcParams['pdf.fonttype']=42
    pkentries={};pklabels={}
    if a.add_ready_pks:
        pkrun,pkentries,pklabels,metapath=available_pk()
    audit={'caller':str(caller),'caller_sha256':digest(caller),'fit_summary':str(a.fit_summary),'fit_summary_sha256':digest(a.fit_summary),'parameters':dict(zip(keys,pars)),'parameters_refitted':False,'fit_transfer':'saved optix nuisance parameters applied to new Yb counts on both assemblies; not a Yb/reference-specific fit','prior':'equal over nine diploid genotype pairs','posterior_cutoff':.95,'min_count_depth_per_merged_side':2,'aggregation':'T2 and T3 counts summed within side','count_method':'existing MAPQ30/baseQ30 fragment counts; overlapping mates counted once per site','count_method_differs_from_original_VCF_AD':True,'selection_uses_focal':False,'other_alleles':'retained in table; excluded from biallelic likelihood denominator','intervals':'legacy annotations; Luca primer span an anchor, not exact functional boundaries','inputs':[],'panels':[]}
    if a.add_ready_pks:audit['PK']={'run':str(pkrun),'metadata':str(metapath),'metadata_sha256':digest(metapath),'phenotypes':pklabels,'selection_changed':False,'inference':'unmixed diploid control; original error/bias/concentration transferred; posterior95'}
    pdfpath=out/('Yb_d_v_posterior95_corrected_genotypes_extra_race_strips.pdf' if a.add_race_strips else ('Yb_d_v_posterior95_corrected_genotypes_with_ready_PK_phenotypes.pdf' if a.add_ready_pks else 'Yb_d_v_posterior95_corrected_genotypes_both_references.pdf'))
    with PdfPages(pdfpath) as pdf:
        for ref in REFS:
            fasta=a.focal_run/ref/'reference/genome.fa'
            with pysam.FastaFile(str(fasta)) as genome:
                for idx,interval in enumerate(('d','v')):
                    path=a.counts_run/f'{ref}_Yb{interval}_SNP_counts.tsv.gz';audit['inputs'].append({'path':str(path),'sha256':digest(path),'reference_fasta':str(fasta)})
                    df=pd.read_csv(path,sep='\t');records=[];ctrl=[]
                    for pos,group in df.groupby('pos1',sort=True):
                        first=group.iloc[0];rb=genome.fetch(first.contig,int(pos)-1,int(pos)).upper();ros=first.rosina_base;mel=first.melpomene_base
                        if rb not in (ros,mel):
                            audit.setdefault('excluded_reference_mismatch',[]).append({'reference':ref,'interval':interval,'pos1':int(pos),'reference_base':rb});continue
                        totals=[]
                        for names in (['TN2L_S1','TN3L_S3'],['TN2R_S2','TN3R_S4']):
                            g=group.loc[group['sample'].isin(names)]
                            if len(g)!=2:raise RuntimeError('Missing/duplicate focal counts')
                            totals.append((int(g.rosina_count.sum()),int(g.melpomene_count.sum()),int(g.other_count.sum())))
                        l,r=totals;nl=l[0]+l[1];nr=r[0]+r[1];alt=mel if rb==ros else ros;la=l[1] if rb==ros else l[0];ra=r[1] if rb==ros else r[0]
                        post=np.full(9,np.nan);best=0;confidence=np.nan;ld=rd=np.nan
                        if min(nl,nr)>=model.MIN_CALL_DP:
                            ll=model.state_log_likelihoods(np.array([la]),np.array([nl]),np.array([ra]),np.array([nr]),*pars)[0];post=np.exp(ll-logsumexp(ll));best=int(np.argmax(post));confidence=float(post[best])
                            if confidence>=model.POSTERIOR_CUTOFF:
                                ld=float(model.GL[best]);rd=float(model.GR[best])
                                if ros==rb:ld=1-ld;rd=1-rd
                        rec=dict(reference=ref,interval=interval,contig=first.contig,pos1=int(pos),ref=rb,alt=alt,rosina_base=ros,left_alt=la,left_total=nl,right_alt=ra,right_total=nr,left_other=l[2],right_other=r[2],best_left_GT=str(model.GTL[best]) if np.isfinite(confidence) else 'missing',best_right_GT=str(model.GTR[best]) if np.isfinite(confidence) else 'missing',pair_posterior=confidence,left_rosina_dosage=ld,right_rosina_dosage=rd)
                        rec.update({'posterior_'+str(i):float(v) for i,v in enumerate(post)});records.append(rec)
                        ctrl.append([group.loc[(group.role=='control')&(group.race==race)].rosina_frequency.mean() for race in ('rosina','melpomene')])
                    table=pd.DataFrame(records);table.to_csv(out/f'{ref}_Yb{interval}_genotype_pairs.tsv.gz',sep='\t',index=False)
                    pkvals=[];pknames=[]
                    if a.add_ready_pks and records:
                        pkvals,pknames,pkaudit=pk_matrix(pkentries,pklabels,ref,records,model,pars,out,interval);audit.setdefault('PK_alignments',[]).extend(pkaudit)
                    racevals=[];racenames=[]
                    if a.add_race_strips and records:
                        source=a.counts_run.parent.parent
                        racevals,racenames,raceaudit=race_matrix(source,ref,records,out,interval);audit.setdefault('extra_race_alignments',[]).extend(raceaudit)
                        audit['extra_race_method']='mean of callable individual read frequencies; individual biallelic depth >=8; no marker reselection; n per site may be 1 or 2'
                    fig=plt.figure(figsize=(15,5.5+(len(pknames)+len(racenames))*.25));heat=fig.add_axes([.18,.40,.70,.36]);locus=fig.add_axes([.18,.16,.70,.12]);chrom,ranges=REGIONS[ref];lo,hi=ranges[idx]
                    if records:
                        mat=np.vstack([np.array(ctrl).T,table[['left_rosina_dosage','right_rosina_dosage']].to_numpy().T]+([np.array(pkvals)] if pkvals else [])+([np.array(racevals)] if racevals else []));n=len(records);height=mat.shape[0]
                        image=heat.imshow(mat,cmap=cmap,vmin=0,vmax=1,aspect='auto',interpolation='nearest',extent=(-.5,n-.5,height-.5,-.5));heat.set_yticks(range(height),['rosina controls','melpomene controls','TN Left adjusted','TN Right adjusted']+pknames+racenames);heat.set_xticks([])
                        for j,row in enumerate(records):fig.add_artist(ConnectionPatch(xyA=(j,height-.5),coordsA=heat.transData,xyB=(row['pos1'],.5),coordsB=locus.transData,color='#89939D',lw=.35,alpha=.2))
                        bar=fig.colorbar(image,cax=fig.add_axes([.91,.40,.012,.36]));bar.set_ticks([0,.5,1]);bar.set_label('Rosina-associated frequency / genotype dosage')
                    else:heat.text(.5,.5,'No compatible selected SNPs',ha='center',transform=heat.transAxes);heat.axis('off')
                    locus.axhline(0,color='grey');locus.axvspan(lo,hi,color='grey',alpha=.25,label='Legacy interval')
                    if interval=='d':
                        aa,bb=PRIMERS[ref];locus.plot([aa,bb],[.3,.3],color='#662A83',lw=4);locus.text((aa+bb)/2,.4,'Luca primer span',ha='center',fontsize=8)
                    locus.set_xlim(lo-1000,hi+1000);locus.set_ylim(-.1,.65);locus.set_yticks([]);locus.ticklabel_format(axis='x',style='plain',useOffset=False);locus.set_xlabel(chrom+' position (bp)')
                    called=int(sum(np.isfinite(row['left_rosina_dosage']) for row in records));fig.suptitle(f'{ref} Yb({interval}): {called}/{len(records)} posterior95 genotype pairs')
                    fig.text(.18,.035,('Extra race strips: mean measured read frequency; n per SNP may be 1 or 2 callable controls.\n' if a.add_race_strips else '')+('PK rows: unmixed diploid posterior95 calls; phenotype labels from metadata; only completed alignments included.\n' if a.add_ready_pks else '')+'Controls: measured read frequencies. Focal: merged T2/T3 diploid dosage; pair posterior >=0.95; grey = unavailable.\nOriginal archived beta-binomial likelihood; saved optix fit transferred without refitting. This is not a Yb-calibrated fit.\nRace-associated SNPs; phase and causal Yb genotype not established. Existing fragment counts differ from original VCF AD.',fontsize=8)
                    pdf.savefig(fig);plt.close(fig);audit['panels'].append({'reference':ref,'interval':interval,'selected_compatible':len(records),'posterior95':called})
    (out/'audit.json').write_text(json.dumps(audit,indent=2)+'\n');print('BUNDLE_PDF',pdfpath,sep='\t');print('AUDIT',out/'audit.json',sep='\t')
if __name__=='__main__':main()
