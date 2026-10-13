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
def main():
    p=argparse.ArgumentParser();p.add_argument('--counts-run',required=True,type=Path);p.add_argument('--fit-summary',required=True,type=Path);p.add_argument('--focal-run',required=True,type=Path);a=p.parse_args()
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
    audit={'caller':str(caller),'caller_sha256':digest(caller),'fit_summary':str(a.fit_summary),'fit_summary_sha256':digest(a.fit_summary),'parameters':dict(zip(keys,pars)),'parameters_refitted':False,'fit_transfer':'saved optix nuisance parameters applied to new Yb counts on both assemblies; not a Yb/reference-specific fit','prior':'equal over nine diploid genotype pairs','posterior_cutoff':.95,'min_count_depth_per_merged_side':2,'aggregation':'T2 and T3 counts summed within side','count_method':'existing MAPQ30/baseQ30 fragment counts; overlapping mates counted once per site','count_method_differs_from_original_VCF_AD':True,'selection_uses_focal':False,'other_alleles':'retained in table; excluded from biallelic likelihood denominator','intervals':'legacy annotations; Luca primer span an anchor, not exact functional boundaries','inputs':[],'panels':[]}
    pdfpath=out/'Yb_d_v_posterior95_corrected_genotypes_both_references.pdf'
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
                    fig=plt.figure(figsize=(12,5.5));heat=fig.add_axes([.18,.40,.70,.36]);locus=fig.add_axes([.18,.16,.70,.12]);chrom,ranges=REGIONS[ref];lo,hi=ranges[idx]
                    if records:
                        mat=np.vstack([np.array(ctrl).T,table[['left_rosina_dosage','right_rosina_dosage']].to_numpy().T]);n=len(records)
                        image=heat.imshow(mat,cmap=cmap,vmin=0,vmax=1,aspect='auto',interpolation='nearest',extent=(-.5,n-.5,3.5,-.5));heat.set_yticks(range(4),['rosina controls','melpomene controls','TN Left adjusted','TN Right adjusted']);heat.set_xticks([])
                        for j,row in enumerate(records):fig.add_artist(ConnectionPatch(xyA=(j,3.5),coordsA=heat.transData,xyB=(row['pos1'],.5),coordsB=locus.transData,color='#89939D',lw=.35,alpha=.2))
                        bar=fig.colorbar(image,cax=fig.add_axes([.91,.40,.012,.36]));bar.set_ticks([0,.5,1]);bar.set_label('Rosina-associated frequency / genotype dosage')
                    else:heat.text(.5,.5,'No compatible selected SNPs',ha='center',transform=heat.transAxes);heat.axis('off')
                    locus.axhline(0,color='grey');locus.axvspan(lo,hi,color='grey',alpha=.25,label='Legacy interval')
                    if interval=='d':
                        aa,bb=PRIMERS[ref];locus.plot([aa,bb],[.3,.3],color='#662A83',lw=4);locus.text((aa+bb)/2,.4,'Luca primer span',ha='center',fontsize=8)
                    locus.set_xlim(lo-1000,hi+1000);locus.set_ylim(-.1,.65);locus.set_yticks([]);locus.ticklabel_format(axis='x',style='plain',useOffset=False);locus.set_xlabel(chrom+' position (bp)')
                    called=sum(np.isfinite(row['left_rosina_dosage']) for row in records);fig.suptitle(f'{ref} Yb({interval}): {called}/{len(records)} posterior95 genotype pairs')
                    fig.text(.18,.035,'Controls: measured read frequencies. Focal: merged T2/T3 diploid dosage; pair posterior >=0.95; grey = unavailable.\nOriginal archived beta-binomial likelihood; saved optix fit transferred without refitting. This is not a Yb-calibrated fit.\nRace-associated SNPs; phase and causal Yb genotype not established. Existing fragment counts differ from original VCF AD.',fontsize=8)
                    pdf.savefig(fig);plt.close(fig);audit['panels'].append({'reference':ref,'interval':interval,'selected_compatible':len(records),'posterior95':called})
    (out/'audit.json').write_text(json.dumps(audit,indent=2)+'\n');print('BUNDLE_PDF',pdfpath,sep='\t');print('AUDIT',out/'audit.json',sep='\t')
if __name__=='__main__':main()
