#!/usr/bin/env python3
"""Geography-explicit Yb(v) supplements from existing, reviewed qualified calls.

Run on Pegasus. This diagnostic does not change the main Figure 3 cohort,
markers, calls, mixture estimates, or any existing Optix supplement.
"""
import argparse, getpass, hashlib, importlib.util, json, shlex, shutil, subprocess, sys
from pathlib import Path

ROOT = Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
SOURCE = ROOT/'analysis/Figure3_for_Arnaud'

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def save(p, value):
    Path(p).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')

def load(p, name):
    spec = importlib.util.spec_from_file_location(name, p)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

def prepare(a):
    out = a.output
    if (out/'complete.json').exists() and not a.overwrite:
        raise RuntimeError('Use --overwrite to replace this diagnostic in place')
    if (out/'jobs.json').is_file():
        job = json.loads((out/'jobs.json').read_text())['plot']
        active = subprocess.check_output(['squeue','-h','-u',getpass.getuser(),'-o','%A'], text=True).split()
        if job in active:
            raise RuntimeError('Existing diagnostic job is still active: '+job)
    repo = Path(__file__).resolve().parents[3]
    inputs = [a.source/'complete.json', a.source/'run.json',
              a.source/'reviewed_specimens.tsv', a.source/'code/worker.py',
              a.source/'tables/Ybv_focus_markers_and_updated_BA.tsv']
    inputs += [a.source/'data/Ybv'/n for n in
               ['complete.json','membership.tsv','all_SNP_BA.tsv.gz','all_ALT_dosages.tsv.gz']]
    inputs += [repo/'scripts/figures/figure3'/n for n in
               ['build_Figure3_morph_supplements.py','compare_Ybv_geography_BA.py']]
    for p in inputs:
        if not p.is_file():
            raise FileNotFoundError(p)
    cfg = json.loads((a.source/'run.json').read_text())
    if cfg['metadata_sha256'] != sha(a.source/'reviewed_specimens.tsv'):
        raise RuntimeError('Reviewed snapshot differs from the completed Figure 3 inputs')
    if json.loads((a.source/'complete.json').read_text())['status'] != 'COMPLETE':
        raise RuntimeError('Source Figure 3 is incomplete')
    out.mkdir(parents=True, exist_ok=True)
    (out/'logs').mkdir(exist_ok=True)
    for src, name in [(Path(__file__),'worker.py'),
                      (a.source/'code/worker.py','figure3_helper.py'),
                      (inputs[-2],'morph_helper.py'),(inputs[-1],'geography_helper.py')]:
        shutil.copy2(src, out/name)
    helper = load(out/'figure3_helper.py','source_figure3')
    inputs += [helper.BROAD/n for n in ['membership.tsv','SNPs.tsv','individual_ALT_dosages.tsv.gz']]
    inputs += [helper.YF/'individual_ALT_dosages.tsv.gz']
    for p in inputs:
        if not p.is_file():
            raise FileNotFoundError(p)
    save(out/'run.json', dict(source=str(a.source), cutoff=a.cutoff,
         hashes={str(p):sha(p) for p in inputs}, workflow='Ybv_geography_supplements',
         metadata_authority='Completed Figure3_for_Arnaud/reviewed_specimens.tsv',
         cohort='All available reviewed comparative specimens with known ventral phenotype; global exclusions and PK exclusion retained',
         region='Recorded reviewed West/East assignments; unresolved displayed separately, never assigned by genotype',
         BA='Homozygote rule, >=80% callability separately in each phenotype class and region',
         selection='Original 27 markers plus union of sites passing the declared threshold in West or East independently',
         orientation='West-best direction held fixed for both region tracks; East-best fallback only where West is unscoreable',
         main_figures_changed=False, downloads=0, alignments=0, calling=0))
    (out/'complete.json').unlink(missing_ok=True)
    command = [sys.executable,str(out/'worker.py'),'--run','--output',str(out)]
    job = subprocess.check_output(['sbatch','--parsable','--partition=nano',
        '--time=00:30:00','--cpus-per-task=1','--mem=8G','--job-name=Ybv_geoSupp',
        '--output='+str(out/'logs/plot.%j.out'),'--error='+str(out/'logs/plot.%j.err'),
        '--wrap='+shlex.join(command)], text=True).strip().split(';')[0]
    save(out/'jobs.json',dict(plot=job))
    print('JOB\t'+job+'\nOUTPUT\t'+str(out))

def regional_scores(sites, values, meta, geo):
    import numpy as np
    result = sites[['chrom','pos','ref','alt']].copy()
    regions = meta.geography_used_Fig3b.to_numpy()
    y = meta.phenotype.to_numpy(int)
    raw = {r:geo.raw_scores(values,y,regions==r,'homozygote') for r in ['West','East']}
    for region, rr in raw.items():
        for name in ['ALT','REF','n_total_0','n_total_1','n_called_0','n_called_1','rate_0','rate_1']:
            result[region+'_'+name] = rr[name]
        result[region+'_best_BA'] = np.maximum(rr['ALT'],rr['REF'])
        result[region+'_best_orientation'] = np.where(~np.isfinite(rr['ALT']),'unscored',
            np.where(np.isclose(rr['ALT'],rr['REF'],atol=1e-12,rtol=0),'tie',
                     np.where(rr['ALT']>rr['REF'],'ALT','REF')))
    west = raw['West']; east = raw['East']
    west_good = np.isfinite(west['ALT']) & ~np.isclose(west['ALT'],west['REF'],atol=1e-12,rtol=0)
    direction = np.where(west_good,west['ALT']>=west['REF'],east['ALT']>=east['REF'])
    either = west_good | np.isfinite(east['ALT'])
    result['track_orientation'] = np.where(either,np.where(direction,'ALT','REF'),'unscored')
    result['orientation_source'] = np.where(west_good,'West',np.where(either,'East fallback','unscored'))
    for region, rr in raw.items():
        result[region+'_fixed_direction_BA'] = np.where(either,np.where(direction,rr['ALT'],rr['REF']),np.nan)
    return result

def grouped_rows(meta, values, sites, helper):
    import numpy as np, pandas as pd
    rows, labels, exports, blocks = [], [], [], []
    regions = ['West','East'] + sorted(set(meta.geography_used_Fig3b)-{'West','East'})
    for region in regions:
        for state in [1,0]:
            group = meta.loc[meta.geography_used_Fig3b.eq(region) & meta.phenotype.eq(state)]
            if group.empty:
                continue
            first = len(rows)
            def append(g, name, kind, population='', morph=''):
                x = values[g.index]; n = np.isfinite(x).sum(axis=0)
                freq = np.divide(np.nansum(x,axis=0),n,out=np.full(len(sites),np.nan),where=n>0)
                rows.append(freq); labels.append(name+f' | n={len(g)}')
                for j,r in enumerate(sites.itertuples()):
                    exports.append(dict(row=len(rows)-1,row_type=kind,geography=region,
                        phenotype=state,population=population,morph=morph,chrom=r.chrom,
                        pos=r.pos,ref=r.ref,alt=r.alt,n_total=len(g),n_called=int(n[j]),
                        ALT_frequency=freq[j],
                        REF_REF_fraction=np.sum(x[:,j]==0)/n[j] if n[j] else np.nan,
                        REF_ALT_fraction=np.sum(x[:,j]==.5)/n[j] if n[j] else np.nan,
                        ALT_ALT_fraction=np.sum(x[:,j]==1)/n[j] if n[j] else np.nan))
            append(group,f'All {region} v{"+" if state else "−"}','region_phenotype')
            for (pop,morph),g in group.groupby(['population','morph'],sort=True):
                append(g,'    '+helper.short(pop)+' | '+morph,'population_morph',pop,morph)
            blocks.append((first,len(rows),region,state))
    return np.asarray(rows),labels,pd.DataFrame(exports),blocks

def draw(path, selected, all_scores, values, meta, helper, title, cutoff):
    import numpy as np, matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.patches import ConnectionPatch, Rectangle
    rows,labels,exports,blocks = grouped_rows(meta,values,selected,helper)
    cmap = LinearSegmentedColormap.from_list('previous_frequency_palette',helper.SUPP,N=257)
    cmap.set_bad('#C8C8C8')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'pdf.fonttype':42})
    lo,hi = helper.SPANS['Ybv'][1:]
    # Exhaustive pagination only when necessary; never thin or rank-cap SNPs.
    pages = (len(selected)+159)//160
    with PdfPages(path) as pdf:
        for first in range(0,len(selected),160):
            ss = selected.iloc[first:first+160]; width = len(ss); nr = len(labels)
            fig = plt.figure(figsize=(17,max(10,5.4+.23*nr)))
            left,ww = .34,.64
            heat = fig.add_axes([left,.16,ww,.53])
            heat.pcolormesh(np.arange(width+1)-.5,np.arange(nr+1)-.5,
                np.ma.masked_invalid(rows[:,first:first+width]),cmap=cmap,vmin=0,vmax=1)
            heat.set(xlim=(-.5,width-.5),ylim=(nr-.5,-.5),xticks=[])
            heat.set_yticks(np.arange(nr),labels,fontsize=7); heat.tick_params(length=0,pad=6)
            for start,end,region,state in blocks:
                heat.axhline(start-.5,c='white',lw=2)
                heat.get_yticklabels()[start].set_weight('bold')
            ann = fig.add_axes([left,.735,ww,.04])
            ann.set(xlim=(lo,hi),ylim=(0,1.6),yticks=[])
            ann.axhline(.1,c='#888888',lw=.5)
            for j,(name,a,b) in enumerate([('cortex',1642343,1653300),('Yb(v)',1658079,1659072),('ivory (partial)',1659008,1741475)]):
                a,b=max(a,lo),min(b,hi); y=.4+(j%2)*.45
                ann.add_patch(Rectangle((a,y),b-a,.13,color='#777777'))
                ann.text((a+b)/2,y+.18,name,ha='center',fontsize=7)
            ann.set_xticks(np.linspace(lo,hi,5),[f'{x/1e6:.3f}' for x in np.linspace(lo,hi,5)])
            ann.set_xlabel('Rosina chromosome 15 (Mb)')
            for j,r in enumerate(ss.itertuples()):
                fig.add_artist(ConnectionPatch((r.pos,.1),(j,-.5),coordsA=ann.transData,
                    coordsB=heat.transData,color='#BBBBBB',lw=.25,alpha=.5))
                if r.original_marker:
                    heat.plot(j,-.5,marker='s',ms=2,c='#111111',clip_on=False)
            for region,bottom in [('West',.875),('East',.805)]:
                ax = fig.add_axes([left,bottom,ww,.045])
                y = all_scores[region+'_fixed_direction_BA']
                ax.scatter(all_scores.pos,y,s=3,c='#555555',alpha=.55,lw=0)
                ax.scatter(ss.pos,ss[region+'_fixed_direction_BA'],s=6,c='#1764AB',lw=0)
                ax.axhline(.5,c='#999999',ls=':',lw=.5);ax.axhline(cutoff,c='#999999',ls='--',lw=.5)
                ax.set(xlim=(lo,hi),ylim=(0,1.02),yticks=[0,.5,1])
                ax.set_ylabel(region+' BA');ax.tick_params(axis='x',bottom=False,labelbottom=False)
                for sp in ['top','right']:ax.spines[sp].set_visible(False)
            for ax in [heat,ann]:
                for sp in ax.spines.values():sp.set_visible(False)
            fig.text(.025,.98,title,fontsize=13,weight='bold',va='top')
            fig.text(.025,.95,f'{len(selected)} retained SNPs | columns {first+1}–{first+width} | page {first//160+1}/{pages}',fontsize=8)
            cb = fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(0,1),cmap=cmap),
                cax=fig.add_axes([left,.11,ww,.012]),orientation='horizontal',ticks=[0,.5,1])
            cb.set_label('ALT allele frequency (mean individual genotype dosage / 2)')
            fig.text(.025,.025,
                'West v+, West v−, East v+, East v−: regional summary followed by all available reviewed population/morph rows. Unresolved geography, if present, is separate.\n'
                'BA: homozygote rule; ≥80% called in each phenotype class within each region. Same allele direction across regions: West-best, with labelled East fallback in exported tables.\n'
                'Wide selection: original 27 plus every site passing the declared BA threshold in either region independently; BA is descriptive, not independent validation.\n'
                'Black squares mark original sites. Grey cells = missing; grey BA background = all regional scores; blue BA points = displayed sites. Missing BA is not zero.\n'
                'n = specimens; per-site called denominators and genotype fractions exported. Frequency 0.5 does not imply all specimens are heterozygous. No focal genotypes extrapolated.',fontsize=7)
            pdf.savefig(fig);plt.close(fig)
    return exports,pages

def run(out):
    import numpy as np, pandas as pd
    cfg = json.loads((out/'run.json').read_text()); source = Path(cfg['source'])
    for p,h in cfg['hashes'].items():
        if sha(p)!=h:
            raise RuntimeError('Input changed after submission: '+p)
    helper = load(out/'figure3_helper.py','frozen_figure3')
    morph = load(out/'morph_helper.py','morph_helper')
    geo = load(out/'geography_helper.py','geo_helper')
    folder = source/'data/Ybv'
    if json.loads((folder/'complete.json').read_text())['status'] != 'COMPLETE':
        raise RuntimeError('Yb source incomplete')
    scores = pd.read_csv(folder/'all_SNP_BA.tsv.gz',sep='\t')
    selection = pd.read_csv(folder/'membership.tsv',sep='\t',dtype=str,keep_default_na=False)
    frozen = morph.matrix(folder/'all_ALT_dosages.tsv.gz',selection,scores)
    meta,values,conflicts = morph.yb_display(helper,source,scores,selection,frozen)
    meta = meta.reset_index(drop=True)
    if meta.specimen_id.duplicated().any() or any(helper.blocked(r) for r in meta.to_dict('records')):
        raise RuntimeError('Duplicate or globally excluded specimen')
    if np.any(np.isfinite(values)&~np.isin(values,[0,.5,1])):
        raise RuntimeError('Invalid diploid dosage')
    meta['population_original']=meta.population
    meta['population']=meta.population.map(morph.population_name)
    meta['morph']=['d'+str(r.Yb_dorsal)+' v'+str(r.Yb_ventral) for r in meta.itertuples()]
    meta['used_in_regional_BA']=meta.geography_used_Fig3b.isin(['West','East'])
    meta.to_csv(out/'membership.tsv',sep='\t',index=False)
    pd.DataFrame(conflicts).to_csv(out/'source_conflicts.tsv',sep='\t',index=False)
    counts=meta.groupby(['geography_used_Fig3b','phenotype']).size().reset_index(name='specimens')
    counts.to_csv(out/'cohort_counts.tsv',sep='\t',index=False)
    regional = regional_scores(scores,values,meta,geo)
    regional['_index']=np.arange(len(regional))
    keys=['chrom','pos','ref','alt']
    if regional.duplicated(keys).any():raise RuntimeError('Duplicate allele identity')
    focus = pd.read_csv(source/'tables/Ybv_focus_markers_and_updated_BA.tsv',sep='\t')
    focus_keys={helper.key(r.chrom,r.pos,r.ref,r.alt) for r in focus.itertuples()}
    if len(focus_keys)!=27:raise RuntimeError('Expected 27 exact original marker identities')
    regional['original_marker']=[helper.key(r.chrom,r.pos,r.ref,r.alt) in focus_keys for r in regional.itertuples()]
    if regional.original_marker.sum()!=27:raise RuntimeError('Original markers absent from completed full call matrix')
    regional=regional.sort_values(keys).reset_index(drop=True)
    regional['passes_West']=regional.West_best_BA.ge(cfg['cutoff'])
    regional['passes_East']=regional.East_best_BA.ge(cfg['cutoff'])
    regional.to_csv(out/'all_regional_BA_and_denominators.tsv.gz',sep='\t',index=False)
    threshold_rows=[]
    for t in [.80,.85,.90,.93,.95]:
        w,e=regional.West_best_BA.ge(t),regional.East_best_BA.ge(t)
        threshold_rows.append(dict(threshold=t,West=int(w.sum()),East=int(e.sum()),union=int((w|e).sum()),
            union_plus_original=int((w|e|regional.original_marker).sum())))
    pd.DataFrame(threshold_rows).to_csv(out/'threshold_counts.tsv',sep='\t',index=False)
    outputs=[]
    for filename,mask,title in [
        ('Ybv_geography_27_original.pdf',regional.original_marker,'Yb(v): original 27 SNPs by geography and morph'),
        ('Ybv_geography_regional_union.pdf',regional.original_marker|regional.passes_West|regional.passes_East,
         f'Yb(v): original markers + regional BA ≥{cfg["cutoff"]:.2f} union')]:
        sites=regional.loc[mask].copy()
        v=values[:,sites['_index'].to_numpy(int)]
        exports,pages=draw(out/filename,sites,regional,v,meta,helper,title,cfg['cutoff'])
        stem=Path(filename).stem
        sites.to_csv(out/(stem+'_sites.tsv'),sep='\t',index=False)
        exports.to_csv(out/(stem+'_frequencies_and_callability.tsv.gz'),sep='\t',index=False)
        outputs.append(dict(PDF=str(out/filename),SNPs=len(sites),pages=pages))
    save(out/'complete.json',dict(status='COMPLETE',outputs=outputs,individuals=len(meta),
        cohorts=counts.to_dict('records'),main_figures_changed=False,main_BA_cohort_changed=False,
        diagnostic_cohort='All available reviewed comparative specimens; original restricted main BA flags are not used to discard geographic exceptions',
        source_config_sha256=sha(out/'run.json'),downloads=0,alignments=0,calling=0))
    for row in outputs:print('PDF\t'+row['PDF']+'\tSNPs='+str(row['SNPs']))
    print('STATUS\tCOMPLETE')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=SOURCE)
    p.add_argument('--output',type=Path,default=SOURCE/'Ybv_geography_supplements')
    p.add_argument('--cutoff',type=float,default=.90)
    p.add_argument('--submit',action='store_true');p.add_argument('--run',action='store_true')
    p.add_argument('--overwrite',action='store_true');a=p.parse_args()
    if not .5<=a.cutoff<=1:p.error('--cutoff must be between 0.5 and 1')
    if a.submit:prepare(a)
    elif a.run:run(a.output)
    else:p.error('Use --submit')

if __name__=='__main__':main()
