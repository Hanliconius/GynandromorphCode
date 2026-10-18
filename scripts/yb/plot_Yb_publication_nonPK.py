#!/usr/bin/env python3
"""Compact Yb(d/v) panels from restored counts and unchanged focal calls."""
import argparse,hashlib,json,re
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
INTERVALS={'HmRos':[(1633041,1635646),(1658079,1659072)],'Hmel2.5':[(1403500,1405500),(1429500,1430500)]}
PRIMERS={'HmRos':(1633657,1635225),'Hmel2.5':(1404112,1405075)}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def binary(x):
    try:return int(float(x)) if float(x) in (0,1) else None
    except (ValueError,TypeError):return None

def find_source(root):
    candidates=[]
    for p in (root/'analysis/yb_previous_calls_new_samples').glob('*/audit.json'):
        a=json.loads(p.read_text())
        groups={str(r['group']).lower() for r in a.get('new_samples',[])}
        if a.get('status')=='COMPLETE' and a.get('original_SNP_selection_unchanged') is True and {'cythera','vulcanus'}<=groups:
            candidates.append(p.parent)
    if not candidates:raise RuntimeError('No completed restored Yb run including cythera/vulcanus. Supply --source to a completed run.')
    return sorted(candidates)[-1]

def phenotype(sample,group,lookup,interval):
    g=group.lower()
    if g=='vulcanus':return 0 if interval=='d' else 1
    if g=='cythera':return 0
    col='Yb-dorsal' if interval=='d' else 'Yb-ventral'
    vals={binary(r.get(col)) for r in lookup.get(sample,[])}-{None}
    if len(vals)>1:raise RuntimeError('Conflicting phenotype: '+sample+' '+col)
    return next(iter(vals)) if vals else None

def pooled(values,samples,phen):
    rows=[];membership=[]
    for state in (1,0):
        ids=[s for s in samples if phen.get(s)==state]
        if not ids:raise RuntimeError('Empty non-PK phenotype class: '+str(state))
        mat=np.array([values[s] for s in ids]);n=np.isfinite(mat).sum(axis=0)
        mean=np.divide(np.nansum(mat,axis=0),n,out=np.full(mat.shape[1],np.nan),where=n>0)
        rows.append((mean,n));membership.append(dict(phenotype=state,samples=ids,n=len(ids)))
    return rows,membership

def panel(table,matrix,labels,ref,interval,pdf):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.patches import Rectangle,ConnectionPatch
    fig=plt.figure(figsize=(7.2,2.45));fig.text(.025,.94,'A' if interval=='d' else 'B',weight='bold',fontsize=10)
    fig.text(.23,.94,f'Yb({interval})',weight='bold',fontsize=9);fig.text(.94,.94,ref,ha='right',fontsize=6.5)
    locus=fig.add_axes([.23,.73,.70,.09]);heat=fig.add_axes([.23,.25,.70,.39])
    cmap=LinearSegmentedColormap.from_list('rosina',['#D7EDF5','#ECCB63','#662A83']);cmap.set_bad('#C8C8C8')
    n=len(table);im=heat.pcolormesh(np.arange(n+1)-.5,np.arange(5)-.5,matrix,cmap=cmap,vmin=0,vmax=1,shading='flat')
    heat.set_xlim(-.5,n-.5);heat.set_ylim(3.5,-.5);heat.set_yticks(range(4),labels,fontsize=7);heat.set_xticks([]);heat.tick_params(length=0)
    lo,hi=INTERVALS[ref][0 if interval=='d' else 1]
    locus.set_xlim(min(lo-1000,table.pos1.min()-100),max(hi+1000,table.pos1.max()+100));locus.set_ylim(-.1,1)
    locus.axhline(0,color='#888',lw=.45);locus.add_patch(Rectangle((lo,.05),hi-lo,.25,color='#777',lw=0))
    locus.text((lo+hi)/2,.48,f'Yb({interval}) interval',ha='center',fontsize=6)
    if interval=='d':
        x,y=PRIMERS[ref];locus.plot([x,y],[.8,.8],color='#662A83',lw=1.8)
    for j,pos in enumerate(table.pos1):
        fig.add_artist(ConnectionPatch(xyA=(pos,0),coordsA=locus.transData,xyB=(j,-.5),coordsB=heat.transData,color='#9ca2a8',alpha=.22,lw=.25,clip_on=False))
    locus.set_yticks([]);locus.ticklabel_format(axis='x',style='plain',useOffset=False);locus.tick_params(axis='x',labelsize=6,length=2,pad=1)
    cb=fig.colorbar(im,cax=fig.add_axes([.69,.105,.23,.025]),orientation='horizontal');cb.set_ticks([0,.5,1]);cb.ax.tick_params(labelsize=6,length=2,pad=1);cb.outline.set_visible(False)
    cb.set_label('Rosina allele frequency / focal dosage',fontsize=6,labelpad=1)
    fig.text(.23,.025,f'{table.contig.iloc[0]} (bp) · {n} original markers · grey: unavailable',fontsize=6)
    for ax in (locus,heat):
        for sp in ax.spines.values():sp.set_visible(False)
    pdf.savefig(fig);plt.close(fig)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,default=ROOT)
    p.add_argument('--source',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    source=a.source or find_source(a.root);audit=json.loads((source/'audit.json').read_text())
    if audit.get('status')!='COMPLETE' or audit.get('original_SNP_selection_unchanged') is not True:raise RuntimeError('Source not audited restored calls')
    cfg=json.loads((source/'run.json').read_text());out=a.output;out.mkdir(parents=True,exist_ok=False)
    metadata_path=a.root/'data/derived/region_panel_metadata.tsv'
    meta=pd.read_csv(metadata_path,sep='\t',dtype=str).fillna('');meta.columns=[c.strip() for c in meta.columns]
    lookup={};pk=set()
    for r in meta.to_dict('records'):
        tokens=set()
        for k in ('Accession','SequenceID','EarthCapeID'):
            v=r.get(k,'').strip()
            if v:tokens.add(v);tokens.update(re.findall(r'(?:SRR|ERR|DRR|SAMN|SAMEA|ERS)\d+',v))
        for s in tokens:lookup.setdefault(s,[]).append(r)
        if re.search(r'piano|BWPK',r.get('Taxon Name',''),re.I):pk.update(tokens)
    entries={r['sample']:r for r in cfg['entries']}
    pk.update(s for s,r in entries.items() if str(r['group']).startswith('PK_'))
    import matplotlib
    matplotlib.use('Agg');matplotlib.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':7,'pdf.fonttype':42,'ps.fonttype':42})
    from matplotlib.backends.backend_pdf import PdfPages
    pages=[];unknown=[];hashes={};members={}
    with PdfPages(out/'Yb_d_v_nonPK_compact_Rosina.pdf') as mainpdf, PdfPages(out/'Yb_d_v_nonPK_compact_Hmel25.pdf') as supp:
        for ref in ('HmRos','Hmel2.5'):
            for interval in ('d','v'):
                stem=f'{ref}_Yb{interval}';focal=source/(stem+'_unchanged_focal_calls.tsv.gz')
                table=pd.read_csv(focal,sep='\t').sort_values('pos1');positions=table.pos1.tolist()
                if table.empty or table.pos1.duplicated().any():raise RuntimeError('Empty/duplicate markers')
                parts=[]
                rawpath=source/'baseline'/(stem+'_SNP_counts.tsv.gz')
                raw=pd.read_csv(rawpath,sep='\t');raw=raw.loc[raw.role.eq('control')]
                for sample,t in raw.groupby('sample'):
                    parts.append((sample,str(t.race.iloc[0]),t.set_index('pos1').rosina_frequency.reindex(positions).to_numpy()))
                extra=source/'baseline'/(stem+'_extra_race_frequencies.tsv.gz')
                if extra.is_file():
                    d=pd.read_csv(extra,sep='\t')
                    for sample,t in d.groupby('sample'):parts.append((sample,str(t.race.iloc[0]),t.set_index('pos1').frequency.reindex(positions).to_numpy()))
                added=source/(stem+'_new_sample_counts.tsv.gz')
                if added.is_file() and added.stat().st_size>30:
                    try:d=pd.read_csv(added,sep='\t')
                    except pd.errors.EmptyDataError:d=pd.DataFrame(columns=['sample','group','pos1','display_value'])
                    for sample,t in d.groupby('sample'):parts.append((sample,str(t.group.iloc[0]),t.set_index('pos1').display_value.reindex(positions).to_numpy()))
                vals={};phen={}
                for sample,group,v in parts:
                    if sample in pk or group.startswith('PK_'):continue
                    if sample in vals:
                        np.testing.assert_equal(vals[sample],v,err_msg='Conflicting duplicate sample '+sample);continue
                    vals[sample]=v;phen[sample]=phenotype(sample,group,lookup,interval)
                    if phen[sample] is None:unknown.append(dict(reference=ref,interval=interval,sample=sample,group=group))
                rows,membership=pooled(vals,list(vals),phen);members[stem]=membership
                matrix=np.vstack([rows[0][0],table.left_rosina_dosage,table.right_rosina_dosage,rows[1][0]])
                np.testing.assert_equal(matrix[1:3],table[['left_rosina_dosage','right_rosina_dosage']].to_numpy().T)
                labels=[f'Yb({interval})+ (n = {membership[0]["n"]})','TNLeft','TNRight',f'Yb({interval})− (n = {membership[1]["n"]})']
                panel(table,matrix,labels,ref,interval,mainpdf if ref=='HmRos' else supp)
                with PdfPages(out/(stem+'_nonPK_compact.pdf')) as single:panel(table,matrix,labels,ref,interval,single)
                pd.DataFrame(matrix,index=labels,columns=positions).to_csv(out/(stem+'_plotted_values.tsv'),sep='\t',na_rep='NA')
                pd.DataFrame(dict(pos1=positions,positive_n_called=rows[0][1],negative_n_called=rows[1][1])).to_csv(out/(stem+'_callability.tsv'),sep='\t',index=False)
                hashes[str(focal)]=sha(focal);hashes[str(rawpath)]=sha(rawpath)
                pages.append(dict(reference=ref,interval=interval,markers=len(table),positive_n=membership[0]['n'],negative_n=membership[1]['n'],focal_identical=True))
    pd.DataFrame(unknown,columns=['reference','interval','sample','group']).to_csv(out/'unknown_phenotypes_excluded.tsv',sep='\t',index=False)
    (out/'control_membership.json').write_text(json.dumps(members,indent=2)+'\n')
    result=dict(status='COMPLETE',source=str(source),source_audit_sha256=sha(source/'audit.json'),input_hashes=hashes,metadata_sha256=sha(metadata_path),
        SNP_selection_changed=False,calling_changed=False,parameters_refitted=False,PK_in_controls=False,
        control_statistic='Equal individual mean stored read frequency; original per-sample rules and available-value averaging preserved',
        polarity='Rosina allele, not asserted causal Yb allele',phenotype_overrides={'cythera':{'d':0,'v':0},'vulcanus':{'d':0,'v':1}},
        unknown_phenotypes=unknown,pages=pages)
    (out/'audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print('SOURCE\t'+str(source));print('PDF\t'+str(out/'Yb_d_v_nonPK_compact_Rosina.pdf'))
    for page in pages:print('PANEL\t'+json.dumps(page))
    print('UNKNOWN_EXCLUDED\t'+str(len(unknown)))
if __name__=='__main__':main()
