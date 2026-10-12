#!/usr/bin/env python3
"""Three focused alternatives for centroid clashes; no simulations/refitting."""
import json
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, Normalize
from scipy.stats import betabinom
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'figure2'))
from common import set_publication_style
ROOT=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
STATES=np.array([(l,r) for l in (0.,.5,1.) for r in (0.,.5,1.)])
GRID=np.arange(41)/40

def distributions(cl,cr):
    result=[]
    centres=[]
    for l,r in STATES:
        ql=(1-cl)*l+cl*r; qr=cr*l+(1-cr)*r
        def pmf(q):
            if q<=0 or q>=1:
                p=np.zeros(41);p[0 if q<=0 else 40]=1;return p
            p=betabinom.pmf(np.arange(41),40,q*22,(1-q)*22)
            return p/p.sum()
        p=np.outer(pmf(qr),pmf(ql))
        ranked=np.sort(p.ravel())[::-1]
        cut=ranked[np.searchsorted(np.cumsum(ranked),.5)]
        result.append((p,p>=cut));centres.append((ql,qr))
    return result,np.array(centres)

def overlap(a,b):
    p,ma=a;q,mb=b
    return float(np.minimum(p,q)[ma&mb].sum())

def strongest(ds):
    return max(overlap(ds[i],ds[j]) for i in range(9) for j in range(i))

def main():
    out=ROOT/'analysis/arnaud_mixing_supplement/clash_options'
    out.mkdir(parents=True,exist_ok=True)
    cmap=ListedColormap(json.loads((ROOT/'palettes/barfbow/luka_bleute.json').read_text()))
    set_publication_style(9)
    fig,axes=plt.subplots(2,3,figsize=(9.4,7.0),sharex=True,sharey=True)
    fig.subplots_adjust(left=.08,right=.97,bottom=.16,top=.90,hspace=.24,wspace=.13)
    for ax,c in zip(axes.flat,(0.,.1,.2,.3,.4,.5)):
        ds,centres=distributions(c,c)
        # Count distinct underlying states, even when their predicted centres coincide.
        multiplicity=np.sum([mask for _,mask in ds],axis=0)
        if multiplicity.max():
            ax.contourf(GRID,GRID,(multiplicity>0).astype(float),levels=[.5,1.5],colors=['#dfe5e9'],alpha=.6)
        if multiplicity.max()>=2:
            ax.contourf(GRID,GRID,(multiplicity>=2).astype(float),levels=[.5,1.5],colors=['#B61F72'],alpha=.45)
        ax.scatter(centres[:,0],centres[:,1],s=13,c='#333333',zorder=5,clip_on=False)
        ax.set(xlim=(-.025,1.025),ylim=(-.025,1.025),aspect='equal',xticks=[0,.5,1],yticks=[0,.5,1])
        ax.set_title(f'{c:.0%} mixing per side')
        ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
    fig.supxlabel('Left ALT-read fraction',y=.10);fig.supylabel('Right ALT-read fraction',x=.01)
    fig.text(.08,.05,'Grey: conditional 50% predictive regions. Magenta: intersection of two or more genotype-state regions.',fontsize=8)
    fig.text(.08,.022,'40 reads/side; concentration 22; no sequencing error or bias. Overlap describes ambiguity, not a significance test.',fontsize=8)
    fig.savefig(out/'A_clash_regions_by_mixing.pdf');plt.close(fig)
    parameter_grid=np.linspace(0,.5,21)
    heat=np.zeros((21,21))
    for ri,cr in enumerate(parameter_grid):
        for li,cl in enumerate(parameter_grid):
            ds,_=distributions(cl,cr);heat[ri,li]=strongest(ds)
    fig,ax=plt.subplots(figsize=(6.4,5.8),layout='constrained')
    image=ax.pcolormesh(parameter_grid*100,parameter_grid*100,heat,cmap=cmap,norm=Normalize(0,1),shading='nearest',rasterized=False)
    cb=fig.colorbar(image,ax=ax);cb.set_label('Largest shared predictive probability mass')
    ax.set(xlabel='Reads in Left contributed by Right (%)',ylabel='Reads in Right contributed by Left (%)',title='Where genotype predictions overlap',aspect='equal')
    fit=json.loads((ROOT/'analysis/genomewide_mixing_corrected_v1/expanded_calibration/observed_expanded.json').read_text())['full']['alternative']
    cl,cr=fit['c_left'],fit['c_right']
    ax.plot(cl*100,cr*100,'o',mfc='white',mec='#E600A9',mew=2,ms=8)
    ax.annotate('Saved fit',(cl*100,cr*100),xytext=(8,-15),textcoords='offset points',fontsize=9,color='#E600A9')
    fig.savefig(out/'B_clash_heatmap_left_right.pdf');plt.close(fig)
    mixing=np.linspace(0,.5,101)
    pairs=[(1,3,'Opposite sides: 0 / 0.5 vs 0.5 / 0'),(2,6,'Opposite sides: 0 / 1 vs 1 / 0'),(1,4,'Different vs equal: 0 / 0.5 vs 0.5 / 0.5')]
    fig,ax=plt.subplots(figsize=(7.8,5.0),layout='constrained')
    curves=[]
    for index,(i,j,label) in enumerate(pairs):
        values=[]
        for c in mixing:
            ds,_=distributions(c,c);values.append(overlap(ds[i],ds[j]))
        curves.append(values)
        ax.plot(mixing*100,values,label=label,color=['#5A4FA3','#E56B00','#B61F72'][index],lw=1.8)
    ax.set(xlabel='Symmetric mixing per side (%)',ylabel='Shared probability mass within both 50% regions',xlim=(0,50),ylim=(0,1),title='How genotype-region clashes develop')
    ax.legend(frameon=False,fontsize=8,loc='upper left')
    ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
    fig.savefig(out/'C_clash_progression_curves.pdf');plt.close(fig)
    (out/'README.txt').write_text('All alternatives compare genotype states at the SAME mixing parameters. Conditional discrete 50% highest-density regions include threshold ties. Assumptions: beta-binomial depth 40 per side, concentration 22, independent side sampling, no error or bias. Grey and magenta regions in A indicate coverage and overlap; they do not test significance. In B and C shared mass is sum(min(P_i,P_j)) over the intersection of their 50% regions. B takes the largest value over all distinct genotype-state pairs. The colour scale is linear 0-1, not site density. Shared mass may exceed 0.5 because discrete threshold ties are included. C shows selected pairs rather than a universal detection threshold; jaggedness reflects discrete highest-density boundaries. Fixed and coincident states remain distinct biological states in overlap calculations. No observations are simulated and no parameters are refitted. The saved-fit marker in B identifies the verified million-site autosomal alternative estimate.\n')
    np.savez_compressed(out/'clash_metrics.npz',mixing=mixing,curves=np.array(curves),parameter_grid=parameter_grid,heat=heat)
    print(f'WROTE three alternatives: {out}')
if __name__=='__main__':main()
