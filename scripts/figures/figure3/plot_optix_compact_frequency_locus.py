#!/usr/bin/env python3
"""Compact existing optix frequency values and connect them to locus coordinates."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import ConnectionPatch,Rectangle
plt.rcParams.update({"pdf.fonttype":42,"font.size":8})
KEY=["chrom","pos","ref","alt"]
MODULES=("Band","Dennis","Ray")
# Existing OPTIX_POSITION_INTERVALS analysis windows (90b workflow).
# These are analysis windows, not SNP extrema or exact regulatory boundaries.
WINDOWS=(("B",1111500,1117000,"#087DBB"),
         ("R1",1156000,1161000,"#CC70A7"),
         ("R2",1174000,1189000,"#CC70A7"),
         ("D1",1215000,1217000,"#D79500"),
         ("D2",1218000,1223000,"#D79500"))

def layout(values):
    blocks=[];matrices=[];offset=0;mapping=[]
    for module in MODULES:
        table=values.loc[values.module.eq(module)].copy()
        sites=table[KEY].drop_duplicates().sort_values(KEY).reset_index(drop=True)
        rows=["absent","TN_Left","TN_Right","present"] if module=="Band" else ["present","TN_Left","TN_Right","absent"]
        if sites.empty:continue
        sites["column"]=np.arange(offset,offset+len(sites))
        table=table.merge(sites,on=KEY,validate="many_to_one")
        if table.duplicated(["row","column"]).any():raise RuntimeError("Duplicate genotype cell")
        matrix=table.pivot(index="row",columns="column",values="frequency").reindex(index=rows,columns=sites.column)
        matrices.append(matrix.to_numpy(dtype=float))
        blocks.append((module,offset,len(sites),rows))
        mapping.append(sites.assign(module=module))
        offset+=len(sites)
    if not matrices:raise RuntimeError("No frequency values")
    return np.concatenate(matrices,axis=1),blocks,pd.concat(mapping,ignore_index=True)

def render(values,input_dir,zygosity=False):
    matrix,blocks,mapping=layout(values)
    if zygosity:
        # Orient each trait-associated marker toward the DR phenotype endpoint.
        # B absence is the DR endpoint; D/R presence is the DR endpoint.
        for module,start,count,rows in blocks:
            if module=="Band":
                matrix[1:3,start:start+count]=1-matrix[1:3,start:start+count]
        matrix[0,:]=1
        matrix[3,:]=0
    row_labels=(["DR+/DR+ reference","TN_Left adjusted","TN_Right adjusted","B+/B+ reference"]
                if zygosity else
                ["B absent / D+ / R+","TN_Left adjusted","TN_Right adjusted","B present / D− / R−"])
    cmap=LinearSegmentedColormap.from_list("pattern_frequency",[(0,"#D7EDF5"),(.5,"#ECCB63"),(1,"#662A83")])
    cmap.set_bad("#C8C8C8")
    fig=plt.figure(figsize=(13,5.2))
    locus=fig.add_axes([.17,.70,.75,.15])
    heat=fig.add_axes([.17,.20,.75,.32])
    locus.set_xlim(1.06,1.24);locus.set_ylim(-.4,1.4)
    locus.hlines(0,1.06,1.24,color="#616B74",lw=1.1)
    for name,start,end,color in WINDOWS:
        locus.add_patch(Rectangle((start/1e6,.15),(end-start)/1e6,.25,facecolor=color,edgecolor="none"))
        y=.65 if name!="D2" else 1.03
        locus.text((start+end)/2e6,y,name,ha="center",color=color,fontweight="bold")
    locus.set_yticks([])
    locus.set_xticks(np.arange(1.06,1.241,.02))
    locus.set_xlabel("Hmr18 genomic position (Mb)",labelpad=3)
    locus.set_title("optix locus: annotated analysis windows",loc="left",fontsize=11,pad=16)
    for spine in locus.spines.values():spine.set_visible(False)
    n=matrix.shape[1]
    for row in mapping.itertuples():
        connection=ConnectionPatch(xyA=(row.pos/1e6,-.35),coordsA=locus.transData,
                    xyB=(row.column,-.5),coordsB=heat.transData,
                    color="#89939D",lw=.35,alpha=.20,clip_on=False,zorder=0)
        fig.add_artist(connection)
    image=heat.imshow(matrix,cmap=cmap,vmin=0,vmax=1,aspect="auto",interpolation="nearest",
                      extent=(-.5,n-.5,3.5,-.5))
    heat.set_yticks(range(4),row_labels)
    heat.set_xticks([]);heat.tick_params(axis="y",length=0)
    for module,start,count,rows in blocks:
        if start:heat.axvline(start-.5,color="white",lw=2)
        label={"Band":"B","Dennis":"D","Ray":"R"}[module]
        heat.text(start+(count-1)/2,1.04,label,transform=heat.get_xaxis_transform(),
                  ha="center",va="bottom",fontweight="bold",fontsize=10)
    for spine in heat.spines.values():spine.set_visible(False)
    bar=fig.colorbar(image,cax=fig.add_axes([.945,.20,.012,.32]))
    bar.set_ticks([0,.5,1])
    bar.set_label("DR-oriented adjusted allele dosage" if zygosity else "Allele frequency / adjusted dosage",fontsize=8)
    if zygosity:bar.set_ticklabels(["0","0.5","1"])
    note=("Top/bottom are schematic homozygous phenotype endpoints, not measured control genotypes. B markers are inverted; D/R retain presence orientation.\n"
          "Left/Right are unphased adjusted dosages: 0 = B-like homozygote, 0.5 = heterozygote, 1 = DR-like homozygote; grey = unavailable.\n"
          "Endpoint labels assume complementary B versus DR alleles at each marker; linkage/phase is not inferred."
          if zygosity else
          "Top row: Band absent, Dennis present, Ray present. Bottom row: Band present, Dennis absent, Ray absent.\n"
          "Comparative rows show measured population frequencies; Left/Right show adjusted dosage. Grey = unavailable adjusted call.")
    fig.text(.17,.065,note+"\nConnections map compact SNP columns to genomic positions; annotated windows are analysis intervals.",fontsize=7)
    out=input_dir/("optix_compact_DR_B_zygosity_linked_locus.pdf" if zygosity else "optix_compact_frequency_linked_locus.pdf")
    fig.savefig(out);fig.savefig(out.with_suffix(".png"),dpi=250);plt.close(fig)
    mapping.to_csv(input_dir/"compact_column_to_genomic_position.tsv",sep="\t",index=False)
    (input_dir/("compact_DR_B_zygosity.audit.json" if zygosity else "compact_layout.audit.json")).write_text(json.dumps({
        "input":str(input_dir/"frequency_plot_values.tsv.gz"),
        "selection_changed":False,"schematic_references":zygosity,
        "orientation":"DR endpoint; Band inverted" if zygosity else "trait-presence allele",
        "phase_inferred":False,"layout":"consecutive Band, Dennis, Ray columns; separate memberships retained",
        "windows":[{"label":n,"start":s,"end":e} for n,s,e,c in WINDOWS],
        "window_definition":"90b OPTIX_POSITION_INTERVALS analysis windows; not exact regulatory boundaries",
        "counts":{m:n for m,start,n,rows in blocks},"pdf":str(out)},indent=2)+"\n")
    print(f"PDF\t{out}")
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--input-dir",type=Path,required=True)
    args=p.parse_args()
    values=pd.read_csv(args.input_dir/"frequency_plot_values.tsv.gz",sep="\t")
    render(values,args.input_dir)
    render(values,args.input_dir,zygosity=True)
if __name__=="__main__":main()
