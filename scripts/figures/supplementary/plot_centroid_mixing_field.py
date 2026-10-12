#!/usr/bin/env python3
"""Centroid mixing field: analytic geometry only; run on Pegasus."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, Normalize
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figure2"))
from common import set_publication_style

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project-root", default="/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling")
    args = ap.parse_args()
    root = Path(args.project_root)
    out = root / "analysis/arnaud_mixing_supplement"
    out.mkdir(parents=True, exist_ok=True)
    palette = root / "palettes/barfbow/luka_bleute.json"
    cmap = ListedColormap(json.loads(palette.read_text()), name="luka_bleute")
    set_publication_style(9)
    fig, ax = plt.subplots(figsize=(7.6,7.6))
    fig.subplots_adjust(left=.12,right=.82,bottom=.17,top=.90)
    grid = np.arange(10) * .05
    cl, cr = np.meshgrid(grid,grid)
    labels = ["0/0","0/1","1/1"]
    norm = Normalize(0,45)
    artist = None
    for i, left in enumerate([0.,.5,1.]):
        for j, right in enumerate([0.,.5,1.]):
            ax.annotate(labels[i]+" | "+labels[j],(left,right),xytext=(4,6),
                        textcoords="offset points",fontsize=8,clip_on=False)
            if left == right:
                ax.plot(left,right,"s",color="#777777",ms=6,clip_on=False,zorder=5)
                continue
            x = (1-cl)*left + cl*right
            y = cr*left + (1-cr)*right
            artist = ax.quiver(x,y,(x-left)/12,(y-right)/12,
                        (cl+cr)*50,cmap=cmap,norm=norm,angles="xy",
                        scale_units="xy",scale=1,pivot="mid",width=.0035,
                        headwidth=3.5,headlength=4.5,zorder=3)
            ax.plot(left,right,"o",mfc="white",mec="#333333",ms=6,zorder=5,clip_on=False)
    ax.set(xlim=(-.035,1.065),ylim=(-.035,1.065),aspect="equal",
           xticks=[0,.25,.5,.75,1],yticks=[0,.25,.5,.75,1],
           xlabel="Left expected ALT fraction",ylabel="Right expected ALT fraction")
    ax.spines["top"].set_visible(False);ax.spines["right"].set_visible(False)
    ax.set_title("Centroid shifts under reciprocal mixing",weight="bold",pad=14)
    cb = fig.colorbar(artist,cax=fig.add_axes([.87,.29,.025,.48]))
    cb.set_label("Mean of Left and Right mixing (%)")
    fig.text(.12,.11,"Independent Left/Right mixing: 0-45%, in 5-percentage-point steps.",fontsize=9)
    fig.text(.12,.075,"Arrow centre: candidate centroid; direction: shift from its unmixed origin.",fontsize=8)
    fig.text(.12,.045,"Arrow length: displacement / 12. Circles: unmixed; squares: stationary.",fontsize=8)
    fig.text(.12,.017,"Analytic geometry before sequencing error or bias; no simulations or refitting.",fontsize=8)
    fig.savefig(out / "04_centroid_mixing_parameter_field.pdf")
    plt.close(fig)
    (out / "04_centroid_mixing_parameter_field.legend.txt").write_text(
        "Centroid positions for nine Left/Right diploid genotype pairs. "
        "Left and Right reciprocal read mixing vary independently from 0 to 45% "
        "in five-percentage-point steps. Colour gives their arithmetic mean. "
        "Arrow centres give candidate centroids, directions give displacement from "
        "the unmixed origin, and lengths are displacement divided by 12. "
        "Arrows are not optimizer trajectories or a local vector flow. "
        "Centroids are calculated before sequencing error and allelic bias. "
        "The linear 0-45% colour scale uses the original Luka-Bleute palette; "
        "it represents mixing parameters rather than site density.\n")
    print(f"WROTE {out / '04_centroid_mixing_parameter_field.pdf'}")
if __name__ == "__main__":
    main()
