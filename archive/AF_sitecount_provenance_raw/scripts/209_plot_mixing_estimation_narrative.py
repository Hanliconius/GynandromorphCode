#!/usr/bin/env python3
import argparse, importlib.util, json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm

ROOT=Path("analysis/genomewide_mixing_corrected_v1")
OUT=ROOT/"mixing_estimation_narrative_v2"
PDF=OUT/"figures/mixing_estimation_from_counts_to_simulations.pdf"
GRID=OUT/"tables/conditional_likelihood_grid.tsv.gz"
AUDIT=OUT/"mixing_estimation_narrative.audit.json"
NSITES,NGRID,SEED=25000,25,209001

def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def header(fig,n,title,subtitle):
    fig.text(.055,.965,f"STEP {n}",va="top",fontsize=9,
             weight="bold",color="#666666")
    fig.text(.055,.925,title,va="top",fontsize=19,
             weight="bold",color="#222222")
    fig.text(.055,.875,subtitle,va="top",fontsize=9.2,
             color="#666666")

def arrays(inputs):
    d=inputs["observed"]
    la=np.asarray(d["left_alt"],float)
    ln=np.asarray(d["left_total"],float)
    ra=np.asarray(d["right_alt"],float)
    rn=np.asarray(d["right_total"],float)
    keep=np.isfinite(la+ln+ra+rn)&(ln>0)&(rn>0)
    return (
        la[keep],ln[keep],ra[keep],rn[keep],
        la[keep]/ln[keep],ra[keep]/rn[keep],
    )

def page_locus(pdf, x, s):
    c_left = x["c_left"]
    c_right = x["c_right"]

    labels = [
        "0/0",
        "0/1",
        "1/1",
    ]

    fractions = [
        0.0,
        0.5,
        1.0,
    ]

    fig, axes = plt.subplots(
        3,
        3,
        figsize=(10.8, 8.2),
        sharex=True,
        sharey=True,
    )

    fig.subplots_adjust(
        left=0.10,
        right=0.96,
        bottom=0.10,
        top=0.69,
        wspace=0.17,
        hspace=0.25,
    )

    header(
        fig,
        2,
        "Reciprocal mixing moves every possible genotype pair",
        (
            "Rows are Left genotypes; columns are Right genotypes. "
            "Each arrow runs from the unmixed allele fractions to "
            "the fitted expectation."
        ),
    )

    equation = fig.add_axes([
        0.08,
        0.72,
        0.84,
        0.10,
    ])

    equation.axis("off")

    equation.text(
        0.5,
        0.66,
        r"$q_L=(1-c_L)p_L+c_Lp_R$"
        "                 "
        r"$q_R=(1-c_R)p_R+c_Rp_L$",
        ha="center",
        va="center",
        fontsize=17,
        color=s.DARK,
    )

    equation.text(
        0.5,
        0.12,
        (
            f"Fitted reciprocal mixing:  "
            f"cL = {c_left:.3f}   |   "
            f"cR = {c_right:.3f}"
        ),
        ha="center",
        va="center",
        fontsize=10,
        color=s.GREEN,
    )

    for row, (
        left_label,
        p_left,
    ) in enumerate(
        zip(
            labels,
            fractions,
        )
    ):
        for column, (
            right_label,
            p_right,
        ) in enumerate(
            zip(
                labels,
                fractions,
            )
        ):
            axis = axes[
                row,
                column,
            ]

            q_left = (
                (1.0 - c_left) * p_left
                + c_left * p_right
            )

            q_right = (
                (1.0 - c_right) * p_right
                + c_right * p_left
            )

            axis.plot(
                [0, 1],
                [0, 1],
                color="#E1E1E1",
                linewidth=0.7,
                zorder=0,
            )

            displacement = np.hypot(
                q_left - p_left,
                q_right - p_right,
            )

            if displacement > 1e-10:
                axis.annotate(
                    "",
                    xy=(
                        q_left,
                        q_right,
                    ),
                    xytext=(
                        p_left,
                        p_right,
                    ),
                    arrowprops=dict(
                        arrowstyle="-|>",
                        color=s.GREEN,
                        linewidth=1.6,
                        mutation_scale=10,
                    ),
                    zorder=3,
                )

            axis.scatter(
                [p_left],
                [p_right],
                s=34,
                facecolors="white",
                edgecolors=s.MID,
                linewidths=1.1,
                zorder=4,
            )

            axis.scatter(
                [q_left],
                [q_right],
                s=44,
                marker="D",
                color=s.ORANGE,
                edgecolors="white",
                linewidths=0.5,
                zorder=5,
            )

            axis.text(
                0.04,
                0.94,
                (
                    f"q = "
                    f"({q_left:.2f}, {q_right:.2f})"
                ),
                transform=axis.transAxes,
                ha="left",
                va="top",
                fontsize=7.2,
                color=s.DARK,
            )

            axis.set_xlim(
                -0.06,
                1.06,
            )

            axis.set_ylim(
                -0.06,
                1.06,
            )

            axis.set_aspect("equal")

            axis.set_xticks([
                0,
                0.5,
                1,
            ])

            axis.set_yticks([
                0,
                0.5,
                1,
            ])

            s.style_axis(
                axis
            )

            if row == 0:
                axis.set_title(
                    f"Right {right_label}",
                    fontsize=9,
                    pad=5,
                )

            if column == 0:
                axis.set_ylabel(
                    (
                        f"Left {left_label}\n"
                        "Right expected AF"
                    ),
                    fontsize=8,
                )

            if row == 2:
                axis.set_xlabel(
                    "Left expected AF",
                    fontsize=8,
                )

    fig.text(
        0.50,
        0.045,
        (
            "Open circle = unmixed (pL, pR); orange diamond = "
            "fitted mixed expectation (qL, qR). This is the complete "
            "diploid genotype-pair space."
        ),
        ha="center",
        fontsize=8.3,
        color=s.MID,
    )

    pdf.savefig(
        fig
    )

    plt.close(
        fig
    )


def page_ladder(pdf,x,s):
    cl=x["c_left"]
    cr=x["c_right"]

    settings=[
        (0.,0.,"No mixing"),
        (.1,.1,"10% each way"),
        (cl,cr,f"Fitted {100*cl:.1f}% / {100*cr:.1f}%"),
        (.3,.3,"30% each way"),
    ]

    fig,axs=plt.subplots(1,4,figsize=(13.5,4.3))
    fig.subplots_adjust(
        left=.06,right=.98,bottom=.17,top=.68,wspace=.30
    )

    header(
        fig,3,
        "Mixing moves the complete nine-centroid genotype pattern",
        "Opposing-homozygote centroids are diamonds; "
        "heterozygous combinations move simultaneously.",
    )

    for ax,(a,b,title) in zip(axs,settings):
        c=s.genotype_centroids(a,b)
        ax.scatter(
            c[:,2],c[:,3],s=45,facecolors="white",
            edgecolors=s.BLUE,lw=1.3,
        )

        opposing=(
            (c[:,0]!=c[:,1])
            & np.isin(c[:,0],[0,1])
            & np.isin(c[:,1],[0,1])
        )

        ax.scatter(
            c[opposing,2],c[opposing,3],
            s=55,marker="D",color=s.ORANGE,
            edgecolors="white",lw=.5,
        )
        ax.plot([0,1],[0,1],color="#DDDDDD",lw=.8)
        ax.set(
            xlim=(-.04,1.04),ylim=(-.04,1.04),
            xlabel="Expected Left AF",
            ylabel="Expected Right AF",
            title=title,
        )
        ax.set_aspect("equal")
        s.style_axis(ax)

    pdf.savefig(fig)
    plt.close(fig)

def page_candidates(pdf, x, s):
    *_, left_af, right_af = arrays(x)

    c_left = x["c_left"]
    c_right = x["c_right"]

    candidate_values = np.linspace(
        0.0,
        0.45,
        10,
    )

    candidate_rows = []

    for candidate_left in candidate_values:
        for candidate_right in candidate_values:
            table = s.genotype_centroids(
                candidate_left,
                candidate_right,
            )

            mean_mixing = 0.5 * (
                candidate_left
                + candidate_right
            )

            for row in table:
                candidate_rows.append((
                    row[2],
                    row[3],
                    mean_mixing,
                ))

    candidates = np.asarray(
        candidate_rows,
        dtype=float,
    )

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(12.9, 5.1),
    )

    fig.subplots_adjust(
        left=0.055,
        right=0.93,
        bottom=0.15,
        top=0.69,
        wspace=0.25,
    )

    header(
        fig,
        4,
        (
            "The fitted centroids emerge from the "
            "observed AF pattern"
        ),
        (
            "The observed million-site hexbin is unchanged: "
            "first naked, then with the candidate model space, "
            "then with only the fitted solution emphasized."
        ),
    )

    titles = (
        "Observed AF density only",
        "Candidate centroid space",
        (
            f"Fitted: {100*c_left:.1f}% Left / "
            f"{100*c_right:.1f}% Right"
        ),
    )

    for axis, title in zip(
        axes,
        titles,
    ):
        axis.hexbin(
            left_af,
            right_af,
            gridsize=70,
            extent=(
                0,
                1,
                0,
                1,
            ),
            mincnt=1,
            cmap=s.GREY_CMAP,
            norm=LogNorm(),
            linewidths=0,
        )

        axis.plot(
            [0, 1],
            [0, 1],
            color="#CCCCCC",
            linewidth=0.8,
        )

        axis.set(
            xlim=(-0.02, 1.02),
            ylim=(-0.02, 1.02),
            xlabel="Left ALT-read fraction",
            ylabel="Right ALT-read fraction",
            title=title,
        )

        axis.set_aspect(
            "equal"
        )

        s.style_axis(
            axis
        )

    candidate_artist = axes[1].scatter(
        candidates[:, 0],
        candidates[:, 1],
        c=candidates[:, 2],
        cmap="viridis",
        vmin=0.0,
        vmax=0.45,
        s=8,
        alpha=0.42,
        linewidths=0,
        zorder=5,
    )

    axes[1].text(
        0.03,
        0.04,
        (
            "cL and cR: 0–45%\n"
            "10 × 10 display grid"
        ),
        transform=axes[1].transAxes,
        fontsize=7.2,
        color=s.DARK,
        ha="left",
        va="bottom",
        bbox=dict(
            boxstyle="round,pad=0.25",
            facecolor="white",
            edgecolor="#CCCCCC",
            alpha=0.88,
        ),
        zorder=8,
    )

    axes[2].scatter(
        candidates[:, 0],
        candidates[:, 1],
        s=7,
        color="#BDBDBD",
        alpha=0.28,
        linewidths=0,
        zorder=4,
    )

    fitted = s.genotype_centroids(
        c_left,
        c_right,
    )

    axes[2].scatter(
        fitted[:, 2],
        fitted[:, 3],
        s=42,
        facecolors="white",
        edgecolors=s.BLUE,
        linewidths=1.3,
        label="Fitted nine centroids",
        zorder=6,
    )

    opposing = (
        (fitted[:, 0] != fitted[:, 1])
        & np.isin(
            fitted[:, 0],
            [0, 1],
        )
        & np.isin(
            fitted[:, 1],
            [0, 1],
        )
    )

    axes[2].scatter(
        fitted[
            opposing,
            2,
        ],
        fitted[
            opposing,
            3,
        ],
        s=50,
        marker="D",
        color=s.ORANGE,
        edgecolors="white",
        linewidths=0.5,
        zorder=7,
    )

    axes[2].legend(
        frameon=True,
        fontsize=7.2,
        loc="lower left",
    )

    color_axis = fig.add_axes([
        0.945,
        0.20,
        0.012,
        0.42,
    ])

    colorbar = fig.colorbar(
        candidate_artist,
        cax=color_axis,
    )

    colorbar.set_label(
        "Mean candidate mixing",
        fontsize=8,
    )

    colorbar.ax.tick_params(
        labelsize=7,
    )

    fig.text(
        0.055,
        0.065,
        (
            "Grey points in the final panel are the same candidate "
            "cloud shown in colour in the middle panel; coloured "
            "outlines mark the fitted nine-state pattern."
        ),
        ha="left",
        fontsize=8.2,
        color=s.MID,
    )

    pdf.savefig(
        fig
    )

    plt.close(
        fig
    )


def page_anchors(pdf,x,s):
    *_,lf,rf=arrays(x)
    cl=x["c_left"]
    cr=x["c_right"]

    da=np.hypot(lf-cl,rf-(1-cr))
    db=np.hypot(lf-(1-cl),rf-cr)
    orientation_a=da<=db
    keep=np.minimum(da,db)<=.12

    left=np.where(orientation_a,lf,1-lf)[keep]
    right=np.where(orientation_a,1-rf,rf)[keep]

    if len(left)==0 or len(right)==0:
        raise RuntimeError(
            "No opposing-homozygote anchor observations found"
        )

    fig,axs=plt.subplots(1,3,figsize=(12.8,5.7))
    fig.subplots_adjust(
        left=.07,right=.97,bottom=.16,top=.70,wspace=.30
    )

    header(
        fig,5,
        "The strongest genotype contrast gives an intuitive "
        "near-20% signal",
        "Post-fit anchor assignments explain the signal but do not "
        "replace the full likelihood estimator.",
    )

    axs[0].hexbin(
        lf[keep],rf[keep],
        gridsize=45,extent=(0,1,0,1),mincnt=1,
        cmap="Greys",norm=LogNorm(),linewidths=0,
    )

    c=s.genotype_centroids(cl,cr)
    axs[0].scatter(
        c[:,2],c[:,3],s=30,
        facecolors="white",edgecolors=s.BLUE,
    )
    axs[0].set(
        xlim=(-.02,1.02),ylim=(-.02,1.02),
        xlabel="Left ALT-read fraction",
        ylabel="Right ALT-read fraction",
        title=f"Opposing anchors (n={keep.sum():,})",
    )
    axs[0].set_aspect("equal")

    bins=np.linspace(0,.45,80)

    for ax,value,estimate,color,side in (
        (axs[1],left,cl,s.BLUE,"Left"),
        (axs[2],right,cr,s.ORANGE,"Right"),
    ):
        ax.hist(
            value,bins=bins,density=True,
            color=color,alpha=.70,
        )
        ax.axvline(estimate,color="black",lw=1.6)
        ax.set(
            xlabel=f"Per-site implied {side} mixing",
            ylabel="Density",
            title=f"{side} anchor median = {np.median(value):.3f}",
        )

    for ax in axs:
        s.style_axis(ax)

    fig.text(
        .5,.065,
        "Black lines are the official million-site MLEs; "
        "these histograms are explanatory diagnostics.",
        ha="center",fontsize=8.5,color=s.MID,
    )

    pdf.savefig(fig)
    plt.close(fig)

def subset(data,index):
    n=len(data["left_alt"])
    result={}

    for key,value in data.items():
        array=np.asarray(value)

        if array.ndim>0 and len(array)==n:
            result[key]=array[index]
        else:
            result[key]=value

    return result

def make_grid(x,m):
    fit=x["observed_fit"]["full"]["alternative"]
    epsilon=float(fit.get("epsilon",0.0))
    n=len(x["observed"]["left_alt"])

    index=np.sort(
        np.random.default_rng(SEED).choice(
            n,size=min(NSITES,n),replace=False
        )
    )

    data=subset(x["observed"],index)

    left_grid=np.linspace(
        max(0,x["c_left"]-.08),
        min(.45,x["c_left"]+.08),
        NGRID,
    )

    right_grid=np.linspace(
        max(0,x["c_right"]-.08),
        min(.45,x["c_right"]+.08),
        NGRID,
    )

    rows=[]

    for c_left in left_grid:
        for c_right in right_grid:
            states=m.model.state_log_likelihoods(
                data,
                float(c_left),
                float(c_right),
                x["emission"],
            )

            log_likelihood=float(
                m.model_evidence(
                    states,
                    epsilon,
                ).sum()
            )

            rows.append((
                c_left,c_right,epsilon,len(index),
                log_likelihood,
                x["c_left"],x["c_right"],SEED,
            ))

    table=pd.DataFrame(
        rows,
        columns=(
            "c_left c_right epsilon sites logL "
            "reference_left reference_right seed"
        ).split(),
    )

    table["relative_logL"]=(
        table["logL"]-table["logL"].max()
    )

    return table

def likelihood_grid(x,m):
    required=set(
        "c_left c_right epsilon sites logL reference_left "
        "reference_right seed relative_logL".split()
    )

    if GRID.is_file() and GRID.stat().st_size>0:
        table=pd.read_csv(GRID,sep="\t")

        reusable=(
            required<=set(table.columns)
            and len(table)==NGRID*NGRID
            and int(table["sites"].iloc[0])==NSITES
            and int(table["seed"].iloc[0])==SEED
            and np.isclose(
                table["reference_left"].iloc[0],
                x["c_left"],
            )
            and np.isclose(
                table["reference_right"].iloc[0],
                x["c_right"],
            )
        )

        if reusable:
            return table,True

    table=make_grid(x,m)
    GRID.parent.mkdir(parents=True,exist_ok=True)
    table.to_csv(GRID,sep="\t",index=False)

    return table,False

def page_likelihood(pdf,x,table,s):
    left_values=np.sort(table["c_left"].unique())
    right_values=np.sort(table["c_right"].unique())

    matrix=(
        table.pivot(
            index="c_right",
            columns="c_left",
            values="relative_logL",
        )
        .reindex(
            index=right_values,
            columns=left_values,
        )
        .to_numpy()
    )

    best=table.loc[table["logL"].idxmax()]

    nearest_right=right_values[
        np.argmin(abs(right_values-x["c_right"]))
    ]

    nearest_left=left_values[
        np.argmin(abs(left_values-x["c_left"]))
    ]

    left_slice=(
        table[np.isclose(table["c_right"],nearest_right)]
        .sort_values("c_left")
    )

    right_slice=(
        table[np.isclose(table["c_left"],nearest_left)]
        .sort_values("c_right")
    )

    fig=plt.figure(figsize=(12.7,8.2))

    header(
        fig,6,
        "The numeric estimate comes from the likelihood maximum",
        "The real restored-model objective is evaluated on a "
        "deterministic 25,000-site subset; epsilon is fixed at "
        "its full-fit value.",
    )

    ax=fig.add_axes([.07,.13,.43,.62])

    image=ax.imshow(
        matrix,origin="lower",aspect="auto",
        extent=(
            left_values.min(),left_values.max(),
            right_values.min(),right_values.max(),
        ),
        cmap="viridis",
    )

    ax.scatter(
        [x["c_left"]],[x["c_right"]],
        marker="*",s=180,color="white",
        edgecolors="black",lw=.8,
        label="Official 1M-site MLE",
    )

    ax.scatter(
        [best["c_left"]],[best["c_right"]],
        marker="x",s=70,color="#F0E442",lw=2,
        label="25k-site grid maximum",
    )

    ax.set(
        xlabel="Candidate Left mixing cL",
        ylabel="Candidate Right mixing cR",
        title="Conditional log-likelihood surface",
    )
    ax.legend(fontsize=7.5,loc="lower right")

    color_axis=fig.add_axes([.515,.18,.014,.52])
    fig.colorbar(
        image,cax=color_axis
    ).set_label(
        "Log likelihood relative to maximum"
    )

    slices=(
        (
            [.62,.47,.32,.27],
            left_slice,
            "c_left",
            s.BLUE,
            "Slice near fitted cR",
        ),
        (
            [.62,.13,.32,.27],
            right_slice,
            "c_right",
            s.ORANGE,
            "Slice near fitted cL",
        ),
    )

    for position,data,key,color,title in slices:
        z=fig.add_axes(position)
        z.plot(
            data[key],data["relative_logL"],
            color=color,lw=2,
        )
        z.axvline(
            x[key],color="black",ls="--"
        )
        z.set(
            xlabel=f"Candidate {key.replace('_',' ')}",
            ylabel="Relative log likelihood",
            title=title,
        )
        s.style_axis(z,grid=True)

    fig.text(
        .07,.055,
        "The subset verifies the likelihood shape. Reported values "
        "remain the authoritative full one-million-site MLEs.",
        fontsize=8.5,color=s.MID,
    )

    pdf.savefig(fig)
    plt.close(fig)

def page_handoff(pdf,x,s):
    table=x["simulations"].copy()

    table=table[
        table["true_c_left"].between(.14,.26)
        & table["true_c_right"].between(.14,.26)
    ]

    if table.empty:
        raise RuntimeError(
            "No focused simulation rows were available"
        )

    fig=plt.figure(figsize=(12.8,7.3))

    header(
        fig,7,
        "The fitted values become inputs to observed-matched "
        "simulations",
        "Simulations validate recovery; they do not generate "
        "the original mixing estimates.",
    )

    ax=fig.add_axes([.055,.54,.89,.25])
    ax.axis("off")

    from matplotlib.patches import FancyBboxPatch

    items=[
        (0.,"Observed AD counts",
         "1,000,000 sites",s.BLUE),
        (.26,"Full likelihood fit",
         "latent genotype states",s.PURPLE),
        (.52,"Official MLE",
         f"Left {100*x['c_left']:.2f}%\n"
         f"Right {100*x['c_right']:.2f}%",s.GREEN),
        (.78,"Simulation inputs",
         "known cL and cR\nrecovery checked",s.ORANGE),
    ]

    for px,title,body,color in items:
        ax.add_patch(
            FancyBboxPatch(
                (px,.18),.20,.60,
                boxstyle=(
                    "round,pad=.018,"
                    "rounding_size=.02"
                ),
                facecolor="white",
                edgecolor=color,
                lw=1.5,
                transform=ax.transAxes,
            )
        )
        ax.text(
            px+.10,.61,title,
            ha="center",fontsize=10,
            weight="bold",color=color,
        )
        ax.text(
            px+.10,.35,body,
            ha="center",fontsize=8.5,
            linespacing=1.4,
        )

    for start,end in ((.21,.25),(.47,.51),(.73,.77)):
        ax.annotate(
            "",xy=(end,.48),xytext=(start,.48),
            xycoords="axes fraction",
            arrowprops=dict(
                arrowstyle="->",
                color="#666666",
                lw=1.4,
            ),
        )

    panels=(
        ([.08,.10,.37,.32],"left",s.BLUE),
        ([.57,.10,.37,.32],"right",s.ORANGE),
    )

    for position,side,color in panels:
        z=fig.add_axes(position)
        z.scatter(
            table[f"true_c_{side}"],
            table[f"estimated_c_{side}_alt"],
            s=8,alpha=.18,color=color,lw=0,
        )
        z.plot(
            [.13,.27],[.13,.27],
            color="black",lw=1,
        )
        z.axhline(
            x[f"c_{side}"],
            color=s.GREEN,ls="--",lw=1.4,
        )
        z.set(
            xlim=(.13,.27),ylim=(.13,.27),
            xlabel=f"Known simulated {side} mixing",
            ylabel=f"Recovered {side} mixing",
            title=(
                f"{side.capitalize()} recovery "
                "near observed value"
            ),
        )
        s.style_axis(z)

    pdf.savefig(fig)
    plt.close(fig)

def validate(x,m):
    required=(
        "left_alt left_total right_alt "
        "right_total controls"
    ).split()

    missing=[
        key for key in required
        if key not in x["observed"]
    ]

    if missing:
        raise RuntimeError(
            f"Observed data missing {missing}"
        )

    observed_sites=len(
        x["observed"]["left_alt"]
    )

    if observed_sites!=1_000_000:
        raise RuntimeError(
            "Expected exactly 1,000,000 sites; "
            f"found {observed_sites}"
        )

    test_data=subset(
        x["observed"],
        np.arange(100),
    )

    states=m.model.state_log_likelihoods(
        test_data,
        x["c_left"],
        x["c_right"],
        x["emission"],
    )

    epsilon=float(
        x["observed_fit"]["full"]["alternative"]
        .get("epsilon",0.0)
    )

    evidence=m.model_evidence(
        states,
        epsilon,
    )

    if (
        evidence.shape!=(100,)
        or not np.isfinite(evidence).all()
    ):
        raise RuntimeError(
            "Likelihood smoke test failed"
        )

    print("VALIDATION_ONLY=PASS")
    print(f"observed_sites={observed_sites}")
    print(f"c_left={x['c_left']:.8f}")
    print(f"c_right={x['c_right']:.8f}")
    print("likelihood_smoke_sites=100")

def main(validate_only=False):
    stage208=load(
        "scripts/"
        "208_plot_mixing_number_mechanism.py",
        "stage208_for_209",
    )

    stage185=load(
        "scripts/185_restored_model_pilot.py",
        "stage185_for_209",
    )

    inputs=stage208.load_inputs(
        load_counts=True
    )

    validate(inputs,stage185)

    if validate_only:
        return

    PDF.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    GRID.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    mpl.rcParams.update({
        "font.family":"DejaVu Sans",
        "axes.titleweight":"bold",
        "figure.facecolor":"white",
    })

    likelihood,reused=likelihood_grid(
        inputs,
        stage185,
    )

    with PdfPages(PDF) as pdf:
        stage208.page_inputs(pdf,inputs)
        page_locus(pdf,inputs,stage208)
        page_ladder(pdf,inputs,stage208)
        page_candidates(pdf,inputs,stage208)
        page_anchors(pdf,inputs,stage208)
        page_likelihood(
            pdf,inputs,likelihood,stage208
        )
        page_handoff(pdf,inputs,stage208)

    best=likelihood.loc[
        likelihood["logL"].idxmax()
    ]

    audit={
        "status":"PASS",
        "pages":7,
        "official_full_data_c_left":
            inputs["c_left"],
        "official_full_data_c_right":
            inputs["c_right"],
        "observed_sites":1_000_000,
        "conditional_grid_sites":
            int(best["sites"]),
        "conditional_grid_best_c_left":
            float(best["c_left"]),
        "conditional_grid_best_c_right":
            float(best["c_right"]),
        "conditional_grid_reused":reused,
        "likelihood":(
            "stage185 state_log_likelihoods "
            "+ model_evidence"
        ),
        "focal_GT_used":False,
        "model_refitted":False,
        "fundamental_method_changed":False,
        "pdf":str(PDF),
        "grid":str(GRID),
    }

    AUDIT.write_text(
        json.dumps(audit,indent=2)+"\n"
    )

    for path in (PDF,GRID,AUDIT):
        if (
            not path.is_file()
            or path.stat().st_size==0
        ):
            raise RuntimeError(
                f"Missing output {path}"
            )

    print(
        "MIXING ESTIMATION NARRATIVE COMPLETE"
    )
    print(f"PDF={PDF.resolve()}")
    print(f"AUDIT={AUDIT.resolve()}")

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument(
        "--validate-only",
        action="store_true",
    )
    arguments=parser.parse_args()
    main(arguments.validate_only)
