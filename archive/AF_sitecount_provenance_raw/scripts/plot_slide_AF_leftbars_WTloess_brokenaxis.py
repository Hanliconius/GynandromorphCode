from pathlib import Path
import json
import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter


ROOT = Path("analysis/figures/pairwise_AF_cli")
OUT = ROOT / "AF_narrative_frames/original_GT_histograms"
OUT.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------
# Input audit information
# ---------------------------------------------------------------------

control_audit = json.loads((
    ROOT / "control_benchmarks/"
    "best_coverage_control_self_hex45_batlowW.AF.audit.json"
).read_text())

target_audit = json.loads((
    ROOT /
    "TN_merged_control_sites_hex40_batlowW.AF.audit.json"
).read_text())


# Keep the sample definitions exactly as in the previous plotting work.
cp = control_audit
tp = target_audit


# ---------------------------------------------------------------------
# Histogram cache produced by the original-GT histogram analysis
#
# sample order:
#   0 = wild-type/control
#   1 = gynandromorph Left
#   2 = gynandromorph Right
#
# genotype order:
#   0 = 0/0
#   1 = 0/1
#   2 = 1/1
#   3 = missing/other
# ---------------------------------------------------------------------

CACHE = OUT / "original_GT_counts.npz"

if not CACHE.exists():
    raise FileNotFoundError(
        f"Missing histogram cache:\n{CACHE}\n\n"
        "Run scripts/plot_slide_AF_original_GT_loess.py first."
    )

dat = np.load(CACHE)

# Allow for the key names used by the existing cache.
if "hist" in dat:
    hist = dat["hist"]
elif "counts" in dat:
    hist = dat["counts"]
else:
    raise RuntimeError(
        f"Could not find hist/counts array in {CACHE}. "
        f"Available keys: {list(dat.keys())}"
    )

if "edges" in dat:
    EDGES = dat["edges"]
else:
    EDGES = np.linspace(0.0, 1.0, 51)

CENTRES = 0.5 * (EDGES[:-1] + EDGES[1:])


# ---------------------------------------------------------------------
# Colours retained from the original figure
# ---------------------------------------------------------------------

COLOURS = [
    "#2C9BCB",  # 0/0
    "#F28E2B",  # 0/1
    "#B61F72",  # 1/1
    "#AAAAAA",  # missing
]

LABELS = [
    "0/0",
    "0/1",
    "1/1",
    "missing",
]

PLOT_GT = [0, 1, 2]


# ---------------------------------------------------------------------
# Simple LOESS-style local polynomial smoother
#
# This reproduces the local smoothing logic used for the previous
# plotting version without requiring an additional package.
# ---------------------------------------------------------------------

def loess(x, y, span=0.20):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    n = len(x)

    if n < 3:
        return y.copy()

    k = max(3, int(np.ceil(span * n)))

    out = np.full(n, np.nan, dtype=float)

    for i, x0 in enumerate(x):
        distance = np.abs(x - x0)

        order = np.argsort(distance)
        use = order[:k]

        xu = x[use]
        yu = y[use]

        dmax = distance[use].max()

        if dmax == 0:
            out[i] = np.mean(yu)
            continue

        u = distance[use] / dmax

        # tricube weights
        w = (1.0 - u**3)**3

        # locally weighted linear fit centred on x0
        X = np.column_stack([
            np.ones(len(xu)),
            xu - x0,
        ])

        sw = np.sqrt(w)

        Xw = X * sw[:, None]
        yw = yu * sw

        try:
            beta, _, _, _ = np.linalg.lstsq(
                Xw,
                yw,
                rcond=None,
            )

            out[i] = beta[0]

        except np.linalg.LinAlgError:
            out[i] = np.average(
                yu,
                weights=w,
            )

    return out


# ---------------------------------------------------------------------
# Convert counts to fractions
# ---------------------------------------------------------------------

def genotype_fractions(counts):
    counts = np.asarray(counts, dtype=float)

    total = counts.sum()

    if total <= 0:
        return np.zeros_like(
            counts,
            dtype=float,
        )

    return counts / total


# ---------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------

BAR_WIDTH = (EDGES[1] - EDGES[0]) * 0.92

LOW_Y_MAX = 0.015


def draw_left_bars(ax, fraction):
    """
    Left gynandromorph shown as stacked genotype-specific histogram bars.
    """

    bottom = np.zeros(
        len(CENTRES),
        dtype=float,
    )

    for gt_i in PLOT_GT:
        vals = fraction[gt_i]

        ax.bar(
            CENTRES,
            vals,
            width=BAR_WIDTH,
            bottom=bottom,
            color=COLOURS[gt_i],
            edgecolor="none",
            alpha=0.67,
            align="center",
            zorder=2,
        )

        bottom += vals


def draw_wildtype_loess(ax, fraction):
    """
    Wild-type/control shown as genotype-specific LOESS curves.

    The extreme AF=0 and AF=1 bins are not included in the smoother.
    Their actual observed fractions are retained as points.
    """

    x = CENTRES[1:-1]

    for gt_i in PLOT_GT:
        y = fraction[gt_i, 1:-1]

        smoothed = loess(
            x,
            y,
            span=0.20,
        )

        ax.plot(
            x,
            smoothed,
            color=COLOURS[gt_i],
            linewidth=2.3,
            alpha=1.0,
            zorder=6,
        )

        # Actual endpoint observations
        ax.scatter(
            [CENTRES[0], CENTRES[-1]],
            [
                fraction[gt_i, 0],
                fraction[gt_i, -1],
            ],
            s=31,
            facecolor=COLOURS[gt_i],
            edgecolor="white",
            linewidth=0.7,
            zorder=8,
        )


def draw_panel(
    ax_low,
    ax_high,
    left_counts,
    wt_counts,
):
    left_fraction = genotype_fractions(
        left_counts
    )

    wt_fraction = genotype_fractions(
        wt_counts
    )

    # Draw the same objects on each axis.
    # Different ylim values create the genuine axis break.
    for ax in (ax_low, ax_high):

        draw_left_bars(
            ax,
            left_fraction,
        )

        draw_wildtype_loess(
            ax,
            wt_fraction,
        )

        ax.axvline(
            0.5,
            color="#555555",
            linestyle="--",
            linewidth=0.9,
            alpha=0.75,
            zorder=1,
        )

        ax.set_xlim(
            0.0,
            1.0,
        )

        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)

        ax.tick_params(
            axis="both",
            direction="out",
            width=0.8,
            length=3,
            labelsize=9,
        )


    # -------------------------------------------------------------
    # Work out where the upper section should begin.
    # -------------------------------------------------------------

    left_stack = left_fraction[
        PLOT_GT,
        :
    ].sum(axis=0)

    high_values = []

    high_values.extend(
        left_stack[
            left_stack > LOW_Y_MAX
        ].tolist()
    )

    for gt_i in PLOT_GT:

        left_endpoint = wt_fraction[
            gt_i,
            0
        ]

        right_endpoint = wt_fraction[
            gt_i,
            -1
        ]

        if left_endpoint > LOW_Y_MAX:
            high_values.append(
                float(left_endpoint)
            )

        if right_endpoint > LOW_Y_MAX:
            high_values.append(
                float(right_endpoint)
            )


    if len(high_values) > 0:

        high_values = np.asarray(
            high_values,
            dtype=float,
        )

        smallest_high = float(
            np.min(high_values)
        )

        largest_high = float(
            np.max(high_values)
        )

        high_y_min = max(
            LOW_Y_MAX * 1.7,
            smallest_high * 0.88,
        )

        high_y_max = (
            largest_high * 1.08
        )

        if high_y_max <= high_y_min:
            high_y_max = (
                high_y_min * 1.20
            )

    else:

        high_y_min = (
            LOW_Y_MAX * 2.0
        )

        high_y_max = (
            LOW_Y_MAX * 4.0
        )


    # -------------------------------------------------------------
    # Axis ranges
    # -------------------------------------------------------------

    ax_low.set_ylim(
        0.0,
        LOW_Y_MAX,
    )

    ax_high.set_ylim(
        high_y_min,
        high_y_max,
    )


    # -------------------------------------------------------------
    # Percentage formatting
    # -------------------------------------------------------------

    ax_low.yaxis.set_major_formatter(
        PercentFormatter(
            xmax=1.0,
            decimals=1,
        )
    )

    ax_high.yaxis.set_major_formatter(
        PercentFormatter(
            xmax=1.0,
            decimals=0,
        )
    )


    # -------------------------------------------------------------
    # Remove duplicated x axis from upper region
    # -------------------------------------------------------------

    ax_high.tick_params(
        axis="x",
        which="both",
        bottom=False,
        labelbottom=False,
    )

    ax_high.spines[
        "bottom"
    ].set_visible(False)

    ax_low.spines[
        "top"
    ].set_visible(False)


    # -------------------------------------------------------------
    # Genuine broken-axis diagonal marks
    # -------------------------------------------------------------

    d = 0.010

    kwargs = dict(
        color="black",
        clip_on=False,
        linewidth=0.9,
    )

    ax_high.plot(
        (-d, +d),
        (-d, +d),
        transform=ax_high.transAxes,
        **kwargs,
    )

    ax_high.plot(
        (1 - d, 1 + d),
        (-d, +d),
        transform=ax_high.transAxes,
        **kwargs,
    )

    ax_low.plot(
        (-d, +d),
        (1 - d, 1 + d),
        transform=ax_low.transAxes,
        **kwargs,
    )

    ax_low.plot(
        (1 - d, 1 + d),
        (1 - d, 1 + d),
        transform=ax_low.transAxes,
        **kwargs,
    )


# ---------------------------------------------------------------------
# Extract samples
#
# Existing ordering from original_GT_counts.npz:
#
#   0 = wild type / control
#   1 = Left gynandromorph
#   2 = Right gynandromorph
#
# We deliberately omit Right here.
# ---------------------------------------------------------------------

if hist.ndim != 3:
    raise RuntimeError(
        f"Expected histogram array with three dimensions; "
        f"found shape {hist.shape}"
    )

if hist.shape[0] < 2:
    raise RuntimeError(
        f"Expected at least WT + Left samples; "
        f"found histogram shape {hist.shape}"
    )

wt_counts = hist[0]
left_counts = hist[1]


# ---------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------

fig = plt.figure(
    figsize=(7.2, 5.7),
)

gs = fig.add_gridspec(
    nrows=2,
    ncols=1,
    height_ratios=[
        1.05,
        3.9,
    ],
    hspace=0.055,
)

ax_high = fig.add_subplot(
    gs[0]
)

ax_low = fig.add_subplot(
    gs[1],
    sharex=ax_high,
)


draw_panel(
    ax_low=ax_low,
    ax_high=ax_high,
    left_counts=left_counts,
    wt_counts=wt_counts,
)


# ---------------------------------------------------------------------
# Labels
# ---------------------------------------------------------------------

ax_low.set_xlabel(
    "Observed ALT-allele fraction",
    fontsize=11,
)

ax_low.set_ylabel(
    "Fraction of SNPs per AF bin",
    fontsize=11,
)

# Position label so it visually spans the broken y-axis.
ax_low.yaxis.set_label_coords(
    -0.095,
    0.70,
)


fig.suptitle(
    "Observed allele-fraction distributions",
    x=0.12,
    y=0.985,
    ha="left",
    va="top",
    fontsize=14,
    fontweight="bold",
)


# ---------------------------------------------------------------------
# Legend
#
# Colour = original genotype
# Geometry = sample
# ---------------------------------------------------------------------

legend_handles = [
    Patch(
        facecolor="#777777",
        alpha=0.67,
        edgecolor="none",
        label="Left gynandromorph: histogram",
    ),

    Line2D(
        [0],
        [0],
        color="#333333",
        linewidth=2.3,
        label="Wild type: LOESS",
    ),

    Line2D(
        [0],
        [0],
        marker="s",
        linestyle="none",
        markersize=7,
        markerfacecolor=COLOURS[0],
        markeredgecolor="none",
        label="0/0",
    ),

    Line2D(
        [0],
        [0],
        marker="s",
        linestyle="none",
        markersize=7,
        markerfacecolor=COLOURS[1],
        markeredgecolor="none",
        label="0/1",
    ),

    Line2D(
        [0],
        [0],
        marker="s",
        linestyle="none",
        markersize=7,
        markerfacecolor=COLOURS[2],
        markeredgecolor="none",
        label="1/1",
    ),
]


fig.legend(
    handles=legend_handles,
    loc="upper center",
    bbox_to_anchor=(
        0.54,
        0.925,
    ),
    ncol=5,
    frameon=False,
    fontsize=8,
    handlelength=2.0,
    columnspacing=1.25,
)


fig.subplots_adjust(
    left=0.14,
    right=0.98,
    bottom=0.12,
    top=0.80,
)


# ---------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------

pdf_file = (
    OUT /
    "03_left_gyn_bars_WT_LOESS_broken_y.pdf"
)

png_file = (
    OUT /
    "03_left_gyn_bars_WT_LOESS_broken_y.png"
)


fig.savefig(
    pdf_file,
    bbox_inches="tight",
)

fig.savefig(
    png_file,
    dpi=300,
    bbox_inches="tight",
)

plt.close(fig)


print(
    f"Wrote: {pdf_file}"
)

print(
    f"Wrote: {png_file}"
)
