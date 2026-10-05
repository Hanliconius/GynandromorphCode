from pathlib import Path
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter


OUT = Path(
    "analysis/figures/pairwise_AF_cli/"
    "AF_narrative_frames/original_GT_histograms"
)

CACHE = OUT / "ERR2298208_vs_TNLeft_original_GT_counts.npz"

if not CACHE.exists():
    raise FileNotFoundError(f"Missing cache: {CACHE}")

dat = np.load(CACHE)

print("CACHE:", CACHE)
print("KEYS:", list(dat.keys()))

if "hist" in dat:
    hist = dat["hist"]
elif "counts" in dat:
    hist = dat["counts"]
else:
    raise RuntimeError(
        f"Could not find hist/counts in cache. Keys: {list(dat.keys())}"
    )

if "edges" in dat:
    edges = dat["edges"]
else:
    edges = np.linspace(0.0, 1.0, hist.shape[-1] + 1)

centres = 0.5 * (edges[:-1] + edges[1:])

print("HIST SHAPE:", hist.shape)


# ------------------------------------------------------------
# Cache ordering from the ERR2298208 fork:
#   0 = ERR2298208
#   1 = TN Left
#
# Genotypes:
#   0 = 0/0
#   1 = 0/1
#   2 = 1/1
# ------------------------------------------------------------

wt = hist[0, :3, :].astype(float)
left = hist[1, :3, :].astype(float)


def fractions(x):
    total = x.sum()
    if total == 0:
        return np.zeros_like(x)
    return x / total


wt = fractions(wt)
left = fractions(left)


# ------------------------------------------------------------
# Colours
# ------------------------------------------------------------

colors = {
    0: "#2C9BCB",  # 0/0
    1: "#F28E2B",  # 0/1
    2: "#B61F72",  # 1/1
}

labels = {
    0: "0/0",
    1: "0/1",
    2: "1/1",
}


# ------------------------------------------------------------
# Simple smoother for ERR2298208
# ------------------------------------------------------------

def smooth(x, y, window=4):
    out = np.zeros_like(y, dtype=float)

    for i in range(len(y)):
        lo = max(0, i - window)
        hi = min(len(y), i + window + 1)

        xx = x[lo:hi]
        yy = y[lo:hi]

        if len(xx) < 2:
            out[i] = y[i]
        else:
            coef = np.polyfit(xx, yy, 1)
            out[i] = np.polyval(coef, x[i])

    return out


# ------------------------------------------------------------
# Broken-axis limits
# ------------------------------------------------------------

LOW_MAX = 0.015

high_values = []

for gt in range(3):
    vals = left[gt]
    high_values.extend(vals[vals > LOW_MAX].tolist())

    for idx in (0, -1):
        if wt[gt, idx] > LOW_MAX:
            high_values.append(float(wt[gt, idx]))

if len(high_values) == 0:
    high_min = 0.03
    high_max = 0.10
else:
    high_values = np.asarray(high_values)
    high_min = max(LOW_MAX * 1.7, high_values.min() * 0.88)
    high_max = high_values.max() * 1.08


# ------------------------------------------------------------
# Figure
# ------------------------------------------------------------

fig = plt.figure(figsize=(7.4, 5.8))

gs = fig.add_gridspec(
    2,
    1,
    height_ratios=[1.0, 3.8],
    hspace=0.055
)

ax_high = fig.add_subplot(gs[0])
ax_low = fig.add_subplot(gs[1], sharex=ax_high)


# ------------------------------------------------------------
# Grouped bars for TN Left
#
# Each AF bin gets three narrow bars:
#   blue   = 0/0
#   magenta= 1/1
#   orange = 0/1
#
# Orange is drawn last / visually foremost.
# ------------------------------------------------------------

bin_width = edges[1] - edges[0]

group_width = bin_width * 0.92
bar_width = group_width / 3.0

offsets = {
    0: -bar_width,
    2:  0.0,
    1: +bar_width,
}

plot_order = [0, 2, 1]


def draw_grouped_bars(ax):
    for gt in plot_order:
        ax.bar(
            centres + offsets[gt],
            left[gt],
            width=bar_width * 0.99,
            color=colors[gt],
            edgecolor="none",
            alpha=0.50,
            align="center",
            zorder=3 if gt == 1 else 2,
        )


# ------------------------------------------------------------
# ERR2298208 as LOESS-like curves
# ------------------------------------------------------------

def draw_wt_curves(ax):
    x_mid = centres[1:-1]

    for gt in plot_order:
        y_mid = wt[gt, 1:-1]

        ys_mid = smooth(
            x_mid,
            y_mid,
            window=6,
        )

        # Build smooth shoulders that connect the interior curve
        # continuously to the exact endpoint values at AF=0 and AF=1.

        # LEFT: exact AF=0 value -> first smoothed interior value
        x_left = np.linspace(
            0.0,
            x_mid[0],
            80
        )

        t_left = x_left / x_mid[0]

        # Steep near the endpoint, flattening toward the interior curve.
        s_left = 1.0 - (1.0 - t_left) ** 4

        y_left = (
            wt[gt, 0]
            + (ys_mid[0] - wt[gt, 0]) * s_left
        )

        # RIGHT: last smoothed interior value -> exact AF=1 value
        x_right = np.linspace(
            x_mid[-1],
            1.0,
            80
        )

        t_right = (
            (x_right - x_mid[-1])
            / (1.0 - x_mid[-1])
        )

        # Stay near the interior curve at first, then rise steeply
        # toward the exact endpoint value.
        s_right = t_right ** 4

        y_right = (
            ys_mid[-1]
            + (wt[gt, -1] - ys_mid[-1]) * s_right
        )

        x_full = np.concatenate([
            x_left[:-1],
            x_mid,
            x_right[1:],
        ])

        y_full = np.concatenate([
            y_left[:-1],
            ys_mid,
            y_right[1:],
        ])

        # Black outline / underlay for contrast against the bars
        ax.plot(
            x_full,
            y_full,
            color="black",
            linewidth=3.6 if gt == 1 else 3.1,
            alpha=0.95,
            zorder=8 if gt == 1 else 7,
        )

        # Coloured WT line on top
        ax.plot(
            x_full,
            y_full,
            color=colors[gt],
            linewidth=2.5 if gt == 1 else 2.1,
            alpha=1.0,
            zorder=9 if gt == 1 else 8,
        )


# ------------------------------------------------------------
# Draw both broken-axis panels
# ------------------------------------------------------------

for ax in (ax_low, ax_high):

    draw_grouped_bars(ax)
    draw_wt_curves(ax)

    ax.axvline(
        0.5,
        color="#555555",
        linestyle="--",
        linewidth=0.9,
        alpha=0.75,
        zorder=1,
    )

    ax.set_xlim(0.0, 1.0)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    ax.tick_params(
        axis="both",
        direction="out",
        width=0.8,
        length=3,
        labelsize=9,
    )


# ------------------------------------------------------------
# Broken y axis
# ------------------------------------------------------------

ax_low.set_ylim(0.0, LOW_MAX)
ax_high.set_ylim(high_min, high_max)

ax_low.yaxis.set_major_formatter(
    PercentFormatter(xmax=1.0, decimals=1)
)

ax_high.yaxis.set_major_formatter(
    PercentFormatter(xmax=1.0, decimals=0)
)

ax_high.tick_params(
    axis="x",
    which="both",
    bottom=False,
    labelbottom=False,
)

ax_high.spines["bottom"].set_visible(False)
ax_low.spines["top"].set_visible(False)


# diagonal break marks
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
    **kwargs
)

ax_high.plot(
    (1 - d, 1 + d),
    (-d, +d),
    transform=ax_high.transAxes,
    **kwargs
)

ax_low.plot(
    (-d, +d),
    (1 - d, 1 + d),
    transform=ax_low.transAxes,
    **kwargs
)

ax_low.plot(
    (1 - d, 1 + d),
    (1 - d, 1 + d),
    transform=ax_low.transAxes,
    **kwargs
)


# ------------------------------------------------------------
# Labels
# ------------------------------------------------------------

ax_low.set_xlabel(
    "Observed ALT-allele fraction",
    fontsize=11,
)

ax_low.set_ylabel(
    "Fraction of SNPs per AF bin",
    fontsize=11,
)

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


# ------------------------------------------------------------
# Legend
# ------------------------------------------------------------

legend_handles = [
    Patch(
        facecolor="#777777",
        edgecolor="none",
        alpha=0.50,
        label="TN Left: grouped bars",
    ),

    Line2D(
        [0], [0],
        color="#333333",
        linewidth=2.2,
        label="ERR2298208: smoothed curve",
    ),

    Patch(
        facecolor=colors[0],
        edgecolor="none",
        label="0/0",
    ),

    Patch(
        facecolor=colors[1],
        edgecolor="none",
        label="0/1",
    ),

    Patch(
        facecolor=colors[2],
        edgecolor="none",
        label="1/1",
    ),
]

fig.legend(
    handles=legend_handles,
    loc="upper center",
    bbox_to_anchor=(0.54, 0.925),
    ncol=5,
    frameon=False,
    fontsize=8,
    columnspacing=1.2,
)


fig.subplots_adjust(
    left=0.14,
    right=0.98,
    bottom=0.12,
    top=0.80,
)


# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

pdf_file = OUT / "ERR2298208_vs_TNLeft_groupedbars_LOESS_broken_y.pdf"
png_file = OUT / "ERR2298208_vs_TNLeft_groupedbars_LOESS_broken_y.png"

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

print("WROTE:", pdf_file)
print("WROTE:", png_file)
