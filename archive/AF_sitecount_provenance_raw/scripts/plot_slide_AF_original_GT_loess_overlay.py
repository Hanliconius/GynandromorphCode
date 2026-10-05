#!/usr/bin/env python3

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42

# ============================================================
# Load the original script preamble exactly as written, but
# stop before it creates the two existing figures.
#
# This preserves:
#   - the cached histogram extraction
#   - all sample definitions
#   - the original LOESS helper
#   - colours, labels, CAP, CENTRES, OUT, etc.
# ============================================================

SOURCE = Path(
    "scripts/plot_slide_AF_original_GT_loess.py"
)

text = SOURCE.read_text()

marker = '\nfig1 = frame("Allele fractions in an unmixed control")'
cut = text.find(marker)

if cut < 0:
    raise RuntimeError(
        "Could not find the figure-creation boundary in "
        f"{SOURCE}"
    )

prefix = text[:cut]

ns = {
    "__file__": str(SOURCE),
    "__name__": "__overlay_loader__",
}

exec(compile(prefix, str(SOURCE), "exec"), ns)

hist = ns["hist"]
samples = ns["samples"]
COLOURS = ns["COLOURS"]
LABELS = ns["LABELS"]
CENTRES = ns["CENTRES"]
CAP = ns["CAP"]
OUT = ns["OUT"]
loess = ns["loess"]

if hist.shape != (3, 4, 50):
    raise RuntimeError(
        f"Unexpected histogram shape: {hist.shape}"
    )

control_counts = hist[0]
left_counts = hist[1]
right_counts = hist[2]

control_name = samples[0]

OUT.mkdir(parents=True, exist_ok=True)

# ============================================================
# Shared helpers
# ============================================================

GENOTYPE_INDICES = [0, 1, 2, 3]


def setup_axis(ax, label_text):
    ax.axvline(
        0.5,
        color="#555555",
        ls="--",
        lw=0.9,
        zorder=1,
    )
    ax.set_xlim(-0.015, 1.015)
    ax.set_ylim(0, CAP)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticks([0, 0.005, 0.010, 0.015])
    ax.yaxis.set_major_formatter(
        PercentFormatter(1, decimals=1)
    )
    ax.set_xlabel(
        "Observed ALT-read fraction",
        fontsize=12.5,
    )
    ax.set_ylabel(
        "Sites per 0.02 AF bin",
        fontsize=12,
    )
    ax.text(
        0.5,
        0.92,
        label_text,
        transform=ax.transAxes,
        ha="center",
        fontsize=11,
        bbox=dict(
            facecolor="white",
            edgecolor="none",
            alpha=0.85,
        ),
    )
    ax.grid(axis="y", color="#E5E5E5", lw=0.7)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=11)


def draw_profiles(
    ax,
    counts,
    *,
    linestyle="-",
    linewidth=2.0,
    alpha=0.90,
    zorder=3,
    show_clipped_labels=True,
    annotation_y_points=5,
    categories=GENOTYPE_INDICES,
):
    """
    Draw genotype-coloured LOESS profiles from one 4x50 histogram.
    Endpoint bins are kept unsmoothed, exactly as in the source.

    show_clipped_labels controls whether clipped-bin percentage
    labels are shown for bins above CAP.
    """
    counts = np.asarray(counts, dtype=float)
    n = float(counts.sum())

    if n <= 0:
        raise RuntimeError("Empty histogram provided.")

    fraction = counts / n

    for index in categories:
        colour = COLOURS[index]

        fitted = loess(
            CENTRES[1:-1],
            fraction[index, 1:-1],
        )

        ax.plot(
            CENTRES[1:-1],
            fitted,
            color=colour,
            alpha=alpha,
            linewidth=linewidth,
            linestyle=linestyle,
            zorder=zorder,
        )

        # Endpoint bins remain unsmoothed.
        for endpoint in (0, len(CENTRES) - 1):
            value = fraction[index, endpoint]

            if value <= 0:
                continue

            ax.plot(
                [CENTRES[endpoint], CENTRES[endpoint]],
                [0, min(value, CAP)],
                color=colour,
                alpha=alpha,
                linewidth=linewidth,
                linestyle=linestyle,
                zorder=zorder,
            )

        # If a bin exceeds the plotted cap, show a truncation marker
        # and, optionally, the actual clipped percentage.
        for bin_id in np.flatnonzero(fraction[index] > CAP):
            value = fraction[index, bin_id]

            ax.plot(
                CENTRES[bin_id],
                CAP * 0.98,
                "^",
                color=colour,
                markersize=5,
                zorder=zorder + 0.2,
                clip_on=False,
            )

            if show_clipped_labels:
                ax.annotate(
                    f"{100 * value:.2f}%",
                    (CENTRES[bin_id], CAP),
                    xytext=((index - 1.5) * 9, annotation_y_points),
                    textcoords="offset points",
                    rotation=90,
                    ha="center",
                    va="bottom",
                    fontsize=9,
                    color=colour,
                    annotation_clip=False,
                )


def add_panel_letter(ax, letter):
    ax.text(
        -0.14,
        1.08,
        letter,
        transform=ax.transAxes,
        fontsize=18,
        fontweight="bold",
        ha="left",
        va="bottom",
    )


# ============================================================
# Build the new triptych figure
# ============================================================

fig = plt.figure(figsize=(16.5, 7.8))

fig.text(
    0.055,
    0.94,
    "How allele mixing distorts wild-type and mutant classes",
    fontsize=28,
    weight="bold",
)

fig.text(
    0.055,
    0.885,
    "Reference control at left; dashed control overlays show how the gynandromorph deviates from an unmixed sample.",
    fontsize=15,
    color="#555555",
)

# Three science-style horizontal panels.
RECTS = [
    [0.065, 0.20, 0.255, 0.52],  # A
    [0.375, 0.20, 0.255, 0.52],  # B
    [0.685, 0.20, 0.255, 0.52],  # C
]

ax_a = fig.add_axes(RECTS[0])
ax_b = fig.add_axes(RECTS[1])
ax_c = fig.add_axes(RECTS[2])

# ------------------------------------------------------------
# A. Control reference
# ------------------------------------------------------------

add_panel_letter(ax_a, "A")
setup_axis(
    ax_a,
    f"{control_name}  |  n = {int(control_counts.sum()):,}",
)
draw_profiles(
    ax_a,
    control_counts,
    linestyle="-",
    linewidth=2.2,
    alpha=0.88,
    zorder=3,
    show_clipped_labels=True,
    annotation_y_points=5,
)

ax_a.set_title(
    "Unmixed control reference",
    fontsize=13,
    pad=14,
)

# ------------------------------------------------------------
# B. Left versus control
# ------------------------------------------------------------

add_panel_letter(ax_b, "B")
setup_axis(
    ax_b,
    f"Left  |  n = {int(left_counts.sum()):,}",
)

# Dashed control reference underneath.
draw_profiles(
    ax_b,
    control_counts,
    linestyle=(0, (4, 2)),
    linewidth=1.5,
    alpha=0.42,
    zorder=2,
    show_clipped_labels=False,
    annotation_y_points=-12,
)

# Solid gynandromorph left on top.
draw_profiles(
    ax_b,
    left_counts,
    linestyle="-",
    linewidth=2.2,
    alpha=0.92,
    zorder=4,
    show_clipped_labels=True,
    annotation_y_points=5,
)

ax_b.set_title(
    "Gynandromorph Left overlaid on control",
    fontsize=13,
    pad=14,
)

# ------------------------------------------------------------
# C. Right versus control
# ------------------------------------------------------------

add_panel_letter(ax_c, "C")
setup_axis(
    ax_c,
    f"Right  |  n = {int(right_counts.sum()):,}",
)

draw_profiles(
    ax_c,
    control_counts,
    linestyle=(0, (4, 2)),
    linewidth=1.5,
    alpha=0.42,
    zorder=2,
    show_clipped_labels=False,
    annotation_y_points=-12,
)

draw_profiles(
    ax_c,
    right_counts,
    linestyle="-",
    linewidth=2.2,
    alpha=0.92,
    zorder=4,
    show_clipped_labels=True,
    annotation_y_points=5,
)

ax_c.set_title(
    "Gynandromorph Right overlaid on control",
    fontsize=13,
    pad=14,
)

# ============================================================
# Legends
# ============================================================

genotype_handles = [
    Patch(facecolor=c, edgecolor="none", label=l)
    for c, l in zip(COLOURS, LABELS)
]

style_handles = [
    Line2D(
        [0], [0],
        color="#444444",
        lw=1.5,
        alpha=0.55,
        linestyle=(0, (4, 2)),
        label=f"Dashed = control ({control_name})",
    ),
    Line2D(
        [0], [0],
        color="#444444",
        lw=2.2,
        alpha=0.95,
        linestyle="-",
        label="Solid = focal sample",
    ),
]

fig.legend(
    handles=genotype_handles,
    loc="upper right",
    bbox_to_anchor=(0.965, 0.82),
    frameon=False,
    ncol=4,
    fontsize=11.5,
    title="Original VCF genotype",
    title_fontsize=12.5,
)

fig.legend(
    handles=style_handles,
    loc="upper left",
    bbox_to_anchor=(0.67, 0.84),
    frameon=False,
    ncol=1,
    fontsize=11,
    title="Line style",
    title_fontsize=12,
)

# ============================================================
# Caption
# ============================================================

fig.text(
    0.065,
    0.085,
    "LOESS curves of observed ALT-read fractions, coloured by original VCF genotype. "
    "Endpoint bins (0 and 1) are not smoothed; triangles and labels indicate clipped percentages above the plotted 1.5% cap. "
    "Panel A shows the unmixed control reference; panels B and C overlay the control (dashed) beneath the gynandromorph Left and Right distributions (solid).",
    fontsize=11,
    color="#555555",
)

# ============================================================
# Save
# ============================================================

stem = "control_vs_gynandromorph_original_GT_LOESS_overlay_triptych"

pdf_file = OUT / f"{stem}.pdf"
png_file = OUT / f"{stem}.png"

fig.savefig(pdf_file)
fig.savefig(png_file, dpi=200)

print(f"WROTE\t{pdf_file}", flush=True)
print(f"WROTE\t{png_file}", flush=True)
