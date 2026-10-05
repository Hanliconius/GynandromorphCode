from pathlib import Path
import ast
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe
from matplotlib.collections import PolyCollection
from matplotlib.colors import ListedColormap, LogNorm
from matplotlib.ticker import PercentFormatter
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "analysis/figures/pairwise_AF_cli"
OUT = BASE / "TN_main_figure"
OUT.mkdir(parents=True, exist_ok=True)

BINS = BASE / "TN_merged_control_sites_hex40_batlowW.AF_bins.tsv.gz"
HIST_ROOT = BASE / "AF_narrative_frames/original_GT_histograms"
COUNTS = HIST_ROOT / "original_GT_counts.npz"
HIST_AUDIT = HIST_ROOT / "original_GT_counts.audit.json"
PALETTE = ROOT / "palettes/barfbow/luka_bleute.json"

ORIGINAL_MODEL = ROOT / "scripts/plot_simulated_founder_slides.py"
CORRECTED_MODEL = ROOT / "scripts/plot_founder_slides_empirical_depth.py"
LOESS_SOURCE = ROOT / "scripts/plot_slide_AF_original_GT_loess.py"

for path in (
    BINS, COUNTS, HIST_AUDIT, PALETTE,
    ORIGINAL_MODEL, CORRECTED_MODEL, LOESS_SOURCE,
):
    if not path.is_file():
        raise FileNotFoundError(path)

# Reuse the existing LOESS implementation.
tree = ast.parse(LOESS_SOURCE.read_text())
functions = [
    node for node in tree.body
    if isinstance(node, ast.FunctionDef) and node.name == "loess"
]
if len(functions) != 1:
    raise RuntimeError("Expected exactly one existing loess function.")
loess_ns = {"np": np}
exec(
    compile(ast.Module(body=functions, type_ignores=[]),
            str(LOESS_SOURCE), "exec"),
    loess_ns,
)
loess = loess_ns["loess"]

# Load definitions without executing either simulation script.
model_ns = {
    "__name__": "model_overlay_definitions",
    "__file__": str(ORIGINAL_MODEL),
}
for path in (ORIGINAL_MODEL, CORRECTED_MODEL):
    tree = ast.parse(path.read_text())
    for node in tree.body:
        if isinstance(node, (ast.Import, ast.ImportFrom, ast.FunctionDef)):
            exec(
                compile(ast.Module(body=[node], type_ignores=[]),
                        str(path), "exec"),
                model_ns,
            )
        elif isinstance(node, ast.Assign):
            try:
                value = ast.literal_eval(node.value)
            except (ValueError, TypeError, SyntaxError):
                continue
            for target in node.targets:
                if isinstance(target, ast.Name):
                    model_ns[target.id] = value

MODEL_COLOUR = "#E600A9"
CENTROID_SIZE = 14

# Full-strength colour is retained for centroid points.
# Only the predictive-region overlay is made translucent.
MODEL_OVERLAY_COLOUR = mcolors.to_rgba(
    MODEL_COLOUR,
    0.42,
)

model_ns["MODEL_COLOUR"] = MODEL_OVERLAY_COLOUR

# Thin, partially transparent white separation from the heatmap.
# This is deliberately much lighter than the earlier 3.2-pt halo.
model_ns["HALO"] = [
    pe.Stroke(
        linewidth=1.8,
        foreground=(1.0, 1.0, 1.0, 0.70),
    ),
    pe.Normal(),
]

model_ns["ns"] = model_ns

if "predictive_region" not in model_ns:
    raise RuntimeError("Existing predictive_region function not found.")

_original_predictive_region = model_ns["predictive_region"]

def predictive_region_without_diagonal_corner_boxes(
    ax,
    q_left,
    q_right,
):
    ql = float(q_left)
    qr = float(q_right)

    is_00 = abs(ql) < 1e-12 and abs(qr) < 1e-12
    is_11 = abs(ql - 1.0) < 1e-12 and abs(qr - 1.0) < 1e-12

    if is_00 or is_11:
        return

    return _original_predictive_region(
        ax,
        q_left,
        q_right,
    )

model_ns["predictive_region"] = (
    predictive_region_without_diagonal_corner_boxes
)

# Read and check the two independently cached representations.
bins = pd.read_csv(BINS, sep="\t")
required = {
    "left_centre", "right_centre", "width", "height", "count"
}
if not required.issubset(bins.columns):
    raise RuntimeError(f"Hex table missing columns: {required - set(bins.columns)}")

with np.load(COUNTS, allow_pickle=False) as cache:
    hist = cache["hist"].copy()

metadata = json.loads(HIST_AUDIT.read_text())
samples = metadata["signature"]["samples"]

if hist.shape != (3, 4, 50):
    raise RuntimeError(f"Unexpected histogram shape: {hist.shape}")
if len(set(samples)) != len(samples):
    raise RuntimeError("Duplicate sample names in histogram audit.")

left_index = samples.index("TN_Left_merged")
right_index = samples.index("TN_Right_merged")
left_counts = hist[left_index]
right_counts = hist[right_index]
n = int(bins["count"].sum())

if n != int(left_counts.sum()) or n != int(right_counts.sum()):
    raise RuntimeError(
        "Hex plot and marginal histograms have different site totals: "
        f"hex={n}, left={left_counts.sum()}, right={right_counts.sum()}"
    )

colours = json.loads(PALETTE.read_text())
if not isinstance(colours, list) or len(colours) < 2:
    raise RuntimeError("Expected palette JSON containing a list of colours.")
cmap = ListedColormap(colours)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 7,
    "axes.labelsize": 8,
    "axes.linewidth": 0.6,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.facecolor": "white",
})

# Fixed dimensions: square main plot and narrow aligned marginals.
fig = plt.figure(figsize=(5.3, 5.3))
ax = fig.add_axes([0.275, 0.275, 0.60, 0.60])
bottom = fig.add_axes([0.275, 0.885, 0.60, 0.065])
side = fig.add_axes([0.120, 0.275, 0.065, 0.60])
cax = fig.add_axes([0.905, 0.405, 0.022, 0.34])

# Reconstruct the original hexagons directly from the saved geometry.
occupied = bins.loc[bins["count"] > 0].copy()
unit_hex = np.array([
    [0.0, -0.5],
    [0.5, -0.25],
    [0.5, 0.25],
    [0.0, 0.5],
    [-0.5, 0.25],
    [-0.5, -0.25],
])
centres = occupied[["left_centre", "right_centre"]].to_numpy()
sizes = occupied[["width", "height"]].to_numpy()
vertices = centres[:, None, :] + unit_hex[None, :, :] * sizes[:, None, :]

collection = PolyCollection(
    vertices,
    array=occupied["count"].to_numpy(dtype=float),
    cmap=cmap,
    norm=LogNorm(vmin=1, vmax=float(occupied["count"].max())),
    edgecolors="none",
    linewidths=0,
    rasterized=False,
)
ax.add_collection(collection)

# Preserve complete edge hexagons.
padding = max(
    float(bins["width"].max()),
    float(bins["height"].max()),
) / 2 + 0.003
limits = (-padding, 1 + padding)
ax.set(xlim=limits, ylim=limits, aspect="equal")
ax.plot([0, 1], [0, 1], color="white", lw=0.5, alpha=0.45)

# Same four-founder, symmetric 20% overlay as the earlier cartoons.
mixing = 0.20
error = 0.001
centroid_rows = []
for p_left in (0.0, 0.5, 1.0):
    for p_right in (0.0, 0.5, 1.0):
        q_left = error + (1 - 2 * error) * (
            (1 - mixing) * p_left + mixing * p_right
        )
        q_right = error + (1 - 2 * error) * (
            mixing * p_left + (1 - mixing) * p_right
        )

        previous_artists = {id(artist) for artist in ax.get_children()}
        model_ns["predictive_region"](ax, q_left, q_right)

        # Keep the existing contour geometry; refine only its styling.
        for artist in ax.get_children():
            if id(artist) not in previous_artists:
                if hasattr(artist, "set_alpha"):
                    artist.set_alpha(0.80)
                if hasattr(artist, "set_linewidth"):
                    artist.set_linewidth(0.70)

        ax.scatter(
            [q_left], [q_right],
            s=CENTROID_SIZE, color=MODEL_COLOUR,
            alpha=0.68, edgecolors="white",
            linewidths=0.45, zorder=10,
        )
        centroid_rows.append((p_left, p_right, q_left, q_right))

ticks = [0, 0.25, 0.5, 0.75, 1]
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_xlabel("Left allele fraction", labelpad=3)
ax.set_ylabel("Right allele fraction", labelpad=3)
ax.tick_params(length=2.5, pad=2)
for spine in ("top", "right"):
    ax.spines[spine].set_visible(False)

cb = fig.colorbar(collection, cax=cax)
cb.set_label("Sites per hexagon", fontsize=7, labelpad=5)
cb.ax.tick_params(labelsize=6, length=2)
cb.outline.set_linewidth(0.5)

# Marginals: all genotype categories combined, same 0.02 AF bins.
# Smooth interior bins only; retain endpoint bins as observed.
x = (np.arange(50) + 0.5) / 50
cap = 0.015

def marginal(axis, counts, vertical=False):
    fraction = counts.astype(float) / counts.sum()
    smooth = loess(x[1:-1], fraction[1:-1])
    colour = "#555555"

    if vertical:
        axis.fill_betweenx(
            x[1:-1], 0, smooth, color=colour, alpha=0.17
        )
        axis.plot(smooth, x[1:-1], color=colour, lw=0.9)
        for index in (0, 49):
            value = min(fraction[index], cap)
            axis.plot([0, value], [x[index], x[index]],
                      color=colour, lw=0.9)
            if fraction[index] > cap:
                axis.plot(value, x[index], marker="<", ms=3,
                          color=colour, clip_on=False)
                axis.annotate(
                    f"{fraction[index]:.1%}",
                    (value, x[index]), xytext=(-3, 0),
                    textcoords="offset points",
                    ha="right", va="center", fontsize=6,
                    clip_on=False,
                )
        axis.set(ylim=limits, xlim=(cap * 1.05, 0))
        axis.set_xticks([0, cap])
        axis.xaxis.set_major_formatter(PercentFormatter(1, decimals=1))
        axis.set_yticks([])
        axis.set_xlabel("Sites / bin", fontsize=6, labelpad=2)
    else:
        axis.fill_between(
            x[1:-1], 0, smooth, color=colour, alpha=0.17
        )
        axis.plot(x[1:-1], smooth, color=colour, lw=0.9)
        for index in (0, 49):
            value = min(fraction[index], cap)
            axis.plot([x[index], x[index]], [0, value],
                      color=colour, lw=0.9)
            if fraction[index] > cap:
                axis.plot(x[index], value, marker="v", ms=3,
                          color=colour, clip_on=False)
                axis.annotate(
                    f"{fraction[index]:.1%}",
                    (x[index], value), xytext=(0, -4),
                    textcoords="offset points",
                    ha="center", va="top", fontsize=6,
                    clip_on=False,
                )
        axis.set(xlim=limits, ylim=(cap * 1.05, 0))
        axis.set_yticks([0, cap])
        axis.yaxis.set_major_formatter(PercentFormatter(1, decimals=1))
        axis.set_xticks([])
        axis.set_ylabel("Sites / bin", fontsize=6, labelpad=2)

    axis.tick_params(labelsize=6, length=2, pad=2)
    for spine in axis.spines.values():
        spine.set_visible(False)

GT_COLOURS = (
    "#2C9BCB",  # 0/0 REF
    "#F28E2B",  # 0/1 heterozygote
    "#B61F72",  # 1/1 ALT
)


def coloured_marginal(
    axis,
    counts_by_genotype,
    vertical=False,
):
    """
    Genotype-split marginal AF density using the same LOESS
    logic as the existing marginal() function, while restoring
    endpoint truncation markers at 0 and 1.

    Missing/non-diploid rows are omitted from this compact
    publication annotation.
    """

    counts_by_genotype = np.asarray(
        counts_by_genotype,
        dtype=float,
    )

    if (
        counts_by_genotype.ndim != 2
        or counts_by_genotype.shape[0] < 3
    ):
        raise RuntimeError(
            "Expected genotype-by-AF-bin histogram; "
            f"observed shape {counts_by_genotype.shape}"
        )

    callable_counts = counts_by_genotype[:3]

    total = float(
        callable_counts.sum()
    )

    if total <= 0:
        raise RuntimeError(
            "Marginal genotype histogram is empty."
        )

    # Let the existing helper set up the publication styling.
    pre_collection_ids = {
        id(artist)
        for artist in axis.collections
    }

    pre_line_ids = {
        id(artist)
        for artist in axis.lines
    }

    marginal(
        axis,
        callable_counts.sum(axis=0),
        vertical=vertical,
    )

    # Remove only the neutral-grey artists produced by the
    # temporary summed marginal.
    for artist in list(axis.collections):
        if id(artist) not in pre_collection_ids:
            artist.remove()

    for artist in list(axis.lines):
        if id(artist) not in pre_line_ids:
            artist.remove()

    upper_limit = (
        axis.get_xlim()[1]
        if vertical
        else axis.get_ylim()[1]
    )

    endpoint_marker_level = 0.985 * upper_limit
    endpoint_trigger = 0.92 * upper_limit

    for genotype_index, colour in enumerate(
        GT_COLOURS
    ):
        fraction = (
            callable_counts[genotype_index]
            / total
        )

        smooth = loess(
            x[1:-1],
            fraction[1:-1],
        )

        if vertical:
            axis.fill_betweenx(
                x[1:-1],
                0,
                smooth,
                color=colour,
                alpha=0.15,
                linewidth=0,
                zorder=2,
            )

            axis.plot(
                smooth,
                x[1:-1],
                color=colour,
                linewidth=0.85,
                alpha=0.95,
                zorder=3,
            )

        else:
            axis.fill_between(
                x[1:-1],
                0,
                smooth,
                color=colour,
                alpha=0.15,
                linewidth=0,
                zorder=2,
            )

            axis.plot(
                x[1:-1],
                smooth,
                color=colour,
                linewidth=0.85,
                alpha=0.95,
                zorder=3,
            )

        # Endpoint truncation markers: show when the mass at
        # AF = 0 or AF = 1 is very high relative to the plotted
        # y/x limit, so the reader can see that the curve is cut.
        endpoint0 = float(fraction[0])
        endpoint1 = float(fraction[-1])

        if vertical:
            if endpoint0 >= endpoint_trigger:
                axis.scatter(
                    endpoint_marker_level,
                    0.0,
                    marker=">",
                    s=22,
                    color=colour,
                    linewidths=0,
                    clip_on=False,
                    zorder=5,
                )

            if endpoint1 >= endpoint_trigger:
                axis.scatter(
                    endpoint_marker_level,
                    1.0,
                    marker=">",
                    s=22,
                    color=colour,
                    linewidths=0,
                    clip_on=False,
                    zorder=5,
                )

        else:
            if endpoint0 >= endpoint_trigger:
                axis.scatter(
                    0.0,
                    endpoint_marker_level,
                    marker="^",
                    s=22,
                    color=colour,
                    linewidths=0,
                    clip_on=False,
                    zorder=5,
                )

            if endpoint1 >= endpoint_trigger:
                axis.scatter(
                    1.0,
                    endpoint_marker_level,
                    marker="^",
                    s=22,
                    color=colour,
                    linewidths=0,
                    clip_on=False,
                    zorder=5,
                )




coloured_marginal(
    bottom,
    left_counts,
)

coloured_marginal(
    side,
    right_counts,
    vertical=True,
)

fig.text(0.275, 0.905, f"{n:,} sites", fontsize=7)
fig.text(
    0.275, 0.035,
    "20% mixing model  •  Centroids and 50% predictive regions",
    fontsize=6.5,
)

prefix = OUT / "TN_AF_with_marginal_LOESS"
fig.savefig(str(prefix) + ".pdf", dpi=600)
fig.savefig(str(prefix) + ".png", dpi=600)
plt.close(fig)

pd.DataFrame(
    centroid_rows,
    columns=["unmixed_left", "unmixed_right", "expected_left", "expected_right"],
).to_csv(OUT / "centroids.tsv", sep="\t", index=False)

audit = {
    "hex_input": str(BINS),
    "histogram_input": str(COUNTS),
    "histogram_audit": str(HIST_AUDIT),
    "samples": [samples[left_index], samples[right_index]],
    "sites": n,
    "palette": str(PALETTE),
    "marginals": "LOESS of total AF counts, all original-GT categories combined",
    "marginal_bin_width": 0.02,
    "marginal_display_cap": cap,
    "endpoints": "unsmoothed; clipped labels give actual endpoint-bin percentages",
    "model": "four-founder illustrative overlay, symmetric 20% read mixing",
    "predictive_region_source": str(CORRECTED_MODEL),
    "predictive_region_note": (
        "Existing 50% predictive contour implementation; "
        "not confidence intervals for fitted parameters"
    ),
    "new_fitting_performed": False,
}
(OUT / "figure.audit.json").write_text(json.dumps(audit, indent=2) + "\n")

print("PDF\t" + str(prefix) + ".pdf")
print("PNG\t" + str(prefix) + ".png")
