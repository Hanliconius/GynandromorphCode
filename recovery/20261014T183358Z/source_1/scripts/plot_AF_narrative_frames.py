from pathlib import Path
from itertools import product
import csv
import gzip
import importlib.util
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from scipy.stats import betabinom

ROOT = Path("analysis/figures/pairwise_AF_cli")
PREFIX = ROOT / "TN_merged_control_sites_hex40_batlowW"
OUT = ROOT / "AF_narrative_frames"
OUT.mkdir(parents=True, exist_ok=True)

# Illustrative model settings, unchanged across every frame.
DEPTH = 40
CONCENTRATION = 22.
ERROR = .001
P_ALT = .5
K = np.arange(DEPTH + 1)
AF = K / DEPTH
STATES = list(product((0., .5, 1.), repeat=2))
BLUE = "#0072B2"
ORANGE = "#D55E00"
MODEL = "#F000B8"
HALO = [pe.Stroke(linewidth=3, foreground="white"), pe.Normal()]

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 12,
    "pdf.fonttype": 42,
})

# Load actual observations once, before writing any figures.
audit = json.loads(Path(str(PREFIX) + ".AF.audit.json").read_text())
with gzip.open(Path(str(PREFIX) + ".AF_bins.tsv.gz"), "rt") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

polygons, observed_counts = [], []
for row in rows:
    x, y = float(row["left_centre"]), float(row["right_centre"])
    w, h = float(row["width"]), float(row["height"])
    n = int(row["count"])
    if not np.isfinite([x, y, w, h]).all() or min(w, h) <= 0 or n < 0:
        raise RuntimeError("Invalid saved hexagon.")
    polygons.append([
        (x+w/2, y-h/4), (x+w/2, y+h/4), (x, y+h/2),
        (x-w/2, y+h/4), (x-w/2, y-h/4), (x, y-h/2),
    ])
    observed_counts.append(n)

polygons = np.asarray(polygons)
observed_counts = np.asarray(observed_counts, dtype=np.int64)
if not len(observed_counts) or not (observed_counts > 0).any():
    raise RuntimeError("No observed counts.")
if int(observed_counts.sum()) != audit["records_retained"]:
    raise RuntimeError("Observed bin total differs from audit.")

spec = importlib.util.spec_from_file_location(
    "af_tool", "pairwise-af-hex20/plot_pairwise_AF.py"
)
tool = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tool)
background_cmap = tool.colourmap("batlowW_r")

parameters = audit["parameters"]
vmin = parameters.get("vmin")
vmax = parameters.get("vmax")
vmin = 1 if vmin is None else float(vmin)
vmax = (
    max(float(observed_counts.max()), vmin * 10)
    if vmax is None else float(vmax)
)
background_norm = LogNorm(vmin=vmin, vmax=vmax)

def allele_probability(bits):
    return float(np.prod([P_ALT if b else 1-P_ALT for b in bits]))

def weights(founders):
    result = {state: 0. for state in STATES}
    for bits in product((0, 1), repeat=founders):
        if founders == 2:
            a, b = bits
            state = ((a+b)/2, (a+b)/2)
        elif founders == 3:
            a, b, c = bits
            state = ((a+b)/2, (a+c)/2)
        else:
            a, b, c, d = bits
            state = ((a+b)/2, (c+d)/2)
        result[state] += allele_probability(bits)
    assert np.isclose(sum(result.values()), 1)
    return result

WEIGHTS = {n: weights(n) for n in (2, 3, 4)}

def centroid(state, mixing):
    left, right = state
    return (
        (1-mixing)*left + mixing*right,
        mixing*left + (1-mixing)*right,
    )

def pmf(q):
    mean = ERROR + (1-2*ERROR)*q
    values = betabinom.pmf(
        K, DEPTH,
        mean*CONCENTRATION,
        (1-mean)*CONCENTRATION,
    )
    return values / values.sum()

def threshold(density, mass=.5):
    ordered = np.sort(density.ravel())[::-1]
    index = np.searchsorted(np.cumsum(ordered), mass)
    return ordered[min(index, len(ordered)-1)]

# Fixed frame dimensions and fixed positions; never use tight_layout
# or bbox_inches="tight", which would shift successive frames.
POSITIONS = [
    [.065, .27, .265, .471111],
    [.367, .27, .265, .471111],
    [.669, .27, .265, .471111],
]
TITLES = [
    "Two-founder\nSame two haplotypes",
    "Three-founder\nOne shared haplotype",
    "Four-founder\nFour independent haplotypes",
]

def canvas(title, subtitle):
    fig = plt.figure(figsize=(16, 9))
    fig.text(.065, .93, title, fontsize=24, weight="bold")
    fig.text(.065, .865, subtitle, fontsize=14, color="#444444")
    return fig

def footer(fig, extra):
    fig.text(.065, .15, extra, fontsize=12)
    fig.text(
        .065, .055,
        "Illustrative model: allele frequency 0.5; 40 reads/side; "
        "beta-binomial concentration 22; error 0.001; no bias.\n"
        "Founders refer to haplotype origins. Predictive contours describe "
        "read-fraction spread, not confidence in fitted centroids.",
        fontsize=10, color="#555555", linespacing=1.5,
    )

def style_af(axis):
    axis.set(
        xlim=(-.03, 1.03), ylim=(-.03, 1.03),
        aspect="equal",
        xlabel="Left ALT-read fraction",
        ylabel="Right ALT-read fraction",
    )
    axis.set_xticks([0, .25, .5, .75, 1])
    axis.set_yticks([0, .25, .5, .75, 1])
    axis.spines[["top", "right"]].set_visible(False)
    axis.plot([0, 1], [0, 1], color="#CCCCCC", lw=.7, zorder=1)

def background(axis):
    keep = observed_counts > 0
    collection = PolyCollection(
        polygons[keep],
        array=observed_counts[keep].astype(float),
        cmap=background_cmap, norm=background_norm,
        edgecolors="none", zorder=0,
    )
    axis.add_collection(collection)
    return collection

def model_frame(mixing, contours, show_data=False):
    percentage = round(100*mixing)
    title = (
        "The same model predictions over the observed gynandromorph"
        if show_data else
        f"Two samples under {percentage}% symmetric mixing"
    )
    subtitle = (
        f"Identical TN data in all three panels: "
        f"{observed_counts.sum():,} retained sites; no model-based site selection."
        if show_data else
        "Compare two, three and four founding haplotypes."
    )
    fig = canvas(title, subtitle)
    image = None

    for position, title, founders in zip(POSITIONS, TITLES, (2, 3, 4)):
        axis = fig.add_axes(position)
        style_af(axis)
        axis.set_title(title, fontsize=15, pad=15)

        if show_data:
            image = background(axis)

        for state, weight in WEIGHTS[founders].items():
            if weight == 0:
                continue
            ql, qr = centroid(state, mixing)

            if contours:
                density = np.outer(pmf(qr), pmf(ql))
                level = threshold(density)
                if density.min() < level < density.max():
                    result = axis.contour(
                        AF, AF, density, levels=[level],
                        colors=[MODEL], linewidths=1.3, zorder=4,
                    )
                    if hasattr(result, "set_path_effects"):
                        result.set_path_effects(HALO)
                    else:
                        for collection in result.collections:
                            collection.set_path_effects(HALO)

            if mixing > 0 and state[0] != state[1]:
                arrow = axis.annotate(
                    "", xy=(ql, qr), xytext=state,
                    arrowprops=dict(
                        arrowstyle="-|>", color="#333333",
                        lw=1.3, mutation_scale=13,
                        shrinkA=4, shrinkB=5,
                    ),
                    zorder=5,
                )
                arrow.arrow_patch.set_path_effects(HALO)
                axis.scatter(
                    *state, s=30, facecolors="white",
                    edgecolors="#333333", zorder=6,
                )

            axis.scatter(
                ql, qr, s=48, color=MODEL,
                edgecolors="white", linewidths=.8, zorder=7,
            )

    if show_data:
        cax = fig.add_axes([.951, .30, .010, .39])
        cb = fig.colorbar(image, cax=cax)
        cb.set_label("Observed sites / hexagon", fontsize=10)
        cb.ax.tick_params(labelsize=9)

    description = (
        "Filled dots: predicted centroids."
        if not contours else
        "Filled dots: predicted centroids   |   Rings: 50% predictive regions"
        "   |   Arrows: displacement from unmixed states"
    )
    footer(fig, description)
    return fig

# Histograms use full model probabilities, including boundary mass.
single = sum(
    probability * pmf(state[0])
    for state, probability in WEIGHTS[2].items()
)

mixed_left = np.zeros(DEPTH+1)
mixed_right = np.zeros(DEPTH+1)
for state, probability in WEIGHTS[4].items():
    ql, qr = centroid(state, .2)
    mixed_left += probability * pmf(ql)
    mixed_right += probability * pmf(qr)

hist_max = 1.08 * max(single.max(), mixed_left.max(), mixed_right.max())
hist_edges = np.r_[AF - .5/DEPTH, 1 + .5/DEPTH]

def histogram_frame(mixed=False):
    fig = canvas(
        "Mixing redistributes the allele-read fractions" if mixed else
        "A single unmixed diploid sample",
        "Illustrative marginal distributions from the same read-count model.",
    )
    axis = fig.add_axes([.11, .27, .56, .49])
    if mixed:
        axis.stairs(
            mixed_left, hist_edges, color=BLUE, lw=2,
            alpha=.8, label="Left: 20% mixing",
        )
        axis.stairs(
            mixed_right, hist_edges, color=ORANGE,
            linestyle="--", lw=2, alpha=.8, label="Right: 20% mixing",
        )
        note = (
            "Two independent diploid sources.\n"
            "Each library contains 20% from the other.\n\n"
            "Under these symmetric assumptions,\n"
            "Left and Right have identical marginal\n"
            "histograms, so the curves overlap."
        )
    else:
        axis.stairs(
            single, hist_edges, color=BLUE,
            fill=True, alpha=.65, label="Unmixed sample",
        )
        note = (
            "Diploid genotype fractions:\n"
            "0/0 → 0\n"
            "0/1 → 0.5\n"
            "1/1 → 1\n\n"
            "Read sampling and emission noise\n"
            "broaden these peaks."
        )
    axis.set(
        xlim=(-.03, 1.03), ylim=(0, hist_max),
        xlabel="ALT-read fraction",
        ylabel="Probability per read-fraction bin",
    )
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(frameon=False)
    fig.text(.72, .69, note, va="top", fontsize=14, linespacing=1.7)
    footer(fig, "All read-fraction outcomes retained; linear probability scale.")
    return fig

def save_frame(pdf, fig, name):
    fig.savefig(OUT / f"{name}.pdf")
    fig.savefig(OUT / f"{name}.png", dpi=160)
    pdf.savefig(fig)
    plt.close(fig)
    print(f"FRAME\t{name}", flush=True)

with PdfPages(OUT / "AF_narrative_all_frames.pdf") as pdf:
    save_frame(pdf, histogram_frame(False), "01_single_sample_histogram")
    save_frame(pdf, histogram_frame(True), "02_mixed_LR_histograms")
    save_frame(pdf, model_frame(0., False), "03_unmixed_founder_spaces")
    for number, mixing in enumerate((.05, .10, .15, .20), start=4):
        save_frame(
            pdf, model_frame(mixing, True),
            f"{number:02d}_mixing_{round(100*mixing):02d}pct",
        )
    save_frame(
        pdf, model_frame(.20, True, True),
        "08_mixing_20pct_observed_TN",
    )

print(f"\nOUTPUT_DIRECTORY\t{OUT.resolve()}")
print("PNG_DIMENSIONS\t2560 x 1440")
print("PDF_DIMENSIONS\t16 x 9 inches")
