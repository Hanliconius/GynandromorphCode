from pathlib import Path
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Keep text editable/vector in PDF and SVG output.
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"

from matplotlib.colors import ListedColormap, LogNorm
from matplotlib.patches import FancyBboxPatch
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path("analysis/figures/pairwise_AF_cli/AF_narrative_frames")
OUT = ROOT / "supplementary_haplotype_scenarios"
OUT.mkdir(parents=True, exist_ok=True)

FILES = {
    "standard": ROOT / "simulated_founder_empirical_depth/simulated_AF.npz",
    "AB_A": ROOT / "AB_A_digametic/AB_A_simulated_AF.npz",
    "AA_BC": ROOT / "gyno_three_founder/gyno_simulated_AF.npz",
}
PALETTE = Path("palettes/barfbow/luka_bleute.json")

for path in [*FILES.values(), PALETTE]:
    if not path.is_file():
        raise FileNotFoundError(f"Required simulation or palette missing: {path}")

palette = json.loads(PALETTE.read_text())
if not isinstance(palette, list):
    raise ValueError("Expected a list of palette colours.")
cmap = ListedColormap(palette, name="luka_bleute")

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.linewidth": 0.7,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


# Load the corrected analytical contour code without running simulations.
import ast

original_path = Path("scripts/plot_simulated_founder_slides.py")
corrected_path = Path("scripts/plot_founder_slides_empirical_depth.py")

original_tree = ast.parse(original_path.read_text())
nodes = []
for node in original_tree.body:
    if isinstance(node, ast.Assign) and any(
        isinstance(t, ast.Name) and t.id == "rng"
        for t in node.targets
    ):
        break
    if isinstance(node, (ast.Import, ast.ImportFrom, ast.Assign)):
        nodes.append(node)
nodes.extend(
    node for node in original_tree.body
    if isinstance(node, ast.FunctionDef)
)

overlay_ns = {"__name__": "supplementary_overlay"}
exec(
    compile(ast.Module(body=nodes, type_ignores=[]),
            str(original_path), "exec"),
    overlay_ns,
)

corrected_source = corrected_path.read_text()
corrected_tree = ast.parse(corrected_source)
functions = [
    node for node in corrected_tree.body
    if isinstance(node, ast.FunctionDef)
]
matches = [node for node in functions if node.name == "predictive_region"]
if len(matches) != 1:
    raise RuntimeError("Expected one corrected predictive_region function.")

overlay_source = ast.get_source_segment(corrected_source, matches[0])
if "histogram2d" in overlay_source or "region_cache" in overlay_source:
    raise RuntimeError("Found the old binned contours, not the corrected overlay.")

nodes = [
    node for node in corrected_tree.body
    if isinstance(node, (ast.Import, ast.ImportFrom))
]
for node in corrected_tree.body:
    if isinstance(node, ast.Assign):
        try:
            ast.literal_eval(node.value)
        except (ValueError, TypeError, SyntaxError):
            continue
        nodes.append(node)
nodes.extend(functions)

exec(
    compile(ast.Module(body=nodes, type_ignores=[]),
            str(corrected_path), "exec"),
    overlay_ns,
)
overlay_ns["MODEL_COLOUR"] = "#F000B8"
overlay_ns["ns"] = overlay_ns

# Check all expected centres before generating the final figure.
temporary, test_ax = plt.subplots()
try:
    for mixing in (0.0, 0.2):
        for left in (0.0, 0.5, 1.0):
            for right in (0.0, 0.5, 1.0):
                ql = 0.001 + 0.998 * ((1-mixing)*left + mixing*right)
                qr = 0.001 + 0.998 * (mixing*left + (1-mixing)*right)
                overlay_ns["predictive_region"](test_ax, ql, qr)
finally:
    plt.close(temporary)


MODELS = [
    {
        "name": "AB/AB",
        "title": "Two haplotypes",
        "description": "Same two haplotypes",
        "left": "AB", "right": "AB",
        "states": [(0., 0.), (.5, .5), (1., 1.)],
    },
    {
        "name": "AB/A",
        "title": "Two haplotypes",
        "description": "Diploid / haploid",
        "left": "AB", "right": "A",
        "states": [(0., 0.), (.5, 0.), (.5, 1.), (1., 1.)],
    },
    {
        "name": "AB/AC",
        "title": "Three haplotypes",
        "description": "One shared",
        "left": "AB", "right": "AC",
        "states": [
            (l, r) for l in (0., .5, 1.) for r in (0., .5, 1.)
            if (l, r) not in [(0., 1.), (1., 0.)]
        ],
    },
    {
        "name": "AA/BC",
        "title": "Three haplotypes",
        "description": "One duplicated",
        "left": "AA", "right": "BC",
        "states": [(l, r) for l in (0., 1.) for r in (0., .5, 1.)],
    },
    {
        "name": "AB/CD",
        "title": "Four haplotypes",
        "description": "No shared haplotype",
        "left": "AB", "right": "CD",
        "states": [(l, r) for l in (0., .5, 1.) for r in (0., .5, 1.)],
    },
]

data = {}
with (
    np.load(FILES["standard"], allow_pickle=False) as standard,
    np.load(FILES["AB_A"], allow_pickle=False) as haploid,
    np.load(FILES["AA_BC"], allow_pickle=False) as duplicated,
):
    def pair(cache, prefix):
        return (
            cache[prefix + "_left"].copy(),
            cache[prefix + "_right"].copy(),
        )

    for mixing in (0, 20):
        data[mixing] = [
            pair(standard, f"mix{mixing}_model2"),
            pair(haploid, f"mix{mixing}"),
            pair(standard, f"mix{mixing}_model3"),
            pair(duplicated, f"mix{mixing}"),
            pair(standard, f"mix{mixing}_model4"),
        ]

sizes = set()
for panels in data.values():
    for left, right in panels:
        if left.ndim != 1 or left.shape != right.shape or not len(left):
            raise ValueError("Invalid simulation array dimensions.")
        for values in (left, right):
            if (
                not np.isfinite(values).all()
                or ((values < 0) | (values > 1)).any()
            ):
                raise ValueError("Invalid simulated AF values.")
        sizes.add(len(left))

if len(sizes) != 1:
    raise ValueError("Panel site counts differ; count colours would not match.")
n = sizes.pop()

BINS = 40
ERROR = 0.001
HAP_COLOURS = {
    "A": "#0072B2",
    "B": "#E69F00",
    "C": "#009E73",
    "D": "#CC79A7",
}

# Fixed physical dimensions keep every AF plot square.
# Five model columns: A-E.
FIG_W, FIG_H = 13.8, 7.3

LEFT = 0.105
RIGHT = 0.925
COLUMN_GAP = 0.015

WIDTH = (RIGHT - LEFT - 4 * COLUMN_GAP) / 5
HEIGHT = WIDTH * FIG_W / FIG_H

XS = [
    LEFT + i * (WIDTH + COLUMN_GAP)
    for i in range(5)
]

ROW_Y = {0: 0.410, 20: 0.085}

fig = plt.figure(figsize=(FIG_W, FIG_H))
collections = []

def cell(ax, centre, haplotypes, label):
    width, height = 0.34, 0.39
    bottom = 0.15
    outline = FancyBboxPatch(
        (centre - width / 2, bottom), width, height,
        boxstyle="round,pad=0.015,rounding_size=0.07",
        facecolor="#F7F7F7", edgecolor="#777777", linewidth=0.8,
    )
    ax.add_patch(outline)
    positions = (
        [centre] if len(haplotypes) == 1
        else [centre - 0.075, centre + 0.075]
    )
    for x, hap in zip(positions, haplotypes):
        bar = FancyBboxPatch(
            (x - 0.026, bottom + 0.09), 0.052, 0.22,
            boxstyle="round,pad=0,rounding_size=0.02",
            facecolor=HAP_COLOURS[hap], edgecolor="none",
        )
        ax.add_patch(bar)
        ax.text(
            x, bottom + 0.20, hap,
            ha="center", va="center",
            fontsize=8, weight="bold", color="white",
        )
    ax.text(centre, 0.065, label, ha="center", va="center", fontsize=8)

for column, model in enumerate(MODELS):
    cartoon = fig.add_axes([XS[column], 0.755, WIDTH, 0.205])
    cartoon.set_xlim(0, 1)
    cartoon.set_ylim(0, 1)
    cartoon.axis("off")
    cartoon.text(
        0.0, 1.02, chr(ord("A") + column),
        fontsize=12, weight="bold", va="top",
    )
    cartoon.text(
        .5, .98, model["name"].replace("/", " / "),
        fontsize=11, weight="bold", ha="center", va="top",
    )
    cartoon.text(
        .5, .79, model["title"],
        fontsize=8.5, ha="center", va="top",
    )
    cartoon.text(
        .5, .65, model["description"],
        fontsize=8, ha="center", va="top", color="#555555",
    )
    cell(cartoon, .26, model["left"], "Left")
    cell(cartoon, .74, model["right"], "Right")

for mixing in (0, 20):
    for column, model in enumerate(MODELS):
        ax = fig.add_axes([XS[column], ROW_Y[mixing], WIDTH, HEIGHT])
        left, right = data[mixing][column]
        collection = ax.hexbin(
            left, right,
            gridsize=BINS, extent=(0, 1, 0, 1),
            mincnt=1, cmap=cmap,
            linewidths=0, rasterized=False,
        )
        collections.append(collection)

        # Symmetric read mixing, followed by the simulation's error mapping.
        c = mixing / 100
        for p_left, p_right in model["states"]:
            q_left = (1 - c) * p_left + c * p_right
            q_right = c * p_left + (1 - c) * p_right
            q_left = ERROR + (1 - 2 * ERROR) * q_left
            q_right = ERROR + (1 - 2 * ERROR) * q_right
            overlay_ns["predictive_region"](ax, q_left, q_right)
            ax.plot(
                q_left, q_right,
                marker="o", linestyle="none",
                markersize=3.8,
                markerfacecolor="#F000B8",
                markeredgecolor="white", markeredgewidth=.65,
                zorder=5,
            )

        # Equal padding and limits preserve square axes and complete edge bins.
        ax.set_xlim(-.025, 1.025)
        ax.set_ylim(-.025, 1.025)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([0, .5, 1])
        ax.set_yticks([0, .5, 1])
        ax.set_xticklabels(["0", "0.5", "1"])
        ax.set_yticklabels(["0", "0.5", "1"])
        ax.tick_params(length=3, pad=2, labelsize=8)
        if column > 0:
            ax.tick_params(labelleft=False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

maximum = max(float(item.get_array().max()) for item in collections)
norm = LogNorm(vmin=1, vmax=maximum)
for item in collections:
    item.set_norm(norm)

cax = fig.add_axes([0.945, 0.31, 0.012, 0.31])
colourbar = fig.colorbar(collections[0], cax=cax)
colourbar.set_label("Simulated sites per hexagon", fontsize=8)
colourbar.ax.tick_params(labelsize=7, length=2)

for mixing, label in [
    (0, "Unmixed"),
    (20, "20% reciprocal read mixing"),
]:
    fig.text(
        0.023, ROW_Y[mixing] + HEIGHT / 2,
        label, rotation=90,
        va="center", ha="center", fontsize=10, weight="bold",
    )

fig.text(
    0.068, 0.385, "Right ALT-read fraction",
    rotation=90, ha="center", va="center", fontsize=10,
)
fig.text(
    0.505, 0.027, "Left ALT-read fraction",
    ha="center", va="center", fontsize=10,
)

prefix = OUT / "Fig_Sx_alternative_haplotype_compositions_50pct_horizontal_vector"


# ============================================================
# HORIZONTAL LAYOUT RETAINED
#
# The finished pre-transpose figure is saved directly:
#     columns = haplotype models A-E
#     rows    = Unmixed / 20% reciprocal read mixing
# ============================================================

fig.savefig(prefix.with_suffix(".pdf"))
fig.savefig(prefix.with_suffix(".svg"))

plt.close(fig)

legend = f"""Fig. Sx. Predicted allele-fraction patterns under alternative haplotype compositions.

(A) AB/AB: two diploid lineages with the same two founding haplotypes.
(B) AB/A: a diploid lineage and a haploid lineage sharing haplotype A.
(C) AB/AC: two diploid lineages sharing one founding haplotype.
(D) AA/BC: one diploid lineage containing two copies of A and another containing B and C.
(E) AB/CD: two diploid lineages with no shared founding haplotype.

Letters identify founding haplotypes; identical letters denote copies of the same haplotype.
The cartoons depict haplotype composition, not individual chromosomes or established developmental mechanisms.
The left and right AF columns show 0% and 20% reciprocal read mixing, respectively.
Each plot contains {n:,} simulated sites. Hexagon colours show site counts on a common logarithmic scale;
empty bins are white. Magenta dots indicate expected read-fraction centroids after the sequencing-error mapping.
Outlines show the corrected analytical 50% predictive regions at 40 reads per side and concentration 22.
These illustrative overlays use a fixed depth; the simulated backgrounds use empirical paired depths.
The outlines describe read-fraction variation, not confidence intervals for estimated centroids.

Simulations assume independent founding alleles with ALT probability 0.5, empirical paired TN
allele-assigned read depths, beta-binomial concentration 22, sequencing error 0.001 and no allelic bias.
The same underlying genotypes and depth pairs are retained between mixing scenarios within each model.
Genotypes for different model families were generated in separate simulation runs.
No genotype, posterior or centroid-distance selection was applied to the plotted simulations.

Mixing denotes the proportion of reads contributed by the opposite lineage, not the proportion of cells.
Empirical sequencing depths describe read-sampling variation; these plots do not predict absolute DNA
content or the coverage differences caused by ploidy. The doubled composition AABB/AA has the same
expected AF centroids as AB/A. These simulations illustrate predicted patterns rather than formal
model-selection results.
"""
(prefix.parent / (prefix.name + ".legend.txt")).write_text(legend)

(OUT / "figure_50pct.audit.json").write_text(json.dumps({
    "models": [model["name"] for model in MODELS],
    "simulation_sources": {
        key: str(path.resolve()) for key, path in FILES.items()
    },
    "sites_per_panel": n,
    "mixing": [0, .2],
    "hex_gridsize": BINS,
    "palette": str(PALETTE.resolve()),
    "shared_colour_limits": [1, maximum],
    "centroids": True,
    "predictive_contours": True,
    "predictive_mass": 0.5,
    "overlay_depth_per_side": 40,
    "overlay_source": str(corrected_path.resolve()),
    "predictive_region_source": overlay_source,
    "resimulation": False,
    "figure_inches": [FIG_W, FIG_H],
    "png_dpi": 300,
}, indent=2) + "\n")

print(f"PDF\t{prefix.with_suffix('.pdf').resolve()}")
print(f"PNG\t{prefix.with_suffix('.png').resolve()}")
print(f"LEGEND\t{(OUT / (prefix.name + '.legend.txt')).resolve()}")
