from pathlib import Path
import ast
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import ListedColormap, LogNorm
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import betabinom

ROOT = Path("analysis/figures/pairwise_AF_cli/AF_narrative_frames")
OUT = ROOT / "all_founder_scenarios"
OUT.mkdir(parents=True, exist_ok=True)

original_script = Path("scripts/plot_simulated_founder_slides.py")
corrected_script = Path("scripts/plot_founder_slides_empirical_depth.py")
palette_path = Path("palettes/barfbow/luka_bleute.json")

inputs = {
    "standard": ROOT / "simulated_founder_empirical_depth/simulated_AF.npz",
    "AA_BC": ROOT / "gyno_three_founder/gyno_simulated_AF.npz",
    "AB_A": ROOT / "AB_A_digametic/AB_A_simulated_AF.npz",
}
for path in [
    original_script, corrected_script, palette_path, *inputs.values()
]:
    if not path.is_file():
        raise FileNotFoundError(f"Required file is not ready: {path}")

# Load original drawing definitions without running its simulation.
tree = ast.parse(original_script.read_text())
nodes = []
for node in tree.body:
    if isinstance(node, ast.Assign) and any(
        isinstance(target, ast.Name) and target.id == "rng"
        for target in node.targets
    ):
        break
    if isinstance(node, (ast.Import, ast.ImportFrom, ast.Assign)):
        nodes.append(node)
nodes.extend(
    node for node in tree.body if isinstance(node, ast.FunctionDef)
)
ns = {"__name__": "combined_scenario_plot"}
exec(
    compile(ast.Module(body=nodes, type_ignores=[]),
            str(original_script), "exec"),
    ns,
)

# Reuse the corrected contour implementation, including helper functions
# and literal constants, without executing its extraction/simulation.
corrected_source = corrected_script.read_text()
corrected_tree = ast.parse(corrected_source)
functions = [
    node for node in corrected_tree.body
    if isinstance(node, ast.FunctionDef)
]
predictive = [
    node for node in functions if node.name == "predictive_region"
]
if len(predictive) != 1:
    raise RuntimeError("Expected one corrected predictive_region function.")

predictive_text = ast.get_source_segment(
    corrected_source, predictive[0]
)
if "histogram2d" in predictive_text or "region_cache" in predictive_text:
    raise RuntimeError(
        "This script still contains the earlier binned contour method. "
        "Use the corrected analytical-overlay script before plotting."
    )

safe_nodes = [
    node for node in corrected_tree.body
    if isinstance(node, (ast.Import, ast.ImportFrom))
]
for node in corrected_tree.body:
    if isinstance(node, ast.Assign):
        try:
            ast.literal_eval(node.value)
        except (ValueError, TypeError, SyntaxError):
            continue
        safe_nodes.append(node)
safe_nodes.extend(functions)
exec(
    compile(ast.Module(body=safe_nodes, type_ignores=[]),
            str(corrected_script), "exec"),
    ns,
)

palette = json.loads(palette_path.read_text())
if not isinstance(palette, list):
    raise ValueError("Expected a palette colour list.")
cmap = ListedColormap(palette, name="luka_bleute")

names = ["AB/AB", "AB/A", "AB/AC", "AA/BC", "AB/CD"]
titles = [
    "AB / AB\nSame diploid",
    "AB / A\nDiploid / haploid",
    "AB / AC\nShared haplotype",
    "AA / BC\nDuplicated haplotype",
    "AB / CD\nFour haplotypes",
]
dosages = (0.0, 0.5, 1.0)
all_pairs = [(l, r) for l in dosages for r in dosages]
states = [
    [(p, p) for p in dosages],
    [(0., 0.), (.5, 0.), (.5, 1.), (1., 1.)],
    [pair for pair in all_pairs
     if pair not in [(0., 1.), (1., 0.)]],
    [(l, r) for l in (0., 1.) for r in dosages],
    all_pairs,
]

data = {}
with (
    np.load(inputs["standard"], allow_pickle=False) as standard,
    np.load(inputs["AA_BC"], allow_pickle=False) as duplicated,
    np.load(inputs["AB_A"], allow_pickle=False) as haploid,
):
    for mixing in (0, 20):
        def pair(cache, prefix):
            return (
                cache[prefix + "_left"].copy(),
                cache[prefix + "_right"].copy(),
            )
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
        if (
            not np.isfinite(left).all()
            or not np.isfinite(right).all()
            or ((left < 0) | (left > 1)).any()
            or ((right < 0) | (right > 1)).any()
        ):
            raise ValueError("Invalid simulated allele fractions.")
        sizes.add(len(left))
if len(sizes) != 1:
    raise ValueError(
        "Simulation sizes differ; cannot use directly comparable count colours."
    )
n = sizes.pop()

BINS = 40
temporary, ax = plt.subplots()
maximum = 1
for panels in data.values():
    for left, right in panels:
        artist = ax.hexbin(
            left, right, gridsize=BINS,
            extent=(0, 1, 0, 1), mincnt=1,
        )
        maximum = max(maximum, int(artist.get_array().max()))
        artist.remove()
plt.close(temporary)
norm = LogNorm(vmin=1, vmax=maximum)

# Connect the corrected overlay's nested namespace reference.
ns["MODEL_COLOUR"] = "#F000B8"
ns["ns"] = ns

# Test the corrected overlay before writing any final figures.
temporary, ax = plt.subplots()
try:
    ns["predictive_region"](ax, 0.5, 0.5)
finally:
    plt.close(temporary)

def make_frame(mixing, show_overlay):
    fig = plt.figure(figsize=(16, 9))
    fig.text(
        0.05, 0.93,
        f"Alternative haplotype compositions — {mixing}% read mixing",
        fontsize=25, weight="bold",
    )
    fig.text(
        0.05, 0.865,
        f"{n:,} simulated sites per model; empirical paired read depths; "
        "all simulated sites retained.",
        fontsize=12, color="#555555",
    )

    for index, (left, right) in enumerate(data[mixing]):
        ax = fig.add_axes([
            0.045 + index * 0.188, 0.31, 0.164, 0.30
        ])
        artist = ax.hexbin(
            left, right, gridsize=BINS,
            extent=(0, 1, 0, 1), mincnt=1,
            cmap=cmap, norm=norm,
            linewidths=0, rasterized=True,
        )
        ax.set_xlim(-0.018, 1.018)
        ax.set_ylim(-0.022, 1.022)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(titles[index], fontsize=12, pad=12)
        ax.set_xlabel("Left ALT-read fraction", fontsize=10)
        if index == 0:
            ax.set_ylabel("Right ALT-read fraction", fontsize=10)
        ax.set_xticks([0, .5, 1])
        ax.set_yticks([0, .5, 1])
        ax.tick_params(labelsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.plot([0, 1], [0, 1], color="white", alpha=.35, lw=.7)

        if show_overlay:
            for p_left, p_right in states[index]:
                q_left, q_right = ns["expected"](
                    p_left, p_right, mixing / 100
                )
                ns["predictive_region"](ax, q_left, q_right)
                ax.plot(
                    q_left, q_right, "o",
                    color="#F000B8", markersize=4,
                    markeredgecolor="white", markeredgewidth=.6,
                    zorder=8, clip_on=False,
                )

    cax = fig.add_axes([0.971, 0.32, 0.008, 0.27])
    colourbar = fig.colorbar(artist, cax=cax)
    colourbar.ax.tick_params(labelsize=8)
    colourbar.set_label("Sites / hexagon", fontsize=9)

    fig.text(
        0.05, 0.205,
        (
            "Dots: predicted centroids | Outlines: original analytical "
            "50% predictive overlays"
            if show_overlay else
            "Identical axes, binning and colour scale across all four slides."
        ),
        fontsize=12,
    )
    fig.text(
        0.05, 0.085,
        "Illustrative founder allele frequency 0.5; concentration 22; "
        "sequencing error 0.001; no bias. Mixing denotes read contribution.\n"
        "Backgrounds use empirical depths. Analytical overlays retain the "
        "illustrative 40 reads/side setting; they are not parameter confidence intervals.",
        fontsize=9, color="#555555",
    )
    return fig

frames = [
    ("01_all_models_0pct", 0, False),
    ("02_all_models_0pct_overlay", 0, True),
    ("03_all_models_20pct", 20, False),
    ("04_all_models_20pct_overlay", 20, True),
]
with PdfPages(OUT / "all_founder_scenarios.pdf") as pdf:
    for name, mixing, show_overlay in frames:
        fig = make_frame(mixing, show_overlay)
        fig.savefig(OUT / f"{name}.pdf")
        fig.savefig(OUT / f"{name}.png", dpi=160)
        pdf.savefig(fig)
        plt.close(fig)

(OUT / "plot.audit.json").write_text(json.dumps({
    "models": names,
    "sites_per_panel": n,
    "simulation_sources": {
        key: str(path.resolve()) for key, path in inputs.items()
    },
    "overlay_source": str(corrected_script.resolve()),
    "predictive_region_source": predictive_text,
    "palette": str(palette_path.resolve()),
    "colour_limits": [1, maximum],
    "mixing": [0, 0.2],
    "figure_inches": [16, 9],
    "png_pixels": [2560, 1440],
    "resimulation": False,
}, indent=2) + "\n")

print(f"COMPLETE\t{OUT.resolve()}")
