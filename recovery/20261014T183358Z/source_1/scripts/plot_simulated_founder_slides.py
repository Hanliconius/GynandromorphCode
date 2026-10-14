from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import betabinom

OUT = Path(
    "analysis/figures/pairwise_AF_cli/"
    "AF_narrative_frames/simulated_founder_slides"
)
OUT.mkdir(parents=True, exist_ok=True)

SEED = 20260921
SITES = 200000
DEPTH = 40
CONCENTRATION = 22.0
ERROR = 0.001
BINS = 40

RECTS = [
    [0.065, 0.27, 0.265, 0.471111],
    [0.367, 0.27, 0.265, 0.471111],
    [0.669, 0.27, 0.265, 0.471111],
]
TITLES = [
    "Diploid / two-founder\nSame two haplotypes",
    "Three-founder\nOne shared haplotype",
    "Four-founder\nFour independent haplotypes",
]
MODEL_COLOUR = "#F000B8"
HALO = [pe.Stroke(linewidth=3.2, foreground="white"), pe.Normal()]

# Same founder draws underpin all models and mixing scenarios.
rng = np.random.default_rng(SEED)
haplotypes = rng.binomial(1, 0.5, size=(SITES, 4))
a, b, c, d = haplotypes.T
pairs = [
    ((a + b) / 2, (a + b) / 2),
    ((a + b) / 2, (a + c) / 2),
    ((a + b) / 2, (c + d) / 2),
]

def expected(p_left, p_right, mixing):
    biological_left = (1 - mixing) * p_left + mixing * p_right
    biological_right = mixing * p_left + (1 - mixing) * p_right
    # Same sequencing-error mapping for simulations and overlays.
    return (
        ERROR + (1 - 2 * ERROR) * biological_left,
        ERROR + (1 - 2 * ERROR) * biological_right,
    )

def draw_reads(q, generator):
    latent = generator.beta(
        q * CONCENTRATION, (1 - q) * CONCENTRATION
    )
    return generator.binomial(DEPTH, latent) / DEPTH

simulations = {}
for scenario, mixing in enumerate((0.0, 0.2)):
    panels = []
    for model_index, (p_left, p_right) in enumerate(pairs):
        generator = np.random.default_rng(
            SEED + 1000 + scenario * 10 + model_index
        )
        q_left, q_right = expected(p_left, p_right, mixing)
        panels.append((
            draw_reads(q_left, generator),
            draw_reads(q_right, generator),
        ))
    simulations[mixing] = panels

states = [
    sorted(set(zip(left, right)))
    for left, right in pairs
]
if [len(value) for value in states] != [3, 7, 9]:
    raise RuntimeError("Unexpected founder-state support.")

# Derive one shared colour normalization across all six datasets.
temporary = plt.figure()
axis = temporary.add_subplot(111)
maximum = 1
for panels in simulations.values():
    for left, right in panels:
        artist = axis.hexbin(
            left, right, gridsize=BINS,
            extent=(0, 1, 0, 1), mincnt=1,
        )
        maximum = max(maximum, int(artist.get_array().max()))
        artist.remove()
plt.close(temporary)
norm = LogNorm(vmin=1, vmax=maximum)

grid = np.arange(DEPTH + 1) / DEPTH
xx, yy = np.meshgrid(grid, grid)

def predictive_region(ax, q_left, q_right):
    left_pmf = betabinom.pmf(
        np.arange(DEPTH + 1), DEPTH,
        q_left * CONCENTRATION,
        (1 - q_left) * CONCENTRATION,
    )
    right_pmf = betabinom.pmf(
        np.arange(DEPTH + 1), DEPTH,
        q_right * CONCENTRATION,
        (1 - q_right) * CONCENTRATION,
    )
    probability = np.outer(right_pmf, left_pmf)
    ranked = np.sort(probability.ravel())[::-1]
    threshold = ranked[
        np.searchsorted(np.cumsum(ranked), 0.5 * probability.sum())
    ]

    # Discrete endpoint states can have >50% mass in one cell.
    # Their predictive region is represented as a small square.
    if threshold >= probability.max() * (1 - 1e-12):
        iy, ix = np.unravel_index(
            probability.argmax(), probability.shape
        )
        ax.plot(
            grid[ix], grid[iy], marker="s", markersize=7,
            markerfacecolor="none", markeredgecolor=MODEL_COLOUR,
            markeredgewidth=1.2, zorder=5, clip_on=False,
            path_effects=HALO,
        )
    else:
        contour = ax.contour(
            xx, yy, probability, levels=[threshold],
            colors=[MODEL_COLOUR], linewidths=1.3, zorder=5,
        )
        contour.set_path_effects(HALO)

def overlay(ax, model_index, mixing):
    for p_left, p_right in states[model_index]:
        q_left, q_right = expected(p_left, p_right, mixing)
        predictive_region(ax, q_left, q_right)

        if mixing > 0 and p_left != p_right:
            origin_left, origin_right = expected(p_left, p_right, 0)
            arrow = ax.annotate(
                "", xy=(q_left, q_right),
                xytext=(origin_left, origin_right),
                arrowprops={
                    "arrowstyle": "->",
                    "color": "#555555",
                    "lw": 1.0,
                    "shrinkA": 2,
                    "shrinkB": 4,
                },
                zorder=6,
            )
            arrow.arrow_patch.set_path_effects([
                pe.Stroke(linewidth=2.6, foreground="white"),
                pe.Normal(),
            ])

        ax.plot(
            q_left, q_right, "o",
            color=MODEL_COLOUR, markersize=4.5,
            markeredgecolor="white", markeredgewidth=0.7,
            zorder=7, clip_on=False,
        )

def make_frame(mixing, show_model):
    fig = plt.figure(figsize=(16, 9))
    title = (
        f"Simulated founder states — {int(mixing * 100)}% mixing"
        + (" with model predictions" if show_model else "")
    )
    fig.text(0.065, 0.935, title, fontsize=24, weight="bold")
    fig.text(
        0.065, 0.875,
        f"{SITES:,} simulated sites per panel; "
        "all sites shown; independent read sampling on Left and Right.",
        fontsize=12, color="#555555",
    )

    for index, rect in enumerate(RECTS):
        ax = fig.add_axes(rect)
        left, right = simulations[mixing][index]
        artist = ax.hexbin(
            left, right,
            gridsize=BINS, extent=(0, 1, 0, 1),
            mincnt=1, cmap="viridis", norm=norm,
            linewidths=0, rasterized=True,
        )
        # Padding displays boundary hexagons in full.
        ax.set_xlim(-0.018, 1.018)
        ax.set_ylim(-0.022, 1.022)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(TITLES[index], fontsize=13, pad=12)
        ax.set_xlabel("Left ALT-read fraction", fontsize=11)
        ax.set_ylabel("Right ALT-read fraction", fontsize=11)
        ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
        ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
        ax.tick_params(labelsize=10)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.plot(
            [0, 1], [0, 1],
            color="white", alpha=0.35, lw=0.7,
        )
        if show_model:
            overlay(ax, index, mixing)

    colour_axis = fig.add_axes([0.957, 0.30, 0.010, 0.39])
    bar = fig.colorbar(artist, cax=colour_axis)
    bar.set_label("Simulated sites / hexagon", fontsize=10)
    bar.ax.tick_params(labelsize=9)

    fig.text(
        0.065, 0.165,
        (
            "Dots: expected read fractions  |  "
            "Contours: 50% predictive regions  |  "
            "Arrows: displacement from unmixed states"
            if show_model else
            "White: empty bins  |  Identical binning and colour scale "
            "across all four slides"
        ),
        fontsize=11,
    )
    fig.text(
        0.065, 0.075,
        "Illustrative simulation: independent founder alleles at frequency 0.5; "
        "40 reads/side; beta-binomial concentration 22; error 0.001; no bias.\n"
        "Founder counts refer to haplotype origins. Predictive regions describe "
        "read variation, not confidence in fitted parameters.",
        fontsize=9, color="#555555",
    )
    return fig

frames = [
    ("01_unmixed_simulations", 0.0, False),
    ("02_unmixed_simulations_model", 0.0, True),
    ("03_mixed20_simulations", 0.2, False),
    ("04_mixed20_simulations_model", 0.2, True),
]

with PdfPages(OUT / "simulated_founder_slides.pdf") as pdf:
    for name, mixing, show_model in frames:
        fig = make_frame(mixing, show_model)
        fig.savefig(OUT / f"{name}.pdf")
        fig.savefig(OUT / f"{name}.png", dpi=160)
        pdf.savefig(fig)
        plt.close(fig)

np.savez_compressed(
    OUT / "simulated_AF.npz",
    **{
        f"mix{int(mixing * 100)}_model{index + 2}_{side}": values
        for mixing, panels in simulations.items()
        for index, pair in enumerate(panels)
        for side, values in zip(("left", "right"), pair)
    },
)
(OUT / "simulation.audit.json").write_text(json.dumps({
    "seed": SEED,
    "sites_per_panel": SITES,
    "depth_per_side": DEPTH,
    "concentration": CONCENTRATION,
    "error_rate": ERROR,
    "founder_allele_frequency": 0.5,
    "symmetric_mixing": [0.0, 0.2],
    "state_counts": [len(value) for value in states],
    "genotype_selection": "none",
    "full_sib_structure": False,
    "predictive_mass": 0.5,
    "colourmap": "viridis",
    "shared_count_colour_limits": [1, maximum],
    "figure_inches": [16, 9],
    "png_pixels": [2560, 1440],
    "interpretation": "illustrative simulations, not fitted biological data",
}, indent=2) + "\n")

print(f"COMPLETE\t{OUT.resolve()}")
