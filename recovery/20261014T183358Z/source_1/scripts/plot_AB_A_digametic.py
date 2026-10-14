from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LogNorm

ROOT = Path("analysis/figures/pairwise_AF_cli/AF_narrative_frames")
DEPTH_FILE = (
    ROOT / "simulated_founder_empirical_depth/paired_depth_sample.npz"
)
PALETTE_FILE = Path("palettes/barfbow/luka_bleute.json")
OUT = ROOT / "AB_A_digametic"
OUT.mkdir(parents=True, exist_ok=True)

SEED = 20260923
CONCENTRATION = 22.0
ERROR = 0.001
BINS = 40

with np.load(DEPTH_FILE, allow_pickle=False) as cache:
    depths = cache["depths"]

if (
    depths.ndim != 2
    or depths.shape[1] != 2
    or len(depths) == 0
    or not np.isfinite(depths).all()
    or (depths < 10).any()
    or not np.equal(depths, np.floor(depths)).all()
):
    raise ValueError("Invalid empirical paired-depth cache.")
depths = depths.astype(np.int64)
n = len(depths)

palette = json.loads(PALETTE_FILE.read_text())
if not isinstance(palette, list):
    raise ValueError("Expected a list of palette colours.")
cmap = ListedColormap(palette, name="luka_bleute")

# Two independent founding haplotypes: Left AB, Right A.
rng = np.random.default_rng(SEED)
a, b = rng.binomial(1, 0.5, size=(n, 2)).T
p_left = (a + b) / 2.0
p_right = a.astype(float)

def expected(mixing):
    left = (1 - mixing) * p_left + mixing * p_right
    right = mixing * p_left + (1 - mixing) * p_right
    return (
        ERROR + (1 - 2 * ERROR) * left,
        ERROR + (1 - 2 * ERROR) * right,
    )

def reads(q, depth, generator):
    latent = generator.beta(
        q * CONCENTRATION,
        (1 - q) * CONCENTRATION,
    )
    return generator.binomial(depth, latent) / depth

simulations = {}
for index, mixing in enumerate((0.0, 0.2)):
    generator = np.random.default_rng(SEED + 10 + index)
    q_left, q_right = expected(mixing)
    simulations[mixing] = (
        reads(q_left, depths[:, 0], generator),
        reads(q_right, depths[:, 1], generator),
    )

# Same colour normalization for both scenarios.
temporary, ax = plt.subplots()
maximum = 1
for left, right in simulations.values():
    artist = ax.hexbin(
        left, right, gridsize=BINS,
        extent=(0, 1, 0, 1), mincnt=1,
    )
    maximum = max(maximum, int(artist.get_array().max()))
    artist.remove()
plt.close(temporary)
norm = LogNorm(vmin=1, vmax=maximum)

fig = plt.figure(figsize=(16, 9))
fig.text(
    0.065, 0.935,
    "Digametic scenario: Left AB (2n), Right A (n)",
    fontsize=24, weight="bold",
)
fig.text(
    0.065, 0.875,
    f"{n:,} simulated sites per panel; identical founding genotypes "
    "and empirical paired depths in both scenarios.",
    fontsize=12, color="#555555",
)

# Same panel dimensions as the previous three-column slides,
# with the pair centred on the slide.
rectangles = [
    [0.216, 0.27, 0.265, 0.471111],
    [0.518, 0.27, 0.265, 0.471111],
]

for rect, mixing in zip(rectangles, (0.0, 0.2)):
    ax = fig.add_axes(rect)
    left, right = simulations[mixing]
    artist = ax.hexbin(
        left, right, gridsize=BINS,
        extent=(0, 1, 0, 1),
        mincnt=1, cmap=cmap, norm=norm,
        linewidths=0, rasterized=True,
    )
    ax.set_xlim(-0.018, 1.018)
    ax.set_ylim(-0.022, 1.022)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(
        f"{int(mixing * 100)}% reciprocal read mixing\nLeft AB | Right A",
        fontsize=13, pad=12,
    )
    ax.set_xlabel("Left ALT-read fraction", fontsize=11)
    ax.set_ylabel("Right ALT-read fraction", fontsize=11)
    ax.set_xticks([0, .25, .5, .75, 1])
    ax.set_yticks([0, .25, .5, .75, 1])
    ax.tick_params(labelsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.plot([0, 1], [0, 1], color="white", alpha=.35, lw=.7)

colour_axis = fig.add_axes([0.817, 0.30, 0.010, 0.39])
bar = fig.colorbar(artist, cax=colour_axis)
bar.set_label("Simulated sites / hexagon", fontsize=10)
bar.ax.tick_params(labelsize=9)

fig.text(
    0.065, 0.165,
    "Left contains A and B; Right contains the shared haplotype A only.",
    fontsize=12,
)
fig.text(
    0.065, 0.075,
    "Illustrative founder allele frequency 0.5; empirical paired TN depths; "
    "beta-binomial concentration 22; error 0.001; no bias.\n"
    "All simulated sites retained. Shared colour scale between panels; "
    "white denotes empty bins. No model contours.",
    fontsize=9, color="#555555",
)

prefix = OUT / "AB_A_digametic_0_and_20pct"
fig.savefig(prefix.with_suffix(".pdf"))
fig.savefig(prefix.with_suffix(".png"), dpi=160)
plt.close(fig)

np.savez_compressed(
    OUT / "AB_A_simulated_AF.npz",
    true_left=p_left,
    true_right=p_right,
    **{
        f"mix{int(mixing*100)}_{side}": values
        for mixing, pair in simulations.items()
        for side, values in zip(("left", "right"), pair)
    },
)
(OUT / "AB_A_simulation.audit.json").write_text(json.dumps({
    "seed": SEED,
    "sites": n,
    "model": "Left=AB (2n); Right=A (n)",
    "founder_haplotypes": 2,
    "independent_founder_allele_frequency": 0.5,
    "depth_source": str(DEPTH_FILE.resolve()),
    "same_genotypes_and_depth_pairs_across_scenarios": True,
    "symmetric_mixing": [0.0, 0.2],
    "mixing_definition": "fraction of reads from the other lineage, not cells",
    "concentration": CONCENTRATION,
    "error_rate": ERROR,
    "palette": "luka_bleute",
    "hex_gridsize": BINS,
    "colour_limits": [1, maximum],
    "contours": False,
    "interpretation": "illustrative simulation, not a fit to observed data",
}, indent=2) + "\n")

print(f"PDF\t{prefix.with_suffix('.pdf').resolve()}")
print(f"PNG\t{prefix.with_suffix('.png').resolve()}")
