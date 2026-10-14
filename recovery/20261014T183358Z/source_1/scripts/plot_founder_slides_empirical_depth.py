from pathlib import Path
import ast
import json
import os

import numpy as np
import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LogNorm
from matplotlib.backends.backend_pdf import PdfPages

BASE = Path("analysis/figures/pairwise_AF_cli")
OUT = BASE / "AF_narrative_frames/simulated_founder_empirical_depth"
OUT.mkdir(parents=True, exist_ok=True)

SOURCE = Path("scripts/plot_simulated_founder_slides.py")
PALETTE = Path("palettes/barfbow/luka_bleute.json")
INPUT_AUDIT = BASE / "TN_merged_control_sites_hex40_batlowW.AF.audit.json"
DEPTH_CACHE = OUT / "paired_depth_sample.npz"
DEPTH_AUDIT = OUT / "paired_depth_sample.audit.json"

SEED = 20260922
N = 200000
CONTOUR_N = 200000

audit = json.loads(INPUT_AUDIT.read_text())
params = audit["parameters"]
vcf = Path(params["vcf"])
samples = [params["left_sample"], params["right_sample"]]

if (
    params["min_depth"] != 10
    or params.get("max_depth") is not None
    or params["filter"] != "all"
    or params.get("sites")
    or params.get("contig")
):
    raise ValueError("Unexpected source selection settings.")

stat = vcf.stat()
fingerprint = {
    "path": str(vcf.resolve()),
    "bytes": stat.st_size,
    "mtime_ns": stat.st_mtime_ns,
}
for key in ("bytes", "mtime_ns"):
    if fingerprint[key] != audit["input"][key]:
        raise ValueError(f"Source BCF changed: {key}")

signature = {
    "input": fingerprint,
    "samples": samples,
    "seed": SEED,
    "sample_size": N,
    "minimum_AD_depth_both_sides": 10,
    "selection": "uniform random eligible sites; paired depths preserved",
}

if DEPTH_CACHE.exists() and DEPTH_AUDIT.exists():
    saved = json.loads(DEPTH_AUDIT.read_text())
    if saved["signature"] != signature:
        raise ValueError("Depth cache does not match current inputs.")
    with np.load(DEPTH_CACHE, allow_pickle=False) as cache:
        depths = cache["depths"]
    print("REUSING_PAIRED_DEPTH_SAMPLE", flush=True)
else:
    # Random-priority sampling gives a uniform sample without replacement,
    # without storing all 31 million eligible sites.
    generator = np.random.default_rng(SEED)
    kept_depths = np.empty((0, 2), dtype=np.int64)
    kept_keys = np.empty(0)
    buffer = []
    read = 0
    eligible = 0

    def absorb():
        global kept_depths, kept_keys
        if not buffer:
            return
        values = np.asarray(buffer, dtype=np.int64)
        keys = generator.random(len(values))
        values = np.concatenate([kept_depths, values])
        keys = np.concatenate([kept_keys, keys])
        if len(keys) > N:
            selected = np.argpartition(keys, N - 1)[:N]
            values = values[selected]
            keys = keys[selected]
        kept_depths, kept_keys = values, keys
        buffer.clear()

    def depth(call):
        ad = call.get("AD")
        if ad is None or len(ad) != 2:
            return None
        if any(value is None or value < 0 for value in ad):
            return None
        total = int(ad[0]) + int(ad[1])
        return total if total >= 10 else None

    with pysam.VariantFile(str(vcf)) as handle:
        if not set(samples).issubset(handle.header.samples):
            raise ValueError("TN samples missing from BCF.")
        handle.subset_samples(samples)

        for record in handle:
            read += 1
            if read % 1000000 == 0:
                print(f"READ\t{read:,}\tELIGIBLE\t{eligible:,}", flush=True)
            if (
                record.ref not in ("A", "C", "G", "T")
                or not record.alts
                or len(record.alts) != 1
                or record.alts[0] not in ("A", "C", "G", "T")
            ):
                continue
            left = depth(record.samples[samples[0]])
            right = depth(record.samples[samples[1]])
            if left is None or right is None:
                continue
            eligible += 1
            buffer.append((left, right))
            if len(buffer) >= 100000:
                absorb()
        absorb()

    if read != audit["records_read"] or eligible != audit["records_retained"]:
        raise ValueError(
            f"Selection mismatch: read={read}, eligible={eligible}"
        )
    after = vcf.stat()
    if (
        after.st_size != fingerprint["bytes"]
        or after.st_mtime_ns != fingerprint["mtime_ns"]
    ):
        raise RuntimeError("BCF changed during extraction.")

    depths = kept_depths
    temporary = OUT / "paired_depth_sample.partial.npz"
    np.savez_compressed(temporary, depths=depths)
    os.replace(temporary, DEPTH_CACHE)
    saved = {
        "signature": signature,
        "records_read": read,
        "eligible_sites": eligible,
        "depth_quantiles": {
            side: np.quantile(depths[:, i], [0, .05, .5, .95, 1]).tolist()
            for i, side in enumerate(("left", "right"))
        },
    }
    temporary = OUT / "paired_depth_sample.audit.partial.json"
    temporary.write_text(json.dumps(saved, indent=2) + "\n")
    os.replace(temporary, DEPTH_AUDIT)

if depths.shape != (N, 2) or (depths < 10).any():
    raise ValueError("Invalid paired-depth cache.")

# Reuse the existing layout and model definitions without running
# the old fixed-depth simulation.
tree = ast.parse(SOURCE.read_text())
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
ns = {"__name__": "empirical_depth_slides"}
exec(
    compile(ast.Module(body=nodes, type_ignores=[]), str(SOURCE), "exec"),
    ns,
)

palette = json.loads(PALETTE.read_text())
if not isinstance(palette, list):
    raise ValueError("Expected a list of colours in luka_bleute.json.")
cmap = ListedColormap(palette, name="luka_bleute")

def sample_af(q, depth_values, rng):
    latent = rng.beta(
        np.asarray(q) * ns["CONCENTRATION"],
        (1 - np.asarray(q)) * ns["CONCENTRATION"],
        size=len(depth_values),
    )
    return rng.binomial(depth_values, latent) / depth_values

rng = np.random.default_rng(SEED + 1)
a, b, c, d = rng.binomial(1, .5, size=(N, 4)).T
pairs = [
    ((a + b) / 2, (a + b) / 2),
    ((a + b) / 2, (a + c) / 2),
    ((a + b) / 2, (c + d) / 2),
]
simulations = {}
for scenario, mixing in enumerate((0.0, 0.2)):
    panels = []
    for index, (left, right) in enumerate(pairs):
        rng = np.random.default_rng(SEED + 100 + scenario * 10 + index)
        ql, qr = ns["expected"](left, right, mixing)
        panels.append((
            sample_af(ql, depths[:, 0], rng),
            sample_af(qr, depths[:, 1], rng),
        ))
    simulations[mixing] = panels

states = [sorted(set(zip(left, right))) for left, right in pairs]
if [len(value) for value in states] != [3, 7, 9]:
    raise ValueError("Unexpected genotype-state support.")

# Shared colour scale across all six simulated datasets.
fig, ax = plt.subplots()
maximum = 1
for panels in simulations.values():
    for left, right in panels:
        collection = ax.hexbin(
            left, right, gridsize=ns["BINS"],
            extent=(0, 1, 0, 1), mincnt=1,
        )
        maximum = max(maximum, int(collection.get_array().max()))
        collection.remove()
plt.close(fig)

# Monte Carlo predictive regions averaged over the SAME paired-depth
# distribution. Separate draws from those used for displayed simulations.
# Equal-area 0.025 AF cells; no KDE, smoothing, or jitter.
region_cache = {}
region_edges = np.linspace(0, 1, 41)
region_centres = (region_edges[:-1] + region_edges[1:]) / 2
padded_centres = np.r_[
    region_centres[0] - .025,
    region_centres,
    region_centres[-1] + .025,
]

def predictive_region(ax, q_left, q_right):
    """
    Original analytical model overlay.

    The hexbin background continues to use empirical paired TN depths.
    The model overlay is calculated at the original illustrative depth
    of 40 reads/side, matching the earlier clean predictive regions.
    """
    from scipy.stats import betabinom

    depth = 40
    concentration = 22.0

    k = np.arange(depth + 1)
    grid = k / depth

    def pmf(q):
        q = float(np.clip(q, 1e-12, 1.0 - 1e-12))

        p = betabinom.pmf(
            k,
            depth,
            q * concentration,
            (1.0 - q) * concentration,
        )

        p = np.asarray(p, dtype=float)
        p /= p.sum()

        return p

    density = np.outer(
        pmf(q_right),
        pmf(q_left),
    )

    ranked = np.sort(density.ravel())[::-1]
    cumulative = np.cumsum(ranked)

    threshold_index = int(
        np.searchsorted(cumulative, 0.50)
    )
    threshold = ranked[
        min(threshold_index, len(ranked) - 1)
    ]

    if threshold >= density.max() * (1.0 - 1e-12):
        iy, ix = np.unravel_index(
            np.argmax(density),
            density.shape,
        )

        marker, = ax.plot(
            grid[ix],
            grid[iy],
            marker="s",
            markersize=5.5,
            markerfacecolor="none",
            markeredgecolor=ns["MODEL_COLOUR"],
            markeredgewidth=1.3,
            zorder=5,
        )
        marker.set_path_effects(ns["HALO"])
        return

    contour = ax.contour(
        grid,
        grid,
        density,
        levels=[threshold],
        colors=[ns["MODEL_COLOUR"]],
        linewidths=1.3,
        zorder=5,
    )

    contour.set_path_effects(ns["HALO"])



ns.update({
    "simulations": simulations,
    "states": states,
    "norm": LogNorm(vmin=1, vmax=maximum),
    "predictive_region": predictive_region,
    "SITES": N,
})

frames = [
    ("01_unmixed_simulations", 0.0, False),
    ("02_unmixed_simulations_model", 0.0, True),
    ("03_mixed20_simulations", 0.2, False),
    ("04_mixed20_simulations_model", 0.2, True),
]
with PdfPages(OUT / "simulated_founder_empirical_depth_slides.pdf") as pdf:
    for name, mixing, show_model in frames:
        fig = ns["make_frame"](mixing, show_model)

        for ax in fig.axes[:3]:
            for collection in ax.collections:
                if collection.get_array() is not None:
                    if collection.get_cmap().name == "viridis":
                        collection.set_cmap(cmap)

        for text in fig.texts:
            if text.get_text().startswith("Illustrative simulation:"):
                text.set_text(
                    "Illustrative founder alleles at frequency 0.5; "
                    "empirical paired TN read depths; concentration 22; "
                    "error 0.001; no bias.\n"
                    "Contours: analytical 50% beta-binomial predictive regions at "
                    "40 reads/side; not fitted-parameter confidence."
                )

        fig.canvas.draw()
        fig.savefig(OUT / f"{name}.pdf")
        fig.savefig(OUT / f"{name}.png", dpi=160)
        pdf.savefig(fig)
        plt.close(fig)

np.savez_compressed(
    OUT / "simulated_AF.npz",
    **{
        f"mix{int(mixing*100)}_model{index+2}_{side}": values
        for mixing, panels in simulations.items()
        for index, pair in enumerate(panels)
        for side, values in zip(("left", "right"), pair)
    },
)
(OUT / "simulation.audit.json").write_text(json.dumps({
    "seed": SEED,
    "sites_per_panel": N,
    "depth_source": str(DEPTH_CACHE.resolve()),
    "depth_sampling": signature["selection"],
    "same_depth_pairs_across_models": True,
    "founder_allele_frequency": .5,
    "concentration": ns["CONCENTRATION"],
    "error_rate": ns["ERROR"],
    "symmetric_mixing": [0, .2],
    "palette": "luka_bleute",
    "shared_colour_limits": [1, maximum],
    "predictive_draws_per_state": CONTOUR_N,
    "predictive_AF_cell_width": .025,
    "predictive_regions": [
        {
            "expected_left_AF": key[0],
            "expected_right_AF": key[1],
            "enclosed_probability": value[2],
        }
        for key, value in region_cache.items()
    ],
    "interpretation": "illustrative; not empirical founder-support calibration",
}, indent=2) + "\n")

print(f"COMPLETE\t{OUT.resolve()}", flush=True)
