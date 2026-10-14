from pathlib import Path
import csv
import gzip
import json
import importlib.util
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from scipy.stats import betabinom

root = Path("analysis/figures/pairwise_AF_cli")
prefix = root / "TN_merged_control_sites_hex40_batlowW"
bins_path = Path(str(prefix) + ".AF_bins.tsv.gz")
audit_path = Path(str(prefix) + ".AF.audit.json")
output = root / "TN_AF_hex40_batlowW_overlay_0_20_50pct.pdf"

spec = importlib.util.spec_from_file_location(
    "af_tool", "pairwise-af-hex20/plot_pairwise_AF.py"
)
tool = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tool)

audit = json.loads(audit_path.read_text())
with gzip.open(bins_path, "rt") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

polygons, counts = [], []
for row in rows:
    x, y = float(row["left_centre"]), float(row["right_centre"])
    w, h = float(row["width"]), float(row["height"])
    count = int(row["count"])
    if not np.isfinite([x, y, w, h]).all() or min(w, h) <= 0 or count < 0:
        raise RuntimeError("Invalid hexagon record.")
    polygons.append([
        (x + w/2, y - h/4), (x + w/2, y + h/4),
        (x, y + h/2),
        (x - w/2, y + h/4), (x - w/2, y - h/4),
        (x, y - h/2),
    ])
    counts.append(count)

polygons = np.asarray(polygons)
counts = np.asarray(counts, dtype=np.int64)
if not len(counts) or not (counts > 0).any():
    raise RuntimeError("No observed counts.")
if int(counts.sum()) != audit["records_retained"]:
    raise RuntimeError("Bin total does not match the audit.")

params = audit["parameters"]
vmin = params.get("vmin")
vmax = params.get("vmax")
vmin = 1 if vmin is None else float(vmin)
vmax = max(float(counts.max()), vmin * 10) if vmax is None else float(vmax)

depth, concentration, error = 40, 22., .001
k = np.arange(depth + 1)
af = k / depth
scenarios = [(0., "#F000B8"), (.2, "#00D9EF")]
states = [(l, r) for l in (0., .5, 1.) for r in (0., .5, 1.)]
halo = [pe.Stroke(linewidth=3, foreground="white"), pe.Normal()]

def centroid(left, right, mixing):
    return (
        (1 - mixing) * left + mixing * right,
        mixing * left + (1 - mixing) * right,
    )

def probabilities(q):
    mean = error + (1 - 2 * error) * q
    values = betabinom.pmf(
        k, depth, mean * concentration, (1 - mean) * concentration
    )
    return values / values.sum()

def threshold(pmf, mass):
    ordered = np.sort(pmf.ravel())[::-1]
    index = np.searchsorted(np.cumsum(ordered), mass)
    return ordered[min(index, len(ordered) - 1)]

plt.rcParams.update({"font.size": 11, "pdf.fonttype": 42})
fig = plt.figure(figsize=(11, 10))
ax = fig.add_axes([.10, .20, .70, .70])
cax = fig.add_axes([.85, .32, .025, .46])

occupied = counts > 0
hexes = PolyCollection(
    polygons[occupied],
    array=counts[occupied].astype(float),
    cmap=tool.colourmap("batlowW_r"),
    norm=LogNorm(vmin=vmin, vmax=vmax),
    edgecolors="none", zorder=0,
)
ax.add_collection(hexes)
ax.set_facecolor("white")
fig.colorbar(hexes, cax=cax, label="Observed sites per hexagon")
cax.text(.5, -.09, "White = 0", transform=cax.transAxes,
         ha="center", fontsize=9)

for left, right in states:
    first = centroid(left, right, 0.)
    second = centroid(left, right, .2)
    if left != right:
        line, = ax.plot(
            [first[0], second[0]], [first[1], second[1]],
            color="#333333", lw=1.1, zorder=3,
        )
        line.set_path_effects(halo)

    for mixing, color in scenarios:
        # Identical genotype pairs have unchanged centroids and spread.
        if left == right and mixing != 0:
            continue
        if left == right:
            color = "#333333"

        ql, qr = centroid(left, right, mixing)
        pmf = np.outer(probabilities(qr), probabilities(ql))
        assert np.isclose(pmf.sum(), 1)

        for mass, style in [(0.50, "-")]:
            level = threshold(pmf, mass)
            if pmf.min() < level < pmf.max():
                contours = ax.contour(
                    af, af, pmf, levels=[level],
                    colors=[color], linestyles=[style],
                    linewidths=1.35, zorder=4,
                )
                # Support both older and newer Matplotlib.
                if hasattr(contours, "set_path_effects"):
                    contours.set_path_effects(halo)
                else:
                    for collection in contours.collections:
                        collection.set_path_effects(halo)

        ax.scatter(
            ql, qr, s=55, color=color, edgecolors="white",
            linewidths=1.1, zorder=6,
        )

vertices = polygons.reshape(-1, 2)
ax.set(
    xlim=(vertices[:, 0].min() - .004, vertices[:, 0].max() + .004),
    ylim=(vertices[:, 1].min() - .004, vertices[:, 1].max() + .004),
    aspect="equal",
    xlabel=f"{params['left_sample']}\nALT-read fraction",
    ylabel=f"{params['right_sample']}\nALT-read fraction",
)
ax.set_xticks(np.linspace(0, 1, 6))
ax.set_yticks(np.linspace(0, 1, 6))
ax.spines[["top", "right"]].set_visible(False)

fig.suptitle(
    "TN gynandromorph: observed AF and mixing predictions",
    fontsize=18, weight="bold", y=.965,
)
fig.text(.10, .92, f"{counts.sum():,} observed sites | 40 hexagonal columns")

fig.legend(
    handles=[
        Line2D([0], [0], marker="o", color="#F000B8", label="0% mixing"),
        Line2D([0], [0], marker="o", color="#00D9EF", label="20% mixing"),
        Line2D([0], [0], color="#333333", label="50% predictive region"),

    ],
    loc="lower center", bbox_to_anchor=(.5, .10),
    ncol=2, frameon=False,
)
fig.text(
    .10, .035,
    "Illustrative contours: depth 40/side; beta-binomial concentration 22; "
    "error 0.001; no bias.\n"
    "Symmetric mixing. Stationary states shown once in grey. "
    "Contours describe read-fraction spread, not fitted-centroid confidence.",
    fontsize=9, color="#555555", linespacing=1.5,
)

fig.savefig(output)
plt.close(fig)
print(output.resolve())
