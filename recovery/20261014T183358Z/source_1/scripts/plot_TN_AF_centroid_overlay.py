from pathlib import Path
import csv
import gzip
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib.colors import LogNorm, Normalize, LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
from scipy.stats import betabinom

root = Path("analysis/figures/pairwise_AF_cli")
bins = root / "TN_merged_control_sites_hex20_batlowW.AF_bins.tsv.gz"
output = root / "TN_AF_with_centroid_contours_0_10_20_30_40.pdf"

with gzip.open(bins, "rt") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

required = {"left_centre", "right_centre", "width", "height", "count"}
if not rows or not required.issubset(rows[0]):
    raise RuntimeError("Missing rows or required hexagon columns.")

# Reconstruct the original hexagons directly; no rebinning.
polygons = []
counts = []
for row in rows:
    x = float(row["left_centre"])
    y = float(row["right_centre"])
    w = float(row["width"])
    h = float(row["height"])
    count = int(row["count"])
    if not np.isfinite([x, y, w, h]).all() or w <= 0 or h <= 0:
        raise RuntimeError("Invalid hexagon geometry.")
    if count < 0:
        raise RuntimeError("Negative bin count.")
    polygons.append([
        (x + w/2, y - h/4),
        (x + w/2, y + h/4),
        (x,       y + h/2),
        (x - w/2, y + h/4),
        (x - w/2, y - h/4),
        (x,       y - h/2),
    ])
    counts.append(count)

polygons = np.asarray(polygons)
counts = np.asarray(counts, dtype=np.int64)
occupied = counts > 0
if not occupied.any():
    raise RuntimeError("No occupied bins.")

# Same illustrative contour assumptions as the cartoon.
depth = 40
concentration = 22.
error = .001
k = np.arange(depth + 1)
af = k / depth
mixing = np.array([0., .1, .2, .3, .4])
states = [(l, r) for l in (0., .5, 1.) for r in (0., .5, 1.)]
cmap = plt.get_cmap("turbo")
mix_norm = Normalize(0, 40)

def probabilities(q):
    mean = error + (1 - 2 * error) * q
    values = betabinom.pmf(
        k, depth,
        mean * concentration,
        (1 - mean) * concentration,
    )
    return values / values.sum()

def threshold(pmf, mass):
    ordered = np.sort(pmf.ravel())[::-1]
    index = np.searchsorted(np.cumsum(ordered), mass)
    return ordered[min(index, len(ordered) - 1)]

plt.rcParams.update({"font.size": 11, "pdf.fonttype": 42})
fig = plt.figure(figsize=(12, 10))
ax = fig.add_axes([.09, .19, .67, .70])

# Neutral grey keeps the coloured model contours distinguishable.
grey = LinearSegmentedColormap.from_list(
    "observed_grey", ["#F7F7F7", "#646464"]
)
hexes = PolyCollection(
    polygons[occupied],
    array=counts[occupied].astype(float),
    cmap=grey,
    norm=LogNorm(vmin=1, vmax=max(2, int(counts.max()))),
    edgecolors="none",
    zorder=0,
)
ax.add_collection(hexes)

for left, right in states:
    x = (1 - mixing) * left + mixing * right
    y = mixing * left + (1 - mixing) * right

    for c, ql, qr in zip(mixing, x, y):
        if left == right and c != 0:
            continue

        pmf = np.outer(probabilities(qr), probabilities(ql))
        assert np.isclose(pmf.sum(), 1)
        color = "#333333" if left == right else cmap(mix_norm(c * 100))

        for mass, style, opacity in [
            (.25, "-", .95),
            (.75, "--", .65),
        ]:
            level = threshold(pmf, mass)
            if pmf.min() < level < pmf.max():
                ax.contour(
                    af, af, pmf,
                    levels=[level], colors=[color],
                    linestyles=[style], linewidths=1.1,
                    alpha=opacity, zorder=3,
                )

    if left != right:
        points = np.column_stack([x, y])
        segments = np.stack([points[:-1], points[1:]], axis=1)
        lines = LineCollection(
            segments, cmap=cmap, norm=mix_norm,
            linewidths=1.5, zorder=4,
        )
        lines.set_array((mixing[:-1] + mixing[1:]) * 50)
        ax.add_collection(lines)
        ax.scatter(
            x, y, c=mixing * 100,
            cmap=cmap, norm=mix_norm,
            s=35, edgecolors="white", linewidths=.6, zorder=5,
        )
    else:
        ax.scatter(
            left, right, color="#333333", marker="s",
            s=35, edgecolors="white", linewidths=.6, zorder=5,
        )

# Include complete edge hexagons.
vertices = polygons.reshape(-1, 2)
ax.set(
    xlim=(vertices[:, 0].min() - .004, vertices[:, 0].max() + .004),
    ylim=(vertices[:, 1].min() - .004, vertices[:, 1].max() + .004),
    aspect="equal",
    xlabel="TN Left ALT-read fraction",
    ylabel="TN Right ALT-read fraction",
)
ax.set_xticks(np.linspace(0, 1, 6))
ax.set_yticks(np.linspace(0, 1, 6))
ax.spines[["top", "right"]].set_visible(False)

mix_axis = fig.add_axes([.81, .56, .022, .26])
mix_bar = fig.colorbar(
    ScalarMappable(norm=mix_norm, cmap=cmap), cax=mix_axis
)
mix_bar.set_label("Symmetric mixing (%)")
mix_bar.set_ticks([0, 10, 20, 30, 40])

count_axis = fig.add_axes([.81, .23, .022, .23])
count_bar = fig.colorbar(hexes, cax=count_axis)
count_bar.set_label("Observed sites per hexagon")

fig.suptitle(
    "TN gynandromorph: observed AF and model contours",
    fontsize=19, weight="bold", y=.96,
)
fig.text(
    .09, .915,
    f"Observed sites: {int(counts.sum()):,}  |  Original saved hexagons; no additional filtering",
    fontsize=11,
)
fig.legend(
    handles=[
        Line2D([0], [0], color="#333333", lw=1.2,
               label="25% predictive region"),
        Line2D([0], [0], color="#333333", lw=1.2,
               linestyle="--", label="75% predictive region"),
    ],
    loc="lower left", bbox_to_anchor=(.09, .10),
    frameon=False, ncol=2,
)
fig.text(
    .09, .04,
    "Contours are illustrative, not fitted: 40 reads/side; concentration 22; "
    "error 0.001; no allele bias.\n"
    "Regions are conditional on each genotype pair. They describe read-fraction spread, "
    "not uncertainty in centroid estimates.",
    fontsize=9, color="#555555", linespacing=1.5,
)

fig.savefig(output)
plt.close(fig)
print(f"OBSERVED_SITES\t{int(counts.sum()):,}")
print(f"NEW_OVERLAY_PDF\t{output.resolve()}")
