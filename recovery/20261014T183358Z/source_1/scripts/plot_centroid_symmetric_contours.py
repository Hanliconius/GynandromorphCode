from pathlib import Path
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from scipy.stats import betabinom

parser = argparse.ArgumentParser()
parser.add_argument("--depth", type=int, default=40)
parser.add_argument("--concentration", type=float, default=22)
parser.add_argument("--error", type=float, default=0.001)
parser.add_argument(
    "--output",
    default="analysis/figures/centroid_symmetric_contours.pdf",
)
args = parser.parse_args()

if args.depth < 2 or args.concentration <= 0:
    parser.error("Depth must be >=2; concentration must be positive.")
if not 0 < args.error < 0.5:
    parser.error("Error must lie between 0 and 0.5.")

states = [(l, r) for l in (0., .5, 1.) for r in (0., .5, 1.)]
mixing = np.arange(0, 51, 5) / 100
counts = np.arange(args.depth + 1)
af = counts / args.depth
cmap = plt.get_cmap("turbo")
norm = Normalize(0, 50)

def probabilities(q):
    mean = args.error + (1 - 2 * args.error) * q
    pmf = betabinom.pmf(
        counts, args.depth,
        mean * args.concentration,
        (1 - mean) * args.concentration,
    )
    return pmf / pmf.sum()

def mass_threshold(pmf, mass=.5):
    ordered = np.sort(pmf.ravel())[::-1]
    index = np.searchsorted(np.cumsum(ordered), mass)
    return ordered[min(index, len(ordered) - 1)]

plt.rcParams.update({"font.size": 12, "pdf.fonttype": 42})
fig, ax = plt.subplots(figsize=(10, 10))
fig.subplots_adjust(left=.11, right=.85, bottom=.19, top=.88)

# Avoid drawing identical distributions repeatedly for stationary states
# or for genotype pairs that converge at 50% mixing.
seen_distributions = set()

for left, right in states:
    x = (1 - mixing) * left + mixing * right
    y = mixing * left + (1 - mixing) * right

    for c, ql, qr in zip(mixing, x, y):
        key = (round(float(ql), 10), round(float(qr), 10))
        if key in seen_distributions:
            continue
        seen_distributions.add(key)

        pmf = np.outer(probabilities(qr), probabilities(ql))
        assert np.isclose(pmf.sum(), 1)
        # Inner 25% and outer 75% predictive regions.
        color = "#777777" if left == right else cmap(norm(c * 100))

        for mass, linestyle, opacity in [
            (.25, "-", .80),
            (.75, "--", .45),
        ]:
            threshold = mass_threshold(pmf, mass=mass)
            if pmf.min() < threshold < pmf.max():
                ax.contour(
                    af, af, pmf,
                    levels=[threshold],
                    colors=[color],
                    linestyles=[linestyle],
                    linewidths=.85,
                    alpha=opacity,
                    zorder=1,
                )

    if left != right:
        points = np.column_stack([x, y])
        segments = np.stack([points[:-1], points[1:]], axis=1)
        track = LineCollection(
            segments, cmap=cmap, norm=norm,
            linewidths=2, zorder=3,
        )
        track.set_array((mixing[:-1] + mixing[1:]) * 50)
        ax.add_collection(track)

        ax.scatter(
            x, y, c=mixing * 100, cmap=cmap, norm=norm,
            s=12, edgecolors="none", zorder=4,
        )
        ax.scatter(
            [left], [right], s=60, facecolors="white",
            edgecolors=cmap(norm(0)), linewidths=1.5, zorder=5,
        )
    else:
        ax.scatter(
            [left], [right], s=45, marker="s",
            color="#555555", zorder=5,
        )

ax.plot([0, 1], [0, 1], color="#DDDDDD", lw=.8, zorder=0)
ax.set(
    xlim=(-.025, 1.025), ylim=(-.025, 1.025),
    aspect="equal",
    xlabel="Left ALT-read fraction",
    ylabel="Right ALT-read fraction",
)
ax.spines[["top", "right"]].set_visible(False)

cax = fig.add_axes([.89, .27, .023, .51])
bar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax)
bar.set_label("Symmetric mixing on each side (%)")
bar.set_ticks([0, 10, 20, 30, 40, 50])

fig.suptitle(
    "Centroid movement and predicted spread",
    fontsize=20, weight="bold", y=.96,
)
fig.text(
    .11, .91,
    "Symmetric mixing: 0–50%, in five-percentage-point steps",
    fontsize=12,
)
fig.text(
    .11, .115,
    "Connected points: centroids   |   Solid: 25% region   |   Dashed: 75% region\n"
    "Open circles: unmixed origins   |   Grey squares: stationary centroids",
    fontsize=11, linespacing=1.6,
)
fig.text(
    .11, .045,
    f"Illustrative beta-binomial model: depth={args.depth}/side; "
    f"concentration={args.concentration:g}; error={args.error:g}; no bias.\n"
    "Contours describe read-fraction spread, not centroid uncertainty. "
    "Discrete boundaries are interpolated.\n"
    "Coincident distributions are drawn once; contours are conditional on each genotype pair.",
    fontsize=9, color="#555555", linespacing=1.5,
)

output = Path(args.output)
output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(output)
plt.close(fig)
print(output.resolve())
