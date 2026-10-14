from pathlib import Path
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import betabinom

parser = argparse.ArgumentParser()
parser.add_argument("--depth", type=int, default=40)
parser.add_argument("--concentration", type=float, default=22)
parser.add_argument("--error", type=float, default=0.001)
parser.add_argument("--output", default="centroid_contours_preview.pdf")
args = parser.parse_args()

if args.depth < 2 or args.concentration <= 0:
    parser.error("Depth must be >=2 and concentration >0.")
if not 0 < args.error < 0.5:
    parser.error("Error must be between 0 and 0.5.")

states = [(left, right)
          for left in (0., 0.5, 1.)
          for right in (0., 0.5, 1.)]
scenarios = [
    (0.00, 0.00, "#555555", "No mixing"),
    (0.20, 0.20, "#0072B2", "Left 20%, Right 20%"),
    (0.10, 0.30, "#D55E00", "Left 10%, Right 30%"),
]
fractions = np.arange(args.depth + 1) / args.depth
counts = np.arange(args.depth + 1)

def centroid(p_left, p_right, c_left, c_right):
    return (
        (1 - c_left) * p_left + c_left * p_right,
        c_right * p_left + (1 - c_right) * p_right,
    )

def read_probabilities(q):
    mean = args.error + (1 - 2 * args.error) * q
    result = betabinom.pmf(
        counts, args.depth,
        mean * args.concentration,
        (1 - mean) * args.concentration,
    )
    return result / result.sum()

def prediction(c_left, c_right):
    density = np.zeros((args.depth + 1, args.depth + 1))
    centres = []
    for p_left, p_right in states:
        q_left, q_right = centroid(
            p_left, p_right, c_left, c_right
        )
        centres.append((q_left, q_right))
        density += np.outer(
            read_probabilities(q_right),
            read_probabilities(q_left),
        ) / len(states)
    assert np.isclose(density.sum(), 1)
    return density, np.asarray(centres)

def threshold(density, mass):
    ordered = np.sort(density.ravel())[::-1]
    index = np.searchsorted(np.cumsum(ordered), mass)
    return ordered[min(index, len(ordered) - 1)]

plt.rcParams.update({"font.size": 11, "pdf.fonttype": 42})
fig, ax = plt.subplots(figsize=(10, 9))
fig.subplots_adjust(left=0.10, right=0.97, bottom=0.23, top=0.88)

# All independent mixing combinations, including asymmetric cases.
grid = np.linspace(0, 0.45, 10)
points = [
    centroid(pl, pr, cl, cr)
    for pl, pr in states
    for cl in grid
    for cr in grid
]
points = np.unique(np.asarray(points), axis=0)
ax.scatter(
    points[:, 0], points[:, 1],
    s=9, color="#ADB5BD", alpha=0.45, linewidths=0,
    zorder=1,
)

handles = []
for cl, cr, color, label in scenarios:
    density, centres = prediction(cl, cr)
    # Discrete highest-probability sets; lines interpolate their boundaries.
    for mass, style in [(0.90, "--"), (0.50, "-")]:
        level = threshold(density, mass)
        if density.min() < level < density.max():
            ax.contour(
                fractions, fractions, density,
                levels=[level], colors=[color],
                linestyles=[style], linewidths=1.3,
                alpha=0.85, zorder=2,
            )
    ax.scatter(
        centres[:, 0], centres[:, 1],
        s=45, facecolors="white", edgecolors=color,
        linewidths=1.5, zorder=4,
    )
    handles.append(Line2D(
        [0], [0], color=color, marker="o",
        markerfacecolor="white", label=label,
    ))

ax.plot([0, 1], [0, 1], color="#DDDDDD", lw=0.8, zorder=0)
ax.set(
    xlim=(-0.025, 1.025), ylim=(-0.025, 1.025),
    aspect="equal",
    xlabel="Left ALT-read fraction",
    ylabel="Right ALT-read fraction",
)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(
    handles=handles, loc="upper center",
    bbox_to_anchor=(0.5, -0.10), ncol=3, frameon=False,
)
fig.suptitle(
    "Candidate centroids and predicted read-fraction spread",
    fontsize=17, weight="bold", y=0.96,
)
fig.text(
    0.10, 0.91,
    "Grey points: full mixing grid. Coloured contours: separate mixing scenarios.",
    fontsize=11,
)
fig.text(
    0.10, 0.095,
    "Solid: 50% predictive mass   |   Dashed: 90% predictive mass\n"
    "Contours describe each nine-state mixture, not uncertainty in fitted centroids.",
    fontsize=10, linespacing=1.5,
)
fig.text(
    0.10, 0.035,
    f"ILLUSTRATIVE: depth={args.depth}/side; concentration={args.concentration:g}; "
    f"error={args.error:g}; no bias.\n"
    "Equal genotype-state weights; independent read sampling between sides. "
    "Discrete contour boundaries are interpolated.",
    fontsize=9, color="#555555", linespacing=1.5,
)

output = Path(args.output)
output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(output)
plt.close(fig)
print(output.resolve())
