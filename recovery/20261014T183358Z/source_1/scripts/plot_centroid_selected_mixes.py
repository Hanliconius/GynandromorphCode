from pathlib import Path
import ast
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import betabinom
from types import SimpleNamespace

# Reuse the probability calculations without executing the old plotter.
source = Path("scripts/plot_centroid_symmetric_contours.py")
tree = ast.parse(source.read_text())
wanted = {"probabilities", "mass_threshold"}
functions = [
    node for node in tree.body
    if isinstance(node, ast.FunctionDef) and node.name in wanted
]
if {node.name for node in functions} != wanted:
    raise RuntimeError("Required probability functions not found.")

# Same illustrative assumptions as the previous previews.
args = SimpleNamespace(depth=40, concentration=22, error=0.001)
counts = np.arange(args.depth + 1)
af = counts / args.depth
exec(compile(ast.Module(body=functions, type_ignores=[]),
             str(source), "exec"))

scenarios = [
    (0.00, 0.00, "#555555", "0 / 0%"),
    (0.20, 0.20, "#009E73", "20 / 20%"),
    (0.40, 0.40, "#CC79A7", "40 / 40%"),
    (0.10, 0.30, "#56B4E9", "10 / 30%"),
    (0.30, 0.10, "#0072B2", "30 / 10%"),
    (0.20, 0.40, "#E69F00", "20 / 40%"),
    (0.40, 0.20, "#D55E00", "40 / 20%"),
]
states = [(l, r) for l in (0., .5, 1.) for r in (0., .5, 1.)]

plt.rcParams.update({"font.size": 12, "pdf.fonttype": 42})
fig, ax = plt.subplots(figsize=(11, 9))
fig.subplots_adjust(left=.10, right=.76, bottom=.17, top=.88)

for cl, cr, color, label in scenarios:
    for pl, pr in states:
        # Identical-genotype states are stationary; draw them once.
        if pl == pr and (cl != 0 or cr != 0):
            continue

        ql = (1 - cl) * pl + cl * pr
        qr = cr * pl + (1 - cr) * pr
        pmf = np.outer(probabilities(qr), probabilities(ql))
        threshold = mass_threshold(pmf, mass=.50)

        if pmf.min() < threshold < pmf.max():
            ax.contour(
                af, af, pmf,
                levels=[threshold], colors=[color],
                linewidths=1.15, alpha=.75, zorder=2,
            )
        ax.scatter(
            [ql], [qr], s=35,
            facecolors=color, edgecolors="white",
            linewidths=.6, zorder=4,
        )

ax.plot([0, 1], [0, 1], color="#DDDDDD", lw=.8, zorder=0)
ax.set(
    xlim=(-.025, 1.025), ylim=(-.025, 1.025),
    aspect="equal",
    xlabel="Left ALT-read fraction",
    ylabel="Right ALT-read fraction",
)
ax.spines[["top", "right"]].set_visible(False)

handles = [
    Line2D([0], [0], color=color, marker="o",
           linewidth=1.2, label=label)
    for _, _, color, label in scenarios
]
ax.legend(
    handles=handles, title="Mixing: Left / Right",
    loc="upper left", bbox_to_anchor=(1.02, 1),
    frameon=False, labelspacing=1.2,
)

fig.suptitle(
    "Centroids and predicted spread across mixing scenarios",
    fontsize=18, weight="bold", y=.95,
)
fig.text(
    .10, .105,
    "Points: centroids   |   Contours: 50% predictive region per genotype pair",
    fontsize=11,
)
fig.text(
    .10, .045,
    "Illustrative: 40 reads/side; beta-binomial concentration 22; error 0.001; no bias.\n"
    "Stationary states shown once in grey. Contours are interpolated discrete boundaries,\n"
    "not confidence intervals on estimated centroids. Separate scenarios are not connected.",
    fontsize=9, color="#555555", linespacing=1.4,
)

output = Path("analysis/figures/centroid_selected_mixes_50pct.pdf")
output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(output)
plt.close(fig)
print(output.resolve())
