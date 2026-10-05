from pathlib import Path
import csv
import gzip
import json
import shutil

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import PercentFormatter

ROOT = Path("analysis/figures/pairwise_AF_cli")
OUT = ROOT / "AF_narrative_frames"
OUT.mkdir(parents=True, exist_ok=True)

CONTROL = ROOT / (
    "control_benchmarks/"
    "best_coverage_control_self_hex45_batlowW"
)
TN = ROOT / "TN_merged_control_sites_hex40_batlowW"

def load(prefix):
    audit = json.loads(
        Path(str(prefix) + ".AF.audit.json").read_text()
    )
    with gzip.open(
        Path(str(prefix) + ".AF_bins.tsv.gz"), "rt"
    ) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    left = np.array([float(r["left_centre"]) for r in rows])
    right = np.array([float(r["right_centre"]) for r in rows])
    counts = np.array([int(r["count"]) for r in rows], dtype=np.int64)

    if (
        not len(counts)
        or not np.isfinite(left).all()
        or not np.isfinite(right).all()
        or (counts < 0).any()
    ):
        raise RuntimeError(f"Invalid saved bins: {prefix}")
    if int(counts.sum()) != audit["records_retained"]:
        raise RuntimeError(f"Bin total differs from audit: {prefix}")

    return left, right, counts, audit

control_l, control_r, control_n, control_audit = load(CONTROL)
tn_l, tn_r, tn_n, tn_audit = load(TN)

control_params = control_audit["parameters"]
tn_params = tn_audit["parameters"]
if control_params["left_sample"] != control_params["right_sample"]:
    raise RuntimeError("Expected a control self-plot.")

# Common coarse bins avoid suggesting that cached hexagons retain
# exact per-site allele fractions. Boundaries lie between AF=0.05 steps.
edges = np.r_[0., np.arange(.025, 1., .05), 1.]

def marginal(centres, counts):
    histogram = np.histogram(
        np.clip(centres, 0, 1), bins=edges, weights=counts
    )[0]
    if int(histogram.sum()) != int(counts.sum()):
        raise RuntimeError("Marginal histogram lost counts.")
    return histogram / counts.sum()

# One marginal only for the self-plot; do not count the control twice.
control = marginal(control_l, control_n)
left = marginal(tn_l, tn_n)
right = marginal(tn_r, tn_n)

maximum = 1.15 * max(control.max(), left.max(), right.max())
BLUE, ORANGE = "#0072B2", "#D55E00"
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 12,
    "pdf.fonttype": 42,
})

def make_frame(mixed):
    fig = plt.figure(figsize=(16, 9))
    fig.text(
        .065, .93,
        "Observed Left and Right AF distributions" if mixed else
        "Observed AF distribution of a single control sample",
        fontsize=24, weight="bold",
    )
    fig.text(
        .065, .865,
        "Real read-count data; approximate marginals reconstructed from saved hexagons.",
        fontsize=14, color="#444444",
    )

    # Exactly the same histogram position as the original frames.
    axis = fig.add_axes([.11, .27, .56, .49])
    if mixed:
        axis.stairs(
            left, edges, color=BLUE, lw=2,
            alpha=.75, label=tn_params["left_sample"],
        )
        axis.stairs(
            right, edges, color=ORANGE, lw=2,
            alpha=.75, label=tn_params["right_sample"],
        )
        note = (
            "TN gynandromorph\n\n"
            f"{tn_n.sum():,} retained sites\n"
            "Both sides measured at the same sites.\n\n"
            "These are observed distributions;\n"
            "no mixing value is imposed here."
        )
    else:
        axis.stairs(
            control, edges, color=BLUE, lw=2,
            alpha=.75, label=control_params["left_sample"],
        )
        note = (
            f"Control: {control_params['left_sample']}\n\n"
            f"{control_n.sum():,} retained sites\n\n"
            "One observed sample.\n"
            "The self-plot cache supplies its\n"
            "single-sample AF distribution."
        )

    axis.set_yscale("symlog", linthresh=.001, linscale=1)
    axis.set(
        xlim=(-.03, 1.03), ylim=(0, maximum),
        xlabel="ALT-read fraction",
        ylabel="Fraction of retained sites per AF bin",
    )
    axis.set_yticks([0, .001, .01, .1, 1])
    axis.yaxis.set_major_formatter(PercentFormatter(xmax=1))
    axis.set_ylim(0, maximum)
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(frameon=False)
    fig.text(.72, .69, note, va="top", fontsize=14, linespacing=1.7)

    fig.text(
        .065, .15,
        "Same axes in both frames: linear below 0.1% of sites, logarithmic above.",
        fontsize=12,
    )
    fig.text(
        .065, .055,
        "Approximate AF bins centred at 0, 0.05, ..., 1; endpoint bins are half-width. "
        "All cached counts retained.\n"
        "Control and TN use their respective eligible-site sets. "
        "These panels show read fractions, not recalled genotypes.",
        fontsize=10, color="#555555", linespacing=1.5,
    )
    return fig

names = [
    "01_single_sample_histogram",
    "02_mixed_LR_histograms",
]

# Preserve the original illustrative frames once.
for name in names:
    for extension in ("pdf", "png"):
        path = OUT / f"{name}.{extension}"
        backup = OUT / f"{name}.illustrative_original.{extension}"
        if path.exists() and not backup.exists():
            shutil.copy2(path, backup)

with PdfPages(OUT / "AF_real_histogram_frames_01_02.pdf") as pdf:
    for mixed, name in zip((False, True), names):
        fig = make_frame(mixed)
        fig.savefig(OUT / f"{name}.pdf")
        fig.savefig(OUT / f"{name}.png", dpi=160)
        pdf.savefig(fig)
        plt.close(fig)
        print(f"REPLACED\t{name}.pdf / .png")

print(f"OUTPUT_DIRECTORY\t{OUT.resolve()}")
print("PNG_DIMENSIONS\t2560 x 1440")
print("The original combined all-frames PDF is unchanged.")
