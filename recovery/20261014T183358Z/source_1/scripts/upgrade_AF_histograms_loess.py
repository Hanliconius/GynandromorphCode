from pathlib import Path
import argparse
import importlib.util
import json
import shutil

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path("analysis/figures/pairwise_AF_cli")
OUT = ROOT / "AF_narrative_frames"
CACHE = OUT / "real_histogram_loess_inputs"
CACHE.mkdir(parents=True, exist_ok=True)

AUDITS = [
    ROOT / "control_benchmarks/best_coverage_control_self_hex45_batlowW.AF.audit.json",
    ROOT / "TN_merged_control_sites_hex40_batlowW.AF.audit.json",
]
LABELS = ["control", "TN"]

def extract(task):
    spec = importlib.util.spec_from_file_location(
        "af_tool", "pairwise-af-hex20/plot_pairwise_AF.py"
    )
    tool = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(tool)

    old = json.loads(AUDITS[task].read_text())
    params = old["parameters"]
    vcf = Path(old["input"]["path"])
    before = tool.fingerprint(vcf)
    if before != old["input"]:
        raise RuntimeError("Input VCF changed since the original plot.")

    # Fine internal grid for separating boundary peaks and constructing
    # 100-bin marginals. No new 2D AF plots are generated.
    argv = [
        "--vcf", str(vcf),
        "--output-prefix", str(CACHE / LABELS[task]),
        "--bins", "1000",
        "--bin-shape", "square",
    ]
    for key in (
        "left_sample", "right_sample", "min_depth",
        "max_depth", "filter",
    ):
        value = params.get(key)
        if value is not None:
            argv.extend(["--" + key.replace("_", "-"), str(value)])
    for contig in params.get("contig") or []:
        argv.extend(["--contig", contig])
    if params.get("allow_self"):
        argv.append("--allow-self")
    if params.get("sites"):
        white = old.get("whitelist")
        if not white or tool.sha256(Path(white["path"])) != white["sha256"]:
            raise RuntimeError("Missing or changed site whitelist.")
        argv.extend(["--sites", white["path"]])

    args = tool.parser().parse_args(argv)
    tool.check_args(args)
    histogram, audit = tool.accumulate(args, tool.load_sites(args.sites))

    if tool.fingerprint(vcf) != before:
        raise RuntimeError("VCF changed during extraction.")
    for key in (
        "records_read", "records_retained",
        "excluded_first_reason", "retained_by_contig",
    ):
        if audit[key] != old[key]:
            raise RuntimeError(f"Selection mismatch: {key}")

    # In the module's square histogram: first axis Left, second Right.
    left = histogram.sum(axis=1)
    right = histogram.sum(axis=0)
    if int(left.sum()) != old["records_retained"]:
        raise RuntimeError("Marginal counts do not match retained total.")

    destination = CACHE / f"{LABELS[task]}.npz"
    temporary = CACHE / f"{LABELS[task]}.partial.npz"
    np.savez_compressed(
        temporary,
        left=left,
        right=right,
        samples=np.array([params["left_sample"], params["right_sample"]]),
    )
    temporary.replace(destination)
    (CACHE / f"{LABELS[task]}.audit.json").write_text(
        json.dumps({
            "source_audit": str(AUDITS[task]),
            "source_input": before,
            "records_retained": audit["records_retained"],
            "selection_match": True,
            "internal_bins": 1000,
            "histogram_bins": 100,
        }, indent=2) + "\n"
    )
    print(f"EXTRACTION_PASS\t{LABELS[task]}\t{int(left.sum()):,}", flush=True)

def loess(x, y, fraction=.15):
    """Local-linear LOESS with tricube distance weights; no robust reweighting."""
    grid = np.linspace(x.min(), x.max(), 500)
    fitted = np.empty_like(grid)
    neighbours = min(len(x), max(7, int(np.ceil(fraction * len(x)))))
    for index, centre in enumerate(grid):
        distance = np.abs(x - centre)
        bandwidth = np.partition(distance, neighbours - 1)[neighbours - 1]
        bandwidth = max(bandwidth * 1.001, 1e-12)
        u = np.minimum(distance / bandwidth, 1)
        weight = (1 - u**3)**3
        design = np.column_stack([np.ones_like(x), x - centre])
        root = np.sqrt(weight)
        coefficients = np.linalg.lstsq(
            design * root[:, None], y * root, rcond=None
        )[0]
        fitted[index] = max(0., coefficients[0])
    return grid, fitted

def prepare(values):
    values = values.astype(float)
    # Sum fine counts into exact 1%-wide histogram bins.
    probability = values.reshape(100, 10).sum(axis=1) / values.sum()
    centres = (np.arange(100) + .5) / 100
    # Exclude endpoint bins from smoothing: preserve their high peaks.
    x, smooth = loess(centres[1:-1], probability[1:-1])
    return probability, x, smooth

def plot():
    with np.load(CACHE / "control.npz") as data:
        control = data["left"].copy()
        control_sample = str(data["samples"][0])
    with np.load(CACHE / "TN.npz") as data:
        left, right = data["left"].copy(), data["right"].copy()
        left_name, right_name = map(str, data["samples"])

    prepared = [prepare(values) for values in (control, left, right)]
    maximum = 1.15 * max(
        max(probability.max(), smooth.max())
        for probability, _, smooth in prepared
    )
    edges = np.linspace(0, 1, 101)
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 12,
        "pdf.fonttype": 42,
    })

    names = ["01_single_sample_histogram", "02_mixed_LR_histograms"]
    with PdfPages(OUT / "AF_real_histogram_frames_01_02_loess.pdf") as pdf:
        for frame, name in enumerate(names):
            fig = plt.figure(figsize=(16, 9))
            fig.text(
                .065, .93,
                "Observed AF distribution of a single control sample"
                if frame == 0 else "Observed Left and Right AF distributions",
                fontsize=24, weight="bold",
            )
            fig.text(
                .065, .865,
                "Fresh VCF read counts: 100-bin histograms with LOESS-smoothed interior distributions.",
                fontsize=14, color="#444444",
            )
            axis = fig.add_axes([.11, .27, .56, .49])

            series = (
                [(prepared[0], "#0072B2", control_sample)]
                if frame == 0 else
                [(prepared[1], "#0072B2", left_name),
                 (prepared[2], "#D55E00", right_name)]
            )
            for (probability, x, smooth), color, label in series:
                axis.stairs(
                    probability, edges, color=color,
                    alpha=.25, linewidth=1,
                )
                axis.plot(
                    x, smooth, color=color,
                    alpha=.8, linewidth=2.2, label=label,
                )
                axis.vlines(
                    [.005, .995], 0, probability[[0, -1]],
                    color=color, alpha=.7, linewidth=2,
                )

            axis.set_yscale("symlog", linthresh=.001, linscale=1)
            axis.set(
                xlim=(-.03, 1.03), ylim=(0, maximum),
                xlabel="ALT-read fraction",
                ylabel="Fraction of retained sites per 1%-wide bin",
            )
            axis.set_yticks([0, .001, .01, .1, 1])
            axis.yaxis.set_major_formatter(
                matplotlib.ticker.PercentFormatter(xmax=1)
            )
            axis.set_ylim(0, maximum)
            axis.spines[["top", "right"]].set_visible(False)
            axis.legend(frameon=False)

            count = int(control.sum() if frame == 0 else left.sum())
            note = (
                f"Control: {control_sample}" if frame == 0 else "TN gynandromorph"
            )
            note += (
                f"\n\n{count:,} retained sites\n\n"
                "Faint lines: observed histogram\n"
                "Strong curves: LOESS trend\n"
                "End bins: shown unsmoothed"
            )
            fig.text(.72, .69, note, va="top", fontsize=14, linespacing=1.7)
            fig.text(
                .065, .15,
                "Same axes in both frames: linear below 0.1% of sites, logarithmic above.",
                fontsize=12,
            )
            fig.text(
                .065, .055,
                "Original AF-module site filters retained. No genotype, posterior or anchor selection added.\n"
                "LOESS span: 15% of interior bins; smoothing is descriptive. "
                "The 0–1% and 99–100% bins are excluded from smoothing.",
                fontsize=10, color="#555555", linespacing=1.5,
            )

            for extension in ("pdf", "png"):
                destination = OUT / f"{name}.{extension}"
                backup = OUT / f"{name}.before_fresh_loess.{extension}"
                if destination.exists() and not backup.exists():
                    shutil.copy2(destination, backup)
                fig.savefig(destination, dpi=160)
            pdf.savefig(fig)
            plt.close(fig)
            print(f"UPGRADED\t{name}", flush=True)

    print(f"OUTPUT_DIRECTORY\t{OUT.resolve()}")
    print("Other frames and the original combined all-frames PDF are unchanged.")

parser = argparse.ArgumentParser()
parser.add_argument("mode", choices=["extract", "plot"])
parser.add_argument("--task", type=int, choices=[0, 1])
args = parser.parse_args()
if args.mode == "extract":
    if args.task is None:
        parser.error("extract requires --task")
    extract(args.task)
else:
    plot()
