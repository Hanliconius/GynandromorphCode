from pathlib import Path
import importlib.util
import json
import os

import numpy as np
import pysam
from scipy.special import logsumexp

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter
from matplotlib.backends.backend_pdf import PdfPages

PROJECT = Path.cwd()
ROOT = PROJECT / "analysis/figures/pairwise_AF_cli"
OUT = ROOT / "AF_narrative_frames/genotype_histograms"
OUT.mkdir(parents=True, exist_ok=True)

CONTROL_AUDIT = (
    ROOT / "control_benchmarks/"
    "best_coverage_control_self_hex45_batlowW.AF.audit.json"
)
TN_AUDIT = ROOT / "TN_merged_control_sites_hex40_batlowW.AF.audit.json"
FIT = (
    PROJECT / "analysis/genomewide_mixing_corrected_v1/"
    "restored_pilot/observed.json"
)
MODEL = PROJECT / "scripts/159_fit_autosomal_three_vs_four_founder.py"

BIN_COUNT = 50
CUTOFF = 0.95
CAP = 0.015
CACHE = OUT / "genotype_histogram_counts.npz"
AUDIT = OUT / "genotype_histogram_counts.audit.json"

control_audit = json.loads(CONTROL_AUDIT.read_text())
tn_audit = json.loads(TN_AUDIT.read_text())
fit = json.loads(FIT.read_text())

cp = control_audit["parameters"]
tp = tn_audit["parameters"]
vcf = Path(tp["vcf"]).resolve()

if Path(cp["vcf"]).resolve() != vcf:
    raise ValueError("Control and TN audits refer to different VCFs.")

for params in (cp, tp):
    if (
        params["min_depth"] != 10
        or params.get("max_depth") is not None
        or params["filter"] != "all"
        or params.get("sites")
        or params.get("contig")
    ):
        raise ValueError("Input selection differs from the audited setup.")

samples = [
    cp["left_sample"],
    tp["left_sample"],
    tp["right_sample"],
]
if len(set(samples)) != 3:
    raise ValueError(f"Expected three distinct samples: {samples}")

stat = vcf.stat()
fingerprint = {
    "path": str(vcf),
    "bytes": stat.st_size,
    "mtime_ns": stat.st_mtime_ns,
}
for prior in (control_audit, tn_audit):
    for key in ("bytes", "mtime_ns"):
        if prior["input"][key] != fingerprint[key]:
            raise ValueError(f"BCF changed since source audit: {key}")

c_left = float(fit["alternative"]["c_left"])
c_right = float(fit["alternative"]["c_right"])
emission = fit["emission"]

signature = {
    "input": fingerprint,
    "samples": samples,
    "min_AD_depth": 10,
    "bins": BIN_COUNT,
    "c_left": c_left,
    "c_right": c_right,
    "emission": emission,
    "parameter_source": str(FIT),
    "parameter_fit_sites": fit["sites"],
    "parameter_scope": fit.get("scope"),
    "genotype_prior": "equal weights over all nine diploid pairs",
    "confidence": "maximum joint pair posterior >= 0.95",
    "site_selection": "no genotype or posterior-based exclusion",
    "control_colours": "original VCF GT; no GQ threshold",
}

def ad_counts(sample):
    ad = sample.get("AD")
    if ad is None or len(ad) != 2:
        return None
    if any(v is None or v < 0 for v in ad):
        return None
    ref, alt = int(ad[0]), int(ad[1])
    total = ref + alt
    return (alt, total) if total >= 10 else None

def original_category(sample):
    gt = sample.get("GT")
    if gt is None or len(gt) != 2:
        return 3
    if any(value not in (0, 1) for value in gt):
        return 3
    return sum(gt)

def bin_index(alt, total):
    # Exact integer arithmetic for 0.02-wide bins.
    return min(BIN_COUNT - 1, alt * BIN_COUNT // total)

if CACHE.exists() and AUDIT.exists():
    saved = json.loads(AUDIT.read_text())
    if saved["signature"] != signature:
        raise RuntimeError(
            "Cached extraction uses different inputs or parameters. "
            "Use a new output directory."
        )
    with np.load(CACHE, allow_pickle=False) as data:
        hist = data["hist"]
    if hist.shape != (3, 4, BIN_COUNT):
        raise ValueError("Unexpected cached histogram dimensions.")
    if not np.array_equal(
        hist.sum(axis=(1, 2)), np.asarray(saved["retained"])
    ):
        raise ValueError("Cache totals do not match its audit.")
    print("REUSING_GENOTYPE_HISTOGRAM_CACHE", flush=True)
else:
    spec = importlib.util.spec_from_file_location("existing_model", MODEL)
    model = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(model)

    # Derive categories from the model's own state order.
    gl = np.asarray(model.GL)
    gr = np.asarray(model.GR)
    expected_states = {
        (a, b) for a in (0.0, 0.5, 1.0)
        for b in (0.0, 0.5, 1.0)
    }
    if len(gl) != 9 or set(zip(gl, gr)) != expected_states:
        raise ValueError("Model does not contain the expected nine states.")

    hist = np.zeros((3, 4, BIN_COUNT), dtype=np.int64)
    batch = []
    read = 0

    def flush():
        if not batch:
            return
        array = np.asarray(batch, dtype=np.int64)
        counts = {
            "left_alt": array[:, 0],
            "left_total": array[:, 1],
            "right_alt": array[:, 2],
            "right_total": array[:, 3],
        }
        likelihood = model.state_log_likelihoods(
            counts, c_left, c_right, emission
        )
        if likelihood.shape != (len(array), 9):
            raise ValueError("Unexpected state likelihood dimensions.")
        if not np.isfinite(likelihood).all():
            raise ValueError("Nonfinite state likelihoods.")

        posterior = np.exp(
            likelihood - logsumexp(likelihood, axis=1)[:, None]
        )
        best = posterior.argmax(axis=1)
        confident = posterior[np.arange(len(array)), best] >= CUTOFF

        left_category = np.where(
            confident, np.rint(2 * gl[best]).astype(int), 3
        )
        right_category = np.where(
            confident, np.rint(2 * gr[best]).astype(int), 3
        )
        left_bins = np.minimum(
            BIN_COUNT - 1, array[:, 0] * BIN_COUNT // array[:, 1]
        )
        right_bins = np.minimum(
            BIN_COUNT - 1, array[:, 2] * BIN_COUNT // array[:, 3]
        )
        np.add.at(hist[1], (left_category, left_bins), 1)
        np.add.at(hist[2], (right_category, right_bins), 1)
        batch.clear()

    with pysam.VariantFile(str(vcf)) as handle:
        missing = set(samples) - set(handle.header.samples)
        if missing:
            raise ValueError(f"Samples absent from BCF: {sorted(missing)}")
        handle.subset_samples(samples)

        for record in handle:
            read += 1
            if read % 1_000_000 == 0:
                print(f"READ\t{read:,}", flush=True)

            if (
                len(record.ref) != 1
                or record.ref not in "ACGT"
                or not record.alts
                or len(record.alts) != 1
                or record.alts[0] not in ("A", "C", "G", "T")
            ):
                continue

            control = record.samples[samples[0]]
            control_counts = ad_counts(control)
            if control_counts is not None:
                alt, total = control_counts
                hist[0, original_category(control), bin_index(alt, total)] += 1

            left = ad_counts(record.samples[samples[1]])
            right = ad_counts(record.samples[samples[2]])
            if left is not None and right is not None:
                batch.append((*left, *right))
                if len(batch) >= 20000:
                    flush()
        flush()

    retained = hist.sum(axis=(1, 2))
    expected = np.array([
        control_audit["records_retained"],
        tn_audit["records_retained"],
        tn_audit["records_retained"],
    ])
    if not np.array_equal(retained, expected):
        raise ValueError(
            f"Site totals differ from original plots: "
            f"observed={retained.tolist()}, expected={expected.tolist()}"
        )
    if read != tn_audit["records_read"]:
        raise ValueError("BCF record count differs from original audit.")
    after = vcf.stat()
    if (
        after.st_size != fingerprint["bytes"]
        or after.st_mtime_ns != fingerprint["mtime_ns"]
    ):
        raise RuntimeError("BCF changed during extraction.")

    saved = {
        "signature": signature,
        "records_read": read,
        "retained": retained.tolist(),
        "category_order": ["0/0", "0/1", "1/1", "unresolved"],
        "category_counts": hist.sum(axis=2).tolist(),
    }
    partial = CACHE.with_name(CACHE.stem + ".partial.npz")
    np.savez_compressed(partial, hist=hist)
    os.replace(partial, CACHE)
    partial_audit = AUDIT.with_suffix(".partial.json")
    partial_audit.write_text(json.dumps(saved, indent=2) + "\n")
    os.replace(partial_audit, AUDIT)

colours = ["#2C9BCB", "#F28E2B", "#B61F72", "#AAAAAA"]
labels = ["0/0", "0/1", "1/1", "Unresolved"]
edges = np.linspace(0, 1, BIN_COUNT + 1)
centres = (edges[:-1] + edges[1:]) / 2

def histogram(ax, counts, title):
    n = int(counts.sum())
    fractions = counts / n
    bottom = np.zeros(BIN_COUNT)
    for category in range(4):
        ax.bar(
            edges[:-1], fractions[category],
            bottom=bottom, width=0.02, align="edge",
            color=colours[category], alpha=0.9,
            edgecolor="white", linewidth=0.35,
        )
        bottom += fractions[category]

    for index in np.flatnonzero(bottom > CAP):
        ax.plot(
            centres[index], CAP * 0.98,
            marker="^", color="#333333", markersize=5,
        )
        ax.annotate(
            f"{100 * bottom[index]:.2f}%",
            xy=(centres[index], CAP),
            xytext=(0, 5),
            textcoords="offset points",
            rotation=90, ha="center", va="bottom",
            fontsize=9, annotation_clip=False,
        )

    ax.axvline(0.5, color="#555555", ls="--", lw=0.9)
    ax.set_xlim(-0.015, 1.015)
    ax.set_ylim(0, CAP)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticks([0, 0.005, 0.010, 0.015])
    ax.yaxis.set_major_formatter(PercentFormatter(1, decimals=1))
    ax.set_xlabel("Observed ALT-read fraction", fontsize=13)
    ax.set_ylabel("Sites per 0.02 AF bin", fontsize=12)
    ax.text(
        0.5, 0.92, f"{title}  |  n = {n:,}",
        transform=ax.transAxes, ha="center", fontsize=12,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.85),
    )
    ax.grid(axis="y", color="#E5E5E5", lw=0.7)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=11)

def frame(title, subtitle):
    fig = plt.figure(figsize=(16, 9))
    fig.text(0.08, 0.93, title, fontsize=25, weight="bold")
    fig.text(0.08, 0.875, subtitle, fontsize=13, color="#555555")
    fig.legend(
        handles=[
            Patch(facecolor=c, label=l)
            for c, l in zip(colours, labels)
        ],
        loc="upper right", bbox_to_anchor=(0.94, 0.83),
        frameon=False, ncol=4, fontsize=12,
        title="Assigned genotype",
    )
    return fig

figures = []
fig = frame(
    "Allele fractions in an unmixed control",
    "Observed read fractions, coloured by the original VCF genotype",
)
histogram(
    fig.add_axes([0.11, 0.27, 0.78, 0.43]),
    hist[0], samples[0],
)
fig.text(
    0.11, 0.14,
    "All eligible control-site SNPs with AD depth ≥10; no genotype-based exclusion.\n"
    "Bars exceeding 1.5% are clipped and labelled with their full percentage.",
    fontsize=12, color="#555555",
)
figures.append(("01_control_genotype_histogram", fig))

fig = frame(
    "Allele fractions in the gynandromorph",
    "Observed read fractions, coloured by reciprocal-mixture genotype assignment",
)
histogram(
    fig.add_axes([0.09, 0.27, 0.39, 0.43]),
    hist[1], "Left",
)
histogram(
    fig.add_axes([0.57, 0.27, 0.39, 0.43]),
    hist[2], "Right",
)
fig.text(
    0.09, 0.14,
    f"Fixed mixing: Left {100*c_left:.3f}%; Right {100*c_right:.3f}%. "
    "Nine equally weighted genotype-pair states.\n"
    "Pair posterior <0.95 remains grey. All sites with AD depth ≥10 in both sides are retained.\n"
    "Parameters from the saved 96,756-site bridge fit; applied here genome-wide.",
    fontsize=11, color="#555555",
)
figures.append(("02_gynandromorph_genotype_histograms", fig))

with PdfPages(OUT / "control_and_gynandromorph_genotype_slides.pdf") as pdf:
    for name, fig in figures:
        fig.savefig(OUT / f"{name}.pdf")
        fig.savefig(OUT / f"{name}.png", dpi=160)
        pdf.savefig(fig)
        plt.close(fig)

print("COMPLETE", flush=True)
print(f"OUTPUT_DIRECTORY\t{OUT}", flush=True)
print(
    f"COMBINED_PDF\t{OUT / 'control_and_gynandromorph_genotype_slides.pdf'}",
    flush=True,
)
