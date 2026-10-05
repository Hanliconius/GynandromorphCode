from pathlib import Path
import json
import os

import numpy as np
import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path("analysis/figures/pairwise_AF_cli")
OUT = ROOT / "AF_narrative_frames/original_GT_histograms"
OUT.mkdir(parents=True, exist_ok=True)

control_audit = json.loads((
    ROOT / "control_benchmarks/"
    "best_coverage_control_self_hex45_batlowW.AF.audit.json"
).read_text())
tn_audit = json.loads((
    ROOT / "TN_merged_control_sites_hex40_batlowW.AF.audit.json"
).read_text())

cp = control_audit["parameters"]
tp = tn_audit["parameters"]
vcf = Path(tp["vcf"]).resolve()

if Path(cp["vcf"]).resolve() != vcf:
    raise ValueError("Control and TN source files differ.")

for params in (cp, tp):
    if (
        params["min_depth"] != 10
        or params.get("max_depth") is not None
        or params["filter"] != "all"
        or params.get("sites")
        or params.get("contig")
    ):
        raise ValueError("Unexpected source filtering parameters.")

samples = [
    cp["left_sample"], tp["left_sample"], tp["right_sample"]
]
if len(set(samples)) != 3:
    raise ValueError("Expected three distinct samples.")

stat = vcf.stat()
fingerprint = {
    "path": str(vcf),
    "bytes": stat.st_size,
    "mtime_ns": stat.st_mtime_ns,
}
for audit in (control_audit, tn_audit):
    for key in ("bytes", "mtime_ns"):
        if audit["input"][key] != fingerprint[key]:
            raise ValueError(f"Source BCF changed: {key}")

expected = np.array([
    control_audit["records_retained"],
    tn_audit["records_retained"],
    tn_audit["records_retained"],
])
signature = {
    "input": fingerprint,
    "samples": samples,
    "bins": 50,
    "minimum_AD_depth": 10,
    "genotypes": "original VCF GT for all samples; no GQ filter",
    "selection": "control individually callable; TN both sides callable",
}
cache = OUT / "original_GT_counts.npz"
audit_path = OUT / "original_GT_counts.audit.json"

def get_counts(sample):
    ad = sample.get("AD")
    if ad is None or len(ad) != 2:
        return None
    if any(value is None or value < 0 for value in ad):
        return None
    ref, alt = map(int, ad)
    total = ref + alt
    return (alt, total) if total >= 10 else None

def category(sample):
    gt = sample.get("GT")
    if gt is None or len(gt) != 2:
        return 3
    if any(value not in (0, 1) for value in gt):
        return 3
    return sum(gt)

if cache.exists() and audit_path.exists():
    saved = json.loads(audit_path.read_text())
    if saved["signature"] != signature:
        raise ValueError("Cached inputs differ from current inputs.")
    with np.load(cache, allow_pickle=False) as data:
        hist = data["hist"]
    print("REUSING_ORIGINAL_GT_CACHE", flush=True)
else:
    hist = np.zeros((3, 4, 50), dtype=np.int64)
    read = 0

    with pysam.VariantFile(str(vcf)) as handle:
        missing = set(samples) - set(handle.header.samples)
        if missing:
            raise ValueError(f"Missing samples: {sorted(missing)}")
        handle.subset_samples(samples)

        for record in handle:
            read += 1
            if read % 1_000_000 == 0:
                print(f"READ\t{read:,}", flush=True)

            if (
                record.ref not in ("A", "C", "G", "T")
                or not record.alts
                or len(record.alts) != 1
                or record.alts[0] not in ("A", "C", "G", "T")
            ):
                continue

            calls = [record.samples[name] for name in samples]
            counts = [get_counts(call) for call in calls]
            keep = [
                counts[0] is not None,
                counts[1] is not None and counts[2] is not None,
                counts[1] is not None and counts[2] is not None,
            ]

            for index in range(3):
                if keep[index]:
                    alt, total = counts[index]
                    bin_id = min(49, alt * 50 // total)
                    hist[index, category(calls[index]), bin_id] += 1

    if read != tn_audit["records_read"]:
        raise ValueError("Record count differs from original audit.")
    if not np.array_equal(hist.sum(axis=(1, 2)), expected):
        raise ValueError(
            f"Retained totals differ: {hist.sum(axis=(1, 2))}; "
            f"expected {expected}"
        )

    after = vcf.stat()
    if (
        after.st_size != fingerprint["bytes"]
        or after.st_mtime_ns != fingerprint["mtime_ns"]
    ):
        raise RuntimeError("BCF changed during extraction.")

    partial = OUT / "original_GT_counts.partial.npz"
    np.savez_compressed(partial, hist=hist)
    os.replace(partial, cache)
    saved = {
        "signature": signature,
        "records_read": read,
        "retained": hist.sum(axis=(1, 2)).tolist(),
        "categories": ["0/0", "0/1", "1/1", "missing_or_non_diploid"],
        "category_counts": hist.sum(axis=2).tolist(),
    }
    partial_audit = OUT / "original_GT_counts.audit.partial.json"
    partial_audit.write_text(json.dumps(saved, indent=2) + "\n")
    os.replace(partial_audit, audit_path)

if hist.shape != (3, 4, 50):
    raise ValueError(f"Unexpected histogram shape: {hist.shape}")
if not np.array_equal(hist.sum(axis=(1, 2)), expected):
    raise ValueError("Cached counts do not match source audits.")

COLOURS = ["#2C9BCB", "#F28E2B", "#B61F72", "#AAAAAA"]
LABELS = ["0/0", "0/1", "1/1", "Missing/non-diploid GT"]
EDGES = np.linspace(0, 1, 51)
CENTRES = (EDGES[:-1] + EDGES[1:]) / 2
CAP = 0.015

# Exactly the same dimensions on both slides.
LEFT_RECT = [0.09, 0.27, 0.39, 0.43]
RIGHT_RECT = [0.57, 0.27, 0.39, 0.43]

def loess(x, y, span=0.20):
    """Local linear smoothing with tricube weights; no extra dependency."""
    neighbours = min(len(x), max(5, int(np.ceil(span * len(x)))))
    fitted = np.empty_like(y, dtype=float)

    for index, centre in enumerate(x):
        distance = np.abs(x - centre)
        bandwidth = np.partition(distance, neighbours - 1)[neighbours - 1]
        bandwidth = max(bandwidth * 1.000001, 1e-12)
        weights = np.maximum(0.0, 1.0 - (distance / bandwidth)**3)**3
        design = np.column_stack([np.ones(len(x)), x - centre])
        root_weights = np.sqrt(weights)
        coefficients = np.linalg.lstsq(
            design * root_weights[:, None],
            y * root_weights,
            rcond=None,
        )[0]
        fitted[index] = coefficients[0]

    return np.maximum(fitted, 0.0)


def draw_histogram(ax, counts, label):
    n = int(counts.sum())
    fraction = counts / n

    # Each genotype uses the same denominator: all eligible sites.
    # Smooth only the interior bins, leaving endpoint peaks separate.
    for index, colour in enumerate(COLOURS):
        fitted = loess(CENTRES[1:-1], fraction[index, 1:-1])
        ax.plot(
            CENTRES[1:-1], fitted,
            color=colour, alpha=0.85, linewidth=2.1,
        )

        for endpoint in (0, len(CENTRES) - 1):
            value = fraction[index, endpoint]
            if value <= 0:
                continue

            ax.vlines(
                CENTRES[endpoint], 0, min(value, CAP),
                color=colour, alpha=0.85, linewidth=2.1,
            )

        # Label actual genotype-bin percentages where they exceed the cap.
        for bin_id in np.flatnonzero(fraction[index] > CAP):
            value = fraction[index, bin_id]
            ax.plot(
                CENTRES[bin_id], CAP * 0.98,
                "^", color=colour, markersize=5,
            )
            ax.annotate(
                f"{100 * value:.2f}%",
                (CENTRES[bin_id], CAP),
                xytext=((index - 1.5) * 9, 5),
                textcoords="offset points",
                rotation=90, ha="center", va="bottom",
                fontsize=9, color=colour,
                annotation_clip=False,
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
        0.5, 0.92, f"{label}  |  n = {n:,}",
        transform=ax.transAxes, ha="center", fontsize=11,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.85),
    )
    ax.grid(axis="y", color="#E5E5E5", lw=0.7)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=11)

def frame(title):
    fig = plt.figure(figsize=(16, 9))
    fig.text(0.08, 0.93, title, fontsize=25, weight="bold")
    fig.text(
        0.08, 0.875,
        "Observed read fractions, coloured by original VCF genotype",
        fontsize=14, color="#555555",
    )
    fig.legend(
        handles=[
            Patch(facecolor=c, label=l)
            for c, l in zip(COLOURS, LABELS)
        ],
        loc="upper right", bbox_to_anchor=(0.96, 0.83),
        frameon=False, ncol=4, fontsize=11,
        title="Original VCF genotype",
    )
    fig.text(
        0.09, 0.12,
        "All depth-qualified sites retained, including missing GT. "
        "No reciprocal-mixture genotype assignment.\n"
        "LOESS curves; endpoint bins unsmoothed. Labels give actual clipped genotype-bin percentages.",
        fontsize=11, color="#555555",
    )
    return fig

fig1 = frame("Allele fractions in an unmixed control")
draw_histogram(fig1.add_axes(LEFT_RECT), hist[0], samples[0])

fig2 = frame("Allele fractions in the gynandromorph")
draw_histogram(fig2.add_axes(LEFT_RECT), hist[1], "Left")
draw_histogram(fig2.add_axes(RIGHT_RECT), hist[2], "Right")

combined = OUT / "control_and_gynandromorph_original_GT_LOESS_slides.pdf"
with PdfPages(combined) as pdf:
    for name, fig in [
        ("01_control_original_GT_LOESS", fig1),
        ("02_gynandromorph_original_GT_LOESS", fig2),
    ]:
        fig.savefig(OUT / f"{name}.pdf")
        fig.savefig(OUT / f"{name}.png", dpi=160)
        pdf.savefig(fig)
        plt.close(fig)

print(f"COMPLETE\t{combined.resolve()}", flush=True)
