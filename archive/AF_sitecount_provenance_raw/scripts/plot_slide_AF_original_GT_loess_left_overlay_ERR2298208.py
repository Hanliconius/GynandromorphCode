#!/usr/bin/env python3

from pathlib import Path
import json
import os

import numpy as np
import pysam

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import PercentFormatter


# ============================================================
# PATHS / INPUTS
# ============================================================

ROOT = Path(
    "analysis/figures/pairwise_AF_cli"
)

OUT = (
    ROOT
    / "AF_narrative_frames"
    / "original_GT_histograms"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)

tn_audit = json.loads(
    (
        ROOT
        / "TN_merged_control_sites_hex40_batlowW.AF.audit.json"
    ).read_text()
)

tp = tn_audit["parameters"]

vcf = Path(
    tp["vcf"]
).resolve()


# ============================================================
# SAMPLES
# ============================================================

CONTROL = "ERR2298208"
TN_LEFT = tp["left_sample"]
TN_RIGHT = tp["right_sample"]

samples = [
    CONTROL,
    TN_LEFT,
    TN_RIGHT,
]

if len(set(samples)) != 3:
    raise RuntimeError(
        f"Expected three distinct samples; observed {samples}"
    )

print(
    "SAMPLES",
    f"control={CONTROL}",
    f"left={TN_LEFT}",
    f"right={TN_RIGHT}",
    sep="\t",
    flush=True,
)


# ============================================================
# VERIFY SOURCE BCF
# ============================================================

if (
    tp["min_depth"] != 10
    or tp.get("max_depth") is not None
    or tp["filter"] != "all"
    or tp.get("sites")
    or tp.get("contig")
):
    raise RuntimeError(
        "Unexpected TN source filtering parameters."
    )

stat = vcf.stat()

fingerprint = {
    "path": str(vcf),
    "bytes": stat.st_size,
    "mtime_ns": stat.st_mtime_ns,
}

for key in (
    "bytes",
    "mtime_ns",
):
    if (
        tn_audit["input"][key]
        != fingerprint[key]
    ):
        raise RuntimeError(
            f"Source BCF changed: {key}"
        )


# ============================================================
# CACHE
#
# Separate from every previous control extraction.
# ============================================================

cache = (
    OUT
    / "ERR2298208_vs_TNLeft_original_GT_counts.npz"
)

audit_path = (
    OUT
    / "ERR2298208_vs_TNLeft_original_GT_counts.audit.json"
)

signature = {
    "input": fingerprint,
    "control": CONTROL,
    "TN_left": TN_LEFT,
    "TN_right_for_joint_callability": TN_RIGHT,
    "bins": 50,
    "minimum_AD_depth": 10,
    "genotypes": (
        "original VCF GT; no GQ filter"
    ),
    "selection": {
        "control": (
            "ERR2298208 individually AD-callable "
            "at depth >=10"
        ),
        "TN_left": (
            "TN Left plotted only where both TN Left "
            "and TN Right are AD-callable at depth >=10"
        ),
    },
}


# ============================================================
# ORIGINAL EXTRACTION HELPERS
# ============================================================

def get_counts(sample):
    ad = sample.get("AD")

    if (
        ad is None
        or len(ad) != 2
    ):
        return None

    if any(
        value is None
        or value < 0
        for value in ad
    ):
        return None

    ref, alt = map(
        int,
        ad,
    )

    total = ref + alt

    if total < 10:
        return None

    return alt, total


def category(sample):
    gt = sample.get("GT")

    if (
        gt is None
        or len(gt) != 2
    ):
        return 3

    if any(
        value not in (0, 1)
        for value in gt
    ):
        return 3

    return sum(gt)


# ============================================================
# EXTRACT / REUSE HISTOGRAMS
#
# hist dimensions:
#   sample x genotype x AF-bin
#
# sample 0 = ERR2298208
# sample 1 = TN Left
#
# We still read TN Right because the original TN figure retained
# a site for Left only when BOTH Left and Right were callable.
# ============================================================

if (
    cache.exists()
    and audit_path.exists()
):
    saved = json.loads(
        audit_path.read_text()
    )

    if (
        saved.get("signature")
        != signature
    ):
        raise RuntimeError(
            "Cached inputs differ from current inputs."
        )

    with np.load(
        cache,
        allow_pickle=False,
    ) as data:
        hist = data["hist"]

    print(
        "REUSING_ERR2298208_TNLEFT_CACHE",
        flush=True,
    )

else:
    hist = np.zeros(
        (2, 4, 50),
        dtype=np.int64,
    )

    read = 0

    with pysam.VariantFile(
        str(vcf)
    ) as handle:

        missing = (
            set(samples)
            - set(handle.header.samples)
        )

        if missing:
            raise RuntimeError(
                "Samples absent from source BCF: "
                + ", ".join(
                    sorted(missing)
                )
            )

        handle.subset_samples(
            samples
        )

        for record in handle:
            read += 1

            if (
                read % 1_000_000
                == 0
            ):
                print(
                    f"READ\t{read:,}",
                    flush=True,
                )

            if (
                record.ref
                not in ("A", "C", "G", "T")
                or not record.alts
                or len(record.alts) != 1
                or record.alts[0]
                not in ("A", "C", "G", "T")
            ):
                continue

            calls = [
                record.samples[name]
                for name in samples
            ]

            counts = [
                get_counts(call)
                for call in calls
            ]

            # ------------------------------------------------
            # ERR2298208:
            # individually callable, matching the original
            # control-side criterion.
            # ------------------------------------------------

            if counts[0] is not None:
                alt, total = counts[0]

                bin_id = min(
                    49,
                    alt * 50 // total,
                )

                hist[
                    0,
                    category(calls[0]),
                    bin_id,
                ] += 1

            # ------------------------------------------------
            # TN Left:
            # preserve original requirement that BOTH TN sides
            # are callable at the site.
            # ------------------------------------------------

            if (
                counts[1] is not None
                and counts[2] is not None
            ):
                alt, total = counts[1]

                bin_id = min(
                    49,
                    alt * 50 // total,
                )

                hist[
                    1,
                    category(calls[1]),
                    bin_id,
                ] += 1

    if (
        read
        != tn_audit["records_read"]
    ):
        raise RuntimeError(
            "Record count differs from original TN audit: "
            f"{read:,} versus "
            f"{tn_audit['records_read']:,}"
        )

    # Re-check that the BCF did not change during extraction.
    after = vcf.stat()

    if (
        after.st_size
        != fingerprint["bytes"]
        or after.st_mtime_ns
        != fingerprint["mtime_ns"]
    ):
        raise RuntimeError(
            "Source BCF changed during extraction."
        )

    partial = (
        OUT
        / "ERR2298208_vs_TNLeft_original_GT_counts.partial.npz"
    )

    np.savez_compressed(
        partial,
        hist=hist,
    )

    os.replace(
        partial,
        cache,
    )

    saved = {
        "signature": signature,
        "records_read": read,
        "retained": {
            "ERR2298208": int(
                hist[0].sum()
            ),
            "TN_Left": int(
                hist[1].sum()
            ),
        },
        "categories": [
            "0/0",
            "0/1",
            "1/1",
            "missing_or_non_diploid",
        ],
        "category_counts": {
            "ERR2298208": (
                hist[0]
                .sum(axis=1)
                .tolist()
            ),
            "TN_Left": (
                hist[1]
                .sum(axis=1)
                .tolist()
            ),
        },
    }

    partial_audit = (
        OUT
        / "ERR2298208_vs_TNLeft_original_GT_counts.audit.partial.json"
    )

    partial_audit.write_text(
        json.dumps(
            saved,
            indent=2,
        )
        + "\n"
    )

    os.replace(
        partial_audit,
        audit_path,
    )


if hist.shape != (2, 4, 50):
    raise RuntimeError(
        f"Unexpected histogram shape: {hist.shape}"
    )


# ============================================================
# ORIGINAL PLOT CONSTANTS
# ============================================================

COLOURS = [
    "#2C9BCB",
    "#F28E2B",
    "#B61F72",
    "#AAAAAA",
]

LABELS = [
    "0/0",
    "0/1",
    "1/1",
    "Missing/non-diploid GT",
]

EDGES = np.linspace(
    0,
    1,
    51,
)

CENTRES = (
    EDGES[:-1]
    + EDGES[1:]
) / 2

CAP = 0.015


# ============================================================
# EXACT ORIGINAL LOESS IMPLEMENTATION
# ============================================================

def loess(
    x,
    y,
    span=0.20,
):
    """
    Local linear smoothing with tricube weights;
    same implementation as the original figure.
    """

    neighbours = min(
        len(x),
        max(
            5,
            int(
                np.ceil(
                    span * len(x)
                )
            ),
        ),
    )

    fitted = np.empty_like(
        y,
        dtype=float,
    )

    for index, centre in enumerate(x):
        distance = np.abs(
            x - centre
        )

        bandwidth = np.partition(
            distance,
            neighbours - 1,
        )[neighbours - 1]

        bandwidth = max(
            bandwidth * 1.000001,
            1e-12,
        )

        weights = (
            np.maximum(
                0.0,
                1.0
                - (
                    distance
                    / bandwidth
                )**3,
            )**3
        )

        design = np.column_stack(
            [
                np.ones(len(x)),
                x - centre,
            ]
        )

        root_weights = np.sqrt(
            weights
        )

        coefficients = np.linalg.lstsq(
            (
                design
                * root_weights[:, None]
            ),
            y * root_weights,
            rcond=None,
        )[0]

        fitted[index] = (
            coefficients[0]
        )

    return np.maximum(
        fitted,
        0.0,
    )


# ============================================================
# PROFILE DRAWING
# ============================================================

def draw_profiles(
    ax,
    counts,
    *,
    linestyle,
    linewidth,
    alpha,
    zorder,
    show_clipped_labels,
):
    counts = np.asarray(
        counts,
        dtype=float,
    )

    n = float(
        counts.sum()
    )

    if n <= 0:
        raise RuntimeError(
            "Empty histogram."
        )

    fraction = (
        counts / n
    )

    for index, colour in enumerate(
        COLOURS
    ):
        fitted = loess(
            CENTRES[1:-1],
            fraction[
                index,
                1:-1,
            ],
        )

        ax.plot(
            CENTRES[1:-1],
            fitted,
            color=colour,
            alpha=alpha,
            linewidth=linewidth,
            linestyle=linestyle,
            zorder=zorder,
        )

        # Endpoint bins are unsmoothed, matching the source.
        for endpoint in (
            0,
            len(CENTRES) - 1,
        ):
            value = fraction[
                index,
                endpoint,
            ]

            if value <= 0:
                continue

            ax.plot(
                [
                    CENTRES[endpoint],
                    CENTRES[endpoint],
                ],
                [
                    0,
                    min(
                        value,
                        CAP,
                    ),
                ],
                color=colour,
                alpha=alpha,
                linewidth=linewidth,
                linestyle=linestyle,
                zorder=zorder,
            )

        # Mark every clipped bin. Only the focal TN curve gets
        # numerical labels, keeping the overlay readable.
        for bin_id in np.flatnonzero(
            fraction[index] > CAP
        ):
            value = fraction[
                index,
                bin_id,
            ]

            ax.plot(
                CENTRES[bin_id],
                CAP * 0.98,
                "^",
                color=colour,
                markersize=5,
                zorder=zorder + 0.2,
                clip_on=False,
            )

            if show_clipped_labels:
                ax.annotate(
                    f"{100 * value:.2f}%",
                    (
                        CENTRES[bin_id],
                        CAP,
                    ),
                    xytext=(
                        (index - 1.5) * 9,
                        5,
                    ),
                    textcoords="offset points",
                    rotation=90,
                    ha="center",
                    va="bottom",
                    fontsize=9,
                    color=colour,
                    annotation_clip=False,
                )


# ============================================================
# SINGLE-PANEL SCIENCE-STYLE FIGURE
# ============================================================

control_counts = hist[0]
left_counts = hist[1]

fig = plt.figure(
    figsize=(8.2, 6.2)
)

ax = fig.add_axes(
    [
        0.14,
        0.19,
        0.78,
        0.66,
    ]
)

# Control reference first.
draw_profiles(
    ax,
    control_counts,
    linestyle=(0, (4, 2)),
    linewidth=1.45,
    alpha=0.43,
    zorder=2,
    show_clipped_labels=False,
)

# TN Left on top.
draw_profiles(
    ax,
    left_counts,
    linestyle="-",
    linewidth=2.15,
    alpha=0.92,
    zorder=4,
    show_clipped_labels=True,
)


# ============================================================
# AXIS FORMATTING
# ============================================================

ax.axvline(
    0.5,
    color="#555555",
    ls="--",
    lw=0.8,
    alpha=0.7,
    zorder=1,
)

ax.set_xlim(
    -0.015,
    1.015,
)

ax.set_ylim(
    0,
    CAP,
)

ax.set_xticks(
    [
        0,
        0.25,
        0.5,
        0.75,
        1,
    ]
)

ax.set_yticks(
    [
        0,
        0.005,
        0.010,
        0.015,
    ]
)

ax.yaxis.set_major_formatter(
    PercentFormatter(
        1,
        decimals=1,
    )
)

ax.set_xlabel(
    "Observed ALT-read fraction",
    fontsize=12,
)

ax.set_ylabel(
    "Sites per 0.02 AF bin",
    fontsize=12,
)

ax.tick_params(
    labelsize=10.5,
)

ax.grid(
    axis="y",
    color="#E5E5E5",
    lw=0.65,
)

ax.set_axisbelow(
    True
)

ax.spines["top"].set_visible(
    False
)

ax.spines["right"].set_visible(
    False
)


# ============================================================
# LEGENDS
# ============================================================

genotype_handles = [
    Patch(
        facecolor=colour,
        edgecolor="none",
        label=label,
    )
    for colour, label
    in zip(
        COLOURS,
        LABELS,
    )
]

genotype_legend = ax.legend(
    handles=genotype_handles,
    loc="upper center",
    bbox_to_anchor=(
        0.5,
        1.20,
    ),
    frameon=False,
    ncol=4,
    fontsize=9.5,
    title="Original VCF genotype",
    title_fontsize=10,
)

ax.add_artist(
    genotype_legend
)

style_handles = [
    Line2D(
        [0],
        [0],
        color="#555555",
        lw=1.45,
        alpha=0.50,
        linestyle=(0, (4, 2)),
        label=f"Control: {CONTROL}",
    ),
    Line2D(
        [0],
        [0],
        color="#555555",
        lw=2.15,
        alpha=0.95,
        linestyle="-",
        label="Gynandromorph Left",
    ),
]

ax.legend(
    handles=style_handles,
    loc="upper right",
    frameon=False,
    fontsize=9.5,
)


# ============================================================
# PANEL LABEL / SAMPLE COUNTS
# ============================================================

ax.text(
    -0.12,
    1.08,
    "A",
    transform=ax.transAxes,
    fontsize=17,
    fontweight="bold",
    ha="left",
    va="bottom",
)

ax.text(
    0.02,
    0.94,
    (
        f"Control n = {int(control_counts.sum()):,}"
        f"   |   Left n = {int(left_counts.sum()):,}"
    ),
    transform=ax.transAxes,
    fontsize=9.5,
    color="#444444",
    ha="left",
    va="top",
    bbox=dict(
        facecolor="white",
        edgecolor="none",
        alpha=0.80,
    ),
)


# ============================================================
# VECTOR OUTPUT
# ============================================================

stem = (
    "ERR2298208_vs_TNLeft_"
    "original_GT_LOESS_overlay"
)

pdf_file = (
    OUT
    / f"{stem}.pdf"
)

png_file = (
    OUT
    / f"{stem}.png"
)

fig.savefig(
    pdf_file,
)

fig.savefig(
    png_file,
    dpi=220,
)

plt.close(fig)


# ============================================================
# REPORT
# ============================================================

print(
    "CONTROL_RETAINED",
    int(control_counts.sum()),
    sep="\t",
)

print(
    "TN_LEFT_RETAINED",
    int(left_counts.sum()),
    sep="\t",
)

print(
    "CONTROL_CATEGORY_COUNTS",
    *(
        control_counts
        .sum(axis=1)
        .astype(int)
        .tolist()
    ),
    sep="\t",
)

print(
    "TN_LEFT_CATEGORY_COUNTS",
    *(
        left_counts
        .sum(axis=1)
        .astype(int)
        .tolist()
    ),
    sep="\t",
)

print(
    "WROTE",
    pdf_file.resolve(),
    sep="\t",
)

print(
    "WROTE",
    png_file.resolve(),
    sep="\t",
)
