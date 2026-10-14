#!/usr/bin/env python3

from pathlib import Path
import csv
import math

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# ============================================================
# OUTPUT
# ============================================================

OUT = Path(
    "analysis/figures/"
    "depth_dispersion_correct_founder_models"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# EMPIRICAL TN DEPTH SOURCE
#
# This is the real TN Left / Right depth audit previously used.
# ============================================================

DEPTH_FILE = Path(
    "analysis/original_502/checkpoints/"
    "optix_R_BA90_focal_depth_audit.tsv"
)


# ============================================================
# MODEL PARAMETERS
#
# These match the previous founder-model simulations:
#
# beta-binomial concentration = 22
# sequencing error            = 0.001
# reciprocal mixing           = 0 or 0.20
# ============================================================

CONCENTRATION = 22.0
ERROR = 0.001
PREDICTIVE_MASS = 0.50

DEPTHS = [
    10,
    20,
    30,
    50,
    100,
    200,
]


# ============================================================
# CORRECT FOUNDER STATES
# ============================================================

# Two-founder:
# Both sides contain the SAME two haplotypes A and B.
#
# Therefore at a biallelic site the two sides have the
# same underlying dosage:
#
#   A=B=REF       -> (0, 0)
#   A!=B          -> (0.5, 0.5)
#   A=B=ALT       -> (1, 1)
#
TWO_FOUNDER_STATES = [
    (0.0, 0.0),
    (0.5, 0.5),
    (1.0, 1.0),
]


# Four-founder:
# Left = AB
# Right = CD
#
# The four founders are independent, so each side can
# independently have dosage 0, 0.5, or 1.
#
FOUR_FOUNDER_STATES = [
    (left, right)
    for left in (0.0, 0.5, 1.0)
    for right in (0.0, 0.5, 1.0)
]


# ============================================================
# READ REAL PAIRED TN DEPTHS
# ============================================================

def read_empirical_depths(path):

    if not path.exists():
        raise FileNotFoundError(
            f"Missing empirical depth table: {path}"
        )

    paired = {}

    with open(path, "rt") as fh:

        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )

        required = {
            "chrom",
            "pos",
            "side",
            "DP",
        }

        missing = required - set(reader.fieldnames)

        if missing:
            raise RuntimeError(
                f"{path} is missing columns: "
                + ", ".join(sorted(missing))
            )

        for row in reader:

            side = row["side"].strip()

            if side not in {
                "Left",
                "Right",
            }:
                continue

            try:
                dp = float(row["DP"])
            except (ValueError, TypeError):
                continue

            if not np.isfinite(dp):
                continue

            if dp <= 0:
                continue

            key = (
                row["chrom"],
                row["pos"],
            )

            paired.setdefault(
                key,
                {}
            )[side] = dp


    left = []
    right = []

    for value in paired.values():

        if (
            "Left" not in value
            or "Right" not in value
        ):
            continue

        left.append(
            value["Left"]
        )

        right.append(
            value["Right"]
        )


    left = np.asarray(
        left,
        dtype=float
    )

    right = np.asarray(
        right,
        dtype=float
    )

    if len(left) == 0:
        raise RuntimeError(
            "No paired empirical Left/Right depths found."
        )

    return left, right


empirical_left, empirical_right = (
    read_empirical_depths(
        DEPTH_FILE
    )
)

empirical_mean = (
    empirical_left
    + empirical_right
) / 2.0


MED_LEFT = float(
    np.median(empirical_left)
)

MED_RIGHT = float(
    np.median(empirical_right)
)

MED_MEAN = float(
    np.median(empirical_mean)
)

Q25_MEAN = float(
    np.percentile(
        empirical_mean,
        25
    )
)

Q75_MEAN = float(
    np.percentile(
        empirical_mean,
        75
    )
)


print(
    f"PAIRED SITES: {len(empirical_left)}"
)

print(
    f"LEFT MEDIAN DP: {MED_LEFT:.1f}"
)

print(
    f"RIGHT MEDIAN DP: {MED_RIGHT:.1f}"
)

print(
    f"PAIRED-MEAN MEDIAN DP: {MED_MEAN:.1f}"
)

print(
    f"PAIRED-MEAN IQR: "
    f"{Q25_MEAN:.1f} - {Q75_MEAN:.1f}"
)


# ============================================================
# MODEL FUNCTIONS
# ============================================================

def mix_centroid(
    p_left,
    p_right,
    mixing
):
    """
    Symmetric reciprocal read mixing.

    qL = (1-m) pL + m pR
    qR = m pL + (1-m) pR

    Then apply sequencing-error mapping used previously.
    """

    mixed_left = (
        (1.0 - mixing) * p_left
        + mixing * p_right
    )

    mixed_right = (
        mixing * p_left
        + (1.0 - mixing) * p_right
    )

    q_left = (
        ERROR
        + (1.0 - 2.0 * ERROR)
        * mixed_left
    )

    q_right = (
        ERROR
        + (1.0 - 2.0 * ERROR)
        * mixed_right
    )

    return q_left, q_right


def beta_binomial_pmf(
    n,
    p,
    concentration
):
    """
    Exact beta-binomial predictive PMF.
    """

    p = min(
        max(
            float(p),
            1e-10
        ),
        1.0 - 1e-10
    )

    alpha = (
        p
        * concentration
    )

    beta = (
        (1.0 - p)
        * concentration
    )

    logp = np.empty(
        n + 1,
        dtype=float
    )

    for k in range(
        n + 1
    ):

        logp[k] = (
            math.lgamma(n + 1)
            - math.lgamma(k + 1)
            - math.lgamma(n - k + 1)

            + math.lgamma(
                k + alpha
            )

            + math.lgamma(
                n - k + beta
            )

            - math.lgamma(
                n + alpha + beta
            )

            + math.lgamma(
                alpha + beta
            )

            - math.lgamma(alpha)
            - math.lgamma(beta)
        )


    logp -= np.max(logp)

    pmf = np.exp(logp)

    pmf /= pmf.sum()

    return pmf


def beta_binomial_af_sd(
    p,
    depth,
    concentration
):
    """
    Predictive SD of read AF under the beta-binomial model.
    """

    n = float(depth)

    variance = (
        p
        * (1.0 - p)
        * (n + concentration)
        / (
            n
            * (concentration + 1.0)
        )
    )

    return math.sqrt(
        variance
    )


def binomial_af_sd(
    p,
    depth
):
    return math.sqrt(
        p
        * (1.0 - p)
        / float(depth)
    )


# ============================================================
# 50% HIGHEST-PREDICTIVE-MASS RING
#
# IMPORTANT:
#   outline only
#   no filled probability blocks
# ============================================================

def draw_predictive_ring(
    ax,
    q_left,
    q_right,
    depth_left,
    depth_right,
    colour,
    linewidth=0.65,
    alpha=0.90,
    zorder=3,
):
    """
    Smooth vector 50% predictive ellipse.

    Variance follows the same beta-binomial model used in the
    founder simulations, but the display is deliberately a smooth
    continuous ellipse rather than a discrete read-count contour.
    """

    # Predictive AF SD under beta-binomial sampling.
    sx = beta_binomial_af_sd(
        q_left,
        depth_left,
        CONCENTRATION
    )

    sy = beta_binomial_af_sd(
        q_right,
        depth_right,
        CONCENTRATION
    )

    # 50% probability contour of a 2-D normal approximation:
    #
    # chi-square(df=2) median = 1.38629436112
    #
    # radius in SD units = sqrt(1.38629436112)
    radius = math.sqrt(
        1.38629436112
    )

    # Handle effectively fixed states at 0 / 1 cleanly.
    # Sequencing error means these are not mathematically zero,
    # but they should remain very narrow in the cartoon.
    min_half_width = 0.004

    half_x = max(
        radius * sx,
        min_half_width
    )

    half_y = max(
        radius * sy,
        min_half_width
    )

    theta = np.linspace(
        0.0,
        2.0 * np.pi,
        240
    )

    x = (
        q_left
        + half_x * np.cos(theta)
    )

    y = (
        q_right
        + half_y * np.sin(theta)
    )

    ax.plot(
        x,
        y,
        color=colour,
        linewidth=linewidth,
        alpha=alpha,
        zorder=zorder,
        solid_joinstyle="round",
        solid_capstyle="round",
    )


# ============================================================
# COLOURS
#
# Restrained viridis depth sequence.
# Empirical depth keeps the old model-magenta identity.
# ============================================================

viridis = plt.get_cmap(
    "viridis"
)

depth_colours = {
    depth: viridis(value)
    for depth, value
    in zip(
        DEPTHS,
        np.linspace(
            0.10,
            0.82,
            len(DEPTHS)
        )
    )
}

EMPIRICAL_COLOUR = "#F000B8"


# ============================================================
# FIGURE 1
#
# DEPTH vs HET AF SD
# ============================================================

depth_grid = np.arange(
    5,
    221,
    1
)


binomial_sd = np.asarray([
    binomial_af_sd(
        0.5,
        depth
    )
    for depth
    in depth_grid
])


model_sd = np.asarray([
    beta_binomial_af_sd(
        0.5,
        depth,
        CONCENTRATION
    )
    for depth
    in depth_grid
])


fig, ax = plt.subplots(
    figsize=(
        7.0,
        5.0
    )
)


# Pure read-sampling expectation
ax.plot(
    depth_grid,
    binomial_sd,
    color="#777777",
    linewidth=1.4,
    linestyle="--",
    label="Binomial read sampling",
    zorder=1,
)


# Model actually used for founder simulations
ax.plot(
    depth_grid,
    model_sd,
    color="#222222",
    linewidth=2.2,
    label=(
        "Beta-binomial model "
        "(concentration = 22)"
    ),
    zorder=2,
)


# Reference depths
for depth in DEPTHS:

    ax.scatter(
        depth,
        beta_binomial_af_sd(
            0.5,
            depth,
            CONCENTRATION
        ),
        s=42,
        facecolor=depth_colours[depth],
        edgecolor="black",
        linewidth=0.5,
        zorder=4,
    )


# Real-sample depth
empirical_sd = (
    beta_binomial_af_sd(
        0.5,
        MED_MEAN,
        CONCENTRATION
    )
)


ax.axvline(
    MED_MEAN,
    color=EMPIRICAL_COLOUR,
    linewidth=1.4,
    alpha=0.85,
    zorder=1,
)


# IQR shown as lines, NOT a filled band.
ax.axvline(
    Q25_MEAN,
    color=EMPIRICAL_COLOUR,
    linewidth=0.8,
    linestyle=":",
    alpha=0.65,
    zorder=1,
)

ax.axvline(
    Q75_MEAN,
    color=EMPIRICAL_COLOUR,
    linewidth=0.8,
    linestyle=":",
    alpha=0.65,
    zorder=1,
)


ax.scatter(
    MED_MEAN,
    empirical_sd,
    s=80,
    facecolor=EMPIRICAL_COLOUR,
    edgecolor="black",
    linewidth=0.7,
    zorder=5,
)


ax.annotate(
    (
        f"TN empirical median = "
        f"{MED_MEAN:.1f}×\n"
        f"Left {MED_LEFT:.0f}×; "
        f"Right {MED_RIGHT:.0f}×"
    ),
    xy=(
        MED_MEAN,
        empirical_sd
    ),
    xytext=(
        MED_MEAN + 13,
        empirical_sd + 0.012
    ),
    fontsize=9,
    color=EMPIRICAL_COLOUR,
    arrowprops=dict(
        arrowstyle="-",
        color=EMPIRICAL_COLOUR,
        linewidth=0.8,
    ),
)


ax.set_xlim(
    0,
    220
)

ax.set_ylim(
    0,
    max(binomial_sd) * 1.05
)

ax.set_xlabel(
    "Read depth (×)",
    fontsize=11
)

ax.set_ylabel(
    "SD of heterozygous ALT-read fraction",
    fontsize=11
)

ax.set_title(
    "Expected heterozygote AF dispersion with sequencing depth",
    loc="left",
    fontsize=13,
    fontweight="bold",
)


ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.tick_params(
    direction="out",
    length=3,
    width=0.8,
)


ax.legend(
    frameon=False,
    fontsize=9,
)


fig.tight_layout()


depth_pdf = (
    OUT
    / "01_depth_vs_hetAF_sd.pdf"
)

depth_png = (
    OUT
    / "01_depth_vs_hetAF_sd.png"
)


fig.savefig(
    depth_pdf,
    bbox_inches="tight"
)

fig.savefig(
    depth_png,
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)


# ============================================================
# FIGURE 2
#
# 2 x 2:
#
#     Two-founder        Four-founder
#
#     0% mix             0% mix
#
#     20% mix            20% mix
#
# ============================================================

fig, axes = plt.subplots(
    nrows=2,
    ncols=2,
    figsize=(
        10.8,
        10.0
    ),
    sharex=True,
    sharey=True,
)


models = [
    (
        "Two-founder",
        "Same two haplotypes",
        TWO_FOUNDER_STATES,
    ),
    (
        "Four-founder",
        "Four independent haplotypes",
        FOUR_FOUNDER_STATES,
    ),
]


mixings = [
    0.0,
    0.20,
]


# Store explicit centroid audit.
audit_rows = []


for row_index, mixing in enumerate(
    mixings
):

    for col_index, (
        model_name,
        model_description,
        states,
    ) in enumerate(models):

        ax = axes[
            row_index,
            col_index
        ]


        # ---------------------------------------------
        # Very light reference diagonal.
        # ---------------------------------------------

        ax.plot(
            [0, 1],
            [0, 1],
            color="#D0D0D0",
            linewidth=0.7,
            linestyle="--",
            zorder=0,
        )


        # ---------------------------------------------
        # Draw LOW depth first -> outer rings first.
        # Then progressively higher depth.
        # ---------------------------------------------

        for (
            p_left,
            p_right
        ) in states:

            q_left, q_right = (
                mix_centroid(
                    p_left,
                    p_right,
                    mixing
                )
            )


            audit_rows.append({
                "model": model_name,
                "mixing": mixing,
                "unmixed_left": p_left,
                "unmixed_right": p_right,
                "expected_left": q_left,
                "expected_right": q_right,
            })


            for depth in DEPTHS:

                draw_predictive_ring(
                    ax=ax,
                    q_left=q_left,
                    q_right=q_right,
                    depth_left=depth,
                    depth_right=depth,
                    colour=depth_colours[depth],
                    linewidth=0.55,
                    alpha=0.90,
                    zorder=2,
                )


            # -----------------------------------------
            # REAL empirical paired median depth:
            # Left = real Left median
            # Right = real Right median
            # -----------------------------------------

            draw_predictive_ring(
                ax=ax,
                q_left=q_left,
                q_right=q_right,
                depth_left=MED_LEFT,
                depth_right=MED_RIGHT,
                colour=EMPIRICAL_COLOUR,
                linewidth=1.05,
                alpha=1.0,
                zorder=5,
            )


            # Small centroid only.
            ax.scatter(
                q_left,
                q_right,
                s=13,
                facecolor="black",
                edgecolor="white",
                linewidth=0.35,
                zorder=7,
            )


        ax.set_xlim(
            -0.035,
            1.035
        )

        ax.set_ylim(
            -0.035,
            1.035
        )

        ax.set_aspect(
            "equal",
            adjustable="box"
        )


        ax.set_xticks([
            0.0,
            0.25,
            0.5,
            0.75,
            1.0,
        ])

        ax.set_yticks([
            0.0,
            0.25,
            0.5,
            0.75,
            1.0,
        ])


        ax.tick_params(
            direction="out",
            length=3,
            width=0.8,
            labelsize=9,
        )


        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)


        if row_index == 0:

            ax.set_title(
                (
                    f"{model_name}\n"
                    f"{model_description}"
                ),
                fontsize=12,
                pad=10,
            )


        if row_index == 1:

            ax.set_xlabel(
                "Left ALT-read fraction",
                fontsize=10,
            )


        if col_index == 0:

            ax.set_ylabel(
                "Right ALT-read fraction",
                fontsize=10,
            )


# Row labels
fig.text(
    0.015,
    0.735,
    "0% reciprocal\nread mixing",
    rotation=90,
    va="center",
    ha="left",
    fontsize=11,
    fontweight="bold",
)

fig.text(
    0.015,
    0.285,
    "20% reciprocal\nread mixing",
    rotation=90,
    va="center",
    ha="left",
    fontsize=11,
    fontweight="bold",
)


# ============================================================
# DEPTH LEGEND
# ============================================================

legend_handles = []

for depth in DEPTHS:

    legend_handles.append(
        Line2D(
            [0],
            [0],
            color=depth_colours[depth],
            linewidth=2.0,
            label=f"{depth}×",
        )
    )


legend_handles.append(
    Line2D(
        [0],
        [0],
        color=EMPIRICAL_COLOUR,
        linewidth=2.4,
        label=(
            f"TN empirical "
            f"({MED_LEFT:.0f}× L / "
            f"{MED_RIGHT:.0f}× R)"
        ),
    )
)


fig.legend(
    handles=legend_handles,
    title="50% predictive-region depth",
    frameon=False,
    ncol=7,
    fontsize=8.5,
    title_fontsize=9,
    loc="lower center",
    bbox_to_anchor=(
        0.52,
        0.015
    ),
)


fig.suptitle(
    (
        "Founder-model AF dispersion across sequencing depth"
    ),
    fontsize=15,
    fontweight="bold",
    y=0.985,
)


fig.text(
    0.52,
    0.952,
    (
        "Outline only: 50% beta-binomial predictive regions; "
        "concentration 22; sequencing error 0.001"
    ),
    ha="center",
    va="top",
    fontsize=9,
)


fig.subplots_adjust(
    left=0.09,
    right=0.98,
    top=0.91,
    bottom=0.11,
    wspace=0.18,
    hspace=0.20,
)


model_pdf = (
    OUT
    / "02_two_vs_four_founder_depth_rings_0_20mix.pdf"
)

model_png = (
    OUT
    / "02_two_vs_four_founder_depth_rings_0_20mix.png"
)


fig.savefig(
    model_pdf,
    bbox_inches="tight"
)

fig.savefig(
    model_png,
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)


# ============================================================
# WRITE CENTROID AUDIT
#
# This makes the model geometry explicit so we can inspect it.
# ============================================================

audit_file = (
    OUT
    / "02_model_centroid_audit.tsv"
)


with open(
    audit_file,
    "wt",
    newline=""
) as fh:

    fieldnames = [
        "model",
        "mixing",
        "unmixed_left",
        "unmixed_right",
        "expected_left",
        "expected_right",
    ]

    writer = csv.DictWriter(
        fh,
        fieldnames=fieldnames,
        delimiter="\t",
    )

    writer.writeheader()

    writer.writerows(
        audit_rows
    )


print()
print(
    f"WROTE: {depth_pdf}"
)
print(
    f"WROTE: {depth_png}"
)
print(
    f"WROTE: {model_pdf}"
)
print(
    f"WROTE: {model_png}"
)
print(
    f"WROTE: {audit_file}"
)
