#!/usr/bin/env python3

from pathlib import Path
import json

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(
    "analysis/genomewide_mixing_corrected_v1"
)

FIGURE_DIR = ROOT / "figures"

PDF = (
    FIGURE_DIR
    / "reciprocal_mixture_snp_recaller_narrative.pdf"
)

AUDIT = (
    FIGURE_DIR
    / "reciprocal_mixture_snp_recaller_narrative.audit.json"
)

BLUE = "#0072B2"
ORANGE = "#D55E00"
PURPLE = "#6A3D9A"
GREEN = "#009E73"
INK = "#202124"
GREY = "#6B737B"
PALE = "#F4F6F8"


plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "axes.titlesize": 20,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def new_page(title, subtitle):
    figure = plt.figure(
        figsize=(13.333, 7.5),
        facecolor="white",
    )

    axis = figure.add_axes([0, 0, 1, 1])

    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")

    axis.text(
        0.055,
        0.94,
        title,
        fontsize=24,
        fontweight="bold",
        color=INK,
        va="top",
    )

    axis.text(
        0.057,
        0.875,
        subtitle,
        fontsize=12.5,
        color=GREY,
        va="top",
    )

    return figure, axis


def add_box(
    axis,
    x,
    y,
    width,
    height,
    title,
    body,
    color,
    fill=PALE,
    title_size=14,
    body_size=11,
):
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle=(
            "round,pad=0.012,"
            "rounding_size=0.015"
        ),
        linewidth=1.6,
        edgecolor=color,
        facecolor=fill,
    )

    axis.add_patch(patch)

    axis.text(
        x + 0.025 * width,
        y + height - 0.08 * height,
        title,
        fontsize=title_size,
        fontweight="bold",
        color=color,
        va="top",
    )

    if body:
        axis.text(
            x + 0.025 * width,
            y + height - 0.29 * height,
            body,
            fontsize=body_size,
            color=INK,
            va="top",
            linespacing=1.35,
        )


def add_arrow(
    axis,
    x1,
    y1,
    x2,
    y2,
    color=GREY,
):
    axis.add_patch(
        FancyArrowPatch(
            (x1, y1),
            (x2, y2),
            arrowstyle="-|>",
            mutation_scale=17,
            linewidth=1.8,
            color=color,
        )
    )


def add_footer(axis, page_number):
    axis.text(
        0.055,
        0.035,
        "Reciprocal-mixture SNP recaller",
        fontsize=9,
        color=GREY,
    )

    axis.text(
        0.945,
        0.035,
        str(page_number),
        fontsize=9,
        color=GREY,
        ha="right",
    )


def finish_page(pdf, figure, axis, number):
    add_footer(axis, number)
    pdf.savefig(figure)
    plt.close(figure)


def page_one(pdf):
    figure, axis = new_page(
        "1. Mixing corrupts the ordinary focal GT calls",
        (
            "A biological homozygote can appear heterozygous "
            "because approximately 20% of its reads came from "
            "the other tissue."
        ),
    )

    add_box(
        axis,
        0.06,
        0.57,
        0.36,
        0.22,
        "Underlying biological state",
        "Left:   0/0\nRight:  1/1",
        BLUE,
        "#EAF4FA",
        body_size=17,
    )

    add_arrow(
        axis,
        0.43,
        0.68,
        0.56,
        0.68,
        PURPLE,
    )

    axis.text(
        0.495,
        0.73,
        "~20% reciprocal\nread transfer",
        ha="center",
        va="bottom",
        fontsize=12,
        fontweight="bold",
        color=PURPLE,
    )

    add_box(
        axis,
        0.57,
        0.57,
        0.37,
        0.22,
        "Observed focal VCF fields",
        (
            "Left:   GT=0/1   AD=32,8\n"
            "Right:  GT=0/1   AD=9,31"
        ),
        ORANGE,
        "#FFF1E8",
        body_size=17,
    )

    add_box(
        axis,
        0.12,
        0.20,
        0.76,
        0.20,
        "First methodological decision",
        (
            "Discard the ordinary focal GT as biological truth.\n"
            "Retain REF, ALT and AD_REF,AD_ALT as the primary "
            "observations to be explained."
        ),
        PURPLE,
        "#F3EEFA",
        body_size=15,
    )

    finish_page(pdf, figure, axis, 1)


def page_two(pdf):
    figure, axis = new_page(
        "2. Define usable read evidence and calibrate its noise",
        (
            "The recaller separates genuine reciprocal mixing "
            "from sequencing error, allele bias and read-count "
            "overdispersion."
        ),
    )

    add_box(
        axis,
        0.055,
        0.57,
        0.25,
        0.24,
        "Broad Optix VCF",
        (
            "Biallelic SNPs\n"
            "Left and Right AD counts\n"
            "Unmixed WGS controls"
        ),
        BLUE,
        "#EAF4FA",
    )

    add_arrow(
        axis,
        0.305,
        0.69,
        0.375,
        0.69,
    )

    add_box(
        axis,
        0.375,
        0.57,
        0.25,
        0.24,
        "Fitting-site filter",
        (
            "AD_REF + AD_ALT >= 10\n"
            "in both focal samples\n"
            "Exclude unusable sites"
        ),
        ORANGE,
        "#FFF1E8",
    )

    add_arrow(
        axis,
        0.625,
        0.69,
        0.695,
        0.69,
    )

    add_box(
        axis,
        0.695,
        0.57,
        0.25,
        0.24,
        "Control calibration",
        (
            "Sequencing error\n"
            "REF/ALT bias\n"
            "Beta-binomial dispersion"
        ),
        GREEN,
        "#EAF7F2",
    )

    add_box(
        axis,
        0.11,
        0.20,
        0.78,
        0.21,
        "Why this step is necessary",
        (
            "Even an unmixed 0/1 genotype does not always give "
            "exactly 50% ALT reads. The control-calibrated "
            "emission model defines the amount and direction "
            "of ordinary read-count variation before mixing "
            "is estimated."
        ),
        GREEN,
        "#F2FAF7",
        body_size=13,
    )

    finish_page(pdf, figure, axis, 2)


def page_three(pdf):
    figure, axis = new_page(
        "3. Fit reciprocal mixing while the genotypes remain unknown",
        (
            "Every site is evaluated under all nine possible "
            "diploid Left/Right genotype pairs."
        ),
    )

    states = ["0/0", "0/1", "1/1"]

    x_start = 0.105
    y_start = 0.19
    cell_width = 0.13
    cell_height = 0.14

    axis.text(
        0.37,
        0.72,
        "Right biological genotype",
        fontsize=13,
        fontweight="bold",
        color=ORANGE,
        ha="center",
    )

    axis.text(
        0.055,
        0.42,
        "Left biological\ngenotype",
        fontsize=13,
        fontweight="bold",
        color=BLUE,
        ha="center",
        va="center",
        rotation=90,
    )

    for column, right_state in enumerate(states):
        axis.text(
            x_start
            + (column + 1) * cell_width
            + cell_width / 2,
            0.665,
            right_state,
            fontsize=13,
            fontweight="bold",
            color=ORANGE,
            ha="center",
        )

    for row, left_state in enumerate(states):
        y = (
            y_start
            + (2 - row) * cell_height
        )

        axis.text(
            x_start + cell_width / 2,
            y + cell_height / 2,
            left_state,
            fontsize=13,
            fontweight="bold",
            color=BLUE,
            ha="center",
            va="center",
        )

        for column, right_state in enumerate(states):
            x = (
                x_start
                + (column + 1) * cell_width
            )

            add_box(
                axis,
                x,
                y,
                cell_width - 0.012,
                cell_height - 0.012,
                (
                    f"({left_state}, "
                    f"{right_state})"
                ),
                "",
                PURPLE,
                "#F3EEFA",
                title_size=10.5,
            )

    add_box(
        axis,
        0.62,
        0.48,
        0.32,
        0.21,
        "Expected fractions after mixing",
        (
            r"$q_L=(1-c_L)p_L+c_Lp_R$"
            "\n\n"
            r"$q_R=c_Rp_L+(1-c_R)p_R$"
        ),
        PURPLE,
        "#F3EEFA",
        body_size=15,
    )

    add_box(
        axis,
        0.62,
        0.19,
        0.32,
        0.21,
        "Joint likelihood",
        (
            "Use a beta-binomial likelihood for every "
            "genotype pair, then average over all nine "
            "pairs with equal prior weight.\n\n"
            "Optimize c_L and c_R over 0-0.45 using "
            "multiple starting points."
        ),
        ORANGE,
        "#FFF1E8",
        body_size=11.5,
    )

    finish_page(pdf, figure, axis, 3)


def page_four(pdf):
    figure, axis = new_page(
        "4. Recall the SNPs, then validate the fitted mixing",
        (
            "The intuitive opposing-allele calculation becomes "
            "possible only after the underlying genotype pairs "
            "have been recalled."
        ),
    )

    stages = [
        (
            "Fitted mixing",
            r"$\hat c_L,\hat c_R$",
            ORANGE,
        ),
        (
            "Nine posteriors",
            (
                r"$P(g_L,g_R\mid AD,$"
                "\n"
                r"$\hat c_L,\hat c_R)$"
            ),
            PURPLE,
        ),
        (
            "Recall pair",
            (
                "Choose maximum-\n"
                "posterior state"
            ),
            GREEN,
        ),
        (
            "Confidence rule",
            (
                r"$P_{\max}\geq0.95$"
                "\n"
                "No raw-GT fallback"
            ),
            BLUE,
        ),
    ]

    x_positions = [
        0.055,
        0.285,
        0.515,
        0.745,
    ]

    for index, (
        title,
        body,
        color,
    ) in enumerate(stages):
        add_box(
            axis,
            x_positions[index],
            0.59,
            0.19,
            0.20,
            title,
            body,
            color,
            PALE,
            title_size=12,
            body_size=12,
        )

        if index < len(stages) - 1:
            add_arrow(
                axis,
                x_positions[index] + 0.19,
                0.69,
                x_positions[index + 1],
                0.69,
            )

    add_box(
        axis,
        0.07,
        0.25,
        0.39,
        0.20,
        "Recalled (0/0, 1/1)",
        (
            r"Left foreign $=AD_{L,ALT}/D_L$"
            "\n"
            r"Right foreign $=AD_{R,REF}/D_R$"
        ),
        BLUE,
        "#EAF4FA",
        body_size=13,
    )

    add_box(
        axis,
        0.54,
        0.25,
        0.39,
        0.20,
        "Recalled (1/1, 0/0)",
        (
            r"Left foreign $=AD_{L,REF}/D_L$"
            "\n"
            r"Right foreign $=AD_{R,ALT}/D_R$"
        ),
        ORANGE,
        "#FFF1E8",
        body_size=13,
    )

    axis.text(
        0.50,
        0.16,
        (
            "These read fractions validate and explain the "
            "fitted mixing; they did not independently create "
            "the precise estimates."
        ),
        ha="center",
        va="center",
        fontsize=13,
        fontweight="bold",
        color=PURPLE,
    )

    axis.text(
        0.50,
        0.095,
        (
            "Original Optix recaller: "
            "20.87% Right->Left, 20.99% Left->Right"
            "     |     "
            "Later genome-wide fit: "
            "19.70% Right->Left, 20.24% Left->Right"
        ),
        ha="center",
        va="center",
        fontsize=11.5,
        color=INK,
    )

    finish_page(pdf, figure, axis, 4)


def main():
    FIGURE_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    with PdfPages(PDF) as pdf:
        page_one(pdf)
        page_two(pdf)
        page_three(pdf)
        page_four(pdf)

    audit = {
        "pages": 4,
        "purpose": (
            "Narrative of the reciprocal-mixture "
            "SNP-recalling protocol"
        ),
        "focal_raw_GT_used_as_truth": False,
        "minimum_AD_depth_each_side": 10,
        "latent_genotype_pairs": 9,
        "pair_prior": "equal",
        "mixing_bounds": [0.0, 0.45],
        "posterior_threshold": 0.95,
        "raw_GT_fallback": False,
        "pdf": str(PDF.resolve()),
    }

    AUDIT.write_text(
        json.dumps(
            audit,
            indent=2,
        )
        + "\n"
    )

    print(
        "NARRATIVE_PDF\t"
        f"{PDF.resolve()}"
    )

    print(
        "AUDIT\t"
        f"{AUDIT.resolve()}"
    )

    print("PAGES\t4")


if __name__ == "__main__":
    main()
