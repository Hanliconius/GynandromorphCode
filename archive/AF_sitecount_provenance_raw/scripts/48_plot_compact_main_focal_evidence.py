#!/usr/bin/env python3

from __future__ import annotations

import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


PROJECT = Path(
    "/scratch/martinlab/joe/GynandroCalling"
)

CHECKPOINTS = (
    PROJECT
    / "analysis/original_502/checkpoints"
)

VCF_DIR = (
    PROJECT
    / "external/krista_genotype_plots_original/Genotype_Plots"
)

FIGURE_DIR = (
    PROJECT
    / "analysis/figures/main_compact_focal_evidence"
)

OPTIX_SITES = (
    CHECKPOINTS
    / "optix_comparative_informative_sites_BA80.tsv"
)

YB_SITES = (
    CHECKPOINTS
    / "ivory_comparative_informative_sites_BA80.tsv"
)

OUTPUT_LONG = (
    CHECKPOINTS
    / "focal_calls_at_independent_comparative_sites.tsv"
)

OUTPUT_SUMMARY = (
    CHECKPOINTS
    / "focal_calls_at_independent_comparative_sites_summary.tsv"
)

OUTPUT_PDF = (
    FIGURE_DIR
    / "focal_left_right_at_independent_sites_compact.pdf"
)

OUTPUT_SVG = (
    FIGURE_DIR
    / "focal_left_right_at_independent_sites_compact.svg"
)

LEFT_SAMPLE = "TN_Left"
RIGHT_SAMPLE = "TN_Right"

LEFT_COLOUR = "#B12A90"
RIGHT_COLOUR = "#2C9CCB"
HET_LINE_COLOUR = "#FCA636"

BLOCKS = [
    {
        "block": "band",
        "display": "Band",
        "module": "B",
        "width": 0.72,
        "vcf": "optix_band-1112651-1114314.vcf",
    },
    {
        "block": "dennis1",
        "display": "Dennis 1",
        "module": "D",
        "width": 1.00,
        "vcf": "dennis1-1215210-1216124.vcf",
    },
    {
        "block": "dennis2",
        "display": "Dennis 2",
        "module": "D",
        "width": 1.00,
        "vcf": "dennis2-1218557-1221494.vcf",
    },
    {
        "block": "ray1",
        "display": "Ray 1",
        "module": "R",
        "width": 1.00,
        "vcf": "ray1-1156746-1160220.vcf",
    },
    {
        "block": "ray2",
        "display": "Ray 2",
        "module": "R",
        "width": 1.18,
        "vcf": "ray2-1174547-1187833.vcf",
    },
    {
        "block": "Yb",
        "display": "Yb",
        "module": "Yb",
        "width": 1.00,
        "vcf": (
            "Ivory_Hmr15-1627000-"
            "1637000_biallelic.vcf"
        ),
    },
]


def as_float(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None

    if math.isfinite(result):
        return result

    return None


def read_tsv(path):
    with path.open(
        newline=""
    ) as handle:
        return list(
            csv.DictReader(
                handle,
                delimiter="\t",
            )
        )


def write_tsv(
    path,
    rows,
    fieldnames,
):
    with path.open(
        "w",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(rows)


def parse_ad(value):
    if (
        not value
        or value in {
            ".",
            ".,.",
        }
    ):
        return None, None

    parts = value.split(",")

    if len(parts) < 2:
        return None, None

    return (
        as_float(parts[0]),
        as_float(parts[1]),
    )


def normalize_gt(gt):
    if (
        not gt
        or "." in gt
    ):
        return "./."

    alleles = (
        gt.replace("|", "/")
        .split("/")
    )

    if len(alleles) != 2:
        return "./."

    if any(
        allele not in {"0", "1"}
        for allele in alleles
    ):
        return "./."

    return "/".join(
        sorted(alleles)
    )


def associated_fraction(
    ref_depth,
    alt_depth,
    orientation,
):
    if (
        ref_depth is None
        or alt_depth is None
    ):
        return None

    total = (
        ref_depth
        + alt_depth
    )

    if total <= 0:
        return None

    if orientation.startswith(
        "REF_"
    ):
        return (
            ref_depth
            / total
        )

    if orientation.startswith(
        "ALT_"
    ):
        return (
            alt_depth
            / total
        )

    return None


def associated_dosage(
    gt,
    orientation,
):
    if gt == "./.":
        return None

    alt_dosage = sum(
        int(allele)
        for allele in gt.split("/")
    )

    if orientation.startswith(
        "REF_"
    ):
        return (
            2 - alt_dosage
        )

    return alt_dosage


def read_focal_calls(
    vcf_path,
    selected_keys,
):
    calls = {}
    samples = []

    with vcf_path.open() as handle:
        for line in handle:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                samples = (
                    line.rstrip("\n")
                    .split("\t")[9:]
                )

                continue

            if (
                line.startswith("#")
                or not line.strip()
            ):
                continue

            fields = (
                line.rstrip("\n")
                .split("\t")
            )

            chrom = fields[0]
            position = int(fields[1])
            ref = fields[3]
            alt = fields[4]

            if len(ref) != 1:
                continue

            if (
                len(alt) != 1
                or "," in alt
            ):
                continue

            key = (
                chrom,
                position,
                ref,
                alt,
            )

            if key not in selected_keys:
                continue

            format_names = (
                fields[8]
                .split(":")
            )

            format_index = {
                name: index
                for index, name
                in enumerate(format_names)
            }

            sample_index = {
                sample: index
                for index, sample
                in enumerate(samples)
            }

            for focal_sample in (
                LEFT_SAMPLE,
                RIGHT_SAMPLE,
            ):
                if (
                    focal_sample
                    not in sample_index
                ):
                    continue

                values = fields[
                    9
                    + sample_index[
                        focal_sample
                    ]
                ].split(":")

                def get_format(
                    name,
                    default="",
                ):
                    index = (
                        format_index.get(name)
                    )

                    if (
                        index is not None
                        and index < len(values)
                    ):
                        return values[index]

                    return default

                gt = normalize_gt(
                    get_format(
                        "GT",
                        ".",
                    )
                )

                (
                    ref_depth,
                    alt_depth,
                ) = parse_ad(
                    get_format(
                        "AD",
                        "",
                    )
                )

                depth = as_float(
                    get_format(
                        "DP",
                        "",
                    )
                )

                genotype_quality = as_float(
                    get_format(
                        "GQ",
                        "",
                    )
                )

                calls[
                    (
                        key,
                        focal_sample,
                    )
                ] = {
                    "GT": gt,
                    "ref_depth": ref_depth,
                    "alt_depth": alt_depth,
                    "DP": depth,
                    "GQ": genotype_quality,
                    "passes_DP2_GQ20": (
                        gt != "./."
                        and depth is not None
                        and depth >= 2
                        and genotype_quality
                        is not None
                        and genotype_quality >= 20
                    ),
                }

    return calls


def build_selected_sites():
    optix_rows = read_tsv(
        OPTIX_SITES
    )

    yb_rows = read_tsv(
        YB_SITES
    )

    selected_by_block = (
        defaultdict(list)
    )

    for row in optix_rows:
        selected_by_block[
            row["panel"]
        ].append(
            {
                "block": row["panel"],
                "module": row["trait"],
                "support": row["trait"],
                "chrom": row["chrom"],
                "pos": int(row["pos"]),
                "ref": row["ref"],
                "alt": row["alt"],
                "orientation": row[
                    "selected_orientation"
                ],
                "minimum_BA": as_float(
                    row["minimum_BA"]
                ),
                "inheritance_model": row[
                    "inheritance_model"
                ],
            }
        )

    yb_merged = {}

    for row in yb_rows:
        key = (
            row["chrom"],
            int(row["pos"]),
            row["ref"],
            row["alt"],
        )

        if key not in yb_merged:
            yb_merged[key] = {
                "block": "Yb",
                "module": "Yb",
                "support": set(),
                "chrom": row["chrom"],
                "pos": int(row["pos"]),
                "ref": row["ref"],
                "alt": row["alt"],
                "orientation": row[
                    "selected_orientation"
                ],
                "minimum_BA": as_float(
                    row["minimum_BA"]
                ),
                "inheritance_model": (
                    "recessive"
                ),
                "orientation_conflict": (
                    False
                ),
            }

        merged = yb_merged[key]

        merged["support"].add(
            row["trait"]
        )

        merged["minimum_BA"] = max(
            merged["minimum_BA"],
            as_float(
                row["minimum_BA"]
            ),
        )

        if (
            merged["orientation"]
            != row[
                "selected_orientation"
            ]
        ):
            merged[
                "orientation_conflict"
            ] = True

    for row in yb_merged.values():
        if row[
            "orientation_conflict"
        ]:
            print(
                "WARNING: excluding Yb site "
                "with conflicting dorsal/"
                "ventral orientation: "
                f'{row["chrom"]}:'
                f'{row["pos"]}'
            )

            continue

        row["support"] = "+".join(
            sorted(
                row["support"]
            )
        )

        selected_by_block[
            "Yb"
        ].append(row)

    for block in selected_by_block:
        selected_by_block[
            block
        ].sort(
            key=lambda row: row["pos"]
        )

    return selected_by_block


def expected_interpretation(module):
    if module == "B":
        return (
            "left_expected_0;"
            "right_expected_1_or_2;"
            "right_heterozygosity_plausible"
        )

    if module in {
        "D",
        "R",
    }:
        return (
            "both_expected_1_or_2;"
            "left_2_right_1_plausible_"
            "under_DR_vs_BDR_model"
        )

    if module == "Yb":
        return (
            "left_expected_1_heterozygote;"
            "right_expected_2_Yb_homozygote"
        )

    return ""


def main():
    FIGURE_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    CHECKPOINTS.mkdir(
        parents=True,
        exist_ok=True,
    )

    required = [
        OPTIX_SITES,
        YB_SITES,
    ]

    required.extend(
        VCF_DIR / block["vcf"]
        for block in BLOCKS
    )

    missing = [
        str(path)
        for path in required
        if not path.is_file()
    ]

    if missing:
        print(
            "ERROR: required inputs are missing:"
        )

        for path in missing:
            print(path)

        return

    selected_by_block = (
        build_selected_sites()
    )

    evidence_rows = []
    block_layout = []
    current_x = 0.0

    for block_config in BLOCKS:
        block = block_config["block"]
        width = block_config["width"]

        start = current_x
        end = start + width

        sites = selected_by_block.get(
            block,
            [],
        )

        block_layout.append(
            {
                **block_config,
                "start": start,
                "end": end,
                "sites": len(sites),
            }
        )

        if sites:
            margin = min(
                0.09,
                width * 0.10,
            )

            if len(sites) == 1:
                x_values = [
                    (
                        start + end
                    )
                    / 2
                ]

            else:
                step = (
                    width
                    - 2 * margin
                ) / (
                    len(sites) - 1
                )

                x_values = [
                    (
                        start
                        + margin
                        + index * step
                    )
                    for index
                    in range(len(sites))
                ]

            selected_keys = {
                (
                    row["chrom"],
                    row["pos"],
                    row["ref"],
                    row["alt"],
                )
                for row in sites
            }

            focal_calls = read_focal_calls(
                VCF_DIR
                / block_config["vcf"],
                selected_keys,
            )

            for site_index, (
                site,
                x_value,
            ) in enumerate(
                zip(
                    sites,
                    x_values,
                ),
                start=1,
            ):
                key = (
                    site["chrom"],
                    site["pos"],
                    site["ref"],
                    site["alt"],
                )

                focal_sides = (
                    (
                        LEFT_SAMPLE,
                        "Left",
                    ),
                    (
                        RIGHT_SAMPLE,
                        "Right",
                    ),
                )

                for (
                    sample,
                    side,
                ) in focal_sides:
                    call = focal_calls.get(
                        (
                            key,
                            sample,
                        ),
                        {},
                    )

                    gt = call.get(
                        "GT",
                        "./.",
                    )

                    ref_depth = call.get(
                        "ref_depth"
                    )

                    alt_depth = call.get(
                        "alt_depth"
                    )

                    evidence_rows.append(
                        {
                            "block": block,
                            "block_display": (
                                block_config[
                                    "display"
                                ]
                            ),
                            "module": site[
                                "module"
                            ],
                            "support": site[
                                "support"
                            ],
                            "site_order_within_block": (
                                site_index
                            ),
                            "plot_x": x_value,
                            "chrom": site[
                                "chrom"
                            ],
                            "pos": site["pos"],
                            "ref": site["ref"],
                            "alt": site["alt"],
                            "orientation": site[
                                "orientation"
                            ],
                            "minimum_comparative_BA": (
                                site[
                                    "minimum_BA"
                                ]
                            ),
                            "inheritance_model": (
                                site[
                                    "inheritance_model"
                                ]
                            ),
                            "focal_sample": sample,
                            "side": side,
                            "GT": gt,
                            "associated_dosage": (
                                associated_dosage(
                                    gt,
                                    site[
                                        "orientation"
                                    ],
                                )
                            ),
                            "ref_depth": ref_depth,
                            "alt_depth": alt_depth,
                            "associated_allele_fraction": (
                                associated_fraction(
                                    ref_depth,
                                    alt_depth,
                                    site[
                                        "orientation"
                                    ],
                                )
                            ),
                            "DP": call.get("DP"),
                            "GQ": call.get("GQ"),
                            "passes_DP2_GQ20": (
                                call.get(
                                    "passes_DP2_GQ20",
                                    False,
                                )
                            ),
                            "expected_interpretation_not_used_for_selection": (
                                expected_interpretation(
                                    site[
                                        "module"
                                    ]
                                )
                            ),
                        }
                    )

        current_x = end

    fieldnames = [
        "block",
        "block_display",
        "module",
        "support",
        "site_order_within_block",
        "plot_x",
        "chrom",
        "pos",
        "ref",
        "alt",
        "orientation",
        "minimum_comparative_BA",
        "inheritance_model",
        "focal_sample",
        "side",
        "GT",
        "associated_dosage",
        "ref_depth",
        "alt_depth",
        "associated_allele_fraction",
        "DP",
        "GQ",
        "passes_DP2_GQ20",
        "expected_interpretation_not_used_for_selection",
    ]

    write_tsv(
        OUTPUT_LONG,
        evidence_rows,
        fieldnames,
    )

    summary_rows = []

    for layout in block_layout:
        subset = [
            row
            for row in evidence_rows
            if (
                row["block"]
                == layout["block"]
            )
        ]

        for side in (
            "Left",
            "Right",
        ):
            values = [
                row[
                    "associated_allele_fraction"
                ]
                for row in subset
                if (
                    row["side"] == side
                    and row[
                        "associated_allele_fraction"
                    ]
                    is not None
                )
            ]

            sorted_values = sorted(
                values
            )

            summary_rows.append(
                {
                    "block": layout[
                        "block"
                    ],
                    "module": layout[
                        "module"
                    ],
                    "side": side,
                    "independent_sites": (
                        layout["sites"]
                    ),
                    "allele_fraction_observed": (
                        len(values)
                    ),
                    "allele_fraction_median": (
                        sorted_values[
                            len(
                                sorted_values
                            )
                            // 2
                        ]
                        if sorted_values
                        else None
                    ),
                    "allele_fraction_min": (
                        min(values)
                        if values
                        else None
                    ),
                    "allele_fraction_max": (
                        max(values)
                        if values
                        else None
                    ),
                    "QC_failed_or_missing": (
                        sum(
                            row["side"] == side
                            and not row[
                                "passes_DP2_GQ20"
                            ]
                            for row in subset
                        )
                    ),
                }
            )

    write_tsv(
        OUTPUT_SUMMARY,
        summary_rows,
        [
            "block",
            "module",
            "side",
            "independent_sites",
            "allele_fraction_observed",
            "allele_fraction_median",
            "allele_fraction_min",
            "allele_fraction_max",
            "QC_failed_or_missing",
        ],
    )

    figure, axis = plt.subplots(
        figsize=(
            9.2,
            3.15,
        )
    )

    axis.axhline(
        0.0,
        color="#C8C8C8",
        linewidth=0.7,
        zorder=0,
    )

    axis.axhline(
        0.5,
        color=HET_LINE_COLOUR,
        linewidth=0.9,
        linestyle=(
            0,
            (
                3,
                3,
            ),
        ),
        zorder=0,
    )

    axis.axhline(
        1.0,
        color="#C8C8C8",
        linewidth=0.7,
        zorder=0,
    )

    for index, layout in enumerate(
        block_layout
    ):
        start = layout["start"]
        end = layout["end"]

        background = (
            "#F4F4F4"
            if index % 2 == 0
            else "#EAEAEA"
        )

        axis.add_patch(
            Rectangle(
                (
                    start,
                    1.035,
                ),
                end - start,
                0.105,
                facecolor=background,
                edgecolor="white",
                linewidth=1.2,
                clip_on=False,
                zorder=0,
            )
        )

        axis.text(
            (
                start + end
            )
            / 2,
            1.086,
            (
                f'{layout["display"]}\n'
                f'(n={layout["sites"]})'
            ),
            ha="center",
            va="center",
            fontsize=8,
            fontweight="bold",
            clip_on=False,
        )

        if index > 0:
            axis.axvline(
                start,
                color="white",
                linewidth=1.5,
                zorder=0,
            )

        if layout["sites"] == 0:
            axis.text(
                (
                    start + end
                )
                / 2,
                0.5,
                (
                    "no BA>=0.80\n"
                    "comparative SNP"
                ),
                ha="center",
                va="center",
                fontsize=7,
                color="#777777",
            )

    rows_by_site = defaultdict(dict)

    for row in evidence_rows:
        site_key = (
            row["block"],
            row["chrom"],
            row["pos"],
            row["ref"],
            row["alt"],
        )

        rows_by_site[
            site_key
        ][
            row["side"]
        ] = row

    x_offset = 0.012

    for site_rows in (
        rows_by_site.values()
    ):
        left = site_rows.get(
            "Left"
        )

        right = site_rows.get(
            "Right"
        )

        if (
            left is None
            or right is None
        ):
            continue

        left_y = left[
            "associated_allele_fraction"
        ]

        right_y = right[
            "associated_allele_fraction"
        ]

        x_value = left[
            "plot_x"
        ]

        if (
            left_y is not None
            and right_y is not None
        ):
            axis.plot(
                [
                    x_value - x_offset,
                    x_value + x_offset,
                ],
                [
                    left_y,
                    right_y,
                ],
                color="#8F8F8F",
                linewidth=0.45,
                alpha=0.45,
                zorder=1,
            )

        point_rows = (
            (
                left,
                LEFT_COLOUR,
                -x_offset,
            ),
            (
                right,
                RIGHT_COLOUR,
                x_offset,
            ),
        )

        for (
            row,
            colour,
            offset,
        ) in point_rows:
            y_value = row[
                "associated_allele_fraction"
            ]

            if y_value is None:
                continue

            if row[
                "passes_DP2_GQ20"
            ]:
                axis.scatter(
                    x_value + offset,
                    y_value,
                    s=12,
                    facecolor=colour,
                    edgecolor="white",
                    linewidth=0.25,
                    alpha=0.82,
                    zorder=3,
                )

            else:
                axis.scatter(
                    x_value + offset,
                    y_value,
                    s=17,
                    facecolor="none",
                    edgecolor=colour,
                    linewidth=0.75,
                    alpha=0.8,
                    zorder=3,
                )

    axis.set_xlim(
        0,
        current_x,
    )

    axis.set_ylim(
        -0.055,
        1.16,
    )

    axis.set_xticks([])

    axis.set_yticks(
        [
            0.0,
            0.5,
            1.0,
        ]
    )

    axis.set_yticklabels(
        [
            (
                "0\nnon-associated "
                "homozygote"
            ),
            (
                "0.5\nheterozygote "
                "/ balanced"
            ),
            (
                "1\nassociated "
                "homozygote"
            ),
        ]
    )

    axis.set_ylabel(
        "Pattern-associated allele fraction",
        fontsize=9,
    )

    axis.tick_params(
        axis="y",
        labelsize=7.5,
        length=3,
    )

    axis.spines[
        "top"
    ].set_visible(False)

    axis.spines[
        "right"
    ].set_visible(False)

    axis.spines[
        "bottom"
    ].set_visible(False)

    legend = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=LEFT_COLOUR,
            markeredgecolor="white",
            markersize=6,
            label="Left (female)",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=RIGHT_COLOUR,
            markeredgecolor="white",
            markersize=6,
            label="Right (male)",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="none",
            markeredgecolor="#666666",
            markersize=6,
            label=(
                "Focal DP/GQ below threshold"
            ),
        ),
    ]

    axis.legend(
        handles=legend,
        loc="lower center",
        bbox_to_anchor=(
            0.5,
            -0.17,
        ),
        ncol=3,
        frameon=False,
        fontsize=8,
        handletextpad=0.4,
        columnspacing=1.2,
    )

    axis.text(
        current_x,
        -0.045,
        (
            "Independent comparative sites; "
            "equal-width blocks and SNP spacing "
            "are not genomic scale"
        ),
        ha="right",
        va="top",
        fontsize=6.8,
        color="#555555",
    )

    figure.savefig(
        OUTPUT_PDF,
        bbox_inches="tight",
        transparent=True,
    )

    figure.savefig(
        OUTPUT_SVG,
        bbox_inches="tight",
        transparent=True,
    )

    plt.close(figure)

    print(
        "COMPACT MAIN-PANEL FOCAL EVIDENCE"
    )

    print()

    print(
        "Site columns were frozen from "
        "comparative-only BA>=0.80 tables."
    )

    print(
        "Focal genotype, AD, DP, GQ, and "
        "expected zygosity were not used to "
        "select or remove sites."
    )

    print()

    print(
        "block\tmodule\tindependent_sites\t"
        "left_observed\tright_observed\t"
        "left_QC_fail\tright_QC_fail"
    )

    for layout in block_layout:
        left = next(
            row
            for row in summary_rows
            if (
                row["block"]
                == layout["block"]
                and row["side"]
                == "Left"
            )
        )

        right = next(
            row
            for row in summary_rows
            if (
                row["block"]
                == layout["block"]
                and row["side"]
                == "Right"
            )
        )

        print(
            f'{layout["block"]}\t'
            f'{layout["module"]}\t'
            f'{layout["sites"]}\t'
            f'{left["allele_fraction_observed"]}\t'
            f'{right["allele_fraction_observed"]}\t'
            f'{left["QC_failed_or_missing"]}\t'
            f'{right["QC_failed_or_missing"]}'
        )

    print()
    print(f"PDF: {OUTPUT_PDF}")
    print(f"SVG: {OUTPUT_SVG}")
    print(
        f"Long evidence table: {OUTPUT_LONG}"
    )
    print(f"Summary: {OUTPUT_SUMMARY}")


if __name__ == "__main__":
    main()
