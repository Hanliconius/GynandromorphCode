#!/usr/bin/env python3

from pathlib import Path
import re

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")
GENOTYPE_FILE = (
    PROJECT
    / "data/derived/region_genotypes/optix_51sample.genotypes.tsv.gz"
)
METADATA_FILE = (
    PROJECT / "data/derived/region_panel_metadata.tsv"
)
FIGURE_DIR = PROJECT / "analysis/figures"
CHECKPOINT_DIR = PROJECT / "analysis/checkpoints"

BAND_START = 1_112_000
BAND_END = 1_116_000

FOCAL_ORDER = [
    "TN2L_S1",
    "TN3L_S3",
    "TN2R_S2",
    "TN3R_S4",
]

FOCAL_LABELS = {
    "TN2L_S1": "TN2L - left, T2",
    "TN3L_S3": "TN3L - left, T3",
    "TN2R_S2": "TN2R - right, T2",
    "TN3R_S4": "TN3R - right, T3",
}

STATE_COLORS = [
    "#D9D9D9",  # missing
    "#E69F00",  # homozygous REF
    "#CC79A7",  # heterozygous
    "#0072B2",  # homozygous ALT
    "#000000",  # other/multiallelic
]

STATE_LABELS = {
    0: "Missing",
    1: "Homozygous REF",
    2: "Heterozygous",
    3: "Homozygous ALT",
    4: "Other/multiallelic",
}


def normalize_name(value):
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def resolve_column(dataframe, requested):
    target = normalize_name(requested)
    matches = [
        column for column in dataframe.columns
        if normalize_name(column) == target
    ]

    if len(matches) != 1:
        raise ValueError(
            f"Could not resolve {requested!r}; matches={matches}"
        )

    return matches[0]


def numeric(series):
    return pd.to_numeric(
        series.replace({"": np.nan, ".": np.nan}),
        errors="coerce",
    )


def genotype_state(value):
    value = str(value)

    if value in {"", ".", "./.", ".|.", "nan"}:
        return 0

    fields = re.split(r"[/|]", value)

    if len(fields) != 2 or any(field == "." for field in fields):
        return 0

    try:
        alleles = [int(field) for field in fields]
    except ValueError:
        return 4

    if alleles == [0, 0]:
        return 1

    if alleles[0] != alleles[1]:
        return 2

    if alleles[0] == alleles[1] and alleles[0] > 0:
        return 3

    return 4


def shorten(text, maximum=42):
    text = str(text)

    if len(text) <= maximum:
        return text

    return text[:maximum - 3] + "..."


def make_sample_label(sample, metadata_by_sample):
    if sample in FOCAL_LABELS:
        return FOCAL_LABELS[sample]

    if sample not in metadata_by_sample.index:
        return f"{sample} - metadata unmatched"

    taxon = metadata_by_sample.loc[sample, "taxon_qc"]

    if not taxon:
        return sample

    return f"{sample} - {shorten(taxon)}"


def build_order(samples, metadata, ordering):
    focal = [
        sample for sample in FOCAL_ORDER
        if sample in samples
    ]

    remaining = metadata.loc[
        metadata["sample"].isin(set(samples) - set(focal))
    ].copy()

    if ordering == "taxonomy":
        remaining = remaining.sort_values(
            ["taxon_qc", "country_qc", "sample"],
            na_position="last",
        )
    else:
        remaining["B_sort"] = remaining["B_qc"].fillna(-1)
        remaining["R_sort"] = remaining["R_qc"].fillna(-1)
        remaining["D_sort"] = remaining["D_qc"].fillna(-1)

        remaining = remaining.sort_values(
            [
                "B_sort",
                "R_sort",
                "D_sort",
                "taxon_qc",
                "sample",
            ],
            ascending=[False, False, False, True, True],
            na_position="last",
        )

    ordered = focal + remaining["sample"].tolist()

    unmatched = sorted(set(samples) - set(ordered))
    ordered.extend(unmatched)

    return ordered


def phenotype_matrix(order, metadata_by_sample):
    matrix = np.full((len(order), 3), np.nan)

    for row_index, sample in enumerate(order):
        if sample not in metadata_by_sample.index:
            continue

        for column_index, trait in enumerate(
            ["B_qc", "R_qc", "D_qc"]
        ):
            value = metadata_by_sample.loc[sample, trait]

            if pd.notna(value):
                matrix[row_index, column_index] = value

    return matrix


def genotype_matrix(calls, order, variants, state_column):
    pivot = calls.pivot_table(
        index="sample",
        columns="variant_key",
        values=state_column,
        aggfunc="first",
    )

    return (
        pivot
        .reindex(index=order, columns=variants["variant_key"])
        .fillna(0)
        .to_numpy(dtype=float)
    )


def plot_heatmap(
    calls,
    variants,
    metadata_by_sample,
    order,
    ordering,
    state_column,
    state_description,
    output_png,
):
    genotype = genotype_matrix(
        calls,
        order,
        variants,
        state_column,
    )
    phenotypes = phenotype_matrix(order, metadata_by_sample)

    figure_height = max(10, 0.24 * len(order) + 3)

    figure = plt.figure(
        figsize=(17, figure_height),
        constrained_layout=False,
    )

    grid = figure.add_gridspec(
        nrows=1,
        ncols=2,
        width_ratios=[0.55, 6.0],
        left=0.20,
        right=0.98,
        bottom=0.13,
        top=0.90,
        wspace=0.04,
    )

    phenotype_axis = figure.add_subplot(grid[0, 0])
    genotype_axis = figure.add_subplot(grid[0, 1])

    phenotype_cmap = ListedColormap(
        ["#FFFFFF", "#333333"]
    )
    phenotype_cmap.set_bad("#D9D9D9")

    phenotype_axis.imshow(
        np.ma.masked_invalid(phenotypes),
        aspect="auto",
        interpolation="nearest",
        cmap=phenotype_cmap,
        vmin=0,
        vmax=1,
    )

    phenotype_axis.set_xticks([0, 1, 2])
    phenotype_axis.set_xticklabels(
        ["B", "R", "D"],
        fontsize=9,
        fontweight="bold",
    )
    phenotype_axis.xaxis.tick_top()
    phenotype_axis.set_yticks(range(len(order)))
    phenotype_axis.set_yticklabels(
        [
            make_sample_label(sample, metadata_by_sample)
            for sample in order
        ],
        fontsize=6.5,
    )
    phenotype_axis.tick_params(axis="y", length=0)
    phenotype_axis.tick_params(axis="x", length=0)

    genotype_axis.imshow(
        genotype,
        aspect="auto",
        interpolation="nearest",
        cmap=ListedColormap(STATE_COLORS),
        vmin=0,
        vmax=4,
    )

    genotype_axis.set_yticks([])
    genotype_axis.set_xlabel(
        "Variant position on Hmr18",
        fontsize=10,
    )

    tick_count = min(14, len(variants))
    tick_indices = np.unique(
        np.linspace(
            0,
            len(variants) - 1,
            tick_count,
            dtype=int,
        )
    )

    genotype_axis.set_xticks(tick_indices)
    genotype_axis.set_xticklabels(
        [
            f"{int(variants.iloc[index]['pos_numeric']):,}"
            for index in tick_indices
        ],
        rotation=45,
        ha="right",
        fontsize=8,
    )

    for axis in [phenotype_axis, genotype_axis]:
        for spine in axis.spines.values():
            spine.set_linewidth(0.6)
            spine.set_color("#555555")

    focal_count = sum(
        sample in FOCAL_ORDER for sample in order
    )

    if focal_count and focal_count < len(order):
        boundary = focal_count - 0.5
        phenotype_axis.axhline(
            boundary,
            color="black",
            linewidth=1.4,
        )
        genotype_axis.axhline(
            boundary,
            color="black",
            linewidth=1.4,
        )

    legend_handles = [
        Patch(
            facecolor=STATE_COLORS[state],
            edgecolor="#555555",
            linewidth=0.5,
            label=STATE_LABELS[state],
        )
        for state in range(5)
    ]

    figure.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.64, 0.965),
        ncol=5,
        frameon=False,
        fontsize=9,
    )

    figure.suptitle(
        (
            "Optix Band-region genotype reconstruction\n"
            f"Hmr18:{BAND_START:,}-{BAND_END:,}; "
            f"{state_description}; samples ordered by {ordering}"
        ),
        x=0.60,
        y=0.995,
        fontsize=13,
    )

    figure.text(
        0.012,
        0.96,
        "Phenotype",
        fontsize=9,
        fontweight="bold",
    )
    figure.text(
        0.012,
        0.935,
        "white = absent; black = present; gray = unavailable",
        fontsize=7,
    )

    figure.savefig(
        output_png,
        dpi=220,
        bbox_inches="tight",
    )

    return figure


def side_agreement_summary(calls, state_column):
    pivot = calls.pivot_table(
        index="variant_key",
        columns="sample",
        values=state_column,
        aggfunc="first",
    )

    for sample in FOCAL_ORDER:
        if sample not in pivot.columns:
            pivot[sample] = 0

    left_called = (
        pivot["TN2L_S1"].ne(0)
        & pivot["TN3L_S3"].ne(0)
    )
    right_called = (
        pivot["TN2R_S2"].ne(0)
        & pivot["TN3R_S4"].ne(0)
    )

    left_agree = (
        left_called
        & (
            pivot["TN2L_S1"]
            == pivot["TN3L_S3"]
        )
    )
    right_agree = (
        right_called
        & (
            pivot["TN2R_S2"]
            == pivot["TN3R_S4"]
        )
    )

    both_sides_agree = left_agree & right_agree
    sides_differ = (
        both_sides_agree
        & (
            pivot["TN2L_S1"]
            != pivot["TN2R_S2"]
        )
    )

    return {
        "variants": len(pivot),
        "left_both_called": int(left_called.sum()),
        "left_legs_agree": int(left_agree.sum()),
        "right_both_called": int(right_called.sum()),
        "right_legs_agree": int(right_agree.sum()),
        "both_sides_internally_agree": int(
            both_sides_agree.sum()
        ),
        "internally_consistent_left_right_differences": int(
            sides_differ.sum()
        ),
    }


def main():
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(
        METADATA_FILE,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    accession = resolve_column(metadata, "Accession")
    taxon = resolve_column(metadata, "Taxon Name")
    country = resolve_column(metadata, "Country")
    band = resolve_column(metadata, "B")
    ray = resolve_column(metadata, "R")
    dennis = resolve_column(metadata, "D")

    metadata["sample"] = metadata[accession].astype(str)
    metadata["taxon_qc"] = metadata[taxon].astype(str)
    metadata["country_qc"] = metadata[country].astype(str)
    metadata["B_qc"] = numeric(metadata[band])
    metadata["R_qc"] = numeric(metadata[ray])
    metadata["D_qc"] = numeric(metadata[dennis])

    metadata_by_sample = metadata.set_index("sample")

    calls = pd.read_csv(
        GENOTYPE_FILE,
        sep="\t",
        dtype=str,
        compression="gzip",
        keep_default_na=False,
    )

    calls["pos_numeric"] = numeric(calls["pos"])
    calls["DP_numeric"] = numeric(calls["DP"])
    calls["GQ_numeric"] = numeric(calls["GQ"])

    calls = calls.loc[
        calls["pos_numeric"].between(
            BAND_START,
            BAND_END,
        )
    ].copy()

    calls["variant_key"] = (
        calls["chrom"].astype(str)
        + ":"
        + calls["pos"].astype(str)
        + ":"
        + calls["ref"].astype(str)
        + ">"
        + calls["alt"].astype(str)
    )

    calls["state_raw"] = calls["GT"].map(genotype_state)

    calls["state_qc"] = calls["state_raw"]
    calls.loc[
        calls["DP_numeric"].lt(2)
        | calls["GQ_numeric"].lt(20)
        | calls["DP_numeric"].isna()
        | calls["GQ_numeric"].isna(),
        "state_qc",
    ] = 0

    variants = (
        calls[
            [
                "variant_key",
                "chrom",
                "pos",
                "pos_numeric",
                "ref",
                "alt",
            ]
        ]
        .drop_duplicates()
        .sort_values(
            ["pos_numeric", "ref", "alt"]
        )
        .reset_index(drop=True)
    )

    samples = calls["sample"].drop_duplicates().tolist()

    taxonomy_order = build_order(
        samples,
        metadata,
        "taxonomy",
    )
    phenotype_order = build_order(
        samples,
        metadata,
        "phenotype",
    )

    configurations = [
        (
            "raw",
            "state_raw",
            "hard GT without additional masking",
        ),
        (
            "DP2_GQ20",
            "state_qc",
            "hard GT masked unless DP >= 2 and GQ >= 20",
        ),
    ]

    output_rows = []

    for suffix, state_column, description in configurations:
        pdf_file = (
            FIGURE_DIR
            / f"optix_51sample_band_genotypes_{suffix}.pdf"
        )
        taxonomy_png = (
            FIGURE_DIR
            / f"optix_51sample_band_genotypes_{suffix}_taxonomy.png"
        )
        phenotype_png = (
            FIGURE_DIR
            / f"optix_51sample_band_genotypes_{suffix}_phenotype.png"
        )

        with PdfPages(pdf_file) as pdf:
            taxonomy_figure = plot_heatmap(
                calls,
                variants,
                metadata_by_sample,
                taxonomy_order,
                "taxonomy",
                state_column,
                description,
                taxonomy_png,
            )
            pdf.savefig(
                taxonomy_figure,
                bbox_inches="tight",
            )
            plt.close(taxonomy_figure)

            phenotype_figure = plot_heatmap(
                calls,
                variants,
                metadata_by_sample,
                phenotype_order,
                "B/R/D phenotype",
                state_column,
                description,
                phenotype_png,
            )
            pdf.savefig(
                phenotype_figure,
                bbox_inches="tight",
            )
            plt.close(phenotype_figure)

        summary = side_agreement_summary(
            calls,
            state_column,
        )
        summary["genotype_representation"] = suffix
        summary["pdf"] = str(pdf_file)
        output_rows.append(summary)

    summary_table = pd.DataFrame(output_rows)
    summary_file = (
        CHECKPOINT_DIR
        / "optix_51sample_band_plot_summary.tsv"
    )
    summary_table.to_csv(
        summary_file,
        sep="\t",
        index=False,
    )

    print("OPTIX BAND GENOTYPE-PLOT RECONSTRUCTION")
    print(f"Samples: {len(samples)}")
    print(f"Variants in interval: {len(variants)}")
    print(f"Interval: Hmr18:{BAND_START}-{BAND_END}")

    print("\nFOUR-LEG AGREEMENT SUMMARY")
    print(
        summary_table.drop(columns=["pdf"]).to_csv(
            sep="\t",
            index=False,
        ).rstrip()
    )

    print("\nFigures:")
    for row in output_rows:
        print(row["pdf"])

    print(f"\nSummary: {summary_file}")


if __name__ == "__main__":
    main()
