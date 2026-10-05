#!/usr/bin/env python3

from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd


PROJECT = Path(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
)

MODELED = (
    PROJECT
    / "analysis/original_502/checkpoints/"
    "focal_calls_at_frozen_sites_biological."
    "modeled_LR_posterior95.wide.tsv"
)

MIXTURE = (
    PROJECT
    / "analysis/original_502/reciprocal_L_R_mixture/"
    "optix_reciprocal_L_R_genotype_pairs.posterior95.tsv"
)

COMPARATIVE = (
    PROJECT
    / "data/derived/region_genotypes/"
    "optix_51sample.genotypes.tsv.gz"
)

METADATA = (
    PROJECT
    / "data/derived/region_panel_metadata.tsv"
)

OUTDIR = (
    PROJECT
    / "analysis/figures/modeled_LR_supplements"
)

AUDIT = (
    PROJECT
    / "analysis/original_502/checkpoints/"
    "modeled_LR_supplementary_figure.audit.tsv"
)

FROZEN_PDF = (
    OUTDIR
    / "optix_frozen_sites_modeled_LR_genotypes_and_read_fractions.pdf"
)

HEATMAP_PDF = (
    OUTDIR
    / "optix_comparative_plus_modeled_LR_Band_genotypes.pdf"
)

MODULES = [
    "Band",
    "Ray1",
    "Ray2",
    "Dennis1",
    "Dennis2",
    "Yb",
]

DISPLAY = {
    "Band": "Band",
    "Ray1": "Ray 1",
    "Ray2": "Ray 2",
    "Dennis1": "Dennis 1",
    "Dennis2": "Dennis 2",
    "Yb": "Yb",
}

RANGES = {
    "Band": (1112651, 1114314),
    "Ray1": (1156746, 1160220),
    "Ray2": (1174547, 1187833),
    "Dennis1": (1215210, 1216124),
    "Dennis2": (1218557, 1221494),
    "Yb": (1627000, 1637000),
}

BAND_START = 1112000
BAND_END = 1116000

OLD_FOCAL = {
    "TN2L_S1",
    "TN3L_S3",
    "TN2R_S2",
    "TN3R_S4",
}

GT_COLORS = {
    0: "#D9D9D9",
    1: "#E69F00",
    2: "#CC79A7",
    3: "#0072B2",
    4: "#000000",
}

BIO_COLORS = {
    "missing": "#D9D9D9",
    "AA": "#E69F00",
    "PA": "#CC79A7",
    "PP": "#0072B2",
}

LEFT_COLOR = "#0072B2"
RIGHT_COLOR = "#D55E00"
EXACT_COLOR = "#009E73"
OTHER_COLOR = "#A8A8A8"


def as_bool(series):
    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .isin({"true", "t", "1", "yes"})
    )


def numeric(series):
    return pd.to_numeric(
        series.replace(
            {
                "": np.nan,
                ".": np.nan,
                "NA": np.nan,
            }
        ),
        errors="coerce",
    )


def canonical_gt(value):
    text = str(value).strip().replace("|", "/")

    if text == "1/0":
        text = "0/1"

    if text in {"0/0", "0/1", "1/1"}:
        return text

    return None


def genotype_state(value):
    genotype = canonical_gt(value)

    if genotype is not None:
        return {
            "0/0": 1,
            "0/1": 2,
            "1/1": 3,
        }[genotype]

    if str(value).strip() in {
        "",
        ".",
        "./.",
        ".|.",
        "nan",
        "NA",
    }:
        return 0

    return 4


def site_key(frame):
    return (
        frame["chrom"].astype(str)
        + ":"
        + frame["pos"].astype(str)
        + ":"
        + frame["ref"].astype(str)
        + ">"
        + frame["alt"].astype(str)
    )


def resolve_column(frame, requested):
    target = re.sub(
        r"[^a-z0-9]+",
        "",
        requested.lower(),
    )

    matches = [
        column
        for column in frame.columns
        if re.sub(
            r"[^a-z0-9]+",
            "",
            column.lower(),
        )
        == target
    ]

    if len(matches) != 1:
        raise RuntimeError(
            f"Could not resolve metadata column "
            f"{requested!r}; matches={matches}"
        )

    return matches[0]


def load_modeled():
    data = pd.read_csv(
        MODELED,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    required = {
        "site_id",
        "module",
        "chrom",
        "pos",
        "ref",
        "alt",
        "left_biological_genotype",
        "right_biological_genotype",
        "left_presence_allele_fraction",
        "right_presence_allele_fraction",
        "left_expected_genotype",
        "right_expected_genotype",
        "both_callable",
        "pair_matches_expectation",
        "mixture_pair_posterior",
        "focal_call_source",
    }

    missing = required - set(data.columns)

    if missing:
        raise RuntimeError(
            "Modeled table missing columns: "
            f"{sorted(missing)}"
        )

    if len(data) != 117:
        raise RuntimeError(
            f"Expected 117 frozen sites; found {len(data)}"
        )

    if data["site_id"].duplicated().any():
        raise RuntimeError(
            "The modeled table contains duplicate site_id values"
        )

    forbidden = {
        "left_original_GT",
        "right_original_GT",
        "original_pair_GT",
        "left_QC_pass_DP2_GQ20",
        "right_QC_pass_DP2_GQ20",
    }

    present_forbidden = forbidden & set(data.columns)

    if present_forbidden:
        raise RuntimeError(
            "Old focal-GT or focal DP/GQ columns entered "
            f"the modeled table: {sorted(present_forbidden)}"
        )

    data["position"] = numeric(data["pos"])
    data["posterior"] = numeric(
        data["mixture_pair_posterior"]
    )
    data["callable"] = as_bool(
        data["both_callable"]
    )
    data["expected"] = as_bool(
        data["pair_matches_expectation"]
    )
    data["left_fraction"] = numeric(
        data["left_presence_allele_fraction"]
    )
    data["right_fraction"] = numeric(
        data["right_presence_allele_fraction"]
    )

    called = data["callable"]

    if not data.loc[
        called,
        "posterior",
    ].ge(0.95).all():
        raise RuntimeError(
            "A callable focal pair has posterior below 0.95"
        )

    if not data.loc[
        called,
        "focal_call_source",
    ].eq(
        "reciprocal_mixture_posterior95_only"
    ).all():
        raise RuntimeError(
            "A callable row does not use the "
            "reciprocal-mixture model"
        )

    for column in [
        "left_biological_genotype",
        "right_biological_genotype",
    ]:
        valid = data.loc[
            called,
            column,
        ].isin({"AA", "PA", "PP"})

        if not valid.all():
            raise RuntimeError(
                f"Invalid biological genotype in {column}"
            )

        data.loc[
            ~called,
            column,
        ] = "missing"

    data.loc[
        ~called,
        [
            "left_fraction",
            "right_fraction",
        ],
    ] = np.nan

    return data


def frozen_genotype_page(data):
    figure, axes = plt.subplots(
        3,
        2,
        figsize=(11.5, 8.5),
    )

    for axis, module in zip(
        axes.flat,
        MODULES,
    ):
        subset = (
            data.loc[data["module"].eq(module)]
            .sort_values("position")
        )

        for row in subset.itertuples(
            index=False
        ):
            position_mb = row.position / 1e6

            if row.callable:
                connector_color = (
                    EXACT_COLOR
                    if row.expected
                    else OTHER_COLOR
                )

                axis.plot(
                    [position_mb, position_mb],
                    [0, 1],
                    color=connector_color,
                    linewidth=0.8,
                    zorder=1,
                )

            axis.scatter(
                [position_mb],
                [1],
                marker="s",
                s=32,
                color=BIO_COLORS[
                    row.left_biological_genotype
                ],
                edgecolor="white",
                linewidth=0.3,
                zorder=2,
            )

            axis.scatter(
                [position_mb],
                [0],
                marker="s",
                s=32,
                color=BIO_COLORS[
                    row.right_biological_genotype
                ],
                edgecolor="white",
                linewidth=0.3,
                zorder=2,
            )

        minimum, maximum = RANGES[module]

        axis.set_xlim(
            minimum / 1e6,
            maximum / 1e6,
        )
        axis.set_ylim(-0.5, 1.5)
        axis.set_yticks(
            [1, 0],
            ["Left", "Right"],
        )
        axis.set_xlabel(
            "Hmr18 position (Mb)"
        )
        axis.set_title(
            f"{DISPLAY[module]} — "
            f"{int(subset['callable'].sum())}/"
            f"{len(subset)} confident pairs",
            loc="left",
            fontsize=10,
            fontweight="bold",
        )
        axis.grid(
            axis="x",
            color="#EEEEEE",
            linewidth=0.6,
        )
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    handles = [
        Patch(
            facecolor=BIO_COLORS["AA"],
            label="A/A (absence homozygote)",
        ),
        Patch(
            facecolor=BIO_COLORS["PA"],
            label="P/A (heterozygote)",
        ),
        Patch(
            facecolor=BIO_COLORS["PP"],
            label="P/P (presence homozygote)",
        ),
        Patch(
            facecolor=BIO_COLORS["missing"],
            label="No posterior≥0.95 pair call",
        ),
        Line2D(
            [0],
            [0],
            color=EXACT_COLOR,
            label="Exact expected pair",
        ),
        Line2D(
            [0],
            [0],
            color=OTHER_COLOR,
            label="Other confident pair",
        ),
    ]

    figure.legend(
        handles=handles,
        loc="lower center",
        ncol=3,
        frameon=False,
        fontsize=8,
    )

    figure.suptitle(
        "Reciprocal-mixture Left/Right genotypes "
        "at frozen comparative sites",
        x=0.01,
        ha="left",
        fontsize=14,
        fontweight="bold",
    )

    figure.text(
        0.01,
        0.945,
        "Joint pair posterior ≥0.95 only; no original-GT "
        "fallback. Each square is one phenotype-polarized "
        "SNP genotype.",
        ha="left",
        fontsize=9,
    )

    figure.subplots_adjust(
        left=0.09,
        right=0.98,
        top=0.89,
        bottom=0.12,
        hspace=0.55,
        wspace=0.28,
    )

    return figure


def frozen_fraction_page(data):
    figure, axes = plt.subplots(
        3,
        2,
        figsize=(11.5, 8.5),
    )

    expected_dose = {
        "AA": 0.0,
        "PA": 0.5,
        "PP": 1.0,
    }

    for axis, module in zip(
        axes.flat,
        MODULES,
    ):
        all_module = data.loc[
            data["module"].eq(module)
        ]

        subset = (
            all_module.loc[
                all_module["callable"]
            ]
            .sort_values("position")
        )

        background_bands = [
            (0.0, "#F7E7C6"),
            (0.5, "#F3E4EF"),
            (1.0, "#DDECF7"),
        ]

        for center, color in background_bands:
            axis.axhspan(
                max(-0.02, center - 0.055),
                min(1.02, center + 0.055),
                color=color,
                alpha=0.75,
                linewidth=0,
            )

        for row in subset.itertuples(
            index=False
        ):
            position_mb = row.position / 1e6

            axis.plot(
                [position_mb, position_mb],
                [
                    row.left_fraction,
                    row.right_fraction,
                ],
                color="#B8B8B8",
                linewidth=0.7,
                zorder=1,
            )

        axis.scatter(
            subset["position"] / 1e6,
            subset["left_fraction"],
            marker="o",
            s=27,
            color=LEFT_COLOR,
            edgecolor="white",
            linewidth=0.3,
            zorder=3,
        )

        axis.scatter(
            subset["position"] / 1e6,
            subset["right_fraction"],
            marker="^",
            s=30,
            color=RIGHT_COLOR,
            edgecolor="white",
            linewidth=0.3,
            zorder=3,
        )

        if len(all_module):
            left_expected = (
                all_module[
                    "left_expected_genotype"
                ].iloc[0]
            )
            right_expected = (
                all_module[
                    "right_expected_genotype"
                ].iloc[0]
            )

            axis.axhline(
                expected_dose[left_expected],
                color=LEFT_COLOR,
                linestyle="--",
                linewidth=0.8,
            )

            axis.axhline(
                expected_dose[right_expected],
                color=RIGHT_COLOR,
                linestyle=":",
                linewidth=0.9,
            )

        minimum, maximum = RANGES[module]

        axis.set_xlim(
            minimum / 1e6,
            maximum / 1e6,
        )
        axis.set_ylim(-0.03, 1.03)
        axis.set_yticks(
            [0, 0.5, 1],
            ["A/A", "P/A", "P/P"],
        )
        axis.set_xlabel(
            "Hmr18 position (Mb)"
        )
        axis.set_title(
            f"{DISPLAY[module]} — "
            f"{len(subset)} displayed SNP pairs",
            loc="left",
            fontsize=10,
            fontweight="bold",
        )
        axis.grid(
            axis="x",
            color="#EEEEEE",
            linewidth=0.6,
        )
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=LEFT_COLOR,
            markersize=6,
            label="Left observed P-allele read fraction",
        ),
        Line2D(
            [0],
            [0],
            marker="^",
            color="none",
            markerfacecolor=RIGHT_COLOR,
            markersize=6,
            label="Right observed P-allele read fraction",
        ),
        Line2D(
            [0],
            [0],
            color="#B8B8B8",
            label="Same SNP",
        ),
        Line2D(
            [0],
            [0],
            color=LEFT_COLOR,
            linestyle="--",
            label="Expected Left fraction",
        ),
        Line2D(
            [0],
            [0],
            color=RIGHT_COLOR,
            linestyle=":",
            label="Expected Right fraction",
        ),
    ]

    figure.legend(
        handles=handles,
        loc="lower center",
        ncol=3,
        frameon=False,
        fontsize=8,
    )

    figure.suptitle(
        "Observed P-allele read fractions for "
        "posterior≥0.95 modeled pairs",
        x=0.01,
        ha="left",
        fontsize=14,
        fontweight="bold",
    )

    figure.text(
        0.01,
        0.945,
        "Each gray segment is one frozen SNP: blue circle "
        "= Left; orange triangle = Right. Pairs without a "
        "posterior≥0.95 model call are omitted.",
        ha="left",
        fontsize=9,
    )

    figure.subplots_adjust(
        left=0.09,
        right=0.98,
        top=0.89,
        bottom=0.12,
        hspace=0.55,
        wspace=0.28,
    )

    return figure


def load_metadata():
    raw = pd.read_csv(
        METADATA,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    requested = {
        "sample": "Accession",
        "taxon": "Taxon Name",
        "country": "Country",
        "B": "B",
        "R": "R",
        "D": "D",
    }

    columns = {
        name: resolve_column(raw, value)
        for name, value in requested.items()
    }

    data = pd.DataFrame(
        {
            name: raw[column]
            for name, column in columns.items()
        }
    )

    for trait in ["B", "R", "D"]:
        data[trait] = numeric(data[trait])

    return data.drop_duplicates(
        "sample",
        keep="first",
    )


def prepare_heatmap():
    metadata = load_metadata()

    calls = pd.read_csv(
        COMPARATIVE,
        sep="\t",
        dtype=str,
        compression="gzip",
        keep_default_na=False,
    )

    required = {
        "sample",
        "chrom",
        "pos",
        "ref",
        "alt",
        "GT",
        "DP",
        "GQ",
    }

    missing = required - set(calls.columns)

    if missing:
        raise RuntimeError(
            "Comparative genotype columns missing: "
            f"{sorted(missing)}"
        )

    calls["position"] = numeric(calls["pos"])

    calls = calls.loc[
        calls["position"].between(
            BAND_START,
            BAND_END,
        )
    ].copy()

    removed = sorted(
        set(calls["sample"]) & OLD_FOCAL
    )

    calls = calls.loc[
        ~calls["sample"].isin(OLD_FOCAL)
    ].copy()

    duplicate_columns = [
        "sample",
        "chrom",
        "pos",
        "ref",
        "alt",
    ]

    if calls.duplicated(
        duplicate_columns
    ).any():
        raise RuntimeError(
            "Duplicate comparative sample/site rows"
        )

    calls["key"] = site_key(calls)
    calls["state"] = calls["GT"].map(
        genotype_state
    )

    depth = numeric(calls["DP"])
    quality = numeric(calls["GQ"])

    fails_qc = (
        depth.lt(2)
        | quality.lt(20)
        | depth.isna()
        | quality.isna()
    )

    calls.loc[
        fails_qc,
        "state",
    ] = 0

    variants = (
        calls[
            [
                "key",
                "chrom",
                "pos",
                "position",
                "ref",
                "alt",
            ]
        ]
        .drop_duplicates("key")
        .sort_values(
            [
                "position",
                "ref",
                "alt",
            ]
        )
        .reset_index(drop=True)
    )

    samples = sorted(
        calls["sample"].unique()
    )

    matrix = (
        calls.pivot(
            index="sample",
            columns="key",
            values="state",
        )
        .reindex(
            index=samples,
            columns=variants["key"],
        )
        .fillna(0)
        .astype(int)
    )

    mixture = pd.read_csv(
        MIXTURE,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    mixture["position"] = numeric(
        mixture["pos"]
    )

    mixture = mixture.loc[
        mixture["position"].between(
            BAND_START,
            BAND_END,
        )
    ].copy()

    if mixture.duplicated(
        [
            "chrom",
            "pos",
            "ref",
            "alt",
        ]
    ).any():
        raise RuntimeError(
            "Duplicate mixture-model sites"
        )

    mixture["key"] = site_key(mixture)

    confident = (
        as_bool(mixture["posterior95"])
        & numeric(
            mixture["pair_posterior"]
        ).ge(0.95)
    )

    mixture["left_state"] = 0
    mixture["right_state"] = 0

    mixture.loc[
        confident,
        "left_state",
    ] = (
        mixture.loc[
            confident,
            "mixture_left_GT",
        ]
        .map(genotype_state)
    )

    mixture.loc[
        confident,
        "right_state",
    ] = (
        mixture.loc[
            confident,
            "mixture_right_GT",
        ]
        .map(genotype_state)
    )

    lookup = mixture.set_index("key")

    matrix.loc["TN_Left"] = (
        variants["key"]
        .map(lookup["left_state"])
        .fillna(0)
        .astype(int)
        .to_numpy()
    )

    matrix.loc["TN_Right"] = (
        variants["key"]
        .map(lookup["right_state"])
        .fillna(0)
        .astype(int)
        .to_numpy()
    )

    order_table = metadata.loc[
        metadata["sample"].isin(samples)
    ].copy()

    missing_samples = sorted(
        set(samples)
        - set(order_table["sample"])
    )

    if missing_samples:
        unmatched = pd.DataFrame(
            {
                "sample": missing_samples,
                "taxon": "metadata unmatched",
                "country": "",
                "B": np.nan,
                "R": np.nan,
                "D": np.nan,
            }
        )

        order_table = pd.concat(
            [
                order_table,
                unmatched,
            ],
            ignore_index=True,
        )

    taxonomy_order = (
        order_table.sort_values(
            [
                "taxon",
                "country",
                "sample",
            ]
        )["sample"]
        .tolist()
    )

    phenotype_table = order_table.copy()

    phenotype_table[
        [
            "B_sort",
            "R_sort",
            "D_sort",
        ]
    ] = (
        phenotype_table[
            [
                "B",
                "R",
                "D",
            ]
        ]
        .fillna(-1)
    )

    phenotype_order = (
        phenotype_table.sort_values(
            [
                "B_sort",
                "R_sort",
                "D_sort",
                "taxon",
                "country",
                "sample",
            ],
            ascending=[
                False,
                False,
                False,
                True,
                True,
                True,
            ],
        )["sample"]
        .tolist()
    )

    taxonomy_order = [
        "TN_Left",
        "TN_Right",
    ] + taxonomy_order

    phenotype_order = [
        "TN_Left",
        "TN_Right",
    ] + phenotype_order

    return (
        matrix,
        variants,
        metadata.set_index("sample"),
        taxonomy_order,
        phenotype_order,
        removed,
        int(confident.sum()),
    )


def sample_label(sample, metadata):
    if sample == "TN_Left":
        return "Focal Left (mixture model)"

    if sample == "TN_Right":
        return "Focal Right (mixture model)"

    if sample not in metadata.index:
        return (
            f"{sample} — metadata unmatched"
        )

    row = metadata.loc[sample]

    context = ", ".join(
        value
        for value in [
            str(row["taxon"]),
            str(row["country"]),
        ]
        if value
    )

    return f"{sample} — {context}"


def phenotype_values(order, metadata):
    rows = []

    for sample in order:
        if sample == "TN_Left":
            rows.append([0, 1, 1])
        elif sample == "TN_Right":
            rows.append([1, 1, 1])
        elif sample in metadata.index:
            rows.append(
                [
                    metadata.loc[
                        sample,
                        trait,
                    ]
                    for trait in [
                        "B",
                        "R",
                        "D",
                    ]
                ]
            )
        else:
            rows.append(
                [
                    np.nan,
                    np.nan,
                    np.nan,
                ]
            )

    return np.asarray(
        rows,
        dtype=float,
    )


def heatmap_page(
    matrix,
    variants,
    metadata,
    order,
    order_name,
):
    height = max(
        10,
        0.225 * len(order) + 2.5,
    )

    figure = plt.figure(
        figsize=(14.2, height)
    )

    grid = figure.add_gridspec(
        1,
        2,
        width_ratios=[
            0.75,
            11,
        ],
        left=0.24,
        right=0.98,
        top=0.92,
        bottom=0.08,
        wspace=0.03,
    )

    phenotype_axis = figure.add_subplot(
        grid[0, 0]
    )
    genotype_axis = figure.add_subplot(
        grid[0, 1]
    )

    phenotype = np.ma.masked_invalid(
        phenotype_values(
            order,
            metadata,
        )
    )

    phenotype_cmap = ListedColormap(
        [
            "#FFFFFF",
            "#333333",
        ]
    )
    phenotype_cmap.set_bad("#D9D9D9")

    phenotype_axis.imshow(
        phenotype,
        aspect="auto",
        interpolation="nearest",
        cmap=phenotype_cmap,
        vmin=0,
        vmax=1,
    )

    phenotype_axis.set_xticks(
        [0, 1, 2],
        ["B", "R", "D"],
    )

    phenotype_axis.tick_params(
        axis="x",
        top=True,
        labeltop=True,
        bottom=False,
        labelbottom=False,
        length=0,
    )

    phenotype_axis.set_yticks([])

    values = (
        matrix.reindex(
            index=order,
            columns=variants["key"],
        )
        .fillna(0)
        .to_numpy()
    )

    genotype_axis.imshow(
        values,
        aspect="auto",
        interpolation="nearest",
        cmap=ListedColormap(
            [
                GT_COLORS[index]
                for index in range(5)
            ]
        ),
        vmin=-0.5,
        vmax=4.5,
    )

    genotype_axis.set_yticks(
        range(len(order)),
        [
            sample_label(
                sample,
                metadata,
            )
            for sample in order
        ],
        fontsize=6.5,
    )

    tick_indices = np.unique(
        np.linspace(
            0,
            len(variants) - 1,
            6,
        )
        .round()
        .astype(int)
    )

    genotype_axis.set_xticks(
        tick_indices,
        [
            (
                f"{variants.iloc[index]['position'] / 1e6:.4f}"
            )
            for index in tick_indices
        ],
    )

    genotype_axis.set_xlabel(
        "Hmr18 position (Mb); columns are "
        "comparative Band-region variants"
    )

    phenotype_axis.axhline(
        1.5,
        color="black",
        linewidth=1.3,
    )

    genotype_axis.axhline(
        1.5,
        color="black",
        linewidth=1.3,
    )

    handles = [
        Patch(
            facecolor=GT_COLORS[state],
            edgecolor="#555555",
            label=label,
        )
        for state, label in [
            (0, "Unavailable"),
            (1, "Homozygous REF"),
            (2, "Heterozygous"),
            (3, "Homozygous ALT"),
            (4, "Other/multiallelic"),
        ]
    ]

    figure.legend(
        handles=handles,
        loc="lower center",
        ncol=5,
        frameon=False,
        fontsize=8,
    )

    figure.suptitle(
        f"Optix Band-region genotypes — "
        f"{order_name}",
        x=0.01,
        ha="left",
        fontsize=14,
        fontweight="bold",
    )

    figure.text(
        0.01,
        0.955,
        "Comparative rows: DP≥2/GQ≥20 hard GT. "
        "Focal rows: reciprocal-mixture posterior≥0.95 "
        "only; no original-GT fallback. "
        "T2/T3 samples are excluded.",
        ha="left",
        fontsize=9,
    )

    return figure


def main():
    input_paths = [
        MODELED,
        MIXTURE,
        COMPARATIVE,
        METADATA,
    ]

    for path in input_paths:
        if not path.exists():
            raise FileNotFoundError(path)

    OUTDIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    modeled = load_modeled()

    with PdfPages(FROZEN_PDF) as pdf:
        for page_maker in [
            frozen_genotype_page,
            frozen_fraction_page,
        ]:
            figure = page_maker(modeled)
            pdf.savefig(figure)
            plt.close(figure)

    (
        matrix,
        variants,
        metadata,
        taxonomy_order,
        phenotype_order,
        removed,
        band_confident,
    ) = prepare_heatmap()

    displayed_samples = set(
        taxonomy_order + phenotype_order
    )

    forbidden_displayed = (
        displayed_samples & OLD_FOCAL
    )

    if forbidden_displayed:
        raise RuntimeError(
            "T2/T3 focal samples reached the heatmap: "
            f"{sorted(forbidden_displayed)}"
        )

    with PdfPages(HEATMAP_PDF) as pdf:
        page_definitions = [
            (
                taxonomy_order,
                "taxonomy order",
            ),
            (
                phenotype_order,
                "B/R/D phenotype order",
            ),
        ]

        for order, label in page_definitions:
            figure = heatmap_page(
                matrix,
                variants,
                metadata,
                order,
                label,
            )
            pdf.savefig(figure)
            plt.close(figure)

    audit_rows = []

    for module in MODULES:
        subset = modeled.loc[
            modeled["module"].eq(module)
        ]

        audit_rows.append(
            {
                "section": "frozen_sites",
                "group": module,
                "rows": len(subset),
                "posterior95_pairs": int(
                    subset["callable"].sum()
                ),
                "exact_expected_pairs": int(
                    subset["expected"].sum()
                ),
            }
        )

    audit_rows.append(
        {
            "section": "Band_heatmap",
            "group": "all",
            "rows": len(taxonomy_order),
            "posterior95_pairs": (
                band_confident
            ),
            "exact_expected_pairs": np.nan,
        }
    )

    pd.DataFrame(
        audit_rows
    ).to_csv(
        AUDIT,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    expected_outputs = [
        FROZEN_PDF,
        HEATMAP_PDF,
        AUDIT,
    ]

    for path in expected_outputs:
        if (
            not path.exists()
            or path.stat().st_size == 0
        ):
            raise RuntimeError(
                "Expected output is missing or empty: "
                f"{path}"
            )

    print(
        "MODELED LEFT/RIGHT SUPPLEMENTARY "
        "FIGURES COMPLETE"
    )
    print(
        f"Frozen-site PDF: {FROZEN_PDF}"
    )
    print(
        f"Band heatmap PDF: {HEATMAP_PDF}"
    )
    print(
        f"Audit: {AUDIT}"
    )
    print(
        f"Frozen sites: {len(modeled)}"
    )
    print(
        "Frozen posterior95 pairs: "
        f"{int(modeled['callable'].sum())}"
    )
    print(
        f"Band variants: {len(variants)}"
    )
    print(
        "Displayed heatmap rows: "
        f"{len(taxonomy_order)}"
    )
    print(
        "Old T2/T3 focal rows removed: "
        + (
            ", ".join(removed)
            if removed
            else "none present"
        )
    )
    print(
        "Displayed focal rows: "
        "TN_Left and TN_Right"
    )
    print(
        "Original focal GT fallback: disabled"
    )


if __name__ == "__main__":
    main()
