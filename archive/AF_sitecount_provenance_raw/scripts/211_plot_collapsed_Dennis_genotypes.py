#!/usr/bin/env python3

from pathlib import Path
import json

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd


PROJECT = Path(
    "/gpfs/automountdir/gpfs/scratch/"
    "martinlab/joe/GynandroCalling"
)

SOURCE_DIR = (
    PROJECT
    / "analysis/original_502/checkpoints/"
    "krista_DR_regions_reciprocal_genoplots"
)

MANIFEST = (
    PROJECT
    / "analysis/figures/"
    "key_regions_population_context_genoplots/"
    "exact_BA_comparative_sample_manifest.tsv"
)

LR_TABLE = (
    PROJECT
    / "analysis/original_502/"
    "reciprocal_L_R_mixture/"
    "optix_reciprocal_L_R_genotype_pairs.posterior95.tsv"
)

OUTPUT_DIR = (
    PROJECT
    / "analysis/figures/"
    "collapsed_Dennis_vs_nonDennis_plus_recalled_LR_old_palette"
)

COMBINED_PDF = (
    OUTPUT_DIR
    / "optix_Dennis1_Dennis2_"
    "collapsed_D_vs_nonD_plus_recalled_LR.pdf"
)

MATRIX_OUTPUT = (
    OUTPUT_DIR
    / "collapsed_Dennis_genotype_plot_values.tsv.gz"
)

AUDIT_OUTPUT = (
    OUTPUT_DIR
    / "collapsed_Dennis_genotype_plots.audit.json"
)

MIN_GROUP_CALL_RATE = 0.80

MODULES = {
    "Dennis1": (
        SOURCE_DIR
        / "optix_Dennis1_original_region_"
        "reciprocal_genoplot.long.tsv.gz"
    ),
    "Dennis2": (
        SOURCE_DIR
        / "optix_Dennis2_original_region_"
        "reciprocal_genoplot.long.tsv.gz"
    ),
}

GROUP_ORDER = [
    "Dennis-positive",
    "Dennis-negative",
    "TN_Left",
    "TN_Right",
]

GT_TO_AF = {
    "0/0": 0.0,
    "0|0": 0.0,
    "0/1": 0.5,
    "1/0": 0.5,
    "0|1": 0.5,
    "1|0": 0.5,
    "1/1": 1.0,
    "1|1": 1.0,
}


def require_columns(table, required, label):
    missing = sorted(
        set(required) - set(table.columns)
    )

    if missing:
        raise RuntimeError(
            f"{label}: missing columns: {missing}"
        )


def parse_boolean(value):
    return str(value).strip().lower() in {
        "true",
        "t",
        "1",
        "yes",
    }


def load_manifest():
    table = pd.read_csv(
        MANIFEST,
        sep="\t",
        dtype=str,
    )

    require_columns(
        table,
        [
            "BA_pattern",
            "sample",
            "phenotype",
        ],
        "Dennis manifest",
    )

    table = table.loc[
        table["BA_pattern"].eq("Dennis")
    ].copy()

    table["phenotype"] = pd.to_numeric(
        table["phenotype"],
        errors="coerce",
    )

    table = table.loc[
        table["phenotype"].isin([0, 1])
    ].copy()

    duplicate_conflicts = (
        table.groupby("sample")["phenotype"]
        .nunique()
        .gt(1)
    )

    if duplicate_conflicts.any():
        bad = duplicate_conflicts[
            duplicate_conflicts
        ].index.tolist()

        raise RuntimeError(
            "Samples assigned conflicting Dennis "
            f"phenotypes: {bad[:10]}"
        )

    table = (
        table.drop_duplicates(
            subset=["sample"],
            keep="first",
        )
        .reset_index(drop=True)
    )

    table["collapsed_group"] = np.where(
        table["phenotype"].eq(1),
        "Dennis-positive",
        "Dennis-negative",
    )

    counts = (
        table.groupby("collapsed_group")[
            "sample"
        ]
        .nunique()
        .to_dict()
    )

    for group in (
        "Dennis-positive",
        "Dennis-negative",
    ):
        if counts.get(group, 0) == 0:
            raise RuntimeError(
                f"No samples assigned to {group}"
            )

    return table, counts


def load_lr_calls():
    table = pd.read_csv(
        LR_TABLE,
        sep="\t",
        low_memory=False,
    )

    require_columns(
        table,
        [
            "chrom",
            "pos",
            "ref",
            "alt",
            "mixture_left_GT",
            "mixture_right_GT",
            "pair_posterior",
            "posterior95",
        ],
        "Reciprocal L/R table",
    )

    table["pos"] = pd.to_numeric(
        table["pos"],
        errors="raise",
    ).astype(int)

    table["posterior95"] = (
        table["posterior95"]
        .map(parse_boolean)
    )

    if table.duplicated(
        ["chrom", "pos", "ref", "alt"]
    ).any():
        raise RuntimeError(
            "Reciprocal L/R table contains "
            "duplicate sites"
        )

    return table


def prepare_module(
    module,
    source_path,
    manifest,
    group_sizes,
    lr_calls,
):
    source = pd.read_csv(
        source_path,
        sep="\t",
        compression="gzip",
        low_memory=False,
    )

    require_columns(
        source,
        [
            "panel",
            "chrom",
            "pos",
            "ref",
            "alt",
            "site_label",
            "sample",
            "raw_GT",
        ],
        f"{module} source",
    )

    source = source.loc[
        source["panel"].eq(module)
    ].copy()

    source["pos"] = pd.to_numeric(
        source["pos"],
        errors="raise",
    ).astype(int)

    keys = [
        "chrom",
        "pos",
        "ref",
        "alt",
    ]

    sites = (
        source[
            keys + ["site_label"]
        ]
        .drop_duplicates()
        .sort_values(
            ["chrom", "pos", "ref", "alt"],
            kind="mergesort",
        )
        .reset_index(drop=True)
    )

    if sites.duplicated(keys).any():
        raise RuntimeError(
            f"{module}: site labels are not unique"
        )

    sites["site_index"] = np.arange(
        len(sites),
        dtype=int,
    )

    comparative = source.merge(
        manifest[
            [
                "sample",
                "collapsed_group",
            ]
        ],
        on="sample",
        how="inner",
        validate="many_to_one",
    )

    if comparative.duplicated(
        keys + ["sample"]
    ).any():
        raise RuntimeError(
            f"{module}: duplicate comparative "
            "sample/site rows"
        )

    comparative["allele_fraction"] = (
        comparative["raw_GT"]
        .astype(str)
        .map(GT_TO_AF)
    )

    records = []

    for group in (
        "Dennis-positive",
        "Dennis-negative",
    ):
        subset = comparative.loc[
            comparative["collapsed_group"].eq(
                group
            )
        ]

        summary = (
            subset.groupby(
                keys,
                as_index=False,
                sort=False,
            )["allele_fraction"]
            .agg(
                mean_AF="mean",
                called_samples="count",
            )
        )

        summary = sites.merge(
            summary,
            on=keys,
            how="left",
            validate="one_to_one",
        )

        summary["called_samples"] = (
            summary["called_samples"]
            .fillna(0)
            .astype(int)
        )

        summary["total_samples"] = int(
            group_sizes[group]
        )

        summary["call_rate"] = (
            summary["called_samples"]
            / summary["total_samples"]
        )

        summary.loc[
            summary["call_rate"]
            < MIN_GROUP_CALL_RATE,
            "mean_AF",
        ] = np.nan

        summary["module"] = module
        summary["plot_group"] = group
        summary["source"] = (
            "equal-individual mean raw comparative GT"
        )
        summary["pair_posterior"] = np.nan

        records.append(summary)

    focal = sites.merge(
        lr_calls[
            keys
            + [
                "mixture_left_GT",
                "mixture_right_GT",
                "pair_posterior",
                "posterior95",
            ]
        ],
        on=keys,
        how="left",
        validate="one_to_one",
    )

    for group, genotype_column in (
        ("TN_Left", "mixture_left_GT"),
        ("TN_Right", "mixture_right_GT"),
    ):
        summary = focal[
            sites.columns.tolist()
            + [
                genotype_column,
                "pair_posterior",
                "posterior95",
            ]
        ].copy()

        summary["mean_AF"] = (
            summary[genotype_column]
            .astype(str)
            .map(GT_TO_AF)
        )

        confident = (
            summary["posterior95"].fillna(False)
            & summary["pair_posterior"].ge(0.95)
        )

        summary.loc[
            ~confident,
            "mean_AF",
        ] = np.nan

        summary["called_samples"] = (
            confident.astype(int)
        )
        summary["total_samples"] = 1
        summary["call_rate"] = (
            confident.astype(float)
        )
        summary["module"] = module
        summary["plot_group"] = group
        summary["source"] = (
            "reciprocal-mixture recalled genotype"
        )

        records.append(summary)

    output = pd.concat(
        records,
        ignore_index=True,
        sort=False,
    )

    output["plot_group"] = pd.Categorical(
        output["plot_group"],
        categories=GROUP_ORDER,
        ordered=True,
    )

    output = output.sort_values(
        [
            "plot_group",
            "site_index",
        ],
        kind="mergesort",
    ).reset_index(drop=True)

    return sites, output


def make_colormap():
    colourmap = (
        LinearSegmentedColormap.from_list(
            "genotype_dosage",
            [
                (0.0, "#2C9BCB"),
                (0.5, "#F28E2B"),
                (1.0, "#B61F72"),
            ],
        )
    )

    colourmap.set_bad("#F2F2F2")

    return colourmap


def draw_module(module, sites, values):
    matrix = (
        values.pivot(
            index="plot_group",
            columns="site_index",
            values="mean_AF",
        )
        .reindex(GROUP_ORDER)
        .reindex(
            columns=sites["site_index"],
        )
    )

    figure_width = max(
        11.0,
        min(
            20.0,
            7.5 + len(sites) / 45.0,
        ),
    )

    figure, axis = plt.subplots(
        figsize=(figure_width, 3.9),
        constrained_layout=True,
    )

    image = axis.imshow(
        matrix.to_numpy(dtype=float),
        cmap=make_colormap(),
        vmin=0.0,
        vmax=1.0,
        aspect="auto",
        interpolation="nearest",
    )

    axis.axhline(
        1.5,
        color="#303030",
        linewidth=1.6,
    )

    group_counts = {
        group: int(
            values.loc[
                values["plot_group"].eq(group),
                "total_samples",
            ].iloc[0]
        )
        for group in (
            "Dennis-positive",
            "Dennis-negative",
        )
    }

    axis.set_yticks(
        np.arange(len(GROUP_ORDER))
    )

    axis.set_yticklabels(
        [
            (
                "Dennis-positive mean "
                f"(n={group_counts['Dennis-positive']})"
            ),
            (
                "Dennis-negative mean "
                f"(n={group_counts['Dennis-negative']})"
            ),
            "TN_Left recalled",
            "TN_Right recalled",
        ]
    )

    tick_count = min(
        10,
        len(sites),
    )

    tick_indices = np.unique(
        np.linspace(
            0,
            len(sites) - 1,
            tick_count,
        ).round().astype(int)
    )

    axis.set_xticks(tick_indices)

    axis.set_xticklabels(
        [
            f"{sites.iloc[index]['pos'] / 1_000_000:.4f}"
            for index in tick_indices
        ],
        rotation=45,
        ha="right",
    )

    axis.set_xlabel(
        "Position on Hmr18 (Mb)"
    )

    axis.set_title(
        (
            f"{module}: Dennis-positive versus "
            "Dennis-negative genotype frequency"
        ),
        loc="left",
        fontweight="bold",
        fontsize=15,
    )

    axis.text(
        0.0,
        1.025,
        (
            "Comparative rows are equal-individual "
            "mean ALT-allele dosage; grey cells fail "
            "call-rate/confidence requirements."
        ),
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=9.5,
        color="#555555",
    )

    colourbar = figure.colorbar(
        image,
        ax=axis,
        fraction=0.025,
        pad=0.018,
    )

    colourbar.set_ticks(
        [0.0, 0.5, 1.0]
    )

    colourbar.set_ticklabels(
        [
            "0/0",
            "0/1",
            "1/1",
        ]
    )

    colourbar.set_label(
        "Mean ALT-allele dosage"
    )

    for spine in axis.spines.values():
        spine.set_visible(False)

    return figure


def main():
    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    for required in (
        MANIFEST,
        LR_TABLE,
        *MODULES.values(),
    ):
        if (
            not required.is_file()
            or required.stat().st_size == 0
        ):
            raise RuntimeError(
                f"Missing or empty input: {required}"
            )

    manifest, group_sizes = load_manifest()
    lr_calls = load_lr_calls()

    all_values = []
    audit_modules = {}

    with PdfPages(COMBINED_PDF) as combined:
        for module, source_path in MODULES.items():
            sites, values = prepare_module(
                module,
                source_path,
                manifest,
                group_sizes,
                lr_calls,
            )

            all_values.append(values)

            figure = draw_module(
                module,
                sites,
                values,
            )

            module_pdf = (
                OUTPUT_DIR
                / (
                    f"optix_{module}_collapsed_"
                    "D_vs_nonD_plus_recalled_LR.pdf"
                )
            )

            module_png = module_pdf.with_suffix(
                ".png"
            )

            figure.savefig(
                module_pdf,
                bbox_inches="tight",
            )

            figure.savefig(
                module_png,
                dpi=300,
                bbox_inches="tight",
            )

            combined.savefig(
                figure,
                bbox_inches="tight",
            )

            plt.close(figure)

            audit_modules[module] = {
                "sites": int(len(sites)),
                "comparative_cells_retained": {
                    group: int(
                        values.loc[
                            values[
                                "plot_group"
                            ].eq(group),
                            "mean_AF",
                        ].notna().sum()
                    )
                    for group in (
                        "Dennis-positive",
                        "Dennis-negative",
                    )
                },
                "recalled_LR_cells_retained": {
                    group: int(
                        values.loc[
                            values[
                                "plot_group"
                            ].eq(group),
                            "mean_AF",
                        ].notna().sum()
                    )
                    for group in (
                        "TN_Left",
                        "TN_Right",
                    )
                },
                "pdf": str(
                    module_pdf.resolve()
                ),
                "png": str(
                    module_png.resolve()
                ),
            }

    combined_values = pd.concat(
        all_values,
        ignore_index=True,
        sort=False,
    )

    output_columns = [
        "module",
        "chrom",
        "pos",
        "ref",
        "alt",
        "site_label",
        "site_index",
        "plot_group",
        "mean_AF",
        "called_samples",
        "total_samples",
        "call_rate",
        "pair_posterior",
        "source",
    ]

    combined_values[
        output_columns
    ].to_csv(
        MATRIX_OUTPUT,
        sep="\t",
        index=False,
        compression="gzip",
    )

    finite_values = combined_values[
        "mean_AF"
    ].dropna()

    if not finite_values.between(
        0.0,
        1.0,
        inclusive="both",
    ).all():
        raise RuntimeError(
            "An allele-frequency value lies "
            "outside 0-1"
        )

    audit = {
        "purpose": (
            "Collapse comparative samples into "
            "Dennis-positive and Dennis-negative "
            "rows while retaining recalled L/R"
        ),
        "comparative_aggregation": (
            "equal-individual mean diploid "
            "ALT-allele dosage"
        ),
        "Dennis_positive_samples": int(
            group_sizes["Dennis-positive"]
        ),
        "Dennis_negative_samples": int(
            group_sizes["Dennis-negative"]
        ),
        "minimum_comparative_call_rate": (
            MIN_GROUP_CALL_RATE
        ),
        "focal_call_source": (
            "reciprocal-mixture posterior calls"
        ),
        "focal_minimum_pair_posterior": 0.95,
        "focal_raw_GT_fallback": False,
        "row_order": GROUP_ORDER,
        "modules": audit_modules,
        "combined_pdf": str(
            COMBINED_PDF.resolve()
        ),
        "plot_values": str(
            MATRIX_OUTPUT.resolve()
        ),
    }

    AUDIT_OUTPUT.write_text(
        json.dumps(
            audit,
            indent=2,
        )
        + "\n"
    )

    required_outputs = [
        COMBINED_PDF,
        MATRIX_OUTPUT,
        AUDIT_OUTPUT,
    ]

    for module in MODULES:
        required_outputs.extend(
            [
                (
                    OUTPUT_DIR
                    / (
                        f"optix_{module}_collapsed_"
                        "D_vs_nonD_plus_recalled_LR.pdf"
                    )
                ),
                (
                    OUTPUT_DIR
                    / (
                        f"optix_{module}_collapsed_"
                        "D_vs_nonD_plus_recalled_LR.png"
                    )
                ),
            ]
        )

    for output in required_outputs:
        if (
            not output.is_file()
            or output.stat().st_size == 0
        ):
            raise RuntimeError(
                f"Missing or empty output: {output}"
            )

    print(
        "COLLAPSED DENNIS GENOTYPE "
        "PLOTS COMPLETE"
    )

    print(
        "Dennis-positive samples\t"
        f"{group_sizes['Dennis-positive']}"
    )

    print(
        "Dennis-negative samples\t"
        f"{group_sizes['Dennis-negative']}"
    )

    for module, details in (
        audit_modules.items()
    ):
        print(
            f"{module}\tsites="
            f"{details['sites']}"
        )

    print(
        "COMBINED_PDF\t"
        f"{COMBINED_PDF.resolve()}"
    )

    print(
        "AUDIT\t"
        f"{AUDIT_OUTPUT.resolve()}"
    )


if __name__ == "__main__":
    main()
