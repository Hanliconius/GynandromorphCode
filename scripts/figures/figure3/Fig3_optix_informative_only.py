#!/usr/bin/env python3

from pathlib import Path
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams.update({
    "font.family": "DejaVu Sans",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "text.usetex": False,
})
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd


PROJECT = Path(
    os.environ.get(
        "GYNANDRO_ROOT",
        "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling",
    )
).expanduser().resolve()

SOURCE_DIR = (
    PROJECT
    / "analysis/original_502/checkpoints/"
    "krista_DR_regions_reciprocal_genoplots"
)

MANIFEST_PATH = (
    PROJECT
    / "analysis/figures/"
    "key_regions_population_context_genoplots/"
    "exact_BA_comparative_sample_manifest.tsv"
)

LR_PATH = (
    PROJECT
    / "analysis/original_502/"
    "reciprocal_L_R_mixture/"
    "optix_reciprocal_L_R_genotype_pairs.posterior95.tsv"
)

INFORMATIVE_PATH = (
    PROJECT
    / "analysis/original_502/"
    "optix_deltaAF_haplotype_analysis/"
    "optix_B_D_R_deltaAF_informative_SNPs.tsv"
)

OUTPUT_DIR = (
    PROJECT
    / "analysis/figure3/"
    "optix_informative_pattern_dosage"
)

COMBINED_PDF = (
    OUTPUT_DIR
    / "optix_all_Krista_loci_"
    "informative_only_collapsed_exact_BA_groups_plus_recalled_LR.pdf"
)

MATRIX_OUTPUT = (
    OUTPUT_DIR
    / "informative_only_plot_values.tsv.gz"
)

AUDIT_OUTPUT = (
    OUTPUT_DIR
    / "informative_only.audit.json"
)

MIN_CALL_RATE = 0.80

PANEL_CONFIG = {
    "Band": {
        "BA_pattern": "Band",
        "trait": "Band",
        "top_phenotype": 0,
        "labels": {
            0: "Band absent",
            1: "Band present",
        },
    },
    "Ray1": {
        "BA_pattern": "Ray",
        "trait": "Ray",
        "top_phenotype": 1,
        "labels": {
            0: "Ray absent",
            1: "Ray present",
        },
    },
    "Ray2": {
        "BA_pattern": "Ray",
        "trait": "Ray",
        "top_phenotype": 1,
        "labels": {
            0: "Ray absent",
            1: "Ray present",
        },
    },
    "Dennis1": {
        "BA_pattern": "Dennis",
        "trait": "Dennis",
        "top_phenotype": 1,
        "labels": {
            0: "Dennis absent",
            1: "Dennis present",
        },
    },
    "Dennis2": {
        "BA_pattern": "Dennis",
        "trait": "Dennis",
        "top_phenotype": 1,
        "labels": {
            0: "Dennis absent",
            1: "Dennis present",
        },
    },
}

DISPLAY_GT_TO_DOSAGE = {
    "A/A": 0.0,
    "A/P": 0.5,
    "P/A": 0.5,
    "P/P": 1.0,
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


def source_path(module):
    if module == "Band":
        return (
            PROJECT
            / "analysis/original_502/checkpoints/"
            "optix_band_original_region_"
            "reciprocal_genoplot.long.tsv.gz"
        )

    return (
        SOURCE_DIR
        / (
            f"optix_{module}_original_region_"
            "reciprocal_genoplot.long.tsv.gz"
        )
    )


def load_manifest():
    table = pd.read_csv(
        MANIFEST_PATH,
        sep="\t",
        dtype=str,
    )

    require_columns(
        table,
        [
            "BA_pattern",
            "sample",
            "phenotype",
            "population",
        ],
        "Exact-BA manifest",
    )

    table["phenotype"] = pd.to_numeric(
        table["phenotype"],
        errors="coerce",
    )

    return table


def load_informative_sites():
    table = pd.read_csv(
        INFORMATIVE_PATH,
        sep="\t",
        low_memory=False,
    )

    required = [
        "module",
        "chrom",
        "pos",
        "ref",
        "alt",
        "abs_delta_AF",
        "positive_call_rate",
        "negative_call_rate",
    ]

    require_columns(
        table,
        required,
        "Stage-117a informative-site table",
    )

    table["pos"] = pd.to_numeric(
        table["pos"],
        errors="raise",
    ).astype(int)

    for column in (
        "abs_delta_AF",
        "positive_call_rate",
        "negative_call_rate",
    ):
        table[column] = pd.to_numeric(
            table[column],
            errors="raise",
        )

    # Stage 117a's first informative-site definition:
    #   |delta AF| >= 0.80
    #   both comparative groups >= 0.80 callable
    #   calls scored with DP >= 2 and GQ >= 20
    # TN_Left/TN_Right were not used to select or polarize these sites.
    bad = table.loc[
        (table["abs_delta_AF"] < 0.80)
        | (table["positive_call_rate"] < 0.80)
        | (table["negative_call_rate"] < 0.80)
    ]

    if not bad.empty:
        raise RuntimeError(
            "Stage-117a informative-site table contains rows "
            "that do not satisfy the original selection rule"
        )

    keys = [
        "module",
        "chrom",
        "pos",
        "ref",
        "alt",
    ]

    table = table.drop_duplicates(
        subset=keys,
        keep="first",
    ).copy()

    return table


def load_lr_confidence():
    table = pd.read_csv(
        LR_PATH,
        sep="\t",
        low_memory=False,
    )

    keys = [
        "chrom",
        "pos",
        "ref",
        "alt",
    ]

    require_columns(
        table,
        keys
        + [
            "pair_posterior",
            "posterior95",
        ],
        "Reciprocal L/R table",
    )

    table["pos"] = pd.to_numeric(
        table["pos"],
        errors="raise",
    ).astype(int)

    table["pair_posterior"] = pd.to_numeric(
        table["pair_posterior"],
        errors="coerce",
    )

    table["posterior95"] = (
        table["posterior95"]
        .map(parse_boolean)
    )

    table["confident_reciprocal_call"] = (
        table["posterior95"]
        & table["pair_posterior"].ge(0.95)
    )

    if table.duplicated(keys).any():
        raise RuntimeError(
            "Duplicate sites in reciprocal L/R table"
        )

    return table[
        keys
        + [
            "pair_posterior",
            "confident_reciprocal_call",
        ]
    ]


def prepare_panel(
    module,
    config,
    manifest,
    lr_confidence,
    informative_sites,
):
    path = source_path(module)

    table = pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        low_memory=False,
    )

    # The dedicated Band table contains only the Band panel,
    # so its older schema does not include a panel column.
    if module == "Band" and "panel" not in table.columns:
        table["panel"] = "Band"

    keys = [
        "chrom",
        "pos",
        "ref",
        "alt",
    ]

    require_columns(
        table,
        [
            "panel",
            "chrom",
            "pos",
            "ref",
            "alt",
            "site_label",
            "sample",
            "call_source",
            "display_GT",
        ],
        f"{module} long table",
    )

    table = table.loc[
        table["panel"].eq(module)
    ].copy()

    table["pos"] = pd.to_numeric(
        table["pos"],
        errors="raise",
    ).astype(int)

    all_site_count = int(
        table[keys]
        .drop_duplicates()
        .shape[0]
    )

    broad_module = config["BA_pattern"]

    informative_keys = (
        informative_sites.loc[
            informative_sites["module"].eq(
                broad_module
            ),
            keys,
        ]
        .drop_duplicates()
    )

    table = table.merge(
        informative_keys,
        on=keys,
        how="inner",
        validate="many_to_one",
    )

    if table.empty:
        raise RuntimeError(
            f"{module}: no stage-117a informative sites "
            "overlap this locus"
        )

    sites = (
        table[
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
            f"{module}: duplicate site definitions"
        )

    sites["site_index"] = np.arange(
        len(sites),
        dtype=int,
    )

    cohort = manifest.loc[
        manifest["BA_pattern"].eq(
            config["BA_pattern"]
        )
    ].copy()

    # Dennis and Ray phenotypes are stored directly in the
    # exact-BA manifest. Band membership is stored in the
    # established individual-level Band plot_group field,
    # while the manifest identifies the restricted BA cohort.
    if module == "Band":
        band_status = (
            table[
                [
                    "sample",
                    "plot_group",
                ]
            ]
            .drop_duplicates()
            .copy()
        )

        status_map = {
            "Band absent (B=0)": 0,
            "Band present (B=1)": 1,
        }

        band_status["band_phenotype"] = (
            band_status["plot_group"]
            .map(status_map)
        )

        band_status = band_status.loc[
            band_status[
                "band_phenotype"
            ].isin([0, 1])
        ][
            [
                "sample",
                "band_phenotype",
            ]
        ].copy()

        conflicts = (
            band_status.groupby("sample")[
                "band_phenotype"
            ]
            .nunique()
            .gt(1)
        )

        if conflicts.any():
            raise RuntimeError(
                "Band table contains conflicting "
                "phenotype assignments"
            )

        band_status = (
            band_status.drop_duplicates(
                subset=["sample"],
                keep="first",
            )
        )

        cohort = (
            cohort.drop(
                columns=["phenotype"]
            )
            .merge(
                band_status,
                on="sample",
                how="left",
                validate="one_to_one",
            )
            .rename(
                columns={
                    "band_phenotype":
                        "phenotype",
                }
            )
        )

    cohort = cohort.loc[
        cohort["phenotype"].isin([0, 1])
    ].copy()

    conflicts = (
        cohort.groupby("sample")["phenotype"]
        .nunique()
        .gt(1)
    )

    if conflicts.any():
        raise RuntimeError(
            f"{module}: conflicting phenotype "
            "assignments in exact-BA manifest"
        )

    cohort = cohort.drop_duplicates(
        subset=["sample"],
        keep="first",
    )

    group_sizes = (
        cohort.groupby("phenotype")["sample"]
        .nunique()
        .to_dict()
    )

    if (
        group_sizes.get(0, 0) == 0
        or group_sizes.get(1, 0) == 0
    ):
        raise RuntimeError(
            f"{module}: missing phenotype group"
        )

    comparative = table.merge(
        cohort[
            [
                "sample",
                "phenotype",
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

    comparative["pattern_dosage"] = (
        comparative["display_GT"]
        .astype(str)
        .map(DISPLAY_GT_TO_DOSAGE)
    )

    records = []

    for phenotype in (0, 1):
        label = config["labels"][phenotype]

        subset = comparative.loc[
            comparative["phenotype"].eq(
                phenotype
            )
        ]

        summary = (
            subset.groupby(
                keys,
                as_index=False,
                sort=False,
            )["pattern_dosage"]
            .agg(
                mean_pattern_dosage="mean",
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
            group_sizes[phenotype]
        )

        summary["call_rate"] = (
            summary["called_samples"]
            / summary["total_samples"]
        )

        summary.loc[
            summary["call_rate"]
            < MIN_CALL_RATE,
            "mean_pattern_dosage",
        ] = np.nan

        summary["plot_group"] = label
        summary["phenotype"] = phenotype
        summary["source"] = (
            "equal-individual mean display_GT "
            "within exact-BA cohort"
        )
        summary["pair_posterior"] = np.nan
        summary["module"] = module

        records.append(summary)

    # The Band source contains both the original mixed-sample
    # VCF call and the posterior-filtered reciprocal-mixture
    # recall for each focal sample/site. The original GT must
    # never be used for TN_Left or TN_Right.
    reciprocal_source = (
        "reciprocal mixture "
        "(posterior >= 0.95)"
    )

    focal = table.loc[
        table["sample"].isin(
            ["TN_Left", "TN_Right"]
        )
        & table["call_source"].eq(
            reciprocal_source
        )
    ].copy()

    if focal.empty:
        raise RuntimeError(
            f"{module}: no posterior-filtered "
            "reciprocal-mixture focal calls"
        )

    unexpected_sources = sorted(
        focal["call_source"]
        .dropna()
        .unique()
        .tolist()
    )

    if unexpected_sources != [
        reciprocal_source
    ]:
        raise RuntimeError(
            f"{module}: unexpected focal call "
            f"sources: {unexpected_sources}"
        )

    if focal.duplicated(
        keys + ["sample"]
    ).any():
        raise RuntimeError(
            f"{module}: duplicate focal sample/site rows"
        )

    focal = focal.merge(
        lr_confidence,
        on=keys,
        how="left",
        validate="many_to_one",
    )

    focal["pattern_dosage"] = (
        focal["display_GT"]
        .astype(str)
        .map(DISPLAY_GT_TO_DOSAGE)
    )

    focal["confident_reciprocal_call"] = (
        focal["confident_reciprocal_call"]
        .fillna(False)
    )

    focal["confident_reciprocal_call"] = (

        focal["confident_reciprocal_call"]

        .fillna(False)

        .astype(bool)

    )


    focal.loc[

        ~focal["confident_reciprocal_call"],
        "pattern_dosage",
    ] = np.nan

    for sample in (
        "TN_Left",
        "TN_Right",
    ):
        sample_rows = focal.loc[
            focal["sample"].eq(sample)
        ][
            keys
            + [
                "pattern_dosage",
                "pair_posterior",
                "confident_reciprocal_call",
                "call_source",
            ]
        ].copy()

        summary = sites.merge(
            sample_rows,
            on=keys,
            how="left",
            validate="one_to_one",
        )

        confident = (
            summary[
                "confident_reciprocal_call"
            ]
            .fillna(False)
        )

        summary["mean_pattern_dosage"] = (
            summary["pattern_dosage"]
        )

        summary.loc[
            ~confident,
            "mean_pattern_dosage",
        ] = np.nan

        summary["called_samples"] = (
            confident.astype(int)
        )
        summary["total_samples"] = 1
        summary["call_rate"] = (
            confident.astype(float)
        )
        summary["plot_group"] = sample
        summary["phenotype"] = np.nan
        summary["source"] = (
            "posterior>=0.95 reciprocal-mixture "
            "display_GT; no raw-GT fallback"
        )
        summary["module"] = module

        records.append(summary)

    values = pd.concat(
        records,
        ignore_index=True,
        sort=False,
    )

    top_phenotype = int(
        config["top_phenotype"]
    )

    bottom_phenotype = 1 - top_phenotype

    row_order = [
        config["labels"][top_phenotype],
        "TN_Left",
        "TN_Right",
        config["labels"][bottom_phenotype],
    ]

    values["plot_group"] = pd.Categorical(
        values["plot_group"],
        categories=row_order,
        ordered=True,
    )

    values = values.sort_values(
        [
            "plot_group",
            "site_index",
        ],
        kind="mergesort",
    ).reset_index(drop=True)

    return (
        sites,
        values,
        row_order,
        group_sizes,
        all_site_count,
    )


def make_colormap():
    # Phase-free dosage scale.  The population rows can take any value
    # between 0 and 1 because they are equal-individual means; TN rows
    # remain discrete at 0, 0.5, or 1 because they are recalled genotypes.
    colourmap = (
        LinearSegmentedColormap.from_list(
            "pattern_dosage_unphased",
            [
                (0.0, "#2C9BCB"),
                (0.5, "#F3E6C8"),
                (1.0, "#B61F72"),
            ],
        )
    )

    colourmap.set_bad("#FFFFFF")

    return colourmap


def draw_panel(
    module,
    config,
    sites,
    values,
    row_order,
    group_sizes,
):
    matrix = (
        values.pivot(
            index="plot_group",
            columns="site_index",
            values="mean_pattern_dosage",
        )
        .reindex(row_order)
        .reindex(
            columns=sites["site_index"]
        )
    )

    width = max(
        11.0,
        min(
            21.0,
            7.5 + len(sites) / 45.0,
        ),
    )

    figure, axis = plt.subplots(
        figsize=(width, 4.1),
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

    # The population contrasts sandwich Left and Right.
    axis.axhline(
        0.5,
        color="#303030",
        linewidth=1.5,
    )

    axis.axhline(
        2.5,
        color="#303030",
        linewidth=1.5,
    )

    phenotype_by_label = {
        label: phenotype
        for phenotype, label
        in config["labels"].items()
    }

    y_labels = []

    for label in row_order:
        if label in phenotype_by_label:
            phenotype = phenotype_by_label[label]

            y_labels.append(
                f"{label} "
                f"(n={group_sizes[phenotype]})"
            )
        elif label == "TN_Left":
            y_labels.append(
                "TN_Left recalled"
            )
        else:
            y_labels.append(
                "TN_Right recalled"
            )

    axis.set_yticks(
        np.arange(4)
    )

    axis.set_yticklabels(
        y_labels
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

    axis.set_xticks(
        tick_indices
    )

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

    axis.text(
        0.0,
        1.105,
        (
            f"{module}: "
            f"{config['trait']} phenotype contrast"
        ),
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=15,
        fontweight="bold",
    )

    axis.text(
        0.0,
        1.045,
        (
            "Exact-BA comparison cohort; population rows are equal-individual "
            "means; TN rows are posterior>=0.95 recalled genotypes."
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
            "0  A/A",
            "0.5  P/A",
            "1  P/P",
        ]
    )

    colourbar.set_label(
        "Pattern-associated allele dosage"
    )

    for spine in axis.spines.values():
        spine.set_visible(False)

    return figure


def main():
    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    required_inputs = [
        MANIFEST_PATH,
        LR_PATH,
        INFORMATIVE_PATH,
    ] + [
        source_path(module)
        for module in PANEL_CONFIG
    ]

    for path in required_inputs:
        if (
            not path.is_file()
            or path.stat().st_size == 0
        ):
            raise RuntimeError(
                f"Missing or empty input: {path}"
            )

    manifest = load_manifest()
    lr_confidence = load_lr_confidence()
    informative_sites = load_informative_sites()

    print(
        "INFORMATIVE_SITE_SOURCE\t"
        f"{INFORMATIVE_PATH.resolve()}"
    )
    print(
        "INFORMATIVE_SELECTION\t"
        "|delta AF|>=0.80; both groups call-rate>=0.80; "
        "DP>=2; GQ>=20; TN excluded from selection/polarization"
    )

    all_values = []
    audit_panels = {}

    with PdfPages(COMBINED_PDF) as combined:
        for module, config in (
            PANEL_CONFIG.items()
        ):
            (
                sites,
                values,
                row_order,
                group_sizes,
                all_site_count,
            ) = prepare_panel(
                module,
                config,
                manifest,
                lr_confidence,
                informative_sites,
            )

            figure = draw_panel(
                module,
                config,
                sites,
                values,
                row_order,
                group_sizes,
            )

            panel_pdf = (
                OUTPUT_DIR
                / (
                    f"optix_{module}_informative_only_"
                    "collapsed_exact_BA_groups_plus_recalled_LR.pdf"
                )
            )

            panel_png = panel_pdf.with_suffix(
                ".png"
            )

            figure.savefig(
                panel_pdf,
                bbox_inches="tight",
            )

            figure.savefig(
                panel_png,
                dpi=300,
                bbox_inches="tight",
            )

            combined.savefig(
                figure,
                bbox_inches="tight",
            )

            plt.close(figure)

            all_values.append(values)

            audit_panels[module] = {
                "BA_pattern": config[
                    "BA_pattern"
                ],
                "all_locus_sites_before_informative_filter": int(
                    all_site_count
                ),
                "informative_sites": int(len(sites)),
                "phenotype_0_samples": int(
                    group_sizes[0]
                ),
                "phenotype_1_samples": int(
                    group_sizes[1]
                ),
                "row_order": row_order,
                "retained_cells": {
                    group: int(
                        values.loc[
                            values[
                                "plot_group"
                            ].eq(group),
                            "mean_pattern_dosage",
                        ].notna().sum()
                    )
                    for group in row_order
                },
                "pdf": str(
                    panel_pdf.resolve()
                ),
                "png": str(
                    panel_png.resolve()
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
        "phenotype",
        "mean_pattern_dosage",
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

    finite = combined_values[
        "mean_pattern_dosage"
    ].dropna()

    if not finite.between(
        0.0,
        1.0,
        inclusive="both",
    ).all():
        raise RuntimeError(
            "Pattern dosage outside 0-1"
        )

    audit = {
        "panels": list(
            PANEL_CONFIG
        ),
        "site_filter": (
            "original stage-117a comparative informative sites only"
        ),
        "site_filter_source": str(
            INFORMATIVE_PATH.resolve()
        ),
        "site_filter_rule": {
            "minimum_abs_delta_AF": 0.80,
            "minimum_group_call_rate": 0.80,
            "minimum_DP": 2,
            "minimum_GQ": 20,
            "focal_samples_used_for_selection": False,
            "focal_samples_used_for_polarization": False,
        },
        "comparative_cohorts": (
            "phenotype-specific exact-BA cohorts"
        ),
        "aggregation": (
            "equal-individual mean pattern-associated "
            "allele dosage"
        ),
        "minimum_comparative_call_rate": (
            MIN_CALL_RATE
        ),
        "genotype_orientation": (
            "phase-free pattern-associated dosage from existing display_GT: "
            "A/A=0, P/A=0.5, P/P=1"
        ),
        "palette": {
            "A/A_0": "#2C9BCB",
            "P/A_0.5": "#F3E6C8",
            "P/P_1": "#B61F72",
            "missing": "#FFFFFF",
        },
        "population_rows": (
            "continuous equal-individual mean dosage on 0-1 scale"
        ),
        "focal_rows": (
            "discrete recalled genotype dosage 0/0.5/1; no phasing required"
        ),
        "focal_calls": (
            "reciprocal-mixture posterior>=0.95"
        ),
        "focal_raw_GT_fallback": False,
        "panels_detail": audit_panels,
        "combined_pdf": str(
            COMBINED_PDF.resolve()
        ),
        "values": str(
            MATRIX_OUTPUT.resolve()
        ),
        "supersedes_stage_211": (
            "Stage 211 used raw ALT orientation; "
            "Stage 212 uses pattern-associated "
            "display_GT orientation"
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

    for output in required_outputs:
        if (
            not output.is_file()
            or output.stat().st_size == 0
        ):
            raise RuntimeError(
                f"Missing or empty output: {output}"
            )

    print(
        "FIGURE 3 INFORMATIVE-ONLY OPTIX PLOTS COMPLETE"
    )

    for module, details in (
        audit_panels.items()
    ):
        print(
            f"{module}"
            f"\tall_sites={details['all_locus_sites_before_informative_filter']}"
            f"\tinformative_sites={details['informative_sites']}"
            f"\trows={' | '.join(details['row_order'])}"
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
