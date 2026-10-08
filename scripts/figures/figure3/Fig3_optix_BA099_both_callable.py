#!/usr/bin/env python3
"""Figure 3 optix genoplots using the BA>0.99 site set from the 170A/170B work.

Selection is independent of the focal genotype state:
  1. comparative BA > 0.99 for the relevant trait
  2. both TN_Left and TN_Right have a reciprocal-mixture pair call with
     posterior >= 0.95
  3. no requirement that the recalled TN genotypes match any expected pattern

The five current Krista loci are retained (Band, Ray1, Ray2, Dennis1, Dennis2).
Population rows are equal-individual mean pattern-associated allele dosage.
TN rows are the recalled genotypes themselves. No phasing is used.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd

matplotlib.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "text.usetex": False,
    "savefig.facecolor": "white",
})

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

BAND_BA_PATH = (
    PROJECT
    / "analysis/B_western_postman_vs_western_Amazon_DR_reanalysis/"
    "checkpoints/B_western_BA_and_adjusted_LR.tsv"
)

DR_BA_PATH = (
    PROJECT
    / "analysis/original_502/Wallbank_Dennis_Ray_BA/"
    "Wallbank_Dennis_Ray_dominant_BA.tsv"
)

OUTPUT_DIR = (
    PROJECT
    / "analysis/figure3/"
    "optix_BA099_callable_pattern_dosage"
)

COMBINED_PDF = (
    OUTPUT_DIR
    / "optix_BAgt099_both_recalled_callable_five_loci.pdf"
)

SITE_OUTPUT = OUTPUT_DIR / "BAgt099_both_recalled_callable_sites.tsv"
MATRIX_OUTPUT = OUTPUT_DIR / "BAgt099_both_recalled_callable_plot_values.tsv.gz"
AUDIT_OUTPUT = OUTPUT_DIR / "BAgt099_both_recalled_callable.audit.json"

MIN_CALL_RATE = 0.80
MIN_BA = 0.99

PANEL_CONFIG = {
    "Band": {
        "BA_pattern": "Band",
        "trait": "Band",
        "top_phenotype": 0,
        "labels": {0: "Band absent", 1: "Band present"},
    },
    "Ray1": {
        "BA_pattern": "Ray",
        "trait": "Ray",
        "top_phenotype": 1,
        "labels": {0: "Ray absent", 1: "Ray present"},
    },
    "Ray2": {
        "BA_pattern": "Ray",
        "trait": "Ray",
        "top_phenotype": 1,
        "labels": {0: "Ray absent", 1: "Ray present"},
    },
    "Dennis1": {
        "BA_pattern": "Dennis",
        "trait": "Dennis",
        "top_phenotype": 1,
        "labels": {0: "Dennis absent", 1: "Dennis present"},
    },
    "Dennis2": {
        "BA_pattern": "Dennis",
        "trait": "Dennis",
        "top_phenotype": 1,
        "labels": {0: "Dennis absent", 1: "Dennis present"},
    },
}

KEY = ["chrom", "pos", "ref", "alt"]


def require_columns(table, required, label):
    missing = sorted(set(required) - set(table.columns))
    if missing:
        raise RuntimeError(f"{label}: missing columns: {missing}")


def as_bool(value):
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def canonical_gt(value):
    if value is None or pd.isna(value):
        return None

    text = str(value).strip().replace("|", "/")
    if text == "1/0":
        text = "0/1"

    return text if text in {"0/0", "0/1", "1/1"} else None


def alt_dosage(value):
    gt = canonical_gt(value)
    return {
        "0/0": 0.0,
        "0/1": 0.5,
        "1/1": 1.0,
    }.get(gt, np.nan)


def pattern_is_alt(orientation):
    orientation = str(orientation).strip()

    if orientation == "ALT_carrier_predicts_presence":
        return True

    if orientation == "REF_carrier_predicts_presence":
        return False

    raise RuntimeError(f"Unrecognized BA orientation: {orientation!r}")


def pattern_dosage(gt, is_alt):
    dose = alt_dosage(gt)

    if not np.isfinite(dose):
        return np.nan

    return dose if bool(is_alt) else 1.0 - dose


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
        / f"optix_{module}_original_region_reciprocal_genoplot.long.tsv.gz"
    )


def load_manifest():
    table = pd.read_csv(MANIFEST_PATH, sep="\t", dtype=str)

    require_columns(
        table,
        ["BA_pattern", "sample", "phenotype", "population"],
        "Exact-BA manifest",
    )

    table["phenotype"] = pd.to_numeric(table["phenotype"], errors="coerce")
    return table


def load_recalled_pairs():
    table = pd.read_csv(LR_PATH, sep="\t", low_memory=False)

    require_columns(
        table,
        KEY + ["pair_posterior", "posterior95"],
        "Reciprocal L/R table",
    )

    table["pos"] = pd.to_numeric(table["pos"], errors="raise").astype(int)
    table["pair_posterior"] = pd.to_numeric(
        table["pair_posterior"], errors="coerce"
    )
    table["posterior95"] = table["posterior95"].map(as_bool)

    pair_column = next(
        (
            column
            for column in (
                "mixture_pair_canonical",
                "mixture_pair_state",
                "pair_state",
            )
            if column in table.columns
        ),
        None,
    )

    if pair_column is None:
        candidates = [
            column
            for column in table.columns
            if "pair" in column.lower()
            and "original" not in column.lower()
            and column != "pair_posterior"
        ]

        if len(candidates) != 1:
            raise RuntimeError(
                "Could not identify one reciprocal pair-state column. "
                "Candidates: " + ",".join(candidates)
            )

        pair_column = candidates[0]

    parts = table[pair_column].astype(str).str.split(
        "|", n=1, expand=True, regex=False
    )

    if parts.shape[1] != 2:
        raise RuntimeError(f"Could not split reciprocal pair states in {pair_column}")

    table["TN_Left_GT"] = parts[0].map(canonical_gt)
    table["TN_Right_GT"] = parts[1].map(canonical_gt)

    table["both_recalled_callable"] = (
        table["posterior95"]
        & table["pair_posterior"].ge(0.95)
        & table["TN_Left_GT"].notna()
        & table["TN_Right_GT"].notna()
    )

    if table.duplicated(KEY).any():
        raise RuntimeError("Duplicate sites in reciprocal L/R table")

    return table[
        KEY
        + [
            "pair_posterior",
            "TN_Left_GT",
            "TN_Right_GT",
            "both_recalled_callable",
        ]
    ].copy()


def load_ba_sites():
    band = pd.read_csv(BAND_BA_PATH, sep="\t", low_memory=False)

    require_columns(
        band,
        KEY
        + [
            "analysis_cohort",
            "scheme",
            "balanced_accuracy",
            "selected_orientation",
            "score_eligible",
        ],
        "Analysis-142 Band BA table",
    )

    band["pos"] = pd.to_numeric(band["pos"], errors="raise").astype(int)
    band["balanced_accuracy"] = pd.to_numeric(
        band["balanced_accuracy"], errors="coerce"
    )

    band = band.loc[
        band["analysis_cohort"].eq("WGS_plus_CS")
        & band["scheme"].eq("DP2_GQ20")
        & band["score_eligible"].map(as_bool)
        & band["balanced_accuracy"].gt(MIN_BA)
    ].copy()

    band["module"] = "Band"
    band["BA_source"] = "analysis142_western_postman_vs_western_Amazon_DR"

    dr = pd.read_csv(DR_BA_PATH, sep="\t", low_memory=False)

    require_columns(
        dr,
        KEY
        + [
            "module",
            "balanced_accuracy",
            "selected_orientation",
        ],
        "Wallbank Dennis/Ray BA table",
    )

    dr["pos"] = pd.to_numeric(dr["pos"], errors="raise").astype(int)
    dr["balanced_accuracy"] = pd.to_numeric(
        dr["balanced_accuracy"], errors="coerce"
    )

    dr = dr.loc[
        dr["module"].isin(["Dennis", "Ray"])
        & dr["balanced_accuracy"].gt(MIN_BA)
    ].copy()

    dr["BA_source"] = "Wallbank_Dennis_Ray_dominant_BA"

    columns = (
        ["module"]
        + KEY
        + [
            "balanced_accuracy",
            "selected_orientation",
            "BA_source",
        ]
    )

    selected = pd.concat(
        [band[columns], dr[columns]],
        ignore_index=True,
        sort=False,
    )

    if selected.duplicated(["module"] + KEY).any():
        dup = selected.loc[
            selected.duplicated(["module"] + KEY, keep=False),
            ["module"] + KEY,
        ].head(20)
        raise RuntimeError(
            "Duplicate BA-selected module/site rows:\n"
            + dup.to_string(index=False)
        )

    selected["pattern_is_alt"] = selected["selected_orientation"].map(
        pattern_is_alt
    )

    return selected.sort_values(
        ["module", "pos", "ref", "alt"],
        kind="mergesort",
    ).reset_index(drop=True)


def make_colormap():
    cmap = LinearSegmentedColormap.from_list(
        "phase_free_pattern_dosage",
        [
            (0.0, "#2C9BCB"),
            (0.5, "#F3E6C8"),
            (1.0, "#B61F72"),
        ],
    )
    cmap.set_bad("#FFFFFF")
    return cmap


def prepare_panel(
    module,
    config,
    manifest,
    recalled,
    ba_sites,
):
    table = pd.read_csv(
        source_path(module),
        sep="\t",
        compression="gzip",
        low_memory=False,
    )

    if module == "Band" and "panel" not in table.columns:
        table["panel"] = "Band"

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
            "plot_group",
            "call_source",
            "raw_GT",
            "display_GT",
        ],
        f"{module} long table",
    )

    table = table.loc[table["panel"].eq(module)].copy()
    table["pos"] = pd.to_numeric(table["pos"], errors="raise").astype(int)

    all_sites = (
        table[KEY + ["site_label"]]
        .drop_duplicates()
        .sort_values(KEY, kind="mergesort")
        .reset_index(drop=True)
    )

    broad_module = config["BA_pattern"]

    candidates = ba_sites.loc[
        ba_sites["module"].eq(broad_module)
    ].copy()

    candidates = candidates.merge(
        recalled,
        on=KEY,
        how="inner",
        validate="one_to_one",
    )

    candidates = candidates.loc[
        candidates["both_recalled_callable"]
    ].copy()

    selected_sites = all_sites.merge(
        candidates[
            KEY
            + [
                "balanced_accuracy",
                "selected_orientation",
                "pattern_is_alt",
                "BA_source",
                "pair_posterior",
                "TN_Left_GT",
                "TN_Right_GT",
            ]
        ],
        on=KEY,
        how="inner",
        validate="one_to_one",
    )

    if selected_sites.empty:
        raise RuntimeError(
            f"{module}: zero BA>0.99 sites with confident recalled calls "
            "overlap this locus"
        )

    selected_sites = selected_sites.sort_values(
        KEY,
        kind="mergesort",
    ).reset_index(drop=True)

    selected_sites["site_index"] = np.arange(
        len(selected_sites), dtype=int
    )

    table = table.merge(
        selected_sites[KEY],
        on=KEY,
        how="inner",
        validate="many_to_one",
    )

    cohort = manifest.loc[
        manifest["BA_pattern"].eq(broad_module)
    ].copy()

    if module == "Band":
        band_status = (
            table[["sample", "plot_group"]]
            .drop_duplicates()
            .copy()
        )

        status_map = {
            "Band absent (B=0)": 0,
            "Band present (B=1)": 1,
        }

        band_status["band_phenotype"] = band_status["plot_group"].map(
            status_map
        )

        band_status = (
            band_status.loc[
                band_status["band_phenotype"].isin([0, 1]),
                ["sample", "band_phenotype"],
            ]
            .drop_duplicates("sample")
        )

        cohort = (
            cohort.drop(columns=["phenotype"])
            .merge(
                band_status,
                on="sample",
                how="left",
                validate="one_to_one",
            )
            .rename(columns={"band_phenotype": "phenotype"})
        )

    cohort = cohort.loc[cohort["phenotype"].isin([0, 1])].copy()
    cohort = cohort.drop_duplicates("sample")

    group_sizes = (
        cohort.groupby("phenotype")["sample"].nunique().to_dict()
    )

    if group_sizes.get(0, 0) == 0 or group_sizes.get(1, 0) == 0:
        raise RuntimeError(f"{module}: missing phenotype group")

    comparative = table.merge(
        cohort[["sample", "phenotype"]],
        on="sample",
        how="inner",
        validate="many_to_one",
    )

    if comparative.duplicated(KEY + ["sample"]).any():
        raise RuntimeError(f"{module}: duplicate comparative sample/site rows")

    comparative = comparative.merge(
        selected_sites[KEY + ["pattern_is_alt"]],
        on=KEY,
        how="left",
        validate="many_to_one",
    )

    comparative["pattern_dosage"] = [
        pattern_dosage(gt, is_alt)
        for gt, is_alt in zip(
            comparative["raw_GT"],
            comparative["pattern_is_alt"],
        )
    ]

    # Preserve the existing regional genoplot QC: cells shown as N/Uncertain
    # in the source table remain missing in the collapsed population mean.
    source_called = comparative["display_GT"].isin(
        ["A/A", "A/P", "P/A", "P/P"]
    )
    comparative.loc[~source_called, "pattern_dosage"] = np.nan

    records = []

    for phenotype in (0, 1):
        label = config["labels"][phenotype]

        subset = comparative.loc[
            comparative["phenotype"].eq(phenotype)
        ]

        summary = (
            subset.groupby(KEY, as_index=False, sort=False)["pattern_dosage"]
            .agg(
                mean_pattern_dosage="mean",
                called_samples="count",
            )
        )

        summary = selected_sites.merge(
            summary,
            on=KEY,
            how="left",
            validate="one_to_one",
        )

        summary["called_samples"] = (
            summary["called_samples"].fillna(0).astype(int)
        )
        summary["total_samples"] = int(group_sizes[phenotype])
        summary["call_rate"] = (
            summary["called_samples"] / summary["total_samples"]
        )

        summary.loc[
            summary["call_rate"] < MIN_CALL_RATE,
            "mean_pattern_dosage",
        ] = np.nan

        summary["plot_group"] = label
        summary["phenotype"] = phenotype
        summary["source"] = (
            "equal-individual mean genotype dosage in exact-BA display cohort"
        )
        summary["module"] = module

        records.append(summary)

    for sample in ("TN_Left", "TN_Right"):
        summary = selected_sites.copy()
        gt_column = f"{sample}_GT"

        summary["mean_pattern_dosage"] = [
            pattern_dosage(gt, is_alt)
            for gt, is_alt in zip(
                summary[gt_column],
                summary["pattern_is_alt"],
            )
        ]

        if summary["mean_pattern_dosage"].isna().any():
            raise RuntimeError(
                f"{module}: selected sites unexpectedly contain missing "
                f"{sample} recalled genotype"
            )

        summary["called_samples"] = 1
        summary["total_samples"] = 1
        summary["call_rate"] = 1.0
        summary["plot_group"] = sample
        summary["phenotype"] = np.nan
        summary["source"] = (
            "reciprocal-mixture recalled genotype; pair posterior>=0.95"
        )
        summary["module"] = module

        records.append(summary)

    values = pd.concat(records, ignore_index=True, sort=False)

    top = int(config["top_phenotype"])
    bottom = 1 - top

    row_order = [
        config["labels"][top],
        "TN_Left",
        "TN_Right",
        config["labels"][bottom],
    ]

    values["plot_group"] = pd.Categorical(
        values["plot_group"],
        categories=row_order,
        ordered=True,
    )

    values = values.sort_values(
        ["plot_group", "site_index"],
        kind="mergesort",
    ).reset_index(drop=True)

    return (
        selected_sites,
        values,
        row_order,
        group_sizes,
        len(all_sites),
        len(candidates),
    )


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
        .reindex(columns=sites["site_index"])
    )

    width = max(7.5, min(15.0, 6.0 + len(sites) / 15.0))

    fig, ax = plt.subplots(
        figsize=(width, 4.4),
    )

    # Explicit margins are more reliable than constrained_layout for
    # these very sparse 4-row panels with a colourbar and long labels.
    fig.subplots_adjust(
        left=0.20,
        right=0.88,
        bottom=0.22,
        top=0.78,
    )

    im = ax.imshow(
        matrix.to_numpy(dtype=float),
        cmap=make_colormap(),
        vmin=0.0,
        vmax=1.0,
        aspect="auto",
        interpolation="nearest",
    )

    ax.axhline(0.5, color="#303030", linewidth=1.2)
    ax.axhline(2.5, color="#303030", linewidth=1.2)

    phenotype_by_label = {
        label: phenotype
        for phenotype, label in config["labels"].items()
    }

    y_labels = []

    for label in row_order:
        if label in phenotype_by_label:
            phenotype = phenotype_by_label[label]
            y_labels.append(
                f"{label} (n={group_sizes[phenotype]})"
            )
        elif label == "TN_Left":
            y_labels.append("TN_Left recalled")
        else:
            y_labels.append("TN_Right recalled")

    ax.set_yticks(np.arange(4))
    ax.set_yticklabels(y_labels)

    tick_count = min(10, len(sites))
    tick_indices = np.unique(
        np.linspace(0, len(sites) - 1, tick_count)
        .round()
        .astype(int)
    )

    ax.set_xticks(tick_indices)
    ax.set_xticklabels(
        [
            f"{sites.iloc[index]['pos'] / 1_000_000:.4f}"
            for index in tick_indices
        ],
        rotation=45,
        ha="right",
    )

    ax.set_xlabel("Position on Hmr18 (Mb)")

    ax.text(
        0.0,
        1.105,
        f"{module}: BA > 0.99, recalled on both sides",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=14,
        fontweight="bold",
    )

    ax.text(
        0.0,
        1.045,
        (
            "Comparative BA selection is independent of TN genotype state; "
            "sites lacking a posterior>=0.95 reciprocal pair call are omitted."
        ),
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.5,
        color="#555555",
    )

    cbar = fig.colorbar(
        im,
        ax=ax,
        fraction=0.028,
        pad=0.018,
    )

    cbar.set_ticks([0.0, 0.5, 1.0])
    cbar.set_ticklabels(["0  A/A", "0.5  P/A", "1  P/P"])
    cbar.set_label("Pattern-associated allele dosage")

    for spine in ax.spines.values():
        spine.set_visible(False)

    return fig


def main():
    required = [
        MANIFEST_PATH,
        LR_PATH,
        BAND_BA_PATH,
        DR_BA_PATH,
    ] + [
        source_path(module)
        for module in PANEL_CONFIG
    ]

    missing = [
        str(path)
        for path in required
        if not path.is_file() or path.stat().st_size == 0
    ]

    if missing:
        raise RuntimeError(
            "Missing/empty required inputs:\n  "
            + "\n  ".join(missing)
        )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    manifest = load_manifest()
    recalled = load_recalled_pairs()
    ba_sites = load_ba_sites()

    print("SITE_SELECTION")
    print("  comparative BA > 0.99")
    print("  both TN sides reciprocal-mixture pair posterior >= 0.95")
    print("  TN genotype state NOT used for selection")
    print()

    print("BA_GT_099_BEFORE_TN_CALLABILITY")
    print(
        ba_sites.groupby("module")
        .size()
        .to_string()
    )
    print()

    all_values = []
    all_sites = []
    details = {}

    with PdfPages(COMBINED_PDF) as combined:
        for module, config in PANEL_CONFIG.items():
            (
                sites,
                values,
                row_order,
                group_sizes,
                locus_sites_before,
                broad_candidates_after_callability,
            ) = prepare_panel(
                module,
                config,
                manifest,
                recalled,
                ba_sites,
            )

            fig = draw_panel(
                module,
                config,
                sites,
                values,
                row_order,
                group_sizes,
            )

            panel_pdf = (
                OUTPUT_DIR
                / f"optix_{module}_BAgt099_both_recalled_callable.pdf"
            )

            fig.savefig(panel_pdf, bbox_inches="tight")
            combined.savefig(fig, bbox_inches="tight")
            plt.close(fig)

            all_values.append(values)
            all_sites.append(
                sites.assign(
                    panel=module,
                    module=config["BA_pattern"],
                )
            )

            details[module] = {
                "all_sites_in_locus": int(locus_sites_before),
                "broad_BAgt099_both_recalled_callable": int(
                    broad_candidates_after_callability
                ),
                "selected_sites_in_locus": int(len(sites)),
                "phenotype_0_samples": int(group_sizes[0]),
                "phenotype_1_samples": int(group_sizes[1]),
                "position_min": int(sites["pos"].min()),
                "position_max": int(sites["pos"].max()),
                "pdf": str(panel_pdf.resolve()),
            }

            print(
                f"{module}"
                f"\tall_locus_sites={locus_sites_before}"
                f"\tselected={len(sites)}"
                f"\trange={int(sites['pos'].min())}-{int(sites['pos'].max())}"
            )

    site_table = pd.concat(
        all_sites,
        ignore_index=True,
        sort=False,
    )

    site_table[
        [
            "panel",
            "module",
            "chrom",
            "pos",
            "ref",
            "alt",
            "site_label",
            "balanced_accuracy",
            "selected_orientation",
            "pattern_is_alt",
            "BA_source",
            "pair_posterior",
            "TN_Left_GT",
            "TN_Right_GT",
        ]
    ].to_csv(
        SITE_OUTPUT,
        sep="\t",
        index=False,
    )

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
        "balanced_accuracy",
        "selected_orientation",
        "pattern_is_alt",
        "BA_source",
        "pair_posterior",
        "source",
    ]

    combined_values[output_columns].to_csv(
        MATRIX_OUTPUT,
        sep="\t",
        index=False,
        compression="gzip",
    )

    focal_rows = combined_values.loc[
        combined_values["plot_group"].astype(str).isin(
            ["TN_Left", "TN_Right"]
        )
    ]

    if focal_rows["mean_pattern_dosage"].isna().any():
        raise RuntimeError(
            "White-space guard failed: a selected focal cell is missing"
        )

    audit = {
        "selection": {
            "minimum_comparative_BA": ">0.99",
            "both_recalled_pair_posterior": ">=0.95",
            "recalled_genotype_state_used_for_selection": False,
            "missing_recalled_sites_retained": False,
        },
        "BA_sources": {
            "Band": str(BAND_BA_PATH.resolve()),
            "Dennis_Ray": str(DR_BA_PATH.resolve()),
        },
        "population_display": (
            "equal-individual mean BA-oriented pattern-associated allele dosage"
        ),
        "focal_display": (
            "reciprocal-mixture recalled genotype dosage; both sides required"
        ),
        "phase_required": False,
        "panels": details,
        "combined_pdf": str(COMBINED_PDF.resolve()),
        "selected_sites": str(SITE_OUTPUT.resolve()),
        "values": str(MATRIX_OUTPUT.resolve()),
    }

    AUDIT_OUTPUT.write_text(
        json.dumps(audit, indent=2) + "\n"
    )

    print()
    print(f"COMBINED_PDF\t{COMBINED_PDF.resolve()}")
    print(f"SELECTED_SITES\t{SITE_OUTPUT.resolve()}")
    print(f"AUDIT\t{AUDIT_OUTPUT.resolve()}")


if __name__ == "__main__":
    main()
