#!/usr/bin/env python3
"""Figure 3 trial: resurrect the original stage-117a informative-site view.

Uses the *existing* stage-117a informative SNP table directly rather than
intersecting it with the later Krista five-locus tables.

Site selection is therefore exactly the original first-pass rule:
  |delta AF| >= 0.80
  both comparative groups >= 0.80 callable
  DP >= 2 and GQ >= 20 for the calls used in the comparison
  TN_Left/TN_Right excluded from site selection and polarization

Population rows are plotted from the stage-117a group allele frequencies.
TN rows use only the reciprocal-mixture posterior>=0.95 recalled genotypes.
No phasing is used.
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

INFORMATIVE = (
    PROJECT
    / "analysis/original_502/optix_deltaAF_haplotype_analysis"
    / "optix_B_D_R_deltaAF_informative_SNPs.tsv"
)

LR = (
    PROJECT
    / "analysis/original_502/reciprocal_L_R_mixture"
    / "optix_reciprocal_L_R_genotype_pairs.posterior95.tsv"
)

OUTDIR = (
    PROJECT
    / "analysis/figure3/optix_first_informative_sites"
)

COMBINED = OUTDIR / "optix_first_informative_B_D_R_collapsed_plus_recalled_LR.pdf"
VALUES = OUTDIR / "optix_first_informative_plot_values.tsv.gz"
AUDIT = OUTDIR / "optix_first_informative.audit.json"

MODULES = ("Band", "Ray", "Dennis")

ROW_ORDER = {
    "Band": ("absent", "TN_Left", "TN_Right", "present"),
    "Ray": ("present", "TN_Left", "TN_Right", "absent"),
    "Dennis": ("present", "TN_Left", "TN_Right", "absent"),
}


def as_bool(value):
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def canonical_gt(value):
    if value is None or pd.isna(value):
        return None
    text = str(value).strip().replace("|", "/")
    if text == "1/0":
        text = "0/1"
    return text if text in {"0/0", "0/1", "1/1"} else None


def alt_dosage(gt):
    gt = canonical_gt(gt)
    return {
        "0/0": 0.0,
        "0/1": 0.5,
        "1/1": 1.0,
    }.get(gt, np.nan)


def pattern_dosage_from_gt(gt, pattern_is_alt):
    dose = alt_dosage(gt)
    if not np.isfinite(dose):
        return np.nan
    return dose if pattern_is_alt else 1.0 - dose


def load_sites():
    table = pd.read_csv(INFORMATIVE, sep="\t", low_memory=False)

    required = {
        "module", "chrom", "pos", "ref", "alt",
        "positive_label", "negative_label",
        "positive_n", "negative_n",
        "positive_call_rate", "negative_call_rate",
        "positive_ALT_AF", "negative_ALT_AF",
        "abs_delta_AF",
    }
    missing = required - set(table.columns)
    if missing:
        raise RuntimeError(
            "Informative-site table missing columns: "
            + ", ".join(sorted(missing))
        )

    table["pos"] = pd.to_numeric(table["pos"], errors="raise").astype(int)

    numeric = [
        "positive_n", "negative_n",
        "positive_call_rate", "negative_call_rate",
        "positive_ALT_AF", "negative_ALT_AF",
        "abs_delta_AF",
    ]
    for column in numeric:
        table[column] = pd.to_numeric(table[column], errors="raise")

    table = table.loc[table["module"].isin(MODULES)].copy()

    bad = table.loc[
        (table["abs_delta_AF"] < 0.80)
        | (table["positive_call_rate"] < 0.80)
        | (table["negative_call_rate"] < 0.80)
    ]
    if not bad.empty:
        raise RuntimeError(
            "The stored stage-117a table contains rows outside "
            "the original informative-site thresholds."
        )

    keys = ["module", "chrom", "pos", "ref", "alt"]
    if table.duplicated(keys).any():
        raise RuntimeError("Duplicate module/site rows in informative table")

    # In stage 117a, the positive group is always the focal pattern-present
    # group (Band present, Dennis present, or Ray present). Therefore the
    # allele enriched in positive vs negative is the pattern-associated allele.
    table["pattern_is_alt"] = (
        table["positive_ALT_AF"] > table["negative_ALT_AF"]
    )

    table["present_mean_dosage"] = np.where(
        table["pattern_is_alt"],
        table["positive_ALT_AF"],
        1.0 - table["positive_ALT_AF"],
    )
    table["absent_mean_dosage"] = np.where(
        table["pattern_is_alt"],
        table["negative_ALT_AF"],
        1.0 - table["negative_ALT_AF"],
    )

    return table.sort_values(
        ["module", "pos", "ref", "alt"],
        kind="mergesort",
    ).reset_index(drop=True)


def load_recalled_pairs():
    table = pd.read_csv(LR, sep="\t", low_memory=False)

    required = {"chrom", "pos", "ref", "alt", "pair_posterior", "posterior95"}
    missing = required - set(table.columns)
    if missing:
        raise RuntimeError(
            "Reciprocal L/R table missing columns: "
            + ", ".join(sorted(missing))
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
            column for column in table.columns
            if "pair" in column.lower()
            and "original" not in column.lower()
            and column not in {"pair_posterior"}
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
        raise RuntimeError(
            f"Could not split reciprocal pair states in {pair_column}"
        )

    table["TN_Left_GT"] = parts[0].map(canonical_gt)
    table["TN_Right_GT"] = parts[1].map(canonical_gt)

    table["confident"] = (
        table["posterior95"]
        & table["pair_posterior"].ge(0.95)
        & table["TN_Left_GT"].notna()
        & table["TN_Right_GT"].notna()
    )

    keys = ["chrom", "pos", "ref", "alt"]
    if table.duplicated(keys).any():
        raise RuntimeError("Duplicate sites in reciprocal L/R table")

    return table[
        keys + ["pair_posterior", "confident", "TN_Left_GT", "TN_Right_GT"]
    ].copy()


def make_cmap():
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


def build_module(module, sites, recalled):
    d = sites.loc[sites["module"].eq(module)].copy()
    if d.empty:
        raise RuntimeError(f"No stage-117a informative sites for {module}")

    keys = ["chrom", "pos", "ref", "alt"]
    d = d.merge(recalled, on=keys, how="left", validate="one_to_one")

    d["TN_Left"] = np.nan
    d["TN_Right"] = np.nan

    confident = d["confident"].fillna(False).astype(bool)

    for sample in ("TN_Left", "TN_Right"):
        gt_column = f"{sample}_GT"
        d.loc[confident, sample] = [
            pattern_dosage_from_gt(gt, is_alt)
            for gt, is_alt in zip(
                d.loc[confident, gt_column],
                d.loc[confident, "pattern_is_alt"],
            )
        ]

    matrix_by_name = {
        "present": d["present_mean_dosage"].to_numpy(dtype=float),
        "absent": d["absent_mean_dosage"].to_numpy(dtype=float),
        "TN_Left": d["TN_Left"].to_numpy(dtype=float),
        "TN_Right": d["TN_Right"].to_numpy(dtype=float),
    }

    order = ROW_ORDER[module]
    matrix = np.vstack([matrix_by_name[name] for name in order])

    pos_n = sorted(d["positive_n"].dropna().astype(int).unique().tolist())
    neg_n = sorted(d["negative_n"].dropna().astype(int).unique().tolist())

    if len(pos_n) != 1 or len(neg_n) != 1:
        raise RuntimeError(
            f"{module}: comparative group sizes vary across selected sites: "
            f"positive={pos_n}, negative={neg_n}"
        )

    labels = {
        "present": f"{module} present (n={pos_n[0]})",
        "absent": f"{module} absent (n={neg_n[0]})",
        "TN_Left": "TN_Left recalled",
        "TN_Right": "TN_Right recalled",
    }

    width = max(7.5, min(15.0, 5.8 + len(d) / 18.0))
    fig, ax = plt.subplots(figsize=(width, 4.0), constrained_layout=True)

    im = ax.imshow(
        matrix,
        cmap=make_cmap(),
        vmin=0.0,
        vmax=1.0,
        aspect="auto",
        interpolation="nearest",
    )

    ax.axhline(0.5, color="#303030", linewidth=1.2)
    ax.axhline(2.5, color="#303030", linewidth=1.2)

    ax.set_yticks(np.arange(4))
    ax.set_yticklabels([labels[name] for name in order])

    tick_count = min(10, len(d))
    tick_indices = np.unique(
        np.linspace(0, len(d) - 1, tick_count).round().astype(int)
    )
    ax.set_xticks(tick_indices)
    ax.set_xticklabels(
        [f"{d.iloc[i]['pos'] / 1_000_000:.4f}" for i in tick_indices],
        rotation=45,
        ha="right",
    )
    ax.set_xlabel("Position on Hmr18 (Mb)")

    ax.text(
        0.0, 1.105,
        f"{module}: original informative sites",
        transform=ax.transAxes,
        ha="left", va="bottom",
        fontsize=14, fontweight="bold",
    )
    ax.text(
        0.0, 1.045,
        (
            "Stage 117a: |ΔAF|≥0.80, both comparative groups ≥80% callable; "
            "TN rows are independent posterior≥0.95 reciprocal-mixture recalls."
        ),
        transform=ax.transAxes,
        ha="left", va="bottom",
        fontsize=8.5, color="#555555",
    )

    cbar = fig.colorbar(im, ax=ax, fraction=0.028, pad=0.018)
    cbar.set_ticks([0.0, 0.5, 1.0])
    cbar.set_ticklabels(["0  A/A", "0.5  P/A", "1  P/P"])
    cbar.set_label("Pattern-associated allele dosage")

    for spine in ax.spines.values():
        spine.set_visible(False)

    return fig, d, pos_n[0], neg_n[0]


def main():
    for path in (INFORMATIVE, LR):
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing or empty input: {path}")

    OUTDIR.mkdir(parents=True, exist_ok=True)

    sites = load_sites()
    recalled = load_recalled_pairs()

    values = []
    details = {}

    with PdfPages(COMBINED) as combined:
        for module in MODULES:
            fig, d, present_n, absent_n = build_module(module, sites, recalled)

            panel_pdf = OUTDIR / f"optix_{module}_first_informative_only.pdf"
            fig.savefig(panel_pdf, bbox_inches="tight")
            combined.savefig(fig, bbox_inches="tight")
            plt.close(fig)

            values.append(d.assign(panel=module))
            details[module] = {
                "informative_sites": int(len(d)),
                "present_n": int(present_n),
                "absent_n": int(absent_n),
                "posterior95_LR_pairs": int(
                    d["confident"].fillna(False).astype(bool).sum()
                ),
                "position_min": int(d["pos"].min()),
                "position_max": int(d["pos"].max()),
                "pdf": str(panel_pdf.resolve()),
            }

            print(
                f"{module}\tinformative_sites={len(d)}"
                f"\trange={int(d['pos'].min())}-{int(d['pos'].max())}"
                f"\tposterior95_pairs={details[module]['posterior95_LR_pairs']}"
            )

    pd.concat(values, ignore_index=True).to_csv(
        VALUES,
        sep="\t",
        index=False,
        compression="gzip",
    )

    audit = {
        "selection_source": str(INFORMATIVE.resolve()),
        "selection": {
            "minimum_abs_delta_AF": 0.80,
            "minimum_group_call_rate": 0.80,
            "minimum_DP": 2,
            "minimum_GQ": 20,
            "TN_used_for_selection": False,
            "TN_used_for_polarization": False,
        },
        "display": {
            "population_rows": (
                "stage-117a comparative group pattern-allele frequencies"
            ),
            "TN_rows": (
                "reciprocal-mixture genotype calls, pair posterior>=0.95 only"
            ),
            "phase_required": False,
            "dosage_scale": {
                "0": "A/A",
                "0.5": "P/A",
                "1": "P/P",
            },
        },
        "modules": details,
        "combined_pdf": str(COMBINED.resolve()),
        "values": str(VALUES.resolve()),
    }

    AUDIT.write_text(json.dumps(audit, indent=2) + "\n")

    print(f"COMBINED_PDF\t{COMBINED.resolve()}")
    print(f"AUDIT\t{AUDIT.resolve()}")


if __name__ == "__main__":
    main()
