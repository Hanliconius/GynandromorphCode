#!/usr/bin/env python3

import csv
import glob
import gzip
import itertools
import json
import os
import random
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap


PROJECT = Path(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/"
    "GynandroCalling"
)

VALIDATION_ROOT = (
    PROJECT
    / "analysis/genomewide_mixing_corrected_v1/"
    "full_sib_validation"
)

IBD_ROOT = (
    VALIDATION_ROOT
    / "ibd_transmission"
)

INPUT_ROOT = (
    IBD_ROOT
    / "by_chromosome"
)

OUTPUT_ROOT = (
    IBD_ROOT
    / "genomewide"
)

FIGURE_ROOT = (
    IBD_ROOT
    / "figures"
)

MENDELIAN_ROOT = (
    VALIDATION_ROOT
    / "mendelian_qc/genomewide"
)

CALL_ROOT = (
    VALIDATION_ROOT
    / "variant_calls/by_chromosome"
)

SAMPLE_FILE = (
    VALIDATION_ROOT
    / "metadata/bam_set/full_sib_14_samples.samples"
)

BCFTOOLS = Path(
    "/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/"
    ".local/share/mamba/envs/colias_w_popgen/bin/bcftools"
)

MAX_AF_POINTS_PER_PAIR = 500_000
RANDOM_SEED = 199
MINIMUM_DEPTH = 5
MINIMUM_PL_GAP = 20

BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
PURPLE = "#CC79A7"
YELLOW = "#E69F00"
GREY = "#6B7280"
LIGHT_GREY = "#E5E7EB"
DARK = "#1F2937"

IBD_COLORS = {
    0: ORANGE,
    1: BLUE,
    2: GREEN,
}

HEX_CMAP = LinearSegmentedColormap.from_list(
    "clean_blue_hex",
    ["#F8FAFC", "#D7EAF5", "#74ADD1", "#0072B2", "#003B63"],
)


def configure_plotting():
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "axes.edgecolor": DARK,
            "axes.linewidth": 0.8,
            "axes.facecolor": "white",
            "figure.facecolor": "white",
            "xtick.color": DARK,
            "ytick.color": DARK,
            "text.color": DARK,
            "axes.titlecolor": DARK,
            "axes.labelcolor": DARK,
            "legend.frameon": False,
            "savefig.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def panel_label(axis, label):
    axis.text(
        -0.12,
        1.08,
        label,
        transform=axis.transAxes,
        fontsize=14,
        fontweight="bold",
        va="top",
        ha="left",
    )


def read_tsv(path):
    return pd.read_csv(path, sep="\t")


def write_tsv_atomic(table, path):
    temporary = Path(str(path) + ".tmp")
    table.to_csv(temporary, sep="\t", index=False)
    os.replace(temporary, path)


def write_gzip_tsv_atomic(table, path):
    temporary = Path(str(path) + ".tmp")

    with gzip.open(temporary, "wt", newline="") as handle:
        table.to_csv(handle, sep="\t", index=False)

    os.replace(temporary, path)


def parse_sample_field(field):
    parts = field.split(":", 3)

    while len(parts) < 4:
        parts.append(".")

    gt_text, dp_text, ad_text, pl_text = parts

    if gt_text in {".", "./.", ".|."}:
        return None

    alleles = gt_text.replace("|", "/").split("/")

    if len(alleles) != 2:
        return None

    try:
        allele_numbers = [int(value) for value in alleles]
        depth = int(dp_text)
        allele_depths = [
            int(value)
            for value in ad_text.split(",")
        ]
        likelihoods = [
            int(value)
            for value in pl_text.split(",")
        ]
    except ValueError:
        return None

    if any(value not in {0, 1} for value in allele_numbers):
        return None

    if len(allele_depths) != 2 or len(likelihoods) < 3:
        return None

    dosage = sum(allele_numbers)
    ordered_pl = sorted(
        (value, index)
        for index, value in enumerate(likelihoods[:3])
    )

    if (
        depth < MINIMUM_DEPTH
        or sum(allele_depths) <= 0
        or ordered_pl[0][1] != dosage
        or ordered_pl[1][0] - ordered_pl[0][0]
            < MINIMUM_PL_GAP
    ):
        return None

    return allele_depths[1] / sum(allele_depths)


def select_representative_pairs(pair_qc):
    eligible = pair_qc.loc[
        pair_qc["evaluated_sites"]
        >= pair_qc["evaluated_sites"].median()
    ].copy()

    if len(eligible) < 3:
        eligible = pair_qc.copy()

    eligible = eligible.sort_values(
        [
            "opposing_homozygote_rate",
            "evaluated_sites",
        ],
        ascending=[True, False],
    ).reset_index(drop=True)

    target_indices = sorted(
        {
            0,
            len(eligible) // 2,
            len(eligible) - 1,
        }
    )

    representatives = eligible.iloc[target_indices].copy()

    representatives["selection_label"] = [
        "lower genotype contrast",
        "median genotype contrast",
        "higher genotype contrast",
    ][:len(representatives)]

    return representatives


def reservoir_add(storage, seen, value, rng):
    seen += 1

    if len(storage) < MAX_AF_POINTS_PER_PAIR:
        storage.append(value)
    else:
        replacement = rng.randrange(seen)

        if replacement < MAX_AF_POINTS_PER_PAIR:
            storage[replacement] = value

    return seen


def collect_af_data(representatives):
    samples = [
        line.strip()
        for line in SAMPLE_FILE.read_text().splitlines()
        if line.strip()
    ]

    selected_pairs = [
        (row.offspring_a, row.offspring_b)
        for row in representatives.itertuples()
    ]

    selected_sample_set = {
        sample
        for pair in selected_pairs
        for sample in pair
    }

    selected_samples = [
        sample
        for sample in samples
        if sample in selected_sample_set
    ]

    selected_index = {
        sample: index
        for index, sample in enumerate(selected_samples)
    }

    pair_points = {
        pair: []
        for pair in selected_pairs
    }

    pair_seen = {
        pair: 0
        for pair in selected_pairs
    }

    rng = random.Random(RANDOM_SEED)

    for chromosome in range(1, 21):
        pattern = str(
            CALL_ROOT
            / f"chr{chromosome}_*"
            / "*.bcf"
        )

        candidates = [
            Path(path)
            for path in glob.glob(pattern)
            if not path.endswith(".bcf.new")
        ]

        if len(candidates) != 1:
            raise RuntimeError(
                "Expected one family BCF for chromosome "
                f"{chromosome}; found {candidates}"
            )

        bcf = candidates[0]

        command = [
            str(BCFTOOLS),
            "query",
            "-s",
            ",".join(selected_samples),
            "-f",
            "%REF\\t%ALT[\\t%GT:%DP:%AD:%PL]\\n",
            str(bcf),
        ]

        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1024 * 1024,
        )

        assert process.stdout is not None

        for line in process.stdout:
            fields = line.rstrip("\n").split("\t")

            if len(fields) != 2 + len(selected_samples):
                continue

            ref = fields[0]
            alt = fields[1]

            if (
                len(ref) != 1
                or len(alt) != 1
                or ref not in "ACGT"
                or alt not in "ACGT"
                or "," in alt
            ):
                continue

            af_values = [
                parse_sample_field(value)
                for value in fields[2:]
            ]

            for pair in selected_pairs:
                first_af = af_values[
                    selected_index[pair[0]]
                ]
                second_af = af_values[
                    selected_index[pair[1]]
                ]

                if first_af is None or second_af is None:
                    continue

                pair_seen[pair] = reservoir_add(
                    pair_points[pair],
                    pair_seen[pair],
                    (first_af, second_af),
                    rng,
                )

        stderr_text = ""

        if process.stderr is not None:
            stderr_text = process.stderr.read()

        return_code = process.wait()

        if return_code != 0:
            if stderr_text:
                print(stderr_text, file=sys.stderr)

            raise RuntimeError(
                "bcftools AF extraction failed for "
                f"chromosome {chromosome}"
            )

        print(
            "AF_SCAN_COMPLETE"
            f"\tchromosome={chromosome}"
        )

    return pair_points, pair_seen


def summarize_ibd(windows):
    rows = []

    for pair, group in windows.groupby(
        ["offspring_a", "offspring_b"],
        sort=True,
    ):
        resolved = pd.to_numeric(
            group["provisional_IBD"],
            errors="coerce",
        ).dropna()

        counts = resolved.astype(int).value_counts()

        total = len(resolved)

        rows.append(
            {
                "offspring_a": pair[0],
                "offspring_b": pair[1],
                "total_windows": len(group),
                "resolved_windows": total,
                "unresolved_windows":
                    len(group) - total,
                "IBD0_windows": int(counts.get(0, 0)),
                "IBD1_windows": int(counts.get(1, 0)),
                "IBD2_windows": int(counts.get(2, 0)),
                "IBD0_fraction":
                    counts.get(0, 0) / total
                    if total else np.nan,
                "IBD1_fraction":
                    counts.get(1, 0) / total
                    if total else np.nan,
                "IBD2_fraction":
                    counts.get(2, 0) / total
                    if total else np.nan,
            }
        )

    return pd.DataFrame(rows)


def transition_summary(windows):
    rows = []

    for (
        offspring_a,
        offspring_b,
        chromosome,
    ), group in windows.groupby(
        [
            "offspring_a",
            "offspring_b",
            "chromosome",
        ],
        sort=True,
    ):
        group = group.sort_values("window_start")

        maternal = [
            state
            for state in group["maternal_state"]
            if state in {"shared", "not_shared"}
        ]

        paternal = [
            state
            for state in group["paternal_state"]
            if state in {"shared", "not_shared"}
        ]

        maternal_transitions = sum(
            first != second
            for first, second in zip(
                maternal,
                maternal[1:],
            )
        )

        paternal_transitions = sum(
            first != second
            for first, second in zip(
                paternal,
                paternal[1:],
            )
        )

        maternal_counts = Counter(maternal)
        maternal_resolved = len(maternal)

        maternal_dominant_fraction = (
            max(maternal_counts.values())
            / maternal_resolved
            if maternal_resolved else np.nan
        )

        rows.append(
            {
                "offspring_a": offspring_a,
                "offspring_b": offspring_b,
                "chromosome": chromosome,
                "maternal_resolved_windows":
                    maternal_resolved,
                "maternal_transitions":
                    maternal_transitions,
                "maternal_dominant_state_fraction":
                    maternal_dominant_fraction,
                "paternal_resolved_windows":
                    len(paternal),
                "paternal_transitions":
                    paternal_transitions,
            }
        )

    return pd.DataFrame(rows)


def make_qc_page(pdf, offspring_qc, pair_qc):
    figure, axes = plt.subplots(
        1,
        2,
        figsize=(12, 5.8),
        constrained_layout=True,
    )

    ordered = offspring_qc.sort_values(
        "mendelian_error_rate"
    )

    colors = [
        ORANGE if flag == "REVIEW" else BLUE
        for flag in ordered["review_flag"]
    ]

    axes[0].barh(
        ordered["offspring"],
        ordered["mendelian_error_rate"] * 100,
        color=colors,
        edgecolor="white",
    )
    axes[0].axvline(
        ordered["review_threshold"].iloc[0] * 100,
        color=ORANGE,
        linestyle="--",
        linewidth=1.4,
        label="Diagnostic review threshold",
    )
    axes[0].set_xlabel("High-confidence Mendelian errors (%)")
    axes[0].set_title("Pedigree consistency by offspring")
    axes[0].legend(loc="lower right")
    axes[0].grid(axis="x", color=LIGHT_GREY, linewidth=0.7)
    panel_label(axes[0], "A")

    scatter = axes[1].scatter(
        pair_qc["evaluated_sites"] / 1_000_000,
        pair_qc["opposing_homozygote_rate"] * 100,
        c=pair_qc["combined_information_score"],
        cmap=HEX_CMAP,
        s=55,
        edgecolor="white",
        linewidth=0.5,
    )
    axes[1].set_xlabel("Evaluated sites (millions)")
    axes[1].set_ylabel("Opposing homozygotes (%)")
    axes[1].set_title("Information available from 66 sibling pairs")
    axes[1].grid(color=LIGHT_GREY, linewidth=0.7)
    panel_label(axes[1], "B")

    colorbar = figure.colorbar(
        scatter,
        ax=axes[1],
        fraction=0.046,
        pad=0.04,
    )
    colorbar.set_label("Coverage × contrast score")

    figure.suptitle(
        "Step 1 — Confirm pedigree quality and empirical genotype contrast",
        fontsize=15,
        fontweight="bold",
    )

    pdf.savefig(figure, bbox_inches="tight")
    plt.close(figure)


def make_af_page(pdf, representatives, pair_points, pair_seen):
    number = len(representatives)

    figure, axes = plt.subplots(
        1,
        number,
        figsize=(5.0 * number, 5.3),
        constrained_layout=True,
        squeeze=False,
    )
    axes = axes[0]

    centroids = np.array(
        list(itertools.product([0.0, 0.5, 1.0], repeat=2))
    )

    for index, row in enumerate(
        representatives.itertuples()
    ):
        axis = axes[index]
        pair = (row.offspring_a, row.offspring_b)
        points = np.asarray(pair_points[pair])

        if len(points) == 0:
            raise RuntimeError(
                f"No AF points retained for pair {pair}"
            )

        hexbin = axis.hexbin(
            points[:, 0],
            points[:, 1],
            gridsize=65,
            extent=(0, 1, 0, 1),
            mincnt=1,
            bins="log",
            cmap=HEX_CMAP,
            linewidths=0,
        )

        axis.scatter(
            centroids[:, 0],
            centroids[:, 1],
            marker="o",
            s=46,
            facecolor="white",
            edgecolor=ORANGE,
            linewidth=1.3,
            zorder=4,
            label="Nine diploid centroids",
        )

        axis.plot(
            [0, 1],
            [0, 1],
            color=GREY,
            linestyle="--",
            linewidth=0.8,
            alpha=0.7,
        )

        axis.set_xlim(-0.02, 1.02)
        axis.set_ylim(-0.02, 1.02)
        axis.set_aspect("equal")
        axis.set_xlabel(f"Alternate-allele fraction\n{pair[0]}")
        axis.set_ylabel(f"Alternate-allele fraction\n{pair[1]}")
        axis.set_title(
            f"{row.selection_label}\n"
            f"{pair[0]} × {pair[1]}"
        )
        axis.grid(color=LIGHT_GREY, linewidth=0.5)
        panel_label(axis, chr(ord("A") + index))

        axis.text(
            0.03,
            0.97,
            (
                f"Plotted: {len(points):,}\n"
                f"Eligible: {pair_seen[pair]:,}"
            ),
            transform=axis.transAxes,
            va="top",
            ha="left",
            fontsize=8.5,
            bbox={
                "facecolor": "white",
                "edgecolor": LIGHT_GREY,
                "alpha": 0.9,
            },
        )

        if index == number - 1:
            axis.legend(loc="lower right", fontsize=8)

        figure.colorbar(
            hexbin,
            ax=axis,
            fraction=0.045,
            pad=0.03,
            label="log10 hexagon count",
        )

    figure.suptitle(
        "Step 2 — Empirical full siblings occupy the nine-genotype space",
        fontsize=15,
        fontweight="bold",
    )

    pdf.savefig(figure, bbox_inches="tight")
    plt.close(figure)


def make_composition_page(pdf, ibd_summary):
    ordered = ibd_summary.sort_values(
        ["IBD0_fraction", "IBD2_fraction"]
    ).reset_index(drop=True)

    x = np.arange(len(ordered))

    figure, axes = plt.subplots(
        1,
        2,
        figsize=(12, 5.8),
        constrained_layout=True,
    )

    bottom = np.zeros(len(ordered))

    for state in [0, 1, 2]:
        values = ordered[f"IBD{state}_fraction"].fillna(0)
        axes[0].bar(
            x,
            values,
            bottom=bottom,
            width=0.9,
            color=IBD_COLORS[state],
            label=f"IBD{state}",
        )
        bottom += values.to_numpy()

    axes[0].set_xlim(-1, len(ordered))
    axes[0].set_ylim(0, 1)
    axes[0].set_xlabel("Sibling pairs, ordered by IBD composition")
    axes[0].set_ylabel("Fraction of resolved 1-Mb windows")
    axes[0].set_title("Genome-wide sharing mosaic")
    axes[0].legend(ncol=3, loc="upper center")
    axes[0].grid(axis="y", color=LIGHT_GREY, linewidth=0.7)
    panel_label(axes[0], "A")

    observed = [
        ordered[f"IBD{state}_fraction"].mean()
        for state in [0, 1, 2]
    ]
    expected = [0.25, 0.50, 0.25]

    bar_x = np.arange(3)
    width = 0.36

    axes[1].bar(
        bar_x - width / 2,
        expected,
        width,
        color=LIGHT_GREY,
        edgecolor=GREY,
        label="Full-sib expectation",
    )
    axes[1].bar(
        bar_x + width / 2,
        observed,
        width,
        color=[IBD_COLORS[state] for state in [0, 1, 2]],
        label="Observed window mean",
    )
    axes[1].set_xticks(bar_x)
    axes[1].set_xticklabels(["IBD0", "IBD1", "IBD2"])
    axes[1].set_ylim(0, max(0.65, max(observed) + 0.08))
    axes[1].set_ylabel("Fraction")
    axes[1].set_title("Observed versus idealized full-sib expectation")
    axes[1].legend()
    axes[1].grid(axis="y", color=LIGHT_GREY, linewidth=0.7)
    panel_label(axes[1], "B")

    figure.suptitle(
        "Step 3 — Reconstruct the full-sibling IBD state space",
        fontsize=15,
        fontweight="bold",
    )

    pdf.savefig(figure, bbox_inches="tight")
    plt.close(figure)


def make_transmission_page(pdf, transitions):
    figure, axes = plt.subplots(
        1,
        2,
        figsize=(12, 5.8),
        constrained_layout=True,
    )

    valid = transitions[
        transitions["maternal_resolved_windows"] > 1
    ]

    axes[0].hist(
        valid["maternal_dominant_state_fraction"],
        bins=np.linspace(0.5, 1.0, 21),
        color=PURPLE,
        edgecolor="white",
    )
    axes[0].axvline(
        1.0,
        color=DARK,
        linestyle="--",
        linewidth=1.2,
    )
    axes[0].set_xlabel(
        "Fraction of resolved windows in dominant maternal state"
    )
    axes[0].set_ylabel("Sibling-pair × chromosome combinations")
    axes[0].set_title("Maternal chromosome-wide consistency")
    axes[0].grid(axis="y", color=LIGHT_GREY, linewidth=0.7)
    panel_label(axes[0], "A")

    maximum = max(
        transitions["maternal_transitions"].max(),
        transitions["paternal_transitions"].max(),
    )

    bins = np.arange(-0.5, maximum + 1.5, 1)

    axes[1].hist(
        transitions["maternal_transitions"],
        bins=bins,
        color=PURPLE,
        alpha=0.75,
        label="Maternal",
    )
    axes[1].hist(
        transitions["paternal_transitions"],
        bins=bins,
        histtype="step",
        color=BLUE,
        linewidth=2,
        label="Paternal",
    )
    axes[1].set_xlabel("Resolved sharing-state transitions")
    axes[1].set_ylabel("Sibling-pair × chromosome combinations")
    axes[1].set_title("Maternal achiasmy versus paternal recombination")
    axes[1].legend()
    axes[1].grid(axis="y", color=LIGHT_GREY, linewidth=0.7)
    panel_label(axes[1], "B")

    figure.suptitle(
        "Step 4 — Test whether transmission patterns behave biologically",
        fontsize=15,
        fontweight="bold",
    )

    pdf.savefig(figure, bbox_inches="tight")
    plt.close(figure)


def add_cumulative_coordinates(windows):
    chromosome_max = (
        windows.groupby("chromosome")["window_end"]
        .max()
        .sort_index()
    )

    offsets = {}
    current = 0

    for chromosome, length in chromosome_max.items():
        offsets[int(chromosome)] = current
        current += int(length)

    windows = windows.copy()
    windows["genome_position"] = (
        windows.apply(
            lambda row: (
                offsets[int(row["chromosome"])]
                + (
                    float(row["window_start"])
                    + float(row["window_end"])
                )
                / 2
            ),
            axis=1,
        )
    )

    centers = {
        chromosome: offset
        + chromosome_max.loc[chromosome] / 2
        for chromosome, offset in offsets.items()
    }

    return windows, offsets, centers


def make_tracks_page(pdf, windows, representatives):
    figure, axes = plt.subplots(
        len(representatives),
        1,
        figsize=(14, 7.5),
        constrained_layout=True,
        sharex=True,
    )

    if len(representatives) == 1:
        axes = [axes]

    for index, row in enumerate(
        representatives.itertuples()
    ):
        axis = axes[index]

        selected = windows.loc[
            (
                windows["offspring_a"]
                == row.offspring_a
            )
            & (
                windows["offspring_b"]
                == row.offspring_b
            )
        ].copy()

        selected, offsets, centers = (
            add_cumulative_coordinates(selected)
        )

        resolved = selected.dropna(
            subset=["provisional_IBD"]
        ).copy()

        resolved["provisional_IBD"] = (
            resolved["provisional_IBD"].astype(int)
        )

        for state in [0, 1, 2]:
            state_rows = resolved.loc[
                resolved["provisional_IBD"] == state
            ]

            axis.scatter(
                state_rows["genome_position"],
                state_rows["provisional_IBD"],
                s=13,
                color=IBD_COLORS[state],
                label=f"IBD{state}",
                alpha=0.9,
            )

        for offset in list(offsets.values())[1:]:
            axis.axvline(
                offset,
                color=LIGHT_GREY,
                linewidth=0.6,
            )

        axis.set_ylim(-0.35, 2.35)
        axis.set_yticks([0, 1, 2])
        axis.set_ylabel("IBD state")
        axis.set_title(
            f"{row.selection_label}: "
            f"{row.offspring_a} × {row.offspring_b}",
            loc="left",
        )
        axis.grid(axis="y", color=LIGHT_GREY, linewidth=0.6)
        panel_label(axis, chr(ord("A") + index))

        if index == 0:
            axis.legend(
                ncol=3,
                loc="upper right",
            )

    final_selected = windows.loc[
        (
            windows["offspring_a"]
            == representatives.iloc[0]["offspring_a"]
        )
        & (
            windows["offspring_b"]
            == representatives.iloc[0]["offspring_b"]
        )
    ].copy()

    _, _, centers = add_cumulative_coordinates(final_selected)

    axes[-1].set_xticks(list(centers.values()))
    axes[-1].set_xticklabels(
        [str(chromosome) for chromosome in centers],
        fontsize=8,
    )
    axes[-1].set_xlabel("Autosome")

    figure.suptitle(
        "Step 5 — Genome-wide IBD mosaics reveal regional variation",
        fontsize=15,
        fontweight="bold",
    )

    pdf.savefig(figure, bbox_inches="tight")
    plt.close(figure)


def make_summary_figure(
    output_path,
    offspring_qc,
    ibd_summary,
    transitions,
):
    figure, axes = plt.subplots(
        1,
        3,
        figsize=(14, 4.6),
        constrained_layout=True,
    )

    ordered = offspring_qc.sort_values(
        "mendelian_error_rate"
    )

    axes[0].barh(
        ordered["offspring"],
        ordered["mendelian_error_rate"] * 100,
        color=BLUE,
    )
    axes[0].set_xlabel("Mendelian errors (%)")
    axes[0].set_title("Pedigree consistency")
    panel_label(axes[0], "A")

    means = [
        ibd_summary[f"IBD{state}_fraction"].mean()
        for state in [0, 1, 2]
    ]

    axes[1].bar(
        ["IBD0", "IBD1", "IBD2"],
        means,
        color=[ORANGE, BLUE, GREEN],
    )
    axes[1].plot(
        [0, 1, 2],
        [0.25, 0.50, 0.25],
        color=DARK,
        marker="o",
        linestyle="--",
        label="0.25 : 0.50 : 0.25",
    )
    axes[1].set_ylabel("Mean fraction of resolved windows")
    axes[1].set_title("Full-sibling IBD space")
    axes[1].legend(fontsize=8)
    panel_label(axes[1], "B")

    axes[2].scatter(
        transitions["maternal_transitions"],
        transitions["paternal_transitions"],
        s=13,
        alpha=0.45,
        color=PURPLE,
        edgecolor="none",
    )
    axes[2].set_xlabel("Maternal transitions")
    axes[2].set_ylabel("Paternal transitions")
    axes[2].set_title("Achiasmy/recombination diagnostic")
    panel_label(axes[2], "C")

    for axis in axes:
        axis.grid(color=LIGHT_GREY, linewidth=0.6)

    figure.suptitle(
        "Empirical full-sibling validation",
        fontsize=15,
        fontweight="bold",
    )

    figure.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(figure)


def main():
    started = time.time()
    configure_plotting()

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)

    completion_files = glob.glob(
        str(
            INPUT_ROOT
            / "chr*_*"
            / "sibling_transmission_IBD.complete.tsv"
        )
    )
    window_files = glob.glob(
        str(
            INPUT_ROOT
            / "chr*_*"
            / "sibling_transmission_IBD_windows.tsv"
        )
    )
    chromosome_files = glob.glob(
        str(
            INPUT_ROOT
            / "chr*_*"
            / "sibling_transmission_IBD_chromosome.tsv"
        )
    )

    if (
        len(completion_files) != 20
        or len(window_files) != 20
        or len(chromosome_files) != 20
    ):
        raise RuntimeError(
            "Expected 20 complete chromosome result sets; "
            f"found completion={len(completion_files)}, "
            f"windows={len(window_files)}, "
            f"chromosomes={len(chromosome_files)}"
        )

    completion = pd.concat(
        [read_tsv(path) for path in completion_files],
        ignore_index=True,
    )

    if (
        set(completion["chromosome"]) != set(range(1, 21))
        or set(completion["status"]) != {"PASS"}
    ):
        raise RuntimeError(
            "Stage 198 chromosome completion validation failed."
        )

    windows = pd.concat(
        [read_tsv(path) for path in window_files],
        ignore_index=True,
    )

    chromosomes = pd.concat(
        [read_tsv(path) for path in chromosome_files],
        ignore_index=True,
    )

    offspring_qc = read_tsv(
        MENDELIAN_ROOT
        / "genomewide_offspring_mendelian_qc.tsv"
    )

    pair_qc = read_tsv(
        MENDELIAN_ROOT
        / "genomewide_sibling_pair_genotype_states.tsv"
    )

    representatives = select_representative_pairs(
        pair_qc
    )

    ibd_summary = summarize_ibd(windows)
    transitions = transition_summary(windows)

    pair_points, pair_seen = collect_af_data(
        representatives
    )

    window_output = (
        OUTPUT_ROOT
        / "all_sibling_IBD_windows.tsv.gz"
    )
    chromosome_output = (
        OUTPUT_ROOT
        / "all_sibling_IBD_chromosomes.tsv"
    )
    summary_output = (
        OUTPUT_ROOT
        / "sibling_pair_IBD_summary.tsv"
    )
    transition_output = (
        OUTPUT_ROOT
        / "sibling_transmission_transition_summary.tsv"
    )
    representative_output = (
        OUTPUT_ROOT
        / "representative_sibling_pairs.tsv"
    )

    write_gzip_tsv_atomic(windows, window_output)
    write_tsv_atomic(chromosomes, chromosome_output)
    write_tsv_atomic(ibd_summary, summary_output)
    write_tsv_atomic(transitions, transition_output)

    representative_table = representatives[
        [
            "selection_label",
            "offspring_a",
            "offspring_b",
            "evaluated_sites",
            "opposing_homozygote_rate",
            "combined_information_score",
        ]
    ].copy()

    representative_table["AF_eligible_sites"] = [
        pair_seen[
            (row.offspring_a, row.offspring_b)
        ]
        for row in representative_table.itertuples()
    ]

    representative_table["AF_plotted_sites"] = [
        len(
            pair_points[
                (row.offspring_a, row.offspring_b)
            ]
        )
        for row in representative_table.itertuples()
    ]

    write_tsv_atomic(
        representative_table,
        representative_output,
    )

    narrative_pdf = (
        FIGURE_ROOT
        / "full_sib_IBD_narrative.pdf"
    )

    with PdfPages(narrative_pdf) as pdf:
        make_qc_page(
            pdf,
            offspring_qc,
            pair_qc,
        )
        make_af_page(
            pdf,
            representatives,
            pair_points,
            pair_seen,
        )
        make_composition_page(
            pdf,
            ibd_summary,
        )
        make_transmission_page(
            pdf,
            transitions,
        )
        make_tracks_page(
            pdf,
            windows,
            representatives,
        )

    summary_png = (
        FIGURE_ROOT
        / "full_sib_IBD_summary.png"
    )

    make_summary_figure(
        summary_png,
        offspring_qc,
        ibd_summary,
        transitions,
    )

    mean_ibd = {
        f"IBD{state}":
            float(
                ibd_summary[
                    f"IBD{state}_fraction"
                ].mean()
            )
        for state in [0, 1, 2]
    }

    audit = {
        "chromosomes": 20,
        "offspring": int(len(offspring_qc)),
        "sibling_pairs": int(len(ibd_summary)),
        "IBD_window_rows": int(len(windows)),
        "chromosome_pair_rows":
            int(len(chromosomes)),
        "representative_pairs":
            representative_table[
                [
                    "selection_label",
                    "offspring_a",
                    "offspring_b",
                    "AF_eligible_sites",
                    "AF_plotted_sites",
                ]
            ].to_dict(orient="records"),
        "mean_resolved_IBD_fractions":
            mean_ibd,
        "AF_plot_type": "hexbin",
        "maximum_AF_points_per_pair":
            MAX_AF_POINTS_PER_PAIR,
        "nine_diploid_centroids_overlaid":
            True,
        "automatic_exclusions": False,
        "founder_model_changed": False,
        "IBD_classifications_smoothed":
            False,
        "narrative_pdf": str(narrative_pdf),
        "summary_png": str(summary_png),
        "elapsed_seconds":
            round(time.time() - started, 3),
        "next_stage": (
            "fit chromosome-aware IBD smoothing and "
            "construct empirical full-sibling mixing simulations"
        ),
    }

    audit_output = (
        FIGURE_ROOT
        / "full_sib_IBD_narrative.audit.json"
    )
    temporary_audit = Path(str(audit_output) + ".tmp")
    temporary_audit.write_text(
        json.dumps(audit, indent=2) + "\n"
    )
    os.replace(temporary_audit, audit_output)

    completion_output = (
        OUTPUT_ROOT
        / "full_sib_IBD_collection.complete.tsv"
    )

    completion_table = pd.DataFrame(
        [
            {
                "status": "PASS",
                "chromosomes": 20,
                "offspring": len(offspring_qc),
                "sibling_pairs": len(ibd_summary),
                "IBD_window_rows": len(windows),
                "narrative_pdf": narrative_pdf,
                "elapsed_seconds":
                    round(time.time() - started, 3),
            }
        ]
    )

    write_tsv_atomic(
        completion_table,
        completion_output,
    )

    print("FULL-SIB IBD NARRATIVE COMPLETE")
    print(
        representative_table.to_string(index=False)
    )
    print(f"NARRATIVE_PDF\t{narrative_pdf}")
    print(f"SUMMARY_PNG\t{summary_png}")
    print(
        "ENDPOINT="
        "FULL_SIB_IBD_DIAGNOSTICS_READY"
    )


if __name__ == "__main__":
    main()
