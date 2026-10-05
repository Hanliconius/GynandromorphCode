#!/usr/bin/env python3

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
import numpy as np


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

VCF_DIR = (
    PROJECT
    / "external/krista_genotype_plots_original/Genotype_Plots"
)

SCORE_FILE = (
    PROJECT
    / "analysis/original_502/checkpoints"
    / "optix_original_502_comparative_scores.tsv"
)

CANDIDATE_FILE = (
    PROJECT
    / "analysis/original_502/checkpoints"
    / "optix_comparative_informative_sites_BA80.tsv"
)

METADATA_FILE = (
    PROJECT
    / "data/derived/region_panel_metadata.tsv"
)

OUTPUT_DIR = (
    PROJECT
    / "analysis/figures/supplementary_aligned_genoplots"
)

CHECKPOINT_DIR = (
    PROJECT
    / "analysis/original_502/checkpoints"
)

PANELS = {
    "band": {
        "trait": "B",
        "title": "Band (B)",
        "vcf": "optix_band-1112651-1114314.vcf",
    },
    "dennis1": {
        "trait": "D",
        "title": "Dennis interval 1 (D)",
        "vcf": "dennis1-1215210-1216124.vcf",
    },
    "dennis2": {
        "trait": "D",
        "title": "Dennis interval 2 (D)",
        "vcf": "dennis2-1218557-1221494.vcf",
    },
    "ray1": {
        "trait": "R",
        "title": "Ray interval 1 (R)",
        "vcf": "ray1-1156746-1160220.vcf",
    },
    "ray2": {
        "trait": "R",
        "title": "Ray interval 2 (R)",
        "vcf": "ray2-1174547-1187833.vcf",
    },
}

FOCAL = {
    "TN2L_S1",
    "TN3L_S3",
    "TN2R_S2",
    "TN3R_S4",
    "TN_Left",
    "TN_Right",
}


def as_float(value):
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None

    if math.isfinite(result):
        return result

    return None


def normalize_gt(gt):
    if not gt or "." in gt:
        return None

    alleles = gt.replace("|", "/").split("/")

    if len(alleles) != 2:
        return None

    if any(allele not in {"0", "1"} for allele in alleles):
        return None

    return sum(int(allele) for allele in alleles)


def read_metadata():
    metadata = {}

    with METADATA_FILE.open(newline="") as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        for row in reader:
            sample = (
                row.get("Accession")
                or ""
            ).strip()

            if sample:
                metadata[sample] = row

    return metadata


def read_scores():
    scores = {}

    with SCORE_FILE.open(newline="") as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        for row in reader:
            panel = row["panel"]

            if panel not in PANELS:
                continue

            key = (
                panel,
                row["chrom"],
                int(row["pos"]),
                row["ref"],
                row["alt"],
                row["scheme"],
            )

            scores[key] = row

    return scores


def read_inclusive_candidates():
    candidates = set()

    if not CANDIDATE_FILE.is_file():
        return candidates

    with CANDIDATE_FILE.open(newline="") as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        for row in reader:
            candidates.add(
                (
                    row["panel"],
                    row["chrom"],
                    int(row["pos"]),
                    row["ref"],
                    row["alt"],
                )
            )

    return candidates


def read_vcf(path):
    samples = []
    variants = []
    genotype_rows = []

    with path.open() as handle:
        for line in handle:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                fields = line.rstrip("\n").split("\t")
                samples = fields[9:]
                continue

            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")

            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]

            if len(ref) != 1:
                continue

            if len(alt) != 1 or "," in alt:
                continue

            format_fields = fields[8].split(":")

            format_index = {
                name: index
                for index, name in enumerate(format_fields)
            }

            if "GT" not in format_index:
                continue

            calls = []

            for sample_field in fields[9:]:
                values = sample_field.split(":")

                gt_index = format_index["GT"]

                gt = (
                    values[gt_index]
                    if gt_index < len(values)
                    else "."
                )

                dosage = normalize_gt(gt)

                dp = None

                if "DP" in format_index:
                    dp_index = format_index["DP"]

                    if dp_index < len(values):
                        dp = as_float(values[dp_index])

                gq = None

                if "GQ" in format_index:
                    gq_index = format_index["GQ"]

                    if gq_index < len(values):
                        gq = as_float(values[gq_index])

                calls.append(
                    (
                        dosage,
                        dp,
                        gq,
                    )
                )

            variants.append(
                (
                    chrom,
                    pos,
                    ref,
                    alt,
                )
            )

            genotype_rows.append(calls)

    if not samples:
        raise RuntimeError(
            f"No VCF header found in {path}"
        )

    return samples, variants, genotype_rows


def choose_orientation(
    scores,
    panel,
    variant,
):
    chrom, pos, ref, alt = variant

    for scheme in (
        "DP2_GQ20",
        "called_GT",
    ):
        row = scores.get(
            (
                panel,
                chrom,
                pos,
                ref,
                alt,
                scheme,
            )
        )

        if row and row.get("selected_orientation"):
            return row["selected_orientation"]

    return None


def phenotype_value(
    row,
    trait,
):
    value = as_float(
        row.get(trait)
    )

    if value in {0.0, 1.0}:
        return int(value)

    return None


def build_panel(
    panel,
    config,
    metadata,
    scores,
    inclusive,
):
    vcf_path = (
        VCF_DIR
        / config["vcf"]
    )

    samples, variants, genotype_rows = read_vcf(
        vcf_path
    )

    sample_index = {
        sample: index
        for index, sample in enumerate(samples)
    }

    cohort = []

    for sample in samples:
        row = metadata.get(sample)

        if row is None:
            continue

        if sample in FOCAL:
            continue

        taxon = (
            row.get("Taxon Name")
            or ""
        ).lower()

        if "piano keys" in taxon:
            continue

        phenotype = phenotype_value(
            row,
            config["trait"],
        )

        if phenotype is None:
            continue

        cohort.append(
            (
                phenotype,
                row.get("Taxon Name") or "",
                row.get("Data Type") or "",
                sample,
            )
        )

    cohort.sort(
        key=lambda item: (
            item[0],
            item[1],
            item[2],
            item[3],
        )
    )

    ordered_samples = [
        item[3]
        for item in cohort
    ]

    phenotype_vector = np.array(
        [
            item[0]
            for item in cohort
        ],
        dtype=int,
    )

    kept_variants = []
    kept_calls = []
    orientations = []

    for variant, calls in zip(
        variants,
        genotype_rows,
    ):
        orientation = choose_orientation(
            scores,
            panel,
            variant,
        )

        if orientation is None:
            continue

        kept_variants.append(variant)
        kept_calls.append(calls)
        orientations.append(orientation)

    n_samples = len(ordered_samples)
    n_sites = len(kept_variants)

    matrix = np.full(
        (
            n_samples,
            n_sites,
        ),
        np.nan,
    )

    for column, (
        calls,
        orientation,
    ) in enumerate(
        zip(
            kept_calls,
            orientations,
        )
    ):
        for row_number, sample in enumerate(
            ordered_samples
        ):
            dosage, dp, gq = calls[
                sample_index[sample]
            ]

            if dosage is None:
                continue

            if dp is None or gq is None:
                continue

            if dp < 2 or gq < 20:
                continue

            associated_dosage = dosage

            if orientation.startswith("REF_"):
                associated_dosage = (
                    2 - dosage
                )

            matrix[
                row_number,
                column,
            ] = associated_dosage

    x = np.arange(n_sites)

    raw_ba = np.full(
        n_sites,
        np.nan,
    )

    qc_ba = np.full(
        n_sites,
        np.nan,
    )

    raw_eligible = np.zeros(
        n_sites,
        dtype=bool,
    )

    qc_eligible = np.zeros(
        n_sites,
        dtype=bool,
    )

    inclusive_mask = np.zeros(
        n_sites,
        dtype=bool,
    )

    for index, variant in enumerate(
        kept_variants
    ):
        chrom, pos, ref, alt = variant

        score_targets = (
            (
                "called_GT",
                raw_ba,
                raw_eligible,
            ),
            (
                "DP2_GQ20",
                qc_ba,
                qc_eligible,
            ),
        )

        for (
            scheme,
            target,
            eligibility,
        ) in score_targets:
            row = scores.get(
                (
                    panel,
                    chrom,
                    pos,
                    ref,
                    alt,
                    scheme,
                )
            )

            if row:
                value = as_float(
                    row.get(
                        "balanced_accuracy"
                    )
                )

                if value is not None:
                    target[index] = value

                eligibility[index] = (
                    str(
                        row.get(
                            "score_eligible",
                            "",
                        )
                    ).lower()
                    == "true"
                )

        inclusive_mask[index] = (
            (
                panel,
                chrom,
                pos,
                ref,
                alt,
            )
            in inclusive
        )

    absent_n = int(
        np.sum(
            phenotype_vector == 0
        )
    )

    present_n = int(
        np.sum(
            phenotype_vector == 1
        )
    )

    boundary = (
        absent_n - 0.5
    )

    figure = plt.figure(
        figsize=(
            11.2,
            8.2,
        )
    )

    grid = figure.add_gridspec(
        2,
        2,
        width_ratios=[
            0.18,
            1.0,
        ],
        height_ratios=[
            4.3,
            1.5,
        ],
        hspace=0.05,
        wspace=0.025,
    )

    phenotype_axis = figure.add_subplot(
        grid[0, 0]
    )

    genotype_axis = figure.add_subplot(
        grid[0, 1]
    )

    score_axis = figure.add_subplot(
        grid[1, 1],
        sharex=genotype_axis,
    )

    genotype_colours = ListedColormap(
        [
            "#B12A90",
            "#FCA636",
            "#2C9CCB",
        ]
    )

    genotype_colours.set_bad(
        "#E6E6E6"
    )

    genotype_norm = BoundaryNorm(
        [
            -0.5,
            0.5,
            1.5,
            2.5,
        ],
        genotype_colours.N,
    )

    genotype_axis.imshow(
        matrix,
        aspect="auto",
        interpolation="nearest",
        origin="upper",
        cmap=genotype_colours,
        norm=genotype_norm,
        rasterized=True,
    )

    genotype_axis.axhline(
        boundary,
        color="black",
        linewidth=0.65,
    )

    genotype_axis.set_xlim(
        -0.5,
        n_sites - 0.5,
    )

    genotype_axis.set_yticks([])

    genotype_axis.tick_params(
        axis="x",
        which="both",
        bottom=False,
        labelbottom=False,
    )

    genotype_axis.set_title(
        (
            f'{config["title"]}: '
            "genotype plot and aligned comparative score"
        ),
        loc="left",
        fontsize=12,
        fontweight="bold",
        pad=8,
    )

    for spine in genotype_axis.spines.values():
        spine.set_linewidth(0.6)

    phenotype_colours = ListedColormap(
        [
            "#E6E6E6",
            "#D55E00",
        ]
    )

    phenotype_axis.imshow(
        phenotype_vector[:, None],
        aspect="auto",
        interpolation="nearest",
        origin="upper",
        cmap=phenotype_colours,
        vmin=0,
        vmax=1,
        rasterized=True,
    )

    phenotype_axis.axhline(
        boundary,
        color="black",
        linewidth=0.65,
    )

    phenotype_axis.set_xticks([])
    phenotype_axis.set_yticks([])

    phenotype_axis.text(
        -0.15,
        (
            1
            - absent_n
            / 2
            / n_samples
        ),
        (
            "Absent\n"
            f"(n={absent_n})"
        ),
        transform=phenotype_axis.transAxes,
        ha="right",
        va="center",
        fontsize=8,
    )

    phenotype_axis.text(
        -0.15,
        (
            1
            - (
                absent_n
                + present_n / 2
            )
            / n_samples
        ),
        (
            "Present\n"
            f"(n={present_n})"
        ),
        transform=phenotype_axis.transAxes,
        ha="right",
        va="center",
        fontsize=8,
        color="#A33F00",
    )

    for spine in phenotype_axis.spines.values():
        spine.set_visible(False)

    score_axis.axhline(
        0.5,
        color="#6F6F6F",
        linewidth=0.7,
        linestyle=(
            0,
            (
                3,
                3,
            ),
        ),
    )

    score_axis.scatter(
        x[~raw_eligible],
        raw_ba[~raw_eligible],
        s=5,
        color="#A6A6A6",
        alpha=0.22,
        linewidths=0,
    )

    score_axis.scatter(
        x[raw_eligible],
        raw_ba[raw_eligible],
        s=8,
        color="#777777",
        alpha=0.55,
        linewidths=0,
    )

    score_axis.scatter(
        x[~qc_eligible],
        qc_ba[~qc_eligible],
        s=5,
        color="#FCA636",
        alpha=0.22,
        linewidths=0,
    )

    score_axis.scatter(
        x[qc_eligible],
        qc_ba[qc_eligible],
        s=9,
        color="#7A0177",
        alpha=0.72,
        linewidths=0,
    )

    if np.any(inclusive_mask):
        score_axis.scatter(
            x[inclusive_mask],
            qc_ba[inclusive_mask],
            marker="v",
            s=36,
            facecolor="#D55E00",
            edgecolor="black",
            linewidth=0.45,
            zorder=5,
        )

    score_axis.set_ylim(
        0.47,
        1.02,
    )

    score_axis.set_ylabel(
        "Balanced accuracy",
        fontsize=9,
    )

    score_axis.set_xlabel(
        "Hmr18 position (Mb)",
        fontsize=9,
    )

    score_axis.tick_params(
        labelsize=8
    )

    score_axis.spines[
        "top"
    ].set_visible(False)

    score_axis.spines[
        "right"
    ].set_visible(False)

    tick_count = min(
        6,
        n_sites,
    )

    tick_indices = np.unique(
        np.linspace(
            0,
            n_sites - 1,
            tick_count,
        )
        .round()
        .astype(int)
    )

    tick_labels = [
        (
            f"{kept_variants[index][1] / 1_000_000:.4f}"
        )
        for index in tick_indices
    ]

    score_axis.set_xticks(
        tick_indices,
        tick_labels,
    )

    score_legend = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#777777",
            markeredgecolor="none",
            markersize=5,
            label="All called GT",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor="#7A0177",
            markeredgecolor="none",
            markersize=5,
            label="DP≥2; GQ≥20",
        ),
        Line2D(
            [0],
            [0],
            marker="v",
            color="none",
            markerfacecolor="#D55E00",
            markeredgecolor="black",
            markeredgewidth=0.45,
            markersize=6,
            label="Comparative-informative (min BA >= 0.80)",
        ),
    ]

    score_axis.legend(
        handles=score_legend,
        loc="lower right",
        frameon=False,
        fontsize=7.5,
        ncol=3,
        handletextpad=0.4,
        columnspacing=0.9,
    )

    genotype_legend = [
        Line2D(
            [0],
            [0],
            marker="s",
            color="none",
            markerfacecolor="#B12A90",
            markeredgecolor="#999999",
            markersize=6,
            label="0 associated-allele copies",
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            color="none",
            markerfacecolor="#FCA636",
            markeredgecolor="none",
            markersize=6,
            label="1 copy (heterozygote)",
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            color="none",
            markerfacecolor="#2C9CCB",
            markeredgecolor="none",
            markersize=6,
            label="2 copies",
        ),
        Line2D(
            [0],
            [0],
            marker="s",
            color="none",
            markerfacecolor="#E6E6E6",
            markeredgecolor="none",
            markersize=6,
            label="Missing/failed QC",
        ),
    ]

    genotype_axis.legend(
        handles=genotype_legend,
        loc="upper center",
        bbox_to_anchor=(
            0.5,
            1.01,
        ),
        frameon=False,
        fontsize=7.3,
        ncol=4,
        handletextpad=0.35,
        columnspacing=0.85,
    )

    figure.text(
        0.995,
        0.005,
        (
            "Rows: comparative butterflies only, sorted by phenotype. "
            "Columns: biallelic SNPs in genomic order. "
            "Genotypes are polarized using the DP≥2/GQ≥20 "
            "selected allele orientation; focal and piano-key "
            "samples are excluded."
        ),
        ha="right",
        va="bottom",
        fontsize=6.8,
        color="#444444",
    )

    output = (
        OUTPUT_DIR
        / (
            f"optix_{panel}_aligned_"
            "genoplot_score_inclusive.pdf"
        )
    )

    figure.savefig(
        output,
        bbox_inches="tight",
    )

    plt.close(figure)

    return {
        "panel": panel,
        "trait": config["trait"],
        "samples": n_samples,
        "phenotype_0": absent_n,
        "phenotype_1": present_n,
        "SNP_columns": n_sites,
        "inclusive_sites": int(
            np.sum(inclusive_mask)
        ),
        "output": str(output),
    }


def main():
    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    CHECKPOINT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    required = [
        METADATA_FILE,
        SCORE_FILE,
        CANDIDATE_FILE,
    ]

    required.extend(
        VCF_DIR / config["vcf"]
        for config in PANELS.values()
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

    metadata = read_metadata()
    scores = read_scores()

    inclusive = (
        read_inclusive_candidates()
    )

    summaries = []

    for panel, config in PANELS.items():
        summaries.append(
            build_panel(
                panel,
                config,
                metadata,
                scores,
                inclusive,
            )
        )

    summary_path = (
        CHECKPOINT_DIR
        / "optix_aligned_genoplot_score_figure_summary.tsv"
    )

    with summary_path.open(
        "w",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                summaries[0]
            ),
            delimiter="\t",
        )

        writer.writeheader()
        writer.writerows(summaries)

    print(
        "ALIGNED OPTIX GENOTYPE-PLOT "
        "AND SCORE PANELS"
    )

    print()

    print(
        "panel\ttrait\tsamples\t"
        "phenotype_0\tphenotype_1\t"
        "SNP_columns\tinclusive_sites"
    )

    for row in summaries:
        print(
            f'{row["panel"]}\t'
            f'{row["trait"]}\t'
            f'{row["samples"]}\t'
            f'{row["phenotype_0"]}\t'
            f'{row["phenotype_1"]}\t'
            f'{row["SNP_columns"]}\t'
            f'{row["inclusive_sites"]}'
        )

    print()
    print(f"Figures: {OUTPUT_DIR}")
    print(f"Summary: {summary_path}")


if __name__ == "__main__":
    main()
