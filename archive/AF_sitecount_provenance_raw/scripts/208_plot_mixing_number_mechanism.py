#!/usr/bin/env python3
import argparse
import importlib.util
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap, LogNorm
import numpy as np
import pandas as pd


ROOT = Path(
    "analysis/genomewide_mixing_corrected_v1"
)
EXPANDED = ROOT / "expanded_calibration_v2"
SIBLING = ROOT / "empirical_full_sib_calibration_v1"
OUT = ROOT / "mixing_number_narrative_v1"
FIGURES = OUT / "figures"

PDF = (
    FIGURES
    / "mixing_number_mechanism_narrative.pdf"
)
PNG = (
    FIGURES
    / "mixing_number_mechanism_summary.png"
)
AUDIT = (
    OUT
    / "mixing_number_mechanism.audit.json"
)

BLUE = "#0072B2"
ORANGE = "#D55E00"
GREEN = "#009E73"
PURPLE = "#7B3294"
DARK = "#222222"
MID = "#666666"
LIGHT = "#E8E8E8"

GREY_CMAP = LinearSegmentedColormap.from_list(
    "mixing_grey",
    [
        "#FFFFFF",
        "#E5E5E5",
        "#A8A8A8",
        "#555555",
        "#181818",
    ],
)


def load_module(path, name):
    specification = (
        importlib.util.spec_from_file_location(
            name,
            path,
        )
    )

    module = (
        importlib.util.module_from_spec(
            specification
        )
    )

    specification.loader.exec_module(
        module
    )

    return module


def first_number(mapping, names):
    for name in names:
        if name in mapping:
            return float(
                mapping[name]
            )

    raise KeyError(
        "None of these keys were found: "
        f"{names}"
    )


def load_observed_fit():
    plotter = load_module(
        (
            "scripts/"
            "188_plot_corrected_genomewide_narrative.py"
        ),
        "mixing_number_observed_fit_source",
    )

    path = (
        Path(plotter.BASE)
        / "observed_expanded.json"
    )

    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(
            "Missing or empty authoritative observed-fit "
            f"file: {path}"
        )

    fit = json.loads(path.read_text())
    full = fit["full"]
    alternative = full["alternative"]
    null = full["null"]

    c_left = first_number(
        alternative,
        ("c_left", "left_mixing"),
    )

    c_right = first_number(
        alternative,
        ("c_right", "right_mixing"),
    )

    if "LR" in full:
        lr = float(full["LR"])
    else:
        alt_logl = first_number(
            alternative,
            ("logL", "log_likelihood"),
        )

        null_logl = first_number(
            null,
            ("logL", "log_likelihood"),
        )

        lr = 2.0 * (alt_logl - null_logl)

    return (
        path,
        fit,
        c_left,
        c_right,
        lr,
    )

def load_emission():
    path = (
        ROOT
        / "control_qc"
        / "emission_calibration.diagnostic.json"
    )

    records = json.loads(
        path.read_text()
    )

    selected = None

    for record in records:
        if record.get("scope") == "snps_only":
            selected = record
            break

    if selected is None:
        selected = records[0]

    return (
        path,
        selected,
    )


def require_columns(
    table,
    required,
    label,
):
    missing = (
        set(required)
        - set(table.columns)
    )

    if missing:
        raise RuntimeError(
            f"{label} missing columns: "
            f"{sorted(missing)}"
        )


def load_inputs(load_counts=True):
    (
        observed_path,
        observed_fit,
        c_left,
        c_right,
        likelihood_ratio,
    ) = load_observed_fit()

    (
        emission_path,
        emission,
    ) = load_emission()

    simulations_path = (
        EXPANDED
        / "simulations.tsv.gz"
    )

    sibling_path = (
        SIBLING
        / "final"
        / "all_empirical_sibling_simulations.tsv.gz"
    )

    simulations = pd.read_csv(
        simulations_path,
        sep="\t",
    )

    siblings = pd.read_csv(
        sibling_path,
        sep="\t",
    )

    simulation_columns = (
        "family",
        "generating_model",
        "true_c_left",
        "true_c_right",
        "genotype_contrast",
        "estimated_c_left_alt",
        "estimated_c_right_alt",
    )

    sibling_columns = (
        "scenario",
        "generating_model",
        "true_c_left",
        "true_c_right",
        "genotype_contrast",
        "IBD0_fraction",
        "IBD1_fraction",
        "IBD2_fraction",
        "estimated_c_left_alt",
        "estimated_c_right_alt",
    )

    require_columns(
        simulations,
        simulation_columns,
        "expanded simulations",
    )

    require_columns(
        siblings,
        sibling_columns,
        "sibling simulations",
    )

    observed = None

    if load_counts:
        plotter = load_module(
            (
                "scripts/"
                "188_plot_corrected_genomewide_narrative.py"
            ),
            "mixing_narrative_stage188",
        )

        observed = (
            plotter.stage.load_data(
                calibration=False
            )
        )

        required_arrays = (
            "left_alt",
            "left_total",
            "right_alt",
            "right_total",
        )

        for key in required_arrays:
            if key not in observed:
                raise RuntimeError(
                    f"Observed data missing {key}"
                )

        sizes = {
            len(observed[key])
            for key
            in required_arrays
        }

        if (
            len(sizes) != 1
            or min(sizes) != 1000000
        ):
            raise RuntimeError(
                "Expected one million equal-scope "
                "sites; observed lengths="
                f"{sizes}"
            )

    return {
        "observed_path": observed_path,
        "observed_fit": observed_fit,
        "c_left": c_left,
        "c_right": c_right,
        "lr": likelihood_ratio,
        "emission_path": emission_path,
        "emission": emission,
        "simulations_path": simulations_path,
        "simulations": simulations,
        "sibling_path": sibling_path,
        "siblings": siblings,
        "observed": observed,
    }


def style_axis(axis, grid=False):
    axis.spines[
        [
            "top",
            "right",
        ]
    ].set_visible(False)

    axis.tick_params(
        labelsize=8,
        colors="#444444",
    )

    if grid:
        axis.grid(
            axis="y",
            color=LIGHT,
            linewidth=0.7,
            zorder=0,
        )


def page_header(
    figure,
    number,
    title,
    subtitle,
):
    figure.text(
        0.055,
        0.965,
        f"STEP {number}",
        ha="left",
        va="top",
        fontsize=9,
        fontweight="bold",
        color=MID,
    )

    figure.text(
        0.055,
        0.925,
        title,
        ha="left",
        va="top",
        fontsize=19,
        fontweight="bold",
        color=DARK,
    )

    figure.text(
        0.055,
        0.875,
        subtitle,
        ha="left",
        va="top",
        fontsize=9.5,
        color=MID,
    )


def flow_box(
    axis,
    x_value,
    y_value,
    width,
    height,
    title,
    body,
    color,
):
    from matplotlib.patches import FancyBboxPatch

    patch = FancyBboxPatch(
        (
            x_value,
            y_value,
        ),
        width,
        height,
        boxstyle=(
            "round,pad=0.018,"
            "rounding_size=0.02"
        ),
        facecolor="white",
        edgecolor=color,
        linewidth=1.5,
        transform=axis.transAxes,
    )

    axis.add_patch(patch)

    axis.text(
        x_value + 0.025,
        y_value + height - 0.045,
        title,
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        fontweight="bold",
        color=color,
    )

    axis.text(
        x_value + 0.025,
        y_value + height - 0.105,
        body,
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=8.2,
        color=DARK,
        linespacing=1.35,
    )


def observed_fractions(observed):
    left_alt = np.asarray(
        observed["left_alt"],
        dtype=float,
    )

    left_total = np.asarray(
        observed["left_total"],
        dtype=float,
    )

    right_alt = np.asarray(
        observed["right_alt"],
        dtype=float,
    )

    right_total = np.asarray(
        observed["right_total"],
        dtype=float,
    )

    keep = (
        np.isfinite(
            left_alt
            + left_total
            + right_alt
            + right_total
        )
        & (left_total > 0)
        & (right_total > 0)
    )

    return (
        left_alt[keep]
        / left_total[keep],
        right_alt[keep]
        / right_total[keep],
        left_total[keep],
        right_total[keep],
    )


def load_unmixed_control_af():
    import subprocess

    sample = "ERR1307111"

    cache = (
        ROOT
        / "mixing_estimation_narrative_v2"
        / "tables"
        / f"{sample}.chr1.empirical_unmixed_AF.npz"
    )

    if cache.is_file() and cache.stat().st_size > 0:
        with np.load(
            cache,
            allow_pickle=False,
        ) as stored:
            allele_fraction = np.asarray(
                stored["allele_fraction"],
                dtype=float,
            )

            source = str(
                stored["source"].item()
            )

        cache_valid = (
            len(allele_fraction) >= 10_000
            and np.isfinite(allele_fraction).all()
            and (allele_fraction >= 0).all()
            and (allele_fraction <= 1).all()
        )

        if cache_valid:
            print(
                "UNMIXED_CONTROL_AF"
                f"\tsample={sample}"
                f"\tsites={len(allele_fraction)}"
                "\tcache=REUSED",
                flush=True,
            )

            return (
                allele_fraction,
                sample,
                source,
            )

    call_root = Path(
        "analysis/genomewide_mixing_input/"
        "variant_calls/joint_24/by_chromosome"
    )

    candidates = sorted(
        call_root.glob(
            "chr1_*/"
            "all_24_samples.joint_at_control_sites.bcf"
        )
    )

    if len(candidates) != 1:
        raise RuntimeError(
            "Expected exactly one chromosome-1 joint BCF; "
            f"found {len(candidates)}: {candidates}"
        )

    bcf = candidates[0]

    bcftools = Path(
        "/gpfs/automountdir/gpfs/homes/CCAS/home/"
        "jjhanly/.local/share/mamba/envs/"
        "colias_w_popgen/bin/bcftools"
    )

    if not bcftools.is_file():
        raise RuntimeError(
            f"bcftools not found: {bcftools}"
        )

    process = subprocess.Popen(
        [
            str(bcftools),
            "query",
            "-s",
            sample,
            "-f",
            "[%AD]\\n",
            str(bcf),
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    values = []

    for line in process.stdout:
        fields = (
            line.rstrip("\n")
            .split(",")
        )

        if (
            len(fields) < 2
            or fields[0] == "."
            or fields[1] == "."
        ):
            continue

        try:
            ref_count = int(
                fields[0]
            )

            alt_count = int(
                fields[1]
            )
        except ValueError:
            continue

        total = (
            ref_count
            + alt_count
        )

        if total > 0:
            values.append(
                alt_count / total
            )

    stderr = process.stderr.read()
    status = process.wait()

    if status != 0:
        raise RuntimeError(
            "bcftools query failed with "
            f"status {status}:\n"
            f"{stderr[-2000:]}"
        )

    allele_fraction = np.asarray(
        values,
        dtype=float,
    )

    if len(allele_fraction) < 10_000:
        raise RuntimeError(
            "Only "
            f"{len(allele_fraction)} "
            "usable control-AF observations"
        )

    cache.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    temporary = cache.with_suffix(
        ".npz.tmp"
    )

    with temporary.open("wb") as handle:
        np.savez_compressed(
            handle,
            allele_fraction=allele_fraction,
            sample=np.asarray(sample),
            source=np.asarray(str(bcf)),
        )

    temporary.replace(
        cache
    )

    print(
        "UNMIXED_CONTROL_AF"
        f"\tsample={sample}"
        f"\tsites={len(allele_fraction)}"
        f"\tsource={bcf}"
        "\tcache=CREATED",
        flush=True,
    )

    return (
        allele_fraction,
        sample,
        str(bcf),
    )


def page_inputs(pdf, inputs):
    observed = inputs["observed"]

    (
        left_af,
        right_af,
        left_dp,
        right_dp,
    ) = observed_fractions(
        observed
    )

    (
        control_af,
        control_sample,
        control_source,
    ) = load_unmixed_control_af()

    figure = plt.figure(
        figsize=(12, 8.2)
    )

    page_header(
        figure,
        1,
        (
            "From joint VCF records to "
            "mixing-model observations"
        ),
        (
            "Pooled Left/Right AD fractions are contrasted "
            "with one empirical unmixed control individual; "
            "focal GT calls are not used."
        ),
    )

    flow = figure.add_axes([
        0.055,
        0.52,
        0.89,
        0.28,
    ])

    flow.axis("off")

    boxes = [
        (
            0.00,
            "Joint chromosome BCFs",
            (
                "24 samples\n"
                "18 controls + 6 focal\n"
                "GT:AD:DP:PL fields"
            ),
            BLUE,
        ),
        (
            0.255,
            "QC and site extraction",
            (
                "biallelic autosomal sites\n"
                "control DP/PL support\n"
                "Left and Right AD retained"
            ),
            GREEN,
        ),
        (
            0.510,
            "Equal-scope dataset",
            (
                "1,000,000 autosomal sites\n"
                "Left ALT, Left AD sum\n"
                "Right ALT, Right AD sum"
            ),
            PURPLE,
        ),
        (
            0.765,
            "Likelihood input",
            (
                "kL | nL and kR | nR\n"
                "control genotype states\n"
                "emission calibration"
            ),
            ORANGE,
        ),
    ]

    for (
        xpos,
        title,
        body,
        color,
    ) in boxes:
        flow_box(
            flow,
            xpos,
            0.17,
            0.205,
            0.66,
            title,
            body,
            color,
        )

    for xpos in (
        0.225,
        0.480,
        0.735,
    ):
        flow.annotate(
            "",
            xy=(
                xpos + 0.025,
                0.50,
            ),
            xytext=(
                xpos,
                0.50,
            ),
            xycoords="axes fraction",
            arrowprops=dict(
                arrowstyle="->",
                color="#777777",
                linewidth=1.4,
            ),
        )

    axes = [
        figure.add_axes([
            0.075,
            0.10,
            0.25,
            0.30,
        ]),
        figure.add_axes([
            0.385,
            0.10,
            0.25,
            0.30,
        ]),
        figure.add_axes([
            0.695,
            0.10,
            0.25,
            0.30,
        ]),
    ]

    positive_depth = np.concatenate([
        left_dp[
            left_dp > 0
        ],
        right_dp[
            right_dp > 0
        ],
    ])

    maximum = max(
        2.0,
        float(
            positive_depth.max()
        ),
    )

    depth_bins = np.geomspace(
        1.0,
        maximum + 1.0,
        55,
    )

    axes[0].hist(
        left_dp,
        bins=depth_bins,
        histtype="step",
        linewidth=1.8,
        color=BLUE,
        alpha=0.72,
        label="Left",
    )

    axes[0].hist(
        right_dp,
        bins=depth_bins,
        histtype="step",
        linewidth=1.8,
        color=ORANGE,
        alpha=0.72,
        label="Right",
    )

    axes[0].set_xscale(
        "log"
    )

    axes[0].set_yscale(
        "log"
    )

    axes[0].set(
        xlabel="AD_REF + AD_ALT",
        ylabel="Sites",
        title="Observed effective depth",
    )

    axes[0].legend(
        frameon=False,
        fontsize=8,
    )

    af_bins = np.linspace(
        0,
        1,
        101,
    )

    series = (
        (
            left_af,
            BLUE,
            0.72,
            "Pooled Left",
        ),
        (
            right_af,
            ORANGE,
            0.72,
            "Pooled Right",
        ),
        (
            control_af,
            GREEN,
            0.88,
            f"Unmixed {control_sample}",
        ),
    )

    for (
        values,
        color,
        alpha,
        label,
    ) in series:
        weights = np.full(
            len(values),
            1.0 / len(values),
        )

        axes[1].hist(
            values,
            bins=af_bins,
            weights=weights,
            histtype="step",
            linewidth=1.7,
            color=color,
            alpha=alpha,
            label=label,
        )

    axes[1].set_yscale(
        "symlog",
        linthresh=1e-5,
        linscale=0.7,
    )

    axes[1].set(
        xlabel="ALT-read fraction",
        ylabel=(
            "Fraction of sites "
            "(compressed scale)"
        ),
        title=(
            "Pooled sides versus "
            "one unmixed individual"
        ),
    )

    axes[1].text(
        0.98,
        0.73,
        "symmetric-log y scale",
        transform=axes[1].transAxes,
        ha="right",
        va="top",
        fontsize=7.1,
        color=MID,
    )

    axes[1].legend(
        frameon=False,
        fontsize=6.8,
        loc="upper center",
    )

    emission = inputs[
        "emission"
    ]

    axes[2].axis(
        "off"
    )

    axes[2].text(
        0,
        1,
        "Control-calibrated emission",
        ha="left",
        va="top",
        fontsize=11,
        fontweight="bold",
    )

    axes[2].text(
        0,
        0.83,
        (
            "Sequencing error  ε = "
            f"{float(emission['error_rate']):.5f}\n"
            "Allelic-bias log odds  b = "
            f"{float(emission['bias_log_odds']):.3f}\n"
            "Beta-binomial concentration  φ = "
            f"{float(emission['concentration']):.2f}\n\n"
            "Green: empirical, unpooled control individual\n"
            "from chromosome 1 of the same joint-call set.\n\n"
            "Emission parameters are learned from all 18 controls,\n"
            "then held fixed while cL and cR are estimated."
        ),
        ha="left",
        va="top",
        fontsize=8.8,
        color=DARK,
        linespacing=1.42,
    )

    for axis in axes[:2]:
        style_axis(
            axis
        )

    figure.text(
        0.945,
        0.025,
        (
            f"Pooled observations: {len(left_af):,} sites; "
            f"unmixed control: {len(control_af):,} "
            "callable chr1 sites"
        ),
        ha="right",
        fontsize=8.2,
        color=MID,
    )

    pdf.savefig(
        figure
    )

    plt.close(
        figure
    )


def genotype_centroids(
    c_left,
    c_right,
):
    rows = []

    for left_genotype in (
        0.0,
        0.5,
        1.0,
    ):
        for right_genotype in (
            0.0,
            0.5,
            1.0,
        ):
            rows.append(
                (
                    left_genotype,
                    right_genotype,
                    (
                        (1.0 - c_left)
                        * left_genotype
                        + c_left
                        * right_genotype
                    ),
                    (
                        (1.0 - c_right)
                        * right_genotype
                        + c_right
                        * left_genotype
                    ),
                )
            )

    return np.asarray(
        rows,
        dtype=float,
    )


def observed_hex(
    axis,
    left_af,
    right_af,
    centroids,
    title,
):
    artist = axis.hexbin(
        left_af,
        right_af,
        gridsize=70,
        extent=(
            0.0,
            1.0,
            0.0,
            1.0,
        ),
        mincnt=1,
        cmap=GREY_CMAP,
        norm=LogNorm(),
        linewidths=0,
    )

    axis.plot(
        [
            0.0,
            1.0,
        ],
        [
            0.0,
            1.0,
        ],
        color="#CCCCCC",
        linewidth=0.8,
        zorder=0,
    )

    axis.scatter(
        centroids[:, 2],
        centroids[:, 3],
        s=36,
        facecolors="white",
        edgecolors=BLUE,
        linewidths=1.1,
        zorder=4,
    )

    opposing = (
        (
            (centroids[:, 0] == 0)
            & (centroids[:, 1] == 1)
        )
        | (
            (centroids[:, 0] == 1)
            & (centroids[:, 1] == 0)
        )
    )

    axis.scatter(
        centroids[
            opposing,
            2,
        ],
        centroids[
            opposing,
            3,
        ],
        s=40,
        marker="D",
        color=ORANGE,
        edgecolors="white",
        linewidths=0.5,
        zorder=5,
    )

    axis.set(
        xlim=(-0.02, 1.02),
        ylim=(-0.02, 1.02),
        xlabel=(
            "Left ALT-read fraction"
        ),
        ylabel=(
            "Right ALT-read fraction"
        ),
        title=title,
    )

    axis.set_aspect("equal")
    style_axis(axis)

    return artist


def page_mechanism(pdf, inputs):
    (
        left_af,
        right_af,
        _,
        _,
    ) = observed_fractions(
        inputs["observed"]
    )

    c_left = inputs["c_left"]
    c_right = inputs["c_right"]

    unmixed = genotype_centroids(
        0.0,
        0.0,
    )

    fitted = genotype_centroids(
        c_left,
        c_right,
    )

    figure = plt.figure(
        figsize=(12, 8.2)
    )

    page_header(
        figure,
        2,
        "The reciprocal-mixture mechanism",
        (
            "Each side retains most of its own "
            "molecules and receives a fitted "
            "fraction from the other side."
        ),
    )

    equation_axis = figure.add_axes(
        [
            0.06,
            0.66,
            0.88,
            0.14,
        ]
    )

    equation_axis.axis("off")

    equation_axis.text(
        0.5,
        0.75,
        (
            r"$q_L=(1-c_L)p_L+c_Lp_R$"
            "        "
            r"$q_R=(1-c_R)p_R+c_Rp_L$"
        ),
        ha="center",
        va="center",
        fontsize=18,
        color=DARK,
    )

    equation_axis.text(
        0.5,
        0.24,
        (
            r"Observed counts:  "
            r"$k_L\sim\mathrm{BetaBinomial}"
            r"(n_L,q_L;\epsilon,b,\phi)$"
            "        "
            r"$k_R\sim\mathrm{BetaBinomial}"
            r"(n_R,q_R;\epsilon,b,\phi)$"
        ),
        ha="center",
        va="center",
        fontsize=11,
        color=MID,
    )

    axes = [
        figure.add_axes(
            [
                0.07,
                0.10,
                0.38,
                0.49,
            ]
        ),
        figure.add_axes(
            [
                0.55,
                0.10,
                0.38,
                0.49,
            ]
        ),
    ]

    axes[0].scatter(
        unmixed[:, 0],
        unmixed[:, 1],
        s=38,
        facecolors="white",
        edgecolors=BLUE,
        zorder=3,
    )

    for source, target in zip(
        unmixed,
        fitted,
    ):
        axes[0].annotate(
            "",
            xy=(
                target[2],
                target[3],
            ),
            xytext=(
                source[0],
                source[1],
            ),
            arrowprops={
                "arrowstyle": "->",
                "color": "#999999",
                "linewidth": 1.0,
            },
        )

    axes[0].scatter(
        fitted[:, 2],
        fitted[:, 3],
        s=40,
        color=ORANGE,
        edgecolors="white",
        linewidths=0.5,
        zorder=4,
    )

    axes[0].plot(
        [
            0.0,
            1.0,
        ],
        [
            0.0,
            1.0,
        ],
        color="#DDDDDD",
        linewidth=0.8,
    )

    axes[0].set(
        xlim=(-0.04, 1.04),
        ylim=(-0.04, 1.04),
        aspect="equal",
        xlabel=(
            "Expected Left ALT fraction"
        ),
        ylabel=(
            "Expected Right ALT fraction"
        ),
        title=(
            "Centroid displacement at "
            f"cL={c_left:.3f}, "
            f"cR={c_right:.3f}"
        ),
    )

    style_axis(axes[0])

    artist = observed_hex(
        axes[1],
        left_af,
        right_af,
        fitted,
        (
            "One million observed sites "
            "with fitted centroids"
        ),
    )

    color_axis = figure.add_axes(
        [
            0.945,
            0.17,
            0.014,
            0.34,
        ]
    )

    figure.colorbar(
        artist,
        cax=color_axis,
    ).set_label(
        "Sites per hexagon"
    )

    pdf.savefig(figure)
    plt.close(figure)


def nearest_distance(
    left_af,
    right_af,
    centroids,
):
    points = np.column_stack(
        [
            left_af,
            right_af,
        ]
    )

    targets = centroids[:, 2:4]

    minimum = np.full(
        points.shape[0],
        np.inf,
    )

    for target in targets:
        distance = np.sqrt(
            (
                (
                    points
                    - target
                )
                ** 2
            ).sum(axis=1)
        )

        minimum = np.minimum(
            minimum,
            distance,
        )

    return minimum


def page_value(pdf, inputs):
    (
        left_af,
        right_af,
        _,
        _,
    ) = observed_fractions(
        inputs["observed"]
    )

    rng = np.random.default_rng(
        208001
    )

    if len(left_af) > 200000:
        indices = rng.choice(
            len(left_af),
            size=200000,
            replace=False,
        )

        x_values = left_af[indices]
        y_values = right_af[indices]
    else:
        x_values = left_af
        y_values = right_af

    zero_distance = nearest_distance(
        x_values,
        y_values,
        genotype_centroids(
            0.0,
            0.0,
        ),
    )

    fitted_distance = nearest_distance(
        x_values,
        y_values,
        genotype_centroids(
            inputs["c_left"],
            inputs["c_right"],
        ),
    )

    figure = plt.figure(
        figsize=(12, 8.2)
    )

    page_header(
        figure,
        3,
        (
            "From all sites to the "
            "two fitted values"
        ),
        (
            "The optimizer evaluates candidate "
            "cL and cR values across sites while "
            "summing over allowed latent "
            "genotype states."
        ),
    )

    left_axis = figure.add_axes(
        [
            0.07,
            0.13,
            0.39,
            0.58,
        ]
    )

    bins = np.linspace(
        0.0,
        0.45,
        80,
    )

    left_axis.hist(
        zero_distance,
        bins=bins,
        density=True,
        histtype="step",
        linewidth=1.8,
        color=MID,
        label="No-mixing centroids",
    )

    left_axis.hist(
        fitted_distance,
        bins=bins,
        density=True,
        histtype="step",
        linewidth=2.0,
        color=ORANGE,
        label="Fitted-mixing centroids",
    )

    left_axis.set(
        xlabel=(
            "Distance to nearest expected "
            "AF centroid"
        ),
        ylabel="Site density",
        title=(
            "Geometric diagnostic of "
            "centroid alignment"
        ),
    )

    left_axis.legend(
        frameon=False,
        fontsize=8,
    )

    left_axis.text(
        0.02,
        0.02,
        (
            "Diagnostic only; the fitted "
            "objective is the\nbeta-binomial "
            "likelihood, not this distance."
        ),
        transform=left_axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=7.5,
        color=MID,
    )

    style_axis(left_axis)

    right_axis = figure.add_axes(
        [
            0.53,
            0.13,
            0.40,
            0.58,
        ]
    )

    right_axis.axis("off")

    right_axis.text(
        0.0,
        0.98,
        "Likelihood calculation",
        ha="left",
        va="top",
        fontsize=12,
        fontweight="bold",
    )

    right_axis.text(
        0.0,
        0.86,
        (
            "1. Propose cL and cR.\n"
            "2. Move each latent genotype-pair "
            "centroid using the reciprocal "
            "equations.\n"
            "3. Evaluate Left and Right AD counts "
            "with the calibrated emission.\n"
            "4. Weight/sum over the permitted "
            "genotype states at every site.\n"
            "5. Sum log likelihood across "
            "1,000,000 sites.\n"
            "6. Optimize cL and cR jointly."
        ),
        ha="left",
        va="top",
        fontsize=9.5,
        color=DARK,
        linespacing=1.55,
    )

    right_axis.text(
        0.02,
        0.39,
        "Maximum-likelihood mixing",
        fontsize=11,
        fontweight="bold",
        color=MID,
    )

    right_axis.text(
        0.02,
        0.25,
        (
            "Left   "
            f"{100 * inputs['c_left']:.2f}%"
        ),
        fontsize=24,
        fontweight="bold",
        color=BLUE,
    )

    right_axis.text(
        0.02,
        0.11,
        (
            "Right  "
            f"{100 * inputs['c_right']:.2f}%"
        ),
        fontsize=24,
        fontweight="bold",
        color=ORANGE,
    )

    right_axis.text(
        0.98,
        0.02,
        (
            "Founder-state LR in the same "
            "expanded fit: "
            f"{inputs['lr']:.1f}"
        ),
        ha="right",
        va="bottom",
        fontsize=7.5,
        color=MID,
    )

    pdf.savefig(figure)
    plt.close(figure)


def founder_color(values):
    return np.where(
        values.astype(str).str.contains(
            "three"
        ),
        BLUE,
        ORANGE,
    )


def recovery_metrics(table, side):
    truth = table[
        f"true_c_{side}"
    ].to_numpy(
        dtype=float
    )

    estimate = table[
        f"estimated_c_{side}_alt"
    ].to_numpy(
        dtype=float
    )

    keep = np.isfinite(
        truth + estimate
    )

    error = (
        estimate[keep]
        - truth[keep]
    )

    return {
        "n": int(
            keep.sum()
        ),
        "bias": float(
            error.mean()
        ),
        "mae": float(
            np.abs(error).mean()
        ),
        "rmse": float(
            np.sqrt(
                np.mean(
                    error ** 2
                )
            )
        ),
    }


def page_simulations(pdf, inputs):
    table = inputs[
        "simulations"
    ].copy()

    figure, axes = plt.subplots(
        2,
        2,
        figsize=(12, 8.2),
    )

    figure.subplots_adjust(
        left=0.08,
        right=0.96,
        bottom=0.09,
        top=0.77,
        hspace=0.32,
        wspace=0.24,
    )

    page_header(
        figure,
        4,
        (
            "Verification with known "
            "simulated mixing"
        ),
        (
            "Recovery is checked across the "
            "global and focused parameter grids; "
            "the exact-20% repeats test the "
            "observed region densely."
        ),
    )

    colors = founder_color(
        table["generating_model"]
    )

    for column, side in enumerate(
        (
            "left",
            "right",
        )
    ):
        truth = table[
            f"true_c_{side}"
        ].to_numpy(
            dtype=float
        )

        estimate = table[
            f"estimated_c_{side}_alt"
        ].to_numpy(
            dtype=float
        )

        axes[
            0,
            column,
        ].scatter(
            truth,
            estimate,
            c=colors,
            s=7,
            alpha=0.16,
            linewidths=0,
        )

        axes[
            0,
            column,
        ].plot(
            [
                0.0,
                0.45,
            ],
            [
                0.0,
                0.45,
            ],
            color="black",
            linewidth=1.0,
        )

        axes[
            0,
            column,
        ].axhline(
            inputs[
                f"c_{side}"
            ],
            color=GREEN,
            linestyle="--",
            linewidth=1.3,
            label="Observed estimate",
        )

        metrics = recovery_metrics(
            table,
            side,
        )

        axes[
            0,
            column,
        ].set(
            xlim=(-0.01, 0.46),
            ylim=(-0.01, 0.46),
            xlabel=(
                f"True {side} mixing"
            ),
            ylabel=(
                f"Estimated {side} mixing"
            ),
            title=(
                f"{side.capitalize()} recovery "
                "across parameter space"
            ),
        )

        axes[
            0,
            column,
        ].text(
            0.03,
            0.96,
            (
                f"n={metrics['n']:,}\n"
                f"bias={metrics['bias']:+.4f}\n"
                f"RMSE={metrics['rmse']:.4f}"
            ),
            transform=axes[
                0,
                column,
            ].transAxes,
            ha="left",
            va="top",
            fontsize=8,
        )

        style_axis(
            axes[
                0,
                column,
            ]
        )

    exact = table[
        np.isclose(
            table["true_c_left"],
            0.2,
        )
        & np.isclose(
            table["true_c_right"],
            0.2,
        )
    ].copy()

    if (
        exact["family"]
        == "repeat_20"
    ).any():
        exact = exact[
            exact["family"]
            == "repeat_20"
        ].copy()

    bins = np.linspace(
        0.12,
        0.28,
        65,
    )

    for column, side in enumerate(
        (
            "left",
            "right",
        )
    ):
        for (
            model,
            color,
            label,
        ) in (
            (
                "three",
                BLUE,
                "Three-founder",
            ),
            (
                "four",
                ORANGE,
                "Four-founder",
            ),
        ):
            subset = exact[
                exact[
                    "generating_model"
                ]
                .astype(str)
                .str.contains(model)
            ]

            axes[
                1,
                column,
            ].hist(
                subset[
                    f"estimated_c_{side}_alt"
                ],
                bins=bins,
                density=True,
                histtype="step",
                linewidth=1.8,
                color=color,
                label=label,
            )

        axes[
            1,
            column,
        ].axvline(
            0.2,
            color="black",
            linewidth=1.1,
            label="True 20%",
        )

        axes[
            1,
            column,
        ].axvline(
            inputs[
                f"c_{side}"
            ],
            color=GREEN,
            linestyle="--",
            linewidth=1.5,
            label="Observed estimate",
        )

        axes[
            1,
            column,
        ].set(
            xlabel=(
                f"Recovered {side} mixing"
            ),
            ylabel=(
                "Simulation density"
            ),
            title=(
                "Exact-20% recovery: "
                f"{side}"
            ),
        )

        style_axis(
            axes[
                1,
                column,
            ]
        )

    axes[
        1,
        1,
    ].legend(
        frameon=False,
        fontsize=7.5,
    )

    pdf.savefig(figure)
    plt.close(figure)


def page_siblings(pdf, inputs):
    table = inputs[
        "siblings"
    ].copy()

    exact = table[
        table["scenario"]
        == "exact_20"
    ].copy()

    matched = table[
        table["scenario"]
        == "observed_matched"
    ].copy()

    figure, axes = plt.subplots(
        1,
        3,
        figsize=(12.6, 5.8),
    )

    figure.subplots_adjust(
        left=0.07,
        right=0.97,
        bottom=0.18,
        top=0.70,
        wspace=0.28,
    )

    page_header(
        figure,
        5,
        (
            "Robustness to empirical full-sibling "
            "genotype structure"
        ),
        (
            "The same mixing estimator is challenged "
            "with real sibling IBD and heterozygosity "
            "patterns rather than unrelated "
            "idealized donors."
        ),
    )

    bins = np.linspace(
        0.12,
        0.28,
        60,
    )

    for side, color in (
        (
            "left",
            BLUE,
        ),
        (
            "right",
            ORANGE,
        ),
    ):
        values = exact[
            f"estimated_c_{side}_alt"
        ].to_numpy(
            dtype=float
        )

        axes[0].hist(
            values,
            bins=bins,
            density=True,
            histtype="step",
            linewidth=1.9,
            color=color,
            label=side.capitalize(),
        )

        axes[0].axvline(
            inputs[
                f"c_{side}"
            ],
            color=color,
            linestyle="--",
            linewidth=1.0,
        )

    axes[0].axvline(
        0.2,
        color="black",
        linewidth=1.1,
        label="True 20%",
    )

    axes[0].set(
        xlabel="Recovered mixing",
        ylabel="Density",
        title=(
            "Exact-20% sibling recovery"
        ),
    )

    axes[0].legend(
        frameon=False,
        fontsize=7.5,
    )

    axes[1].scatter(
        matched[
            "estimated_c_left_alt"
        ],
        matched[
            "estimated_c_right_alt"
        ],
        c=founder_color(
            matched["generating_model"]
        ),
        s=10,
        alpha=0.28,
        linewidths=0,
    )

    axes[1].scatter(
        [
            inputs["c_left"],
        ],
        [
            inputs["c_right"],
        ],
        marker="*",
        s=150,
        color=GREEN,
        edgecolors="white",
        linewidths=0.7,
        label="Observed fit",
        zorder=5,
    )

    axes[1].plot(
        [
            0.12,
            0.28,
        ],
        [
            0.12,
            0.28,
        ],
        color="#BBBBBB",
        linewidth=0.8,
    )

    axes[1].set(
        xlim=(0.12, 0.28),
        ylim=(0.12, 0.28),
        xlabel=(
            "Recovered Left mixing"
        ),
        ylabel=(
            "Recovered Right mixing"
        ),
        title=(
            "Observed-matched sibling "
            "simulations"
        ),
    )

    axes[1].legend(
        frameon=False,
        fontsize=7.5,
    )

    sibling_only = exact[
        exact["generating_model"]
        == "empirical_full_sib"
    ].copy()

    error = 0.5 * (
        (
            sibling_only[
                "estimated_c_left_alt"
            ]
            - sibling_only[
                "true_c_left"
            ]
        )
        + (
            sibling_only[
                "estimated_c_right_alt"
            ]
            - sibling_only[
                "true_c_right"
            ]
        )
    )

    axes[2].scatter(
        sibling_only[
            "IBD0_fraction"
        ],
        error,
        s=11,
        alpha=0.35,
        color=PURPLE,
        linewidths=0,
    )

    axes[2].axhline(
        0.0,
        color="black",
        linewidth=0.9,
    )

    axes[2].set(
        xlabel=(
            "Sibling IBD0 fraction"
        ),
        ylabel=(
            "Mean recovered-minus-true mixing"
        ),
        title=(
            "Recovery error versus "
            "relatedness realization"
        ),
    )

    for axis in axes:
        style_axis(axis)

    figure.text(
        0.07,
        0.07,
        (
            "Interpretation: the green star is "
            "the animal's fitted value; the "
            "distributions test whether the same "
            "estimator recovers known mixing under "
            "empirical sibling structure."
        ),
        ha="left",
        fontsize=8.5,
        color=MID,
    )

    pdf.savefig(figure)
    plt.close(figure)


def summary_figure(inputs, path):
    (
        left_af,
        right_af,
        _,
        _,
    ) = observed_fractions(
        inputs["observed"]
    )

    fitted = genotype_centroids(
        inputs["c_left"],
        inputs["c_right"],
    )

    table = inputs["simulations"]

    figure = plt.figure(
        figsize=(14, 4.7)
    )

    figure.text(
        0.035,
        0.965,
        (
            "How reciprocal tissue mixing "
            "was estimated"
        ),
        ha="left",
        va="top",
        fontsize=19,
        fontweight="bold",
        color=DARK,
    )

    flow = figure.add_axes(
        [
            0.035,
            0.13,
            0.25,
            0.70,
        ]
    )

    flow.axis("off")

    flow_box(
        flow,
        0.02,
        0.70,
        0.92,
        0.22,
        "1  Joint VCFs",
        (
            "AD/DP/PL from 18 controls and "
            "pooled Left/Right samples"
        ),
        BLUE,
    )

    flow_box(
        flow,
        0.02,
        0.39,
        0.92,
        0.22,
        "2  Reciprocal model",
        (
            "qL=(1-cL)pL+cLpR\n"
            "qR=(1-cR)pR+cRpL"
        ),
        PURPLE,
    )

    flow_box(
        flow,
        0.02,
        0.08,
        0.92,
        0.22,
        "3  Joint likelihood",
        (
            "Beta-binomial emissions; latent "
            "genotype states summed across "
            "1M sites"
        ),
        ORANGE,
    )

    hex_axis = figure.add_axes(
        [
            0.33,
            0.16,
            0.28,
            0.62,
        ]
    )

    observed_hex(
        hex_axis,
        left_af,
        right_af,
        fitted,
        (
            "Observed data and "
            "fitted centroids"
        ),
    )

    recovery_axis = figure.add_axes(
        [
            0.67,
            0.17,
            0.20,
            0.60,
        ]
    )

    recovery_axis.scatter(
        table["true_c_left"],
        table["estimated_c_left_alt"],
        s=5,
        alpha=0.12,
        color=BLUE,
        linewidths=0,
    )

    recovery_axis.plot(
        [
            0.0,
            0.45,
        ],
        [
            0.0,
            0.45,
        ],
        color="black",
        linewidth=0.9,
    )

    recovery_axis.axhline(
        inputs["c_left"],
        color=GREEN,
        linestyle="--",
        linewidth=1.3,
    )

    recovery_axis.set(
        xlim=(-0.01, 0.46),
        ylim=(-0.01, 0.46),
        xlabel="True Left mixing",
        ylabel="Recovered Left mixing",
        title="Simulation verification",
    )

    style_axis(recovery_axis)

    value_axis = figure.add_axes(
        [
            0.89,
            0.16,
            0.10,
            0.61,
        ]
    )

    value_axis.axis("off")

    value_axis.text(
        0.5,
        0.88,
        "Final fit",
        ha="center",
        fontsize=11,
        fontweight="bold",
        color=MID,
    )

    value_axis.text(
        0.5,
        0.61,
        (
            f"{100 * inputs['c_left']:.2f}%"
        ),
        ha="center",
        fontsize=23,
        fontweight="bold",
        color=BLUE,
    )

    value_axis.text(
        0.5,
        0.51,
        "Left",
        ha="center",
        fontsize=9,
        color=MID,
    )

    value_axis.text(
        0.5,
        0.28,
        (
            f"{100 * inputs['c_right']:.2f}%"
        ),
        ha="center",
        fontsize=23,
        fontweight="bold",
        color=ORANGE,
    )

    value_axis.text(
        0.5,
        0.18,
        "Right",
        ha="center",
        fontsize=9,
        color=MID,
    )

    figure.savefig(
        path,
        dpi=300,
        bbox_inches="tight",
        facecolor="white",
    )

    plt.close(figure)


def validate_only():
    inputs = load_inputs(
        load_counts=False
    )

    print(
        "VALIDATION_ONLY=PASS"
    )

    print(
        f"observed_fit="
        f"{inputs['observed_path']}"
    )

    print(
        "expanded_simulations="
        f"{len(inputs['simulations'])}"
    )

    print(
        "sibling_simulations="
        f"{len(inputs['siblings'])}"
    )

    print(
        f"c_left="
        f"{inputs['c_left']:.8f}"
    )

    print(
        f"c_right="
        f"{inputs['c_right']:.8f}"
    )


def main():
    FIGURES.mkdir(
        parents=True,
        exist_ok=True,
    )

    inputs = load_inputs(
        load_counts=True
    )

    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.titleweight": "bold",
            "axes.edgecolor": "#777777",
            "axes.linewidth": 0.7,
            "figure.facecolor": "white",
        }
    )

    with PdfPages(PDF) as pdf:
        page_inputs(
            pdf,
            inputs,
        )

        page_mechanism(
            pdf,
            inputs,
        )

        page_value(
            pdf,
            inputs,
        )

        page_simulations(
            pdf,
            inputs,
        )

        page_siblings(
            pdf,
            inputs,
        )

    summary_figure(
        inputs,
        PNG,
    )

    expanded_metrics = {
        side: recovery_metrics(
            inputs["simulations"],
            side,
        )
        for side in (
            "left",
            "right",
        )
    }

    sibling_exact = (
        inputs["siblings"].loc[
            inputs["siblings"]["scenario"]
            == "exact_20"
        ]
    )

    sibling_metrics = {
        side: recovery_metrics(
            sibling_exact,
            side,
        )
        for side in (
            "left",
            "right",
        )
    }

    audit = {
        "status": "PASS",
        "scope": (
            "mixing-number mechanism only"
        ),
        "observed_sites": 1000000,
        "observed_fit": str(
            inputs["observed_path"]
        ),
        "emission_fit": str(
            inputs["emission_path"]
        ),
        "expanded_simulations": str(
            inputs["simulations_path"]
        ),
        "empirical_sibling_simulations": str(
            inputs["sibling_path"]
        ),
        "c_left": inputs["c_left"],
        "c_right": inputs["c_right"],
        "expanded_recovery": (
            expanded_metrics
        ),
        "sibling_exact_20_recovery": (
            sibling_metrics
        ),
        "pages": 5,
        "pdf": str(PDF),
        "summary_png": str(PNG),
        "focal_GT_used": False,
        "model_refitted": False,
        "fundamental_method_changed": False,
        "note": (
            "Read-only narrative generated "
            "from existing validated fits "
            "and simulations."
        ),
    }

    AUDIT.write_text(
        json.dumps(
            audit,
            indent=2,
        )
        + "\n"
    )

    print(
        "MIXING NUMBER NARRATIVE COMPLETE"
    )

    print(
        "LEFT_MIXING="
        f"{inputs['c_left']:.8f}"
    )

    print(
        "RIGHT_MIXING="
        f"{inputs['c_right']:.8f}"
    )

    print(
        f"PDF={PDF.resolve()}"
    )

    print(
        f"PNG={PNG.resolve()}"
    )

    print(
        f"AUDIT={AUDIT.resolve()}"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--validate-only",
        action="store_true",
    )

    arguments = parser.parse_args()

    if arguments.validate_only:
        validate_only()
    else:
        main()
