#!/usr/bin/env python3
"""Shared plotting utilities for manuscript Figure 2."""

from __future__ import annotations

import json
import math
import os
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LogNorm
from matplotlib.patches import Ellipse
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter

DEFAULT_PROJECT_ROOT = Path(
    os.environ.get(
        "GYNANDRO_ROOT",
        "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling",
    )
)

FIGURE2_RELATIVE_OUTDIR = Path("analysis/figure2/panels")

LEFT_T2 = "#5A4FA3"
LEFT_T3 = "#8A78C2"
RIGHT_T2 = "#E56B00"
RIGHT_T3 = "#F3A038"
GT_COLOURS = {
    0: "#2C9BCB",
    1: "#F28E2B",
    2: "#B61F72",
}
MODEL_COLOUR = "#F000B8"
HAPLOTYPE_COLOURS = {
    "A": "#0072B2",
    "B": "#E69F00",
    "C": "#CC79A7",
    "D": "#009E73",
}


def set_publication_style(font_size: float = 8.0) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": font_size,
            "axes.labelsize": font_size,
            "axes.titlesize": font_size + 1,
            "axes.linewidth": 0.7,
            "xtick.labelsize": font_size - 0.5,
            "ytick.labelsize": font_size - 0.5,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
        }
    )


def project_root(value: str | Path | None = None) -> Path:
    root = DEFAULT_PROJECT_ROOT if value is None else Path(value)
    return root.expanduser().resolve()


def panel_output(root: Path, panel: str, explicit: str | Path | None = None) -> Path:
    if explicit is not None:
        path = Path(explicit)
        if not path.is_absolute():
            path = root / path
    else:
        path = root / FIGURE2_RELATIVE_OUTDIR / f"Fig2{panel}.pdf"
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def add_panel_letter(fig, letter: str, x: float = 0.01, y: float = 0.99, size: float = 14) -> None:
    fig.text(
        x,
        y,
        letter,
        ha="left",
        va="top",
        fontsize=size,
        fontweight="bold",
    )


def load_luka_bleute(root: Path) -> ListedColormap:
    palette_path = root / "palettes/barfbow/luka_bleute.json"
    if palette_path.is_file():
        colours = json.loads(palette_path.read_text())
        if not isinstance(colours, list) or len(colours) < 2:
            raise RuntimeError(f"Invalid Luka-Bleute palette: {palette_path}")
        return ListedColormap(colours, name="luka_bleute")

    fallback = [
        "#f7d4d9",
        "#d9c56c",
        "#9da75b",
        "#5b9668",
        "#218b78",
        "#126b78",
        "#124d70",
        "#172d59",
    ]
    print(
        "WARNING\tLuka-Bleute palette JSON was not found; using bundled fallback",
        flush=True,
    )
    return mcolors.LinearSegmentedColormap.from_list(
        "luka_bleute_fallback", fallback, N=256
    )


def loess_like(x: np.ndarray, y: np.ndarray, span: float = 0.20) -> np.ndarray:
    """Small dependency-free local-linear smoother used in the existing AF plots."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or x.shape != y.shape:
        raise ValueError("x and y must be equal-length one-dimensional arrays")
    if len(x) < 3:
        return np.maximum(y, 0.0)

    neighbours = min(len(x), max(5, int(np.ceil(span * len(x)))))
    fitted = np.empty_like(y, dtype=float)

    for index, centre in enumerate(x):
        distance = np.abs(x - centre)
        bandwidth = np.partition(distance, neighbours - 1)[neighbours - 1]
        bandwidth = max(float(bandwidth) * 1.000001, 1e-12)
        weights = np.maximum(0.0, 1.0 - (distance / bandwidth) ** 3) ** 3
        design = np.column_stack([np.ones(len(x)), x - centre])
        root_weights = np.sqrt(weights)
        coefficients = np.linalg.lstsq(
            design * root_weights[:, None],
            y * root_weights,
            rcond=None,
        )[0]
        fitted[index] = coefficients[0]

    return np.maximum(fitted, 0.0)


def expected_mixed_af(
    p_left: float,
    p_right: float,
    c_left: float = 0.20,
    c_right: float | None = None,
    error: float = 0.001,
) -> tuple[float, float]:
    if c_right is None:
        c_right = c_left
    q_left = (1.0 - c_left) * p_left + c_left * p_right
    q_right = c_right * p_left + (1.0 - c_right) * p_right
    q_left = error + (1.0 - 2.0 * error) * q_left
    q_right = error + (1.0 - 2.0 * error) * q_right
    return float(q_left), float(q_right)


def _beta_binomial_fraction_variance(q: float, depth: float, concentration: float) -> float:
    """Variance of K/n for beta-binomial mean q and alpha+beta=concentration."""
    q = float(np.clip(q, 0.0, 1.0))
    depth = float(depth)
    concentration = float(concentration)
    if depth <= 0 or concentration <= 0:
        raise ValueError("depth and concentration must be positive")
    return q * (1.0 - q) * (concentration + depth) / (
        depth * (concentration + 1.0)
    )


def add_predictive_region(
    ax,
    q_left: float,
    q_right: float,
    *,
    depth_left: float = 40.0,
    depth_right: float = 40.0,
    concentration: float = 22.0,
    mass: float = 0.50,
    colour: str = MODEL_COLOUR,
    linewidth: float = 1.35,
    alpha: float = 0.78,
    skip_diagonal_corners: bool = True,
) -> None:
    """
    Draw the smooth predictive-ring convention used for Figure 2.

    Spread uses beta-binomial sampling variance. At a fixed boundary one
    dimension collapses to a line; at (0,0)/(1,1) no ring is drawn.
    """
    ql = float(q_left)
    qr = float(q_right)

    if skip_diagonal_corners:
        at_00 = ql <= 0.002 and qr <= 0.002
        at_11 = ql >= 0.998 and qr >= 0.998
        if at_00 or at_11:
            return

    var_l = _beta_binomial_fraction_variance(ql, depth_left, concentration)
    var_r = _beta_binomial_fraction_variance(qr, depth_right, concentration)
    sd_l = math.sqrt(max(var_l, 0.0))
    sd_r = math.sqrt(max(var_r, 0.0))

    r2 = -2.0 * math.log(max(1e-12, 1.0 - mass))
    radial = math.sqrt(r2)
    z1 = 0.6744897501960817

    eps = 2e-3
    fixed_l = ql <= eps or ql >= 1.0 - eps or sd_l < 1e-8
    fixed_r = qr <= eps or qr >= 1.0 - eps or sd_r < 1e-8

    if fixed_l and fixed_r:
        return

    if fixed_l:
        half = z1 * sd_r
        ax.plot(
            [ql, ql],
            [max(0.0, qr - half), min(1.0, qr + half)],
            color=colour,
            linewidth=linewidth,
            alpha=alpha,
            zorder=8,
        )
        return

    if fixed_r:
        half = z1 * sd_l
        ax.plot(
            [max(0.0, ql - half), min(1.0, ql + half)],
            [qr, qr],
            color=colour,
            linewidth=linewidth,
            alpha=alpha,
            zorder=8,
        )
        return

    ellipse = Ellipse(
        (ql, qr),
        width=2.0 * radial * sd_l,
        height=2.0 * radial * sd_r,
        facecolor="none",
        edgecolor=colour,
        linewidth=linewidth,
        alpha=alpha,
        zorder=8,
    )
    ax.add_patch(ellipse)


def read_smu(path: Path) -> dict[str, np.ndarray | int]:
    if not path.is_file():
        raise FileNotFoundError(path)
    table = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["coverage_A", "coverage_B", "pairs"],
    )
    for column in table.columns:
        table[column] = pd.to_numeric(table[column], errors="raise")
    table = table.loc[
        (table["coverage_A"] >= 0)
        & (table["coverage_B"] >= 0)
        & (table["pairs"] > 0)
    ].copy()

    a = table["coverage_A"].to_numpy(dtype=float)
    b = table["coverage_B"].to_numpy(dtype=float)
    weight = table["pairs"].to_numpy(dtype=float)
    total = a + b
    keep = np.isfinite(total) & np.isfinite(weight) & (total > 0) & (weight > 0)
    a = a[keep]
    b = b[keep]
    weight = weight[keep]
    total = total[keep]
    x = np.minimum(a, b) / total
    y = total
    return {
        "x": x,
        "y": y,
        "weights": weight,
        "rows": int(len(table)),
        "total_pairs": int(weight.sum()),
    }


def smudge_density(
    data: dict[str, np.ndarray | int],
    *,
    ymax: float = 120.0,
    x_bins: int = 180,
    y_bins: int = 220,
    sigma: float = 1.15,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    density, x_edges, y_edges = np.histogram2d(
        np.asarray(data["x"], dtype=float),
        np.asarray(data["y"], dtype=float),
        bins=[x_bins, y_bins],
        range=[[0.0, 0.5], [0.0, ymax]],
        weights=np.asarray(data["weights"], dtype=float),
    )
    density = gaussian_filter(density.T, sigma=sigma)
    return density, x_edges, y_edges


def shared_smudge_scale(densities: Sequence[np.ndarray]) -> tuple[float, float, np.ndarray]:
    positives = [d[d > 0] for d in densities if np.any(d > 0)]
    if not positives:
        raise RuntimeError("No positive Smudgeplot density")
    values = np.concatenate(positives)
    vmax = float(values.max())
    vmin = max(1.0, float(np.quantile(values, 0.04)))
    if vmin >= vmax:
        vmin = max(1.0, vmax / 100.0)
    levels = np.geomspace(vmin, vmax, 32)
    return vmin, vmax, levels


def draw_smudge_axis(
    ax,
    density: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    *,
    vmin: float,
    vmax: float,
    levels: np.ndarray,
    ymax: float,
    title: str,
    cmap=None,
):
    if cmap is None:
        cmap = plt.get_cmap("viridis").copy()
        cmap.set_under("white")
    x_centres = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centres = 0.5 * (y_edges[:-1] + y_edges[1:])
    artist = ax.contourf(
        x_centres,
        y_centres,
        density,
        levels=levels,
        norm=LogNorm(vmin=vmin, vmax=vmax),
        cmap=cmap,
        extend="max",
        antialiased=True,
    )
    ax.set_xlim(0, 0.5)
    ax.set_ylim(0, ymax)
    ax.set_title(title, fontweight="bold", pad=5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(length=2.5, pad=2)
    return artist


def load_simulation_pair(path: Path, prefix: str) -> tuple[np.ndarray, np.ndarray]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with np.load(path, allow_pickle=False) as data:
        left_key = f"{prefix}_left"
        right_key = f"{prefix}_right"
        if left_key not in data or right_key not in data:
            raise KeyError(f"{path}: missing {left_key}/{right_key}")
        left = np.asarray(data[left_key], dtype=float).copy()
        right = np.asarray(data[right_key], dtype=float).copy()
    if left.shape != right.shape or left.ndim != 1 or len(left) == 0:
        raise RuntimeError(f"Invalid simulation pair in {path}: {left.shape}/{right.shape}")
    return left, right


def validate_fraction_arrays(*arrays: Iterable[float]) -> None:
    for values in arrays:
        arr = np.asarray(values, dtype=float)
        if not np.isfinite(arr).all() or ((arr < 0) | (arr > 1)).any():
            raise RuntimeError("Allele-fraction array contains invalid values")
