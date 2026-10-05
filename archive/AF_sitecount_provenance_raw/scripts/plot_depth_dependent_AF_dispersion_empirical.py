from pathlib import Path
import argparse
import gzip
import json
import math
import csv

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D


# ---------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------

Z50_1D = 0.67448975
CHI2_50_DF2 = 1.38629436
R50_2D = math.sqrt(CHI2_50_DF2)


# ---------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------

def open_text(path):
    path = Path(path)
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def extract_existing_paths(obj, base_dir):
    """
    Recursively scan JSON-like objects for string values that resolve
    to existing files.
    """
    found = []

    def rec(x):
        if isinstance(x, dict):
            for v in x.values():
                rec(v)
        elif isinstance(x, list):
            for v in x:
                rec(v)
        elif isinstance(x, str):
            candidates = []

            p = Path(x)
            candidates.append(p)

            if not p.is_absolute():
                candidates.append(base_dir / p)

            for c in candidates:
                if c.exists():
                    found.append(c)

    rec(obj)
    return found


def choose_table_from_audit(audit_path):
    audit_path = Path(audit_path)
    if not audit_path.exists():
        raise FileNotFoundError(f"Audit JSON not found: {audit_path}")

    with open(audit_path, "rt") as fh:
        audit = json.load(fh)

    candidates = extract_existing_paths(audit, audit_path.parent)

    # Also add a few obvious neighbours
    for patt in [
        "*.tsv",
        "*.tsv.gz",
        "*.csv",
        "*.csv.gz",
    ]:
        candidates.extend(audit_path.parent.glob(patt))
        candidates.extend(audit_path.parent.parent.glob(patt))

    # unique / preserve order
    seen = set()
    uniq = []
    for c in candidates:
        c = Path(c)
        if c not in seen and c.is_file():
            uniq.append(c)
            seen.add(c)

    # Prefer AF-like tables
    score = []
    for c in uniq:
        name = c.name.lower()
        s = 0
        if "af" in name:
            s += 3
        if "merged" in name:
            s += 2
        if "tn" in name:
            s += 2
        if "site" in name:
            s += 1
        if name.endswith(".tsv") or name.endswith(".tsv.gz"):
            s += 1
        score.append((s, c))

    score.sort(reverse=True)

    for _, c in score:
        try:
            info = inspect_table_header(c)
            if info["has_left"] and info["has_right"]:
                return c, info
        except Exception:
            pass

    raise RuntimeError(
        "Could not auto-detect an empirical AF table with left/right depth "
        f"columns from {audit_path}. Please rerun with --table <path>."
    )


def sniff_delimiter(header_line):
    if header_line.count("\t") >= header_line.count(","):
        return "\t"
    return ","


def inspect_table_header(path):
    with open_text(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("#"):
                continue
            delim = sniff_delimiter(line)
            header = [x.strip() for x in line.rstrip("\n").split(delim)]
            break
        else:
            raise RuntimeError(f"No header found in {path}")

    info = detect_depth_columns(header)
    info["header"] = header
    info["delimiter"] = delim
    return info


def find_first(header, names):
    header_set = set(header)
    for n in names:
        if n in header_set:
            return n
    return None


def find_token_combo(header, must_have):
    out = []
    for h in header:
        hs = h.lower()
        if all(tok in hs for tok in must_have):
            out.append(h)
    if len(out) == 1:
        return out[0]
    return None


def detect_depth_columns(header):
    # Direct depth columns
    left_dp = find_first(header, [
        "left_DP", "left_dp", "Left_DP", "L_DP",
        "TN_Left_DP", "left_depth", "left_total_depth",
        "leftDepth", "LeftDepth",
    ])
    right_dp = find_first(header, [
        "right_DP", "right_dp", "Right_DP", "R_DP",
        "TN_Right_DP", "right_depth", "right_total_depth",
        "rightDepth", "RightDepth",
    ])

    if left_dp is None:
        left_dp = find_token_combo(header, ["left", "dp"])
    if right_dp is None:
        right_dp = find_token_combo(header, ["right", "dp"])

    # Ref/alt fallback
    left_ref = find_first(header, [
        "left_ref_depth", "left_ref_DP", "TN_Left_ref_depth",
        "Left_ref_depth", "leftRefDepth",
    ])
    left_alt = find_first(header, [
        "left_alt_depth", "left_alt_DP", "TN_Left_alt_depth",
        "Left_alt_depth", "leftAltDepth",
    ])
    right_ref = find_first(header, [
        "right_ref_depth", "right_ref_DP", "TN_Right_ref_depth",
        "Right_ref_depth", "rightRefDepth",
    ])
    right_alt = find_first(header, [
        "right_alt_depth", "right_alt_DP", "TN_Right_alt_depth",
        "Right_alt_depth", "rightAltDepth",
    ])

    if left_ref is None:
        left_ref = find_token_combo(header, ["left", "ref", "depth"])
    if left_alt is None:
        left_alt = find_token_combo(header, ["left", "alt", "depth"])
    if right_ref is None:
        right_ref = find_token_combo(header, ["right", "ref", "depth"])
    if right_alt is None:
        right_alt = find_token_combo(header, ["right", "alt", "depth"])

    has_left = (left_dp is not None) or (left_ref is not None and left_alt is not None)
    has_right = (right_dp is not None) or (right_ref is not None and right_alt is not None)

    return {
        "left_dp": left_dp,
        "right_dp": right_dp,
        "left_ref": left_ref,
        "left_alt": left_alt,
        "right_ref": right_ref,
        "right_alt": right_alt,
        "has_left": has_left,
        "has_right": has_right,
    }


def safe_float(x):
    try:
        if x is None:
            return np.nan
        x = str(x).strip()
        if x == "" or x.upper() == "NA":
            return np.nan
        return float(x)
    except Exception:
        return np.nan


def read_empirical_depth_pairs(table_path):
    info = inspect_table_header(table_path)
    delim = info["delimiter"]

    left_vals = []
    right_vals = []

    with open_text(table_path) as fh:
        reader = csv.DictReader(
            (line for line in fh if line.strip() and not line.startswith("#")),
            delimiter=delim
        )

        for row in reader:
            # left
            if info["left_dp"] is not None:
                left = safe_float(row.get(info["left_dp"]))
            else:
                left = (
                    safe_float(row.get(info["left_ref"])) +
                    safe_float(row.get(info["left_alt"]))
                )

            # right
            if info["right_dp"] is not None:
                right = safe_float(row.get(info["right_dp"]))
            else:
                right = (
                    safe_float(row.get(info["right_ref"])) +
                    safe_float(row.get(info["right_alt"]))
                )

            if not np.isfinite(left) or not np.isfinite(right):
                continue

            if left <= 0 or right <= 0:
                continue

            left_vals.append(left)
            right_vals.append(right)

    if len(left_vals) == 0:
        raise RuntimeError(f"No usable left/right depth pairs found in {table_path}")

    return np.asarray(left_vals, dtype=float), np.asarray(right_vals, dtype=float), info


# ---------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------

def af_sd(p, depth):
    if depth <= 0:
        return 0.0
    return math.sqrt(p * (1.0 - p) / depth)


def het_sd(depth):
    return af_sd(0.5, depth)


def discrete_depth_colors(depths, cmap_name="viridis"):
    cmap = plt.get_cmap(cmap_name)
    n = len(depths)
    if n == 1:
        return {depths[0]: cmap(0.7)}
    vals = np.linspace(0.15, 0.90, n)
    return {d: cmap(v) for d, v in zip(depths, vals)}


# ---------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------

def draw_50ci_region(ax, mx, my, depth_x, depth_y, color, lw=0.8, alpha=0.18, zorder=2):
    sx = af_sd(mx, depth_x)
    sy = af_sd(my, depth_y)

    if sx == 0 and sy == 0:
        ax.scatter([mx], [my], s=9, color=color, edgecolor="none", zorder=zorder + 1)
        return

    if sx == 0 and sy > 0:
        half = Z50_1D * sy
        ax.plot([mx, mx], [my - half, my + half],
                color=color, linewidth=lw + 0.5, alpha=min(1.0, alpha + 0.25), zorder=zorder + 1)
        return

    if sy == 0 and sx > 0:
        half = Z50_1D * sx
        ax.plot([mx - half, mx + half], [my, my],
                color=color, linewidth=lw + 0.5, alpha=min(1.0, alpha + 0.25), zorder=zorder + 1)
        return

    width = 2.0 * R50_2D * sx
    height = 2.0 * R50_2D * sy

    e = Ellipse(
        xy=(mx, my),
        width=width,
        height=height,
        angle=0.0,
        facecolor=color,
        edgecolor=color,
        linewidth=lw,
        alpha=alpha,
        zorder=zorder,
    )
    ax.add_patch(e)


# ---------------------------------------------------------------------
# Plot 1
# ---------------------------------------------------------------------

def plot_sd_curve(ax, depths, empirical_site_mean, empirical_color, color_map):
    q1 = np.percentile(empirical_site_mean, 25)
    med = np.percentile(empirical_site_mean, 50)
    q3 = np.percentile(empirical_site_mean, 75)

    max_depth = max(max(depths), q3, med) * 1.15
    max_depth = max(max_depth, 50)
    x = np.arange(5, int(math.ceil(max_depth)) + 5, 1)
    y = np.array([het_sd(d) for d in x])

    ax.plot(x, y, color="black", linewidth=2.0, zorder=2)

    # reference depths
    for d in depths:
        ax.scatter(d, het_sd(d),
                   s=55, color=color_map[d], edgecolor="black", linewidth=0.6, zorder=3)
        ax.text(d, het_sd(d) + 0.003, f"{int(d)}",
                ha="center", va="bottom", fontsize=8, color=color_map[d])

    # empirical IQR band
    ax.axvspan(q1, q3, color=empirical_color, alpha=0.14, zorder=1)

    # empirical median point
    ax.scatter(med, het_sd(med),
               s=95, color=empirical_color, edgecolor="black", linewidth=0.9, zorder=4)

    ax.text(
        med,
        het_sd(med) + 0.006,
        f"empirical median = {med:.1f}×\nIQR = {q1:.1f}–{q3:.1f}×",
        ha="center",
        va="bottom",
        fontsize=9,
        color=empirical_color,
        fontweight="bold",
    )

    ax.set_xlabel("Depth (×)")
    ax.set_ylabel("Expected SD of het AF")
    ax.set_title("Expected heterozygote AF dispersion vs depth", loc="left", fontsize=12)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=3, width=0.8)

    ax.set_xlim(0, max(x) * 1.02)
    ax.set_ylim(0, max(y) * 1.10)


# ---------------------------------------------------------------------
# Plot 2 / 3
# ---------------------------------------------------------------------

def plot_founder_panel(ax, n_founders, depths, empirical_left, empirical_right,
                       empirical_color, color_map, empirical_n_overlay=180):
    vals = np.arange(n_founders + 1) / n_founders

    # base centroids
    for x in vals:
        for y in vals:
            ax.scatter(x, y, s=8, color="black", alpha=0.65, zorder=6)

    # reference fixed depths
    for d in depths:
        for x in vals:
            for y in vals:
                draw_50ci_region(
                    ax=ax,
                    mx=x,
                    my=y,
                    depth_x=d,
                    depth_y=d,
                    color=color_map[d],
                    lw=0.8,
                    alpha=0.20,
                    zorder=1,
                )

    # empirical sample of real paired depths
    rng = np.random.default_rng(1)
    n_sites = len(empirical_left)
    use_n = min(empirical_n_overlay, n_sites)
    idx = rng.choice(n_sites, size=use_n, replace=False)

    for i in idx:
        dx = empirical_left[i]
        dy = empirical_right[i]
        for x in vals:
            for y in vals:
                draw_50ci_region(
                    ax=ax,
                    mx=x,
                    my=y,
                    depth_x=dx,
                    depth_y=dy,
                    color=empirical_color,
                    lw=0.5,
                    alpha=0.035,
                    zorder=3,
                )

    # empirical median pair
    med_left = float(np.median(empirical_left))
    med_right = float(np.median(empirical_right))

    for x in vals:
        for y in vals:
            draw_50ci_region(
                ax=ax,
                mx=x,
                my=y,
                depth_x=med_left,
                depth_y=med_right,
                color=empirical_color,
                lw=1.0,
                alpha=0.13,
                zorder=4,
            )

    ax.set_xlim(-0.04, 1.04)
    ax.set_ylim(-0.04, 1.04)
    ax.set_aspect("equal", adjustable="box")

    ax.set_xlabel("Left ALT-allele fraction")
    ax.set_ylabel("Right ALT-allele fraction")

    ax.set_title(
        f"{n_founders}-founder model: 50% sampling intervals across depth",
        loc="left",
        fontsize=12,
    )

    ax.axvline(0.5, color="#BBBBBB", linestyle="--", linewidth=0.8, zorder=0)
    ax.axhline(0.5, color="#BBBBBB", linestyle="--", linewidth=0.8, zorder=0)

    ticks = [0, 0.25, 0.5, 0.75, 1.0]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=3, width=0.8)

    handles = []
    for d in depths:
        handles.append(
            Line2D([0], [0], color=color_map[d], linewidth=6, alpha=0.8, label=f"{int(d)}×")
        )

    handles.append(
        Line2D([0], [0], color=empirical_color, linewidth=6, alpha=0.85,
               label="empirical paired depths")
    )

    ax.legend(
        handles=handles,
        title="Depth",
        frameon=False,
        fontsize=8,
        title_fontsize=9,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.0,
    )


# ---------------------------------------------------------------------
# Save helpers
# ---------------------------------------------------------------------

def finish_and_save(fig, outbase):
    fig.savefig(f"{outbase}.pdf", bbox_inches="tight")
    fig.savefig(f"{outbase}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"WROTE: {outbase}.pdf")
    print(f"WROTE: {outbase}.png")


def write_summary(outdir, table_path, left_depth, right_depth, info):
    site_mean = 0.5 * (left_depth + right_depth)

    out = outdir / "00_empirical_depth_summary.tsv"
    with open(out, "wt") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"table\t{table_path}\n")
        fh.write(f"n_sites\t{len(site_mean)}\n")
        fh.write(f"left_median\t{np.median(left_depth):.6f}\n")
        fh.write(f"right_median\t{np.median(right_depth):.6f}\n")
        fh.write(f"paired_mean_median\t{np.median(site_mean):.6f}\n")
        fh.write(f"paired_mean_mean\t{np.mean(site_mean):.6f}\n")
        fh.write(f"paired_mean_q25\t{np.percentile(site_mean, 25):.6f}\n")
        fh.write(f"paired_mean_q75\t{np.percentile(site_mean, 75):.6f}\n")
        fh.write(f"left_column_mode\t{'direct_DP' if info['left_dp'] else 'ref_plus_alt'}\n")
        fh.write(f"right_column_mode\t{'direct_DP' if info['right_dp'] else 'ref_plus_alt'}\n")
        fh.write(f"left_dp_column\t{info['left_dp']}\n")
        fh.write(f"right_dp_column\t{info['right_dp']}\n")
        fh.write(f"left_ref_column\t{info['left_ref']}\n")
        fh.write(f"left_alt_column\t{info['left_alt']}\n")
        fh.write(f"right_ref_column\t{info['right_ref']}\n")
        fh.write(f"right_alt_column\t{info['right_alt']}\n")

    print(f"WROTE: {out}")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--audit",
        default="analysis/figures/pairwise_AF_cli/TN_merged_control_sites_hex40_batlowW.AF.audit.json",
        help="Audit JSON from which to auto-detect the empirical TN AF table",
    )
    ap.add_argument(
        "--table",
        default=None,
        help="Optional explicit AF table to use instead of auto-detection",
    )
    ap.add_argument(
        "--depths",
        type=float,
        nargs="+",
        default=[10, 20, 30, 50, 100, 200],
        help="Reference depths to display",
    )
    ap.add_argument(
        "--outdir",
        default="analysis/figures/depth_dispersion_model_empirical",
        help="Output directory",
    )
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if args.table is not None:
        table_path = Path(args.table)
        if not table_path.exists():
            raise FileNotFoundError(f"--table not found: {table_path}")
        info = inspect_table_header(table_path)
        if not (info["has_left"] and info["has_right"]):
            raise RuntimeError(
                f"Could not detect usable left/right depth columns in {table_path}"
            )
    else:
        table_path, info = choose_table_from_audit(args.audit)

    print(f"EMPIRICAL TABLE: {table_path}")

    left_depth, right_depth, info = read_empirical_depth_pairs(table_path)

    site_mean = 0.5 * (left_depth + right_depth)

    print(f"N SITES: {len(site_mean)}")
    print(f"LEFT median depth:  {np.median(left_depth):.3f}")
    print(f"RIGHT median depth: {np.median(right_depth):.3f}")
    print(f"PAIRED mean-depth median: {np.median(site_mean):.3f}")
    print(f"PAIRED mean-depth IQR: {np.percentile(site_mean,25):.3f} - {np.percentile(site_mean,75):.3f}")

    write_summary(outdir, table_path, left_depth, right_depth, info)

    depths = [float(d) for d in args.depths]
    color_map = discrete_depth_colors(depths, cmap_name="viridis")
    empirical_color = "#D62728"

    # 1. SD curve
    fig1, ax1 = plt.subplots(figsize=(6.8, 4.8))
    plot_sd_curve(
        ax=ax1,
        depths=depths,
        empirical_site_mean=site_mean,
        empirical_color=empirical_color,
        color_map=color_map,
    )
    fig1.tight_layout()
    finish_and_save(fig1, outdir / "01_depth_vs_hetAF_sd_empirical")

    # 2. Two-founder
    fig2, ax2 = plt.subplots(figsize=(7.4, 6.0))
    plot_founder_panel(
        ax=ax2,
        n_founders=2,
        depths=depths,
        empirical_left=left_depth,
        empirical_right=right_depth,
        empirical_color=empirical_color,
        color_map=color_map,
    )
    fig2.tight_layout()
    finish_and_save(fig2, outdir / "02_two_founder_depth_50CI_empirical")

    # 3. Four-founder
    fig3, ax3 = plt.subplots(figsize=(7.4, 6.0))
    plot_founder_panel(
        ax=ax3,
        n_founders=4,
        depths=depths,
        empirical_left=left_depth,
        empirical_right=right_depth,
        empirical_color=empirical_color,
        color_map=color_map,
    )
    fig3.tight_layout()
    finish_and_save(fig3, outdir / "03_four_founder_depth_50CI_empirical")

    # 4. Combined triptych
    fig4 = plt.figure(figsize=(17.8, 5.4))
    gs = fig4.add_gridspec(1, 3, width_ratios=[1.05, 1.0, 1.0], wspace=0.42)

    axA = fig4.add_subplot(gs[0, 0])
    axB = fig4.add_subplot(gs[0, 1])
    axC = fig4.add_subplot(gs[0, 2])

    plot_sd_curve(
        ax=axA,
        depths=depths,
        empirical_site_mean=site_mean,
        empirical_color=empirical_color,
        color_map=color_map,
    )

    plot_founder_panel(
        ax=axB,
        n_founders=2,
        depths=depths,
        empirical_left=left_depth,
        empirical_right=right_depth,
        empirical_color=empirical_color,
        color_map=color_map,
    )

    plot_founder_panel(
        ax=axC,
        n_founders=4,
        depths=depths,
        empirical_left=left_depth,
        empirical_right=right_depth,
        empirical_color=empirical_color,
        color_map=color_map,
    )

    fig4.tight_layout()
    finish_and_save(fig4, outdir / "00_depth_dispersion_triptych_empirical")


if __name__ == "__main__":
    main()
