#!/usr/bin/env python3
"""Generate Arnaud's separate supplementary PDFs on Pegasus.

Require the original empirical paired AD depths; never infer 10% simulations
by interpolating noisy AF arrays or substitute FORMAT/DP.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import ListedColormap, LogNorm
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figure2"))
from common import (set_publication_style, expected_mixed_af,
                    add_historical_predictive_region, informative_density_vmax)

MODELS = [("2n", "AB", "AB"), ("2n/n", "AB", "A"),
          ("3n (AB/CC)", "AB", "CC"), ("3n (AB/AC)", "AB", "AC"),
          ("4n", "AB", "CD")]
BASE = Path("analysis/figures/pairwise_AF_cli")

def states(left, right):
    letters = sorted(set(left + right))
    pairs = []
    for values in itertools.product((0., 1.), repeat=len(letters)):
        alleles = dict(zip(letters, values))
        pairs.append((sum(alleles[h] for h in left)/len(left),
                      sum(alleles[h] for h in right)/len(right)))
    return np.array(pairs)

def axis_style(ax):
    ax.set(xlim=(-.018, 1.018), ylim=(-.022, 1.022), aspect="equal",
           xticks=[0, .5, 1], yticks=[0, .5, 1])
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(length=3)

def overlays(ax, pairs, mixing, colour="#F000B8"):
    for l, r in np.unique(pairs, axis=0):
        ql, qr = expected_mixed_af(l, r, mixing, error=.001)
        add_historical_predictive_region(ax, ql, qr, colour=colour,
                                        linewidth=.70, mass=.5)
        ax.plot(ql, qr, "o", color=colour, ms=3.5, mec="white", mew=.5,
                zorder=9, clip_on=False)

def scale(counts, x, y):
    compressed = informative_density_vmax(counts, x, y)
    return LogNorm(1, float(np.sqrt(compressed * np.max(counts))), clip=True)

def render_fitted_mixing(args, root, out):
    fit_path = Path(args.fit_json)
    if not fit_path.is_absolute(): fit_path = root / fit_path
    fit = json.loads(fit_path.read_text())
    def fitted_value(explicit, aliases):
        if explicit:
            value = fit
            for key in explicit.split("."): value = value[key]
            return float(value), explicit
        matches = []
        def visit(value, path=""):
            if isinstance(value, dict):
                for key, item in value.items():
                    here = (path + "." + key).strip(".")
                    if key in aliases and isinstance(item, (int, float)):
                        matches.append((float(item), here))
                    visit(item, here)
        visit(fit)
        if len(matches) != 1:
            raise ValueError(f"Need explicit fitted parameter JSON path; candidates={matches}")
        return matches[0]
    cl, kl = fitted_value(args.fit_left_key, {"cL", "c_left", "c_left_from_right", "estimated_c_left_from_right"})
    cr, kr = fitted_value(args.fit_right_key, {"cR", "c_right", "c_right_from_left", "estimated_c_right_from_left"})
    if not (0 <= cl < .5 and 0 <= cr < .5):
        raise ValueError("Unexpected reciprocal mixing estimate; verify saved-fit parameter conventions")
    fig = plt.figure(figsize=(9.4,5.2))
    ax = fig.add_axes([.11,.19,.32,.64])
    values = [100*cl,100*cr]
    ax.bar([0,1],values,color=["#5A4FA3","#E56B00"],width=.55)
    ax.set(xticks=[0,1],xticklabels=["Left\nreads from Right","Right\nreads from Left"],
           ylabel="Fitted reciprocal read mixing (%)",ylim=(0,max(values)*1.3))
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    for x,v in enumerate(values): ax.text(x,v+.6,f"{v:.8f}%",ha="center",fontsize=10)
    fig.text(.50,.86,"Saved directional mixing fit",weight="bold",fontsize=12)
    fig.text(.50,.79,f"{fit.get('sites', 'Unknown'):,} selected sites" if isinstance(fit.get('sites'), int) else "Saved fit",fontsize=9)
    fig.text(.50,.69,f"Left:  cL = {cl:.10f}\nRight: cR = {cr:.10f}",fontsize=10,linespacing=1.7)
    fig.text(.50,.48,"Expected allele fractions before sequencing effects:",fontsize=9)
    fig.text(.50,.36,rf"$m_L = {1-cl:.10f}p_L + {cl:.10f}p_R$",fontsize=11)
    fig.text(.50,.25,rf"$m_R = {cr:.10f}p_L + {1-cr:.10f}p_R$",fontsize=11)
    emission = fit.get("emission", {})
    if emission:
        fig.text(.50,.15, "Saved sequencing model: " +
            f"error={emission['error_rate']:.8f}\nbias log-odds={emission['bias_log_odds']:.8f}; concentration={emission['concentration']:.6f}",fontsize=8)
    fig.text(.11,.065,"Scope: " + fit.get("scope", "unspecified"),fontsize=8)
    fig.text(.11,.025,"Saved pilot fit; no refitting. Display precision does not imply estimation precision.",fontsize=8)
    fig.savefig(out / "03_fitted_left_right_mixing_calculation.pdf")
    plt.close(fig)
    (out / "fitted_mixing_values.json").write_text(json.dumps({"source":str(fit_path),
        "source_sha256":hashlib.sha256(fit_path.read_bytes()).hexdigest(),
        "left_parameter_path":kl,"right_parameter_path":kr,"c_left":cl,"c_right":cr, "sites":fit.get("sites"), "scope":fit.get("scope"), "emission":fit.get("emission")},indent=2))
    return fit_path

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--project-root", default="/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling")
    ap.add_argument("--depth-pairs", default=str(BASE / "AF_narrative_frames/simulated_founder_empirical_depth/paired_depth_sample.npz"),
                    help="Original empirical paired AD-depth NPZ or TSV; NOT FORMAT/DP")
    ap.add_argument("--left-depth-key", default=None)
    ap.add_argument("--right-depth-key", default=None)
    ap.add_argument("--only-fitted-mixing", action="store_true", help="Plot only the saved directional mixing fit, without reading simulation inputs")
    ap.add_argument("--fit-json", default="analysis/genomewide_mixing_corrected_v1/restored_pilot/observed.json")
    ap.add_argument("--fit-left-key", default=None, help="Exact dotted JSON path for fitted Left mixing")
    ap.add_argument("--fit-right-key", default=None, help="Exact dotted JSON path for fitted Right mixing")
    ap.add_argument("--generate-missing", action="store_true", help="Explicitly permit simulations absent from existing caches")
    ap.add_argument("--sites", type=int, default=200000)
    ap.add_argument("--seed", type=int, default=20261012)
    ap.add_argument("--output-dir", default="analysis/arnaud_mixing_supplement")
    args = ap.parse_args()
    root = Path(args.project_root)
    if args.only_fitted_mixing:
        out = Path(args.output_dir)
        if not out.is_absolute(): out = root / out
        out.mkdir(parents=True, exist_ok=True)
        set_publication_style(8)
        render_fitted_mixing(args, root, out)
        print(f"WROTE saved-fit calculation PDF: {out}")
        return
    source = Path(args.depth_pairs)
    if not source.is_absolute(): source = root / source
    if source.suffix == ".npz":
        with np.load(source, allow_pickle=False) as d:
            if args.left_depth_key is None and args.right_depth_key is None:
                depths = d["depths"]
                if depths.ndim != 2 or depths.shape[1] != 2:
                    raise ValueError("Expected original depths cache with shape (N, 2)")
                dl, dr = depths[:, 0].copy(), depths[:, 1].copy()
            else:
                dl, dr = d[args.left_depth_key].copy(), d[args.right_depth_key].copy()
    else:
        d = pd.read_csv(source, sep="\t")
        dl, dr = d[args.left_depth_key].to_numpy(), d[args.right_depth_key].to_numpy()
    dl, dr = np.asarray(dl, dtype=float), np.asarray(dr, dtype=float)
    if dl.ndim != 1 or dl.shape != dr.shape or not len(dl):
        raise ValueError("Empirical depth pairs must be nonempty matched vectors")
    if not (np.isfinite(dl).all() and np.isfinite(dr).all()):
        raise ValueError("Nonfinite depth pairs")
    if ((dl < 1) | (dr < 1) | (dl != np.floor(dl)) | (dr != np.floor(dr))).any():
        raise ValueError("Depths must be positive integer allele-assigned read counts")
    if args.sites < 1: raise ValueError("--sites must be positive")
    palette = root / "palettes/barfbow/luka_bleute.json"
    cmap = ListedColormap(json.loads(palette.read_text()), name="luka_bleute")
    bins_path = root / BASE / "TN_merged_control_sites_hex40_batlowW.AF_bins.tsv.gz"
    bins = pd.read_csv(bins_path, sep="\t")
    if int(bins["count"].sum()) != 31126758:
        raise ValueError("Observed AF site scope differs from approved 31,126,758 sites")
    out = Path(args.output_dir)
    if not out.is_absolute(): out = root / out
    out.mkdir(parents=True, exist_ok=True)
    set_publication_style(8)
    rng = np.random.default_rng(args.seed)
    sampled = rng.integers(len(dl), size=args.sites)
    depths = (dl[sampled].astype(int), dr[sampled].astype(int))
    # Load complete existing model panels without re-extracting VCFs or re-simulating.
    narrative = root / BASE / "AF_narrative_frames"
    sources = [
        (narrative / "simulated_founder_empirical_depth/simulated_AF.npz", "model2", False),
        (narrative / "AB_A_digametic/AB_A_simulated_AF.npz", "", False),
        (narrative / "gyno_three_founder/gyno_simulated_AF.npz", "", True),
        (narrative / "simulated_founder_empirical_depth/simulated_AF.npz", "model3", False),
        (narrative / "simulated_founder_empirical_depth/simulated_AF.npz", "model4", False),
    ]
    simulations = []
    panel_provenance = []
    for (label, left, right), (cache_path, suffix, transpose) in zip(MODELS, sources):
        pairs = states(left, right)
        chosen = pairs[rng.integers(len(pairs), size=args.sites)]
        if cache_path.is_file():
            with np.load(cache_path, allow_pickle=False) as cache:
                if "true_left" in cache and "true_right" in cache:
                    chosen = np.column_stack([cache["true_left"], cache["true_right"]])
                    if transpose: chosen = chosen[:, ::-1]
                    if chosen.shape != (args.sites, 2):
                        raise ValueError("Saved founder-state count differs from --sites")
        for mixing in (0., .1, .2):
            prefix = f"mix{round(mixing*100)}" + ("_" + suffix if suffix else "")
            missing_path = out / ("missing_" + left + "_" + right + "_" + prefix + ".npz")
            signature = json.dumps({"version":2,"seed":args.seed,"sites":args.sites,
                "depth_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),
                "model":label,"mixing":mixing},sort_keys=True)
            af = None
            provenance = None
            for candidate in (cache_path, missing_path):
                if candidate.is_file():
                    with np.load(candidate, allow_pickle=False) as cache:
                        if candidate == missing_path and ("signature" not in cache or str(cache["signature"]) != signature):
                            raise ValueError(f"Generated cache provenance differs: {candidate}")
                        keys = (prefix + "_left", prefix + "_right")
                        if all(k in cache for k in keys):
                            af = [cache[k].copy() for k in keys]
                            if transpose and candidate == cache_path: af = af[::-1]
                            provenance = str(candidate)
                            break
            if af is None:
                if not args.generate_missing:
                    available = []
                    if cache_path.is_file():
                        with np.load(cache_path, allow_pickle=False) as cache: available = cache.files
                    raise RuntimeError(f"Missing {label} {mixing:.0%} cache: {cache_path}; wanted {prefix}_left/right; available={available}. No simulations generated. Inspect caches before --generate-missing.")
                means = .001 + .998 * np.column_stack([
                    (1-mixing)*chosen[:, 0] + mixing*chosen[:, 1],
                    mixing*chosen[:, 0] + (1-mixing)*chosen[:, 1]])
                af = [rng.binomial(dep, rng.beta(means[:, i]*22,
                                                (1-means[:, i])*22))/dep
                      for i, dep in enumerate(depths)]
                # The persisted missing-panel cache is already in requested axis orientation.
                np.savez_compressed(missing_path, signature=np.array(signature), **{prefix+"_left":af[0], prefix+"_right":af[1]})
                provenance = str(missing_path)
            if af[0].ndim != 1 or af[0].shape != af[1].shape or len(af[0]) != args.sites:
                raise ValueError(f"{label} {mixing:.0%}: expected {args.sites} matching sites; adjust --sites to existing cache count")
            if any(not np.isfinite(a).all() or ((a<0)|(a>1)).any() for a in af):
                raise ValueError("Invalid cached AF values")
            simulations.append((label, pairs, mixing, af))
            panel_provenance.append({"model":label,"mixing":mixing,"source":provenance})
    temporary, ax = plt.subplots()
    counts, xs, ys = [], [], []
    for _, _, _, af in simulations:
        h = ax.hexbin(*af, gridsize=40, extent=(0,1,0,1), mincnt=1)
        counts.extend(h.get_array()); xs.extend(h.get_offsets()[:,0]); ys.extend(h.get_offsets()[:,1])
        h.remove()
    plt.close(temporary)
    norm = scale(np.array(counts), np.array(xs), np.array(ys))
    fig, axes = plt.subplots(5, 3, figsize=(8.4,12.6), sharex=True, sharey=True)
    fig.subplots_adjust(left=.22, right=.88, bottom=.07, top=.94, hspace=.15, wspace=.12)
    for index, (ax, (label, pairs, mixing, af)) in enumerate(zip(axes.flat, simulations)):
        h = ax.hexbin(*af, gridsize=40, extent=(0,1,0,1), mincnt=1,
                      cmap=cmap, norm=norm, linewidths=0, rasterized=False)
        axis_style(ax); overlays(ax, pairs, mixing)
        if index < 3: ax.set_title(f"{mixing:.0%} mixing", pad=10, weight="bold")
        if index % 3 == 0:
            left, right = MODELS[index//3][1:]
            ax.text(-.38,.5, label+"\n"+left+" / "+right, transform=ax.transAxes,
                    ha="right", va="center", weight="bold", fontsize=10)
    fig.supxlabel("Left ALT-read fraction", y=.04)
    fig.supylabel("Right ALT-read fraction", x=.02)
    cb = fig.colorbar(h, cax=fig.add_axes([.92,.30,.018,.40]), extend="max")
    cb.set_label("Simulated sites / hexagon")
    fig.text(.22,.012,"Dots: predicted centroids; outlines: 50% predictive regions", fontsize=8)
    fig.savefig(out / "01_five_BAF_models_0_10_20.pdf")
    plt.close(fig)
    unit = np.array([[0,-.5],[.5,-.25],[.5,.25],[0,.5],[-.5,.25],[-.5,-.25]])
    occupied = bins[bins["count"]>0]
    vertices = occupied[["left_centre","right_centre"]].to_numpy()[:,None,:] + unit[None,:,:]*occupied[["width","height"]].to_numpy()[:,None,:]
    observed_norm = scale(occupied["count"].to_numpy(), occupied["left_centre"].to_numpy(), occupied["right_centre"].to_numpy())
    fig, axes = plt.subplots(1,5, figsize=(14,3.7), sharex=True, sharey=True)
    fig.subplots_adjust(left=.055,right=.93,bottom=.20,top=.88,wspace=.12)
    for ax, mixing in zip(axes, (0.,.1,.2,.3,.4)):
        h = PolyCollection(vertices, array=occupied["count"].to_numpy(dtype=float),
                           cmap=cmap, norm=observed_norm, edgecolors="none", rasterized=False)
        ax.add_collection(h); axis_style(ax)
        overlays(ax, states("AB","CD"), mixing)
        ax.set_title(f"{mixing:.0%} mixing", weight="bold")
    fig.supxlabel("TN Left ALT-read fraction", y=.10)
    axes[0].set_ylabel("TN Right ALT-read fraction")
    cb=fig.colorbar(h,cax=fig.add_axes([.95,.26,.012,.50]),extend="max")
    cb.set_label("Observed sites / hexagon")
    fig.text(.055,.035,"Same 31,126,758 observed sites in every panel; illustrative centroids and 50% predictive regions.",fontsize=8)
    fig.savefig(out / "02_mixing_sensitivity_centroid_heatmaps.pdf")
    plt.close(fig)
    fit_path = render_fitted_mixing(args, root, out)
    (out / "figure_legends.txt").write_text(
        "Figure 1. Five founder-composition models at 0%, 10%, and 20% symmetric reciprocal read mixing. "
        "Row labels follow the requested shorthand; haplotype strings specify each side. "
        "AB/CC is plotted in the requested left/right orientation (the earlier AA/BC version is transposed). "
        f"Each panel contains {args.sites:,} simulated sites, independent founder ALT alleles at probability 0.5, "
        "empirical paired allele-assigned read depths, beta-binomial concentration 22, sequencing error 0.001, and no allelic bias. "
        "Existing simulation arrays are reused; any missing panels are independently sampled under the same assumptions. "
        "All panels share hexagonal binning and a saturated logarithmic Luka-Bleute scale.\n\n"
        "Figure 2. Sensitivity of predicted centroid geometry to 0-40% symmetric reciprocal read mixing. "
        "The observed density is unchanged across panels, with no additional filtering. "
        "This is a visual sensitivity comparison, not a likelihood scan or model-selection test.\n\n"
        "In both figures, dots indicate expected ALT-read fractions; outline-only 50% predictive regions "
        "use the established beta-binomial convention at 40 reads per side and concentration 22. "
        "Regions describe read-fraction variation, not confidence in centroid estimates. "
        "Figure 3 reports the saved pilot directional mixing estimates on the previous selected-site bridge scope (not calibrated founder support) and substitutes them into the reciprocal mixing equations before sequencing error and allelic bias. No refitting is performed. Displayed decimal places are not a measure of uncertainty. Mixing refers to the fraction of reads contributed by the opposite lineage, not cell proportions. "
        "Text is embedded as TrueType; hexagons and overlays are vector graphics.\n")
    audit = vars(args) | {"depth_source":str(source),
        "depth_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),
        "palette_sha256":hashlib.sha256(palette.read_bytes()).hexdigest(),
        "observed_bins":str(bins_path), "empirical_depth_pairs":len(dl),
        "simulation_display_vmax":norm.vmax, "observed_display_vmax":observed_norm.vmax,
        "panel_provenance":panel_provenance, "fit_source":str(fit_path), "density_blend":.5, "density_quantile":.995, "origin_exclude":.05}
    (out / "generation_audit.json").write_text(json.dumps(audit,indent=2))
    print(f"WROTE three separate PDFs: {out}")

if __name__ == "__main__":
    main()
