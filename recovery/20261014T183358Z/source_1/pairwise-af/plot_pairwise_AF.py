#!/usr/bin/env python3
"""Stream a sorted multisample VCF/BCF into an auditable pairwise AD histogram."""

import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
from importlib.metadata import version as package_version
import json
from pathlib import Path
import platform
import shlex
import sys
import tempfile
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, ListedColormap
import numpy as np
import pysam

__version__ = "0.2.0"
BASES = set("ACGT")


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--version", action="version", version=__version__)
    p.add_argument("--vcf", required=True, type=Path,
                   help="Coordinate-sorted multisample VCF, BGZF VCF or BCF; no index required")
    p.add_argument("--left-sample", required=True)
    p.add_argument("--right-sample", required=True)
    p.add_argument("--output-prefix", required=True, type=Path)
    p.add_argument("--min-depth", type=int, default=10, help="Minimum REF+ALT AD on EACH side (default 10)")
    p.add_argument("--max-depth", type=int, help="Optional maximum REF+ALT AD on EACH side")
    p.add_argument("--filter", choices=["all", "pass"], default="all",
                   help="all: ignore FILTER; pass: require explicit PASS, exclude '.' (default all)")
    p.add_argument("--sites", type=Path,
                   help="Exact allele whitelist: tab-separated CHROM POS REF ALT; 1-based POS; optionally .gz")
    p.add_argument("--contig", action="append", help="Keep this exact contig name; may be repeated")
    p.add_argument("--bins", type=int, default=65,
                   help="Hexagonal columns or square bins per axis (default 65)")
    p.add_argument("--bin-shape", choices=["hex", "square"], default="hex")
    p.add_argument("--cmap", choices=["batlowW_r", "Greys"], default="batlowW_r")
    p.add_argument("--scale", choices=["count", "proportion"], default="count")
    p.add_argument("--vmin", type=float, help="Positive lower colour limit; pair with --vmax")
    p.add_argument("--vmax", type=float, help="Positive upper colour limit; pair with --vmin")
    p.add_argument("--mixing-left", type=float, help="Optional illustrative centroid overlay, NOT fitted")
    p.add_argument("--mixing-right", type=float, help="Must accompany --mixing-left")
    p.add_argument("--title", default="Pairwise allele-read fractions")
    p.add_argument("--validate-only", action="store_true", help="Check the ENTIRE input; no output files")
    p.add_argument("--progress-every", type=int, default=1_000_000, help="Report every N records; 0 disables")
    return p


def fingerprint(path):
    stat = path.stat()
    return {"path": str(path.resolve()), "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def check_args(args):
    if not args.vcf.is_file():
        raise ValueError(f"Input does not exist: {args.vcf}")
    if args.left_sample == args.right_sample:
        raise ValueError("Left and Right must be different samples")
    if args.min_depth < 1 or (args.max_depth is not None and args.max_depth < args.min_depth):
        raise ValueError("Require min-depth >= 1 and max-depth >= min-depth")
    if not 2 <= args.bins <= 1000:
        raise ValueError("bins must be between 2 and 1000")
    if args.progress_every < 0:
        raise ValueError("progress-every must be nonnegative")
    if (args.vmin is None) != (args.vmax is None):
        raise ValueError("Provide both vmin and vmax, or neither")
    if args.vmin is not None and not (0 < args.vmin < args.vmax < float("inf")):
        raise ValueError("Require finite 0 < vmin < vmax")
    if (args.mixing_left is None) != (args.mixing_right is None):
        raise ValueError("Provide both mixing-left and mixing-right, or neither")
    for mixing in [args.mixing_left, args.mixing_right]:
        if mixing is not None and not 0 <= mixing <= 0.45:
            raise ValueError("Illustrative mixing values must lie between 0 and 0.45")


def load_sites(path):
    if path is None:
        return None
    sites = set()
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for number, line in enumerate(handle, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if [s.lower() for s in fields] == ["chrom", "pos", "ref", "alt"]:
                continue
            if len(fields) != 4:
                raise ValueError(f"{path}:{number}: expected CHROM POS REF ALT (4 tab-separated columns); not BED")
            chrom, pos, ref, alt = fields
            try:
                pos = int(pos)
            except ValueError:
                raise ValueError(f"{path}:{number}: POS must be a 1-based integer") from None
            if pos < 1 or ref not in BASES or alt not in BASES or ref == alt:
                raise ValueError(f"{path}:{number}: invalid biallelic SNP")
            key = (chrom, pos, ref, alt)
            if key in sites:
                raise ValueError(f"{path}:{number}: duplicate whitelist variant {key}")
            sites.add(key)
    if not sites:
        raise ValueError("Site whitelist is empty")
    return sites


def valid_ad(value):
    if value is None or len(value) != 2:
        return False
    return all(isinstance(x, int) and x >= 0 for x in value)


def hex_geometry(nx):
    """Matplotlib linear hexbin geometry on extent (0,1,0,1)."""
    ny = int(nx / np.sqrt(3))
    sx, sy = ((1 + 1e-9) - (-1e-9)) / nx, 1 / ny
    centres = np.concatenate([
        np.column_stack([np.repeat(np.arange(nx + 1), ny + 1),
                         np.tile(np.arange(ny + 1), nx + 1)]),
        np.column_stack([np.repeat(np.arange(nx) + .5, ny),
                         np.tile(np.arange(ny) + .5, nx)]),
    ])
    centres *= [sx, sy]
    centres[:, 0] -= 1e-9
    return centres, sx, sy


def hex_index(left, right, nx):
    # Same nearest-centre metric and tie rule as Axes.hexbin. Streaming:
    # each retained SNP contributes exactly once; no square-to-hex conversion.
    ny = int(nx / np.sqrt(3))
    sx, sy = ((1 + 1e-9) - (-1e-9)) / nx, 1 / ny
    x, y = (left - (-1e-9)) / sx, right / sy
    ix, iy = round(x), round(y)
    jx, jy = int(x // 1), int(y // 1)
    d1 = (x - ix)**2 + 3 * (y - iy)**2
    d2 = (x - jx - .5)**2 + 3 * (y - jy - .5)**2
    if d1 < d2:
        return ix * (ny + 1) + iy
    return (nx + 1) * (ny + 1) + jx * ny + jy


def colourmap(name):
    if name == "batlowW_r":
        from cmcrameri import cm
        base = cm.batlowW_r
        # Skip the near-white first 4% for occupied cells. Pure white is zero.
        cmap = ListedColormap(base(np.linspace(.04, 1, 256)))
    else:
        cmap = ListedColormap(plt.get_cmap("Greys")(np.linspace(.18, 1, 256)))
    cmap.set_bad("white")
    return cmap


def accumulate(args, targets=None):
    histogram = np.zeros(
        len(hex_geometry(args.bins)[0]) if args.bin_shape == "hex" else (args.bins, args.bins),
        dtype=np.int64)
    edges = np.linspace(0.0, 1.0, args.bins + 1)
    skipped = Counter()
    by_contig = Counter()
    seen = retained = matched_targets = 0
    current_contig, previous_pos = None, -1
    closed_contigs, alleles_at_position = set(), set()
    wanted_contigs = set(args.contig or [])
    with pysam.VariantFile(str(args.vcf)) as vcf:
        names = list(vcf.header.samples)
        missing = {args.left_sample, args.right_sample} - set(names)
        if missing:
            raise ValueError(f"Missing sample(s): {sorted(missing)}. Available: {names}")
        if "AD" not in vcf.header.formats:
            raise ValueError("VCF has no FORMAT/AD; cannot calculate read fractions from GT or DP")
        ad_header = vcf.header.formats["AD"]
        if ad_header.type != "Integer" or str(ad_header.number) != "R":
            raise ValueError("FORMAT/AD must be declared Number=R,Type=Integer")
        if wanted_contigs - set(vcf.header.contigs):
            raise ValueError(f"Unknown contig(s): {sorted(wanted_contigs - set(vcf.header.contigs))}")
        # Always address the two samples by NAME, regardless of their header order.
        vcf.subset_samples([args.left_sample, args.right_sample])
        header_digest = hashlib.sha256(str(vcf.header).encode()).hexdigest()
        for record in vcf:
            seen += 1
            if args.progress_every and seen % args.progress_every == 0:
                print(f"READ\t{seen:,}\tRETAINED\t{retained:,}", file=sys.stderr, flush=True)
            chrom, pos = record.contig, record.pos
            if chrom != current_contig:
                if chrom in closed_contigs:
                    raise ValueError(f"Unsorted input: contig {chrom} occurs in multiple blocks")
                if current_contig is not None:
                    closed_contigs.add(current_contig)
                current_contig, previous_pos = chrom, -1
            if pos < previous_pos:
                raise ValueError(f"Unsorted input at {chrom}:{pos}")
            if pos != previous_pos:
                alleles_at_position.clear()
            previous_pos = pos
            alts = record.alts or ()
            if any(a in {"<NON_REF>", "<*>"} for a in alts):
                raise ValueError(f"gVCF symbolic allele at {chrom}:{pos}. Jointly genotype or recall at a common site set first; reference blocks are not per-site AD.")
            if wanted_contigs and chrom not in wanted_contigs:
                skipped["outside_contigs"] += 1
                continue
            if len(alts) != 1 or record.ref not in BASES or alts[0] not in BASES or record.ref == alts[0]:
                skipped["not_biallelic_snp"] += 1
                continue
            allele_pair = (record.ref, alts[0])
            if allele_pair in alleles_at_position:
                raise ValueError(f"Duplicate SNP record: {chrom}:{pos}:{record.ref}:{alts[0]}; refusing to double-count")
            alleles_at_position.add(allele_pair)
            key = (chrom, pos, record.ref, alts[0])
            if targets is not None:
                if key not in targets:
                    skipped["outside_sites"] += 1
                    continue
                matched_targets += 1
            if args.filter == "pass" and list(record.filter) != ["PASS"]:
                skipped["filter_not_pass"] += 1
                continue
            left_ad = record.samples[args.left_sample].get("AD")
            right_ad = record.samples[args.right_sample].get("AD")
            if not valid_ad(left_ad) or not valid_ad(right_ad):
                skipped["missing_or_invalid_ad"] += 1
                continue
            left_depth, right_depth = sum(left_ad), sum(right_ad)
            if min(left_depth, right_depth) < args.min_depth:
                skipped["below_min_depth"] += 1
                continue
            if args.max_depth is not None and max(left_depth, right_depth) > args.max_depth:
                skipped["above_max_depth"] += 1
                continue
            left_af, right_af = left_ad[1] / left_depth, right_ad[1] / right_depth
            if args.bin_shape == "hex":
                histogram[hex_index(left_af, right_af, args.bins)] += 1
            else:
                lx = min(int(np.searchsorted(edges, left_af, side="right") - 1), args.bins - 1)
                ry = min(int(np.searchsorted(edges, right_af, side="right") - 1), args.bins - 1)
                histogram[lx, ry] += 1
            retained += 1
            by_contig[chrom] += 1
    if not retained:
        raise ValueError(f"No usable sites. Read {seen}; exclusions: {dict(skipped)}")
    assert retained == int(histogram.sum())
    assert seen == retained + sum(skipped.values())
    audit = {
        "records_read": seen, "records_retained": retained,
        "excluded_first_reason": dict(skipped), "retained_by_contig": dict(by_contig),
        "selected_header_sha256": header_digest,
        "sites_requested": len(targets) if targets is not None else None,
        "sites_matched_before_filter_and_depth": matched_targets if targets is not None else None,
        "sites_not_matched_in_selected_contigs": len(targets) - matched_targets if targets is not None else None,
    }
    return histogram, audit


def draw(histogram, args, directory):
    total = int(histogram.sum())
    values = histogram.astype(float)
    if args.scale == "proportion":
        values /= total
    positives = values[values > 0]
    vmin = args.vmin if args.vmin is not None else (1 if args.scale == "count" else 1 / total)
    vmax = args.vmax if args.vmax is not None else max(float(positives.max()), vmin * 10)
    edges = np.linspace(0, 1, args.bins + 1)
    fig = plt.figure(figsize=(8.5, 8.1))
    grid = fig.add_gridspec(1, 2, width_ratios=[1, 0.045], left=0.12, right=0.88,
                           bottom=0.24, top=0.84, wspace=0.12)
    ax, cax = fig.add_subplot(grid[0]), fig.add_subplot(grid[1])
    cmap = colourmap(args.cmap)
    norm = LogNorm(vmin=vmin, vmax=vmax)
    if args.bin_shape == "hex":
        im = ax.hexbin([], [], gridsize=args.bins, extent=(0, 1, 0, 1),
                       mincnt=0, cmap=cmap, norm=norm, linewidths=0,
                       edgecolors="none", rasterized=True)
        if len(im.get_array()) != len(values):
            raise ValueError("Hexagonal geometry disagrees with stored counts")
        im.set_array(np.ma.masked_equal(values, 0))
    else:
        im = ax.pcolormesh(edges, edges, np.ma.masked_equal(values.T, 0), cmap=cmap,
                          norm=norm, rasterized=True, shading="flat")
    ax.set_facecolor("white")
    ax.set(xlim=(0, 1), ylim=(0, 1), aspect="equal",
           xlabel=f"{args.left_sample}\nALT / (REF + ALT reads)",
           ylabel=f"{args.right_sample}\nALT / (REF + ALT reads)")
    ax.set_xticks(np.linspace(0, 1, 6))
    ax.set_yticks(np.linspace(0, 1, 6))
    ax.text(0.025, 0.975, f"n = {total:,} SNP records", transform=ax.transAxes,
            va="top", fontsize=10, bbox=dict(facecolor="white", edgecolor="none", alpha=0.85))
    clipped = bool(np.any(positives < vmin) or np.any(positives > vmax))
    below, above = np.any(positives < vmin), np.any(positives > vmax)
    extend = "both" if below and above else "min" if below else "max" if above else "neither"
    label = "SNP records per bin (log colour scale)" if args.scale == "count" else "Fraction of retained SNP records per bin (log colour scale)"
    fig.colorbar(im, cax=cax, label=label.replace("per bin", "per hexagon")
                 if args.bin_shape == "hex" else label, extend=extend)
    cax.text(.5, -.075, "White = 0", transform=cax.transAxes,
             ha="center", va="top", fontsize=9)
    overlay = "No genotype calls or model-based site selection."
    if args.mixing_left is not None:
        for pl in [0, 0.5, 1]:
            for pr in [0, 0.5, 1]:
                ql = (1 - args.mixing_left) * pl + args.mixing_left * pr
                qr = args.mixing_right * pl + (1 - args.mixing_right) * pr
                ax.scatter([ql], [qr], s=48, facecolors="none", edgecolors="#0072B2", clip_on=False)
        overlay = (f"Ideal centroids only: supplied cL={args.mixing_left:g}, cR={args.mixing_right:g}. "
                   "No fitting or emission-bias correction.")
    fig.suptitle(textwrap.fill(args.title, 55), y=0.975, fontsize=17, fontweight="bold")
    selection = "Exact allele whitelist" if args.sites else "All eligible SNP records present in input VCF"
    maximum = f"; maximum {args.max_depth}" if args.max_depth is not None else ""
    geometry = (f"{args.bins} hexagonal columns ({histogram.size:,} cells)"
                if args.bin_shape == "hex" else f"{args.bins} x {args.bins} square bins")
    footer = (f"{selection}. Minimum REF+ALT AD {args.min_depth} on each side{maximum}. "
              f"FILTER policy: {args.filter}. {geometry}.\n"
              f"{overlay}\nSite discovery/ascertainment is inherited from the input; no claim of unbiased genome sampling.")
    footer = "\n".join(textwrap.fill(line, 100) for line in footer.splitlines())
    fig.text(0.12, 0.035, footer, ha="left", va="bottom", fontsize=8.5, color="#444444")
    fig.savefig(directory / "AF.pdf", dpi=220)
    fig.savefig(directory / "AF.png", dpi=220)
    plt.close(fig)
    return {"scale": args.scale, "vmin": vmin, "vmax": vmax, "colour_limits_clip_positive_bins": clipped,
            "grid_resolution": args.bins, "bin_shape": args.bin_shape,
            "cells": histogram.size, "zero_bins": "white; not pseudo-counted",
            "colormap": args.cmap + "; low-end truncated for occupied bins; white reserved for zero",
            "orientation": "x=left, y=right",
            "bin_edges": edges.tolist() if args.bin_shape == "square" else None,
            "hex_extent": [0, 1, 0, 1] if args.bin_shape == "hex" else None}


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        check_args(args)
        colourmap(args.cmap)  # Fail before scanning if palette dependency is missing.
        targets = load_sites(args.sites)
        input_before = fingerprint(args.vcf)
        suffixes = ["AF.pdf", "AF.png", "AF_bins.tsv.gz", "AF.audit.json"]
        destinations = [Path(str(args.output_prefix) + "." + suffix) for suffix in suffixes]
        if not args.validate_only and any(p.exists() for p in destinations):
            raise ValueError("Output already exists. Choose a NEW output-prefix; existing results are never overwritten.")
        histogram, audit = accumulate(args, targets)
        if fingerprint(args.vcf) != input_before:
            raise ValueError("Input changed during the run; refusing to publish outputs")
        audit.update({
            "tool": "pairwise-af", "version": __version__,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "command": shlex.join(["plot-pairwise-af"] + (list(argv) if argv is not None else sys.argv[1:])),
            "input": input_before, "script_sha256": sha256(Path(__file__)),
            "parameters": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
            "versions": {"python": platform.python_version(), "numpy": np.__version__,
                         "matplotlib": matplotlib.__version__, "pysam": pysam.__version__,
                         "cmcrameri": package_version("cmcrameri") if args.cmap == "batlowW_r" else None,
                         "samtools_bundled_with_pysam": pysam.__samtools_version__},
            "af_definition": "AD[1] / (AD[0] + AD[1]); not FORMAT/DP; no GT filter",
            "ascertainment": "Inherited from VCF and optional whitelist; no posterior/anchor selection, imputation, or sampling",
            "whitelist": dict(fingerprint(args.sites), sha256=sha256(args.sites)) if args.sites else None,
        })
        if args.validate_only:
            print(json.dumps(audit, indent=2, allow_nan=False))
            print("VALIDATION_PASS (full input scan; no files written)")
            return 0
        args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
        # Stage all artifacts. Audit is published LAST and is the completion marker.
        with tempfile.TemporaryDirectory(prefix=".pairwise-af-", dir=args.output_prefix.parent) as tmp:
            directory = Path(tmp)
            audit["plot"] = draw(histogram, args, directory)
            edges = np.linspace(0, 1, args.bins + 1)
            bins_path = directory / "AF_bins.tsv.gz"
            with gzip.open(bins_path, "wt") as handle:
                if args.bin_shape == "hex":
                    centres, sx, sy = hex_geometry(args.bins)
                    handle.write("hex_id\tleft_centre\tright_centre\twidth\theight\tcount\tproportion\n")
                    for index, (left, right) in enumerate(centres):
                        count = int(histogram[index])
                        fields = [index, left, right, sx, 2 * sy / 3,
                                  count, count / audit["records_retained"]]
                        handle.write("\t".join(map(str, fields)) + "\n")
                else:
                    handle.write("left_bin\tright_bin\tleft_min\tleft_max\tright_min\tright_max\tcount\tproportion\n")
                    for left in range(args.bins):
                        for right in range(args.bins):
                            count = int(histogram[left, right])
                            fields = [left, right, edges[left], edges[left + 1], edges[right], edges[right + 1],
                                      count, count / audit["records_retained"]]
                            handle.write("\t".join(map(str, fields)) + "\n")
            audit["outputs"] = {str(dest.resolve()): {"sha256": sha256(directory / suffix)}
                                for suffix, dest in zip(suffixes[:-1], destinations[:-1])}
            (directory / "AF.audit.json").write_text(json.dumps(audit, indent=2, allow_nan=False) + "\n")
            # Exclusive creation prevents accidental overwrite even with simultaneous runs.
            created = []
            try:
                import shutil
                for suffix, dest in zip(suffixes, destinations):
                    with dest.open("xb") as output:
                        created.append(dest)
                        with (directory / suffix).open("rb") as source:
                            shutil.copyfileobj(source, output)
            except BaseException:
                for dest in created:
                    dest.unlink(missing_ok=True)
                raise
        print(f"COMPLETE\tretained={audit['records_retained']:,}\tread={audit['records_read']:,}")
        for dest in destinations:
            print(dest.resolve())
        return 0
    except (ValueError, OSError, KeyError, ImportError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
