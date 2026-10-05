#!/usr/bin/env python3

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import logsumexp
from scipy.stats import gaussian_kde


PROJECT = Path("/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling")
ROOT = PROJECT / "analysis/genomewide_mixing"
INPUT = ROOT / "inputs"
FIT = ROOT / "observed_fit"
SIM = ROOT / "simulations"
JACK = ROOT / "block_jackknife"
FINAL = ROOT / "final"
FIG = ROOT / "figures"
MANIFEST = (
    PROJECT
    / "analysis/genomewide_mixing_input/metadata/chromosome_calling_manifest.tsv"
)

FOCAL = [
    "TN_Left_T2",
    "TN_Left_T3",
    "TN_Right_T2",
    "TN_Right_T3",
    "TN_Left_merged",
    "TN_Right_merged",
]
LEFT = "TN_Left_merged"
RIGHT = "TN_Right_merged"
FIT_PER_AUTOSOME = 5000
CALIBRATION_SITES = 30000
DEPTH_MIN = 5
DEPTH_MAX = 300
EFFECTIVE_DEPTH_CAP = 30
WINDOW_SIZE = 2_000_000
JACKKNIFE_BLOCK_SIZE = 5_000_000
MIXING_GRID = np.array([0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30])
GENOTYPE_LABELS = np.array(["0/0", "0/1", "1/1"])
G9_LEFT = np.repeat(np.array([0.0, 0.5, 1.0]), 3)
G9_RIGHT = np.tile(np.array([0.0, 0.5, 1.0]), 3)


def ensure_directories():
    for directory in [INPUT, FIT, SIM, JACK, FINAL, FIG]:
        directory.mkdir(parents=True, exist_ok=True)


def read_manifest():
    table = pd.read_csv(MANIFEST, sep="\t", dtype={"contig": str})

    if "molecule" not in table.columns:
        molecule_aliases = [
            column
            for column in (
                "assigned_molecule",
                "chromosome",
                "chromosome_number",
            )
            if column in table.columns
        ]

        if len(molecule_aliases) != 1:
            raise RuntimeError(
                "Could not resolve the chromosome-number column. "
                f"Available columns: {list(table.columns)}"
            )

        table = table.rename(
            columns={
                molecule_aliases[0]: "molecule",
            }
        )
    required = {"array_id", "molecule", "contig", "length", "chromosome_type"}
    missing = sorted(required.difference(table.columns))
    if missing:
        raise RuntimeError("Chromosome manifest lacks columns: " + ", ".join(missing))
    if len(table) != 21:
        raise RuntimeError(f"Expected 21 chromosome rows; found {len(table)}")
    return table.sort_values("array_id").reset_index(drop=True)


def input_path(row):
    return INPUT / f"chr{int(row.molecule)}_{row.contig}.counts.tsv.gz"


def parse_gt(gt):
    gt = gt.replace("|", "/")
    if gt in {"0/0", "0"}:
        return 0
    if gt in {"0/1", "1/0"}:
        return 1
    if gt in {"1/1", "1"}:
        return 2
    return -1


def parse_sample_field(value):
    pieces = value.split(":")
    gt = parse_gt(pieces[0]) if pieces else -1
    alt = -1
    dp = -1
    if len(pieces) > 1 and pieces[1] not in {"", "."}:
        ad = pieces[1].split(",")
        if len(ad) >= 2 and ad[1] not in {"", "."}:
            try:
                alt = int(ad[1])
            except ValueError:
                alt = -1
    if len(pieces) > 2 and pieces[2] not in {"", "."}:
        try:
            dp = int(pieces[2])
        except ValueError:
            dp = -1
    if dp >= 0 and alt > dp:
        alt = -1
    return gt, alt, dp


def read_header(path):
    with gzip.open(path, "rt") as handle:
        fields = handle.readline().rstrip("\n").split("\t")
    if fields[:4] != ["#CHROM", "POS", "REF", "ALT"]:
        raise RuntimeError(f"Unexpected header in {path}: {fields[:4]}")
    names = fields[4:]
    if len(names) != 24 or len(set(names)) != 24:
        raise RuntimeError(f"Expected 24 unique samples in {path}; found {len(names)}")
    absent = sorted(set(FOCAL).difference(names))
    if absent:
        raise RuntimeError(f"Missing focal samples in {path}: {', '.join(absent)}")
    return names


def effective_counts(alt, depth):
    n = np.minimum(np.asarray(depth, dtype=float), EFFECTIVE_DEPTH_CAP).astype(np.int16)
    af = np.asarray(alt, dtype=float) / np.asarray(depth, dtype=float)
    y = np.rint(af * n).astype(np.int16)
    return y, n


def weights_four(q):
    q = np.asarray(q, dtype=float)
    p = np.column_stack(((1.0 - q) ** 2, 2.0 * q * (1.0 - q), q**2))
    return (p[:, :, None] * p[:, None, :]).reshape(len(q), 9)


def weights_three(q):
    q = np.asarray(q, dtype=float)
    result = np.zeros((len(q), 9), dtype=float)
    for maternal in (0, 1):
        for sperm_left in (0, 1):
            for sperm_right in (0, 1):
                ones = maternal + sperm_left + sperm_right
                probability = q**ones * (1.0 - q) ** (3 - ones)
                left_dosage = maternal + sperm_left
                right_dosage = maternal + sperm_right
                index = left_dosage * 3 + right_dosage
                result[:, index] += probability
    return result


def sex_states_and_weights(q, female_side):
    q = np.asarray(q, dtype=float)
    diploid = np.column_stack(((1.0 - q) ** 2, 2.0 * q * (1.0 - q), q**2))
    haploid = np.column_stack((1.0 - q, q))
    if female_side == "left":
        g_left = np.repeat(np.array([0.0, 1.0]), 3)
        g_right = np.tile(np.array([0.0, 0.5, 1.0]), 2)
        weights = (haploid[:, :, None] * diploid[:, None, :]).reshape(len(q), 6)
    else:
        g_left = np.repeat(np.array([0.0, 0.5, 1.0]), 2)
        g_right = np.tile(np.array([0.0, 1.0]), 3)
        weights = (diploid[:, :, None] * haploid[:, None, :]).reshape(len(q), 6)
    return g_left, g_right, weights


def likelihood_batch(y_left, n_left, y_right, n_right, weights, g_left, g_right, c_left, c_right):
    error = 0.005
    p_left = (1.0 - c_left) * g_left + c_left * g_right
    p_right = (1.0 - c_right) * g_right + c_right * g_left
    p_left = np.clip(error + (1.0 - 2.0 * error) * p_left, 1e-7, 1.0 - 1e-7)
    p_right = np.clip(error + (1.0 - 2.0 * error) * p_right, 1e-7, 1.0 - 1e-7)
    log_prior = np.log(np.clip(weights, 1e-300, None))
    components = (
        log_prior
        + y_left[:, None] * np.log(p_left)[None, :]
        + (n_left - y_left)[:, None] * np.log1p(-p_left)[None, :]
        + y_right[:, None] * np.log(p_right)[None, :]
        + (n_right - y_right)[:, None] * np.log1p(-p_right)[None, :]
    )
    site_log_likelihood = logsumexp(components, axis=1)
    return site_log_likelihood, components


def total_log_likelihood(data, model, c_left, c_right):
    q = data["q"]
    weights = weights_three(q) if model == "three_founder" else weights_four(q)
    total = 0.0
    chunk = 20000
    for start in range(0, len(q), chunk):
        stop = min(start + chunk, len(q))
        site_ll, _ = likelihood_batch(
            data["y_left"][start:stop],
            data["n_left"][start:stop],
            data["y_right"][start:stop],
            data["n_right"][start:stop],
            weights[start:stop],
            G9_LEFT,
            G9_RIGHT,
            c_left,
            c_right,
        )
        total += float(site_ll.sum())
    return total


def fit_model(data, model):
    starts = [(0.02, 0.02), (0.10, 0.10), (0.18, 0.18), (0.25, 0.20)]
    best = None

    def objective(parameters):
        return -total_log_likelihood(data, model, float(parameters[0]), float(parameters[1]))

    for start in starts:
        result = minimize(
            objective,
            np.asarray(start),
            method="L-BFGS-B",
            bounds=((0.0, 0.40), (0.0, 0.40)),
            options={"maxiter": 100, "ftol": 1e-8},
        )
        if result.success and np.isfinite(result.fun):
            if best is None or result.fun < best.fun:
                best = result
    if best is None:
        raise RuntimeError(f"Every optimizer start failed for {model}")
    return {
        "model": model,
        "c_left": float(best.x[0]),
        "c_right": float(best.x[1]),
        "log_likelihood": float(-best.fun),
        "optimizer_message": str(best.message),
    }


def reservoir_for_chromosome(path, molecule, rng, maximum):
    names = read_header(path)
    left_index = names.index(LEFT)
    right_index = names.index(RIGHT)
    control_indices = [i for i, name in enumerate(names) if name not in FOCAL]
    if len(control_indices) != 18:
        raise RuntimeError(f"Expected 18 controls in {path}; found {len(control_indices)}")
    reservoir = []
    eligible = 0
    heterozygote_af = []
    with gzip.open(path, "rt") as handle:
        next(handle)
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 28:
                raise RuntimeError(f"Malformed row with {len(fields)} fields in {path}")
            sample_fields = fields[4:]
            left_gt, left_alt, left_dp = parse_sample_field(sample_fields[left_index])
            right_gt, right_alt, right_dp = parse_sample_field(sample_fields[right_index])
            if not (
                DEPTH_MIN <= left_dp <= DEPTH_MAX
                and DEPTH_MIN <= right_dp <= DEPTH_MAX
                and 0 <= left_alt <= left_dp
                and 0 <= right_alt <= right_dp
            ):
                continue
            control_gt = np.full(18, -1, dtype=np.int8)
            control_alt = np.full(18, -1, dtype=np.int16)
            control_dp = np.full(18, -1, dtype=np.int16)
            allele_sum = 0
            called = 0
            for out_index, sample_index in enumerate(control_indices):
                gt, alt, dp = parse_sample_field(sample_fields[sample_index])
                control_gt[out_index] = gt
                control_alt[out_index] = alt
                control_dp[out_index] = dp
                if gt >= 0:
                    allele_sum += gt
                    called += 1
                if gt == 1 and dp >= 8 and alt >= 0 and len(heterozygote_af) < 100000:
                    heterozygote_af.append(alt / dp)
            if called < 13:
                continue
            q = (allele_sum + 1.0) / (2.0 * called + 2.0)
            if not (0.01 <= q <= 0.99):
                continue
            record = (
                int(molecule), int(fields[1]), q,
                left_alt, left_dp, right_alt, right_dp,
                control_gt, control_alt, control_dp,
            )
            eligible += 1
            if len(reservoir) < maximum:
                reservoir.append(record)
            else:
                replacement = int(rng.integers(0, eligible))
                if replacement < maximum:
                    reservoir[replacement] = record
    return names, reservoir, eligible, heterozygote_af


def records_to_data(records):
    if not records:
        raise RuntimeError("No eligible records were available")
    chrom = np.array([row[0] for row in records], dtype=np.int16)
    pos = np.array([row[1] for row in records], dtype=np.int32)
    q = np.array([row[2] for row in records], dtype=np.float64)
    left_alt = np.array([row[3] for row in records], dtype=np.int16)
    left_dp = np.array([row[4] for row in records], dtype=np.int16)
    right_alt = np.array([row[5] for row in records], dtype=np.int16)
    right_dp = np.array([row[6] for row in records], dtype=np.int16)
    controls = np.stack([row[7] for row in records]).astype(np.int8)
    controls_alt = np.stack([row[8] for row in records]).astype(np.int16)
    controls_dp = np.stack([row[9] for row in records]).astype(np.int16)
    y_left, n_left = effective_counts(left_alt, left_dp)
    y_right, n_right = effective_counts(right_alt, right_dp)
    return {
        "chrom": chrom,
        "pos": pos,
        "q": q,
        "y_left": y_left,
        "n_left": n_left,
        "y_right": y_right,
        "n_right": n_right,
        "left_alt": left_alt,
        "left_dp": left_dp,
        "right_alt": right_alt,
        "right_dp": right_dp,
        "controls": controls,
        "controls_alt": controls_alt,
        "controls_dp": controls_dp,
    }


def subset_data(data, indices):
    return {key: value[indices] for key, value in data.items()}


def estimate_phi(heterozygote_af):
    values = np.asarray(heterozygote_af, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) < 100:
        return 80.0
    variance = float(np.var(values))
    rho = max(1e-4, min(0.05, (variance - 0.0125) / 0.2375))
    return float(max(20.0, min(500.0, 1.0 / rho - 1.0)))


def infer_female_side(manifest, autosomal_data):
    z_row = manifest.loc[manifest["molecule"].astype(int) == 21]
    if len(z_row) != 1:
        raise RuntimeError("Expected exactly one chromosome-21/Z row")
    z_path = input_path(z_row.iloc[0])
    names = read_header(z_path)
    left_index = names.index(LEFT)
    right_index = names.index(RIGHT)
    z_left = []
    z_right = []
    with gzip.open(z_path, "rt") as handle:
        next(handle)
        for number, line in enumerate(handle):
            if number % 25 != 0:
                continue
            fields = line.rstrip("\n").split("\t")
            _, _, left_dp = parse_sample_field(fields[4 + left_index])
            _, _, right_dp = parse_sample_field(fields[4 + right_index])
            if left_dp > 0 and right_dp > 0:
                z_left.append(left_dp)
                z_right.append(right_dp)
            if len(z_left) >= 25000:
                break
    auto_left = float(np.median(autosomal_data["left_dp"]))
    auto_right = float(np.median(autosomal_data["right_dp"]))
    ratio_left = float(np.median(z_left)) / auto_left
    ratio_right = float(np.median(z_right)) / auto_right
    female_side = "left" if ratio_left < ratio_right else "right"
    return female_side, ratio_left, ratio_right


def maps_for_batch(y_left, n_left, y_right, n_right, q, fit, model):
    weights = weights_three(q) if model == "three_founder" else weights_four(q)
    site_ll, components = likelihood_batch(
        y_left, n_left, y_right, n_right, weights,
        G9_LEFT, G9_RIGHT, fit["c_left"], fit["c_right"],
    )
    posterior = np.exp(components - site_ll[:, None])
    map_index = np.argmax(posterior, axis=1)
    return (
        site_ll,
        G9_LEFT[map_index],
        G9_RIGHT[map_index],
        posterior[np.arange(len(map_index)), map_index],
    )


def format_genotype(value):
    if value == 0.0:
        return "0/0"
    if value == 0.5:
        return "0/1"
    return "1/1"


def analyse_all_sites(manifest, fit_three, fit_four, female_side):
    window_rows = defaultdict(lambda: {
        "sites": 0, "delta": 0.0, "sex_delta": 0.0,
        "sex_sites": 0, "left_dp": 0.0, "right_dp": 0.0,
    })
    for row in manifest.itertuples(index=False):
        path = input_path(row)
        names = read_header(path)
        left_index = names.index(LEFT)
        right_index = names.index(RIGHT)
        control_indices = [i for i, name in enumerate(names) if name not in FOCAL]
        output = FIT / f"chr{int(row.molecule)}_{row.contig}.adjusted_calls.tsv.gz"
        with gzip.open(path, "rt") as source, gzip.open(output.with_suffix(output.suffix + ".new"), "wt") as destination:
            next(source)
            destination.write(
                "chrom\tpos\tref\talt\tleft_alt\tleft_dp\tright_alt\tright_dp\t"
                "left_af\tright_af\tmap3_left\tmap3_right\tmap3_posterior\t"
                "map4_left\tmap4_right\tmap4_posterior\tdelta_logL_four_minus_three\t"
                "delta_logL_sex_aware_minus_four\n"
            )
            batch = []

            def flush(records):
                if not records:
                    return
                positions = np.array([x[0] for x in records], dtype=np.int32)
                q = np.array([x[3] for x in records], dtype=float)
                left_alt = np.array([x[4] for x in records], dtype=np.int16)
                left_dp = np.array([x[5] for x in records], dtype=np.int16)
                right_alt = np.array([x[6] for x in records], dtype=np.int16)
                right_dp = np.array([x[7] for x in records], dtype=np.int16)
                y_left, n_left = effective_counts(left_alt, left_dp)
                y_right, n_right = effective_counts(right_alt, right_dp)
                ll3, ml3, mr3, pp3 = maps_for_batch(
                    y_left, n_left, y_right, n_right, q, fit_three, "three_founder"
                )
                ll4, ml4, mr4, pp4 = maps_for_batch(
                    y_left, n_left, y_right, n_right, q, fit_four, "four_founder"
                )
                delta = ll4 - ll3
                if int(row.molecule) == 21:
                    sex_left, sex_right, sex_weights = sex_states_and_weights(
                        q, female_side
                    )
                    llsex, _ = likelihood_batch(
                        y_left, n_left, y_right, n_right, sex_weights,
                        sex_left, sex_right,
                        fit_four["c_left"], fit_four["c_right"],
                    )
                    sex_delta = llsex - ll4
                else:
                    sex_delta = np.full(len(records), np.nan)
                for i, record in enumerate(records):
                    position, ref, alt = record[:3]
                    destination.write(
                        f"{row.contig}\t{position}\t{ref}\t{alt}\t"
                        f"{left_alt[i]}\t{left_dp[i]}\t{right_alt[i]}\t{right_dp[i]}\t"
                        f"{left_alt[i] / left_dp[i]:.6f}\t{right_alt[i] / right_dp[i]:.6f}\t"
                        f"{format_genotype(ml3[i])}\t{format_genotype(mr3[i])}\t{pp3[i]:.6g}\t"
                        f"{format_genotype(ml4[i])}\t{format_genotype(mr4[i])}\t{pp4[i]:.6g}\t"
                        f"{delta[i]:.8g}\t{sex_delta[i]:.8g}\n"
                    )
                    window_start = ((int(position) - 1) // WINDOW_SIZE) * WINDOW_SIZE
                    key = (int(row.molecule), str(row.contig), window_start, str(row.chromosome_type))
                    target = window_rows[key]
                    target["sites"] += 1
                    target["delta"] += float(delta[i])
                    if np.isfinite(sex_delta[i]):
                        target["sex_delta"] += float(sex_delta[i])
                        target["sex_sites"] += 1
                    target["left_dp"] += float(left_dp[i])
                    target["right_dp"] += float(right_dp[i])

            for line in source:
                fields = line.rstrip("\n").split("\t")
                sample_fields = fields[4:]
                _, left_alt, left_dp = parse_sample_field(sample_fields[left_index])
                _, right_alt, right_dp = parse_sample_field(sample_fields[right_index])
                if not (
                    DEPTH_MIN <= left_dp <= DEPTH_MAX
                    and DEPTH_MIN <= right_dp <= DEPTH_MAX
                    and 0 <= left_alt <= left_dp
                    and 0 <= right_alt <= right_dp
                ):
                    continue
                allele_sum = 0
                called = 0
                for sample_index in control_indices:
                    gt, _, _ = parse_sample_field(sample_fields[sample_index])
                    if gt >= 0:
                        allele_sum += gt
                        called += 1
                if called < 13:
                    continue
                q = (allele_sum + 1.0) / (2.0 * called + 2.0)
                if not (0.01 <= q <= 0.99):
                    continue
                batch.append((int(fields[1]), fields[2], fields[3], q, left_alt, left_dp, right_alt, right_dp))
                if len(batch) >= 5000:
                    flush(batch)
                    batch = []
            flush(batch)
        os.replace(output.with_suffix(output.suffix + ".new"), output)
        print(f"ADJUSTED_CALLS_COMPLETE chromosome={row.molecule} file={output}", flush=True)

    rows = []
    for key, values in sorted(window_rows.items()):
        molecule, contig, start, chromosome_type = key
        n = values["sites"]
        rows.append({
            "molecule": molecule,
            "contig": contig,
            "window_start": start,
            "window_end": start + WINDOW_SIZE,
            "chromosome_type": chromosome_type,
            "sites": n,
            "delta_logL_four_minus_three": values["delta"],
            "delta_logL_four_minus_three_per_site": values["delta"] / n,
            "delta_logL_sex_aware_minus_four": (
                values["sex_delta"] if values["sex_sites"] else np.nan
            ),
            "delta_logL_sex_aware_minus_four_per_site": (
                values["sex_delta"] / values["sex_sites"]
                if values["sex_sites"] else np.nan
            ),
            "mean_left_dp": values["left_dp"] / n,
            "mean_right_dp": values["right_dp"] / n,
        })
    pd.DataFrame(rows).to_csv(FIT / "genomewide_2Mb_window_evidence.tsv", sep="\t", index=False)


def command_fit(_arguments):
    ensure_directories()
    manifest = read_manifest()
    rng = np.random.default_rng(20260911)
    print("STAGE=observed_fit STATUS=START", flush=True)
    records = []
    audit_rows = []
    heterozygote_af = []
    canonical_names = None
    for row in manifest.itertuples(index=False):
        path = input_path(row)
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing extracted count table: {path}")
        names, selected, eligible, het_af = reservoir_for_chromosome(
            path, int(row.molecule), rng,
            FIT_PER_AUTOSOME if int(row.molecule) <= 20 else 2000,
        )
        print(
            f"INPUT_SCAN chromosome={int(row.molecule)} "
            f"eligible={eligible} retained={len(selected)}",
            flush=True,
        )
        if canonical_names is None:
            canonical_names = names
        elif names != canonical_names:
            raise RuntimeError(f"Sample order differs in {path}")
        audit_rows.append({
            "molecule": int(row.molecule), "contig": row.contig,
            "eligible_sites": eligible, "reservoir_sites": len(selected),
        })
        if int(row.molecule) <= 20:
            records.extend(selected)
            heterozygote_af.extend(het_af)
    data = records_to_data(records)
    order = rng.permutation(len(data["q"]))
    data = subset_data(data, order)
    fit_three = fit_model(data, "three_founder")
    fit_four = fit_model(data, "four_founder")
    calibration_indices = np.arange(min(CALIBRATION_SITES, len(data["q"])))
    calibration_data = subset_data(data, calibration_indices)
    calibration_three = fit_model(calibration_data, "three_founder")
    calibration_four = fit_model(calibration_data, "four_founder")
    fit_rows = []
    for scope, result_set in [
        ("observed_fit", [fit_three, fit_four]),
        ("calibration_subset", [calibration_three, calibration_four]),
    ]:
        for result in result_set:
            fit_rows.append({"scope": scope, **result})
    pd.DataFrame(fit_rows).to_csv(FIT / "observed_founder_model_fits.tsv", sep="\t", index=False)
    phi = estimate_phi(heterozygote_af)
    female_side, z_left_ratio, z_right_ratio = infer_female_side(manifest, data)
    np.savez_compressed(
        FIT / "autosomal_fit_input.npz",
        **data,
        phi=np.array([phi]),
        sample_names=np.array(canonical_names),
    )
    np.savez_compressed(
        FIT / "calibration_fit_input.npz",
        **calibration_data,
        phi=np.array([phi]),
        sample_names=np.array(canonical_names),
    )
    block_rows = []
    block_id = 0
    for molecule in range(1, 21):
        positions = data["pos"][data["chrom"] == molecule]
        if len(positions) == 0:
            continue
        for start in sorted(set(((positions - 1) // JACKKNIFE_BLOCK_SIZE) * JACKKNIFE_BLOCK_SIZE)):
            block_rows.append({
                "block_id": block_id, "molecule": molecule,
                "block_start": int(start), "block_end": int(start + JACKKNIFE_BLOCK_SIZE),
                "fit_sites": int(np.sum((data["chrom"] == molecule) & (data["pos"] > start) & (data["pos"] <= start + JACKKNIFE_BLOCK_SIZE))),
            })
            block_id += 1
    pd.DataFrame(block_rows).to_csv(FIT / "autosomal_5Mb_block_manifest.tsv", sep="\t", index=False)
    pd.DataFrame(audit_rows).to_csv(FIT / "input_site_audit.tsv", sep="\t", index=False)
    metadata = {
        "autosomal_fit_sites": int(len(data["q"])),
        "calibration_sites": int(len(calibration_data["q"])),
        "controls": 18,
        "focal_samples": FOCAL,
        "primary_focal_pair": [LEFT, RIGHT],
        "effective_depth_cap": EFFECTIVE_DEPTH_CAP,
        "estimated_beta_binomial_phi_for_simulation": phi,
        "inferred_lower_Z_depth_side": female_side,
        "left_Z_to_autosome_variant_site_depth_ratio": z_left_ratio,
        "right_Z_to_autosome_variant_site_depth_ratio": z_right_ratio,
        "founder_inference_chromosomes": "1-20",
        "Z_usage": "held-out sex-ploidy diagnostic",
    }
    with open(FIT / "model_metadata.json", "w") as handle:
        json.dump(metadata, handle, indent=2)
    analyse_all_sites(manifest, fit_three, fit_four, female_side)
    (FIT / "observed_fit.complete").write_text("complete\n")
    print("OBSERVED GENOME-WIDE FIT COMPLETE")
    print(pd.DataFrame(fit_rows).to_string(index=False))
    print(json.dumps(metadata, indent=2))


def load_npz(path):
    source = np.load(path, allow_pickle=False)
    return {key: source[key] for key in source.files if key not in {"phi", "sample_names"}}, float(source["phi"][0])


def impute_genotypes(genotypes, q, rng):
    result = genotypes.copy()
    missing = result < 0
    if np.any(missing):
        result[missing] = rng.binomial(2, q[missing])
    return result


def sampled_allele(genotype, rng):
    result = np.zeros(len(genotype), dtype=np.int8)
    result[genotype == 2] = 1
    heterozygous = genotype == 1
    result[heterozygous] = rng.binomial(1, 0.5, int(heterozygous.sum()))
    return result


def simulate_counts(g_left, g_right, n_left, n_right, c_left, c_right, phi, rng):
    error = 0.005
    p_left = error + (1.0 - 2.0 * error) * ((1.0 - c_left) * g_left + c_left * g_right)
    p_right = error + (1.0 - 2.0 * error) * ((1.0 - c_right) * g_right + c_right * g_left)
    p_left = np.clip(p_left, 1e-5, 1.0 - 1e-5)
    p_right = np.clip(p_right, 1e-5, 1.0 - 1e-5)
    draw_left = rng.beta(p_left * phi, (1.0 - p_left) * phi)
    draw_right = rng.beta(p_right * phi, (1.0 - p_right) * phi)
    return rng.binomial(n_left, draw_left), rng.binomial(n_right, draw_right)


def command_simulate(arguments):
    ensure_directories()
    data, phi = load_npz(FIT / "calibration_fit_input.npz")
    task = int(arguments.task)
    if not 0 <= task < 98:
        raise RuntimeError(f"Simulation task must be 0-97; received {task}")
    true_model = "three_founder" if task < 49 else "four_founder"
    grid_index = task % 49
    true_c_left = float(MIXING_GRID[grid_index // 7])
    true_c_right = float(MIXING_GRID[grid_index % 7])
    print(
        f"STAGE=simulation STATUS=START task={task} model={true_model} "
        f"c_left={true_c_left} c_right={true_c_right}",
        flush=True,
    )
    controls = data["controls"]
    q = data["q"]
    rows = []
    for replicate in range(4):
        rng = np.random.default_rng(2026091100 + task * 10 + replicate)
        donors = rng.choice(controls.shape[1], size=3, replace=False)
        donor_a = impute_genotypes(controls[:, donors[0]], q, rng)
        donor_b = impute_genotypes(controls[:, donors[1]], q, rng)
        donor_c = impute_genotypes(controls[:, donors[2]], q, rng)
        if true_model == "three_founder":
            shared = sampled_allele(donor_a, rng)
            g_left = (shared + sampled_allele(donor_b, rng)) / 2.0
            g_right = (shared + sampled_allele(donor_c, rng)) / 2.0
        else:
            g_left = donor_a / 2.0
            g_right = donor_b / 2.0
        y_left, y_right = simulate_counts(
            g_left, g_right, data["n_left"], data["n_right"],
            true_c_left, true_c_right, phi, rng,
        )
        simulated = dict(data)
        simulated["y_left"] = y_left.astype(np.int16)
        simulated["y_right"] = y_right.astype(np.int16)
        fit_three = fit_model(simulated, "three_founder")
        fit_four = fit_model(simulated, "four_founder")
        rows.append({
            "task": task, "replicate": replicate, "true_model": true_model,
            "true_c_left": true_c_left, "true_c_right": true_c_right,
            "donor_a": int(donors[0]), "donor_b": int(donors[1]), "donor_c": int(donors[2]),
            "genotype_contrast": float(np.mean(np.abs(g_left - g_right))),
            "informative_fraction": float(np.mean(g_left != g_right)),
            "estimated_c_left_three": fit_three["c_left"],
            "estimated_c_right_three": fit_three["c_right"],
            "estimated_c_left_four": fit_four["c_left"],
            "estimated_c_right_four": fit_four["c_right"],
            "log_likelihood_three": fit_three["log_likelihood"],
            "log_likelihood_four": fit_four["log_likelihood"],
            "delta_logL_four_minus_three": fit_four["log_likelihood"] - fit_three["log_likelihood"],
        })
        print(
            f"SIMULATION_PROGRESS task={task} replicate={replicate + 1}/4",
            flush=True,
        )
    output = SIM / f"simulation_task_{task:03d}.tsv"
    pd.DataFrame(rows).to_csv(output.with_suffix(".tsv.new"), sep="\t", index=False)
    os.replace(output.with_suffix(".tsv.new"), output)
    print(f"SIMULATION_COMPLETE task={task} model={true_model} c_left={true_c_left} c_right={true_c_right}")


def command_jackknife(arguments):
    ensure_directories()
    blocks = pd.read_csv(FIT / "autosomal_5Mb_block_manifest.tsv", sep="\t")
    task = int(arguments.task)
    if task >= len(blocks):
        output = JACK / f"block_{task:03d}.skip.tsv"
        output.write_text(f"task\treason\n{task}\toutside_block_manifest\n")
        print(f"JACKKNIFE_SKIP task={task} blocks={len(blocks)}")
        return
    block = blocks.iloc[task]
    print(
        f"STAGE=jackknife STATUS=START block={task} "
        f"chromosome={int(block.molecule)}",
        flush=True,
    )
    data, _phi = load_npz(FIT / "autosomal_fit_input.npz")
    keep = ~(
        (data["chrom"] == int(block.molecule))
        & (data["pos"] > int(block.block_start))
        & (data["pos"] <= int(block.block_end))
    )
    reduced = subset_data(data, keep)
    fit_three = fit_model(reduced, "three_founder")
    fit_four = fit_model(reduced, "four_founder")
    row = {
        "block_id": int(block.block_id), "molecule": int(block.molecule),
        "block_start": int(block.block_start), "block_end": int(block.block_end),
        "omitted_fit_sites": int((~keep).sum()), "retained_fit_sites": int(keep.sum()),
        "c_left_three": fit_three["c_left"], "c_right_three": fit_three["c_right"],
        "c_left_four": fit_four["c_left"], "c_right_four": fit_four["c_right"],
        "log_likelihood_three": fit_three["log_likelihood"],
        "log_likelihood_four": fit_four["log_likelihood"],
        "delta_logL_four_minus_three": fit_four["log_likelihood"] - fit_three["log_likelihood"],
    }
    output = JACK / f"block_{task:03d}.tsv"
    pd.DataFrame([row]).to_csv(output.with_suffix(".tsv.new"), sep="\t", index=False)
    os.replace(output.with_suffix(".tsv.new"), output)
    print(f"JACKKNIFE_COMPLETE block={task} chromosome={int(block.molecule)}")


def monotone_calibrator(simulations, model, side):
    estimate_column = f"estimated_c_{side}_{'three' if model == 'three_founder' else 'four'}"
    subset = simulations[simulations.true_model == model]
    trend = (
        subset.groupby(f"true_c_{side}", as_index=False)[estimate_column]
        .median()
        .sort_values(f"true_c_{side}")
    )
    x = np.maximum.accumulate(trend[estimate_column].to_numpy(float))
    x = x + np.arange(len(x)) * 1e-10
    y = trend[f"true_c_{side}"].to_numpy(float)

    def calibrate(value):
        return float(np.interp(value, x, y, left=y[0], right=y[-1]))

    return calibrate, trend


def jackknife_interval(values, full_value):
    values = np.asarray(values, dtype=float)
    count = len(values)
    mean = float(np.mean(values))
    standard_error = math.sqrt((count - 1.0) / count * float(np.sum((values - mean) ** 2)))
    return standard_error, max(0.0, full_value - 1.96 * standard_error), min(0.40, full_value + 1.96 * standard_error)


def command_collect(_arguments):
    ensure_directories()
    simulation_files = sorted(SIM.glob("simulation_task_*.tsv"))
    if len(simulation_files) != 98:
        raise RuntimeError(f"Expected 98 simulation files; found {len(simulation_files)}")
    simulations = pd.concat([pd.read_csv(path, sep="\t") for path in simulation_files], ignore_index=True)
    if len(simulations) != 392:
        raise RuntimeError(f"Expected 392 simulation rows; found {len(simulations)}")
    blocks = pd.read_csv(FIT / "autosomal_5Mb_block_manifest.tsv", sep="\t")
    jackknife_files = [JACK / f"block_{int(block_id):03d}.tsv" for block_id in blocks.block_id]
    missing = [str(path) for path in jackknife_files if not path.is_file()]
    if missing:
        raise RuntimeError("Missing jackknife outputs:\n" + "\n".join(missing))
    jackknife = pd.concat([pd.read_csv(path, sep="\t") for path in jackknife_files], ignore_index=True)
    fits = pd.read_csv(FIT / "observed_founder_model_fits.tsv", sep="\t")
    observed = fits[fits.scope == "observed_fit"].set_index("model")
    observed_cal = fits[fits.scope == "calibration_subset"].set_index("model")
    observed_delta_cal = (
        observed_cal.loc["four_founder", "log_likelihood"]
        - observed_cal.loc["three_founder", "log_likelihood"]
    )
    observed_delta_full = (
        observed.loc["four_founder", "log_likelihood"]
        - observed.loc["three_founder", "log_likelihood"]
    )
    jackknife["removed_block_delta_logL_influence"] = (
        observed_delta_full - jackknife["delta_logL_four_minus_three"]
    )
    jackknife["removed_block_delta_logL_influence_per_site"] = (
        jackknife["removed_block_delta_logL_influence"]
        / jackknife["omitted_fit_sites"]
    )
    three_delta = simulations.loc[
        simulations.true_model == "three_founder",
        "delta_logL_four_minus_three",
    ].to_numpy(float)
    four_delta = simulations.loc[
        simulations.true_model == "four_founder",
        "delta_logL_four_minus_three",
    ].to_numpy(float)

    def log_density(values, point):
        try:
            density = float(gaussian_kde(values)(point)[0])
            if np.isfinite(density) and density > 0.0:
                return float(np.log(density))
        except (ValueError, np.linalg.LinAlgError):
            pass
        scale = max(float(np.std(values, ddof=1)), 1e-6)
        center = float(np.mean(values))
        return float(-0.5 * ((point - center) / scale) ** 2 - np.log(scale))

    log_three = log_density(three_delta, observed_delta_cal)
    log_four = log_density(four_delta, observed_delta_cal)
    posterior_four = float(1.0 / (1.0 + np.exp(np.clip(log_three - log_four, -700, 700))))
    selected_model = "four_founder" if posterior_four >= 0.5 else "three_founder"
    summary_rows = []
    calibration_rows = []
    for model in ["three_founder", "four_founder"]:
        short = "three" if model == "three_founder" else "four"
        for side in ["left", "right"]:
            calibrate, trend = monotone_calibrator(simulations, model, side)
            trend = trend.copy()
            trend["model"] = model
            trend["side"] = side
            calibration_rows.append(trend)
            raw = float(observed.loc[model, f"c_{side}"])
            calibrated = calibrate(raw)
            jk_raw = jackknife[f"c_{side}_{short}"].to_numpy(float)
            jk_calibrated = np.array([calibrate(value) for value in jk_raw])
            se, lower, upper = jackknife_interval(jk_calibrated, calibrated)
            summary_rows.append({
                "model": model, "side": side, "raw_mixing": raw,
                "calibrated_mixing": calibrated, "jackknife_se": se,
                "ci95_lower": lower, "ci95_upper": upper,
            })
    summary = pd.DataFrame(summary_rows)
    summary["selected_model"] = selected_model
    summary["calibrated_posterior_four_founder"] = posterior_four
    summary["observed_calibration_delta_logL_four_minus_three"] = observed_delta_cal
    summary.to_csv(FINAL / "calibrated_founder_and_mixing_summary.tsv", sep="\t", index=False)
    pd.concat(calibration_rows, ignore_index=True).to_csv(
        FINAL / "mixing_calibration_curves.tsv", sep="\t", index=False
    )
    simulations.to_csv(FINAL / "all_simulation_results.tsv.gz", sep="\t", index=False, compression="gzip")
    jackknife.to_csv(FINAL / "all_5Mb_block_jackknife_results.tsv", sep="\t", index=False)
    jackknife[[
        "block_id", "molecule", "block_start", "block_end",
        "omitted_fit_sites", "removed_block_delta_logL_influence",
        "removed_block_delta_logL_influence_per_site",
    ]].to_csv(
        FINAL / "jackknife_block_evidence_manhattan.tsv", sep="\t", index=False
    )
    selection = {
        "selected_founder_model": selected_model,
        "calibrated_posterior_four_founder": posterior_four,
        "observed_calibration_delta_logL_four_minus_three": float(observed_delta_cal),
        "simulation_replicates": int(len(simulations)),
        "jackknife_blocks": int(len(jackknife)),
        "selected_left_calibrated_mixing": float(summary.loc[(summary.model == selected_model) & (summary.side == "left"), "calibrated_mixing"].iloc[0]),
        "selected_right_calibrated_mixing": float(summary.loc[(summary.model == selected_model) & (summary.side == "right"), "calibrated_mixing"].iloc[0]),
    }
    with open(FINAL / "final_model_selection.json", "w") as handle:
        json.dump(selection, handle, indent=2)
    checksums = []
    for path in sorted(FINAL.glob("*")):
        if path.is_file() and path.name != "output_checksums.sha256":
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            checksums.append(f"{digest}  {path.name}")
    (FINAL / "output_checksums.sha256").write_text("\n".join(checksums) + "\n")
    (FINAL / "collection.complete").write_text("complete\n")
    print("CALIBRATION AND JACKKNIFE COLLECTION COMPLETE")
    print(summary.to_string(index=False))
    print(json.dumps(selection, indent=2))


def command_adjust(arguments):
    ensure_directories()
    manifest = read_manifest()
    task = int(arguments.task)
    matches = manifest[manifest.array_id.astype(int) == task]
    if len(matches) != 1:
        raise RuntimeError(
            f"Expected one chromosome row for task {task}; found {len(matches)}"
        )
    row = next(matches.itertuples(index=False))
    path = input_path(row)
    names = read_header(path)
    left_index = names.index(LEFT)
    right_index = names.index(RIGHT)
    control_indices = [i for i, name in enumerate(names) if name not in FOCAL]
    summary = pd.read_csv(
        FINAL / "calibrated_founder_and_mixing_summary.tsv", sep="\t"
    )
    selection = json.loads((FINAL / "final_model_selection.json").read_text())
    metadata = json.loads((FIT / "model_metadata.json").read_text())
    selected_model = selection["selected_founder_model"]
    selected_rows = summary[summary.model == selected_model].set_index("side")
    c_left = float(selected_rows.loc["left", "calibrated_mixing"])
    c_right = float(selected_rows.loc["right", "calibrated_mixing"])
    output_dir = FINAL / "adjusted_calls"
    output_dir.mkdir(parents=True, exist_ok=True)
    output = (
        output_dir
        / f"chr{int(row.molecule)}_{row.contig}.mixing_adjusted_true_states.tsv.gz"
    )
    temporary = output.with_suffix(output.suffix + ".new")

    with gzip.open(path, "rt") as source, gzip.open(
        temporary, "wt"
    ) as destination:
        next(source)
        destination.write(
            "chrom\tpos\tref\talt\tleft_alt\tleft_dp\tright_alt\tright_dp\t"
            "left_af\tright_af\tground_left_GT\tground_right_GT\t"
            "pair_posterior\tmodel_used\tcalibrated_c_left\tcalibrated_c_right\n"
        )
        batch = []

        def flush(records):
            if not records:
                return
            q_values = np.array([record[3] for record in records], dtype=float)
            left_alt_values = np.array(
                [record[4] for record in records], dtype=np.int16
            )
            left_dp_values = np.array(
                [record[5] for record in records], dtype=np.int16
            )
            right_alt_values = np.array(
                [record[6] for record in records], dtype=np.int16
            )
            right_dp_values = np.array(
                [record[7] for record in records], dtype=np.int16
            )
            y_left, n_left = effective_counts(left_alt_values, left_dp_values)
            y_right, n_right = effective_counts(right_alt_values, right_dp_values)

            if int(row.molecule) == 21:
                female_side = metadata["inferred_lower_Z_depth_side"]
                g_left, g_right, weights = sex_states_and_weights(
                    q_values, female_side
                )
                site_ll, components = likelihood_batch(
                    y_left,
                    n_left,
                    y_right,
                    n_right,
                    weights,
                    g_left,
                    g_right,
                    c_left,
                    c_right,
                )
                model_used = f"sex_aware_Z_{female_side}_hemizygous"
            else:
                weights = (
                    weights_three(q_values)
                    if selected_model == "three_founder"
                    else weights_four(q_values)
                )
                g_left = G9_LEFT
                g_right = G9_RIGHT
                female_side = "none"
                site_ll, components = likelihood_batch(
                    y_left,
                    n_left,
                    y_right,
                    n_right,
                    weights,
                    g_left,
                    g_right,
                    c_left,
                    c_right,
                )
                model_used = selected_model

            posterior = np.exp(components - site_ll[:, None])
            map_index = np.argmax(posterior, axis=1)
            map_probability = posterior[np.arange(len(map_index)), map_index]

            for index, record in enumerate(records):
                position, ref, alt = record[:3]
                left_value = float(g_left[map_index[index]])
                right_value = float(g_right[map_index[index]])

                if int(row.molecule) == 21 and female_side == "left":
                    left_label = "0" if left_value == 0.0 else "1"
                else:
                    left_label = format_genotype(left_value)

                if int(row.molecule) == 21 and female_side == "right":
                    right_label = "0" if right_value == 0.0 else "1"
                else:
                    right_label = format_genotype(right_value)

                destination.write(
                    f"{row.contig}\t{position}\t{ref}\t{alt}\t"
                    f"{left_alt_values[index]}\t{left_dp_values[index]}\t"
                    f"{right_alt_values[index]}\t{right_dp_values[index]}\t"
                    f"{left_alt_values[index] / left_dp_values[index]:.6f}\t"
                    f"{right_alt_values[index] / right_dp_values[index]:.6f}\t"
                    f"{left_label}\t{right_label}\t{map_probability[index]:.8g}\t"
                    f"{model_used}\t{c_left:.8g}\t{c_right:.8g}\n"
                )

        for line in source:
            fields = line.rstrip("\n").split("\t")
            sample_fields = fields[4:]
            _, left_alt, left_dp = parse_sample_field(sample_fields[left_index])
            _, right_alt, right_dp = parse_sample_field(sample_fields[right_index])
            if not (
                DEPTH_MIN <= left_dp <= DEPTH_MAX
                and DEPTH_MIN <= right_dp <= DEPTH_MAX
                and 0 <= left_alt <= left_dp
                and 0 <= right_alt <= right_dp
            ):
                continue
            allele_sum = 0
            called = 0
            for sample_index in control_indices:
                genotype, _, _ = parse_sample_field(sample_fields[sample_index])
                if genotype >= 0:
                    allele_sum += genotype
                    called += 1
            if called < 13:
                continue
            q_value = (allele_sum + 1.0) / (2.0 * called + 2.0)
            if not (0.01 <= q_value <= 0.99):
                continue
            batch.append(
                (
                    int(fields[1]),
                    fields[2],
                    fields[3],
                    q_value,
                    left_alt,
                    left_dp,
                    right_alt,
                    right_dp,
                )
            )
            if len(batch) >= 5000:
                flush(batch)
                batch = []
        flush(batch)

    os.replace(temporary, output)
    print(
        f"FINAL_ADJUSTED_CALLS_COMPLETE chromosome={int(row.molecule)} "
        f"model={selected_model} c_left={c_left:.5f} "
        f"c_right={c_right:.5f} file={output}"
    )


def command_plot(_arguments):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    ensure_directories()
    simulations = pd.read_csv(FINAL / "all_simulation_results.tsv.gz", sep="\t")
    summary = pd.read_csv(FINAL / "calibrated_founder_and_mixing_summary.tsv", sep="\t")
    jackknife = pd.read_csv(FINAL / "all_5Mb_block_jackknife_results.tsv", sep="\t")
    windows = pd.read_csv(FIT / "genomewide_2Mb_window_evidence.tsv", sep="\t")
    fits = pd.read_csv(FIT / "observed_founder_model_fits.tsv", sep="\t")
    metadata = json.loads((FIT / "model_metadata.json").read_text())
    selection = json.loads((FINAL / "final_model_selection.json").read_text())
    blue = "#0072B2"
    orange = "#D55E00"
    green = "#009E73"
    purple = "#CC79A7"
    grey = "#666666"
    pdf_path = FIG / "genomewide_three_vs_four_founder_mixing_narrative.pdf"
    with PdfPages(pdf_path) as pdf:
        fig, ax = plt.subplots(figsize=(11, 7.5))
        ax.axis("off")
        text = (
            "Genome-wide reciprocal-mixing inference\n\n"
            f"Autosomal fit sites: {metadata['autosomal_fit_sites']:,}\n"
            f"Controls: {metadata['controls']}\n"
            f"Founder simulations: {selection['simulation_replicates']}\n"
            f"5-Mb delete-one-block fits: {selection['jackknife_blocks']}\n\n"
            f"Selected model: {selection['selected_founder_model'].replace('_', ' ')}\n"
            f"Calibrated P(four-founder): {selection['calibrated_posterior_four_founder']:.3f}\n"
            f"Right material in Left: {100*selection['selected_left_calibrated_mixing']:.1f}%\n"
            f"Left material in Right: {100*selection['selected_right_calibrated_mixing']:.1f}%\n\n"
            f"Lower normalized Z depth: {metadata['inferred_lower_Z_depth_side']} side\n"
            "Chromosomes 1-20 determine founder model and mixing.\n"
            "Chromosome 21/Z is retained as a held-out sex-ploidy diagnostic."
        )
        ax.text(0.05, 0.95, text, va="top", ha="left", fontsize=14, linespacing=1.45)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        calibration_data, phi = load_npz(FIT / "calibration_fit_input.npz")
        rng = np.random.default_rng(20260912)
        observed_left = calibration_data["left_alt"] / calibration_data["left_dp"]
        observed_right = calibration_data["right_alt"] / calibration_data["right_dp"]

        control_called = np.mean(calibration_data["controls"] >= 0, axis=0)
        control_het = np.array([
            np.mean(column[column >= 0] == 1)
            for column in calibration_data["controls"].T
        ])
        eligible_controls = np.where(control_called >= 0.90)[0]
        median_het = np.median(control_het[eligible_controls])
        control_index = int(
            eligible_controls[np.argmin(np.abs(control_het[eligible_controls] - median_het))]
        )
        control_alt = calibration_data["controls_alt"][:, control_index]
        control_dp = calibration_data["controls_dp"][:, control_index]
        usable_control = (control_alt >= 0) & (control_dp >= 10) & (control_alt <= control_dp)
        split_dp_left = control_dp[usable_control] // 2
        split_dp_right = control_dp[usable_control] - split_dp_left
        split_alt_left = np.array([
            rng.hypergeometric(int(alt), int(dp - alt), int(draws))
            for alt, dp, draws in zip(
                control_alt[usable_control], control_dp[usable_control], split_dp_left
            )
        ])
        split_alt_right = control_alt[usable_control] - split_alt_left
        unmixed_left = split_alt_left / split_dp_left
        unmixed_right = split_alt_right / split_dp_right

        controls = calibration_data["controls"]
        q = calibration_data["q"]
        donors = rng.choice(controls.shape[1], size=3, replace=False)
        donor_a = impute_genotypes(controls[:, donors[0]], q, rng)
        donor_b = impute_genotypes(controls[:, donors[1]], q, rng)
        donor_c = impute_genotypes(controls[:, donors[2]], q, rng)
        shared = sampled_allele(donor_a, rng)
        three_left_gt = (shared + sampled_allele(donor_b, rng)) / 2.0
        three_right_gt = (shared + sampled_allele(donor_c, rng)) / 2.0
        four_left_gt = donor_a / 2.0
        four_right_gt = donor_b / 2.0
        three_left_alt, three_right_alt = simulate_counts(
            three_left_gt, three_right_gt,
            calibration_data["n_left"], calibration_data["n_right"],
            0.20, 0.20, phi, rng,
        )
        four_left_alt, four_right_alt = simulate_counts(
            four_left_gt, four_right_gt,
            calibration_data["n_left"], calibration_data["n_right"],
            0.20, 0.20, phi, rng,
        )
        panels = [
            (observed_left, observed_right, "Observed merged Left and Right"),
            (unmixed_left, unmixed_right, f"Unmixed control {control_index + 1}, read-split"),
            (three_left_alt / calibration_data["n_left"], three_right_alt / calibration_data["n_right"], "Simulated three-founder, 20%/20%"),
            (four_left_alt / calibration_data["n_left"], four_right_alt / calibration_data["n_right"], "Simulated four-founder, 20%/20%"),
        ]
        fig, axes = plt.subplots(2, 2, figsize=(11, 10), sharex=True, sharey=True)
        for axis, (x_values, y_values, title) in zip(axes.flat, panels):
            axis.hexbin(
                x_values, y_values, gridsize=45, extent=(0, 1, 0, 1),
                mincnt=1, bins="log", cmap="viridis",
            )
            axis.plot([0, 1], [0, 1], color="white", linewidth=1.0, alpha=0.8)
            axis.set_title(title)
            axis.set_xlabel("Left-side ALT-read fraction")
            axis.set_ylabel("Right-side ALT-read fraction")
        fig.suptitle("Observed and benchmark allele-fraction patterns", fontsize=15)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        fig, axes = plt.subplots(1, 2, figsize=(13, 5.8), sharex=True, sharey=True)
        contrast_codes = pd.qcut(
            simulations.genotype_contrast,
            3,
            labels=False,
            duplicates="drop",
        )
        contrast_names = {0: "low", 1: "middle", 2: "high"}
        simulations = simulations.assign(
            contrast_bin=contrast_codes.map(contrast_names).fillna("middle")
        )
        colors = {"low": blue, "middle": green, "high": orange}
        for axis, model in zip(axes, ["three_founder", "four_founder"]):
            subset = simulations[simulations.true_model == model]
            estimate = "estimated_c_left_three" if model == "three_founder" else "estimated_c_left_four"
            for contrast, group in subset.groupby("contrast_bin"):
                axis.scatter(100 * group.true_c_left, 100 * group[estimate], s=13, alpha=0.28, color=colors.get(contrast, grey), label=f"{contrast} contrast")
                trend = group.groupby("true_c_left", as_index=False)[estimate].median().sort_values("true_c_left")
                axis.plot(100 * trend.true_c_left, 100 * trend[estimate], color=colors.get(contrast, grey), linewidth=2)
            axis.plot([0, 30], [0, 30], linestyle="--", color="black", linewidth=1)
            axis.set_title(model.replace("_", " ").title())
            axis.set_xlabel("True mixing into Left (%)")
            axis.grid(alpha=0.18)
        axes[0].set_ylabel("Recovered mixing (%)")
        axes[1].legend(frameon=False, fontsize=9)
        fig.suptitle("Simulation-calibrated recovery by genotype contrast", fontsize=15)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(10, 6))
        for model, color in [("three_founder", blue), ("four_founder", orange)]:
            values = simulations.loc[simulations.true_model == model, "delta_logL_four_minus_three"]
            ax.hist(values, bins=35, density=True, histtype="step", linewidth=2.2, color=color, label=f"Simulated {model.replace('_', ' ')}")
        observed_cal = fits[fits.scope == "calibration_subset"].set_index("model")
        observed_delta = observed_cal.loc["four_founder", "log_likelihood"] - observed_cal.loc["three_founder", "log_likelihood"]
        ax.axvline(observed_delta, color="black", linewidth=2.3, label="Observed")
        ax.set_xlabel("log L(four-founder) - log L(three-founder)")
        ax.set_ylabel("Simulation density")
        ax.set_title("Calibrated support for the founder model")
        ax.legend(frameon=False)
        ax.grid(alpha=0.18)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
        selected = selection["selected_founder_model"]
        short = "four" if selected == "four_founder" else "three"
        for axis, side in zip(axes, ["left", "right"]):
            values = 100 * jackknife[f"c_{side}_{short}"]
            axis.scatter(jackknife.block_id, values, c=jackknife.molecule, cmap="viridis", s=28)
            raw = 100 * summary.loc[(summary.model == selected) & (summary.side == side), "raw_mixing"].iloc[0]
            axis.axhline(raw, color="black", linestyle="--", linewidth=1.5)
            axis.set_title(f"Mixing into {side.title()}")
            axis.set_xlabel("Omitted 5-Mb block")
            axis.set_ylabel("Leave-one-block-out raw estimate (%)")
            axis.grid(alpha=0.18)
        fig.suptitle("Block-jackknife stability", fontsize=15)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        lengths = windows.groupby("molecule").window_end.max().sort_index()
        offsets = {}
        running = 0
        for molecule, length in lengths.items():
            offsets[int(molecule)] = running
            running += int(length)
        windows["genome_position"] = windows.apply(lambda row: offsets[int(row.molecule)] + (row.window_start + row.window_end) / 2, axis=1)
        fig, axes = plt.subplots(4, 1, figsize=(15, 13), sharex=True)
        alternating = [blue, grey]
        for molecule, group in windows.groupby("molecule"):
            color = purple if int(molecule) == 21 else alternating[(int(molecule) - 1) % 2]
            axes[0].scatter(group.genome_position, group.delta_logL_four_minus_three_per_site, s=28, color=color)
            axes[1].scatter(group.genome_position, group.delta_logL_sex_aware_minus_four_per_site, s=28, color=color)
            axes[2].scatter(group.genome_position, group.mean_left_dp / group.mean_right_dp, s=28, color=color)
        centers = [offsets[int(m)] + lengths.loc[m] / 2 for m in lengths.index]
        block_evidence = pd.read_csv(
            FINAL / "jackknife_block_evidence_manhattan.tsv", sep="\t"
        )
        block_evidence["genome_position"] = block_evidence.apply(
            lambda row: offsets[int(row.molecule)]
            + (row.block_start + row.block_end) / 2,
            axis=1,
        )
        for molecule, group in block_evidence.groupby("molecule"):
            color = alternating[(int(molecule) - 1) % 2]
            axes[3].scatter(
                group.genome_position,
                group.removed_block_delta_logL_influence_per_site,
                s=34,
                color=color,
            )
        for axis in [axes[0], axes[1], axes[3]]:
            axis.axhline(0, color="black", linewidth=0.8)
            axis.grid(axis="y", alpha=0.18)
        axes[2].axhline(1, color="black", linewidth=0.8)
        axes[2].grid(axis="y", alpha=0.18)
        axes[0].set_ylabel("Four - three\nlog L per site")
        axes[1].set_ylabel("Sex-aware - diploid\nlog L per site")
        axes[2].set_ylabel("Mean DP Left / Right")
        axes[3].set_ylabel("5-Mb block influence\non four - three log L/site")
        axes[3].set_xticks(centers)
        axes[3].set_xticklabels([str(int(m)) for m in lengths.index])
        axes[3].set_xlabel("Chromosome (21 = Z; magenta in Z panels)")
        fig.suptitle(
            "Broad genomic evidence, jackknife influence, and held-out Z diagnostic",
            fontsize=15,
        )
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

    png_path = FIG / "genomewide_founder_mixing_summary.png"
    selected_rows = summary[summary.model == selection["selected_founder_model"]]
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    x = np.arange(2)
    estimates = 100 * selected_rows.set_index("side").loc[["left", "right"], "calibrated_mixing"].to_numpy()
    lower = 100 * selected_rows.set_index("side").loc[["left", "right"], "ci95_lower"].to_numpy()
    upper = 100 * selected_rows.set_index("side").loc[["left", "right"], "ci95_upper"].to_numpy()
    ax.errorbar(x, estimates, yerr=np.vstack((estimates - lower, upper - estimates)), fmt="o", markersize=9, capsize=5, color=blue)
    ax.set_xticks(x, ["Right material in Left", "Left material in Right"])
    ax.set_ylabel("Calibrated mixing (%)")
    ax.set_title(selection["selected_founder_model"].replace("_", " ").title())
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(png_path, dpi=300)
    plt.close(fig)
    audit = pd.DataFrame({
        "item": ["narrative_pdf", "summary_png", "selected_model", "simulation_replicates", "jackknife_blocks", "window_count", "focal_GT_used"],
        "value": [str(pdf_path), str(png_path), selection["selected_founder_model"], len(simulations), len(jackknife), len(windows), False],
    })
    audit.to_csv(FIG / "genomewide_narrative.audit.tsv", sep="\t", index=False)
    (FIG / "plotting.complete").write_text("complete\n")
    print("GENOME-WIDE NARRATIVE FIGURES COMPLETE")
    print(audit.to_string(index=False))


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("fit")
    simulation = subparsers.add_parser("simulate")
    simulation.add_argument("--task", type=int, required=True)
    jackknife = subparsers.add_parser("jackknife")
    jackknife.add_argument("--task", type=int, required=True)
    subparsers.add_parser("collect")
    adjustment = subparsers.add_parser("adjust")
    adjustment.add_argument("--task", type=int, required=True)
    subparsers.add_parser("plot")
    arguments = parser.parse_args()
    commands = {
        "fit": command_fit,
        "simulate": command_simulate,
        "jackknife": command_jackknife,
        "collect": command_collect,
        "adjust": command_adjust,
        "plot": command_plot,
    }
    commands[arguments.command](arguments)


if __name__ == "__main__":
    main()
