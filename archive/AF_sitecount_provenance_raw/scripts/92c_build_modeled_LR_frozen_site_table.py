#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path("/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling")
CHECKPOINTS = PROJECT / "analysis/original_502/checkpoints"
SITES_FILE = CHECKPOINTS / "frozen_comparative_sites_reference_polarized.tsv"
MIXTURE_FILE = (
    PROJECT
    / "analysis/original_502/reciprocal_L_R_mixture"
    / "optix_reciprocal_L_R_genotype_pairs.posterior95.tsv"
)
OUTPUT = CHECKPOINTS / "focal_calls_at_frozen_sites_biological.modeled_LR_posterior95.wide.tsv"
AUDIT = CHECKPOINTS / "focal_calls_at_frozen_sites_biological.modeled_LR_posterior95.audit.tsv"

VALID_GT = {"0/0", "0/1", "1/1"}


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "t", "1", "yes", "y"})


def canonical_gt(value: object) -> str | None:
    text = str(value).strip().replace("|", "/")
    if text == "1/0":
        text = "0/1"
    return text if text in VALID_GT else None


def biological_call(gt: str | None, presence_is_ref: bool) -> tuple[str, float]:
    if gt is None:
        return "missing", np.nan
    alt_count = {"0/0": 0, "0/1": 1, "1/1": 2}[gt]
    presence_dosage = (2 - alt_count) / 2 if presence_is_ref else alt_count / 2
    label = {0.0: "AA", 0.5: "PA", 1.0: "PP"}[presence_dosage]
    return label, presence_dosage


def format_depth(ref_depth: object, alt_depth: object) -> str:
    ref_value = pd.to_numeric(pd.Series([ref_depth]), errors="coerce").iloc[0]
    alt_value = pd.to_numeric(pd.Series([alt_depth]), errors="coerce").iloc[0]
    if not np.isfinite(ref_value) or not np.isfinite(alt_value):
        return "missing"
    return f"{int(ref_value)},{int(alt_value)}"


def main() -> None:
    for path in (SITES_FILE, MIXTURE_FILE):
        if not path.exists():
            raise FileNotFoundError(path)

    sites = pd.read_csv(SITES_FILE, sep="\t", dtype=str, keep_default_na=False)
    mixture = pd.read_csv(MIXTURE_FILE, sep="\t", dtype=str, keep_default_na=False)

    site_required = {
        "site_id", "module", "module_display", "block", "chrom", "pos", "ref", "alt",
        "minimum_comparative_BA", "presence_associated_allele", "absence_associated_allele",
        "biological_orientation", "left_expected_genotype", "right_expected_genotype",
    }
    mixture_required = {
        "chrom", "pos", "ref", "alt", "mixture_left_GT", "mixture_right_GT",
        "pair_posterior", "posterior95", "left_ref_depth", "left_alt_depth", "left_DP",
        "left_observed_AF", "right_ref_depth", "right_alt_depth", "right_DP",
        "right_observed_AF",
    }
    missing_sites = site_required - set(sites.columns)
    missing_mixture = mixture_required - set(mixture.columns)
    if missing_sites:
        raise RuntimeError(f"Sites table missing columns: {sorted(missing_sites)}")
    if missing_mixture:
        raise RuntimeError(f"Mixture table missing columns: {sorted(missing_mixture)}")

    key = ["chrom", "pos", "ref", "alt"]
    if sites.duplicated(key).any():
        raise RuntimeError("Frozen-site table contains duplicate physical SNPs")
    if mixture.duplicated(key).any():
        raise RuntimeError("Mixture table contains duplicate physical SNPs")

    mixture_columns = key + [
        "mixture_left_GT", "mixture_right_GT", "pair_posterior", "posterior95",
        "left_ref_depth", "left_alt_depth", "left_DP", "left_observed_AF",
        "right_ref_depth", "right_alt_depth", "right_DP", "right_observed_AF",
    ]
    joined = sites.merge(
        mixture[mixture_columns], on=key, how="left", validate="one_to_one", indicator=True
    )
    if len(joined) != len(sites):
        raise RuntimeError("Mixture join changed the frozen-site row count")
    matched_mixture_rows = int(joined["_merge"].eq("both").sum())

    joined["pair_posterior_numeric"] = pd.to_numeric(joined["pair_posterior"], errors="coerce")
    joined["left_model_GT"] = joined["mixture_left_GT"].map(canonical_gt)
    joined["right_model_GT"] = joined["mixture_right_GT"].map(canonical_gt)
    joined["model_confident"] = (
        as_bool(joined["posterior95"])
        & joined["pair_posterior_numeric"].ge(0.95)
        & joined["left_model_GT"].notna()
        & joined["right_model_GT"].notna()
    )

    records: list[dict[str, object]] = []
    for row in joined.itertuples(index=False):
        presence = str(row.presence_associated_allele)
        ref = str(row.ref)
        alt = str(row.alt)
        if presence not in {ref, alt}:
            raise RuntimeError(
                f"{row.site_id}: presence allele {presence!r} is neither REF {ref!r} nor ALT {alt!r}"
            )
        presence_is_ref = presence == ref
        confident = bool(row.model_confident)
        left_gt = row.left_model_GT if confident else None
        right_gt = row.right_model_GT if confident else None
        left_bio, left_dose = biological_call(left_gt, presence_is_ref)
        right_bio, right_dose = biological_call(right_gt, presence_is_ref)

        left_alt_fraction = pd.to_numeric(pd.Series([row.left_observed_AF]), errors="coerce").iloc[0]
        right_alt_fraction = pd.to_numeric(pd.Series([row.right_observed_AF]), errors="coerce").iloc[0]
        left_presence_fraction = (
            (1 - left_alt_fraction) if presence_is_ref else left_alt_fraction
        ) if confident and np.isfinite(left_alt_fraction) else np.nan
        right_presence_fraction = (
            (1 - right_alt_fraction) if presence_is_ref else right_alt_fraction
        ) if confident and np.isfinite(right_alt_fraction) else np.nan

        left_expected = str(row.left_expected_genotype)
        right_expected = str(row.right_expected_genotype)
        left_match = confident and left_bio == left_expected
        right_match = confident and right_bio == right_expected

        records.append(
            {
                "site_id": row.site_id,
                "module": row.module,
                "module_display": row.module_display,
                "block": row.block,
                "chrom": row.chrom,
                "pos": row.pos,
                "ref": ref,
                "alt": alt,
                "minimum_comparative_BA": row.minimum_comparative_BA,
                "presence_associated_allele": presence,
                "absence_associated_allele": row.absence_associated_allele,
                "biological_orientation": row.biological_orientation,
                "source_vcf": "reciprocal-mixture model from merged TN_Left/TN_Right counts",
                "left_GT_raw": left_gt if left_gt is not None else "missing",
                "left_AD_raw": format_depth(row.left_ref_depth, row.left_alt_depth),
                "left_DP": pd.to_numeric(pd.Series([row.left_DP]), errors="coerce").iloc[0],
                "left_GQ": np.nan,
                "left_model_posterior95": confident,
                "left_biological_genotype": left_bio,
                "left_presence_dosage": left_dose,
                "left_presence_allele_fraction": left_presence_fraction,
                "left_expected_genotype": left_expected,
                "left_matches_expectation": left_match if confident else np.nan,
                "left_focal_call_status": (
                    "reciprocal_mixture_posterior95" if confident else "no_posterior95_model_call"
                ),
                "right_GT_raw": right_gt if right_gt is not None else "missing",
                "right_AD_raw": format_depth(row.right_ref_depth, row.right_alt_depth),
                "right_DP": pd.to_numeric(pd.Series([row.right_DP]), errors="coerce").iloc[0],
                "right_GQ": np.nan,
                "right_model_posterior95": confident,
                "right_biological_genotype": right_bio,
                "right_presence_dosage": right_dose,
                "right_presence_allele_fraction": right_presence_fraction,
                "right_expected_genotype": right_expected,
                "right_matches_expectation": right_match if confident else np.nan,
                "right_focal_call_status": (
                    "reciprocal_mixture_posterior95" if confident else "no_posterior95_model_call"
                ),
                "both_callable": confident,
                "pair_matches_expectation": left_match and right_match if confident else np.nan,
                "mixture_pair_posterior": row.pair_posterior_numeric,
                "focal_call_source": (
                    "reciprocal_mixture_posterior95_only" if confident else "no_posterior95_call"
                ),
            }
        )

    output = pd.DataFrame(records)
    if len(output) != 117 or output["site_id"].duplicated().any():
        raise RuntimeError("Final modeled table must contain 117 unique frozen sites")
    if output.filter(regex=r"T2|T3|TN2|TN3", axis=1).shape[1] != 0:
        raise RuntimeError("Unexpected T2/T3-specific columns in modeled table")
    forbidden = {
        "left_original_GT", "right_original_GT", "original_pair_GT",
        "left_QC_pass_DP2_GQ20", "right_QC_pass_DP2_GQ20",
    }
    present_forbidden = forbidden & set(output.columns)
    if present_forbidden:
        raise RuntimeError(
            f"Forbidden original-GT or legacy focal-QC columns present: {sorted(present_forbidden)}"
        )
    called = output["both_callable"].astype(bool)
    if called.any() and not output.loc[called, "focal_call_source"].eq(
        "reciprocal_mixture_posterior95_only"
    ).all():
        raise RuntimeError("A called focal row does not use the posterior95 mixture model")

    output.to_csv(OUTPUT, sep="\t", index=False, na_rep="NA")
    audit = (
        output.groupby("module", sort=False)
        .agg(
            frozen_sites=("site_id", "size"),
            posterior95_pairs=("both_callable", "sum"),
            expected_pairs=("pair_matches_expectation", lambda x: int(pd.Series(x).fillna(False).sum())),
        )
        .reset_index()
    )
    audit["unavailable_pairs"] = audit["frozen_sites"] - audit["posterior95_pairs"]
    audit["unexpected_pairs"] = audit["posterior95_pairs"] - audit["expected_pairs"]
    audit.to_csv(AUDIT, sep="\t", index=False)

    print("MODELED FROZEN-SITE TABLE COMPLETE")
    print(f"Output: {OUTPUT}")
    print(f"Audit: {AUDIT}")
    print(f"Rows: {len(output)}")
    print(f"Frozen sites found in mixture table: {matched_mixture_rows}")
    print(f"Posterior95 modeled pairs: {int(called.sum())}")
    print("\nBY MODULE")
    print(audit.to_csv(sep="\t", index=False).rstrip())


if __name__ == "__main__":
    main()
