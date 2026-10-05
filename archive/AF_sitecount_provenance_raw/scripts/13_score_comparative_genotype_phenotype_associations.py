#!/usr/bin/env python3

from pathlib import Path
import math
import re

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")
GENOTYPE_DIR = PROJECT / "data/derived/region_genotypes"
METADATA_FILE = PROJECT / "data/derived/region_panel_metadata.tsv"
CHECKPOINT_DIR = PROJECT / "analysis/checkpoints"

FOCAL = {
    "TN2L_S1", "TN2R_S2", "TN3L_S3",
    "TN3R_S4", "TN_Left", "TN_Right",
}

CONTROL_PANELS = [
    "neutral1", "neutral2", "neutral3", "neutral4",
]

TRAITS = {
    "B": {
        "target_panel": "optix",
        "figure": "Optix Figure",
        "phenotype": "B",
        "model": "dominant",
        "reference_expectation": "REF-associated",
    },
    "R": {
        "target_panel": "optix",
        "figure": "Optix Figure",
        "phenotype": "R",
        "model": "dominant",
        "reference_expectation": "ALT-associated",
    },
    "D": {
        "target_panel": "optix",
        "figure": "Optix Figure",
        "phenotype": "D",
        "model": "dominant",
        "reference_expectation": "ALT-associated",
    },
    "Yb_dorsal": {
        "target_panel": "ivory",
        "figure": "Ivory Figure",
        "phenotype": "Yb-dorsal",
        "model": "recessive",
        "reference_expectation": "REF-associated",
    },
    "Yb_ventral": {
        "target_panel": "ivory",
        "figure": "Ivory Figure",
        "phenotype": "Yb-ventral",
        "model": "recessive",
        "reference_expectation": "REF-associated",
    },
}

MIN_DP = 2
MIN_GQ = 20
MIN_STATE_CALL_RATE = 0.50
MIN_STATE_CALLED = 20
MIN_ASSAY_CALL_RATE = 0.50
MIN_PREDICTOR_CLASS = 10


def normalized(value):
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def resolve_column(dataframe, requested):
    target = normalized(requested)
    matches = [
        column for column in dataframe.columns
        if normalized(column) == target
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Could not resolve {requested!r}; matches={matches}"
        )
    return matches[0]


def numeric(series):
    return pd.to_numeric(
        series.replace({"": np.nan, ".": np.nan}),
        errors="coerce",
    )


def alt_dosage(value):
    value = str(value)
    if value in {"", ".", "./.", ".|.", "nan"}:
        return np.nan

    fields = re.split(r"[/|]", value)
    if len(fields) != 2 or any(x == "." for x in fields):
        return np.nan

    try:
        alleles = [int(x) for x in fields]
    except ValueError:
        return np.nan

    if any(x not in {0, 1} for x in alleles):
        return np.nan

    return float(sum(alleles))


def model_metrics(phenotype, predictor):
    phenotype = np.asarray(phenotype, dtype=int)
    predictor = np.asarray(predictor, dtype=int)

    a = int(((phenotype == 1) & (predictor == 1)).sum())
    b = int(((phenotype == 0) & (predictor == 1)).sum())
    c = int(((phenotype == 1) & (predictor == 0)).sum())
    d = int(((phenotype == 0) & (predictor == 0)).sum())

    sensitivity = a / (a + c) if (a + c) else np.nan
    specificity = d / (b + d) if (b + d) else np.nan
    balanced_accuracy = (
        (sensitivity + specificity) / 2
        if np.isfinite(sensitivity) and np.isfinite(specificity)
        else np.nan
    )

    denominator = math.sqrt(
        (a + b) * (a + c) * (d + b) * (d + c)
    )
    mcc = (
        ((a * d) - (b * c)) / denominator
        if denominator else np.nan
    )

    odds_ratio, fisher_p = fisher_exact(
        [[a, b], [c, d]],
        alternative="two-sided",
    )

    corrected_or = (
        ((a + 0.5) * (d + 0.5))
        / ((b + 0.5) * (c + 0.5))
    )

    return {
        "a_true_positive": a,
        "b_false_positive": b,
        "c_false_negative": c,
        "d_true_negative": d,
        "predictor_1_n": a + b,
        "predictor_0_n": c + d,
        "sensitivity": sensitivity,
        "specificity": specificity,
        "balanced_accuracy": balanced_accuracy,
        "MCC": mcc,
        "fisher_odds_ratio": odds_ratio,
        "corrected_odds_ratio": corrected_or,
        "fisher_p": fisher_p,
    }


def predictor_from_dosage(dosage, model, orientation):
    if model == "dominant":
        if orientation == "ALT-associated":
            return dosage >= 1
        return dosage <= 1

    if model == "recessive":
        if orientation == "ALT-associated":
            return dosage == 2
        return dosage == 0

    raise ValueError(model)


def bh_adjust(pvalues):
    pvalues = np.asarray(pvalues, dtype=float)
    result = np.full(len(pvalues), np.nan)

    valid = np.isfinite(pvalues)
    valid_indices = np.where(valid)[0]
    valid_p = pvalues[valid]

    if not len(valid_p):
        return result

    order = np.argsort(valid_p)
    ranked = valid_p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)

    reordered = np.empty(len(valid_p))
    reordered[order] = adjusted
    result[valid_indices] = reordered

    return result


def read_panel(panel, metadata_by_sample):
    path = GENOTYPE_DIR / f"{panel}.genotypes.tsv.gz"
    calls = pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        compression="gzip",
        keep_default_na=False,
    )

    calls["data_type"] = calls["sample"].map(
        metadata_by_sample["data_type_qc"]
    )
    calls["role"] = calls["sample"].map(
        metadata_by_sample["role_qc"]
    )
    calls["DP_numeric"] = numeric(calls["DP"])
    calls["GQ_numeric"] = numeric(calls["GQ"])
    calls["dosage"] = calls["GT"].map(alt_dosage)

    calls["biallelic_SNP"] = (
        ~calls["alt"].str.contains(",", regex=False)
        & calls["ref"].str.len().eq(1)
        & calls["alt"].str.len().eq(1)
    )

    calls["passes_qc"] = (
        calls["dosage"].notna()
        & calls["DP_numeric"].ge(MIN_DP)
        & calls["GQ_numeric"].ge(MIN_GQ)
    )

    return calls


def score_site(site, trait_name, config):
    phenotype_column = f"{trait_name}_qc"
    phenotype = site[phenotype_column].astype(int)

    state_statistics = {}
    for state in [0, 1]:
        state_rows = site.loc[phenotype.eq(state)]
        state_statistics[state] = {
            "total": int(len(state_rows)),
            "called": int(state_rows["passes_qc"].sum()),
            "rate": float(state_rows["passes_qc"].mean()),
        }

    assay_statistics = {}
    for assay in ["CS", "WGS"]:
        assay_rows = site.loc[site["data_type"].eq(assay)]
        assay_statistics[assay] = {
            "total": int(len(assay_rows)),
            "called": int(assay_rows["passes_qc"].sum()),
            "rate": (
                float(assay_rows["passes_qc"].mean())
                if len(assay_rows) else np.nan
            ),
        }

    callability_eligible = (
        all(
            state_statistics[state]["called"] >= MIN_STATE_CALLED
            and state_statistics[state]["rate"] >= MIN_STATE_CALL_RATE
            for state in [0, 1]
        )
        and all(
            assay_statistics[assay]["rate"] >= MIN_ASSAY_CALL_RATE
            for assay in ["CS", "WGS"]
        )
    )

    called = site.loc[site["passes_qc"]].copy()

    missingness_table = [
        [
            state_statistics[1]["called"],
            state_statistics[1]["total"]
            - state_statistics[1]["called"],
        ],
        [
            state_statistics[0]["called"],
            state_statistics[0]["total"]
            - state_statistics[0]["called"],
        ],
    ]
    missingness_p = fisher_exact(
        missingness_table,
        alternative="two-sided",
    ).pvalue

    orientation_results = {}

    for orientation in ["REF-associated", "ALT-associated"]:
        predictor = predictor_from_dosage(
            called["dosage"].to_numpy(),
            config["model"],
            orientation,
        ).astype(int)

        metrics = model_metrics(
            called[phenotype_column].astype(int),
            predictor,
        )
        metrics["orientation"] = orientation
        metrics["valid_predictor_classes"] = (
            metrics["predictor_1_n"] >= MIN_PREDICTOR_CLASS
            and metrics["predictor_0_n"] >= MIN_PREDICTOR_CLASS
        )
        orientation_results[orientation] = metrics

    valid_orientations = [
        result for result in orientation_results.values()
        if result["valid_predictor_classes"]
    ]

    if valid_orientations:
        selected = sorted(
            valid_orientations,
            key=lambda x: (
                -x["balanced_accuracy"],
                x["fisher_p"],
            ),
        )[0]
    else:
        selected = max(
            orientation_results.values(),
            key=lambda x: (
                x["balanced_accuracy"]
                if np.isfinite(x["balanced_accuracy"])
                else -np.inf
            ),
        )

    raw_orientation_p = min(
        result["fisher_p"]
        for result in orientation_results.values()
    )
    orientation_corrected_p = min(1.0, 2 * raw_orientation_p)

    row = {
        "trait": trait_name,
        "model": config["model"],
        "reference_expected_orientation":
            config["reference_expectation"],
        "chrom": site["chrom"].iloc[0],
        "pos": int(site["pos"].iloc[0]),
        "ref": site["ref"].iloc[0],
        "alt": site["alt"].iloc[0],
        "phenotype_0_total": state_statistics[0]["total"],
        "phenotype_1_total": state_statistics[1]["total"],
        "phenotype_0_called": state_statistics[0]["called"],
        "phenotype_1_called": state_statistics[1]["called"],
        "phenotype_0_call_rate": state_statistics[0]["rate"],
        "phenotype_1_call_rate": state_statistics[1]["rate"],
        "CS_call_rate": assay_statistics["CS"]["rate"],
        "WGS_call_rate": assay_statistics["WGS"]["rate"],
        "missingness_fisher_p": missingness_p,
        "callability_eligible": callability_eligible,
        "selected_orientation": selected["orientation"],
        "selected_matches_reference_expectation": (
            selected["orientation"]
            == config["reference_expectation"]
        ),
        "predictor_classes_eligible":
            selected["valid_predictor_classes"],
        "score_eligible": (
            callability_eligible
            and selected["valid_predictor_classes"]
        ),
        "orientation_corrected_fisher_p":
            orientation_corrected_p,
    }

    for key, value in selected.items():
        if key not in {"orientation", "valid_predictor_classes"}:
            row[key] = value

    for orientation, result in orientation_results.items():
        prefix = (
            "REF" if orientation == "REF-associated" else "ALT"
        )
        row[f"{prefix}_balanced_accuracy"] = (
            result["balanced_accuracy"]
        )
        row[f"{prefix}_fisher_p"] = result["fisher_p"]
        row[f"{prefix}_predictor_1_n"] = (
            result["predictor_1_n"]
        )
        row[f"{prefix}_predictor_0_n"] = (
            result["predictor_0_n"]
        )

    for assay in ["CS", "WGS"]:
        assay_called = called.loc[
            called["data_type"].eq(assay)
        ]

        if len(assay_called):
            assay_predictor = predictor_from_dosage(
                assay_called["dosage"].to_numpy(),
                config["model"],
                selected["orientation"],
            ).astype(int)

            assay_metrics = model_metrics(
                assay_called[phenotype_column].astype(int),
                assay_predictor,
            )
            row[f"{assay}_balanced_accuracy"] = (
                assay_metrics["balanced_accuracy"]
            )
            row[f"{assay}_MCC"] = assay_metrics["MCC"]
            row[f"{assay}_fisher_p"] = (
                assay_metrics["fisher_p"]
            )
        else:
            row[f"{assay}_balanced_accuracy"] = np.nan
            row[f"{assay}_MCC"] = np.nan
            row[f"{assay}_fisher_p"] = np.nan

    return row


def main():
    CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(
        METADATA_FILE,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    accession = resolve_column(metadata, "Accession")
    data_type = resolve_column(metadata, "Data Type")
    taxon = resolve_column(metadata, "Taxon Name")

    metadata["sample"] = metadata[accession].astype(str)
    metadata["data_type_qc"] = metadata[data_type].astype(str)
    metadata["taxon_qc"] = metadata[taxon].astype(str)
    metadata["role_qc"] = "comparative"

    metadata.loc[
        metadata["taxon_qc"].str.contains(
            "piano keys", case=False, na=False
        ),
        "role_qc",
    ] = "piano_keys"

    metadata.loc[
        metadata["sample"].isin(FOCAL),
        "role_qc",
    ] = "focal"

    for trait_name, config in TRAITS.items():
        phenotype_column = resolve_column(
            metadata,
            config["phenotype"],
        )
        figure_column = resolve_column(
            metadata,
            config["figure"],
        )

        metadata[f"{trait_name}_qc"] = numeric(
            metadata[phenotype_column]
        )
        metadata[f"{trait_name}_figure_qc"] = numeric(
            metadata[figure_column]
        )

    metadata_by_sample = metadata.set_index("sample")

    needed_panels = sorted(
        set(CONTROL_PANELS)
        | {config["target_panel"] for config in TRAITS.values()}
    )
    panel_calls = {
        panel: read_panel(panel, metadata_by_sample)
        for panel in needed_panels
    }

    results = []
    audit_rows = []

    for trait_name, config in TRAITS.items():
        primary_samples = set(
            metadata.loc[
                metadata["role_qc"].eq("comparative")
                & metadata[f"{trait_name}_figure_qc"].eq(1)
                & metadata[f"{trait_name}_qc"].isin([0, 1]),
                "sample",
            ]
        )

        for panel in [config["target_panel"]] + CONTROL_PANELS:
            calls = panel_calls[panel].copy()
            calls[f"{trait_name}_qc"] = calls["sample"].map(
                metadata_by_sample[f"{trait_name}_qc"]
            )

            cohort = calls.loc[
                calls["sample"].isin(primary_samples)
                & calls["biallelic_SNP"]
            ].copy()

            panel_rows = []

            for _, site in cohort.groupby(
                ["chrom", "pos", "ref", "alt"],
                sort=False,
            ):
                row = score_site(site, trait_name, config)
                row["panel"] = panel
                row["region_role"] = (
                    "target"
                    if panel == config["target_panel"]
                    else "neutral_control"
                )
                panel_rows.append(row)

            results.extend(panel_rows)

            audit_rows.append({
                "trait": trait_name,
                "panel": panel,
                "region_role": (
                    "target"
                    if panel == config["target_panel"]
                    else "neutral_control"
                ),
                "biallelic_SNPs": len(panel_rows),
                "callability_eligible": sum(
                    row["callability_eligible"]
                    for row in panel_rows
                ),
                "score_eligible": sum(
                    row["score_eligible"]
                    for row in panel_rows
                ),
            })

    results = pd.DataFrame(results)

    results["BH_q_within_panel_trait"] = np.nan
    for _, indices in results.groupby(
        ["panel", "trait"]
    ).groups.items():
        indices = list(indices)
        eligible_indices = [
            index for index in indices
            if results.loc[index, "score_eligible"]
        ]

        if eligible_indices:
            results.loc[
                eligible_indices,
                "BH_q_within_panel_trait",
            ] = bh_adjust(
                results.loc[
                    eligible_indices,
                    "orientation_corrected_fisher_p",
                ]
            )

    results["minus_log10_p"] = results[
        "orientation_corrected_fisher_p"
    ].map(
        lambda value: (
            -math.log10(max(float(value), 1e-300))
            if pd.notna(value) else np.nan
        )
    )

    results["rank_within_panel_trait"] = np.nan
    eligible = results["score_eligible"]

    results.loc[eligible, "rank_within_panel_trait"] = (
        results.loc[eligible]
        .groupby(["panel", "trait"])["balanced_accuracy"]
        .rank(method="min", ascending=False)
    )

    results_file = (
        CHECKPOINT_DIR
        / "comparative_genotype_phenotype_scores.tsv"
    )
    audit_file = (
        CHECKPOINT_DIR
        / "comparative_association_eligibility_summary.tsv"
    )
    specification_file = (
        CHECKPOINT_DIR
        / "candidate_discovery_specification.tsv"
    )

    results.to_csv(
        results_file,
        sep="\t",
        index=False,
        float_format="%.8g",
    )

    audit = pd.DataFrame(audit_rows)
    audit.to_csv(audit_file, sep="\t", index=False)

    specification = pd.DataFrame([
        ["candidate_samples", "Figure 1 comparative samples only"],
        ["excluded_samples", "six TN focal samples"],
        ["excluded_context", "five Butterfly World piano-key samples"],
        ["grouping_variables", "not used"],
        ["variant_type", "biallelic SNPs"],
        ["minimum_DP", MIN_DP],
        ["minimum_GQ", MIN_GQ],
        ["minimum_state_call_rate", MIN_STATE_CALL_RATE],
        ["minimum_state_called", MIN_STATE_CALLED],
        ["minimum_CS_call_rate", MIN_ASSAY_CALL_RATE],
        ["minimum_WGS_call_rate", MIN_ASSAY_CALL_RATE],
        ["minimum_predictor_class", MIN_PREDICTOR_CLASS],
        ["optix_model", "dominant allele-carrier"],
        ["ivory_model", "recessive allele-homozygote"],
        ["orientation_testing", "REF and ALT; two-model p correction"],
    ], columns=["parameter", "value"])

    specification.to_csv(
        specification_file,
        sep="\t",
        index=False,
    )

    top = (
        results.loc[
            results["region_role"].eq("target")
            & results["score_eligible"]
        ]
        .sort_values(
            [
                "trait",
                "balanced_accuracy",
                "orientation_corrected_fisher_p",
            ],
            ascending=[True, False, True],
        )
        .groupby("trait", sort=False)
        .head(8)
    )

    top_columns = [
        "trait", "panel", "chrom", "pos", "ref", "alt",
        "selected_orientation",
        "selected_matches_reference_expectation",
        "phenotype_0_called", "phenotype_1_called",
        "balanced_accuracy", "MCC",
        "sensitivity", "specificity",
        "orientation_corrected_fisher_p",
        "BH_q_within_panel_trait",
        "CS_balanced_accuracy", "WGS_balanced_accuracy",
    ]

    print("ASSOCIATION ELIGIBILITY")
    print(audit.to_csv(sep="\t", index=False).rstrip())

    print("\nTOP TARGET-REGION SCORES")
    print(
        top[top_columns].to_csv(
            sep="\t",
            index=False,
            float_format="%.5g",
        ).rstrip()
    )

    print("\nOutputs:")
    print(results_file)
    print(audit_file)
    print(specification_file)
    print(
        "\nFocal and piano-key genotypes were not used "
        "in any score."
    )


if __name__ == "__main__":
    main()
