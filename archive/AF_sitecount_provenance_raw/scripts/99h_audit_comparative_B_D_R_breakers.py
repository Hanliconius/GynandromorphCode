#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

CLASSIFICATION = (
    PROJECT
    / "analysis/original_502/checkpoints"
    / "B_DR_three_way_physical_SNP_classification.tsv"
)

GENOTYPES = (
    PROJECT
    / "analysis/original_502/checkpoints"
    / "core_explanatory_optix_genotypes.long.tsv"
)

METADATA = (
    PROJECT
    / "data/derived/region_panel_metadata.tsv"
)

OUTDIR = PROJECT / "analysis/original_502/checkpoints"

CALL_AUDIT = OUTDIR / "comparative_B_D_R_selected_site_individual_calls.tsv"
SITE_SUMMARY = OUTDIR / "comparative_B_D_R_selected_site_breaker_summary.tsv"
INDIVIDUAL_SUMMARY = OUTDIR / "comparative_B_D_R_individual_breaker_summary.tsv"
POPULATION_SUMMARY = OUTDIR / "comparative_B_D_R_population_breaker_summary.tsv"
GROUP_SUMMARY = OUTDIR / "comparative_B_D_R_diagnostic_group_summary.tsv"
HYPOTHESIS_SUMMARY = OUTDIR / "comparative_B_D_R_hypothesis_breaker_summary.tsv"


HYPOTHESIS_TO_MODULE = {
    "B": "Band",
    "D": "Dennis",
    "R": "Ray",
}

FINAL_CLASS = {
    "B": "B",
    "D": "DR",
    "R": "DR",
}

# Presence/absence expectations for the three diagnostic groups.
# B, D and R are dominant, so P/A and P/P both count as presence.
EXPECTED_PRESENCE = {
    "B": {
        "besckei_B": True,
        "meriana_D": False,
        "contigua_R": False,
    },
    "D": {
        "besckei_B": False,
        "meriana_D": True,
        "contigua_R": False,
    },
    "R": {
        "besckei_B": False,
        "meriana_D": False,
        "contigua_R": True,
    },
}

EXPECTED_GROUP_SIZES = {
    "besckei_B": 11,
    "meriana_D": 7,
    "contigua_R": 10,
}


def truthy(value):
    if pd.isna(value):
        return False
    return str(value).strip().lower() in {
        "true", "t", "1", "yes", "y"
    }


def numeric(value):
    return pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]


def clean_column_name(value):
    return " ".join(str(value).replace("\n", " ").split())


def qc_failure_reason(row):
    reasons = []

    gt = str(row.get("canonical_GT", "")).strip()
    if gt in {"", ".", "./.", ".|.", "nan", "None"}:
        reasons.append("missing_GT")

    dp = numeric(row.get("DP"))
    if pd.isna(dp):
        reasons.append("missing_DP")
    elif dp < 2:
        reasons.append("low_DP")

    gq = numeric(row.get("GQ"))
    if pd.isna(gq):
        reasons.append("missing_GQ")
    elif gq < 20:
        reasons.append("low_GQ")

    if not reasons:
        reasons.append("failed_DP2_GQ20_unspecified")

    return "+".join(reasons)


def joined_unique(values):
    cleaned = sorted({
        str(value)
        for value in values
        if pd.notna(value) and str(value).strip() not in {"", "nan"}
    })
    return ",".join(cleaned) if cleaned else "-"


def safe_fraction(numerator, denominator):
    if denominator == 0:
        return np.nan
    return numerator / denominator


def main():
    required = [CLASSIFICATION, GENOTYPES]

    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Required inputs are missing:\n" + "\n".join(missing)
        )

    classification = pd.read_csv(
        CLASSIFICATION,
        sep="\t",
        low_memory=False,
    )

    genotypes = pd.read_csv(
        GENOTYPES,
        sep="\t",
        low_memory=False,
    )

    # Convert the display labels used in the genotype table to the
    # internal diagnostic-group labels used by this audit.
    diagnostic_group_map = {
        "H. besckei": "besckei_B",
        "H. m. meriana": "meriana_D",
        "H. t. contigua": "contigua_R",
    }

    genotypes["diagnostic_group"] = (
        genotypes["group"].map(diagnostic_group_map)
    )

    # Select SNPs entirely from the comparative classification.
    selected = classification[
        classification["three_way_classification_pass"].map(truthy)
        & classification["best_hypothesis"].isin(["B", "D", "R"])
    ].copy()

    selected = selected.sort_values(
        ["best_hypothesis", "chrom", "pos"]
    ).reset_index(drop=True)

    duplicate_selected = selected.duplicated(
        ["chrom", "pos", "ref", "alt"],
        keep=False,
    )

    if duplicate_selected.any():
        duplicates = selected.loc[
            duplicate_selected,
            ["chrom", "pos", "ref", "alt", "best_hypothesis"],
        ]
        raise RuntimeError(
            "Comparative classification contains duplicated physical SNPs:\n"
            + duplicates.to_string(index=False)
        )

    # Canonical metadata is optional for the genotype audit but useful
    # for taxon/country summaries.
    metadata_lookup = {}

    if METADATA.is_file():
        metadata = pd.read_csv(
            METADATA,
            sep="\t",
            low_memory=False,
        )

        metadata.columns = [
            clean_column_name(column)
            for column in metadata.columns
        ]

        accession_column = next(
            (
                column for column in metadata.columns
                if column.lower() == "accession"
            ),
            None,
        )

        taxon_column = next(
            (
                column for column in metadata.columns
                if column.lower() == "taxon name"
            ),
            None,
        )

        country_column = next(
            (
                column for column in metadata.columns
                if column.lower() == "country"
            ),
            None,
        )

        if accession_column is not None:
            for _, row in metadata.iterrows():
                accession = row.get(accession_column)

                if pd.isna(accession):
                    continue

                metadata_lookup[str(accession)] = {
                    "taxon": (
                        row.get(taxon_column)
                        if taxon_column is not None
                        else np.nan
                    ),
                    "country": (
                        row.get(country_column)
                        if country_column is not None
                        else np.nan
                    ),
                }

    call_rows = []

    for _, site in selected.iterrows():
        hypothesis = str(site["best_hypothesis"])
        module = HYPOTHESIS_TO_MODULE[hypothesis]
        final_class = FINAL_CLASS[hypothesis]

        site_mask = (
            genotypes["module"].astype(str).eq(module)
            & genotypes["chrom"].astype(str).eq(str(site["chrom"]))
            & genotypes["pos"].astype(int).eq(int(site["pos"]))
            & genotypes["ref"].astype(str).eq(str(site["ref"]))
            & genotypes["alt"].astype(str).eq(str(site["alt"]))
            & genotypes["diagnostic_group"].isin(
                EXPECTED_GROUP_SIZES
            )
        )

        site_calls = genotypes.loc[site_mask].copy()

        duplicated_samples = site_calls.duplicated(
            ["sample"],
            keep=False,
        )

        if duplicated_samples.any():
            details = site_calls.loc[
                duplicated_samples,
                [
                    "module", "chrom", "pos", "sample",
                    "canonical_GT", "pattern_dosage",
                ],
            ]
            raise RuntimeError(
                "Duplicated sample calls were found:\n"
                + details.to_string(index=False)
            )

        expected_samples = (
            genotypes[
                genotypes["diagnostic_group"].isin(
                    EXPECTED_GROUP_SIZES
                )
            ][["sample", "diagnostic_group", "sample_label"]]
            .drop_duplicates("sample")
        )

        # This reindex ensures a missing genotype-table record is
        # represented explicitly as N rather than silently disappearing.
        site_calls = expected_samples.merge(
            site_calls,
            on=["sample", "diagnostic_group", "sample_label"],
            how="left",
            suffixes=("", "_genotype"),
        )

        for _, call in site_calls.iterrows():
            sample = str(call["sample"])
            group = str(call["diagnostic_group"])
            expected_present = EXPECTED_PRESENCE[hypothesis][group]

            metadata_values = metadata_lookup.get(
                sample,
                {"taxon": np.nan, "country": np.nan},
            )

            if pd.isna(call.get("chrom")):
                status = "N_missing_site_record"
                callable_call = False
                correct = False
                observed_present = np.nan
                biological_gt = np.nan
            else:
                qc_pass = truthy(call.get("DP2_GQ20_pass"))

                if not qc_pass:
                    status = "N_or_QC_" + qc_failure_reason(call)
                    callable_call = False
                    correct = False
                    observed_present = np.nan
                    biological_gt = call.get("biological_GT")
                else:
                    dosage = numeric(call.get("pattern_dosage"))

                    if pd.isna(dosage):
                        status = "N_missing_pattern_dosage"
                        callable_call = False
                        correct = False
                        observed_present = np.nan
                    else:
                        callable_call = True
                        observed_present = bool(dosage > 0)
                        correct = observed_present == expected_present
                        status = (
                            "correct"
                            if correct
                            else "callable_genotype_discordance"
                        )

                    biological_gt = call.get("biological_GT")

            call_rows.append({
                "final_class": final_class,
                "hypothesis": hypothesis,
                "module": module,
                "chrom": site["chrom"],
                "pos": int(site["pos"]),
                "ref": site["ref"],
                "alt": site["alt"],
                "best_macro_accuracy": numeric(
                    site.get("best_macro_accuracy")
                ),
                "score_margin": numeric(site.get("score_margin")),
                "comparative_BA": numeric(
                    site.get("best_hypothesis_comparative_BA")
                ),
                "sample": sample,
                "diagnostic_group": group,
                "taxon": metadata_values["taxon"],
                "country": metadata_values["country"],
                "expected_pattern_state": (
                    "carrier"
                    if expected_present
                    else "A/A"
                ),
                "observed_pattern_state": (
                    "carrier"
                    if observed_present is True
                    else (
                        "A/A"
                        if observed_present is False
                        else "N"
                    )
                ),
                "biological_GT": biological_gt,
                "canonical_GT": call.get("canonical_GT"),
                "AD": call.get("AD"),
                "DP": numeric(call.get("DP")),
                "GQ": numeric(call.get("GQ")),
                "pattern_dosage": numeric(
                    call.get("pattern_dosage")
                ),
                "callable": callable_call,
                "correct": correct if callable_call else np.nan,
                "status": status,
            })

    calls = pd.DataFrame(call_rows)

    if calls.empty:
        raise RuntimeError(
            "No calls were recovered for comparative-selected SNPs."
        )

    calls.to_csv(CALL_AUDIT, sep="\t", index=False)

    site_rows = []

    site_keys = [
        "final_class", "hypothesis", "module",
        "chrom", "pos", "ref", "alt",
        "best_macro_accuracy", "score_margin", "comparative_BA",
    ]

    for key, group in calls.groupby(site_keys, dropna=False):
        base = dict(zip(site_keys, key))

        callable_rows = group[group["callable"]]
        discordant = callable_rows[
            callable_rows["status"].eq(
                "callable_genotype_discordance"
            )
        ]
        n_or_qc = group[~group["callable"]]

        base.update({
            "expected_individuals": len(group),
            "callable_individuals": len(callable_rows),
            "correct_individuals": int(
                callable_rows["correct"].fillna(False).sum()
            ),
            "callable_discordances": len(discordant),
            "N_or_QC": len(n_or_qc),
            "callable_accuracy": safe_fraction(
                int(callable_rows["correct"].fillna(False).sum()),
                len(callable_rows),
            ),
            "callable_perfect": len(discordant) == 0,
            "discordant_samples": joined_unique(
                discordant["sample"]
            ),
            "N_or_QC_samples": joined_unique(
                n_or_qc["sample"]
            ),
        })

        for diagnostic_group in EXPECTED_GROUP_SIZES:
            subset = group[
                group["diagnostic_group"].eq(diagnostic_group)
            ]
            subset_callable = subset[subset["callable"]]
            subset_discordant = subset_callable[
                subset_callable["status"].eq(
                    "callable_genotype_discordance"
                )
            ]

            base[
                f"{diagnostic_group}_callable"
            ] = len(subset_callable)

            base[
                f"{diagnostic_group}_discordant"
            ] = len(subset_discordant)

            base[
                f"{diagnostic_group}_N_or_QC"
            ] = int((~subset["callable"]).sum())

            base[
                f"{diagnostic_group}_discordant_samples"
            ] = joined_unique(subset_discordant["sample"])

        site_rows.append(base)

    site_summary = pd.DataFrame(site_rows).sort_values(
        ["hypothesis", "pos"]
    )

    site_summary.to_csv(SITE_SUMMARY, sep="\t", index=False)

    individual_rows = []

    for key, group in calls.groupby(
        [
            "final_class", "hypothesis", "sample",
            "diagnostic_group", "taxon", "country",
        ],
        dropna=False,
    ):
        (
            final_class,
            hypothesis,
            sample,
            diagnostic_group,
            taxon,
            country,
        ) = key

        callable_rows = group[group["callable"]]
        discordant = callable_rows[
            callable_rows["status"].eq(
                "callable_genotype_discordance"
            )
        ]

        individual_rows.append({
            "final_class": final_class,
            "hypothesis": hypothesis,
            "sample": sample,
            "diagnostic_group": diagnostic_group,
            "taxon": taxon,
            "country": country,
            "selected_sites": len(group),
            "callable_sites": len(callable_rows),
            "correct_sites": int(
                callable_rows["correct"].fillna(False).sum()
            ),
            "discordant_sites": len(discordant),
            "N_or_QC_sites": int((~group["callable"]).sum()),
            "callable_discordance_fraction": safe_fraction(
                len(discordant),
                len(callable_rows),
            ),
            "discordant_positions": joined_unique(
                discordant["pos"]
            ),
        })

    individual_summary = pd.DataFrame(individual_rows).sort_values(
        [
            "final_class",
            "hypothesis",
            "discordant_sites",
            "sample",
        ],
        ascending=[True, True, False, True],
    )

    individual_summary.to_csv(
        INDIVIDUAL_SUMMARY,
        sep="\t",
        index=False,
    )

    population_rows = []

    for key, group in calls.groupby(
        [
            "final_class", "hypothesis",
            "diagnostic_group", "taxon", "country",
        ],
        dropna=False,
    ):
        (
            final_class,
            hypothesis,
            diagnostic_group,
            taxon,
            country,
        ) = key

        callable_rows = group[group["callable"]]
        discordant = callable_rows[
            callable_rows["status"].eq(
                "callable_genotype_discordance"
            )
        ]

        population_rows.append({
            "final_class": final_class,
            "hypothesis": hypothesis,
            "diagnostic_group": diagnostic_group,
            "taxon": taxon,
            "country": country,
            "individuals": group["sample"].nunique(),
            "selected_site_calls": len(group),
            "callable_site_calls": len(callable_rows),
            "discordant_site_calls": len(discordant),
            "N_or_QC_site_calls": int((~group["callable"]).sum()),
            "callable_discordance_fraction": safe_fraction(
                len(discordant),
                len(callable_rows),
            ),
            "individuals_with_discordance": (
                discordant["sample"].nunique()
            ),
            "discordant_samples": joined_unique(
                discordant["sample"]
            ),
        })

    population_summary = pd.DataFrame(population_rows).sort_values(
        [
            "final_class",
            "hypothesis",
            "callable_discordance_fraction",
        ],
        ascending=[True, True, False],
    )

    population_summary.to_csv(
        POPULATION_SUMMARY,
        sep="\t",
        index=False,
    )

    group_rows = []

    for key, group in calls.groupby(
        ["final_class", "hypothesis", "diagnostic_group"],
        dropna=False,
    ):
        final_class, hypothesis, diagnostic_group = key
        callable_rows = group[group["callable"]]
        discordant = callable_rows[
            callable_rows["status"].eq(
                "callable_genotype_discordance"
            )
        ]

        group_rows.append({
            "final_class": final_class,
            "hypothesis": hypothesis,
            "diagnostic_group": diagnostic_group,
            "individuals": group["sample"].nunique(),
            "selected_sites": group["pos"].nunique(),
            "callable_site_calls": len(callable_rows),
            "discordant_site_calls": len(discordant),
            "N_or_QC_site_calls": int((~group["callable"]).sum()),
            "callable_discordance_fraction": safe_fraction(
                len(discordant),
                len(callable_rows),
            ),
            "individuals_with_discordance": (
                discordant["sample"].nunique()
            ),
        })

    group_summary = pd.DataFrame(group_rows).sort_values(
        ["hypothesis", "diagnostic_group"]
    )

    group_summary.to_csv(
        GROUP_SUMMARY,
        sep="\t",
        index=False,
    )

    hypothesis_rows = []

    for hypothesis, group in site_summary.groupby("hypothesis"):
        hypothesis_rows.append({
            "final_class": FINAL_CLASS[hypothesis],
            "hypothesis": hypothesis,
            "comparative_selected_SNPs": len(group),
            "callable_perfect_SNPs": int(
                group["callable_perfect"].sum()
            ),
            "SNPs_with_callable_discordance": int(
                (~group["callable_perfect"]).sum()
            ),
            "total_callable_discordances": int(
                group["callable_discordances"].sum()
            ),
            "total_N_or_QC_calls": int(group["N_or_QC"].sum()),
            "minimum_macro_accuracy": (
                group["best_macro_accuracy"].min()
            ),
            "median_macro_accuracy": (
                group["best_macro_accuracy"].median()
            ),
            "maximum_macro_accuracy": (
                group["best_macro_accuracy"].max()
            ),
        })

    hypothesis_summary = pd.DataFrame(hypothesis_rows)
    hypothesis_summary.to_csv(
        HYPOTHESIS_SUMMARY,
        sep="\t",
        index=False,
    )

    print("COMPARATIVE-ONLY B / D / R SCORE-BREAKER AUDIT")
    print()
    print(
        "SNPs were selected and classified using comparative "
        "butterflies only."
    )
    print(
        "TN_Left and TN_Right were not read or used in this audit."
    )
    print(
        "Callable genotype discordance reduces accuracy; "
        "N/low-DP/low-GQ reduces callability."
    )

    print()
    print("HYPOTHESIS SUMMARY")
    print(
        hypothesis_summary.to_csv(
            sep="\t",
            index=False,
        ).rstrip()
    )

    print()
    print("DIAGNOSTIC-GROUP SUMMARY")
    print(
        group_summary.to_csv(
            sep="\t",
            index=False,
        ).rstrip()
    )

    print()
    print("ALL FIVE COMPARATIVE-SELECTED BAND SNPS")
    band_sites = site_summary[
        site_summary["hypothesis"].eq("B")
    ][
        [
            "chrom", "pos",
            "best_macro_accuracy",
            "callable_discordances",
            "N_or_QC",
            "besckei_B_discordant_samples",
            "meriana_D_discordant_samples",
            "contigua_R_discordant_samples",
            "N_or_QC_samples",
        ]
    ]

    print(
        band_sites.to_csv(
            sep="\t",
            index=False,
        ).rstrip()
    )

    print()
    print("INDIVIDUALS REPEATEDLY BREAKING COMPARATIVE SNP PATTERNS")
    repeated = individual_summary[
        individual_summary["discordant_sites"] > 0
    ][
        [
            "final_class", "hypothesis", "sample",
            "diagnostic_group", "taxon", "country",
            "selected_sites", "callable_sites",
            "discordant_sites", "N_or_QC_sites",
            "callable_discordance_fraction",
        ]
    ].sort_values(
        [
            "final_class",
            "hypothesis",
            "discordant_sites",
            "callable_discordance_fraction",
        ],
        ascending=[True, True, False, False],
    )

    print(
        repeated.to_csv(
            sep="\t",
            index=False,
        ).rstrip()
    )

    print()
    print("Outputs:")
    for path in [
        CALL_AUDIT,
        SITE_SUMMARY,
        INDIVIDUAL_SUMMARY,
        POPULATION_SUMMARY,
        GROUP_SUMMARY,
        HYPOTHESIS_SUMMARY,
    ]:
        print(path)


if __name__ == "__main__":
    main()
