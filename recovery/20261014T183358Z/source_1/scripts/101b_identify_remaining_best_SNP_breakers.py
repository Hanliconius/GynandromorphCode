#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

CHECKPOINTS = (
    PROJECT
    / "analysis"
    / "drop_problematic_reanalysis"
    / "checkpoints"
)

VCF = (
    PROJECT
    / "external"
    / "krista_genotype_plots_original"
    / "Genotype_Plots"
    / "Optix_Hmr18-938195-1239421_biallelic.vcf"
)

MANIFEST = (
    CHECKPOINTS
    / "original502_filtered_cohort_manifest.WGS_vs_WGSplusCS.tsv"
)

METADATA = (
    PROJECT
    / "data"
    / "derived"
    / "region_panel_metadata.tsv"
)

OUT_BREAKERS = (
    CHECKPOINTS
    / "best_B_D_R_SNP_remaining_breakers.tsv"
)

OUT_ALL_CALLS = (
    CHECKPOINTS
    / "best_B_D_R_SNP_all_comparative_calls.tsv"
)

OUT_SITES = (
    CHECKPOINTS
    / "best_B_D_R_SNP_audit_summary.tsv"
)


MODULES = {
    "Band": {
        "prefix": "longBand",
        "manifest_module": "Band",
        "trait": "B",
    },
    "Dennis": {
        "prefix": "longDennis",
        "manifest_module": "Dennis2",
        "trait": "D",
    },
    "Ray": {
        "prefix": "longRay",
        "manifest_module": "Ray2",
        "trait": "R",
    },
}


def clean(value):

    if value is None:
        return ""

    return str(value).strip()


def as_bool(value):

    return clean(value).lower() in {
        "true",
        "t",
        "1",
        "yes",
        "y",
    }


def canonical_gt(value):

    value = clean(value)

    if value in {
        "",
        ".",
        "./.",
        ".|.",
    }:
        return None

    value = value.replace("|", "/")

    alleles = value.split("/")

    if (
        len(alleles) != 2
        or "." in alleles
    ):
        return None

    try:
        alleles = sorted(
            int(allele)
            for allele in alleles
        )
    except ValueError:
        return None

    return "/".join(
        str(allele)
        for allele in alleles
    )


def numeric(value):

    value = clean(value)

    if value in {
        "",
        ".",
    }:
        return None

    try:
        return float(value)
    except ValueError:
        return None


def find_column(columns, target):

    normalized_target = (
        target.lower()
        .replace(" ", "")
        .replace("_", "")
        .replace("-", "")
    )

    for column in columns:

        normalized_column = (
            clean(column)
            .lower()
            .replace(" ", "")
            .replace("_", "")
            .replace("-", "")
        )

        if normalized_column == normalized_target:
            return column

    raise RuntimeError(
        f"Could not find column corresponding to {target}. "
        f"Available columns: {list(columns)}"
    )


def load_best_sites():

    best_sites = []

    for module, specification in MODULES.items():

        score_file = (
            CHECKPOINTS
            / (
                f'{specification["prefix"]}_original502_'
                "expected_pattern_sites.WGS_vs_WGSplusCS.tsv"
            )
        )

        if not score_file.exists():
            raise RuntimeError(
                f"Missing score table: {score_file}"
            )

        scores = pd.read_csv(
            score_file,
            sep="\t",
            low_memory=False,
        )

        scores["balanced_accuracy"] = pd.to_numeric(
            scores["balanced_accuracy"],
            errors="coerce",
        )

        scores["misclassified_n"] = pd.to_numeric(
            scores["misclassified_n"],
            errors="coerce",
        )

        subset = scores.loc[
            (
                scores["analysis_cohort"]
                == "WGS_plus_CS"
            )
            & (
                scores["scheme"]
                == "DP2_GQ20"
            )
        ].copy()

        if "expected_pattern_eligible" in subset.columns:
            subset = subset.loc[
                subset[
                    "expected_pattern_eligible"
                ].map(as_bool)
            ].copy()

        if subset.empty:
            raise RuntimeError(
                f"No eligible WGS+CS DP2/GQ20 sites "
                f"were found for {module}."
            )

        subset = subset.sort_values(
            [
                "balanced_accuracy",
                "misclassified_n",
                "call_rate",
                "pos",
            ],
            ascending=[
                False,
                True,
                False,
                True,
            ],
        )

        best = subset.iloc[0]

        best_sites.append(
            {
                "module": module,
                "manifest_module": (
                    specification["manifest_module"]
                ),
                "trait": specification["trait"],
                "chrom": clean(best["chrom"]),
                "pos": int(best["pos"]),
                "ref": clean(best["ref"]),
                "alt": clean(best["alt"]),
                "selected_orientation": clean(
                    best["selected_orientation"]
                ),
                "reported_BA": float(
                    best["balanced_accuracy"]
                ),
                "reported_misclassified_n": int(
                    best["misclassified_n"]
                ),
                "reported_phenotype_0_called": int(
                    best["phenotype_0_called"]
                ),
                "reported_phenotype_1_called": int(
                    best["phenotype_1_called"]
                ),
            }
        )

    return best_sites


def load_allowed_samples(best_sites):

    manifest = pd.read_csv(
        MANIFEST,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    sample_sets = {}

    for site in best_sites:

        subset = manifest.loc[
            (
                manifest["cohort"]
                == "WGS_plus_CS"
            )
            & (
                manifest["module"]
                == site["manifest_module"]
            )
            & manifest["keep"].map(as_bool)
        ].copy()

        sample_sets[site["module"]] = set(
            subset["sample"].map(clean)
        )

    return sample_sets


def load_metadata():

    metadata = pd.read_csv(
        METADATA,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    accession_column = find_column(
        metadata.columns,
        "Accession",
    )

    taxon_column = find_column(
        metadata.columns,
        "Taxon Name",
    )

    country_column = find_column(
        metadata.columns,
        "Country",
    )

    trait_columns = {
        trait: find_column(
            metadata.columns,
            trait,
        )
        for trait in (
            "B",
            "D",
            "R",
        )
    }

    records = {}

    for _, row in metadata.iterrows():

        accession = clean(
            row[accession_column]
        )

        if not accession:
            continue

        record = {
            "taxon": clean(
                row[taxon_column]
            ).replace("\n", " "),
            "country": clean(
                row[country_column]
            ).replace("\n", " "),
        }

        for trait, column in trait_columns.items():

            value = numeric(
                row[column]
            )

            record[trait] = (
                int(value)
                if value in {
                    0.0,
                    1.0,
                }
                else None
            )

        records[accession] = record

    return records


def extract_requested_vcf_records(best_sites):

    requested = {
        (
            site["chrom"],
            site["pos"],
            site["ref"],
            site["alt"],
        )
        for site in best_sites
    }

    sample_names = None
    recovered = {}

    with VCF.open(
        "rt",
        encoding="utf-8",
        errors="replace",
    ) as handle:

        for line in handle:

            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):

                fields = line.rstrip("\n").split("\t")
                sample_names = fields[9:]
                continue

            if sample_names is None:
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 10:
                continue

            chrom = fields[0]

            try:
                pos = int(fields[1])
            except ValueError:
                continue

            ref = fields[3]
            alt = fields[4]

            key = (
                chrom,
                pos,
                ref,
                alt,
            )

            if key not in requested:
                continue

            format_fields = fields[8].split(":")

            per_sample = {}

            for sample, sample_value in zip(
                sample_names,
                fields[9:],
            ):

                values = sample_value.split(":")

                record = {
                    field: (
                        values[index]
                        if index < len(values)
                        else ""
                    )
                    for index, field in enumerate(
                        format_fields
                    )
                }

                per_sample[sample] = record

            recovered[key] = per_sample

            if len(recovered) == len(requested):
                break

    missing = requested.difference(
        recovered
    )

    if missing:
        raise RuntimeError(
            "The following best sites were not recovered "
            f"from the VCF: {sorted(missing)}"
        )

    return recovered


def genotype_prediction(
    canonical_genotype,
    orientation,
):

    if canonical_genotype is None:
        return None

    alleles = {
        int(allele)
        for allele in canonical_genotype.split("/")
    }

    if orientation == "ALT_carrier_predicts_presence":
        return int(
            1 in alleles
        )

    if orientation == "REF_carrier_predicts_presence":
        return int(
            0 in alleles
        )

    raise RuntimeError(
        "Unrecognized genotype orientation: "
        f"{orientation}"
    )


def main():

    for required in (
        VCF,
        MANIFEST,
        METADATA,
    ):
        if not required.exists():
            raise RuntimeError(
                f"Missing required input: {required}"
            )

    best_sites = load_best_sites()
    sample_sets = load_allowed_samples(
        best_sites
    )
    metadata = load_metadata()
    vcf_records = extract_requested_vcf_records(
        best_sites
    )

    all_calls = []
    breaker_rows = []
    summary_rows = []

    for site in best_sites:

        key = (
            site["chrom"],
            site["pos"],
            site["ref"],
            site["alt"],
        )

        allowed = sample_sets[
            site["module"]
        ]

        callable_n = 0
        correct_n = 0
        discordant_n = 0
        N_or_QC_n = 0
        phenotype_0_called = 0
        phenotype_1_called = 0

        for sample in sorted(allowed):

            metadata_record = metadata.get(
                sample
            )

            if metadata_record is None:
                continue

            phenotype = metadata_record.get(
                site["trait"]
            )

            if phenotype not in {
                0,
                1,
            }:
                continue

            vcf_call = vcf_records[
                key
            ].get(
                sample,
                {},
            )

            gt = canonical_gt(
                vcf_call.get(
                    "GT",
                    "",
                )
            )

            dp = numeric(
                vcf_call.get(
                    "DP",
                    "",
                )
            )

            gq = numeric(
                vcf_call.get(
                    "GQ",
                    "",
                )
            )

            ad = clean(
                vcf_call.get(
                    "AD",
                    "",
                )
            )

            qc_pass = (
                gt is not None
                and dp is not None
                and dp >= 2
                and gq is not None
                and gq >= 20
            )

            prediction = (
                genotype_prediction(
                    gt,
                    site["selected_orientation"],
                )
                if qc_pass
                else None
            )

            discordant = (
                qc_pass
                and prediction != phenotype
            )

            if qc_pass:

                callable_n += 1

                if phenotype == 0:
                    phenotype_0_called += 1
                else:
                    phenotype_1_called += 1

                if discordant:
                    discordant_n += 1
                else:
                    correct_n += 1

            else:
                N_or_QC_n += 1

            row = {
                "module": site["module"],
                "trait": site["trait"],
                "chrom": site["chrom"],
                "pos": site["pos"],
                "ref": site["ref"],
                "alt": site["alt"],
                "orientation": (
                    site["selected_orientation"]
                ),
                "sample": sample,
                "taxon": metadata_record["taxon"],
                "country": metadata_record["country"],
                "phenotype": phenotype,
                "GT": gt if gt is not None else "",
                "AD": ad,
                "DP": dp if dp is not None else "",
                "GQ": gq if gq is not None else "",
                "DP2_GQ20_pass": qc_pass,
                "predicted_phenotype": (
                    prediction
                    if prediction is not None
                    else ""
                ),
                "callable_discordance": discordant,
            }

            all_calls.append(
                row
            )

            if discordant:
                breaker_rows.append(
                    row
                )

        summary_rows.append(
            {
                **site,
                "allowed_manifest_samples": len(
                    allowed
                ),
                "phenotyped_callable_n": callable_n,
                "phenotype_0_called_recalculated": (
                    phenotype_0_called
                ),
                "phenotype_1_called_recalculated": (
                    phenotype_1_called
                ),
                "correct_callable_n": correct_n,
                "callable_discordances": discordant_n,
                "N_or_QC_calls": N_or_QC_n,
                "reported_counts_reproduced": (
                    phenotype_0_called
                    == site[
                        "reported_phenotype_0_called"
                    ]
                    and phenotype_1_called
                    == site[
                        "reported_phenotype_1_called"
                    ]
                    and discordant_n
                    == site[
                        "reported_misclassified_n"
                    ]
                ),
            }
        )

    all_calls_table = pd.DataFrame(
        all_calls
    )

    breakers_table = pd.DataFrame(
        breaker_rows
    )

    summary_table = pd.DataFrame(
        summary_rows
    )

    all_calls_table.to_csv(
        OUT_ALL_CALLS,
        sep="\t",
        index=False,
    )

    breakers_table.to_csv(
        OUT_BREAKERS,
        sep="\t",
        index=False,
    )

    summary_table.to_csv(
        OUT_SITES,
        sep="\t",
        index=False,
    )

    print("BEST COMPARATIVE SNP BREAKER AUDIT")
    print()
    print("BEST SITES")

    columns = [
        "module",
        "trait",
        "chrom",
        "pos",
        "ref",
        "alt",
        "selected_orientation",
        "reported_BA",
        "reported_misclassified_n",
        "callable_discordances",
        "N_or_QC_calls",
        "reported_counts_reproduced",
    ]

    print(
        summary_table[
            columns
        ].to_string(
            index=False
        )
    )

    print()
    print("CALLABLE INDIVIDUALS BREAKING EACH BEST SNP")

    if breakers_table.empty:
        print("None")
    else:
        print(
            breakers_table[
                [
                    "module",
                    "chrom",
                    "pos",
                    "sample",
                    "taxon",
                    "country",
                    "phenotype",
                    "GT",
                    "AD",
                    "DP",
                    "GQ",
                    "predicted_phenotype",
                ]
            ].to_string(
                index=False
            )
        )

    print()
    print(f"Summary: {OUT_SITES}")
    print(f"Breakers: {OUT_BREAKERS}")
    print(f"All calls: {OUT_ALL_CALLS}")


if __name__ == "__main__":
    main()
