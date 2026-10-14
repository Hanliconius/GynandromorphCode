#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

CHECKPOINTS = (
    PROJECT
    / "analysis/drop_problematic_reanalysis/checkpoints"
)

OLD_COHORT_MANIFEST = (
    PROJECT
    / "analysis/original_502/checkpoints"
    / "original502_filtered_cohort_manifest.WGS_vs_WGSplusCS.tsv"
)

FILTERED_METADATA = (
    PROJECT
    / "data/derived"
    / "region_panel_metadata.drop_problematic.tsv"
)

VCF = (
    PROJECT
    / "external/krista_genotype_plots_original/Genotype_Plots"
    / "Optix_Hmr18-938195-1239421_biallelic.vcf"
)

OUTPUT = (
    CHECKPOINTS
    / "best_B_D_R_WGSplusCS_site_individual_breakers.tsv"
)

MODULES = {
    "Band": {
        "trait": "B",
        "score_file":
            "longBand_original502_expected_pattern_sites."
            "WGS_vs_WGSplusCS.tsv",
    },
    "Dennis": {
        "trait": "D",
        "score_file":
            "longDennis_original502_expected_pattern_sites."
            "WGS_vs_WGSplusCS.tsv",
    },
    "Ray": {
        "trait": "R",
        "score_file":
            "longRay_original502_expected_pattern_sites."
            "WGS_vs_WGSplusCS.tsv",
    },
}

EXPECTED_COHORT_SIZES = {
    "Band": 364,
    "Dennis": 373,
    "Ray": 372,
}


def truthy(value):
    return str(value).strip().lower() in {
        "true", "1", "yes", "t"
    }


def clean(value):
    return " ".join(str(value).replace("\n", " ").split())


def number(value):
    return pd.to_numeric(
        pd.Series([value]),
        errors="coerce",
    ).iloc[0]


def main():
    metadata = pd.read_csv(
        FILTERED_METADATA,
        sep="\t",
        low_memory=False,
    )

    metadata.columns = [
        clean(column)
        for column in metadata.columns
    ]

    metadata = metadata.set_index(
        metadata["Accession"].astype(str),
        drop=False,
    )

    cohort_manifest = pd.read_csv(
        OLD_COHORT_MANIFEST,
        sep="\t",
        low_memory=False,
    )

    retained_accessions = set(metadata.index.astype(str))

    top_sites = {}
    cohort_samples = {}

    for module, specification in MODULES.items():
        score_path = (
            CHECKPOINTS / specification["score_file"]
        )

        scores = pd.read_csv(
            score_path,
            sep="\t",
            low_memory=False,
        )

        subset = scores[
            scores["analysis_cohort"].astype(str).eq(
                "WGS_plus_CS"
            )
            & scores["scheme"].astype(str).eq("DP2_GQ20")
        ].copy()

        if "expected_pattern_eligible" in subset.columns:
            subset = subset[
                subset["expected_pattern_eligible"].map(truthy)
            ]

        subset["balanced_accuracy"] = pd.to_numeric(
            subset["balanced_accuracy"],
            errors="coerce",
        )

        subset = subset.sort_values(
            [
                "balanced_accuracy",
                "call_rate",
                "pos",
            ],
            ascending=[False, False, True],
        )

        if subset.empty:
            raise RuntimeError(
                f"No eligible WGS+CS scores for {module}"
            )

        top = subset.iloc[0]
        key = (
            str(top["chrom"]),
            int(top["pos"]),
            str(top["ref"]),
            str(top["alt"]),
        )

        top_sites[module] = {
            "key": key,
            "balanced_accuracy": top["balanced_accuracy"],
            "orientation": str(top["selected_orientation"]),
            "reported_errors": int(top["misclassified_n"]),
        }

        manifest_subset = cohort_manifest[
            cohort_manifest["cohort"].astype(str).eq(
                "WGS_plus_CS"
            )
            & cohort_manifest["module"].astype(str).eq(module)
            & cohort_manifest["keep"].map(truthy)
            & cohort_manifest["sample"]
                .astype(str)
                .isin(retained_accessions)
        ]

        samples = sorted(
            manifest_subset["sample"]
            .dropna()
            .astype(str)
            .unique()
        )

        expected_n = int(
            number(top["comparative_cohort_n"])
        )

        if len(samples) != expected_n:
            raise RuntimeError(
                f"{module}: recovered {len(samples)} cohort samples; "
                f"expected {expected_n}"
            )

        cohort_samples[module] = samples

    target_keys = {
        values["key"]
        for values in top_sites.values()
    }

    calls = {}
    sample_order = []

    with VCF.open() as handle:
        for line in handle:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                sample_order = header[9:]
                continue

            prefix = line.split("\t", 5)

            if len(prefix) < 5:
                continue

            chrom = prefix[0]
            pos = int(prefix[1])

            possible = [
                key for key in target_keys
                if key[0] == chrom and key[1] == pos
            ]

            if not possible:
                continue

            fields = line.rstrip("\n").split("\t")
            key = (
                fields[0],
                int(fields[1]),
                fields[3],
                fields[4],
            )

            if key not in target_keys:
                continue

            format_names = fields[8].split(":")
            format_index = {
                name: index
                for index, name in enumerate(format_names)
            }

            site_calls = {}

            for sample, value in zip(
                sample_order,
                fields[9:],
            ):
                components = value.split(":")

                def get_format(name):
                    index = format_index.get(name)
                    if index is None or index >= len(components):
                        return None
                    return components[index]

                site_calls[sample] = {
                    "GT": get_format("GT"),
                    "DP": get_format("DP"),
                    "GQ": get_format("GQ"),
                    "AD": get_format("AD"),
                }

            calls[key] = site_calls

    missing_sites = [
        key for key in target_keys
        if key not in calls
    ]

    if missing_sites:
        raise RuntimeError(
            f"Top sites not recovered from VCF: {missing_sites}"
        )

    output_rows = []

    for module, site in top_sites.items():
        trait = MODULES[module]["trait"]
        key = site["key"]
        orientation = site["orientation"]

        for sample in cohort_samples[module]:
            call = calls[key].get(sample, {})

            GT = str(call.get("GT"))
            DP = number(call.get("DP"))
            GQ = number(call.get("GQ"))

            if (
                GT in {"None", ".", "./.", ".|."}
                or pd.isna(DP)
                or pd.isna(GQ)
                or DP < 2
                or GQ < 20
            ):
                continue

            alleles = (
                GT.replace("|", "/").split("/")
            )

            if orientation == "REF_carrier_predicts_presence":
                predicted = int("0" in alleles)
            elif orientation == "ALT_carrier_predicts_presence":
                predicted = int("1" in alleles)
            else:
                raise RuntimeError(
                    f"Unexpected orientation: {orientation}"
                )

            phenotype = number(
                metadata.loc[sample, trait]
            )

            taxon = clean(
                metadata.loc[sample, "Taxon Name"]
            )

            # Preserve the established Band phenotype correction.
            if (
                trait == "B"
                and (
                    "malleti" in taxon.lower()
                    or "ecuadoriensis" in taxon.lower()
                )
            ):
                phenotype = 0

            if pd.isna(phenotype):
                continue

            phenotype = int(phenotype)

            if predicted != phenotype:
                output_rows.append({
                    "module": module,
                    "trait": trait,
                    "chrom": key[0],
                    "pos": key[1],
                    "ref": key[2],
                    "alt": key[3],
                    "balanced_accuracy":
                        site["balanced_accuracy"],
                    "reported_misclassified_n":
                        site["reported_errors"],
                    "sample": sample,
                    "taxon": taxon,
                    "country": clean(
                        metadata.loc[sample, "Country"]
                    ),
                    "phenotype": phenotype,
                    "predicted_phenotype": predicted,
                    "GT": GT,
                    "AD": call.get("AD"),
                    "DP": DP,
                    "GQ": GQ,
                    "orientation": orientation,
                })

    breakers = pd.DataFrame(output_rows)

    breakers.to_csv(
        OUTPUT,
        sep="\t",
        index=False,
    )

    print("BEST B / D / R WGS+CS SITE BREAKERS")
    print()

    for module, site in top_sites.items():
        subset = breakers[
            breakers["module"].eq(module)
        ]

        print(
            f'{module}\t'
            f'{site["key"][0]}:{site["key"][1]}\t'
            f'BA={site["balanced_accuracy"]:.6f}\t'
            f'reported_errors={site["reported_errors"]}\t'
            f'recovered_errors={len(subset)}'
        )

        if subset.empty:
            print("  none")
        else:
            for _, row in subset.iterrows():
                print(
                    f'  {row["sample"]}\t'
                    f'{row["taxon"]}\t'
                    f'phenotype={row["phenotype"]}\t'
                    f'predicted={row["predicted_phenotype"]}\t'
                    f'GT={row["GT"]}\t'
                    f'DP={row["DP"]:.0f}\t'
                    f'GQ={row["GQ"]:.0f}'
                )

    print()
    print(f"Output: {OUTPUT}")


if __name__ == "__main__":
    main()
