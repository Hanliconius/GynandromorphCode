#!/usr/bin/env python3
"""Geographic sensitivity analysis; reuse existing BA code and figure layout."""
import argparse
import hashlib
import importlib.util
import json
import re
import sys
import unicodedata
from pathlib import Path

import numpy as np
import pandas as pd
import scipy

NEW = "B_western_geographic_reanalysis"
OLD = "B_heurippa_besckei_no_DR_positive_reanalysis"
KEY = ["chrom", "pos", "ref", "alt"]
COHORTS = ("WGS_only", "WGS_plus_CS")

SOURCES = {
    "Arias2014": (
        "https://repository.si.edu/bitstreams/"
        "3e39b5e3-03f5-4684-865f-ce904c7bc8f2/download"
    ),
    "Hausdorf2025": (
        "https://onlinelibrary.wiley.com/doi/full/"
        "10.1111%2Fmec.17717"
    ),
    "Salazar2010": (
        "https://journals.plos.org/plosgenetics/"
        "article?id=10.1371/journal.pgen.1000930"
    ),
}

# Strict Pacific-west plus Central America.
# Country consistency is checked. Missing country is accepted only
# for geographically restricted named taxa.
RULES = {
    "melpomene cythera": (
        "Pacific_west",
        ["ecuador"],
        "Hausdorf2025",
    ),
    "melpomene vulcanus": (
        "Pacific_west_or_Panama",
        ["colombia", "panama"],
        "Hausdorf2025",
    ),
    "melpomene rosina": (
        "Central_America",
        ["panama", "costa rica"],
        "Hausdorf2025",
    ),
    "cydno chioneus": (
        "Central_America",
        ["panama", "costa rica"],
        "Arias2014",
    ),
    "cydno zelinde": (
        "Pacific_west",
        ["colombia"],
        "Arias2014",
    ),
    "cydno alithea": (
        "Pacific_west",
        ["ecuador"],
        "Arias2014",
    ),
    "cydno galanthus": (
        "Central_America",
        [
            "mexico",
            "belize",
            "guatemala",
            "honduras",
            "el salvador",
            "nicaragua",
            "costa rica",
        ],
        "Arias2014",
    ),
    "pachinus": (
        "Central_America",
        ["costa rica", "panama"],
        "Arias2014",
    ),
}


def clean(value):
    return re.sub(r"\s+", " ", str(value)).strip()


def taxon_key(value):
    value = (
        clean(value)
        .lower()
        .replace("'", "")
        .replace('"', "")
    )

    return re.sub(
        r"\b(?:heliconius|ssp|spp|subsp)\.?\s*",
        "",
        value,
    ).strip()


def country_key(value):
    value = (
        unicodedata.normalize("NFKD", clean(value))
        .encode("ascii", "ignore")
        .decode()
        .lower()
    )

    if value.startswith("panam"):
        return "panama"

    if value in (
        "",
        "na",
        "nan",
        "none",
        "unknown",
    ):
        return ""

    return value


def geography(taxon, country):
    key = taxon_key(taxon)
    country = country_key(country)

    if key in RULES:
        region, countries, source = RULES[key]

        if country and country not in countries:
            return (
                False,
                "taxon_country_conflict",
                region,
                source,
            )

        decision = (
            "western_taxon"
            if country
            else "western_taxon_country_missing"
        )

        return True, decision, region, source

    if (
        key == "melpomene melpomene"
        and country == "panama"
    ):
        return (
            True,
            "nominal_melpomene_verified_Panama",
            "Central_America",
            "Hausdorf2025",
        )

    if key.startswith(
        (
            "cydno weymeri",
            "cydno cydnides",
        )
    ):
        return (
            False,
            "inter_Andean_Cauca_not_strict_Pacific_west",
            "Cauca_valley",
            "Arias2014",
        )

    if key.startswith("heurippa"):
        return (
            False,
            "eastern_Andes",
            "eastern_Andes",
            "Salazar2010",
        )

    if key.startswith("timareta"):
        return (
            False,
            "eastern_Andes",
            "eastern_Andes",
            "Arias2014",
        )

    return (
        False,
        "not_in_supported_western_inclusion_list",
        "unassigned_in_this_audit",
        "",
    )


def load(path, name):
    specification = importlib.util.spec_from_file_location(
        name,
        path,
    )

    module = importlib.util.module_from_spec(
        specification
    )

    sys.modules[name] = module
    specification.loader.exec_module(module)

    return module


def digest(path):
    result = hashlib.sha256()

    with Path(path).open("rb") as handle:
        for block in iter(
            lambda: handle.read(8388608),
            b"",
        ):
            result.update(block)

    return result.hexdigest()


def save(frame, path):
    frame.to_csv(
        path,
        sep="\t",
        index=False,
        na_rep="NA",
    )


def prepare_plot(source, target):
    text = source.read_text()

    changes = {
        OLD:
            NEW,

        "broad_B_heurippa_besckei_no_DR_positive_BA_manhattan":
            "broad_B_western_BA_manhattan",

        "Band association after excluding D-positive and R-positive samples":
            "Band association in Pacific-west / Central American samples",

        "Pooled genotypes at exact expected left-right sites":
            "Mixture-adjusted genotypes at exact expected left-right sites",

        "SNP selection and balanced accuracy use comparative samples only.":
            "Pair posterior >=0.95; BA is calculated independently of focal calls.",

        "All other genotype combinations and QC failures are omitted.":
            "Other pairs and uncertain adjusted calls are omitted; no original-call fallback.",
    }

    for before, after in changes.items():
        if before not in text:
            raise RuntimeError(
                "Plot source has changed; missing: "
                + before
            )

        text = text.replace(
            before,
            after,
        )

    # Replace the input route without changing figure geometry.
    start = text.index(
        "score_file <- file.path("
    )

    end = text.index(
        "if (!file.exists(score_file))",
        start,
    )

    text = (
        text[:start]
        + (
            'score_file <- file.path('
            'checkpoint_dir, '
            '"B_western_BA_and_adjusted_LR.tsv"'
            ')\n\n'
        )
        + text[end:]
    )

    # The adjusted pair posterior replaces original-GQ filtering.
    start = text.index(
        "d$left_QC <- ("
    )

    end = text.index(
        "d$same_GT <- (",
        start,
    )

    adjusted_qc = '''
if (!("mixture_pair_posterior" %in% names(d))) {
    stop("Missing adjusted-call posterior")
}

d$adjusted_confident <- (
    is.finite(d$mixture_pair_posterior)
    & d$mixture_pair_posterior >= 0.95
)

d$left_QC <- (
    d$adjusted_confident
    & !is.na(d$left_GT)
)

d$right_QC <- (
    d$adjusted_confident
    & !is.na(d$right_GT)
)

d$both_QC <- (
    d$left_QC
    & d$right_QC
)

'''

    text = (
        text[:start]
        + adjusted_qc
        + text[end:]
    )

    # Allow the script to render even if no exact pair exists.
    text = text.replace(
        'side = "TN_Left",\n'
        '        y = 1,',
        'side = rep("TN_Left", nrow(difference_sites)),\n'
        '        y = rep(1, nrow(difference_sites)),',
    )

    text = text.replace(
        'side = "TN_Right",\n'
        '        y = 0,',
        'side = rep("TN_Right", nrow(difference_sites)),\n'
        '        y = rep(0, nrow(difference_sites)),',
    )

    target.write_text(text)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--prepare-only",
        action="store_true",
    )

    parser.add_argument(
        "--project",
        type=Path,
        default=Path(
            "/scratch/martinlab/joe/GynandroCalling"
        ),
    )

    args = parser.parse_args()
    project = args.project.resolve()
    scripts = project / "scripts"

    out = (
        project
        / "analysis"
        / NEW
        / "checkpoints"
    )

    audit = (
        out
        / "scorer_audit_original_fields"
    )

    meta_path = (
        project
        / "data"
        / "derived"
        / "region_panel_metadata.tsv"
    )

    vcf = (
        project
        / "external"
        / "krista_genotype_plots_original"
        / "Genotype_Plots"
        / "Optix_Hmr18-938195-1239421_biallelic.vcf"
    )

    mixture_path = (
        project
        / "analysis"
        / "original_502"
        / "reciprocal_L_R_mixture"
        / "optix_reciprocal_L_R_genotype_pairs.tsv.gz"
    )

    source_paths = [
        scripts
        / "41_score_original_502_optix_callsets.py",

        scripts
        / "110b_base_scorer_B_no_DR_positive.py",

        scripts
        / "110c_scan_B_no_DR_positive.py",

        scripts
        / "110d_plot_B_no_DR_positive_BA_and_focal.R",
    ]

    for path in (
        source_paths
        + [
            meta_path,
            vcf,
            mixture_path,
        ]
    ):
        if not path.is_file():
            raise RuntimeError(
                f"Missing existing input: {path}"
            )

    audit.mkdir(
        parents=True,
        exist_ok=True,
    )

    original = load(
        source_paths[0],
        "west140_original",
    )

    helper = load(
        source_paths[1],
        "west140_helper",
    )

    scanner = load(
        source_paths[2],
        "west140_scanner",
    )

    original.METADATA_FILE = meta_path
    helper.META = meta_path

    metadata = helper.load_metadata(
        original
    )

    source_metadata = pd.read_csv(
        meta_path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    source_metadata.columns = (
        source_metadata.columns.str.strip()
    )

    required_metadata = {
        "Accession",
        "Country",
        "Taxon Name",
    }

    if not required_metadata.issubset(
        source_metadata.columns
    ):
        raise RuntimeError(
            "Missing geographic metadata columns: "
            + str(
                required_metadata
                - set(source_metadata.columns)
            )
        )

    source_metadata["Accession"] = (
        source_metadata["Accession"].map(clean)
    )

    if source_metadata["Accession"].duplicated().any():
        raise RuntimeError(
            "Duplicate accessions in metadata; "
            "resolve before geographic selection"
        )

    locations = (
        source_metadata
        .set_index("Accession")
        .to_dict("index")
    )

    samples = None

    with vcf.open() as handle:
        for line in handle:
            if line.startswith("#CHROM\t"):
                samples = (
                    line.rstrip("\r\n")
                    .split("\t")[9:]
                )
                break

    if (
        samples is None
        or len(samples) != 502
        or len(set(samples)) != 502
    ):
        raise RuntimeError(
            "Expected 502 unique VCF samples"
        )

    if (
        set(samples) - set(metadata)
        or set(samples) - set(locations)
    ):
        raise RuntimeError(
            "VCF samples missing metadata"
        )

    manifest_rows = []

    for cohort in COHORTS:
        for sample in samples:
            metadata_row = metadata[sample]
            location = locations[sample]

            taxon = clean(
                metadata_row["taxon"]
            )

            country = clean(
                location["Country"]
            )

            (
                western,
                geographic_decision,
                geographic_region,
                geographic_source,
            ) = geography(
                taxon,
                country,
            )

            data_type = clean(
                metadata_row["data_type"]
            ).upper()

            reasons = (
                []
                if western
                else [geographic_decision]
            )

            if sample in original.FOCAL_SAMPLES:
                reasons.append("focal")

            if (
                metadata_row.get("piano_key")
                or re.search(
                    r"piano\s*keys?",
                    taxon,
                    re.I,
                )
            ):
                reasons.append("piano_key")

            if metadata_row["B"] not in (0, 1):
                reasons.append("missing_B")

            permitted_types = (
                {"WGS"}
                if cohort == "WGS_only"
                else {"WGS", "CS"}
            )

            if data_type not in permitted_types:
                reasons.append(
                    "sequencing_type"
                )

            manifest_rows.append(
                {
                    "cohort": cohort,
                    "module": "Band",
                    "sample": sample,
                    "taxon_name": taxon,
                    "country": country,
                    "latitude":
                        location.get("Latitude", ""),
                    "longitude":
                        location.get("Longitude", ""),
                    "data_type": data_type,
                    "B": metadata_row["B"],
                    "geographic_region":
                        geographic_region,
                    "geographic_decision":
                        geographic_decision,
                    "source":
                        SOURCES.get(
                            geographic_source,
                            "",
                        ),
                    "keep": not reasons,
                    "keep_bool": not reasons,
                    "reason": (
                        ";".join(reasons)
                        if reasons
                        else "included"
                    ),
                }
            )

    manifest = pd.DataFrame(
        manifest_rows
    )

    manifest_file = (
        out
        / "B_western_cohort_manifest.tsv"
    )

    save(
        manifest,
        manifest_file,
    )

    retained = manifest.loc[
        manifest["keep"]
    ].copy()

    cohort_counts = (
        retained
        .groupby(
            [
                "cohort",
                "B",
            ]
        )
        .size()
        .rename("samples")
        .reset_index()
    )

    save(
        cohort_counts,
        out / "B_western_cohort_counts.tsv",
    )

    taxon_counts = (
        retained
        .groupby(
            [
                "cohort",
                "taxon_name",
                "country",
                "B",
            ]
        )
        .size()
        .rename("samples")
        .reset_index()
    )

    save(
        taxon_counts,
        out / "B_western_retained_taxa.tsv",
    )

    print(
        "WESTERN BAND COHORT\n",
        cohort_counts.to_string(index=False),
        flush=True,
    )

    print(
        taxon_counts.to_string(index=False),
        flush=True,
    )

    print(
        "\nNo D/R filter, no focal selection, "
        "and no individual discordance exclusions.\n"
        "Pachinus included geographically; "
        "Cauca-valley taxa excluded from this "
        "strict definition.\n"
        "Original comparative callability and "
        "minimum class-count thresholds are unchanged.",
        flush=True,
    )

    for cohort in COHORTS:
        counts = (
            retained.loc[
                retained["cohort"] == cohort,
                "B",
            ]
            .value_counts()
        )

        if any(
            counts.get(state, 0)
            < original.MINIMUM_CALLED_PER_PHENOTYPE
            for state in (0, 1)
        ):
            raise RuntimeError(
                f"{cohort}: too few individuals for "
                "the existing per-class eligibility rule; "
                "inspect the manifest"
            )

    plotter = (
        scripts
        / "140_plot_B_western_adjusted_LR.R"
    )

    prepare_plot(
        source_paths[3],
        plotter,
    )

    if args.prepare_only:
        print(
            f"\nPrepared: {plotter}\n"
            f"Manifest: {manifest_file}",
            flush=True,
        )
        return

    provenance = {
        "version": "140.1",

        "definition":
            "Pacific-west plus Central America",

        "western_taxon_rules":
            RULES,

        "sources":
            SOURCES,

        "nominal_melpomene":
            "Panama only",

        "comparative_QC":
            (
                "Original called_GT and DP2_GQ20; "
                "original eligibility thresholds"
            ),

        "focal_calls":
            (
                "Reciprocal mixture only; "
                "pair posterior >=0.95; "
                "no raw-call fallback"
            ),

        "python":
            sys.version,

        "numpy":
            np.__version__,

        "pandas":
            pd.__version__,

        "scipy":
            scipy.__version__,

        "input_SHA256": {
            str(path): digest(path)
            for path in (
                source_paths
                + [
                    meta_path,
                    vcf,
                    mixture_path,
                    Path(__file__).resolve(),
                ]
            )
        },
    }

    signature = hashlib.sha256(
        json.dumps(
            provenance,
            sort_keys=True,
        ).encode()
    ).hexdigest()

    (
        out
        / "B_western_provenance.json"
    ).write_text(
        json.dumps(
            provenance,
            indent=2,
        )
    )

    completion_file = (
        out
        / "scoring_complete.json"
    )

    helper.load_original_module = (
        lambda: original
    )

    helper.load_metadata = (
        lambda scorer: metadata
    )

    helper.load_final_filtered_manifest = (
        lambda: manifest.copy()
    )

    scanner.load_helper = (
        lambda: helper
    )

    scanner.PROJECT = project
    scanner.OUTDIR = audit

    for name, value in list(
        vars(scanner).items()
    ):
        if (
            name.startswith("OUT_")
            and isinstance(value, Path)
        ):
            setattr(
                scanner,
                name,
                audit / value.name,
            )

    cached = (
        json.loads(
            completion_file.read_text()
        )
        if completion_file.exists()
        else {}
    )

    reuse = (
        cached.get("signature") == signature
        and scanner.OUT_RAW.is_file()
        and cached.get("raw_SHA256")
            == digest(scanner.OUT_RAW)
    )

    if not reuse:
        scanner.main()
    else:
        print(
            "Verified existing scores; "
            "rebuilding adjusted-call display input.",
            flush=True,
        )

    raw_scores = pd.read_csv(
        scanner.OUT_RAW,
        sep="\t",
        low_memory=False,
    )

    if (
        raw_scores
        .groupby(
            [
                "analysis_cohort",
                "scheme",
            ]
        )
        .ngroups
        != 4
    ):
        raise RuntimeError(
            "Expected two cohorts and two schemes"
        )

    for (
        cohort,
        scheme,
    ), group in raw_scores.groupby(
        [
            "analysis_cohort",
            "scheme",
        ]
    ):
        expected_samples = int(
            (
                retained["cohort"]
                == cohort
            ).sum()
        )

        if (
            len(group) != 57797
            or group.duplicated(KEY).any()
            or set(
                group["comparative_cohort_n"]
            ) != {expected_samples}
        ):
            raise RuntimeError(
                "Scoring validation failed: "
                f"{cohort}/{scheme}"
            )

    mixture = pd.read_csv(
        mixture_path,
        sep="\t",
        low_memory=False,
    )

    mixture_columns = (
        KEY
        + [
            "mixture_left_GT",
            "mixture_right_GT",
            "pair_posterior",
        ]
    )

    if (
        not set(mixture_columns).issubset(
            mixture.columns
        )
        or mixture.duplicated(KEY).any()
    ):
        raise RuntimeError(
            "Mixture table schema or "
            "site uniqueness failure"
        )

    combined = raw_scores.merge(
        mixture[mixture_columns],
        on=KEY,
        how="left",
        validate="many_to_one",
        sort=False,
    )

    if (
        len(combined) != len(raw_scores)
        or not np.allclose(
            combined["balanced_accuracy"],
            raw_scores["balanced_accuracy"],
            equal_nan=True,
        )
    ):
        raise RuntimeError(
            "Focal join changed comparative BA "
            "or row count"
        )

    confident = (
        pd.to_numeric(
            combined["pair_posterior"],
            errors="coerce",
        )
        .ge(0.95)
    )

    valid_genotypes = {
        "0/0",
        "0/1",
        "1/1",
    }

    if confident.any():
        valid_left = (
            combined.loc[
                confident,
                "mixture_left_GT",
            ]
            .isin(valid_genotypes)
            .all()
        )

        valid_right = (
            combined.loc[
                confident,
                "mixture_right_GT",
            ]
            .isin(valid_genotypes)
            .all()
        )

        if not (
            valid_left
            and valid_right
        ):
            raise RuntimeError(
                "Invalid confident mixture genotype"
            )

    if not confident.any():
        raise RuntimeError(
            "No confident adjusted calls joined; "
            "no original-call fallback"
        )

    for side in (
        "Left",
        "Right",
    ):
        combined[
            f"TN_{side}_GT"
        ] = combined[
            f"mixture_{side.lower()}_GT"
        ].where(
            confident,
            "./.",
        )

    combined[
        "mixture_pair_posterior"
    ] = combined["pair_posterior"]

    # These flags were calculated from raw focal calls.
    # Never carry them into the adjusted-call table.
    combined = combined.drop(
        columns=[
            "pooled_both_DP2_GQ20",
            "pooled_GT_different",
            "left_11_right_01",
            "dominant_focal_model_match",
            "mixture_left_GT",
            "mixture_right_GT",
            "pair_posterior",
        ]
    )

    combined[
        "focal_call_source"
    ] = (
        "reciprocal_mixture_posterior95_only"
    )

    combined[
        "focal_pair_confident"
    ] = confident

    combined[
        "focal_exact_B_pair"
    ] = (
        confident
        & combined["TN_Left_GT"].eq("1/1")
        & combined["TN_Right_GT"].eq("0/1")
    )

    plot_input = (
        out
        / "B_western_BA_and_adjusted_LR.tsv"
    )

    save(
        combined,
        plot_input,
    )

    primary = combined.loc[
        (
            combined["analysis_cohort"]
            == "WGS_plus_CS"
        )
        & (
            combined["scheme"]
            == "DP2_GQ20"
        )
    ].copy()

    primary["eligible"] = (
        primary["score_eligible"]
        .astype(str)
        .str.lower()
        .eq("true")
    )

    regional = primary.loc[
        primary["pos"].between(
            1112651,
            1114314,
        )
    ].copy()

    regional = regional.sort_values(
        [
            "eligible",
            "balanced_accuracy",
            "pos",
        ],
        ascending=[
            False,
            False,
            True,
        ],
    )

    save(
        regional,
        out
        / "Krista_band_western_scores_and_adjusted_LR.tsv",
    )

    display_columns = [
        "pos",
        "selected_orientation",
        "balanced_accuracy",
        "misclassified_n",
        "phenotype_0_called",
        "phenotype_1_called",
        "eligible",
        "TN_Left_GT",
        "TN_Right_GT",
        "focal_pair_confident",
    ]

    print(
        "\nORIGINAL BAND REGION: "
        "WESTERN COHORT, DP2/GQ20",
        flush=True,
    )

    print(
        regional[
            display_columns
        ]
        .head(15)
        .to_string(index=False),
        flush=True,
    )

    summary_rows = []

    regions = (
        (
            "full_VCF",
            primary,
        ),
        (
            "display_1060k_1240k",
            primary.loc[
                primary["pos"].between(
                    1060000,
                    1240000,
                )
            ],
        ),
        (
            "original_Band",
            regional,
        ),
    )

    for region_name, region_data in regions:
        summary_rows.append(
            {
                "region": region_name,
                "SNPs": len(region_data),

                "eligible":
                    int(
                        region_data[
                            "eligible"
                        ].sum()
                    ),

                "eligible_BA90":
                    int(
                        (
                            region_data["eligible"]
                            & region_data[
                                "balanced_accuracy"
                            ].ge(0.90)
                        ).sum()
                    ),

                "eligible_BA1":
                    int(
                        (
                            region_data["eligible"]
                            & np.isclose(
                                region_data[
                                    "balanced_accuracy"
                                ],
                                1,
                            )
                        ).sum()
                    ),

                "adjusted_confident":
                    int(
                        region_data[
                            "focal_pair_confident"
                        ].sum()
                    ),

                "adjusted_exact_B":
                    int(
                        region_data[
                            "focal_exact_B_pair"
                        ].sum()
                    ),

                "adjusted_exact_and_BA90":
                    int(
                        (
                            region_data[
                                "focal_exact_B_pair"
                            ]
                            & region_data["eligible"]
                            & region_data[
                                "balanced_accuracy"
                            ].ge(0.90)
                        ).sum()
                    ),
            }
        )

    summary = pd.DataFrame(
        summary_rows
    )

    save(
        summary,
        out / "B_western_summary.tsv",
    )

    print(
        "\n",
        summary.to_string(index=False),
        flush=True,
    )

    completion_file.write_text(
        json.dumps(
            {
                "signature":
                    signature,

                "raw_SHA256":
                    digest(scanner.OUT_RAW),
            },
            indent=2,
        )
    )

    checksum_rows = []

    for path in (
        plot_input,
        plotter,
        manifest_file,
    ):
        checksum_rows.append(
            {
                "file": str(path),
                "SHA256": digest(path),
            }
        )

    save(
        pd.DataFrame(checksum_rows),
        out / "B_western_output_checksums.tsv",
    )

    print(
        f"\nNew plot input: {plot_input}\n"
        "Original focal fields kept only in: "
        f"{audit}\n"
        f"Plotter: {plotter}",
        flush=True,
    )


if __name__ == "__main__":
    main()
