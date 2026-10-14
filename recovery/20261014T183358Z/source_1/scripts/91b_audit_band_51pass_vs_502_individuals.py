#!/usr/bin/env python3

from pathlib import Path
import csv
import re

import numpy as np
import pandas as pd


PROJECT = Path(
    "/scratch/martinlab/joe/GynandroCalling"
)

CK = PROJECT / "analysis" / "checkpoints"

CK502 = (
    PROJECT
    / "analysis"
    / "original_502"
    / "checkpoints"
)

# ------------------------------------------------------------
# SITE-SCORE SOURCES
# ------------------------------------------------------------

SITES51 = (
    CK
    / "optix_51sample_stable_candidate_sites.tsv"
)

SITES502 = (
    CK502
    / "optix_original_502_paired_site_evidence.tsv"
)

# ------------------------------------------------------------
# ACTUAL GENOTYPE SOURCES
# ------------------------------------------------------------

VCF51 = (
    PROJECT
    / "data"
    / "raw"
    / "Optix_Hmr18-938195-1239421.vcf"
)

VCF502 = (
    PROJECT
    / "external"
    / "krista_genotype_plots_original"
    / "Genotype_Plots"
    / "optix_band-1112651-1114314.vcf"
)

# ------------------------------------------------------------
# EXACT SCORING COHORT / PHENOTYPE SOURCE
# ------------------------------------------------------------

COHORT502 = (
    CK502
    / "optix_original_502_cohort_manifest.tsv"
)

COHORT502_SUMMARY = (
    CK502
    / "optix_original_502_scoring_cohorts.tsv"
)

EXCLUDED51 = (
    CK
    / "optix_51sample_excluded_comparative_samples.tsv"
)

# ------------------------------------------------------------
# METADATA SOURCES
# ------------------------------------------------------------

META_FIGURE = (
    CK
    / "nonfocal_figure2_samples.tsv"
)

META_REGION = (
    PROJECT
    / "data"
    / "derived"
    / "region_panel_metadata.tsv"
)

META_OPTIX = (
    PROJECT
    / "external"
    / "krista_genotype_plots_original"
    / "Genotype_Plots"
    / "Optix_metadata.csv"
)

META_KRISTA = (
    PROJECT
    / "external"
    / "krista_analysis"
    / "metadata.csv"
)

# ------------------------------------------------------------
# OUTPUTS
# ------------------------------------------------------------

OUT_SITES = (
    CK502
    / "band_51pass_vs_502fail_sites.tsv"
)

OUT_ALL = (
    CK502
    / "band_51pass_vs_502fail_individuals.tsv"
)

OUT_ERRORS = (
    CK502
    / "band_51pass_vs_502fail_misclassified.tsv"
)

OUT_INDIVIDUAL = (
    CK502
    / "band_51pass_vs_502fail_individual_summary.tsv"
)

OUT_TAXA = (
    CK502
    / "band_51pass_vs_502fail_error_taxa_summary.tsv"
)


def fail(message):
    raise RuntimeError(message)


def emptyish(x):
    if x is None:
        return True

    return str(x).strip() in {
        "",
        ".",
        "NA",
        "NaN",
        "nan",
        "None",
        "none",
    }


def norm(x):
    if emptyish(x):
        return ""

    return str(x).strip()


def norm_allele(x):
    return norm(x).upper()


def fnum(x):
    if emptyish(x):
        return np.nan

    try:
        return float(x)
    except Exception:
        return np.nan


def truthy(x):
    return (
        str(x).strip().lower()
        in {
            "1",
            "true",
            "t",
            "yes",
            "y",
        }
    )


def phenotype_binary(x):
    if emptyish(x):
        return np.nan

    text = (
        str(x)
        .strip()
        .lower()
    )

    if text in {
        "0",
        "0.0",
        "absent",
        "absence",
        "no",
        "false",
    }:
        return 0

    if text in {
        "1",
        "1.0",
        "present",
        "presence",
        "yes",
        "true",
    }:
        return 1

    return np.nan


def read_table(path, sep="\t"):
    if not path.exists():
        fail(
            f"Missing required file: {path}"
        )

    try:
        return pd.read_csv(
            path,
            sep=sep,
            dtype=str,
            keep_default_na=False,
            engine="python",
        )
    except Exception as exc:
        fail(
            f"Could not read {path}: {exc}"
        )


def vcf_header_samples(path):
    if not path.exists():
        fail(
            f"Missing VCF: {path}"
        )

    with path.open(
        "rt",
        errors="replace",
    ) as handle:

        for line in handle:

            if line.startswith(
                "#CHROM\t"
            ):
                fields = (
                    line.rstrip("\n")
                    .split("\t")
                )

                return fields[9:]

    fail(
        f"No #CHROM header in {path}"
    )


def parse_vcf_selected_sites(
    path,
    wanted_positions,
):
    """
    Direct parser: this avoids requiring an index.

    Returns:
      sample list
      dict[(chrom,pos,ref,alt)] -> {
          sample -> FORMAT dictionary
      }
    """

    samples = None

    sites = {}

    with path.open(
        "rt",
        errors="replace",
    ) as handle:

        for line in handle:

            if line.startswith("##"):
                continue

            if line.startswith(
                "#CHROM\t"
            ):
                header = (
                    line.rstrip("\n")
                    .split("\t")
                )

                samples = header[9:]

                continue

            if line.startswith("#"):
                continue

            fields = (
                line.rstrip("\n")
                .split("\t")
            )

            if len(fields) < 10:
                continue

            chrom = fields[0]

            try:
                pos = int(
                    fields[1]
                )
            except Exception:
                continue

            if pos not in wanted_positions:
                continue

            ref = norm_allele(
                fields[3]
            )

            alt = norm_allele(
                fields[4]
            )

            fmt_names = (
                fields[8]
                .split(":")
            )

            calls = {}

            for sample, sample_field in zip(
                samples,
                fields[9:],
            ):

                values = (
                    sample_field
                    .split(":")
                )

                fmt = {
                    name: (
                        values[index]
                        if index
                        < len(values)
                        else ""
                    )
                    for index, name
                    in enumerate(
                        fmt_names
                    )
                }

                calls[
                    sample
                ] = fmt

            sites[
                (
                    chrom,
                    pos,
                    ref,
                    alt,
                )
            ] = calls

    if samples is None:
        fail(
            f"No VCF header parsed from {path}"
        )

    return (
        samples,
        sites,
    )


def parse_gt(gt):
    if emptyish(gt):
        return None

    text = str(
        gt
    ).strip()

    if text in {
        ".",
        "./.",
        ".|.",
    }:
        return None

    parts = re.split(
        r"[/|]",
        text,
    )

    if len(parts) != 2:
        return None

    if any(
        x not in {"0", "1"}
        for x in parts
    ):
        return None

    return tuple(
        int(x)
        for x in parts
    )


def parse_ad(ad):
    if emptyish(ad):
        return (
            np.nan,
            np.nan,
        )

    parts = (
        str(ad)
        .split(",")
    )

    if len(parts) < 2:
        return (
            np.nan,
            np.nan,
        )

    try:
        return (
            float(parts[0]),
            float(parts[1]),
        )
    except Exception:
        return (
            np.nan,
            np.nan,
        )


def predicted_band_state(
    gt,
    orientation,
):
    """
    Return predicted biological Band phenotype:
      0 = Band absent
      1 = Band present

    This follows the selected dominant orientation,
    not focal expectation.

    For Band under the biological/reference model:
      REF carrier -> Band present
      ALT/ALT     -> Band absent.
    """

    parsed = parse_gt(
        gt
    )

    if parsed is None:
        return np.nan

    alt_copies = sum(
        parsed
    )

    orientation = (
        str(orientation)
        .upper()
    )

    if (
        "REF_CARRIER"
        in orientation
    ):
        return (
            0
            if alt_copies == 2
            else 1
        )

    if (
        "ALT_CARRIER"
        in orientation
    ):
        return (
            0
            if alt_copies == 0
            else 1
        )

    # Band biological fallback:
    # reference = presence.
    return (
        0
        if alt_copies == 2
        else 1
    )


def metadata_lookup(
    df,
    id_columns,
    value_columns,
    prefix,
):
    """
    A sample can appear under SequenceID,
    Accession, EarthCapeID, etc.

    Build one lookup over every available identifier.
    """

    lookup = {}

    available_ids = [
        x
        for x in id_columns
        if x in df.columns
    ]

    available_values = [
        x
        for x in value_columns
        if x in df.columns
    ]

    for _, row in df.iterrows():

        ids = {
            norm(
                row[col]
            )
            for col
            in available_ids
            if not emptyish(
                row[col]
            )
        }

        if not ids:
            continue

        record = {
            f"{prefix}{col}":
                norm(
                    row[col]
                )
            for col
            in available_values
        }

        for identifier in ids:

            if (
                identifier
                not in lookup
            ):
                lookup[
                    identifier
                ] = record.copy()

            else:
                for key, value in (
                    record.items()
                ):

                    if not value:
                        continue

                    old = lookup[
                        identifier
                    ].get(
                        key,
                        "",
                    )

                    if not old:
                        lookup[
                            identifier
                        ][key] = value

                    elif (
                        value
                        not in old.split("|")
                    ):
                        lookup[
                            identifier
                        ][key] = (
                            old
                            + "|"
                            + value
                        )

    return lookup


def add_lookup(
    record,
    lookup,
    sample,
):
    data = lookup.get(
        sample,
        {}
    )

    for key, value in (
        data.items()
    ):
        record[key] = value


def main():

    # ========================================================
    # 1. FIND BAND SITES THAT PASS IN 51-SAMPLE ANALYSIS
    # ========================================================

    s51 = read_table(
        SITES51
    )

    s51["minimum_BA_num"] = (
        pd.to_numeric(
            s51["minimum_BA"],
            errors="coerce",
        )
    )

    band51 = s51.loc[
        (
            s51["panel"]
            .str.lower()
            == "band"
        )
        &
        (
            s51[
                "minimum_BA_num"
            ] >= 0.80
        )
    ].copy()

    if band51.empty:
        fail(
            "No 51-sample Band sites "
            "with minimum BA >= 0.80."
        )

    band51["pos_num"] = (
        pd.to_numeric(
            band51["pos"],
            errors="raise",
        )
        .astype(int)
    )

    # ========================================================
    # 2. FULL-502 SCORE AT THOSE SAME SITES
    # ========================================================

    s502 = read_table(
        SITES502
    )

    band502 = s502.loc[
        (
            s502["panel"]
            .str.lower()
            == "band"
        )
    ].copy()

    band502["pos_num"] = (
        pd.to_numeric(
            band502["pos"],
            errors="coerce",
        )
    )

    band502 = band502.loc[
        band502[
            "pos_num"
        ].notna()
    ].copy()

    band502["pos_num"] = (
        band502[
            "pos_num"
        ].astype(int)
    )

    key_columns = [
        "chrom",
        "pos_num",
        "ref",
        "alt",
    ]

    for df in [
        band51,
        band502,
    ]:
        df["ref"] = (
            df["ref"]
            .str.upper()
        )
        df["alt"] = (
            df["alt"]
            .str.upper()
        )

    # The paired table should already be one row/site.
    merged = band51.merge(
        band502,
        on=key_columns,
        how="left",
        suffixes=(
            "_51",
            "_502",
        ),
        indicator=True,
    )

    merged[
        "shared_exact_site"
    ] = (
        merged["_merge"]
        == "both"
    )

    merged[
        "BA_51"
    ] = pd.to_numeric(
        merged[
            "minimum_BA_num"
        ],
        errors="coerce",
    )

    merged[
        "BA_502"
    ] = pd.to_numeric(
        merged.get(
            "minimum_BA_502",
            np.nan,
        ),
        errors="coerce",
    )

    merged[
        "passes_51_BA80"
    ] = (
        merged[
            "BA_51"
        ] >= 0.80
    )

    merged[
        "passes_502_BA80"
    ] = (
        merged[
            "BA_502"
        ] >= 0.80
    )

    merged[
        "comparison_class"
    ] = np.where(
        ~merged[
            "shared_exact_site"
        ],
        "51_pass_site_absent_from_502_scored_site_universe",
        np.where(
            merged[
                "passes_502_BA80"
            ],
            "passes_BA80_in_both",
            "passes_51_fails_502",
        ),
    )

    # Only these can be interrogated individual-by-individual.
    target_sites = merged.loc[
        (
            merged[
                "comparison_class"
            ]
            == "passes_51_fails_502"
        )
    ].copy()

    if target_sites.empty:
        print(
            "\nIMPORTANT: no exact shared "
            "Band site passes BA80 in 51 "
            "and fails BA80 in 502."
        )

        print(
            "The discrepancy is entirely "
            "due to site-universe differences."
        )

    wanted_positions = set(
        target_sites[
            "pos_num"
        ].astype(int)
    )

    # ========================================================
    # 3. SAMPLE MEMBERSHIP
    # ========================================================

    samples51 = set(
        vcf_header_samples(
            VCF51
        )
    )

    samples502_header = (
        vcf_header_samples(
            VCF502
        )
    )

    # 51-sample exclusions.
    excluded51 = read_table(
        EXCLUDED51
    )

    excluded51_map = {
        norm(row["sample"]):
            norm(row["reason"])
        for _, row
        in excluded51.iterrows()
    }

    # ========================================================
    # 4. EXACT ORIGINAL-502 PHENOTYPE / COHORT MANIFEST
    # ========================================================

    cohort = read_table(
        COHORT502
    )

    band_cohort = cohort.loc[
        (
            cohort["panel"]
            .str.lower()
            == "band"
        )
        |
        (
            cohort["trait"]
            .str.upper()
            == "B"
        )
    ].copy()

    band_cohort[
        "phenotype_binary"
    ] = band_cohort[
        "phenotype"
    ].map(
        phenotype_binary
    )

    band_cohort[
        "reason_lower"
    ] = (
        band_cohort[
            "reason"
        ]
        .astype(str)
        .str.lower()
    )

    # Prefer explicit cohort reason if present.
    explicit_comparative = (
        band_cohort[
            "reason_lower"
        ]
        .str.contains(
            "comparative",
            regex=False,
        )
    )

    if (
        explicit_comparative.sum()
        > 0
    ):
        band_cohort[
            "used_in_502_scoring"
        ] = (
            explicit_comparative
            &
            band_cohort[
                "phenotype_binary"
            ].notna()
        )
    else:
        # Defensive fallback.
        band_cohort[
            "used_in_502_scoring"
        ] = (
            band_cohort[
                "phenotype_binary"
            ].notna()
            &
            ~band_cohort[
                "reason_lower"
            ].str.contains(
                "focal|piano|exclude|missing",
                regex=True,
            )
        )

    cohort_map = {
        norm(
            row["sample"]
        ):
        row.to_dict()
        for _, row
        in band_cohort.iterrows()
    }

    print(
        "\nORIGINAL-502 COHORT MANIFEST "
        "REASON COUNTS"
    )

    print(
        band_cohort[
            "reason"
        ]
        .value_counts(
            dropna=False
        )
        .to_string()
    )

    print(
        "\nINFERRED ORIGINAL-502 "
        "COMPARATIVE SCORING SAMPLES:",
        int(
            band_cohort[
                "used_in_502_scoring"
            ].sum()
        ),
    )

    print(
        "phenotype 0:",
        int(
            (
                band_cohort[
                    "used_in_502_scoring"
                ]
                &
                (
                    band_cohort[
                        "phenotype_binary"
                    ]
                    == 0
                )
            ).sum()
        ),
    )

    print(
        "phenotype 1:",
        int(
            (
                band_cohort[
                    "used_in_502_scoring"
                ]
                &
                (
                    band_cohort[
                        "phenotype_binary"
                    ]
                    == 1
                )
            ).sum()
        ),
    )

    # ========================================================
    # 5. LOAD ALL METADATA SOURCES
    # ========================================================

    fig = read_table(
        META_FIGURE
    )

    fig_lookup = metadata_lookup(
        fig,
        id_columns=[
            "Accession",
        ],
        value_columns=[
            "Taxon Name",
            "Country",
            "Data Type",
            "B",
            "R",
            "D",
            "Yb-dorsal",
            "Yb-ventral",
            "Optix Figure",
            "note",
        ],
        prefix="figure_",
    )

    region = read_table(
        META_REGION
    )

    region_lookup = metadata_lookup(
        region,
        id_columns=[
            "SequenceID",
            "EarthCapeID",
            "Accession",
        ],
        value_columns=[
            "SequenceID",
            "EarthCapeID",
            "Taxon Name",
            "Country",
            "Latitude",
            "Longitude",
            "Accession",
            "Publication",
        ],
        prefix="region_",
    )

    optix_meta = read_table(
        META_OPTIX,
        sep=",",
    )

    optix_lookup = metadata_lookup(
        optix_meta,
        id_columns=[
            "SequenceID",
            "EarthCapeID",
            "Accession",
        ],
        value_columns=[
            "SequenceID",
            "EarthCapeID",
            "Taxon Name",
            "Country",
            "Latitude",
            "Longitude",
            "Accession",
            "Publication",
        ],
        prefix="optixmeta_",
    )

    krista = read_table(
        META_KRISTA,
        sep=",",
    )

    krista_lookup = metadata_lookup(
        krista,
        id_columns=[
            "Sample",
        ],
        value_columns=[
            "Sample",
            "Country",
            "Genus",
            "Species",
            "Subsp",
            "Yellow Bar ",
            "Red band ",
        ],
        prefix="krista_",
    )

    # ========================================================
    # 6. QUERY FULL 502 VCF AT SHARED PASS/FAIL SITES
    # ========================================================

    if wanted_positions:
        (
            parsed_samples502,
            calls502,
        ) = parse_vcf_selected_sites(
            VCF502,
            wanted_positions,
        )

        if (
            parsed_samples502
            != samples502_header
        ):
            fail(
                "502 VCF sample order changed "
                "between header reads."
            )
    else:
        calls502 = {}

    individual_rows = []

    reconstructed_site_rows = []

    for _, site in (
        target_sites.iterrows()
    ):

        chrom = site["chrom"]

        pos = int(
            site["pos_num"]
        )

        ref = norm_allele(
            site["ref"]
        )

        alt = norm_allele(
            site["alt"]
        )

        exact_key = (
            chrom,
            pos,
            ref,
            alt,
        )

        site_calls = (
            calls502.get(
                exact_key
            )
        )

        if site_calls is None:
            print(
                "WARNING: target site not found "
                "as exact REF/ALT in full VCF:",
                exact_key,
            )
            continue

        orientation = norm(
            site.get(
                "selected_orientation_502",
                "",
            )
        )

        if not orientation:
            orientation = norm(
                site.get(
                    "selected_orientation_51",
                    "",
                )
            )

        site_rows = []

        for sample in (
            samples502_header
        ):

            call = site_calls[
                sample
            ]

            gt = norm(
                call.get(
                    "GT",
                    "",
                )
            )

            ad = norm(
                call.get(
                    "AD",
                    "",
                )
            )

            dp = fnum(
                call.get(
                    "DP",
                    "",
                )
            )

            gq = fnum(
                call.get(
                    "GQ",
                    "",
                )
            )

            qc_pass = bool(
                parse_gt(gt)
                is not None
                and np.isfinite(dp)
                and dp >= 2
                and np.isfinite(gq)
                and gq >= 20
            )

            cohort_record = (
                cohort_map.get(
                    sample,
                    {}
                )
            )

            phenotype_raw = (
                cohort_record.get(
                    "phenotype",
                    "",
                )
            )

            phenotype = (
                phenotype_binary(
                    phenotype_raw
                )
            )

            reason = norm(
                cohort_record.get(
                    "reason",
                    "",
                )
            )

            taxon_manifest = norm(
                cohort_record.get(
                    "taxon",
                    "",
                )
            )

            data_type_manifest = norm(
                cohort_record.get(
                    "data_type",
                    "",
                )
            )

            used502 = bool(
                cohort_record.get(
                    "used_in_502_scoring",
                    False,
                )
            )

            pred_raw = (
                predicted_band_state(
                    gt,
                    orientation,
                )
            )

            pred_qc = (
                pred_raw
                if qc_pass
                else np.nan
            )

            if (
                used502
                and qc_pass
                and np.isfinite(
                    phenotype
                )
            ):
                misclassified = bool(
                    int(pred_qc)
                    != int(
                        phenotype
                    )
                )
            else:
                misclassified = np.nan

            if (
                not used502
            ):
                error_type = (
                    "not_in_502_scoring"
                )

            elif (
                not qc_pass
            ):
                error_type = (
                    "QC_missing"
                )

            elif (
                not np.isfinite(
                    phenotype
                )
            ):
                error_type = (
                    "phenotype_missing"
                )

            elif not misclassified:
                error_type = (
                    "concordant"
                )

            elif (
                int(
                    phenotype
                )
                == 0
                and int(
                    pred_qc
                )
                == 1
            ):
                error_type = (
                    "scored_absent_but_genotype_predicts_present"
                )

            elif (
                int(
                    phenotype
                )
                == 1
                and int(
                    pred_qc
                )
                == 0
            ):
                error_type = (
                    "scored_present_but_genotype_predicts_absent"
                )

            else:
                error_type = (
                    "other_mismatch"
                )

            ref_depth, alt_depth = (
                parse_ad(
                    ad
                )
            )

            if (
                np.isfinite(
                    ref_depth
                )
                and np.isfinite(
                    alt_depth
                )
                and (
                    ref_depth
                    + alt_depth
                ) > 0
            ):
                alt_fraction = (
                    alt_depth
                    / (
                        ref_depth
                        + alt_depth
                    )
                )
            else:
                alt_fraction = (
                    np.nan
                )

            in51 = (
                sample
                in samples51
            )

            excluded51_reason = (
                excluded51_map.get(
                    sample,
                    "",
                )
            )

            used51 = bool(
                in51
                and not sample.startswith(
                    "TN"
                )
                and not excluded51_reason
            )

            record = {
                "chrom":
                    chrom,

                "pos":
                    pos,

                "ref":
                    ref,

                "alt":
                    alt,

                "BA_51":
                    fnum(
                        site[
                            "BA_51"
                        ]
                    ),

                "BA_502":
                    fnum(
                        site[
                            "BA_502"
                        ]
                    ),

                "BA_drop_51_to_502":
                    (
                        fnum(
                            site[
                                "BA_51"
                            ]
                        )
                        -
                        fnum(
                            site[
                                "BA_502"
                            ]
                        )
                    ),

                "orientation_502":
                    orientation,

                "sample":
                    sample,

                "in_51sample_VCF":
                    in51,

                "used_in_51sample_scoring":
                    used51,

                "excluded_51_reason":
                    excluded51_reason,

                "used_in_502_scoring":
                    used502,

                "502_cohort_reason":
                    reason,

                "B_phenotype_used_for_502_score":
                    (
                        int(
                            phenotype
                        )
                        if np.isfinite(
                            phenotype
                        )
                        else np.nan
                    ),

                "B_phenotype_label":
                    (
                        "present"
                        if (
                            np.isfinite(
                                phenotype
                            )
                            and int(
                                phenotype
                            )
                            == 1
                        )
                        else (
                            "absent"
                            if (
                                np.isfinite(
                                    phenotype
                                )
                                and int(
                                    phenotype
                                )
                                == 0
                            )
                            else "missing"
                        )
                    ),

                "manifest_taxon":
                    taxon_manifest,

                "manifest_data_type":
                    data_type_manifest,

                "GT":
                    gt,

                "AD":
                    ad,

                "DP":
                    dp,

                "GQ":
                    gq,

                "ALT_fraction":
                    alt_fraction,

                "passes_DP2_GQ20":
                    qc_pass,

                "genotype_predicted_B_state":
                    (
                        int(
                            pred_qc
                        )
                        if np.isfinite(
                            pred_qc
                        )
                        else np.nan
                    ),

                "genotype_predicted_B_label":
                    (
                        "present"
                        if (
                            np.isfinite(
                                pred_qc
                            )
                            and int(
                                pred_qc
                            )
                            == 1
                        )
                        else (
                            "absent"
                            if (
                                np.isfinite(
                                    pred_qc
                                )
                                and int(
                                    pred_qc
                                )
                                == 0
                            )
                            else "missing"
                        )
                    ),

                "misclassified_for_B":
                    misclassified,

                "error_type":
                    error_type,
            }

            add_lookup(
                record,
                fig_lookup,
                sample,
            )

            add_lookup(
                record,
                region_lookup,
                sample,
            )

            add_lookup(
                record,
                optix_lookup,
                sample,
            )

            add_lookup(
                record,
                krista_lookup,
                sample,
            )

            individual_rows.append(
                record
            )

            site_rows.append(
                record
            )

        site_df = pd.DataFrame(
            site_rows
        )

        scoring = site_df.loc[
            site_df[
                "used_in_502_scoring"
            ]
            == True
        ].copy()

        callable_scoring = (
            scoring.loc[
                scoring[
                    "passes_DP2_GQ20"
                ]
                == True
            ]
        )

        error_n = int(
            (
                callable_scoring[
                    "misclassified_for_B"
                ]
                == True
            ).sum()
        )

        expected_qc_errors = fnum(
            site.get(
                "QC_errors_502",
                np.nan,
            )
        )

        reconstructed_site_rows.append(
            {
                "chrom":
                    chrom,

                "pos":
                    pos,

                "ref":
                    ref,

                "alt":
                    alt,

                "BA_51":
                    fnum(
                        site[
                            "BA_51"
                        ]
                    ),

                "BA_502":
                    fnum(
                        site[
                            "BA_502"
                        ]
                    ),

                "BA_drop_51_to_502":
                    (
                        fnum(
                            site[
                                "BA_51"
                            ]
                        )
                        -
                        fnum(
                            site[
                                "BA_502"
                            ]
                        )
                    ),

                "orientation_51":
                    norm(
                        site.get(
                            "orientation_51",
                            "",
                        )
                    ),

                "orientation_502":
                    orientation,

                "scoring_samples_reconstructed":
                    len(
                        scoring
                    ),

                "QC_callable_scoring_samples":
                    len(
                        callable_scoring
                    ),

                "reconstructed_QC_errors":
                    error_n,

                "paired_table_QC_errors":
                    expected_qc_errors,

                "error_count_matches_paired_table":
                    (
                        (
                            error_n
                            == int(
                                expected_qc_errors
                            )
                        )
                        if np.isfinite(
                            expected_qc_errors
                        )
                        else np.nan
                    ),

                "errors_from_samples_not_in_51":
                    int(
                        (
                            (
                                callable_scoring[
                                    "misclassified_for_B"
                                ]
                                == True
                            )
                            &
                            (
                                callable_scoring[
                                    "in_51sample_VCF"
                                ]
                                == False
                            )
                        ).sum()
                    ),

                "errors_from_samples_in_51":
                    int(
                        (
                            (
                                callable_scoring[
                                    "misclassified_for_B"
                                ]
                                == True
                            )
                            &
                            (
                                callable_scoring[
                                    "in_51sample_VCF"
                                ]
                                == True
                            )
                        ).sum()
                    ),
            }
        )

    all_individuals = pd.DataFrame(
        individual_rows
    )

    reconstructed = pd.DataFrame(
        reconstructed_site_rows
    )

    # ========================================================
    # 7. SITE SUMMARY
    # ========================================================

    site_columns = [
        "chrom",
        "pos_num",
        "ref",
        "alt",
        "BA_51",
        "BA_502",
        "comparison_class",
    ]

    site_summary = merged[
        [
            col
            for col
            in site_columns
            if col in merged.columns
        ]
    ].copy()

    site_summary = site_summary.rename(
        columns={
            "pos_num": "pos",
        }
    )

    if not reconstructed.empty:
        site_summary = (
            site_summary.merge(
                reconstructed,
                on=[
                    "chrom",
                    "pos",
                    "ref",
                    "alt",
                ],
                how="left",
                suffixes=(
                    "",
                    "_audit",
                ),
            )
        )

    site_summary.to_csv(
        OUT_SITES,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    # ========================================================
    # 8. LONG ALL-INDIVIDUAL TABLE
    # ========================================================

    if not all_individuals.empty:

        all_individuals = (
            all_individuals.sort_values(
                [
                    "pos",
                    "B_phenotype_used_for_502_score",
                    "misclassified_for_B",
                    "manifest_taxon",
                    "sample",
                ],
                ascending=[
                    True,
                    True,
                    False,
                    True,
                    True,
                ],
            )
        )

        all_individuals.to_csv(
            OUT_ALL,
            sep="\t",
            index=False,
            na_rep="NA",
        )

        errors = (
            all_individuals.loc[
                (
                    all_individuals[
                        "used_in_502_scoring"
                    ]
                    == True
                )
                &
                (
                    all_individuals[
                        "passes_DP2_GQ20"
                    ]
                    == True
                )
                &
                (
                    all_individuals[
                        "misclassified_for_B"
                    ]
                    == True
                )
            ]
            .copy()
        )

        errors.to_csv(
            OUT_ERRORS,
            sep="\t",
            index=False,
            na_rep="NA",
        )

        # ====================================================
        # 9. PER-INDIVIDUAL ERROR SUMMARY
        # ====================================================

        scoring_callable = (
            all_individuals.loc[
                (
                    all_individuals[
                        "used_in_502_scoring"
                    ]
                    == True
                )
                &
                (
                    all_individuals[
                        "passes_DP2_GQ20"
                    ]
                    == True
                )
            ]
            .copy()
        )

        individual = (
            scoring_callable.groupby(
                "sample",
                dropna=False,
            )
            .agg(
                sites_callable=(
                    "pos",
                    "size",
                ),
                sites_misclassified=(
                    "misclassified_for_B",
                    lambda x:
                        int(
                            (
                                x
                                == True
                            ).sum()
                        ),
                ),
                in_51sample_VCF=(
                    "in_51sample_VCF",
                    "first",
                ),
                used_in_51sample_scoring=(
                    "used_in_51sample_scoring",
                    "first",
                ),
                B_phenotype=(
                    "B_phenotype_used_for_502_score",
                    "first",
                ),
                taxon=(
                    "manifest_taxon",
                    "first",
                ),
                data_type=(
                    "manifest_data_type",
                    "first",
                ),
            )
            .reset_index()
        )

        individual[
            "error_fraction"
        ] = (
            individual[
                "sites_misclassified"
            ]
            /
            individual[
                "sites_callable"
            ]
        )

        # Pull useful metadata from the first site row.
        first_metadata = (
            all_individuals
            .drop_duplicates(
                "sample"
            )
            .set_index(
                "sample"
            )
        )

        metadata_columns = [
            x
            for x
            in all_individuals.columns
            if (
                x.startswith(
                    "figure_"
                )
                or x.startswith(
                    "region_"
                )
                or x.startswith(
                    "optixmeta_"
                )
                or x.startswith(
                    "krista_"
                )
                or x
                == "excluded_51_reason"
            )
        ]

        for col in metadata_columns:
            individual[
                col
            ] = individual[
                "sample"
            ].map(
                first_metadata[
                    col
                ]
            )

        individual = (
            individual.sort_values(
                [
                    "sites_misclassified",
                    "error_fraction",
                    "in_51sample_VCF",
                    "taxon",
                    "sample",
                ],
                ascending=[
                    False,
                    False,
                    True,
                    True,
                    True,
                ],
            )
        )

        individual.to_csv(
            OUT_INDIVIDUAL,
            sep="\t",
            index=False,
            na_rep="NA",
        )

        # ====================================================
        # 10. TAXON/COUNTRY ERROR SUMMARY
        # ====================================================

        if not errors.empty:

            errors[
                "best_taxon"
            ] = errors[
                "manifest_taxon"
            ]

            for fallback in [
                "figure_Taxon Name",
                "region_Taxon Name",
                "optixmeta_Taxon Name",
            ]:
                if (
                    fallback
                    not in errors.columns
                ):
                    continue

                blank = (
                    errors[
                        "best_taxon"
                    ]
                    .astype(str)
                    .str.strip()
                    == ""
                )

                errors.loc[
                    blank,
                    "best_taxon",
                ] = errors.loc[
                    blank,
                    fallback,
                ]

            errors[
                "best_country"
            ] = ""

            for source in [
                "figure_Country",
                "region_Country",
                "optixmeta_Country",
                "krista_Country",
            ]:
                if (
                    source
                    not in errors.columns
                ):
                    continue

                blank = (
                    errors[
                        "best_country"
                    ]
                    .astype(str)
                    .str.strip()
                    == ""
                )

                errors.loc[
                    blank,
                    "best_country",
                ] = errors.loc[
                    blank,
                    source,
                ]

            taxa = (
                errors.groupby(
                    [
                        "best_taxon",
                        "best_country",
                        "B_phenotype_label",
                        "in_51sample_VCF",
                        "error_type",
                    ],
                    dropna=False,
                )
                .agg(
                    error_rows=(
                        "sample",
                        "size",
                    ),
                    unique_individuals=(
                        "sample",
                        "nunique",
                    ),
                    unique_sites=(
                        "pos",
                        "nunique",
                    ),
                )
                .reset_index()
                .sort_values(
                    [
                        "error_rows",
                        "unique_individuals",
                    ],
                    ascending=False,
                )
            )

            taxa.to_csv(
                OUT_TAXA,
                sep="\t",
                index=False,
                na_rep="NA",
            )

    # ========================================================
    # 11. PRINT A HUMAN-READABLE AUDIT
    # ========================================================

    print(
        "\n51-SAMPLE BAND BA>=0.80 SITES"
    )

    print(
        site_summary.to_string(
            index=False
        )
    )

    if not all_individuals.empty:

        print(
            "\nTOP INDIVIDUALS DRIVING "
            "FULL-502 BAND ERRORS"
        )

        show = pd.read_csv(
            OUT_INDIVIDUAL,
            sep="\t",
            dtype=str,
        ).head(40)

        preferred = [
            "sample",
            "taxon",
            "B_phenotype",
            "sites_callable",
            "sites_misclassified",
            "error_fraction",
            "in_51sample_VCF",
            "used_in_51sample_scoring",
            "figure_Country",
            "krista_Subsp",
            "krista_Red band ",
            "excluded_51_reason",
        ]

        preferred = [
            x
            for x
            in preferred
            if x in show.columns
        ]

        print(
            show[
                preferred
            ].to_string(
                index=False
            )
        )

        errors = pd.read_csv(
            OUT_ERRORS,
            sep="\t",
            dtype=str,
        )

        print(
            "\nTOTAL MISCLASSIFIED "
            "SITE x INDIVIDUAL ROWS:",
            len(errors),
        )

        if not errors.empty:

            print(
                "\nERRORS BY 51-SAMPLE MEMBERSHIP"
            )

            print(
                errors[
                    "in_51sample_VCF"
                ]
                .value_counts(
                    dropna=False
                )
                .to_string()
            )

            print(
                "\nMOST COMMON ERROR TAXA"
            )

            taxa = pd.read_csv(
                OUT_TAXA,
                sep="\t",
                dtype=str,
            )

            print(
                taxa.head(
                    30
                ).to_string(
                    index=False
                )
            )

    print(
        "\nOutputs:"
    )

    for path in [
        OUT_SITES,
        OUT_ALL,
        OUT_ERRORS,
        OUT_INDIVIDUAL,
        OUT_TAXA,
    ]:
        print(
            f"  {path}"
        )


if __name__ == "__main__":
    main()
