#!/usr/bin/env python3

from pathlib import Path
import re
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

CK = (
    PROJECT
    / "analysis"
    / "original_502"
    / "checkpoints"
)

GENOTYPE_DIR = (
    PROJECT
    / "external"
    / "krista_genotype_plots_original"
    / "Genotype_Plots"
)

SITE_FILE = (
    CK
    / "frozen_comparative_sites_biological_polarity.tsv"
)

EXPECTATION_FILE = (
    CK
    / "focal_module_expectations.tsv"
)

VCF_BY_MODULE = {
    "Band":
        GENOTYPE_DIR
        / "optix_band-1112651-1114314.vcf",

    "Dennis1":
        GENOTYPE_DIR
        / "dennis1-1215210-1216124.vcf",

    "Dennis2":
        GENOTYPE_DIR
        / "dennis2-1218557-1221494.vcf",

    "Ray1":
        GENOTYPE_DIR
        / "ray1-1156746-1160220.vcf",

    "Ray2":
        GENOTYPE_DIR
        / "ray2-1174547-1187833.vcf",

    "Yb":
        GENOTYPE_DIR
        / "Ivory_Hmr15-1627000-1637000_biallelic.vcf",
}

MODULE_ORDER = [
    "Band",
    "Ray1",
    "Ray2",
    "Dennis1",
    "Dennis2",
    "Yb",
]

# ------------------------------------------------------------
# LOCKED BIOLOGICAL POLARITY
#
# This comes from the known state of the reference genome,
# NOT from REF/ALT statistical association orientation.
#
# Band: reference is homozygous presence -> REF=P, ALT=A
# D/R:  reference is homozygous absence  -> REF=A, ALT=P
# Yb:   reference is homozygous presence -> REF=P, ALT=A
# ------------------------------------------------------------

REFERENCE_STATE = {
    "Band": "presence",
    "Dennis1": "absence",
    "Dennis2": "absence",
    "Ray1": "absence",
    "Ray2": "absence",
    "Yb": "presence",
}

OUT_SITES = (
    CK
    / "frozen_comparative_sites_reference_polarized.tsv"
)

OUT_RAW = (
    CK
    / "original502_focal_calls_at_frozen_sites.raw.tsv"
)

OUT_LONG = (
    CK
    / "focal_calls_at_frozen_sites_biological.long.tsv"
)

OUT_WIDE = (
    CK
    / "focal_calls_at_frozen_sites_biological.wide.tsv"
)

OUT_SUMMARY = (
    CK
    / "focal_calls_at_frozen_sites_biological.summary.tsv"
)


def fail(message):
    raise RuntimeError(message)


def emptyish(value):
    if value is None:
        return True

    return str(value).strip() in {
        "",
        ".",
        "NA",
        "NaN",
        "nan",
        "None",
    }


def as_float(value):
    if emptyish(value):
        return np.nan

    try:
        return float(value)
    except Exception:
        return np.nan


def normalise_allele(value):
    if emptyish(value):
        return None

    return str(value).strip().upper()


def parse_gt(gt):
    if emptyish(gt):
        return None

    parts = re.split(
        r"[/|]",
        str(gt).strip(),
    )

    if len(parts) != 2:
        return None

    if any(
        allele not in {"0", "1"}
        for allele in parts
    ):
        return None

    return sum(
        int(allele)
        for allele in parts
    )


def parse_ad(ad):
    if emptyish(ad):
        return (
            np.nan,
            np.nan,
        )

    parts = str(ad).split(",")

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


def genotype_from_presence_copies(copies):
    return {
        0: "AA",
        1: "PA",
        2: "PP",
    }.get(
        copies,
        None,
    )


def dosage_from_genotype(gt):
    return {
        "AA": 0.0,
        "PA": 0.5,
        "PP": 1.0,
    }.get(
        gt,
        np.nan,
    )


def query_vcf(module, path):
    if not path.exists():
        fail(
            f"Missing original 502-sample VCF for "
            f"{module}: {path}"
        )

    fmt = (
        "%CHROM\\t%POS\\t%REF\\t%ALT"
        "[\\t%GT\\t%AD\\t%DP\\t%GQ]\\n"
    )

    command = [
        "bcftools",
        "query",
        "-s",
        "TN_Left,TN_Right",
        "-f",
        fmt,
        str(path),
    ]

    result = subprocess.run(
        command,
        text=True,
        capture_output=True,
    )

    if result.returncode != 0:
        print(
            result.stderr,
            file=sys.stderr,
        )

        fail(
            f"bcftools query failed for {module}"
        )

    rows = []

    for line in result.stdout.splitlines():
        fields = line.split("\t")

        if len(fields) != 12:
            fail(
                f"Unexpected bcftools field count "
                f"({len(fields)}) for {module}: {line}"
            )

        (
            chrom,
            pos,
            ref,
            alt,
            left_gt,
            left_ad,
            left_dp,
            left_gq,
            right_gt,
            right_ad,
            right_dp,
            right_gq,
        ) = fields

        rows.append(
            {
                "module": module,
                "chrom": chrom,
                "pos": int(pos),
                "ref": normalise_allele(ref),
                "alt": normalise_allele(alt),
                "left_GT": left_gt,
                "left_AD": left_ad,
                "left_DP": left_dp,
                "left_GQ": left_gq,
                "right_GT": right_gt,
                "right_AD": right_ad,
                "right_DP": right_dp,
                "right_GQ": right_gq,
                "source_vcf": str(path),
            }
        )

    return pd.DataFrame(rows)


def biological_polarity(module, ref, alt):
    ref = normalise_allele(ref)
    alt = normalise_allele(alt)

    if REFERENCE_STATE[module] == "presence":
        return {
            "presence_allele": ref,
            "absence_allele": alt,
            "biological_orientation":
                "REF=presence;ALT=absence",
        }

    return {
        "presence_allele": alt,
        "absence_allele": ref,
        "biological_orientation":
            "REF=absence;ALT=presence",
    }


def evaluate_side(
    module,
    side,
    raw_row,
    presence_allele,
    expected_gt,
):
    prefix = (
        "left"
        if side == "Left"
        else "right"
    )

    raw_gt = raw_row[
        f"{prefix}_GT"
    ]

    ad = raw_row[
        f"{prefix}_AD"
    ]

    dp = as_float(
        raw_row[
            f"{prefix}_DP"
        ]
    )

    gq = as_float(
        raw_row[
            f"{prefix}_GQ"
        ]
    )

    ref = normalise_allele(
        raw_row["ref"]
    )

    alt = normalise_allele(
        raw_row["alt"]
    )

    alt_copies = parse_gt(
        raw_gt
    )

    ref_depth, alt_depth = parse_ad(
        ad
    )

    presence_copies = None
    presence_fraction = np.nan

    if presence_allele == ref:
        if alt_copies is not None:
            presence_copies = (
                2 - alt_copies
            )

        if (
            np.isfinite(ref_depth)
            and np.isfinite(alt_depth)
            and ref_depth + alt_depth > 0
        ):
            presence_fraction = (
                ref_depth
                / (
                    ref_depth
                    + alt_depth
                )
            )

    elif presence_allele == alt:
        if alt_copies is not None:
            presence_copies = (
                alt_copies
            )

        if (
            np.isfinite(ref_depth)
            and np.isfinite(alt_depth)
            and ref_depth + alt_depth > 0
        ):
            presence_fraction = (
                alt_depth
                / (
                    ref_depth
                    + alt_depth
                )
            )

    biological_gt = (
        genotype_from_presence_copies(
            presence_copies
        )
    )

    qc_pass = bool(
        biological_gt is not None
        and np.isfinite(dp)
        and dp >= 2
        and np.isfinite(gq)
        and gq >= 20
    )

    if not qc_pass:
        display_gt = "missing"
        dosage = np.nan
    else:
        display_gt = biological_gt
        dosage = dosage_from_genotype(
            biological_gt
        )

    match = (
        biological_gt == expected_gt
        if qc_pass
        else np.nan
    )

    return {
        "side": side,
        "GT_raw": raw_gt,
        "AD_raw": ad,
        "DP": dp,
        "GQ": gq,
        "QC_pass_DP2_GQ20": qc_pass,
        "biological_genotype":
            display_gt,
        "presence_dosage":
            dosage,
        "presence_allele_fraction":
            presence_fraction
            if qc_pass
            else np.nan,
        "expected_genotype":
            expected_gt,
        "matches_expectation":
            match,
    }


def main():
    if shutil.which("bcftools") is None:
        fail(
            "bcftools is not available in PATH."
        )

    sites = pd.read_csv(
        SITE_FILE,
        sep="\t",
        dtype=str,
    )

    expectations = pd.read_csv(
        EXPECTATION_FILE,
        sep="\t",
        dtype=str,
    )

    expectation_map = (
        expectations
        .set_index("module")
        .to_dict(
            orient="index"
        )
    )

    # --------------------------------------------------------
    # Query the actual original pooled regional VCFs.
    # --------------------------------------------------------

    raw_frames = []

    for module in MODULE_ORDER:
        raw_frames.append(
            query_vcf(
                module,
                VCF_BY_MODULE[module],
            )
        )

    raw = pd.concat(
        raw_frames,
        ignore_index=True,
    )

    raw.to_csv(
        OUT_RAW,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    # Fast lookup, preserving possible same-position alternate alleles.
    raw_lookup = {}

    raw_position_lookup = {}

    for _, row in raw.iterrows():
        exact_key = (
            row["module"],
            row["chrom"],
            int(row["pos"]),
            row["ref"],
            row["alt"],
        )

        raw_lookup[
            exact_key
        ] = row

        position_key = (
            row["module"],
            row["chrom"],
            int(row["pos"]),
        )

        raw_position_lookup.setdefault(
            position_key,
            []
        ).append(row)

    corrected_sites = []

    long_rows = []

    mismatch_rows = []

    for _, site in sites.iterrows():
        module = site["module"]
        chrom = str(site["chrom"])
        pos = int(float(site["pos"]))

        ref = normalise_allele(
            site["ref"]
        )

        alt = normalise_allele(
            site["alt"]
        )

        polarity = biological_polarity(
            module,
            ref,
            alt,
        )

        presence = polarity[
            "presence_allele"
        ]

        absence = polarity[
            "absence_allele"
        ]

        corrected_site = (
            site.to_dict()
        )

        # Replace the earlier statistical-orientation-based
        # polarity with the LOCKED biological polarity.
        corrected_site[
            "presence_associated_allele"
        ] = presence

        corrected_site[
            "absence_associated_allele"
        ] = absence

        corrected_site[
            "biological_orientation"
        ] = polarity[
            "biological_orientation"
        ]

        corrected_site[
            "biological_polarity_basis"
        ] = (
            "known_reference_phenotype_state"
        )

        corrected_sites.append(
            corrected_site
        )

        exact_key = (
            module,
            chrom,
            pos,
            ref,
            alt,
        )

        raw_row = raw_lookup.get(
            exact_key
        )

        position_key = (
            module,
            chrom,
            pos,
        )

        if raw_row is None:
            position_rows = (
                raw_position_lookup.get(
                    position_key,
                    []
                )
            )

            reason = (
                "site_absent_from_original502_vcf"
                if len(position_rows) == 0
                else "allele_mismatch_at_same_position"
            )

            mismatch_rows.append(
                {
                    "module": module,
                    "chrom": chrom,
                    "pos": pos,
                    "selected_ref": ref,
                    "selected_alt": alt,
                    "reason": reason,
                    "original502_alleles": (
                        ";".join(
                            f"{x['ref']}>{x['alt']}"
                            for x in position_rows
                        )
                        if position_rows
                        else "none"
                    ),
                }
            )

            for side in [
                "Left",
                "Right",
            ]:
                expected_gt = (
                    expectation_map[module][
                        (
                            "left_expected_genotype"
                            if side == "Left"
                            else "right_expected_genotype"
                        )
                    ]
                )

                long_rows.append(
                    {
                        "site_id":
                            site["site_id"],
                        "module":
                            module,
                        "module_display":
                            site["module_display"],
                        "block":
                            site["block"],
                        "chrom":
                            chrom,
                        "pos":
                            pos,
                        "ref":
                            ref,
                        "alt":
                            alt,
                        "minimum_comparative_BA":
                            site[
                                "minimum_comparative_BA"
                            ],
                        "presence_associated_allele":
                            presence,
                        "absence_associated_allele":
                            absence,
                        "biological_orientation":
                            polarity[
                                "biological_orientation"
                            ],
                        "side":
                            side,
                        "GT_raw":
                            "missing",
                        "AD_raw":
                            "missing",
                        "DP":
                            np.nan,
                        "GQ":
                            np.nan,
                        "QC_pass_DP2_GQ20":
                            False,
                        "biological_genotype":
                            "missing",
                        "presence_dosage":
                            np.nan,
                        "presence_allele_fraction":
                            np.nan,
                        "expected_genotype":
                            expected_gt,
                        "matches_expectation":
                            np.nan,
                        "focal_call_status":
                            reason,
                        "source_vcf":
                            str(
                                VCF_BY_MODULE[
                                    module
                                ]
                            ),
                    }
                )

            continue

        for side in [
            "Left",
            "Right",
        ]:
            expected_gt = (
                expectation_map[module][
                    (
                        "left_expected_genotype"
                        if side == "Left"
                        else "right_expected_genotype"
                    )
                ]
            )

            result = evaluate_side(
                module=module,
                side=side,
                raw_row=raw_row,
                presence_allele=presence,
                expected_gt=expected_gt,
            )

            row = {
                "site_id":
                    site["site_id"],
                "module":
                    module,
                "module_display":
                    site["module_display"],
                "block":
                    site["block"],
                "site_order_within_module":
                    site[
                        "site_order_within_module"
                    ],
                "chrom":
                    chrom,
                "pos":
                    pos,
                "ref":
                    ref,
                "alt":
                    alt,
                "minimum_comparative_BA":
                    site[
                        "minimum_comparative_BA"
                    ],
                "presence_associated_allele":
                    presence,
                "absence_associated_allele":
                    absence,
                "biological_orientation":
                    polarity[
                        "biological_orientation"
                    ],
                "focal_call_status":
                    "exact_site_match",
                "source_vcf":
                    str(
                        VCF_BY_MODULE[
                            module
                        ]
                    ),
            }

            row.update(
                result
            )

            long_rows.append(
                row
            )

    corrected_sites = pd.DataFrame(
        corrected_sites
    )

    corrected_sites.to_csv(
        OUT_SITES,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    long = pd.DataFrame(
        long_rows
    )

    long.to_csv(
        OUT_LONG,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    # --------------------------------------------------------
    # Wide table.
    # --------------------------------------------------------

    base_cols = [
        "site_id",
        "module",
        "module_display",
        "block",
        "chrom",
        "pos",
        "ref",
        "alt",
        "minimum_comparative_BA",
        "presence_associated_allele",
        "absence_associated_allele",
        "biological_orientation",
        "source_vcf",
    ]

    value_cols = [
        "GT_raw",
        "AD_raw",
        "DP",
        "GQ",
        "QC_pass_DP2_GQ20",
        "biological_genotype",
        "presence_dosage",
        "presence_allele_fraction",
        "expected_genotype",
        "matches_expectation",
        "focal_call_status",
    ]

    left = (
        long.loc[
            long["side"] == "Left",
            base_cols + value_cols,
        ]
        .copy()
        .rename(
            columns={
                col: f"left_{col}"
                for col in value_cols
            }
        )
    )

    right = (
        long.loc[
            long["side"] == "Right",
            base_cols + value_cols,
        ]
        .copy()
        .rename(
            columns={
                col: f"right_{col}"
                for col in value_cols
            }
        )
    )

    wide = left.merge(
        right,
        on=base_cols,
        how="outer",
        validate="one_to_one",
    )

    wide["both_callable"] = (
        wide[
            "left_QC_pass_DP2_GQ20"
        ].astype(bool)
        &
        wide[
            "right_QC_pass_DP2_GQ20"
        ].astype(bool)
    )

    wide[
        "pair_matches_expectation"
    ] = np.where(
        wide["both_callable"],
        (
            wide[
                "left_matches_expectation"
            ].astype(bool)
            &
            wide[
                "right_matches_expectation"
            ].astype(bool)
        ),
        np.nan,
    )

    wide.to_csv(
        OUT_WIDE,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    # --------------------------------------------------------
    # Summary.
    # --------------------------------------------------------

    summary_rows = []

    for module in MODULE_ORDER:
        d = wide.loc[
            wide["module"] == module
        ].copy()

        if d.empty:
            continue

        left_callable = (
            d[
                "left_QC_pass_DP2_GQ20"
            ].astype(bool)
        )

        right_callable = (
            d[
                "right_QC_pass_DP2_GQ20"
            ].astype(bool)
        )

        both_callable = (
            d["both_callable"].astype(bool)
        )

        pair_callable = d.loc[
            both_callable
        ]

        summary_rows.append(
            {
                "module":
                    module,

                "frozen_comparative_sites":
                    len(d),

                "left_callable":
                    int(
                        left_callable.sum()
                    ),

                "right_callable":
                    int(
                        right_callable.sum()
                    ),

                "both_callable":
                    int(
                        both_callable.sum()
                    ),

                "left_exact_expected_n":
                    int(
                        d.loc[
                            left_callable,
                            "left_matches_expectation",
                        ]
                        .fillna(False)
                        .astype(bool)
                        .sum()
                    ),

                "right_exact_expected_n":
                    int(
                        d.loc[
                            right_callable,
                            "right_matches_expectation",
                        ]
                        .fillna(False)
                        .astype(bool)
                        .sum()
                    ),

                "pair_exact_expected_n":
                    int(
                        pair_callable[
                            "pair_matches_expectation"
                        ]
                        .fillna(False)
                        .astype(bool)
                        .sum()
                    ),

                "pair_exact_expected_fraction":
                    (
                        pair_callable[
                            "pair_matches_expectation"
                        ]
                        .astype(float)
                        .mean()
                        if len(
                            pair_callable
                        )
                        else np.nan
                    ),

                "left_expected_genotype":
                    expectation_map[
                        module
                    ][
                        "left_expected_genotype"
                    ],

                "right_expected_genotype":
                    expectation_map[
                        module
                    ][
                        "right_expected_genotype"
                    ],

                "reference_state":
                    REFERENCE_STATE[
                        module
                    ],

                "presence_allele_is":
                    (
                        "REF"
                        if REFERENCE_STATE[
                            module
                        ] == "presence"
                        else "ALT"
                    ),

                "left_median_presence_fraction":
                    pd.to_numeric(
                        d[
                            "left_presence_allele_fraction"
                        ],
                        errors="coerce",
                    ).median(),

                "right_median_presence_fraction":
                    pd.to_numeric(
                        d[
                            "right_presence_allele_fraction"
                        ],
                        errors="coerce",
                    ).median(),
            }
        )

    summary = pd.DataFrame(
        summary_rows
    )

    summary.to_csv(
        OUT_SUMMARY,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    print(
        "\nLOCKED FROZEN SITE COUNTS"
    )

    print(
        corrected_sites.groupby(
            "module"
        )
        .size()
        .reindex(
            MODULE_ORDER,
            fill_value=0,
        )
        .to_string()
    )

    print(
        "\nREFERENCE-BASED BIOLOGICAL POLARITY"
    )

    for module in MODULE_ORDER:
        state = REFERENCE_STATE[
            module
        ]

        if state == "presence":
            polarity_text = (
                "REF=P; ALT=A"
            )
        else:
            polarity_text = (
                "REF=A; ALT=P"
            )

        exp = expectation_map[
            module
        ]

        print(
            f"{module:8s}  "
            f"{polarity_text:12s}  "
            f"L={exp['left_expected_genotype']}  "
            f"R={exp['right_expected_genotype']}"
        )

    print(
        "\nORIGINAL-502 FOCAL SUMMARY"
    )

    print(
        summary.to_string(
            index=False
        )
    )

    print(
        "\nSITE/ALLELE MISMATCHES:",
        len(mismatch_rows),
    )

    if mismatch_rows:
        print(
            pd.DataFrame(
                mismatch_rows
            ).to_string(
                index=False
            )
        )

    print(
        "\nOutputs:"
    )

    for path in [
        OUT_SITES,
        OUT_RAW,
        OUT_LONG,
        OUT_WIDE,
        OUT_SUMMARY,
    ]:
        print(
            f"  {path}"
        )


if __name__ == "__main__":
    main()
