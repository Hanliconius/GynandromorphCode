#!/usr/bin/env python3

from pathlib import Path
import math
import re

import numpy as np
import pandas as pd


# ============================================================
# PATHS
# ============================================================

PROJECT = Path(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
)

VCF = (
    PROJECT
    / "external"
    / "krista_genotype_plots_original"
    / "Genotype_Plots"
    / "Optix_Hmr18-938195-1239421_biallelic.vcf"
)

MANIFEST = (
    PROJECT
    / "analysis"
    / "B_western_postman_vs_western_Amazon_DR_reanalysis"
    / "checkpoints"
    / "B_western_retained_sample_manifest.tsv"
)

MIXTURE = (
    PROJECT
    / "analysis"
    / "original_502"
    / "reciprocal_L_R_mixture"
    / "optix_reciprocal_L_R_genotype_pairs.posterior95.tsv"
)

OUTDIR = (
    PROJECT
    / "analysis"
    / "B_vs_DR_repolarization"
    / "checkpoints"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True,
)

SCORES_OUT = (
    OUTDIR
    / "B_vs_DR_broad_optix_site_scores.tsv.gz"
)

FOCAL_OUT = (
    OUTDIR
    / "B_vs_DR_supporting_sites_adjusted_LR.tsv"
)

COHORT_OUT = (
    OUTDIR
    / "B_vs_DR_comparative_cohort.tsv"
)

SUMMARY_OUT = (
    OUTDIR
    / "B_vs_DR_repolarization_summary.tsv"
)

CLASS_SUMMARY_OUT = (
    OUTDIR
    / "B_vs_DR_site_class_summary.tsv"
)


# ============================================================
# ANALYSIS SETTINGS
# ============================================================

CHROMOSOME = "Hmr18"

REGION_START = 1_060_000
REGION_END = 1_240_000

DP_MIN = 2
GQ_MIN = 20

MINIMUM_OVERALL_CALL_RATE = 0.70
MINIMUM_STATE_CALL_RATE = 0.50
MINIMUM_CALLED_PER_CLASS = 20
MINIMUM_PREDICTED_CLASS = 10

# Diagnostic D/R coherence.
MINIMUM_DR_SUBGROUP_CALLED = 5

# Initial "supporting" threshold.
SUPPORT_BA = 0.80

FOCAL_POSTERIOR_MIN = 0.95

EXPECTED_COHORT_N = 163


# ============================================================
# HELPERS
# ============================================================

def clean(value):
    if value is None:
        return ""
    return " ".join(
        str(value)
        .replace("\n", " ")
        .replace("\r", " ")
        .split()
    )


def numeric(value):
    try:
        x = float(value)
    except (TypeError, ValueError):
        return np.nan

    if math.isfinite(x):
        return x

    return np.nan


def parse_binary(value):
    x = numeric(value)

    if not math.isfinite(x):
        return np.nan

    if abs(x) < 1e-8:
        return 0

    if abs(x - 1) < 1e-8:
        return 1

    return np.nan


def canonical_gt(value):
    value = clean(value).replace("|", "/")

    if value in {
        "",
        ".",
        "./.",
        ".|.",
    }:
        return None

    alleles = value.split("/")

    if len(alleles) != 2:
        return None

    if any(
        allele not in {"0", "1"}
        for allele in alleles
    ):
        return None

    return "/".join(
        sorted(
            alleles,
            key=int,
        )
    )


def parse_number(value):
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


def parse_call(
    format_ids,
    sample_value,
):
    values = sample_value.split(":")

    fields = dict(
        zip(
            format_ids,
            values,
        )
    )

    gt = canonical_gt(
        fields.get("GT")
    )

    dp = parse_number(
        fields.get("DP")
    )

    gq = parse_number(
        fields.get("GQ")
    )

    passes = (
        gt is not None
        and dp is not None
        and dp >= DP_MIN
        and gq is not None
        and gq >= GQ_MIN
    )

    if not passes:
        return {
            "GT": None,
            "alt_dosage": None,
        }

    alt_dosage = sum(
        int(x)
        for x in gt.split("/")
    )

    return {
        "GT": gt,
        "alt_dosage": alt_dosage,
    }


def safe_fraction(
    numerator,
    denominator,
):
    if denominator == 0:
        return np.nan

    return numerator / denominator


def mcc_from_counts(
    tp,
    tn,
    fp,
    fn,
):
    numerator = (
        tp * tn
        - fp * fn
    )

    denominator = (
        (tp + fp)
        * (tp + fn)
        * (tn + fp)
        * (tn + fn)
    )

    if denominator <= 0:
        return np.nan

    return numerator / math.sqrt(
        denominator
    )


def calculate_orientation_metrics(
    phenotypes,
    genotypes,
    orientation,
):
    """
    phenotype:
        1 = B phenotype
        0 = DR phenotype

    orientation:
        ALT_carrier_predicts_B
        REF_carrier_predicts_B
    """

    predictions = []

    for gt in genotypes:
        a0, a1 = (
            int(x)
            for x in gt.split("/")
        )

        if orientation == "ALT_carrier_predicts_B":
            prediction = int(
                a0 == 1
                or a1 == 1
            )

        elif orientation == "REF_carrier_predicts_B":
            prediction = int(
                a0 == 0
                or a1 == 0
            )

        else:
            raise ValueError(
                orientation
            )

        predictions.append(
            prediction
        )

    tp = sum(
        phenotype == 1
        and prediction == 1
        for phenotype, prediction
        in zip(
            phenotypes,
            predictions,
        )
    )

    tn = sum(
        phenotype == 0
        and prediction == 0
        for phenotype, prediction
        in zip(
            phenotypes,
            predictions,
        )
    )

    fp = sum(
        phenotype == 0
        and prediction == 1
        for phenotype, prediction
        in zip(
            phenotypes,
            predictions,
        )
    )

    fn = sum(
        phenotype == 1
        and prediction == 0
        for phenotype, prediction
        in zip(
            phenotypes,
            predictions,
        )
    )

    sensitivity = safe_fraction(
        tp,
        tp + fn,
    )

    specificity = safe_fraction(
        tn,
        tn + fp,
    )

    if (
        math.isfinite(sensitivity)
        and math.isfinite(specificity)
    ):
        balanced_accuracy = (
            sensitivity
            + specificity
        ) / 2
    else:
        balanced_accuracy = np.nan

    return {
        "orientation": orientation,
        "tp": tp,
        "tn": tn,
        "fp": fp,
        "fn": fn,
        "misclassified_n": (
            fp + fn
        ),
        "sensitivity": sensitivity,
        "specificity": specificity,
        "balanced_accuracy": balanced_accuracy,
        "MCC": mcc_from_counts(
            tp,
            tn,
            fp,
            fn,
        ),
        "predicted_B": sum(
            predictions
        ),
        "predicted_DR": (
            len(predictions)
            - sum(predictions)
        ),
    }


def choose_orientation(
    phenotypes,
    genotypes,
):
    candidates = [
        calculate_orientation_metrics(
            phenotypes,
            genotypes,
            "ALT_carrier_predicts_B",
        ),
        calculate_orientation_metrics(
            phenotypes,
            genotypes,
            "REF_carrier_predicts_B",
        ),
    ]

    def rank(metric):
        ba = metric[
            "balanced_accuracy"
        ]

        mcc = metric[
            "MCC"
        ]

        if not math.isfinite(ba):
            ba = -np.inf

        if not math.isfinite(mcc):
            mcc = -np.inf

        # Deterministic only after BA/MCC/error ties.
        orientation_tiebreak = (
            1
            if metric["orientation"]
            == "ALT_carrier_predicts_B"
            else 0
        )

        return (
            ba,
            mcc,
            -metric[
                "misclassified_n"
            ],
            orientation_tiebreak,
        )

    ordered = sorted(
        candidates,
        key=rank,
        reverse=True,
    )

    best = ordered[0].copy()
    second = ordered[1]

    best[
        "orientation_BA_margin"
    ] = (
        best[
            "balanced_accuracy"
        ]
        - second[
            "balanced_accuracy"
        ]
        if (
            math.isfinite(
                best[
                    "balanced_accuracy"
                ]
            )
            and math.isfinite(
                second[
                    "balanced_accuracy"
                ]
            )
        )
        else np.nan
    )

    best[
        "orientation_tied_BA"
    ] = (
        math.isfinite(
            best[
                "balanced_accuracy"
            ]
        )
        and math.isfinite(
            second[
                "balanced_accuracy"
            ]
        )
        and abs(
            best[
                "balanced_accuracy"
            ]
            - second[
                "balanced_accuracy"
            ]
        ) < 1e-12
    )

    return best


def b_like_dosage(
    alt_dosage,
    orientation,
):
    if alt_dosage is None:
        return None

    if orientation == "ALT_carrier_predicts_B":
        return alt_dosage

    if orientation == "REF_carrier_predicts_B":
        return 2 - alt_dosage

    return None


def b_like_state_from_raw_gt(
    gt,
    orientation,
):
    gt = canonical_gt(
        gt
    )

    if gt is None:
        return "No call"

    alt_dosage = sum(
        int(x)
        for x in gt.split("/")
    )

    dosage = b_like_dosage(
        alt_dosage,
        orientation,
    )

    if dosage == 2:
        return "B/B"

    if dosage == 1:
        return "Het"

    if dosage == 0:
        return "DR/DR"

    return "No call"


# ============================================================
# READ AND VALIDATE COHORT
# ============================================================

for required_file in (
    VCF,
    MANIFEST,
    MIXTURE,
):
    if not required_file.is_file():
        raise RuntimeError(
            f"Missing required input: "
            f"{required_file}"
        )


manifest = pd.read_csv(
    MANIFEST,
    sep="\t",
    low_memory=False,
)

required_manifest_columns = {
    "cohort",
    "sample",
    "taxon_name",
    "country",
    "data_type",
    "B",
    "D",
    "R",
}

missing = (
    required_manifest_columns
    - set(
        manifest.columns
    )
)

if missing:
    raise RuntimeError(
        "Manifest missing columns: "
        + ", ".join(
            sorted(missing)
        )
    )


manifest = manifest.loc[
    manifest[
        "cohort"
    ].astype(str).eq(
        "WGS_plus_CS"
    )
].copy()


for trait in (
    "B",
    "D",
    "R",
):
    manifest[trait] = (
        manifest[trait]
        .map(
            parse_binary
        )
    )


manifest[
    "sample"
] = (
    manifest[
        "sample"
    ]
    .astype(str)
    .map(clean)
)


if (
    manifest[
        "sample"
    ].nunique()
    != EXPECTED_COHORT_N
):
    raise RuntimeError(
        "Expected "
        f"{EXPECTED_COHORT_N} "
        "unique WGS_plus_CS "
        "comparative butterflies; "
        f"found "
        f"{manifest['sample'].nunique()}."
    )


if manifest[
    "sample"
].duplicated().any():
    raise RuntimeError(
        "Duplicate sample rows in "
        "WGS_plus_CS manifest."
    )


manifest[
    "phenotype_class"
] = np.where(
    manifest[
        "B"
    ].eq(1),
    "B",
    np.where(
        manifest[
            "B"
        ].eq(0)
        & (
            manifest[
                "D"
            ].eq(1)
            |
            manifest[
                "R"
            ].eq(1)
        ),
        "DR",
        "EXCLUDE",
    ),
)


bad_B = manifest.loc[
    manifest[
        "phenotype_class"
    ].eq(
        "B"
    )
    & (
        manifest[
            "D"
        ].eq(1)
        |
        manifest[
            "R"
        ].eq(1)
    )
]

if len(bad_B):
    raise RuntimeError(
        "B-class butterflies "
        "with D/R-positive phenotype:\n"
        + bad_B[
            [
                "sample",
                "taxon_name",
                "B",
                "D",
                "R",
            ]
        ].to_string(
            index=False
        )
    )


excluded = manifest.loc[
    manifest[
        "phenotype_class"
    ].eq(
        "EXCLUDE"
    )
]

if len(excluded):
    raise RuntimeError(
        "Final 142 WGS+CS cohort "
        "contains butterflies that are "
        "neither clean B nor clean DR:\n"
        + excluded[
            [
                "sample",
                "taxon_name",
                "B",
                "D",
                "R",
            ]
        ].to_string(
            index=False
        )
    )


manifest[
    "phenotype"
] = (
    manifest[
        "phenotype_class"
    ]
    .map(
        {
            "DR": 0,
            "B": 1,
        }
    )
    .astype(int)
)


manifest[
    "D_subgroup"
] = (
    manifest[
        "phenotype_class"
    ].eq(
        "DR"
    )
    & manifest[
        "D"
    ].eq(1)
)

manifest[
    "R_subgroup"
] = (
    manifest[
        "phenotype_class"
    ].eq(
        "DR"
    )
    & manifest[
        "R"
    ].eq(1)
)


manifest.to_csv(
    COHORT_OUT,
    sep="\t",
    index=False,
)


print(
    "\n========================================"
)
print(
    "UNIFIED B-vs-DR COMPARATIVE COHORT"
)
print(
    "========================================\n"
)

print(
    manifest.groupby(
        [
            "phenotype_class",
            "B",
            "D",
            "R",
        ],
        dropna=False,
    )
    .size()
    .rename(
        "butterflies"
    )
    .reset_index()
    .to_string(
        index=False
    )
)

print()

print(
    "B butterflies:",
    int(
        manifest[
            "phenotype_class"
        ].eq(
            "B"
        ).sum()
    ),
)

print(
    "DR butterflies:",
    int(
        manifest[
            "phenotype_class"
        ].eq(
            "DR"
        ).sum()
    ),
)

print(
    "D-positive within DR:",
    int(
        manifest[
            "D_subgroup"
        ].sum()
    ),
)

print(
    "R-positive within DR:",
    int(
        manifest[
            "R_subgroup"
        ].sum()
    ),
)


# ============================================================
# PREPARE LOOKUPS
# ============================================================

cohort_samples = (
    manifest[
        "sample"
    ]
    .tolist()
)

phenotype_by_sample = dict(
    zip(
        manifest[
            "sample"
        ],
        manifest[
            "phenotype"
        ],
    )
)

class_by_sample = dict(
    zip(
        manifest[
            "sample"
        ],
        manifest[
            "phenotype_class"
        ],
    )
)

D_by_sample = dict(
    zip(
        manifest[
            "sample"
        ],
        manifest[
            "D_subgroup"
        ],
    )
)

R_by_sample = dict(
    zip(
        manifest[
            "sample"
        ],
        manifest[
            "R_subgroup"
        ],
    )
)


n_B_total = int(
    manifest[
        "phenotype_class"
    ].eq(
        "B"
    ).sum()
)

n_DR_total = int(
    manifest[
        "phenotype_class"
    ].eq(
        "DR"
    ).sum()
)


# ============================================================
# SCAN EXISTING BROAD VCF
#
# No VCF extraction or regeneration.
# ============================================================

rows = []

sample_indices = None
header_samples = None


with VCF.open(
    "rt"
) as handle:

    for line in handle:

        if line.startswith(
            "##"
        ):
            continue

        if line.startswith(
            "#CHROM"
        ):
            header = (
                line.rstrip(
                    "\n"
                )
                .split("\t")
            )

            header_samples = (
                header[9:]
            )

            header_index = {
                sample: index
                for index, sample
                in enumerate(
                    header_samples
                )
            }

            missing_samples = [
                sample
                for sample
                in cohort_samples
                if sample
                not in header_index
            ]

            if missing_samples:
                raise RuntimeError(
                    "Comparative samples absent "
                    "from broad VCF:\n"
                    + "\n".join(
                        missing_samples
                    )
                )

            sample_indices = {
                sample: header_index[
                    sample
                ]
                for sample
                in cohort_samples
            }

            continue


        if sample_indices is None:
            continue


        fields = (
            line.rstrip(
                "\n"
            )
            .split("\t")
        )

        if len(fields) < 10:
            continue

        chrom = fields[0]

        try:
            pos = int(
                fields[1]
            )
        except ValueError:
            continue

        if chrom != CHROMOSOME:
            continue

        if (
            pos < REGION_START
            or pos > REGION_END
        ):
            continue


        ref = fields[3]
        alt = fields[4]

        if (
            "," in alt
            or len(ref) != 1
            or len(alt) != 1
        ):
            continue


        format_ids = (
            fields[8]
            .split(":")
        )

        sample_values = fields[9:]


        phenotypes = []
        genotypes = []

        called_by_class = {
            "B": 0,
            "DR": 0,
        }

        missing_by_class = {
            "B": 0,
            "DR": 0,
        }

        b_like_inputs = []

        sample_calls = {}


        for sample in cohort_samples:

            call = parse_call(
                format_ids,
                sample_values[
                    sample_indices[
                        sample
                    ]
                ],
            )

            sample_calls[
                sample
            ] = call

            phenotype_class = (
                class_by_sample[
                    sample
                ]
            )

            if call[
                "GT"
            ] is None:

                missing_by_class[
                    phenotype_class
                ] += 1

            else:

                called_by_class[
                    phenotype_class
                ] += 1

                phenotypes.append(
                    phenotype_by_sample[
                        sample
                    ]
                )

                genotypes.append(
                    call[
                        "GT"
                    ]
                )


        called_n = (
            len(
                genotypes
            )
        )

        call_rate = (
            called_n
            / EXPECTED_COHORT_N
        )

        B_call_rate = safe_fraction(
            called_by_class[
                "B"
            ],
            n_B_total,
        )

        DR_call_rate = safe_fraction(
            called_by_class[
                "DR"
            ],
            n_DR_total,
        )


        orientation = (
            choose_orientation(
                phenotypes,
                genotypes,
            )
        )


        callability_eligible = (
            call_rate
            >= MINIMUM_OVERALL_CALL_RATE
            and B_call_rate
            >= MINIMUM_STATE_CALL_RATE
            and DR_call_rate
            >= MINIMUM_STATE_CALL_RATE
            and called_by_class[
                "B"
            ]
            >= MINIMUM_CALLED_PER_CLASS
            and called_by_class[
                "DR"
            ]
            >= MINIMUM_CALLED_PER_CLASS
        )


        score_eligible = (
            callability_eligible
            and math.isfinite(
                orientation[
                    "MCC"
                ]
            )
            and orientation[
                "predicted_B"
            ]
            >= MINIMUM_PREDICTED_CLASS
            and orientation[
                "predicted_DR"
            ]
            >= MINIMUM_PREDICTED_CLASS
        )


        selected_orientation = (
            orientation[
                "orientation"
            ]
        )


        if (
            selected_orientation
            == "ALT_carrier_predicts_B"
        ):
            B_like_allele = alt
            DR_like_allele = ref

        elif (
            selected_orientation
            == "REF_carrier_predicts_B"
        ):
            B_like_allele = ref
            DR_like_allele = alt

        else:
            B_like_allele = ""
            DR_like_allele = ""


        allele_groups = {
            "B": [],
            "DR": [],
            "D": [],
            "R": [],
        }


        for sample in cohort_samples:

            call = sample_calls[
                sample
            ]

            dosage = b_like_dosage(
                call[
                    "alt_dosage"
                ],
                selected_orientation,
            )

            if dosage is None:
                continue

            allele_groups[
                class_by_sample[
                    sample
                ]
            ].append(
                dosage
            )

            if D_by_sample[
                sample
            ]:
                allele_groups[
                    "D"
                ].append(
                    dosage
                )

            if R_by_sample[
                sample
            ]:
                allele_groups[
                    "R"
                ].append(
                    dosage
                )


        def group_AF(group):
            values = allele_groups[
                group
            ]

            if not values:
                return np.nan

            return (
                sum(values)
                / (
                    2
                    * len(values)
                )
            )


        B_AF = group_AF(
            "B"
        )

        DR_AF = group_AF(
            "DR"
        )

        D_AF = group_AF(
            "D"
        )

        R_AF = group_AF(
            "R"
        )


        B_minus_DR = (
            B_AF - DR_AF
            if (
                math.isfinite(
                    B_AF
                )
                and math.isfinite(
                    DR_AF
                )
            )
            else np.nan
        )

        B_minus_D = (
            B_AF - D_AF
            if (
                math.isfinite(
                    B_AF
                )
                and math.isfinite(
                    D_AF
                )
            )
            else np.nan
        )

        B_minus_R = (
            B_AF - R_AF
            if (
                math.isfinite(
                    B_AF
                )
                and math.isfinite(
                    R_AF
                )
            )
            else np.nan
        )


        D_called = len(
            allele_groups[
                "D"
            ]
        )

        R_called = len(
            allele_groups[
                "R"
            ]
        )


        DR_coherence_testable = (
            D_called
            >= MINIMUM_DR_SUBGROUP_CALLED
            and R_called
            >= MINIMUM_DR_SUBGROUP_CALLED
            and math.isfinite(
                B_minus_D
            )
            and math.isfinite(
                B_minus_R
            )
        )


        DR_direction_coherent = (
            DR_coherence_testable
            and B_minus_D > 0
            and B_minus_R > 0
        )


        DR_strongly_coherent = (
            DR_direction_coherent
            and B_minus_D >= 0.20
            and B_minus_R >= 0.20
        )


        BDR_support = (
            score_eligible
            and math.isfinite(
                orientation[
                    "balanced_accuracy"
                ]
            )
            and orientation[
                "balanced_accuracy"
            ] >= SUPPORT_BA
            and DR_direction_coherent
        )


        if BDR_support:

            if (
                selected_orientation
                == "ALT_carrier_predicts_B"
            ):
                ALT_class = "B"

            elif (
                selected_orientation
                == "REF_carrier_predicts_B"
            ):
                ALT_class = "DR"

            else:
                ALT_class = "neither"

        else:
            ALT_class = "neither"


        rows.append(
            {
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "position_Mb": (
                    pos / 1e6
                ),
                "site_key": (
                    f"{chrom}:{pos}:"
                    f"{ref}:{alt}"
                ),

                "comparative_cohort_n": (
                    EXPECTED_COHORT_N
                ),

                "B_total": (
                    n_B_total
                ),
                "DR_total": (
                    n_DR_total
                ),

                "B_called": (
                    called_by_class[
                        "B"
                    ]
                ),
                "DR_called": (
                    called_by_class[
                        "DR"
                    ]
                ),

                "D_called": D_called,
                "R_called": R_called,

                "overall_call_rate": (
                    call_rate
                ),
                "B_call_rate": (
                    B_call_rate
                ),
                "DR_call_rate": (
                    DR_call_rate
                ),

                "selected_orientation": (
                    selected_orientation
                ),

                "B_like_allele": (
                    B_like_allele
                ),
                "DR_like_allele": (
                    DR_like_allele
                ),

                "ALT_class": ALT_class,

                "balanced_accuracy": (
                    orientation[
                        "balanced_accuracy"
                    ]
                ),
                "sensitivity_B": (
                    orientation[
                        "sensitivity"
                    ]
                ),
                "specificity_DR": (
                    orientation[
                        "specificity"
                    ]
                ),
                "MCC": (
                    orientation[
                        "MCC"
                    ]
                ),

                "tp": orientation[
                    "tp"
                ],
                "tn": orientation[
                    "tn"
                ],
                "fp": orientation[
                    "fp"
                ],
                "fn": orientation[
                    "fn"
                ],

                "misclassified_n": (
                    orientation[
                        "misclassified_n"
                    ]
                ),

                "predicted_B": (
                    orientation[
                        "predicted_B"
                    ]
                ),
                "predicted_DR": (
                    orientation[
                        "predicted_DR"
                    ]
                ),

                "orientation_BA_margin": (
                    orientation[
                        "orientation_BA_margin"
                    ]
                ),

                "orientation_tied_BA": (
                    orientation[
                        "orientation_tied_BA"
                    ]
                ),

                "callability_eligible": (
                    callability_eligible
                ),

                "score_eligible": (
                    score_eligible
                ),

                "B_like_AF_B": (
                    B_AF
                ),
                "B_like_AF_DR": (
                    DR_AF
                ),
                "B_like_AF_D": (
                    D_AF
                ),
                "B_like_AF_R": (
                    R_AF
                ),

                "delta_B_minus_DR": (
                    B_minus_DR
                ),
                "delta_B_minus_D": (
                    B_minus_D
                ),
                "delta_B_minus_R": (
                    B_minus_R
                ),

                "DR_coherence_testable": (
                    DR_coherence_testable
                ),

                "DR_direction_coherent": (
                    DR_direction_coherent
                ),

                "DR_strongly_coherent": (
                    DR_strongly_coherent
                ),

                "BDR_support": (
                    BDR_support
                ),
            }
        )


scores = pd.DataFrame(
    rows
)


if scores.empty:
    raise RuntimeError(
        "No biallelic broad-region "
        "SNPs were scored."
    )


if scores[
    [
        "chrom",
        "pos",
        "ref",
        "alt",
    ]
].duplicated().any():
    raise RuntimeError(
        "Duplicate scored variant keys."
    )


scores.to_csv(
    SCORES_OUT,
    sep="\t",
    index=False,
    compression="gzip",
)


# ============================================================
# JOIN EXISTING RECIPROCAL-MIXTURE FOCAL CALLS
# ============================================================

mixture = pd.read_csv(
    MIXTURE,
    sep="\t",
    low_memory=False,
)


required_mixture_columns = {
    "chrom",
    "pos",
    "ref",
    "alt",
    "mixture_left_GT",
    "mixture_right_GT",
    "pair_posterior",
}

missing = (
    required_mixture_columns
    - set(
        mixture.columns
    )
)

if missing:
    raise RuntimeError(
        "Mixture table missing columns: "
        + ", ".join(
            sorted(missing)
        )
    )


mixture[
    "pos"
] = pd.to_numeric(
    mixture[
        "pos"
    ],
    errors="coerce",
)


mixture = mixture.loc[
    mixture[
        "chrom"
    ].astype(str).eq(
        CHROMOSOME
    )
    & mixture[
        "pos"
    ].between(
        REGION_START,
        REGION_END,
    )
].copy()


if mixture[
    [
        "chrom",
        "pos",
        "ref",
        "alt",
    ]
].duplicated().any():
    raise RuntimeError(
        "Duplicate mixture variant keys."
    )


support = scores.loc[
    scores[
        "BDR_support"
    ]
].copy()


focal = support.merge(
    mixture[
        [
            "chrom",
            "pos",
            "ref",
            "alt",
            "mixture_left_GT",
            "mixture_right_GT",
            "pair_posterior",
        ]
    ],
    on=[
        "chrom",
        "pos",
        "ref",
        "alt",
    ],
    how="left",
    validate="one_to_one",
)


focal[
    "pair_posterior"
] = pd.to_numeric(
    focal[
        "pair_posterior"
    ],
    errors="coerce",
)


focal[
    "focal_pair_confident"
] = (
    focal[
        "pair_posterior"
    ].ge(
        FOCAL_POSTERIOR_MIN
    )
)


focal[
    "TN_Left_raw_GT"
] = np.where(
    focal[
        "focal_pair_confident"
    ],
    focal[
        "mixture_left_GT"
    ],
    "./.",
)


focal[
    "TN_Right_raw_GT"
] = np.where(
    focal[
        "focal_pair_confident"
    ],
    focal[
        "mixture_right_GT"
    ],
    "./.",
)


focal[
    "TN_Left_BDR_GT"
] = [
    b_like_state_from_raw_gt(
        gt,
        orientation,
    )
    for gt, orientation
    in zip(
        focal[
            "TN_Left_raw_GT"
        ],
        focal[
            "selected_orientation"
        ],
    )
]


focal[
    "TN_Right_BDR_GT"
] = [
    b_like_state_from_raw_gt(
        gt,
        orientation,
    )
    for gt, orientation
    in zip(
        focal[
            "TN_Right_raw_GT"
        ],
        focal[
            "selected_orientation"
        ],
    )
]


focal[
    "focal_state_pair"
] = (
    focal[
        "TN_Left_BDR_GT"
    ]
    + " -> "
    + focal[
        "TN_Right_BDR_GT"
    ]
)


focal.to_csv(
    FOCAL_OUT,
    sep="\t",
    index=False,
)


# ============================================================
# SUMMARIES
# ============================================================

classification_summary = (
    scores.groupby(
        [
            "ALT_class",
        ],
        dropna=False,
    )
    .agg(
        SNPs=(
            "site_key",
            "size",
        ),
        eligible=(
            "score_eligible",
            "sum",
        ),
        BDR_support=(
            "BDR_support",
            "sum",
        ),
        mean_BA=(
            "balanced_accuracy",
            "mean",
        ),
        max_BA=(
            "balanced_accuracy",
            "max",
        ),
    )
    .reset_index()
)


classification_summary.to_csv(
    CLASS_SUMMARY_OUT,
    sep="\t",
    index=False,
)


summary = pd.DataFrame(
    [
        {
            "statistic": (
                "comparative_butterflies"
            ),
            "value": (
                EXPECTED_COHORT_N
            ),
        },
        {
            "statistic": (
                "B_butterflies"
            ),
            "value": (
                n_B_total
            ),
        },
        {
            "statistic": (
                "DR_butterflies"
            ),
            "value": (
                n_DR_total
            ),
        },
        {
            "statistic": (
                "scored_biallelic_SNPs"
            ),
            "value": (
                len(scores)
            ),
        },
        {
            "statistic": (
                "score_eligible_SNPs"
            ),
            "value": int(
                scores[
                    "score_eligible"
                ].sum()
            ),
        },
        {
            "statistic": (
                "DR_coherence_testable"
            ),
            "value": int(
                scores[
                    "DR_coherence_testable"
                ].sum()
            ),
        },
        {
            "statistic": (
                "DR_direction_coherent"
            ),
            "value": int(
                scores[
                    "DR_direction_coherent"
                ].sum()
            ),
        },
        {
            "statistic": (
                "BA_ge_0.80_eligible"
            ),
            "value": int(
                (
                    scores[
                        "score_eligible"
                    ]
                    & scores[
                        "balanced_accuracy"
                    ].ge(
                        0.80
                    )
                ).sum()
            ),
        },
        {
            "statistic": (
                "B_DR_supporting_SNPs"
            ),
            "value": int(
                scores[
                    "BDR_support"
                ].sum()
            ),
        },
        {
            "statistic": (
                "ALT_B_like_supporting"
            ),
            "value": int(
                scores[
                    "ALT_class"
                ].eq(
                    "B"
                ).sum()
            ),
        },
        {
            "statistic": (
                "ALT_DR_like_supporting"
            ),
            "value": int(
                scores[
                    "ALT_class"
                ].eq(
                    "DR"
                ).sum()
            ),
        },
        {
            "statistic": (
                "supporting_with_confident_focal_pair"
            ),
            "value": int(
                focal[
                    "focal_pair_confident"
                ].sum()
            ),
        },
    ]
)


summary.to_csv(
    SUMMARY_OUT,
    sep="\t",
    index=False,
)


pair_summary = (
    focal.loc[
        focal[
            "focal_pair_confident"
        ]
    ]
    .groupby(
        [
            "focal_state_pair",
        ],
        dropna=False,
    )
    .size()
    .rename(
        "SNPs"
    )
    .reset_index()
    .sort_values(
        "SNPs",
        ascending=False,
    )
)


print(
    "\n========================================"
)
print(
    "UNIFIED B-vs-DR SITE SUMMARY"
)
print(
    "========================================\n"
)

print(
    summary.to_string(
        index=False
    )
)

print(
    "\nALT-ALLELE CLASSIFICATION\n"
)

print(
    classification_summary.to_string(
        index=False
    )
)

print(
    "\nCONFIDENT FOCAL B/DR GENOTYPE PAIRS\n"
)

print(
    pair_summary.to_string(
        index=False
    )
)

print(
    "\nOutputs:"
)

for path in (
    SCORES_OUT,
    FOCAL_OUT,
    COHORT_OUT,
    SUMMARY_OUT,
    CLASS_SUMMARY_OUT,
):
    print(
        path
    )
