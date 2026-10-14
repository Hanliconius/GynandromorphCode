#!/usr/bin/env python3

from pathlib import Path
import importlib.util
import pandas as pd


PROJECT = Path(
    "/scratch/martinlab/joe/GynandroCalling"
)

ORIGINAL_SCORER = (
    PROJECT
    / "scripts"
    / "41_score_original_502_optix_callsets.py"
)

SCORES = (
    PROJECT
    / "analysis"
    / "original_502"
    / "checkpoints"
    / "original502_allSNP_scores.filtered.WGS_vs_WGSplusCS.DP2_GQ20.tsv"
)

SCORE_MANIFEST = (
    PROJECT
    / "analysis"
    / "original_502"
    / "checkpoints"
    / "original502_allSNP_scores.filtered.WGS_vs_WGSplusCS.score_manifest.tsv"
)

OUTDIR = (
    PROJECT
    / "analysis"
    / "original_502"
    / "checkpoints"
)

OUT_SITES = (
    OUTDIR
    / "original502_recalculated_diagnostic_sites.WGS_vs_WGSplusCS.tsv"
)

OUT_LONG = (
    OUTDIR
    / "original502_recalculated_diagnostic_genotypes.WGS_vs_WGSplusCS.long.tsv"
)

OUT_SUMMARY = (
    OUTDIR
    / "original502_recalculated_diagnostic_genotypes.WGS_vs_WGSplusCS.summary.tsv"
)


MODULES = [
    "Band",
    "Ray2",
    "Dennis2",
]

COHORTS = [
    "WGS_only",
    "WGS_plus_CS",
]

PANEL_FROM_MODULE = {
    "Band": "band",
    "Ray2": "ray2",
    "Dennis2": "dennis2",
}

BA_THRESHOLD = 0.80
SCHEME = "DP2_GQ20"


def fail(message):
    raise RuntimeError(message)


def as_bool(value):

    return str(value).strip().lower() in {
        "true",
        "t",
        "1",
        "yes",
        "y",
    }


def load_original():

    spec = importlib.util.spec_from_file_location(
        "original502_scorer",
        ORIGINAL_SCORER,
    )

    if spec is None or spec.loader is None:
        fail(
            f"Could not import {ORIGINAL_SCORER}"
        )

    module = importlib.util.module_from_spec(
        spec
    )

    spec.loader.exec_module(
        module
    )

    return module


def biological_genotype(
    raw_gt,
    module,
):
    """
    Translate raw REF/ALT GT into phenotype-presence P/A space.

    Band:
        REF = P
        ALT = A

    Ray2 / Dennis2:
        REF = A
        ALT = P
    """

    if raw_gt is None:
        return "missing"

    gt = str(
        raw_gt
    ).replace(
        "|",
        "/",
    )

    alleles = gt.split(
        "/"
    )

    if len(alleles) != 2:
        return "missing"

    if any(
        allele not in {
            "0",
            "1",
        }
        for allele in alleles
    ):
        return "missing"

    alt_dosage = sum(
        allele == "1"
        for allele in alleles
    )

    if module == "Band":

        mapping = {
            0: "PP",
            1: "PA",
            2: "AA",
        }

    else:

        mapping = {
            0: "AA",
            1: "PA",
            2: "PP",
        }

    return mapping[
        alt_dosage
    ]


def select_sites(
    scores,
):

    scores[
        "score_eligible_bool"
    ] = scores[
        "score_eligible"
    ].map(
        as_bool
    )

    scores[
        "matches_reference_expectation_bool"
    ] = scores[
        "matches_reference_expectation"
    ].map(
        as_bool
    )

    scores[
        "balanced_accuracy"
    ] = pd.to_numeric(
        scores[
            "balanced_accuracy"
        ],
        errors="coerce",
    )

    selected = []

    for cohort in COHORTS:

        for module in MODULES:

            d = scores.loc[
                (
                    scores[
                        "analysis_cohort"
                    ]
                    == cohort
                )
                &
                (
                    scores[
                        "module"
                    ]
                    == module
                )
            ].copy()

            eligible = d.loc[
                d[
                    "score_eligible_bool"
                ]
                &
                d[
                    "matches_reference_expectation_bool"
                ]
            ].copy()

            passing = eligible.loc[
                eligible[
                    "balanced_accuracy"
                ]
                >= BA_THRESHOLD
            ].copy()

            if (
                cohort == "WGS_plus_CS"
                and module == "Band"
                and len(passing) == 0
            ):

                use = eligible.copy()

                selection_rule = (
                    "all_expected_pattern_eligible_"
                    "because_none_reached_BA0.80"
                )

            else:

                use = passing.copy()

                selection_rule = (
                    "expected_pattern_eligible_"
                    "and_BA_ge_0.80"
                )

            use[
                "selection_rule"
            ] = selection_rule

            use[
                "passes_BA_0.80"
            ] = (
                use[
                    "balanced_accuracy"
                ]
                >= BA_THRESHOLD
            )

            use = use.sort_values(
                [
                    "chrom",
                    "pos",
                ]
            ).copy()

            use[
                "site_order"
            ] = range(
                1,
                len(use) + 1,
            )

            selected.append(
                use
            )

            print(
                f"{cohort:12s} / {module:8s}: "
                f"eligible={len(eligible):3d}, "
                f"BA>=0.80={len(passing):3d}, "
                f"plotted={len(use):3d}"
            )

    return pd.concat(
        selected,
        ignore_index=True,
    )


def build_included_samples(
    manifest,
):

    manifest[
        "reason"
    ] = manifest[
        "reason"
    ].astype(str)

    included = manifest.loc[
        manifest[
            "reason"
        ]
        == "included"
    ].copy()

    # One row per cohort/module/sample.
    included = included.drop_duplicates(
        subset=[
            "analysis_cohort",
            "module",
            "sample",
        ]
    )

    return included


def parse_selected_vcf(
    original,
    vcf_path,
    module,
    cohort,
    selected_sites,
    included_samples,
):

    wanted_by_key = {}

    for _, row in selected_sites.iterrows():

        key = (
            str(
                row[
                    "chrom"
                ]
            ),
            int(
                row[
                    "pos"
                ]
            ),
            str(
                row[
                    "ref"
                ]
            ),
            str(
                row[
                    "alt"
                ]
            ),
        )

        wanted_by_key[
            key
        ] = row.to_dict()

    sample_meta = {}

    for _, row in included_samples.iterrows():

        sample_meta[
            row[
                "sample"
            ]
        ] = {
            "phenotype":
                row[
                    "phenotype"
                ],

            "taxon_name":
                row[
                    "taxon"
                ],

            "data_type":
                row[
                    "data_type"
                ],
        }

    rows = []

    with Path(
        vcf_path
    ).open(
        errors="replace"
    ) as handle:

        sample_names = []
        sample_indices = {}

        for line in handle:

            if line.startswith(
                "##"
            ):
                continue

            if line.startswith(
                "#CHROM"
            ):

                fields = line.rstrip(
                    "\n"
                ).split(
                    "\t"
                )

                sample_names = fields[
                    9:
                ]

                sample_indices = {
                    sample:
                        index

                    for index, sample
                    in enumerate(
                        sample_names
                    )
                }

                continue

            if line.startswith(
                "#"
            ):
                continue

            fields = line.rstrip(
                "\n"
            ).split(
                "\t"
            )

            if len(fields) < 10:
                continue

            chrom = fields[0]
            pos = int(
                fields[1]
            )
            ref = fields[3]
            alt = fields[4]

            key = (
                chrom,
                pos,
                ref,
                alt,
            )

            if key not in wanted_by_key:
                continue

            site = wanted_by_key[
                key
            ]

            format_ids = fields[
                8
            ].split(
                ":"
            )

            sample_values = fields[
                9:
            ]

            for sample, meta in sample_meta.items():

                index = sample_indices.get(
                    sample
                )

                if index is None:

                    call = {
                        "GT": None,
                        "DP": None,
                        "GQ": None,
                        "AD": None,
                    }

                    passes_qc = False

                else:

                    call = original.parse_call(
                        format_ids,
                        sample_values[
                            index
                        ],
                    )

                    passes_qc = original.call_passes_scheme(
                        call,
                        original.SCHEMES[
                            SCHEME
                        ],
                    )

                if passes_qc:

                    bio_gt = biological_genotype(
                        call[
                            "GT"
                        ],
                        module,
                    )

                else:

                    bio_gt = "missing"

                rows.append(
                    {
                        "cohort":
                            cohort,

                        "module":
                            module,

                        "sample":
                            sample,

                        "taxon_name":
                            meta[
                                "taxon_name"
                            ],

                        "data_type":
                            meta[
                                "data_type"
                            ],

                        "phenotype":
                            meta[
                                "phenotype"
                            ],

                        "chrom":
                            chrom,

                        "pos":
                            pos,

                        "ref":
                            ref,

                        "alt":
                            alt,

                        "site_order":
                            site[
                                "site_order"
                            ],

                        "balanced_accuracy":
                            site[
                                "balanced_accuracy"
                            ],

                        "passes_BA_0.80":
                            site[
                                "passes_BA_0.80"
                            ],

                        "selection_rule":
                            site[
                                "selection_rule"
                            ],

                        "raw_GT":
                            (
                                ""
                                if call[
                                    "GT"
                                ]
                                is None
                                else call[
                                    "GT"
                                ]
                            ),

                        "DP":
                            call.get(
                                "DP"
                            ),

                        "GQ":
                            call.get(
                                "GQ"
                            ),

                        "biological_genotype":
                            bio_gt,

                        "passes_DP2_GQ20":
                            passes_qc,
                    }
                )

    return rows


def main():

    original = load_original()

    if not SCORES.exists():
        fail(
            f"Missing scores: {SCORES}"
        )

    if not SCORE_MANIFEST.exists():
        fail(
            f"Missing score manifest: {SCORE_MANIFEST}"
        )

    scores = pd.read_csv(
        SCORES,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    manifest = pd.read_csv(
        SCORE_MANIFEST,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    selected = select_sites(
        scores
    )

    selected.to_csv(
        OUT_SITES,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    included = build_included_samples(
        manifest
    )

    all_rows = []

    for cohort in COHORTS:

        for module in MODULES:

            panel = PANEL_FROM_MODULE[
                module
            ]

            specification = original.PANELS[
                panel
            ]

            site_subset = selected.loc[
                (
                    selected[
                        "analysis_cohort"
                    ]
                    == cohort
                )
                &
                (
                    selected[
                        "module"
                    ]
                    == module
                )
            ].copy()

            sample_subset = included.loc[
                (
                    included[
                        "analysis_cohort"
                    ]
                    == cohort
                )
                &
                (
                    included[
                        "module"
                    ]
                    == module
                )
            ].copy()

            print(
                f"\nBuilding {cohort} / {module}: "
                f"{len(site_subset)} sites x "
                f"{len(sample_subset)} samples"
            )

            rows = parse_selected_vcf(
                original=original,
                vcf_path=specification[
                    "vcf"
                ],
                module=module,
                cohort=cohort,
                selected_sites=site_subset,
                included_samples=sample_subset,
            )

            all_rows.extend(
                rows
            )

    long = pd.DataFrame(
        all_rows
    )

    if long.empty:
        fail(
            "No genotype rows generated."
        )

    long.to_csv(
        OUT_LONG,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    summary = (
        long.groupby(
            [
                "cohort",
                "module",
            ],
            as_index=False,
        )
        .agg(
            sample_n=(
                "sample",
                "nunique",
            ),
            site_n=(
                "pos",
                "nunique",
            ),
            genotype_rows=(
                "biological_genotype",
                "size",
            ),
            callable_rows=(
                "passes_DP2_GQ20",
                "sum",
            ),
        )
    )

    summary[
        "callable_fraction"
    ] = (
        summary[
            "callable_rows"
        ]
        /
        summary[
            "genotype_rows"
        ]
    )

    summary.to_csv(
        OUT_SUMMARY,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    print(
        "\n=============================================="
    )

    print(
        "RECALCULATED GENOTYPE TABLE SUMMARY"
    )

    print(
        "=============================================="
    )

    print(
        summary.to_string(
            index=False
        )
    )

    print(
        "\nOutputs:"
    )

    for path in [
        OUT_SITES,
        OUT_LONG,
        OUT_SUMMARY,
    ]:
        print(
            f"  {path}"
        )


if __name__ == "__main__":
    main()
