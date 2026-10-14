#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd


PROJECT = Path(
    "/scratch/martinlab/joe/GynandroCalling"
)

CK = (
    PROJECT
    / "analysis"
    / "original_502"
    / "checkpoints"
)

LONG = (
    CK
    / "optix_original502_selected_sites_genotype_plot.long.tsv"
)

SITES = (
    CK
    / "optix_original502_selected_sites_for_plotting.tsv"
)

CALLABILITY = (
    CK
    / "original502_sample_callability_by_module.tsv"
)

META = (
    PROJECT
    / "data"
    / "derived"
    / "region_panel_metadata.drop_problematic.tsv"
)

OUT_SCORES = (
    CK
    / "original502_frozen_site_scores.WGS_vs_WGSplusCS.tsv"
)

OUT_COHORT = (
    CK
    / "original502_filtered_cohort_manifest.WGS_vs_WGSplusCS.tsv"
)

OUT_SUMMARY = (
    CK
    / "original502_frozen_site_scores.WGS_vs_WGSplusCS.summary.tsv"
)


MODULES = [
    "Band",
    "Ray2",
    "Dennis2",
]

CALLABILITY_CUTOFF = 0.25


def balanced_accuracy(obs, pred):

    obs = np.asarray(
        obs,
        dtype=int,
    )

    pred = np.asarray(
        pred,
        dtype=int,
    )

    if len(obs) == 0:
        return np.nan

    pos = obs == 1
    neg = obs == 0

    if (
        pos.sum() == 0
        or neg.sum() == 0
    ):
        return np.nan

    sensitivity = (
        pred[pos] == 1
    ).mean()

    specificity = (
        pred[neg] == 0
    ).mean()

    return (
        sensitivity
        + specificity
    ) / 2


def score_site(d):

    called = d.loc[
        d["biological_genotype"].isin(
            [
                "AA",
                "PA",
                "PP",
            ]
        )
    ].copy()

    if len(called) == 0:
        return {
            "cohort_n": len(d),
            "called_n": 0,
            "missing_n": len(d),
            "call_rate": 0.0,
            "phenotype_0_called": 0,
            "phenotype_1_called": 0,
            "match_n": 0,
            "mismatch_n": 0,
            "match_fraction": np.nan,
            "sensitivity": np.nan,
            "specificity": np.nan,
            "balanced_accuracy": np.nan,
        }

    # Biological genotype is already in P/A space:
    #
    # AA = phenotype absence
    # PA = phenotype presence
    # PP = phenotype presence
    #
    # This is the locked dominant-presence scoring model.

    called["prediction"] = np.where(
        called["biological_genotype"] == "AA",
        0,
        1,
    )

    called["phenotype_value"] = pd.to_numeric(
        called["phenotype_value"],
        errors="raise",
    ).astype(int)

    obs = called[
        "phenotype_value"
    ].to_numpy()

    pred = called[
        "prediction"
    ].to_numpy()

    tp = int(
        (
            (obs == 1)
            &
            (pred == 1)
        ).sum()
    )

    tn = int(
        (
            (obs == 0)
            &
            (pred == 0)
        ).sum()
    )

    fp = int(
        (
            (obs == 0)
            &
            (pred == 1)
        ).sum()
    )

    fn = int(
        (
            (obs == 1)
            &
            (pred == 0)
        ).sum()
    )

    p0 = int(
        (obs == 0).sum()
    )

    p1 = int(
        (obs == 1).sum()
    )

    sensitivity = (
        tp / p1
        if p1 > 0
        else np.nan
    )

    specificity = (
        tn / p0
        if p0 > 0
        else np.nan
    )

    match_n = tp + tn
    mismatch_n = fp + fn

    return {
        "cohort_n":
            len(d),

        "called_n":
            len(called),

        "missing_n":
            len(d) - len(called),

        "call_rate":
            len(called) / len(d)
            if len(d)
            else np.nan,

        "phenotype_0_called":
            p0,

        "phenotype_1_called":
            p1,

        "match_n":
            match_n,

        "mismatch_n":
            mismatch_n,

        "match_fraction":
            match_n / len(called),

        "sensitivity":
            sensitivity,

        "specificity":
            specificity,

        "balanced_accuracy":
            balanced_accuracy(
                obs,
                pred,
            ),
    }


def main():

    long = pd.read_csv(
        LONG,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    sites = pd.read_csv(
        SITES,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    callability = pd.read_csv(
        CALLABILITY,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    meta = pd.read_csv(
        META,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    # --------------------------------------------------------
    # Metadata: sequencing type.
    # --------------------------------------------------------

    seq = (
        meta[
            [
                "Accession",
                "Data Type",
            ]
        ]
        .rename(
            columns={
                "Accession": "sample",
                "Data Type": "data_type",
            }
        )
        .drop_duplicates(
            subset="sample"
        )
    )

    callability = callability.merge(
        seq,
        on="sample",
        how="left",
    )

    callability[
        "callable_fraction"
    ] = pd.to_numeric(
        callability[
            "callable_fraction"
        ],
        errors="coerce",
    )

    # Explicit taxonomic exclusion.
    callability[
        "contains_cydno"
    ] = (
        callability[
            "taxon_name"
        ]
        .astype(str)
        .str.contains(
            "cydno",
            case=False,
            regex=False,
            na=False,
        )
    )

    callability[
        "contains_numata"
    ] = (
        callability[
            "taxon_name"
        ]
        .astype(str)
        .str.contains(
            "numata",
            case=False,
            regex=False,
            na=False,
        )
    )

    # Additional locked taxonomic exclusion.
    callability[
        "contains_pachinus"
    ] = (
        callability[
            "taxon_name"
        ]
        .astype(str)
        .str.contains(
            "pachinus",
            case=False,
            regex=False,
            na=False,
        )
    )

    # Final locked taxonomic exclusion:
    # cydno + numata + pachinus.
    callability[
        "excluded_taxon"
    ] = (
        callability[
            "contains_cydno"
        ]
        |
        callability[
            "contains_numata"
        ]
        |
        callability[
            "contains_pachinus"
        ]
    )

    # Additional locked taxonomic exclusions.
    callability[
        "contains_hecale"
    ] = (
        callability[
            "taxon_name"
        ]
        .astype(str)
        .str.contains(
            "hecale",
            case=False,
            regex=False,
            na=False,
        )
    )

    callability[
        "contains_ismenius"
    ] = (
        callability[
            "taxon_name"
        ]
        .astype(str)
        .str.contains(
            "ismenius",
            case=False,
            regex=False,
            na=False,
        )
    )

    taxon_lower_extra = (
        callability[
            "taxon_name"
        ]
        .astype(str)
        .str.lower()
    )

    callability[
        "contains_timareta_contigua"
    ] = (
        taxon_lower_extra.str.contains(
            "timareta",
            regex=False,
            na=False,
        )
        &
        taxon_lower_extra.str.contains(
            "contigua",
            regex=False,
            na=False,
        )
    )

    callability[
        "contains_timareta_linaresi"
    ] = (
        taxon_lower_extra.str.contains(
            "timareta",
            regex=False,
            na=False,
        )
        &
        taxon_lower_extra.str.contains(
            "linaresi",
            regex=False,
            na=False,
        )
    )

    # Extend the previously locked cydno + numata + pachinus
    # exclusion with the additional requested taxa.
    callability[
        "excluded_taxon"
    ] = (
        callability[
            "excluded_taxon"
        ]
        |
        callability[
            "contains_hecale"
        ]
        |
        callability[
            "contains_ismenius"
        ]
        |
        callability[
            "contains_timareta_contigua"
        ]
        |
        callability[
            "contains_timareta_linaresi"
        ]
    )

    # Technical module-specific callability filter.
    callability[
        "passes_callability"
    ] = (
        callability[
            "callable_fraction"
        ]
        >= CALLABILITY_CUTOFF
    )

    callability[
        "base_keep"
    ] = (
        (~callability["excluded_taxon"])
        &
        callability["passes_callability"]
    )

    # --------------------------------------------------------
    # Define the two sensitivity-analysis cohorts.
    # --------------------------------------------------------

    cohort_rows = []

    for _, row in callability.iterrows():

        base = bool(
            row["base_keep"]
        )

        data_type = str(
            row["data_type"]
        ).strip()

        cohort_rows.append(
            {
                "cohort":
                    "WGS_only",

                "module":
                    row["module"],

                "sample":
                    row["sample"],

                "taxon_name":
                    row["taxon_name"],

                "data_type":
                    data_type,

                "callable_fraction":
                    row["callable_fraction"],

                "contains_cydno":
                    row["contains_cydno"],

                "contains_pachinus":
                    row["contains_pachinus"],

                "excluded_taxon":
                    row["excluded_taxon"],

                "passes_callability":
                    row["passes_callability"],

                "keep":
                    (
                        base
                        and
                        data_type == "WGS"
                    ),
            }
        )

        cohort_rows.append(
            {
                "cohort":
                    "WGS_plus_CS",

                "module":
                    row["module"],

                "sample":
                    row["sample"],

                "taxon_name":
                    row["taxon_name"],

                "data_type":
                    data_type,

                "callable_fraction":
                    row["callable_fraction"],

                "contains_cydno":
                    row["contains_cydno"],

                "contains_pachinus":
                    row["contains_pachinus"],

                "excluded_taxon":
                    row["excluded_taxon"],

                "passes_callability":
                    row["passes_callability"],

                "keep":
                    (
                        base
                        and
                        data_type
                        in {
                            "WGS",
                            "CS",
                        }
                    ),
            }
        )

    cohorts = pd.DataFrame(
        cohort_rows
    )

    cohorts.to_csv(
        OUT_COHORT,
        sep="\t",
        index=False,
    )

    print(
        "\n=============================================="
    )
    print(
        "FILTERED COHORT SIZES"
    )
    print(
        "=============================================="
    )

    retained = cohorts.loc[
        cohorts["keep"]
    ].copy()

    print(
        retained.groupby(
            [
                "cohort",
                "module",
            ]
        )
        .size()
        .to_string()
    )

    print(
        "\nSequencing composition of WGS+CS cohort:"
    )

    print(
        retained.loc[
            retained["cohort"]
            ==
            "WGS_plus_CS"
        ]
        .groupby(
            [
                "module",
                "data_type",
            ]
        )
        .size()
        .to_string()
    )

    # --------------------------------------------------------
    # Score every technically represented frozen SNP.
    # --------------------------------------------------------

    score_rows = []

    for cohort_name in [
        "WGS_only",
        "WGS_plus_CS",
    ]:

        for module in MODULES:

            keep_samples = set(
                cohorts.loc[
                    (
                        cohorts["cohort"]
                        == cohort_name
                    )
                    &
                    (
                        cohorts["module"]
                        == module
                    )
                    &
                    cohorts["keep"],
                    "sample",
                ]
            )

            module_sites = sites.loc[
                sites["module"]
                == module
            ].copy()

            for _, site in module_sites.iterrows():

                site_key = site[
                    "site_key"
                ]

                d = long.loc[
                    (
                        long["module"]
                        == module
                    )
                    &
                    (
                        long["site_key"]
                        == site_key
                    )
                    &
                    (
                        long["sample"]
                        .isin(
                            keep_samples
                        )
                    )
                ].copy()

                metrics = score_site(
                    d
                )

                score_rows.append(
                    {
                        "cohort":
                            cohort_name,

                        "module":
                            module,

                        "site_key":
                            site_key,

                        "chrom":
                            site["chrom"],

                        "pos":
                            site["pos"],

                        "ref":
                            site["ref"],

                        "alt":
                            site["alt"],

                        "source_status":
                            site[
                                "source_status"
                            ],

                        "curated51_BA":
                            site[
                                "plot_BA_curated"
                            ],

                        **metrics,
                    }
                )

    scores = pd.DataFrame(
        score_rows
    )

    scores.to_csv(
        OUT_SCORES,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    # --------------------------------------------------------
    # Threshold summaries.
    # --------------------------------------------------------

    summary_rows = []

    for (
        cohort_name,
        module,
    ), d in scores.groupby(
        [
            "cohort",
            "module",
        ],
        sort=False,
    ):

        represented = d.loc[
            d["source_status"]
            ==
            "exact_REF_ALT_match"
        ].copy()

        ba = pd.to_numeric(
            represented[
                "balanced_accuracy"
            ],
            errors="coerce",
        )

        match_fraction = pd.to_numeric(
            represented[
                "match_fraction"
            ],
            errors="coerce",
        )

        row = {
            "cohort":
                cohort_name,

            "module":
                module,

            "frozen_sites_total":
                len(d),

            "technically_represented_sites":
                len(represented),

            "sites_with_score":
                int(
                    ba.notna().sum()
                ),

            "BA_min":
                ba.min(),

            "BA_median":
                ba.median(),

            "BA_max":
                ba.max(),

            "mean_match_fraction":
                match_fraction.mean(),

            "perfect_match_sites":
                int(
                    (
                        match_fraction
                        == 1.0
                    ).sum()
                ),
        }

        for threshold in [
            0.80,
            0.85,
            0.90,
            0.95,
        ]:

            row[
                f"sites_BA_ge_{threshold:.2f}"
            ] = int(
                (
                    ba
                    >= threshold
                ).sum()
            )

        summary_rows.append(
            row
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
        "\n=============================================="
    )
    print(
        "FROZEN-SITE SCORE SUMMARY"
    )
    print(
        "=============================================="
    )

    cols = [
        "cohort",
        "module",
        "technically_represented_sites",
        "sites_with_score",
        "BA_median",
        "BA_max",
        "sites_BA_ge_0.80",
        "sites_BA_ge_0.85",
        "sites_BA_ge_0.90",
        "sites_BA_ge_0.95",
        "perfect_match_sites",
    ]

    print(
        summary[
            cols
        ].to_string(
            index=False
        )
    )

    print(
        "\n=============================================="
    )
    print(
        "INDIVIDUAL-SNP MATCH COUNTS"
    )
    print(
        "=============================================="
    )

    for cohort_name in [
        "WGS_only",
        "WGS_plus_CS",
    ]:

        for module in MODULES:

            d = scores.loc[
                (
                    scores["cohort"]
                    == cohort_name
                )
                &
                (
                    scores["module"]
                    == module
                )
                &
                (
                    scores["source_status"]
                    ==
                    "exact_REF_ALT_match"
                )
            ]

            print(
                f"\n{cohort_name} / {module}"
            )

            print(
                "  total matching genotype×individual calls:",
                int(
                    pd.to_numeric(
                        d["match_n"]
                    ).sum()
                ),
            )

            print(
                "  total mismatching genotype×individual calls:",
                int(
                    pd.to_numeric(
                        d["mismatch_n"]
                    ).sum()
                ),
            )

    print(
        "\nOutputs:"
    )

    print(
        f"  {OUT_COHORT}"
    )

    print(
        f"  {OUT_SCORES}"
    )

    print(
        f"  {OUT_SUMMARY}"
    )


if __name__ == "__main__":
    main()
