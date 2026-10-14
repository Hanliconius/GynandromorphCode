#!/usr/bin/env python3

from pathlib import Path
import re

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

SITES = (
    CK
    / "frozen_comparative_sites_reference_polarized.tsv"
)

COHORT = (
    CK
    / "optix_original_502_cohort_manifest.tsv"
)

VCFS = {
    "Band": (
        PROJECT
        / "external"
        / "krista_genotype_plots_original"
        / "Genotype_Plots"
        / "optix_band-1112651-1114314.vcf"
    ),
    "Ray2": (
        PROJECT
        / "external"
        / "krista_genotype_plots_original"
        / "Genotype_Plots"
        / "ray2-1174547-1187833.vcf"
    ),
    "Dennis2": (
        PROJECT
        / "external"
        / "krista_genotype_plots_original"
        / "Genotype_Plots"
        / "dennis2-1218557-1221494.vcf"
    ),
}

OUT_LONG = (
    CK
    / "optix_original502_selected_sites_genotype_plot.long.tsv"
)

OUT_SITES = (
    CK
    / "optix_original502_selected_sites_for_plotting.tsv"
)

MODULE_ORDER = [
    "Band",
    "Ray2",
    "Dennis2",
]

PANEL_NAME = {
    "Band": "band",
    "Ray2": "ray2",
    "Dennis2": "dennis2",
}

REFERENCE_STATE = {
    "Band": "presence",
    "Ray2": "absence",
    "Dennis2": "absence",
}


def fail(message):
    raise RuntimeError(message)


def read_tsv(path):
    if not path.exists():
        fail(
            f"Missing required file: {path}"
        )

    return pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )


def emptyish(value):
    return (
        value is None
        or str(value).strip()
        in {
            "",
            ".",
            "NA",
            "NaN",
            "nan",
            "None",
        }
    )


def norm(value):
    if emptyish(value):
        return ""

    return str(value).strip()


def norm_allele(value):
    return norm(
        value
    ).upper()


def fnum(value):
    if emptyish(value):
        return np.nan

    try:
        return float(value)
    except Exception:
        return np.nan


def parse_gt(gt):
    if emptyish(gt):
        return None

    gt = str(gt).strip()

    if gt in {
        ".",
        "./.",
        ".|.",
    }:
        return None

    parts = re.split(
        r"[/|]",
        gt,
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


def biological_gt(
    module,
    raw_gt,
):
    alt_copies = parse_gt(
        raw_gt
    )

    if alt_copies is None:
        return None

    if (
        REFERENCE_STATE[module]
        == "presence"
    ):
        presence_copies = (
            2 - alt_copies
        )
    else:
        presence_copies = (
            alt_copies
        )

    return {
        0: "AA",
        1: "PA",
        2: "PP",
    }.get(
        presence_copies
    )


def genotype_code(gt):
    return {
        "AA": 0,
        "PA": 1,
        "PP": 2,
    }.get(
        gt,
        np.nan,
    )


def predicted_presence(gt):
    """
    Dominant phenotype-presence model in biological
    P/A genotype space:

        AA -> phenotype absent
        PA -> phenotype present
        PP -> phenotype present
    """

    if gt == "AA":
        return 0

    if gt in {
        "PA",
        "PP",
    }:
        return 1

    return np.nan


def balanced_accuracy(
    observed,
    predicted,
):
    observed = np.asarray(
        observed,
        dtype=float,
    )

    predicted = np.asarray(
        predicted,
        dtype=float,
    )

    keep = (
        np.isfinite(
            observed
        )
        &
        np.isfinite(
            predicted
        )
    )

    observed = observed[
        keep
    ]

    predicted = predicted[
        keep
    ]

    if len(
        observed
    ) == 0:
        return np.nan

    present = (
        observed
        == 1
    )

    absent = (
        observed
        == 0
    )

    if (
        present.sum() == 0
        or absent.sum() == 0
    ):
        return np.nan

    sensitivity = (
        predicted[
            present
        ]
        == 1
    ).mean()

    specificity = (
        predicted[
            absent
        ]
        == 0
    ).mean()

    return (
        sensitivity
        + specificity
    ) / 2.0


def parse_vcf(
    path,
    wanted_positions,
):
    if not path.exists():
        fail(
            f"Missing regional VCF: {path}"
        )

    samples = None

    records = {}

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
                fields = (
                    line.rstrip("\n")
                    .split("\t")
                )

                samples = fields[9:]

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

            format_names = (
                fields[8]
                .split(":")
            )

            calls = {}

            for sample, text in zip(
                samples,
                fields[9:],
            ):

                values = (
                    text.split(":")
                )

                calls[
                    sample
                ] = {
                    key: (
                        values[index]
                        if index
                        < len(values)
                        else ""
                    )
                    for index, key
                    in enumerate(
                        format_names
                    )
                }

            records[
                (
                    chrom,
                    pos,
                    ref,
                    alt,
                )
            ] = calls

    if samples is None:
        fail(
            f"Could not parse VCF header: {path}"
        )

    return (
        samples,
        records,
    )


def main():

    sites = read_tsv(
        SITES
    )

    cohort = read_tsv(
        COHORT
    )

    sites = sites.loc[
        sites[
            "module"
        ].isin(
            MODULE_ORDER
        )
    ].copy()

    sites[
        "pos"
    ] = pd.to_numeric(
        sites[
            "pos"
        ],
        errors="raise",
    ).astype(int)

    sites[
        "ref"
    ] = sites[
        "ref"
    ].map(
        norm_allele
    )

    sites[
        "alt"
    ] = sites[
        "alt"
    ].map(
        norm_allele
    )

    # Retain curated-panel BA for comparison.
    if (
        "minimum_BA"
        in sites.columns
    ):
        sites[
            "plot_BA_curated"
        ] = pd.to_numeric(
            sites[
                "minimum_BA"
            ],
            errors="coerce",
        )

    elif (
        "minimum_comparative_BA"
        in sites.columns
    ):
        sites[
            "plot_BA_curated"
        ] = pd.to_numeric(
            sites[
                "minimum_comparative_BA"
            ],
            errors="coerce",
        )

    else:
        sites[
            "plot_BA_curated"
        ] = np.nan

    module_rank = {
        module: index
        for index, module
        in enumerate(
            MODULE_ORDER
        )
    }

    sites[
        "_module_rank"
    ] = sites[
        "module"
    ].map(
        module_rank
    )

    sites = (
        sites.sort_values(
            [
                "_module_rank",
                "pos",
            ]
        )
        .drop(
            columns="_module_rank"
        )
        .reset_index(
            drop=True
        )
    )

    sites[
        "site_order_within_module"
    ] = (
        sites.groupby(
            "module"
        )
        .cumcount()
        + 1
    )

    cohort[
        "sample"
    ] = cohort[
        "sample"
    ].map(
        norm
    )

    cohort[
        "panel_lower"
    ] = (
        cohort[
            "panel"
        ]
        .astype(str)
        .str.strip()
        .str.lower()
    )

    cohort[
        "phenotype_num"
    ] = pd.to_numeric(
        cohort[
            "phenotype"
        ],
        errors="coerce",
    )

    cohort[
        "reason_lower"
    ] = (
        cohort[
            "reason"
        ]
        .astype(str)
        .str.lower()
    )

    cohort[
        "taxon_clean"
    ] = (
        cohort[
            "taxon"
        ]
        .astype(str)
        .str.replace(
            r"\s+",
            " ",
            regex=True,
        )
        .str.strip()
    )

    long_rows = []

    site_rows = []

    for module in MODULE_ORDER:

        panel = PANEL_NAME[
            module
        ]

        module_sites = sites.loc[
            sites[
                "module"
            ]
            == module
        ].copy()

        wanted_positions = set(
            module_sites[
                "pos"
            ].astype(int)
        )

        (
            vcf_samples,
            vcf_records,
        ) = parse_vcf(
            VCFS[
                module
            ],
            wanted_positions,
        )

        vcf_sample_set = set(
            vcf_samples
        )

        module_cohort = cohort.loc[
            cohort[
                "panel_lower"
            ]
            == panel
        ].copy()

        module_cohort = module_cohort.loc[
            module_cohort[
                "phenotype_num"
            ].isin(
                [
                    0,
                    1,
                ]
            )
        ].copy()

        # Prefer the exact comparative-scoring cohort encoded
        # in the manifest.
        explicit_comparative = (
            module_cohort[
                "reason_lower"
            ]
            .str.contains(
                "comparative",
                regex=False,
                na=False,
            )
        )

        if (
            explicit_comparative.any()
        ):

            module_cohort = (
                module_cohort.loc[
                    explicit_comparative
                ]
                .copy()
            )

        else:

            module_cohort = (
                module_cohort.loc[
                    ~module_cohort[
                        "reason_lower"
                    ].str.contains(
                        "focal|piano|exclude|missing",
                        regex=True,
                        na=False,
                    )
                ]
                .copy()
            )

        module_cohort = (
            module_cohort.loc[
                ~module_cohort[
                    "sample"
                ].str.startswith(
                    "TN_",
                    na=False,
                )
            ]
            .copy()
        )

        module_cohort = (
            module_cohort.loc[
                module_cohort[
                    "sample"
                ].isin(
                    vcf_sample_set
                )
            ]
            .copy()
        )

        module_cohort[
            "phenotype_value"
        ] = module_cohort[
            "phenotype_num"
        ].astype(int)

        # ----------------------------------------------------
        # Known Band metadata corrections.
        #
        # malleti and ecuadoriensis are Band-absent.
        # ----------------------------------------------------

        if module == "Band":

            taxon_lower = (
                module_cohort[
                    "taxon_clean"
                ]
                .astype(str)
                .str.lower()
            )

            force_absent = (
                taxon_lower.str.contains(
                    "malleti",
                    regex=False,
                )
                |
                taxon_lower.str.contains(
                    "ecuadoriensis",
                    regex=False,
                )
            )

            originally_present = (
                force_absent
                &
                (
                    module_cohort[
                        "phenotype_value"
                    ]
                    == 1
                )
            )

            module_cohort.loc[
                force_absent,
                "phenotype_value",
            ] = 0

            print(
                "\nBand forced-absent samples "
                "(malleti + ecuadoriensis):",
                int(
                    force_absent.sum()
                ),
            )

            print(
                "Band values actually changed "
                "from present to absent:",
                int(
                    originally_present.sum()
                ),
            )

        module_cohort[
            "phenotype_label"
        ] = np.where(
            module_cohort[
                "phenotype_value"
            ]
            == 1,
            "present",
            "absent",
        )

        module_cohort = (
            module_cohort.sort_values(
                [
                    "phenotype_value",
                    "taxon_clean",
                    "sample",
                ],
                ascending=[
                    True,
                    True,
                    True,
                ],
            )
            .reset_index(
                drop=True
            )
        )

        module_cohort[
            "sample_order"
        ] = np.arange(
            1,
            len(
                module_cohort
            ) + 1,
        )

        print(
            f"\n{module}"
        )

        print(
            "comparative samples:",
            len(
                module_cohort
            ),
        )

        print(
            module_cohort[
                "phenotype_label"
            ]
            .value_counts()
            .to_string()
        )

        for _, site in (
            module_sites.iterrows()
        ):

            chrom = site[
                "chrom"
            ]

            pos = int(
                site[
                    "pos"
                ]
            )

            ref = site[
                "ref"
            ]

            alt = site[
                "alt"
            ]

            site_key = (
                f"{chrom}:"
                f"{pos}:"
                f"{ref}>"
                f"{alt}"
            )

            exact_key = (
                chrom,
                pos,
                ref,
                alt,
            )

            calls = (
                vcf_records.get(
                    exact_key
                )
            )

            if calls is None:

                position_matches = [
                    key
                    for key
                    in vcf_records
                    if (
                        key[0]
                        == chrom
                        and key[1]
                        == pos
                    )
                ]

                if position_matches:
                    source_status = (
                        "allele_mismatch_at_position"
                    )
                else:
                    source_status = (
                        "site_absent_from_regional_VCF"
                    )

            else:
                source_status = (
                    "exact_REF_ALT_match"
                )

            site_observed = []

            site_predicted = []

            site_callable = 0

            for _, individual in (
                module_cohort.iterrows()
            ):

                sample = individual[
                    "sample"
                ]

                phenotype = int(
                    individual[
                        "phenotype_value"
                    ]
                )

                if calls is None:
                    call = None
                else:
                    call = calls.get(
                        sample
                    )

                if call is None:

                    raw_gt = ""
                    raw_ad = ""
                    dp = np.nan
                    gq = np.nan

                else:

                    raw_gt = call.get(
                        "GT",
                        "",
                    )

                    raw_ad = call.get(
                        "AD",
                        "",
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

                bgt = biological_gt(
                    module,
                    raw_gt,
                )

                qc = bool(
                    calls is not None
                    and bgt is not None
                    and np.isfinite(
                        dp
                    )
                    and dp >= 2
                    and np.isfinite(
                        gq
                    )
                    and gq >= 20
                )

                if qc:

                    prediction = (
                        predicted_presence(
                            bgt
                        )
                    )

                    site_observed.append(
                        phenotype
                    )

                    site_predicted.append(
                        prediction
                    )

                    site_callable += 1

                else:

                    prediction = np.nan

                ref_depth, alt_depth = (
                    parse_ad(
                        raw_ad
                    )
                )

                long_rows.append(
                    {
                        "module":
                            module,

                        "site_key":
                            site_key,

                        "site_order_within_module":
                            int(
                                site[
                                    "site_order_within_module"
                                ]
                            ),

                        "chrom":
                            chrom,

                        "pos":
                            pos,

                        "ref":
                            ref,

                        "alt":
                            alt,

                        "source_status":
                            source_status,

                        "sample":
                            sample,

                        "sample_order":
                            int(
                                individual[
                                    "sample_order"
                                ]
                            ),

                        # Display field: one-line taxon only.
                        # Accession remains internal and is
                        # never needed as the plot label.
                        "taxon_name":
                            individual[
                                "taxon_clean"
                            ],

                        "phenotype_value":
                            phenotype,

                        "phenotype_label":
                            individual[
                                "phenotype_label"
                            ],

                        "GT_raw":
                            raw_gt,

                        "AD_raw":
                            raw_ad,

                        "DP":
                            dp,

                        "GQ":
                            gq,

                        "QC_pass_DP2_GQ20":
                            qc,

                        "biological_genotype":
                            (
                                bgt
                                if qc
                                else "missing"
                            ),

                        "genotype_code":
                            (
                                genotype_code(
                                    bgt
                                )
                                if qc
                                else np.nan
                            ),

                        "ref_depth":
                            (
                                ref_depth
                                if qc
                                else np.nan
                            ),

                        "alt_depth":
                            (
                                alt_depth
                                if qc
                                else np.nan
                            ),

                        "predicted_phenotype":
                            prediction,
                    }
                )

            full_ba = balanced_accuracy(
                site_observed,
                site_predicted,
            )

            site_rows.append(
                {
                    "module":
                        module,

                    "site_key":
                        site_key,

                    "site_order_within_module":
                        int(
                            site[
                                "site_order_within_module"
                            ]
                        ),

                    "chrom":
                        chrom,

                    "pos":
                        pos,

                    "ref":
                        ref,

                    "alt":
                        alt,

                    "source_status":
                        source_status,

                    "plot_BA_curated":
                        fnum(
                            site[
                                "plot_BA_curated"
                            ]
                        ),

                    "plot_BA":
                        full_ba,

                    "plot_BA_fullpanel":
                        full_ba,

                    "fullpanel_callable_n":
                        site_callable,

                    "fullpanel_comparative_n":
                        len(
                            module_cohort
                        ),
                }
            )

    long = pd.DataFrame(
        long_rows
    )

    plot_sites = pd.DataFrame(
        site_rows
    )

    if long.empty:
        fail(
            "No full-panel genotype rows built."
        )

    long.to_csv(
        OUT_LONG,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    plot_sites.to_csv(
        OUT_SITES,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    print(
        "\nFULL-PANEL SITE COUNTS"
    )

    print(
        plot_sites.groupby(
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
        "\nSOURCE REPRESENTATION"
    )

    print(
        plot_sites.groupby(
            [
                "module",
                "source_status",
            ]
        )
        .size()
        .to_string()
    )

    print(
        "\nFULL-PANEL BA SUMMARY"
    )

    print(
        plot_sites.groupby(
            "module"
        )[
            "plot_BA"
        ]
        .agg(
            [
                "count",
                "min",
                "median",
                "max",
            ]
        )
        .to_string()
    )

    print(
        "\nOutputs:"
    )

    print(
        f"  {OUT_LONG}"
    )

    print(
        f"  {OUT_SITES}"
    )


if __name__ == "__main__":
    main()
