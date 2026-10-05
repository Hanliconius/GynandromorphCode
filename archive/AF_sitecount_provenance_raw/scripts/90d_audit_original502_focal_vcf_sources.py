#!/usr/bin/env python3

from pathlib import Path
import subprocess
import sys

import pandas as pd


PROJECT = Path(
    "/scratch/martinlab/joe/GynandroCalling"
)

PROVENANCE = (
    PROJECT
    / "analysis"
    / "checkpoints"
    / "complete_vcf_provenance_manifest.tsv"
)

SITES = (
    PROJECT
    / "analysis"
    / "original_502"
    / "checkpoints"
    / "frozen_comparative_sites_reference_polarized.tsv"
)

OUTDIR = (
    PROJECT
    / "analysis"
    / "original_502"
    / "checkpoints"
)

OUT_MATCHES = (
    OUTDIR
    / "frozen_site_vcf_source_matches.tsv"
)

OUT_SUMMARY = (
    OUTDIR
    / "focal_vcf_source_coverage_by_module.tsv"
)

OUT_BEST = (
    OUTDIR
    / "best_focal_vcf_source_by_module.tsv"
)

MODULE_ORDER = [
    "Band",
    "Ray1",
    "Ray2",
    "Dennis1",
    "Dennis2",
    "Yb",
]


def fail(message):
    raise RuntimeError(message)


def clean(value):
    if pd.isna(value):
        return ""

    return str(value).strip()


def resolve_path(relative_path):
    path = Path(
        clean(relative_path)
    )

    if path.is_absolute():
        return path

    return PROJECT / path


def get_samples(vcf):
    result = subprocess.run(
        [
            "bcftools",
            "query",
            "-l",
            str(vcf),
        ],
        text=True,
        capture_output=True,
    )

    if result.returncode != 0:
        return [], result.stderr.strip()

    samples = [
        x.strip()
        for x in result.stdout.splitlines()
        if x.strip()
    ]

    return samples, ""


def query_records(vcf):
    result = subprocess.run(
        [
            "bcftools",
            "query",
            "-f",
            "%CHROM\t%POS\t%REF\t%ALT\n",
            str(vcf),
        ],
        text=True,
        capture_output=True,
    )

    if result.returncode != 0:
        return None, result.stderr.strip()

    records = {}

    for line in result.stdout.splitlines():
        fields = line.split("\t")

        if len(fields) < 4:
            continue

        chrom = fields[0]
        pos = int(fields[1])
        ref = fields[2].upper()
        alt = fields[3].upper()

        records.setdefault(
            (chrom, pos),
            []
        ).append(
            {
                "ref": ref,
                "alt": alt,
            }
        )

    return records, ""


def classify_match(
    selected_ref,
    selected_alt,
    records,
):
    if not records:
        return (
            "absent",
            "",
        )

    selected_ref = (
        str(selected_ref)
        .upper()
    )

    selected_alt = (
        str(selected_alt)
        .upper()
    )

    observed = ";".join(
        (
            f"{record['ref']}>"
            f"{record['alt']}"
        )
        for record in records
    )

    # Exact representation first.
    for record in records:
        if (
            record["ref"]
            == selected_ref
            and record["alt"]
            == selected_alt
        ):
            return (
                "exact_biallelic",
                observed,
            )

    # Same REF, with selected ALT present in a
    # multi-allelic representation.
    for record in records:
        alts = (
            record["alt"]
            .split(",")
        )

        if (
            record["ref"]
            == selected_ref
            and selected_alt
            in alts
        ):
            return (
                "selected_alt_in_multiallelic",
                observed,
            )

    # Swapped simple REF/ALT representation.
    for record in records:
        if (
            record["ref"]
            == selected_alt
            and record["alt"]
            == selected_ref
        ):
            return (
                "ref_alt_swapped",
                observed,
            )

    return (
        "position_only",
        observed,
    )


def main():
    if not PROVENANCE.exists():
        fail(
            f"Missing provenance manifest: "
            f"{PROVENANCE}"
        )

    if not SITES.exists():
        fail(
            f"Missing frozen site table: "
            f"{SITES}"
        )

    provenance = pd.read_csv(
        PROVENANCE,
        sep="\t",
        dtype=str,
    )

    sites = pd.read_csv(
        SITES,
        sep="\t",
        dtype=str,
    )

    required_provenance = {
        "relative_path",
    }

    missing = (
        required_provenance
        - set(
            provenance.columns
        )
    )

    if missing:
        fail(
            "Provenance manifest lacks: "
            + ", ".join(
                sorted(missing)
            )
        )

    sites["pos"] = pd.to_numeric(
        sites["pos"],
        errors="raise",
    ).astype(int)

    # --------------------------------------------------------
    # Work only from provenance-tracked VCFs.
    # Deduplicate paths because the manifest may contain
    # equivalent records.
    # --------------------------------------------------------

    source_rows = []

    seen_paths = set()

    for _, row in provenance.iterrows():
        relative = clean(
            row["relative_path"]
        )

        if not relative:
            continue

        path = resolve_path(
            relative
        )

        path_string = str(path)

        if path_string in seen_paths:
            continue

        seen_paths.add(
            path_string
        )

        if not path.exists():
            source_rows.append(
                {
                    "source_path":
                        path_string,
                    "exists":
                        False,
                    "has_TN_Left":
                        False,
                    "has_TN_Right":
                        False,
                    "has_both_focal":
                        False,
                    "sample_error":
                        "file_missing",
                }
            )

            continue

        samples, sample_error = (
            get_samples(
                path
            )
        )

        has_left = (
            "TN_Left"
            in samples
        )

        has_right = (
            "TN_Right"
            in samples
        )

        source_rows.append(
            {
                "source_path":
                    path_string,
                "exists":
                    True,
                "has_TN_Left":
                    has_left,
                "has_TN_Right":
                    has_right,
                "has_both_focal":
                    (
                        has_left
                        and has_right
                    ),
                "sample_error":
                    sample_error,
            }
        )

    source_table = pd.DataFrame(
        source_rows
    )

    focal_sources = (
        source_table.loc[
            source_table[
                "has_both_focal"
            ]
            == True
        ]
        .copy()
    )

    print(
        "\nPROVENANCE-TRACKED VCF SOURCES"
    )

    print(
        source_table.to_string(
            index=False
        )
    )

    print(
        "\nVCFS CONTAINING BOTH TN_Left "
        "AND TN_Right:",
        len(focal_sources),
    )

    if focal_sources.empty:
        fail(
            "No provenance-tracked VCF "
            "contains both focal samples."
        )

    # --------------------------------------------------------
    # Audit every frozen site against every focal VCF.
    # --------------------------------------------------------

    match_rows = []

    for _, source in (
        focal_sources.iterrows()
    ):
        vcf = Path(
            source["source_path"]
        )

        print(
            f"\nReading site representation: "
            f"{vcf}"
        )

        records, error = (
            query_records(
                vcf
            )
        )

        if records is None:
            print(
                f"WARNING: could not read "
                f"{vcf}: {error}",
                file=sys.stderr,
            )

            continue

        for _, site in sites.iterrows():
            key = (
                str(
                    site["chrom"]
                ),
                int(
                    site["pos"]
                ),
            )

            match_class, observed = (
                classify_match(
                    site["ref"],
                    site["alt"],
                    records.get(
                        key,
                        [],
                    ),
                )
            )

            match_rows.append(
                {
                    "module":
                        site["module"],
                    "site_id":
                        site["site_id"],
                    "chrom":
                        site["chrom"],
                    "pos":
                        int(
                            site["pos"]
                        ),
                    "selected_ref":
                        site["ref"],
                    "selected_alt":
                        site["alt"],
                    "source_vcf":
                        str(vcf),
                    "match_class":
                        match_class,
                    "observed_alleles":
                        observed,
                }
            )

    matches = pd.DataFrame(
        match_rows
    )

    matches.to_csv(
        OUT_MATCHES,
        sep="\t",
        index=False,
    )

    if matches.empty:
        fail(
            "No VCF/site comparisons "
            "were completed."
        )

    # --------------------------------------------------------
    # Source coverage summary.
    #
    # usable representation:
    # exact biallelic or selected ALT contained in a
    # multiallelic record.
    #
    # Source ranking uses representation ONLY.
    # No focal genotype or phenotypic expectation is used.
    # --------------------------------------------------------

    matches["exact"] = (
        matches[
            "match_class"
        ]
        == "exact_biallelic"
    ).astype(int)

    matches["multiallelic_usable"] = (
        matches[
            "match_class"
        ]
        == "selected_alt_in_multiallelic"
    ).astype(int)

    matches["swapped"] = (
        matches[
            "match_class"
        ]
        == "ref_alt_swapped"
    ).astype(int)

    matches["position_only"] = (
        matches[
            "match_class"
        ]
        == "position_only"
    ).astype(int)

    matches["position_present"] = (
        matches[
            "match_class"
        ]
        != "absent"
    ).astype(int)

    matches["usable_selected_allele"] = (
        matches["exact"]
        + matches[
            "multiallelic_usable"
        ]
    )

    summary = (
        matches.groupby(
            [
                "module",
                "source_vcf",
            ],
            dropna=False,
        )
        .agg(
            frozen_sites=(
                "site_id",
                "nunique",
            ),
            exact_biallelic=(
                "exact",
                "sum",
            ),
            selected_alt_in_multiallelic=(
                "multiallelic_usable",
                "sum",
            ),
            usable_selected_allele=(
                "usable_selected_allele",
                "sum",
            ),
            ref_alt_swapped=(
                "swapped",
                "sum",
            ),
            position_only=(
                "position_only",
                "sum",
            ),
            position_present=(
                "position_present",
                "sum",
            ),
        )
        .reset_index()
    )

    summary[
        "usable_fraction"
    ] = (
        summary[
            "usable_selected_allele"
        ]
        / summary[
            "frozen_sites"
        ]
    )

    summary[
        "position_fraction"
    ] = (
        summary[
            "position_present"
        ]
        / summary[
            "frozen_sites"
        ]
    )

    module_rank = {
        module: i
        for i, module in enumerate(
            MODULE_ORDER
        )
    }

    summary[
        "_module_rank"
    ] = summary[
        "module"
    ].map(
        module_rank
    )

    summary = (
        summary.sort_values(
            [
                "_module_rank",
                "usable_selected_allele",
                "position_present",
                "exact_biallelic",
                "source_vcf",
            ],
            ascending=[
                True,
                False,
                False,
                False,
                True,
            ],
        )
        .drop(
            columns=[
                "_module_rank",
            ]
        )
    )

    summary.to_csv(
        OUT_SUMMARY,
        sep="\t",
        index=False,
    )

    # --------------------------------------------------------
    # Best technical source for each module.
    #
    # Again: selected only by frozen-site representation.
    # --------------------------------------------------------

    best_rows = []

    for module in MODULE_ORDER:
        d = summary.loc[
            summary[
                "module"
            ]
            == module
        ].copy()

        if d.empty:
            continue

        d = d.sort_values(
            [
                "usable_selected_allele",
                "position_present",
                "exact_biallelic",
                "source_vcf",
            ],
            ascending=[
                False,
                False,
                False,
                True,
            ],
        )

        best_rows.append(
            d.iloc[0]
        )

    best = pd.DataFrame(
        best_rows
    )

    best.to_csv(
        OUT_BEST,
        sep="\t",
        index=False,
    )

    print(
        "\nBEST REPRESENTATION SOURCE "
        "BY MODULE"
    )

    print(
        best[
            [
                "module",
                "frozen_sites",
                "usable_selected_allele",
                "position_present",
                "exact_biallelic",
                "selected_alt_in_multiallelic",
                "ref_alt_swapped",
                "position_only",
                "usable_fraction",
                "source_vcf",
            ]
        ].to_string(
            index=False
        )
    )

    print(
        "\nTOP THREE SOURCES PER MODULE"
    )

    for module in MODULE_ORDER:
        d = summary.loc[
            summary[
                "module"
            ]
            == module
        ].head(3)

        if d.empty:
            continue

        print(
            f"\n{module}"
        )

        print(
            d[
                [
                    "usable_selected_allele",
                    "position_present",
                    "exact_biallelic",
                    "frozen_sites",
                    "source_vcf",
                ]
            ].to_string(
                index=False
            )
        )

    print(
        "\nOutputs:"
    )

    print(
        f"  {OUT_MATCHES}"
    )

    print(
        f"  {OUT_SUMMARY}"
    )

    print(
        f"  {OUT_BEST}"
    )


if __name__ == "__main__":
    main()
