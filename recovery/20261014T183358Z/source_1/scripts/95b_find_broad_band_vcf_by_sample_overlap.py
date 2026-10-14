#!/usr/bin/env python3

from pathlib import Path
import gzip
import pandas as pd


PROJECT = Path(
    "/scratch/martinlab/joe/GynandroCalling"
)

MANIFEST = (
    PROJECT
    / "analysis"
    / "original_502"
    / "checkpoints"
    / "original502_filtered_cohort_manifest.WGS_vs_WGSplusCS.tsv"
)


def as_bool(x):
    return str(x).strip().lower() in {
        "true",
        "t",
        "1",
        "yes",
        "y",
    }


def read_vcf_samples(path):

    opener = (
        gzip.open
        if str(path).endswith(".gz")
        else open
    )

    try:
        with opener(
            path,
            "rt",
            errors="replace",
        ) as handle:

            for line in handle:

                if line.startswith("#CHROM"):
                    return (
                        line.rstrip("\n")
                        .split("\t")[9:]
                    )

                if not line.startswith("#"):
                    break

    except Exception as exc:
        return None

    return None


def main():

    m = pd.read_csv(
        MANIFEST,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    keep = m["keep"].map(
        as_bool
    )

    band = m.loc[
        keep
        &
        (m["module"] == "Band")
    ].copy()

    targets = {}

    for cohort in [
        "WGS_only",
        "WGS_plus_CS",
    ]:

        targets[cohort] = set(
            band.loc[
                band["cohort"] == cohort,
                "sample",
            ]
        )

    print(
        "Final Band cohorts:"
    )

    for cohort, samples in targets.items():
        print(
            f"  {cohort}: {len(samples)}"
        )

    print(
        "\nScanning VCF headers..."
    )

    vcfs = sorted(
        list(
            PROJECT.rglob(
                "*.vcf"
            )
        )
        +
        list(
            PROJECT.rglob(
                "*.vcf.gz"
            )
        )
    )

    rows = []

    for path in vcfs:

        samples = read_vcf_samples(
            path
        )

        if samples is None:
            continue

        sample_set = set(
            samples
        )

        rows.append(
            {
                "path":
                    str(path),

                "vcf_sample_n":
                    len(samples),

                "WGS_only_overlap":
                    len(
                        sample_set
                        &
                        targets["WGS_only"]
                    ),

                "WGS_plus_CS_overlap":
                    len(
                        sample_set
                        &
                        targets["WGS_plus_CS"]
                    ),

                "WGS_only_fraction":
                    (
                        len(
                            sample_set
                            &
                            targets["WGS_only"]
                        )
                        /
                        len(
                            targets["WGS_only"]
                        )
                    ),

                "WGS_plus_CS_fraction":
                    (
                        len(
                            sample_set
                            &
                            targets["WGS_plus_CS"]
                        )
                        /
                        len(
                            targets["WGS_plus_CS"]
                        )
                    ),
            }
        )

    d = pd.DataFrame(
        rows
    )

    if d.empty:
        print(
            "No readable VCF headers found."
        )
        return

    d = d.sort_values(
        [
            "WGS_plus_CS_overlap",
            "WGS_only_overlap",
            "vcf_sample_n",
        ],
        ascending=[
            False,
            False,
            False,
        ],
    )

    print(
        "\n=============================================="
    )

    print(
        "TOP VCFs BY FINAL BAND COHORT OVERLAP"
    )

    print(
        "=============================================="
    )

    print(
        d.head(
            40
        ).to_string(
            index=False
        )
    )

    out = (
        PROJECT
        / "analysis"
        / "original_502"
        / "checkpoints"
        / "band_vcf_sample_overlap_audit.tsv"
    )

    d.to_csv(
        out,
        sep="\t",
        index=False,
    )

    print(
        "\nOutput:"
    )

    print(
        f"  {out}"
    )


if __name__ == "__main__":
    main()
