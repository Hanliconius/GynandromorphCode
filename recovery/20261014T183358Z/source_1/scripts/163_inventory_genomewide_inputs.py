#!/usr/bin/env python3

import io
import json
import re
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
)

OUT = (
    PROJECT
    / "analysis"
    / "genomewide_mixing_input"
    / "metadata"
)

METRICS = (
    PROJECT
    / "analysis"
    / "autosomal_three_vs_four_founder"
    / "eligible_control_callability_and_heterozygosity.tsv"
)

BIOPROJECT_UID = "1473931"
REFERENCE = "GCA_900068215.1"

TARGET_CONTROLS = 18
MIN_CALL_RATE = 0.65

EUTILS = (
    "https://eutils.ncbi.nlm.nih.gov/"
    "entrez/eutils"
)

FORBIDDEN = (
    "TN2L_S1",
    "TN2R_S2",
    "TN3L_S3",
    "TN3R_S4",
)


def get(endpoint, **params):
    params["tool"] = "gynandrocalling"

    url = (
        f"{EUTILS}/{endpoint}?"
        f"{urllib.parse.urlencode(params)}"
    )

    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "GynandroCalling/1.0"
        },
    )

    try:
        with urllib.request.urlopen(
            request,
            timeout=180,
        ) as response:
            return response.read().decode("utf-8")

    except Exception as exc:
        raise RuntimeError(
            f"NCBI request failed: {endpoint}: {exc}"
        ) from exc


def sra_ids(term):
    result = json.loads(
        get(
            "esearch.fcgi",
            db="sra",
            term=term,
            retmax="100000",
            retmode="json",
        )
    )["esearchresult"]

    ids = result.get("idlist", [])
    count = int(result.get("count", 0))

    if count != len(ids):
        raise RuntimeError(
            "SRA search found "
            f"{count} records but returned {len(ids)}"
        )

    return ids


def runinfo_from_ids(ids):
    if not ids:
        return pd.DataFrame()

    frames = []

    for start in range(0, len(ids), 200):
        text = get(
            "efetch.fcgi",
            db="sra",
            id=",".join(ids[start:start + 200]),
            rettype="runinfo",
            retmode="text",
        )

        frame = pd.read_csv(
            io.StringIO(text),
            dtype=str,
            keep_default_na=False,
        )

        frames.append(frame)

    frame = pd.concat(
        frames,
        ignore_index=True,
    )

    if "Run" not in frame.columns:
        raise RuntimeError(
            "NCBI RunInfo has no Run column. "
            "Returned columns: "
            + ", ".join(frame.columns)
        )

    frame = frame.loc[
        frame["Run"].str.match(
            r"^(?:SRR|ERR|DRR)\d+$"
        )
    ].copy()

    return frame.drop_duplicates("Run")


def runinfo_from_accessions(accessions):
    frames = []

    for start in range(0, len(accessions), 25):
        batch = accessions[start:start + 25]

        term = " OR ".join(
            f"{run}[ACCN]"
            for run in batch
        )

        frames.append(
            runinfo_from_ids(
                sra_ids(term)
            )
        )

    return pd.concat(
        frames,
        ignore_index=True,
    ).drop_duplicates("Run")


def add_missing_columns(frame):
    columns = [
        "Run",
        "BioProject",
        "BioSample",
        "ScientificName",
        "SampleName",
        "Experiment",
        "LibraryName",
        "LibraryStrategy",
        "LibrarySource",
        "LibrarySelection",
        "LibraryLayout",
        "Platform",
        "Model",
        "spots",
        "bases",
        "size_MB",
    ]

    result = frame.copy()

    for column in columns:
        if column not in result.columns:
            result[column] = ""

    return result[columns].copy()


def wgs_mask(frame):
    return (
        frame["LibraryStrategy"]
        .str.upper()
        .eq("WGS")
        &
        frame["LibrarySource"]
        .str.upper()
        .eq("GENOMIC")
        &
        frame["Platform"]
        .str.upper()
        .eq("ILLUMINA")
        &
        frame["LibraryLayout"]
        .str.upper()
        .isin(["PAIRED", "SINGLE"])
    )


def evenly_spaced(frame, number):
    ordered = frame.sort_values(
        [
            "heterozygosity_among_callable",
            "sample",
        ],
        kind="stable",
    ).reset_index(drop=True)

    if len(ordered) <= number:
        return ordered.copy()

    indices = np.rint(
        np.linspace(
            0,
            len(ordered) - 1,
            number,
        )
    ).astype(int)

    if len(np.unique(indices)) != number:
        raise RuntimeError(
            "Could not construct unique, "
            "evenly spaced control ranks"
        )

    return ordered.iloc[indices].copy()


def main():
    OUT.mkdir(
        parents=True,
        exist_ok=True,
    )

    if not METRICS.is_file():
        raise FileNotFoundError(
            f"Missing control metrics: {METRICS}"
        )

    summary_text = get(
        "esummary.fcgi",
        db="bioproject",
        id=BIOPROJECT_UID,
        retmode="json",
    )

    projects = sorted(
        set(
            re.findall(
                r"PRJ(?:NA|EB|DB)\d+",
                summary_text,
            )
        )
    )

    if len(projects) != 1:
        raise RuntimeError(
            f"BioProject UID {BIOPROJECT_UID} "
            f"resolved to {projects}, not one project"
        )

    focal_project = projects[0]

    focal = add_missing_columns(
        runinfo_from_ids(
            sra_ids(
                f"{focal_project}[BioProject]"
            )
        )
    )

    if focal.empty:
        raise RuntimeError(
            f"No SRA runs found for {focal_project}"
        )

    searchable = focal.astype(str).agg(
        " | ".join,
        axis=1,
    )

    explicit_pattern = "|".join(
        re.escape(value)
        for value in FORBIDDEN
    )

    focal["forbidden_old_T2_T3"] = (
        searchable.str.contains(
            explicit_pattern,
            case=False,
            regex=True,
        )
        |
        searchable.str.contains(
            (
                r"(?:^|[^A-Za-z0-9])"
                r"T[23]"
                r"(?:[^A-Za-z0-9]|$)"
            ),
            case=False,
            regex=True,
        )
    )

    focal["ncbi_wgs_ok"] = wgs_mask(focal)

    focal["proposed_keep"] = (
        focal["ncbi_wgs_ok"]
        &
        ~focal["forbidden_old_T2_T3"]
    )

    focal["proposed_reason"] = np.where(
        focal["forbidden_old_T2_T3"],
        "locked_old_T2_T3_exclusion",
        np.where(
            ~focal["ncbi_wgs_ok"],
            "not_Illumina_genomic_WGS",
            "included",
        ),
    )

    metrics = pd.read_csv(
        METRICS,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    required = {
        "sample",
        "call_rate",
        "heterozygosity_among_callable",
    }

    missing = required.difference(
        metrics.columns
    )

    if missing:
        raise RuntimeError(
            "Missing control columns: "
            + ", ".join(sorted(missing))
        )

    metrics = metrics.loc[
        metrics["sample"].str.match(
            r"^(?:SRR|ERR|DRR)\d+$"
        )
    ].copy()

    metrics["call_rate"] = pd.to_numeric(
        metrics["call_rate"],
        errors="coerce",
    )

    metrics[
        "heterozygosity_among_callable"
    ] = pd.to_numeric(
        metrics[
            "heterozygosity_among_callable"
        ],
        errors="coerce",
    )

    metrics = metrics.dropna(
        subset=[
            "call_rate",
            "heterozygosity_among_callable",
        ]
    )

    info = add_missing_columns(
        runinfo_from_accessions(
            sorted(
                metrics["sample"].unique()
            )
        )
    )

    missing_runs = sorted(
        set(metrics["sample"])
        .difference(info["Run"])
    )

    if missing_runs:
        raise RuntimeError(
            "NCBI RunInfo missing: "
            + ", ".join(missing_runs)
        )

    controls = metrics.merge(
        info,
        left_on="sample",
        right_on="Run",
        how="left",
        validate="one_to_one",
    )

    controls["ncbi_wgs_ok"] = wgs_mask(
        controls
    )

    controls["melpomene_ok"] = (
        controls["ScientificName"]
        .str.contains(
            r"\bHeliconius melpomene\b",
            case=False,
            regex=True,
        )
    )

    controls["call_rate_ok"] = (
        controls["call_rate"]
        .ge(MIN_CALL_RATE)
    )

    controls["eligible"] = (
        controls["ncbi_wgs_ok"]
        &
        controls["melpomene_ok"]
        &
        controls["call_rate_ok"]
    )

    eligible = controls.loc[
        controls["eligible"]
    ].copy()

    controls.to_csv(
        OUT / "all_control_run_review.tsv",
        sep="\t",
        index=False,
    )

    if len(eligible) < 12:
        raise RuntimeError(
            f"Only {len(eligible)} H. melpomene "
            "controls passed filters; need at least 12. "
            "The complete review table was written."
        )

    selected = evenly_spaced(
        eligible,
        min(
            TARGET_CONTROLS,
            len(eligible),
        ),
    )

    selected[
        "provisional_het_stratum"
    ] = pd.qcut(
        selected[
            "heterozygosity_among_callable"
        ].rank(method="first"),
        3,
        labels=[
            "low",
            "middle",
            "high",
        ],
    ).astype(str)

    focal.to_csv(
        OUT / "focal_run_review.tsv",
        sep="\t",
        index=False,
    )

    selected.to_csv(
        OUT / "selected_control_run_review.tsv",
        sep="\t",
        index=False,
    )

    pd.DataFrame(
        [
            {
                "accession": REFERENCE,
                "taxon": (
                    "Heliconius melpomene rosina"
                ),
                "assembly_name": (
                    "Heliconius_melpomene_rosina.546"
                ),
            }
        ]
    ).to_csv(
        OUT / "reference_manifest.tsv",
        sep="\t",
        index=False,
    )

    audit = pd.DataFrame(
        [
            (
                "focal_bioproject_uid",
                BIOPROJECT_UID,
            ),
            (
                "focal_bioproject_accession",
                focal_project,
            ),
            (
                "focal_runs_total",
                len(focal),
            ),
            (
                "focal_WGS_runs_proposed",
                int(
                    focal[
                        "proposed_keep"
                    ].sum()
                ),
            ),
            (
                "old_T2_T3_runs_excluded",
                int(
                    focal[
                        "forbidden_old_T2_T3"
                    ].sum()
                ),
            ),
            (
                "control_candidates",
                len(controls),
            ),
            (
                "eligible_H_melpomene_controls",
                len(eligible),
            ),
            (
                "selected_controls",
                len(selected),
            ),
            (
                "reference_accession",
                REFERENCE,
            ),
            (
                "control_selection_uses_focal_data",
                False,
            ),
            (
                "optix_heterozygosity_is_provisional",
                True,
            ),
        ],
        columns=[
            "item",
            "value",
        ],
    )

    audit.to_csv(
        OUT / "inventory.audit.tsv",
        sep="\t",
        index=False,
    )

    print(
        "GENOME-WIDE INVENTORY COMPLETE"
    )

    print(
        audit.to_string(
            index=False
        )
    )

    print(
        "\nFOCAL RUN REVIEW"
    )

    print(
        focal[
            [
                "Run",
                "BioProject",
                "BioSample",
                "SampleName",
                "LibraryName",
                "LibraryStrategy",
                "LibrarySource",
                "LibraryLayout",
                "Platform",
                "bases",
                "proposed_keep",
                "proposed_reason",
            ]
        ].to_string(
            index=False
        )
    )

    print(
        "\nSELECTED CONTROL PROPOSAL"
    )

    print(
        selected[
            [
                "sample",
                "ScientificName",
                "SampleName",
                "LibraryLayout",
                "bases",
                "call_rate",
                (
                    "heterozygosity_"
                    "among_callable"
                ),
                "provisional_het_stratum",
            ]
        ].to_string(
            index=False
        )
    )


if __name__ == "__main__":
    try:
        main()

    except Exception as exc:
        print(
            f"ERROR: {exc}",
            file=sys.stderr,
        )
        raise
