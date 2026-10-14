#!/usr/bin/env python3

from pathlib import Path
import shutil

import pandas as pd


project = Path(
    "/scratch/martinlab/joe/GynandroCalling"
)

checkpoint_dir = (
    project
    / "analysis"
    / "drop_problematic_reanalysis"
    / "checkpoints"
)

manifest_file = (
    checkpoint_dir
    / "original502_filtered_cohort_manifest.WGS_vs_WGSplusCS.tsv"
)

backup_file = (
    checkpoint_dir
    / (
        "original502_filtered_cohort_manifest."
        "WGS_vs_WGSplusCS.before_remaining_breakers.tsv"
    )
)

exclusion_file = (
    checkpoint_dir
    / "excluded_comparative_samples.round2.tsv"
)

new_exclusions = {
    "ERR3650481": (
        "remaining breaker of best Dennis SNP "
        "Hmr18:1115258"
    ),
    "ERR3652402": (
        "remaining breaker of best Ray SNP "
        "Hmr18:1183019"
    ),
}


def as_bool(series):

    return (
        series.astype(str)
        .str.strip()
        .str.lower()
        .isin(
            {
                "true",
                "t",
                "1",
                "yes",
                "y",
            }
        )
    )


if not manifest_file.exists():
    raise RuntimeError(
        f"Missing corrected manifest: {manifest_file}"
    )

manifest = pd.read_csv(
    manifest_file,
    sep="\t",
    dtype=str,
    keep_default_na=False,
)

required = {
    "cohort",
    "module",
    "sample",
    "keep",
}

missing = required.difference(
    manifest.columns
)

if missing:
    raise RuntimeError(
        "Manifest is missing columns: "
        + ", ".join(
            sorted(missing)
        )
    )

before = (
    manifest.loc[
        as_bool(manifest["keep"])
    ]
    .groupby(
        [
            "cohort",
            "module",
        ]
    )["sample"]
    .nunique()
    .reset_index(
        name="retained_samples"
    )
)

if not backup_file.exists():
    shutil.copy2(
        manifest_file,
        backup_file,
    )

change_mask = (
    manifest["sample"].isin(
        new_exclusions
    )
    & as_bool(
        manifest["keep"]
    )
)

changed = manifest.loc[
    change_mask,
    [
        "cohort",
        "module",
        "sample",
        "taxon_name",
        "data_type",
        "keep",
    ],
].copy()

manifest.loc[
    change_mask,
    "keep",
] = "False"

temporary_file = manifest_file.with_suffix(
    ".tsv.part"
)

manifest.to_csv(
    temporary_file,
    sep="\t",
    index=False,
)

temporary_file.replace(
    manifest_file
)

after = (
    manifest.loc[
        as_bool(manifest["keep"])
    ]
    .groupby(
        [
            "cohort",
            "module",
        ]
    )["sample"]
    .nunique()
    .reset_index(
        name="retained_samples"
    )
)

all_exclusions = [
    {
        "sample": "ERR3653118",
        "reason": "previously specified problematic meriana",
    },
    {
        "sample": "ERR3653120",
        "reason": "previously specified problematic meriana",
    },
    {
        "sample": "ERR3653298",
        "reason": "previously specified problematic contigua",
    },
    {
        "sample": "ERR3653306",
        "reason": "previously specified problematic contigua",
    },
    {
        "sample": "ERR3653308",
        "reason": "contigua specimen coded B0D0R0",
    },
    {
        "sample": "ERR3653309",
        "reason": "previously specified problematic contigua",
    },
]

for sample, reason in new_exclusions.items():
    all_exclusions.append(
        {
            "sample": sample,
            "reason": reason,
        }
    )

pd.DataFrame(
    all_exclusions
).to_csv(
    exclusion_file,
    sep="\t",
    index=False,
)

print("ROUND-2 COMPARATIVE EXCLUSIONS")
print()
print("NEWLY CHANGED MANIFEST RECORDS")

if len(changed):
    print(
        changed.to_string(
            index=False
        )
    )
else:
    print(
        "No keep=True records changed; "
        "the samples may already be excluded."
    )

print()
print("RETAINED COUNTS BEFORE")
print(
    before.to_string(
        index=False
    )
)

print()
print("RETAINED COUNTS AFTER")
print(
    after.to_string(
        index=False
    )
)

print()
print(f"Updated manifest: {manifest_file}")
print(f"Backup manifest: {backup_file}")
print(f"Exclusion record: {exclusion_file}")
