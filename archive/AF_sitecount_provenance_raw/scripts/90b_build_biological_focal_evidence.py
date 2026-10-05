#!/usr/bin/env python3

from pathlib import Path
import math
import re
import sys

import numpy as np
import pandas as pd


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

CK = PROJECT / "analysis" / "original_502" / "checkpoints"

EXPECTATION_FILE = CK / "focal_module_expectations.tsv"

# ------------------------------------------------------------
# OPTIX
#
# Site ascertainment:
#   comparative 51-sample stable candidate table ONLY
#
# Focal evaluation:
#   original-502 calls at the 51-sample candidate positions
#
# These two layers remain deliberately separate.
# ------------------------------------------------------------

OPTIX_COMPARATIVE = (
    PROJECT
    / "analysis"
    / "checkpoints"
    / "optix_51sample_stable_candidate_sites.tsv"
)

OPTIX_COMPARATIVE_SCORES = (
    PROJECT
    / "analysis"
    / "checkpoints"
    / "optix_51sample_comparative_scores.tsv"
)

OPTIX_FOCAL_BRIDGE = (
    CK
    / "optix_original_502_comparative_scores.tsv"
)

# ------------------------------------------------------------
# IVORY / Yb
# ------------------------------------------------------------

YB_COMPARATIVE_SELECTED = (
    CK
    / "ivory_comparative_informative_sites_BA80.tsv"
)

YB_COMPARATIVE_SCORES = (
    CK
    / "ivory_original_502_comparative_scores.tsv"
)

YB_FOCAL = (
    CK
    / "ivory_original_502_comparative_scores.tsv"
)

OUT_SITES = CK / "frozen_comparative_sites_biological_polarity.tsv"

OUT_SELECTION_AUDIT = (
    CK
    / "comparative_site_selection_audit.tsv"
)

OUT_WIDE = (
    CK
    / "focal_calls_at_frozen_comparative_sites_biological.wide.tsv"
)

OUT_LONG = (
    CK
    / "focal_calls_at_frozen_comparative_sites_biological.long.tsv"
)

OUT_SUMMARY = (
    CK
    / "focal_calls_at_frozen_comparative_sites_biological.summary.tsv"
)

OUT_SPANS = (
    CK
    / "focal_module_locus_spans.tsv"
)

OUT_SOURCE_SUMMARY = (
    CK
    / "frozen_comparative_site_source_summary.tsv"
)


MODULE_ORDER = [
    "Band",
    "Ray1",
    "Ray2",
    "Dennis1",
    "Dennis2",
    "Yb",
]


BLOCK_MAP = {
    "Band": "B",
    "Dennis1": "D",
    "Dennis2": "D",
    "Ray1": "R",
    "Ray2": "R",
    "Yb": "Yb",
}


DISPLAY_MAP = {
    "Band": "Band",
    "Dennis1": "Dennis 1",
    "Dennis2": "Dennis 2",
    "Ray1": "Ray 1",
    "Ray2": "Ray 2",
    "Yb": "Yb",
}


OPTIX_POSITION_INTERVALS = {
    "Band": (1111500, 1117000),
    "Ray1": (1156000, 1161000),
    "Ray2": (1174000, 1189000),
    "Dennis1": (1215000, 1217000),
    "Dennis2": (1218000, 1223000),
}


def fail(message):
    raise RuntimeError(message)


def read_tsv(path):
    if not path.exists():
        fail(f"Required input is missing: {path}")

    if path.stat().st_size == 0:
        fail(f"Required input is empty: {path}")

    return pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        na_values=[],
    )


def emptyish(value):
    if value is None:
        return True

    text = str(value).strip()

    return text in {
        "",
        ".",
        "NA",
        "NaN",
        "nan",
        "None",
        "none",
    }


def as_float(value):
    if emptyish(value):
        return np.nan

    try:
        return float(value)
    except Exception:
        return np.nan


def as_int(value):
    x = as_float(value)

    if not np.isfinite(x):
        return None

    return int(round(x))


def truthy(value):
    if emptyish(value):
        return False

    return str(value).strip().lower() in {
        "1",
        "true",
        "t",
        "yes",
        "y",
        "pass",
        "passed",
    }


def normalise_allele(value):
    if emptyish(value):
        return None

    return str(value).strip().upper()


def module_from_values(panel="", trait="", pos=None):
    text = " ".join(
        [
            str(panel).lower(),
            str(trait).lower(),
        ]
    )

    compact = re.sub(r"[^a-z0-9]+", "", text)

    if "band" in compact:
        return "Band"

    if (
        "dennis1" in compact
        or "dennisinterval1" in compact
        or "d1" == compact
    ):
        return "Dennis1"

    if (
        "dennis2" in compact
        or "dennisinterval2" in compact
        or "d2" == compact
    ):
        return "Dennis2"

    if (
        "ray1" in compact
        or "rayinterval1" in compact
        or "r1" == compact
    ):
        return "Ray1"

    if (
        "ray2" in compact
        or "rayinterval2" in compact
        or "r2" == compact
    ):
        return "Ray2"

    p = as_int(pos)

    if p is not None:
        for module, (lo, hi) in OPTIX_POSITION_INTERVALS.items():
            if lo <= p <= hi:
                return module

    return None


def choose_orientation(row):
    for column in [
        "stable_orientation",
        "selected_orientation",
        "orientation",
    ]:
        if column not in row.index:
            continue

        value = row[column]

        if emptyish(value):
            continue

        upper = str(value).upper()

        has_ref = "REF" in upper
        has_alt = "ALT" in upper

        if has_ref ^ has_alt:
            return str(value)

    return None


def presence_from_orientation(orientation, ref, alt):
    if orientation is None:
        return None, None

    upper = str(orientation).upper()

    ref = normalise_allele(ref)
    alt = normalise_allele(alt)

    has_ref = "REF" in upper
    has_alt = "ALT" in upper

    if has_ref and not has_alt:
        return ref, alt

    if has_alt and not has_ref:
        return alt, ref

    return None, None


def site_key(chrom, pos, ref, alt):
    return (
        str(chrom),
        int(pos),
        normalise_allele(ref),
        normalise_allele(alt),
    )


def position_key(chrom, pos):
    return (
        str(chrom),
        int(pos),
    )


def parse_gt(gt):
    if emptyish(gt):
        return None

    text = str(gt).strip()

    if text in {
        "./.",
        ".|.",
        ".",
    }:
        return None

    parts = re.split(r"[/|]", text)

    if len(parts) != 2:
        return None

    if any(x not in {"0", "1"} for x in parts):
        return None

    return sum(int(x) for x in parts)


def parse_ad(ad):
    if emptyish(ad):
        return (np.nan, np.nan)

    text = str(ad).strip()

    text = (
        text
        .replace("[", "")
        .replace("]", "")
        .replace("(", "")
        .replace(")", "")
    )

    parts = text.split(",")

    if len(parts) < 2:
        return (np.nan, np.nan)

    try:
        ref_depth = float(parts[0])
        alt_depth = float(parts[1])
    except Exception:
        return (np.nan, np.nan)

    return (ref_depth, alt_depth)


def biological_genotype_from_presence_copies(copies):
    if copies is None:
        return None

    if copies == 0:
        return "AA"

    if copies == 1:
        return "PA"

    if copies == 2:
        return "PP"

    return None


def genotype_to_presence_dosage(genotype):
    return {
        "AA": 0.0,
        "PA": 0.5,
        "PP": 1.0,
    }.get(genotype, np.nan)


def make_selected_optix(expectations):
    df = read_tsv(OPTIX_COMPARATIVE).copy()

    if "minimum_BA" not in df.columns:
        fail(
            "Optix stable candidate table lacks minimum_BA: "
            f"{OPTIX_COMPARATIVE}"
        )

    df["module"] = [
        module_from_values(
            panel=row.get("panel", ""),
            trait=row.get("trait", ""),
            pos=row.get("pos", ""),
        )
        for _, row in df.iterrows()
    ]

    df["minimum_comparative_BA"] = (
        pd.to_numeric(
            df["minimum_BA"],
            errors="coerce",
        )
    )

    selected = df.loc[
        df["module"].isin(MODULE_ORDER)
        & (
            df["minimum_comparative_BA"]
            >= 0.80
        )
    ].copy()

    if selected.empty:
        fail(
            "No Optix comparative sites survived minimum_BA >= 0.80."
        )

    selected["selected_orientation"] = [
        choose_orientation(row)
        for _, row in selected.iterrows()
    ]

    presence = [
        presence_from_orientation(
            row["selected_orientation"],
            row["ref"],
            row["alt"],
        )
        for _, row in selected.iterrows()
    ]

    selected["presence_associated_allele"] = [
        x[0]
        for x in presence
    ]

    selected["absence_associated_allele"] = [
        x[1]
        for x in presence
    ]

    bad = selected[
        selected["presence_associated_allele"].isna()
    ]

    if not bad.empty:
        print(
            "\nERROR: could not polarize these selected Optix rows:",
            file=sys.stderr,
        )

        print(
            bad[
                [
                    "panel",
                    "trait",
                    "chrom",
                    "pos",
                    "ref",
                    "alt",
                    "selected_orientation",
                ]
            ].to_string(index=False),
            file=sys.stderr,
        )

        fail("Selected Optix sites contain unresolved allele orientation.")

    selected["comparative_source"] = (
        "optix_51sample_comparative"
    )

    selected["comparative_source_file"] = str(
        OPTIX_COMPARATIVE
    )

    selected["inheritance_model"] = (
        "dominant_presence"
    )

    selected["selection_rule"] = (
        "51sample_stable_comparative_minimum_BA_ge_0.80"
    )

    return selected


def collapse_yb_selected(expectations):
    raw = read_tsv(YB_COMPARATIVE_SELECTED).copy()

    raw["module"] = "Yb"

    if "minimum_BA" not in raw.columns:
        fail(
            "Yb BA80 table lacks minimum_BA: "
            f"{YB_COMPARATIVE_SELECTED}"
        )

    raw["minimum_comparative_BA"] = (
        pd.to_numeric(
            raw["minimum_BA"],
            errors="coerce",
        )
    )

    raw = raw.loc[
        raw["minimum_comparative_BA"] >= 0.80
    ].copy()

    raw["selected_orientation"] = [
        choose_orientation(row)
        for _, row in raw.iterrows()
    ]

    grouped_rows = []

    group_columns = [
        "chrom",
        "pos",
        "ref",
        "alt",
    ]

    for key, group in raw.groupby(
        group_columns,
        dropna=False,
        sort=True,
    ):
        orientations = sorted(
            {
                str(x)
                for x in group["selected_orientation"]
                if not emptyish(x)
            }
        )

        presence_alleles = set()

        absence_alleles = set()

        for _, row in group.iterrows():
            p, a = presence_from_orientation(
                row["selected_orientation"],
                row["ref"],
                row["alt"],
            )

            if p is not None:
                presence_alleles.add(p)

            if a is not None:
                absence_alleles.add(a)

        if len(presence_alleles) != 1:
            fail(
                "Yb biological polarity disagrees across phenotype "
                f"rows at {key}: {sorted(presence_alleles)}"
            )

        if len(absence_alleles) != 1:
            fail(
                "Yb absence polarity disagrees across phenotype "
                f"rows at {key}: {sorted(absence_alleles)}"
            )

        traits = sorted(
            {
                str(x)
                for x in group.get(
                    "trait",
                    pd.Series(dtype=str),
                )
                if not emptyish(x)
            }
        )

        grouped_rows.append(
            {
                "module": "Yb",
                "panel": "ivory",
                "trait": ",".join(traits),
                "chrom": key[0],
                "pos": key[1],
                "ref": key[2],
                "alt": key[3],
                "selected_orientation": "|".join(
                    orientations
                ),
                "minimum_comparative_BA": float(
                    group[
                        "minimum_comparative_BA"
                    ].max()
                ),
                "presence_associated_allele": next(
                    iter(presence_alleles)
                ),
                "absence_associated_allele": next(
                    iter(absence_alleles)
                ),
                "comparative_source":
                    "ivory_original502_comparative",
                "comparative_source_file": str(
                    YB_COMPARATIVE_SELECTED
                ),
                "inheritance_model":
                    "recessive_presence",
                "selection_rule":
                    "ivory_comparative_minimum_BA_ge_0.80",
            }
        )

    return pd.DataFrame(grouped_rows)


def merge_expectations(selected, expectations):
    keep = expectations[
        [
            "module",
            "left_phenotype",
            "right_phenotype",
            "left_expected_genotype",
            "right_expected_genotype",
            "left_expected_presence_dosage",
            "right_expected_presence_dosage",
            "reference_biological_genotype",
            "reference_presence_dosage",
            "reference_trait_state",
        ]
    ].copy()

    result = selected.merge(
        keep,
        on="module",
        how="left",
        validate="many_to_one",
    )

    if result[
        "left_expected_genotype"
    ].isna().any():
        fail("Expectation-table merge failed for one or more modules.")

    result["block"] = result["module"].map(
        BLOCK_MAP
    )

    result["module_display"] = result["module"].map(
        DISPLAY_MAP
    )

    return result


def deduplicate_selected_sites(selected):
    rows = []

    for (
        module,
        chrom,
        pos,
        ref,
        alt,
    ), group in selected.groupby(
        [
            "module",
            "chrom",
            "pos",
            "ref",
            "alt",
        ],
        dropna=False,
        sort=True,
    ):
        polarities = {
            (
                normalise_allele(x),
                normalise_allele(y),
            )
            for x, y in zip(
                group["presence_associated_allele"],
                group["absence_associated_allele"],
            )
        }

        if len(polarities) != 1:
            fail(
                "Conflicting biological polarity for selected site: "
                f"{module} {chrom}:{pos} {ref}>{alt}"
            )

        row = group.iloc[0].copy()

        row["minimum_comparative_BA"] = float(
            pd.to_numeric(
                group["minimum_comparative_BA"],
                errors="coerce",
            ).max()
        )

        if "trait" in group.columns:
            row["trait"] = ",".join(
                sorted(
                    {
                        str(x)
                        for x in group["trait"]
                        if not emptyish(x)
                    }
                )
            )

        if "panel" in group.columns:
            row["panel"] = ",".join(
                sorted(
                    {
                        str(x)
                        for x in group["panel"]
                        if not emptyish(x)
                    }
                )
            )

        rows.append(row)

    out = pd.DataFrame(rows)

    order_map = {
        module: i
        for i, module in enumerate(
            MODULE_ORDER
        )
    }

    out["_module_order"] = out["module"].map(
        order_map
    )

    out["pos_numeric"] = pd.to_numeric(
        out["pos"],
        errors="raise",
    ).astype(int)

    out = out.sort_values(
        [
            "_module_order",
            "chrom",
            "pos_numeric",
        ]
    ).reset_index(drop=True)

    out["site_order_within_module"] = (
        out.groupby("module").cumcount() + 1
    )

    out["site_id"] = [
        (
            f"{module}:"
            f"{chrom}:"
            f"{int(pos)}:"
            f"{ref}>{alt}"
        )
        for (
            module,
            chrom,
            pos,
            ref,
            alt,
        ) in zip(
            out["module"],
            out["chrom"],
            out["pos_numeric"],
            out["ref"],
            out["alt"],
        )
    ]

    out = out.drop(
        columns=[
            "_module_order",
        ]
    )

    return out


def collapse_ivory_score_trace(selected_site_keys):
    scores = read_tsv(
        YB_COMPARATIVE_SCORES
    ).copy()

    if "balanced_accuracy" not in scores.columns:
        fail(
            "Yb comparative-score table lacks balanced_accuracy."
        )

    scores["selection_BA_value"] = pd.to_numeric(
        scores["balanced_accuracy"],
        errors="coerce",
    )

    rows = []

    group_cols = [
        "trait",
        "chrom",
        "pos",
        "ref",
        "alt",
    ]

    for key, group in scores.groupby(
        group_cols,
        dropna=False,
        sort=True,
    ):
        ba = group[
            "selection_BA_value"
        ].dropna()

        if ba.empty:
            continue

        chrom = key[1]
        pos = as_int(key[2])

        if pos is None:
            continue

        ref = key[3]
        alt = key[4]

        selected = (
            site_key(
                chrom,
                pos,
                ref,
                alt,
            )
            in selected_site_keys
        )

        rows.append(
            {
                "module": "Yb",
                "panel": "ivory",
                "trait": key[0],
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "selection_BA": float(
                    ba.min()
                ),
                "selected": selected,
                "comparative_source":
                    "ivory_original502_comparative",
            }
        )

    return pd.DataFrame(rows)


def make_selection_audit(
    optix_all,
    selected_optix,
    selected_yb,
):
    selected_optix_keys = {
        (
            row["module"],
            site_key(
                row["chrom"],
                int(row["pos_numeric"]),
                row["ref"],
                row["alt"],
            ),
        )
        for _, row in selected_optix.iterrows()
    }

    optix = optix_all.copy()

    optix["module"] = [
        module_from_values(
            panel=row.get("panel", ""),
            trait=row.get("trait", ""),
            pos=row.get("pos", ""),
        )
        for _, row in optix.iterrows()
    ]

    optix = optix.loc[
        optix["module"].isin(
            MODULE_ORDER
        )
    ].copy()

    optix["pos"] = pd.to_numeric(
        optix["pos"],
        errors="coerce",
    )

    optix["selection_BA"] = pd.to_numeric(
        optix["minimum_BA"],
        errors="coerce",
    )

    optix = optix.loc[
        optix["pos"].notna()
        & optix["selection_BA"].notna()
    ].copy()

    optix["pos"] = optix["pos"].astype(int)

    optix["selected"] = [
        (
            module,
            site_key(
                chrom,
                int(pos),
                ref,
                alt,
            ),
        )
        in selected_optix_keys
        for (
            module,
            chrom,
            pos,
            ref,
            alt,
        ) in zip(
            optix["module"],
            optix["chrom"],
            optix["pos"],
            optix["ref"],
            optix["alt"],
        )
    ]

    optix["comparative_source"] = (
        "optix_51sample_comparative"
    )

    selected_yb_keys = {
        site_key(
            row["chrom"],
            int(row["pos_numeric"]),
            row["ref"],
            row["alt"],
        )
        for _, row in selected_yb.iterrows()
    }

    yb = collapse_ivory_score_trace(
        selected_yb_keys
    )

    columns = [
        "module",
        "panel",
        "trait",
        "chrom",
        "pos",
        "ref",
        "alt",
        "selection_BA",
        "selected",
        "comparative_source",
    ]

    result = pd.concat(
        [
            optix.reindex(columns=columns),
            yb.reindex(columns=columns),
        ],
        ignore_index=True,
    )

    order_map = {
        module: i
        for i, module in enumerate(
            MODULE_ORDER
        )
    }

    result["_module_order"] = result[
        "module"
    ].map(order_map)

    result = result.sort_values(
        [
            "_module_order",
            "pos",
        ]
    ).drop(
        columns="_module_order"
    )

    return result


def make_locus_spans(selection_audit):
    rows = []

    for module in MODULE_ORDER:
        d = selection_audit.loc[
            selection_audit["module"]
            == module
        ].copy()

        if d.empty:
            continue

        positions = pd.to_numeric(
            d["pos"],
            errors="coerce",
        ).dropna()

        if positions.empty:
            continue

        start = int(positions.min())
        end = int(positions.max())

        rows.append(
            {
                "module": module,
                "module_display":
                    DISPLAY_MAP[module],
                "block":
                    BLOCK_MAP[module],
                "chrom":
                    str(
                        d[
                            "chrom"
                        ].dropna().iloc[0]
                    ),
                "start": start,
                "end": end,
                "span_bp": end - start + 1,
            }
        )

    return pd.DataFrame(rows)


def choose_focal_row(
    selected_row,
    focal_df,
):
    module = selected_row["module"]
    chrom = str(selected_row["chrom"])
    pos = int(selected_row["pos_numeric"])

    d = focal_df.loc[
        (
            focal_df["_module"]
            == module
        )
        & (
            focal_df["_chrom"]
            == chrom
        )
        & (
            focal_df["_pos"]
            == pos
        )
    ].copy()

    if d.empty:
        return None

    selected_ref = normalise_allele(
        selected_row["ref"]
    )

    selected_alt = normalise_allele(
        selected_row["alt"]
    )

    exact = d.loc[
        (
            d["_ref"]
            == selected_ref
        )
        & (
            d["_alt"]
            == selected_alt
        )
    ]

    if not exact.empty:
        return exact.iloc[0]

    swapped = d.loc[
        (
            d["_ref"]
            == selected_alt
        )
        & (
            d["_alt"]
            == selected_ref
        )
    ]

    if not swapped.empty:
        return swapped.iloc[0]

    return d.iloc[0]


def prepare_optix_focal():
    df = read_tsv(
        OPTIX_FOCAL_BRIDGE
    ).copy()

    df["_module"] = [
        module_from_values(
            panel=row.get("panel", ""),
            trait=row.get("trait", ""),
            pos=row.get("pos", ""),
        )
        for _, row in df.iterrows()
    ]

    df["_chrom"] = df[
        "chrom"
    ].astype(str)

    df["_pos"] = pd.to_numeric(
        df["pos"],
        errors="coerce",
    )

    df = df.loc[
        df["_pos"].notna()
    ].copy()

    df["_pos"] = df["_pos"].astype(int)

    df["_ref"] = df[
        "ref"
    ].map(normalise_allele)

    df["_alt"] = df[
        "alt"
    ].map(normalise_allele)

    return df


def prepare_yb_focal():
    df = read_tsv(
        YB_FOCAL
    ).copy()

    df["_module"] = "Yb"

    df["_chrom"] = df[
        "chrom"
    ].astype(str)

    df["_pos"] = pd.to_numeric(
        df["pos"],
        errors="coerce",
    )

    df = df.loc[
        df["_pos"].notna()
    ].copy()

    df["_pos"] = df["_pos"].astype(int)

    df["_ref"] = df[
        "ref"
    ].map(normalise_allele)

    df["_alt"] = df[
        "alt"
    ].map(normalise_allele)

    return df


def focal_ref_alt(row, selected_row):
    if row is None:
        return (
            normalise_allele(
                selected_row["ref"]
            ),
            normalise_allele(
                selected_row["alt"]
            ),
        )

    for ref_col, alt_col in [
        (
            "original_ref",
            "original_alt",
        ),
        (
            "ref",
            "alt",
        ),
    ]:
        if (
            ref_col in row.index
            and alt_col in row.index
            and not emptyish(row[ref_col])
            and not emptyish(row[alt_col])
        ):
            return (
                normalise_allele(
                    row[ref_col]
                ),
                normalise_allele(
                    row[alt_col]
                ),
            )

    return (
        normalise_allele(
            selected_row["ref"]
        ),
        normalise_allele(
            selected_row["alt"]
        ),
    )


def get_field(row, candidates):
    if row is None:
        return ""

    for column in candidates:
        if column in row.index:
            value = row[column]

            if not emptyish(value):
                return value

    return ""


def evaluate_side(
    selected_row,
    focal_row,
    side,
):
    prefix = (
        "TN_Left"
        if side == "Left"
        else "TN_Right"
    )

    gt = get_field(
        focal_row,
        [
            f"{prefix}_GT",
        ],
    )

    ad = get_field(
        focal_row,
        [
            f"{prefix}_AD",
        ],
    )

    dp = as_float(
        get_field(
            focal_row,
            [
                f"{prefix}_DP",
            ],
        )
    )

    gq = as_float(
        get_field(
            focal_row,
            [
                f"{prefix}_GQ",
            ],
        )
    )

    raw_alt_fraction = as_float(
        get_field(
            focal_row,
            [
                f"{prefix}_ALT_fraction",
            ],
        )
    )

    focal_ref, focal_alt = focal_ref_alt(
        focal_row,
        selected_row,
    )

    presence_allele = normalise_allele(
        selected_row[
            "presence_associated_allele"
        ]
    )

    alt_copies = parse_gt(gt)

    mapping_status = "unresolved"

    presence_copies = None

    if (
        alt_copies is not None
        and presence_allele
        == focal_alt
    ):
        presence_copies = alt_copies
        mapping_status = (
            "presence_is_focal_ALT"
        )

    elif (
        alt_copies is not None
        and presence_allele
        == focal_ref
    ):
        presence_copies = (
            2 - alt_copies
        )

        mapping_status = (
            "presence_is_focal_REF"
        )

    elif presence_allele == focal_alt:
        mapping_status = (
            "presence_is_focal_ALT"
        )

    elif presence_allele == focal_ref:
        mapping_status = (
            "presence_is_focal_REF"
        )

    ref_depth, alt_depth = parse_ad(ad)

    presence_fraction = np.nan

    if (
        np.isfinite(ref_depth)
        and np.isfinite(alt_depth)
        and (
            ref_depth + alt_depth
            > 0
        )
    ):
        if presence_allele == focal_alt:
            presence_fraction = (
                alt_depth
                / (
                    ref_depth
                    + alt_depth
                )
            )

        elif presence_allele == focal_ref:
            presence_fraction = (
                ref_depth
                / (
                    ref_depth
                    + alt_depth
                )
            )

    genotype = (
        biological_genotype_from_presence_copies(
            presence_copies
        )
    )

    genotype_dosage = (
        genotype_to_presence_dosage(
            genotype
        )
    )

    qc_pass = bool(
        genotype is not None
        and np.isfinite(dp)
        and dp >= 2
        and np.isfinite(gq)
        and gq >= 20
        and mapping_status
        != "unresolved"
    )

    expected_genotype = str(
        selected_row[
            (
                "left_expected_genotype"
                if side == "Left"
                else "right_expected_genotype"
            )
        ]
    )

    expected_dosage = float(
        selected_row[
            (
                "left_expected_presence_dosage"
                if side == "Left"
                else "right_expected_presence_dosage"
            )
        ]
    )

    match = (
        bool(
            qc_pass
            and genotype
            == expected_genotype
        )
    )

    return {
        "GT_raw": gt,
        "AD_raw": ad,
        "DP": dp,
        "GQ": gq,
        "ALT_fraction_raw":
            raw_alt_fraction,
        "focal_ref": focal_ref,
        "focal_alt": focal_alt,
        "allele_mapping_status":
            mapping_status,
        "ref_depth": ref_depth,
        "alt_depth": alt_depth,
        "biological_genotype":
            genotype
            if qc_pass
            else "missing",
        "presence_dosage":
            genotype_dosage
            if qc_pass
            else np.nan,
        "presence_allele_fraction":
            presence_fraction
            if qc_pass
            else np.nan,
        "QC_pass_DP2_GQ20":
            qc_pass,
        "expected_genotype":
            expected_genotype,
        "expected_presence_dosage":
            expected_dosage,
        "matches_expectation":
            match
            if qc_pass
            else np.nan,
    }


def build_focal_tables(selected):
    optix_focal = prepare_optix_focal()

    yb_focal = prepare_yb_focal()

    wide_rows = []
    long_rows = []

    for _, site in selected.iterrows():
        if site["module"] == "Yb":
            focal = choose_focal_row(
                site,
                yb_focal,
            )
        else:
            focal = choose_focal_row(
                site,
                optix_focal,
            )

        left = evaluate_side(
            site,
            focal,
            "Left",
        )

        right = evaluate_side(
            site,
            focal,
            "Right",
        )

        base = {
            "site_id":
                site["site_id"],
            "module":
                site["module"],
            "module_display":
                site["module_display"],
            "block":
                site["block"],
            "site_order_within_module":
                int(
                    site[
                        "site_order_within_module"
                    ]
                ),
            "chrom":
                site["chrom"],
            "pos":
                int(
                    site["pos_numeric"]
                ),
            "ref":
                site["ref"],
            "alt":
                site["alt"],
            "selected_orientation":
                site[
                    "selected_orientation"
                ],
            "presence_associated_allele":
                site[
                    "presence_associated_allele"
                ],
            "absence_associated_allele":
                site[
                    "absence_associated_allele"
                ],
            "minimum_comparative_BA":
                float(
                    site[
                        "minimum_comparative_BA"
                    ]
                ),
            "comparative_source":
                site[
                    "comparative_source"
                ],
            "selection_rule":
                site[
                    "selection_rule"
                ],
            "inheritance_model":
                site[
                    "inheritance_model"
                ],
            "reference_biological_genotype":
                site[
                    "reference_biological_genotype"
                ],
            "reference_trait_state":
                site[
                    "reference_trait_state"
                ],
        }

        wide = dict(base)

        for side_label, result in [
            ("left", left),
            ("right", right),
        ]:
            for key, value in result.items():
                wide[
                    f"{side_label}_{key}"
                ] = value

        both_callable = (
            bool(
                left["QC_pass_DP2_GQ20"]
            )
            and bool(
                right["QC_pass_DP2_GQ20"]
            )
        )

        pair_match = (
            bool(
                left["matches_expectation"]
            )
            and bool(
                right["matches_expectation"]
            )
            if both_callable
            else np.nan
        )

        if both_callable:
            total_distance = (
                abs(
                    left[
                        "presence_dosage"
                    ]
                    - left[
                        "expected_presence_dosage"
                    ]
                )
                + abs(
                    right[
                        "presence_dosage"
                    ]
                    - right[
                        "expected_presence_dosage"
                    ]
                )
            )
        else:
            total_distance = np.nan

        wide["both_callable"] = (
            both_callable
        )

        wide["pair_matches_expectation"] = (
            pair_match
        )

        wide[
            "pair_genotype_distance_from_expectation"
        ] = total_distance

        wide_rows.append(wide)

        for side_name, result in [
            ("Left", left),
            ("Right", right),
        ]:
            long = dict(base)

            long["side"] = side_name

            for key, value in result.items():
                long[key] = value

            long_rows.append(long)

    return (
        pd.DataFrame(wide_rows),
        pd.DataFrame(long_rows),
    )


def summarise(wide):
    rows = []

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

        both = (
            d["both_callable"].astype(bool)
        )

        left_match = d.loc[
            left_callable,
            "left_matches_expectation",
        ]

        right_match = d.loc[
            right_callable,
            "right_matches_expectation",
        ]

        pair_match = d.loc[
            both,
            "pair_matches_expectation",
        ]

        rows.append(
            {
                "module": module,
                "module_display":
                    DISPLAY_MAP[module],
                "independent_comparative_sites":
                    len(d),
                "left_callable":
                    int(left_callable.sum()),
                "right_callable":
                    int(right_callable.sum()),
                "both_callable":
                    int(both.sum()),
                "left_exact_expected_n":
                    int(
                        pd.to_numeric(
                            left_match,
                            errors="coerce",
                        ).fillna(0).sum()
                    ),
                "right_exact_expected_n":
                    int(
                        pd.to_numeric(
                            right_match,
                            errors="coerce",
                        ).fillna(0).sum()
                    ),
                "pair_exact_expected_n":
                    int(
                        pd.to_numeric(
                            pair_match,
                            errors="coerce",
                        ).fillna(0).sum()
                    ),
                "pair_exact_expected_fraction":
                    (
                        float(
                            pd.to_numeric(
                                pair_match,
                                errors="coerce",
                            ).mean()
                        )
                        if len(pair_match)
                        else np.nan
                    ),
                "left_expected_genotype":
                    d[
                        "left_expected_genotype"
                    ].iloc[0],
                "right_expected_genotype":
                    d[
                        "right_expected_genotype"
                    ].iloc[0],
                "left_expected_presence_dosage":
                    d[
                        "left_expected_presence_dosage"
                    ].iloc[0],
                "right_expected_presence_dosage":
                    d[
                        "right_expected_presence_dosage"
                    ].iloc[0],
                "left_median_presence_allele_fraction":
                    pd.to_numeric(
                        d[
                            "left_presence_allele_fraction"
                        ],
                        errors="coerce",
                    ).median(),
                "right_median_presence_allele_fraction":
                    pd.to_numeric(
                        d[
                            "right_presence_allele_fraction"
                        ],
                        errors="coerce",
                    ).median(),
            }
        )

    return pd.DataFrame(rows)


def source_summary(selected):
    return (
        selected.groupby(
            [
                "comparative_source",
                "module",
            ],
            dropna=False,
        )
        .agg(
            sites=(
                "site_id",
                "nunique",
            ),
            minimum_BA_min=(
                "minimum_comparative_BA",
                "min",
            ),
            minimum_BA_median=(
                "minimum_comparative_BA",
                "median",
            ),
            minimum_BA_max=(
                "minimum_comparative_BA",
                "max",
            ),
        )
        .reset_index()
    )


def main():
    CK.mkdir(
        parents=True,
        exist_ok=True,
    )

    expectations = read_tsv(
        EXPECTATION_FILE
    )

    numeric_expectation_cols = [
        "left_expected_presence_dosage",
        "right_expected_presence_dosage",
        "reference_presence_dosage",
    ]

    for column in numeric_expectation_cols:
        expectations[column] = (
            pd.to_numeric(
                expectations[column],
                errors="raise",
            )
        )

    optix_all = read_tsv(
        OPTIX_COMPARATIVE
    )

    selected_optix = (
        make_selected_optix(
            expectations
        )
    )

    selected_yb = (
        collapse_yb_selected(
            expectations
        )
    )

    selected = pd.concat(
        [
            selected_optix,
            selected_yb,
        ],
        ignore_index=True,
        sort=False,
    )

    selected = merge_expectations(
        selected,
        expectations,
    )

    selected = (
        deduplicate_selected_sites(
            selected
        )
    )

    counts = (
        selected.groupby("module")
        .size()
        .reindex(
            MODULE_ORDER,
            fill_value=0,
        )
    )

    # --------------------------------------------------------
    # Hard guard: no more silently producing a Band n=0 plot.
    # --------------------------------------------------------

    missing_modules = [
        module
        for module, n in counts.items()
        if int(n) == 0
    ]

    if missing_modules:
        print(
            "\nSELECTED SITE COUNTS\n",
            counts.to_string(),
            sep="",
            file=sys.stderr,
        )

        fail(
            "No frozen comparative sites for: "
            + ", ".join(
                missing_modules
            )
        )

    selection_audit = (
        make_selection_audit(
            optix_all,
            selected.loc[
                selected[
                    "module"
                ] != "Yb"
            ].copy(),
            selected.loc[
                selected[
                    "module"
                ] == "Yb"
            ].copy(),
        )
    )

    spans = make_locus_spans(
        selection_audit
    )

    wide, long = build_focal_tables(
        selected
    )

    summary = summarise(wide)

    sources = source_summary(
        selected
    )

    site_columns = [
        "site_id",
        "module",
        "module_display",
        "block",
        "site_order_within_module",
        "panel",
        "trait",
        "chrom",
        "pos_numeric",
        "ref",
        "alt",
        "selected_orientation",
        "presence_associated_allele",
        "absence_associated_allele",
        "minimum_comparative_BA",
        "inheritance_model",
        "comparative_source",
        "comparative_source_file",
        "selection_rule",
        "left_phenotype",
        "right_phenotype",
        "left_expected_genotype",
        "right_expected_genotype",
        "left_expected_presence_dosage",
        "right_expected_presence_dosage",
        "reference_biological_genotype",
        "reference_presence_dosage",
        "reference_trait_state",
    ]

    sites_out = selected.reindex(
        columns=site_columns
    ).rename(
        columns={
            "pos_numeric": "pos",
        }
    )

    sites_out.to_csv(
        OUT_SITES,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    selection_audit.to_csv(
        OUT_SELECTION_AUDIT,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    wide.to_csv(
        OUT_WIDE,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    long.to_csv(
        OUT_LONG,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    summary.to_csv(
        OUT_SUMMARY,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    spans.to_csv(
        OUT_SPANS,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    sources.to_csv(
        OUT_SOURCE_SUMMARY,
        sep="\t",
        index=False,
        na_rep="NA",
    )

    print(
        "\nFROZEN COMPARATIVE SITE COUNTS"
    )

    print(
        counts.to_string()
    )

    print(
        "\nFOCAL BIOLOGICAL-SPACE SUMMARY"
    )

    print(
        summary.to_string(
            index=False
        )
    )

    unresolved = long.loc[
        long[
            "allele_mapping_status"
        ]
        == "unresolved"
    ]

    print(
        "\nUNRESOLVED FOCAL ALLELE MAPPINGS:",
        len(unresolved),
    )

    if not unresolved.empty:
        print(
            unresolved[
                [
                    "module",
                    "chrom",
                    "pos",
                    "side",
                    "GT_raw",
                    "presence_associated_allele",
                    "focal_ref",
                    "focal_alt",
                ]
            ]
            .head(30)
            .to_string(
                index=False
            )
        )

    print(
        "\nOutputs:"
    )

    for path in [
        OUT_SITES,
        OUT_SELECTION_AUDIT,
        OUT_WIDE,
        OUT_LONG,
        OUT_SUMMARY,
        OUT_SPANS,
        OUT_SOURCE_SUMMARY,
    ]:
        print(
            f"  {path}"
        )


if __name__ == "__main__":
    main()
