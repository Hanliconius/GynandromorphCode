#!/usr/bin/env python3

from pathlib import Path
import re
import numpy as np
import pandas as pd


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

SITES = (
    PROJECT
    / "analysis"
    / "original_502"
    / "checkpoints"
    / "frozen_comparative_sites_reference_polarized.tsv"
)

VCF = (
    PROJECT
    / "data"
    / "raw"
    / "Optix_Hmr18-938195-1239421.vcf"
)

META = (
    PROJECT
    / "data"
    / "derived"
    / "region_panel_metadata.tsv"
)

EXCLUDED = (
    PROJECT
    / "analysis"
    / "checkpoints"
    / "optix_51sample_excluded_comparative_samples.tsv"
)

OUTDIR = (
    PROJECT
    / "analysis"
    / "original_502"
    / "checkpoints"
)

OUT_LONG = (
    OUTDIR
    / "optix_comparative_selected_sites_genotype_plot.long.tsv"
)

OUT_SITES = (
    OUTDIR
    / "optix_comparative_selected_sites_for_plotting.tsv"
)

MODULE_ORDER = [
    "Band",
    "Ray1",
    "Ray2",
    "Dennis1",
    "Dennis2",
]

DISPLAY = {
    "Band": "Band",
    "Ray1": "Ray 1",
    "Ray2": "Ray 2",
    "Dennis1": "Dennis 1",
    "Dennis2": "Dennis 2",
}

PHENO_COL = {
    "Band": "B",
    "Ray1": "R",
    "Ray2": "R",
    "Dennis1": "D",
    "Dennis2": "D",
}

REFERENCE_STATE = {
    "Band": "presence",
    "Ray1": "absence",
    "Ray2": "absence",
    "Dennis1": "absence",
    "Dennis2": "absence",
}


def fail(msg):
    raise RuntimeError(msg)


def read_tsv(path):
    if not path.exists():
        fail(f"Missing required file: {path}")
    return pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )


def emptyish(x):
    return x is None or str(x).strip() in {
        "", ".", "NA", "NaN", "nan", "None"
    }


def fnum(x):
    if emptyish(x):
        return np.nan
    try:
        return float(x)
    except Exception:
        return np.nan


def norm(x):
    if emptyish(x):
        return ""
    return str(x).strip()


def norm_allele(x):
    return norm(x).upper()


def truthy(x):
    return str(x).strip().lower() in {
        "1", "true", "t", "yes", "y"
    }


def parse_gt(gt):
    if emptyish(gt):
        return None

    gt = str(gt).strip()
    if gt in {".", "./.", ".|."}:
        return None

    parts = re.split(r"[/|]", gt)
    if len(parts) != 2:
        return None

    if any(p not in {"0", "1"} for p in parts):
        return None

    return sum(int(p) for p in parts)


def parse_ad(ad):
    if emptyish(ad):
        return (np.nan, np.nan)

    parts = str(ad).split(",")
    if len(parts) < 2:
        return (np.nan, np.nan)

    try:
        return (float(parts[0]), float(parts[1]))
    except Exception:
        return (np.nan, np.nan)


def biological_gt(module, raw_gt):
    alt_copies = parse_gt(raw_gt)
    if alt_copies is None:
        return None

    if REFERENCE_STATE[module] == "presence":
        presence_copies = 2 - alt_copies
    else:
        presence_copies = alt_copies

    return {
        0: "AA",
        1: "PA",
        2: "PP",
    }.get(presence_copies)


def genotype_code(gt):
    return {
        "AA": 0,
        "PA": 1,
        "PP": 2,
    }.get(gt, np.nan)


def find_first_present(df, candidates, required=True):
    for col in candidates:
        if col in df.columns:
            return col
    if required:
        fail(
            "Could not find any of columns: "
            + ", ".join(candidates)
        )
    return None


def parse_vcf_selected_sites(path, wanted_keys):
    samples = None
    rows = []

    wanted_by_pos = {}
    for chrom, pos, ref, alt in wanted_keys:
        wanted_by_pos.setdefault((chrom, pos), []).append((ref, alt))

    with path.open("rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM\t"):
                header = line.rstrip("\n").split("\t")
                samples = header[9:]
                continue

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue

            chrom = fields[0]
            try:
                pos = int(fields[1])
            except Exception:
                continue

            if (chrom, pos) not in wanted_by_pos:
                continue

            ref = norm_allele(fields[3])
            alt = norm_allele(fields[4])

            if (ref, alt) not in wanted_by_pos[(chrom, pos)]:
                continue

            fmt_names = fields[8].split(":")
            sample_calls = fields[9:]

            for sample, text in zip(samples, sample_calls):
                vals = text.split(":")
                fmt = {
                    k: (vals[i] if i < len(vals) else "")
                    for i, k in enumerate(fmt_names)
                }

                rows.append(
                    {
                        "sample": sample,
                        "chrom": chrom,
                        "pos": pos,
                        "ref": ref,
                        "alt": alt,
                        "GT_raw": fmt.get("GT", ""),
                        "AD_raw": fmt.get("AD", ""),
                        "DP": fnum(fmt.get("DP", "")),
                        "GQ": fnum(fmt.get("GQ", "")),
                    }
                )

    if samples is None:
        fail(f"Could not parse VCF header: {path}")

    return pd.DataFrame(rows), samples


def main():
    sites = read_tsv(SITES)
    meta = read_tsv(META)

    sites = sites.loc[
        sites["module"].isin(MODULE_ORDER)
    ].copy()

    sites["pos"] = pd.to_numeric(
        sites["pos"], errors="raise"
    ).astype(int)

    sites["ref"] = sites["ref"].map(norm_allele)
    sites["alt"] = sites["alt"].map(norm_allele)

    ba_col = find_first_present(
        sites,
        [
            "minimum_comparative_BA",
            "minimum_BA",
            "balanced_accuracy",
            "BA",
        ],
        required=False,
    )

    mcc_col = find_first_present(
        sites,
        [
            "minimum_comparative_MCC",
            "minimum_MCC",
            "MCC",
        ],
        required=False,
    )

    if ba_col is None:
        sites["plot_BA"] = np.nan
    else:
        sites["plot_BA"] = pd.to_numeric(
            sites[ba_col], errors="coerce"
        )

    if mcc_col is None:
        sites["plot_MCC"] = np.nan
    else:
        sites["plot_MCC"] = pd.to_numeric(
            sites[mcc_col], errors="coerce"
        )

    sites = sites.sort_values(
        ["module", "pos"]
    ).reset_index(drop=True)

    sites["site_order_within_module"] = (
        sites.groupby("module").cumcount() + 1
    )

    wanted_keys = set(
        zip(
            sites["chrom"],
            sites["pos"],
            sites["ref"],
            sites["alt"],
        )
    )

    calls, vcf_samples = parse_vcf_selected_sites(
        VCF, wanted_keys
    )

    if calls.empty:
        fail("No selected-site calls were recovered from Optix VCF.")

    # metadata columns
    acc_col = find_first_present(meta, ["Accession"], required=True)
    tax_col = find_first_present(meta, ["Taxon Name"], required=True)
    country_col = find_first_present(meta, ["Country"], required=False)

    for req in ["B", "R", "D"]:
        if req not in meta.columns:
            fail(f"Metadata missing required phenotype column: {req}")

    meta["sample"] = meta[acc_col].map(norm)

    meta["taxon_name"] = meta[tax_col].map(norm)

    meta["country_name"] = (
        meta[country_col].map(norm)
        if country_col is not None
        else ""
    )

    # Actual curated-panel membership is defined by presence
    # in the 51-sample Optix VCF.
    meta = meta.loc[
        meta["sample"].isin(vcf_samples)
    ].copy()

    # Focal samples are never comparative rows.
    meta = meta.loc[
        ~meta["sample"].str.startswith("TN_", na=False)
    ].copy()

    # Apply the previously defined comparative exclusions.
    if EXCLUDED.exists():

        excluded = read_tsv(EXCLUDED)

        if "sample" not in excluded.columns:
            fail(
                f"Expected 'sample' column in {EXCLUDED}"
            )

        excluded_samples = set(
            excluded["sample"].map(norm)
        )

        meta = meta.loc[
            ~meta["sample"].isin(excluded_samples)
        ].copy()

        print(
            "\nExplicit 51-panel exclusions:",
            len(excluded_samples),
        )

    else:
        print(
            f"\nWARNING: exclusion file absent: {EXCLUDED}"
        )

    print(
        "Metadata rows retained after VCF membership "
        "and exclusions:",
        len(meta),
    )

    calls["sample"] = calls["sample"].map(norm)
    calls["site_key"] = (
        calls["chrom"]
        + ":"
        + calls["pos"].astype(str)
        + ":"
        + calls["ref"]
        + ">"
        + calls["alt"]
    )

    sites["site_key"] = (
        sites["chrom"]
        + ":"
        + sites["pos"].astype(str)
        + ":"
        + sites["ref"]
        + ">"
        + sites["alt"]
    )

    long_rows = []

    for module in MODULE_ORDER:
        module_sites = sites.loc[
            sites["module"] == module
        ].copy()

        if module_sites.empty:
            continue

        pheno_col = PHENO_COL[module]

        module_meta = meta.copy()

        module_meta["phenotype_value"] = pd.to_numeric(
            module_meta[pheno_col], errors="coerce"
        )

        module_meta = module_meta.loc[
            module_meta["phenotype_value"].isin([0, 1])
        ].copy()

        if module_meta.empty:
            print(f"WARNING: no comparative samples for {module}")
            continue

        module_meta["phenotype_label"] = np.where(
            module_meta["phenotype_value"] == 1,
            "present",
            "absent",
        )

        # stable row order across all sites in module
        module_meta = module_meta.sort_values(
            [
                "phenotype_value",
                "taxon_name",
                "country_name",
                "sample",
            ],
            ascending=[True, True, True, True],
        ).reset_index(drop=True)

        module_meta["sample_order"] = np.arange(1, len(module_meta) + 1)

        for _, site in module_sites.iterrows():
            site_calls = calls.loc[
                (calls["chrom"] == site["chrom"])
                & (calls["pos"] == site["pos"])
                & (calls["ref"] == site["ref"])
                & (calls["alt"] == site["alt"])
            ].copy()

            merged = module_meta.merge(
                site_calls,
                on="sample",
                how="left",
            )

            for _, row in merged.iterrows():
                bgt = biological_gt(module, row.get("GT_raw", ""))

                ref_depth, alt_depth = parse_ad(row.get("AD_raw", ""))

                qc = bool(
                    bgt is not None
                    and np.isfinite(row.get("DP", np.nan))
                    and row.get("DP", np.nan) >= 2
                    and np.isfinite(row.get("GQ", np.nan))
                    and row.get("GQ", np.nan) >= 20
                )

                long_rows.append(
                    {
                        "module": module,
                        "module_display": DISPLAY[module],
                        "site_key": site["site_key"],
                        "site_order_within_module": int(site["site_order_within_module"]),
                        "chrom": site["chrom"],
                        "pos": int(site["pos"]),
                        "ref": site["ref"],
                        "alt": site["alt"],
                        "plot_BA": fnum(site["plot_BA"]),
                        "plot_MCC": fnum(site["plot_MCC"]),
                        "sample": row["sample"],
                        "sample_order": int(row["sample_order"]),
                        "taxon_name": row["taxon_name"],
                        "country_name": row["country_name"],
                        "phenotype_value": int(row["phenotype_value"]),
                        "phenotype_label": row["phenotype_label"],
                        "GT_raw": row.get("GT_raw", ""),
                        "AD_raw": row.get("AD_raw", ""),
                        "DP": row.get("DP", np.nan),
                        "GQ": row.get("GQ", np.nan),
                        "QC_pass_DP2_GQ20": qc,
                        "biological_genotype": bgt if qc else "missing",
                        "genotype_code": genotype_code(bgt) if qc else np.nan,
                        "ref_depth": ref_depth if qc else np.nan,
                        "alt_depth": alt_depth if qc else np.nan,
                    }
                )

    long = pd.DataFrame(long_rows)

    if long.empty:
        fail("No comparative genotype rows were built.")

    plot_sites = (
        sites.loc[
            sites["module"].isin(long["module"].unique())
        ,
            [
                "module",
                "site_key",
                "site_order_within_module",
                "chrom",
                "pos",
                "ref",
                "alt",
                "plot_BA",
                "plot_MCC",
            ]
        ]
        .sort_values(["module", "site_order_within_module"])
        .reset_index(drop=True)
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

    print("\nCOMPARATIVE GENOTYPE-PLOT SITE COUNTS")
    print(
        plot_sites.groupby("module")
        .size()
        .reindex(MODULE_ORDER, fill_value=0)
        .to_string()
    )

    print("\nCOMPARATIVE SAMPLE COUNTS BY MODULE")
    print(
        long[["module", "sample"]]
        .drop_duplicates()
        .groupby("module")
        .size()
        .reindex(MODULE_ORDER, fill_value=0)
        .to_string()
    )

    print("\nPHENOTYPE COUNTS BY MODULE")
    print(
        long[["module", "sample", "phenotype_label"]]
        .drop_duplicates()
        .groupby(["module", "phenotype_label"])
        .size()
        .to_string()
    )

    print("\nOutputs:")
    print(f"  {OUT_LONG}")
    print(f"  {OUT_SITES}")


if __name__ == "__main__":
    main()
