#!/usr/bin/env python3

from pathlib import Path
import csv


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

CHECKPOINT_DIR = (
    PROJECT
    / "analysis"
    / "drop_problematic_reanalysis"
    / "checkpoints"
)

COMPARATIVE_FILE = (
    CHECKPOINT_DIR
    / "perfect_D_R_all_comparative_genotypes.long.tsv"
)

BROAD_VCF = (
    PROJECT
    / "external"
    / "krista_genotype_plots_original"
    / "Genotype_Plots"
    / "Optix_Hmr18-938195-1239421_biallelic.vcf"
)

OUTPUT_FILE = (
    CHECKPOINT_DIR
    / "perfect_D_R_focal_pooled_genotypes.long.tsv"
)

FOCAL_SAMPLES = {
    "TN_Left": {
        "label": "TN_Left | pooled left",
        "taxon": "Gynandromorph pooled left",
        "plot_order": "0",
    },
    "TN_Right": {
        "label": "TN_Right | pooled right",
        "taxon": "Gynandromorph pooled right",
        "plot_order": "1",
    },
}


def clean(value):
    if value is None:
        return ""
    return str(value).strip()


def canonical_gt(value):
    value = clean(value).replace("|", "/")

    if value in {
        "",
        ".",
        "./.",
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


def biological_gt(gt):
    return {
        "0/0": "A/A",
        "0/1": "P/A",
        "1/1": "P/P",
    }.get(
        gt,
        "No call",
    )


def numeric(value):
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


def target_key(row):
    return (
        clean(row.get("chrom")),
        clean(row.get("pos")),
        clean(row.get("ref")),
        clean(row.get("alt")),
    )


if not COMPARATIVE_FILE.exists():
    raise RuntimeError(
        f"Missing comparative table: {COMPARATIVE_FILE}"
    )

if not BROAD_VCF.exists():
    raise RuntimeError(
        f"Missing broad Optix VCF: {BROAD_VCF}"
    )


with COMPARATIVE_FILE.open(newline="") as handle:
    reader = csv.DictReader(
        handle,
        delimiter="\t",
    )

    output_fields = list(
        reader.fieldnames or []
    )

    comparative_rows = list(reader)


if not comparative_rows:
    raise RuntimeError(
        "Comparative genotype table is empty."
    )


target_templates = {}
physical_targets = {}

for row in comparative_rows:
    module = clean(row["module"])
    key = target_key(row)
    module_key = (
        module,
        *key,
    )

    if module_key not in target_templates:
        target_templates[module_key] = row

    physical_targets.setdefault(
        key,
        [],
    ).append(module_key)


for key in physical_targets:
    physical_targets[key] = sorted(
        set(
            physical_targets[key]
        )
    )


sample_column_indices = {}
output_rows = []
observed = set()
vcf_records_scanned = 0
target_records_found = 0


with BROAD_VCF.open(
    "r",
    errors="replace",
) as handle:

    for line in handle:

        if line.startswith("##"):
            continue

        if line.startswith("#CHROM"):
            header = line.rstrip("\n").split("\t")
            vcf_samples = header[9:]

            for sample in FOCAL_SAMPLES:
                if sample not in vcf_samples:
                    raise RuntimeError(
                        f"Focal sample absent from VCF: {sample}"
                    )

                sample_column_indices[sample] = (
                    9 + vcf_samples.index(sample)
                )

            continue

        if line.startswith("#"):
            continue

        vcf_records_scanned += 1

        fields = line.rstrip("\n").split("\t")

        if len(fields) < 10:
            continue

        key = (
            clean(fields[0]),
            clean(fields[1]),
            clean(fields[3]),
            clean(fields[4]),
        )

        if key not in physical_targets:
            continue

        target_records_found += 1

        format_keys = fields[8].split(":")

        for sample, sample_info in FOCAL_SAMPLES.items():
            column_index = sample_column_indices[sample]

            if column_index >= len(fields):
                raise RuntimeError(
                    f"Missing VCF sample field for {sample} "
                    f"at {key[0]}:{key[1]}"
                )

            sample_values = fields[column_index].split(":")

            format_values = {
                name: (
                    sample_values[index]
                    if index < len(sample_values)
                    else ""
                )
                for index, name in enumerate(format_keys)
            }

            gt = canonical_gt(
                format_values.get("GT")
            )

            dp = numeric(
                format_values.get("DP")
            )

            gq = numeric(
                format_values.get("GQ")
            )

            ad = clean(
                format_values.get("AD")
            )

            qc_pass = (
                gt is not None
                and dp is not None
                and gq is not None
                and dp >= 2
                and gq >= 20
            )

            display_biological_gt = (
                biological_gt(gt)
                if qc_pass
                else "No call"
            )

            for module_key in physical_targets[key]:
                template = target_templates[module_key]
                module = clean(template["module"])

                row = {
                    field: ""
                    for field in output_fields
                }

                row.update(
                    {
                        "module": module,
                        "trait": clean(template["trait"]),
                        "chrom": key[0],
                        "pos": key[1],
                        "ref": key[2],
                        "alt": key[3],
                        "site_key": clean(
                            template["site_key"]
                        ),
                        "perfect_both_schemes": clean(
                            template[
                                "perfect_both_schemes"
                            ]
                        ),
                        "sample": sample,
                        "sample_label": sample_info[
                            "label"
                        ],
                        "plot_order": sample_info[
                            "plot_order"
                        ],
                        "phenotype": "focal",
                        "taxon": sample_info["taxon"],
                        "country": "",
                        "data_type": "pooled WGS",
                        "GT": (
                            gt
                            if gt is not None
                            else "./."
                        ),
                        "biological_GT": (
                            display_biological_gt
                        ),
                        "AD": ad,
                        "DP": (
                            ""
                            if dp is None
                            else str(dp)
                        ),
                        "GQ": (
                            ""
                            if gq is None
                            else str(gq)
                        ),
                        "DP2_GQ20_pass": str(
                            qc_pass
                        ),
                        "predicted_presence": (
                            ""
                            if not qc_pass
                            else str(
                                int(
                                    gt in {
                                        "0/1",
                                        "1/1",
                                    }
                                )
                            )
                        ),
                        "callable_discordance": "",
                    }
                )

                output_rows.append(row)

                observed.add(
                    (
                        module,
                        clean(template["site_key"]),
                        sample,
                    )
                )


expected = set()

for module_key, template in target_templates.items():
    module = clean(template["module"])
    site = clean(template["site_key"])

    for sample in FOCAL_SAMPLES:
        expected.add(
            (
                module,
                site,
                sample,
            )
        )


missing = sorted(
    expected - observed
)

if missing:
    preview = "\n".join(
        "\t".join(item)
        for item in missing[:30]
    )

    raise RuntimeError(
        f"Missing {len(missing)} focal genotype cells "
        f"after scanning the broad VCF.\n"
        f"First missing entries:\n{preview}"
    )


output_rows.sort(
    key=lambda row: (
        clean(row["module"]),
        int(clean(row["pos"])),
        int(clean(row["plot_order"])),
    )
)


with OUTPUT_FILE.open(
    "w",
    newline="",
) as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=output_fields,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )

    writer.writeheader()
    writer.writerows(output_rows)


print("FOCAL POOLED GENOTYPES AT PERFECT D/R SNPS")
print()
print(f"VCF records scanned: {vcf_records_scanned}")
print(f"Target physical records found: {target_records_found}")
print()
print(
    "module\tperfect_SNPs\tfocal_samples\t"
    "genotype_cells\tQC_pass_cells\tQC_fail_cells"
)

for module in sorted(
    {
        row["module"]
        for row in output_rows
    }
):
    module_rows = [
        row
        for row in output_rows
        if row["module"] == module
    ]

    site_n = len(
        {
            row["site_key"]
            for row in module_rows
        }
    )

    sample_n = len(
        {
            row["sample"]
            for row in module_rows
        }
    )

    qc_pass_n = sum(
        clean(row["DP2_GQ20_pass"]).lower()
        == "true"
        for row in module_rows
    )

    print(
        f"{module}\t"
        f"{site_n}\t"
        f"{sample_n}\t"
        f"{len(module_rows)}\t"
        f"{qc_pass_n}\t"
        f"{len(module_rows) - qc_pass_n}"
    )

print()
print(f"Output: {OUTPUT_FILE}")
