from pathlib import Path
import csv
import json
import numpy as np
import pysam
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path("analysis/figures/pairwise_AF_cli")
OUT = ROOT / "control_histogram_comparison"
OUT.mkdir(parents=True, exist_ok=True)

audit = json.loads(
    (ROOT / "TN_merged_control_sites_hex40_batlowW.AF.audit.json").read_text()
)
vcf = Path(audit["input"]["path"])
bamlist = Path(
    "analysis/genomewide_mixing_input/metadata/bam_sets_24/controls.bamlist"
)

# Identify genuine controls from their BAM sample names, not filenames.
controls = set()
for line in bamlist.read_text().splitlines():
    if not line.strip():
        continue
    with pysam.AlignmentFile(line.strip(), "rb") as bam:
        names = {
            row.get("SM")
            for row in bam.header.to_dict().get("RG", [])
        }
    if len(names) != 1 or None in names:
        raise RuntimeError(f"Ambiguous sample name: {line}")
    controls.update(names)

reference = "ERR1307111"
if reference not in controls:
    raise RuntimeError("ERR1307111 is absent from the control manifest.")

samples = [reference] + sorted(controls - {reference})[:5]
if len(samples) != 6:
    raise RuntimeError("Could not identify six controls.")

contig = "NC_140389.1"
edges = np.linspace(0, 1, 101)
individual = np.zeros((6, 100), dtype=np.int64)
common = np.zeros_like(individual)
depth_sum = np.zeros(6, dtype=np.int64)
read = 0

with pysam.VariantFile(str(vcf)) as handle:
    missing = set(samples) - set(handle.header.samples)
    if missing:
        raise RuntimeError(f"Controls missing from VCF: {sorted(missing)}")
    if contig not in handle.header.contigs:
        raise RuntimeError(f"Chromosome-1 contig absent: {contig}")
    handle.subset_samples(samples)

    try:
        records = handle.fetch(contig)
    except (ValueError, OSError):
        print("No usable index: streaming until chromosome 1 ends.", flush=True)
        records = iter(handle)

    entered = False
    for record in records:
        if record.contig != contig:
            if entered:
                break
            continue
        entered = True
        read += 1

        if (
            record.ref not in {"A", "C", "G", "T"}
            or not record.alts
            or len(record.alts) != 1
            or record.alts[0] not in {"A", "C", "G", "T"}
        ):
            continue

        indices = []
        for i, sample in enumerate(samples):
            ad = record.samples[sample].get("AD")
            if (
                ad is None or len(ad) != 2
                or any(value is None or value < 0 for value in ad)
                or sum(ad) < 10
            ):
                indices.append(None)
                continue

            depth = sum(ad)
            fraction = ad[1] / depth
            index = min(int(fraction * 100), 99)
            individual[i, index] += 1
            depth_sum[i] += depth
            indices.append(index)

        if all(index is not None for index in indices):
            for i, index in enumerate(indices):
                common[i, index] += 1

        if read % 500000 == 0:
            print(f"READ\t{read:,}", flush=True)

if not read or np.any(individual.sum(axis=1) == 0):
    raise RuntimeError("No usable chromosome-1 data for one or more controls.")
if common.sum() == 0:
    raise RuntimeError("No jointly callable sites.")

np.savez_compressed(
    OUT / "control_AF_counts.npz",
    individual=individual, common=common,
    samples=np.array(samples), edges=edges,
)

with (OUT / "control_AF_summary.tsv").open("w") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow([
        "sample", "callable_sites", "common_sites",
        "mean_AD_depth_at_callable_sites",
    ])
    for i, sample in enumerate(samples):
        n = int(individual[i].sum())
        writer.writerow([
            sample, n, int(common[i].sum()), depth_sum[i] / n,
        ])

plt.rcParams.update({"font.size": 11, "pdf.fonttype": 42})
maximum = max(10, float(individual.max()) * 1.2)

pdf_path = OUT / "six_controls_chr1_raw_AF_histograms.pdf"
with PdfPages(pdf_path) as pdf:
    for label, values, description in [
        (
            "individual_callable",
            individual,
            "Each control's callable sites",
        ),
        (
            "common_callable",
            common,
            "Exactly the same sites: callable in all six controls",
        ),
    ]:
        fig, axes = plt.subplots(2, 3, figsize=(16, 9), sharex=True, sharey=True)
        fig.subplots_adjust(
            left=.07, right=.97, bottom=.12,
            top=.83, hspace=.32, wspace=.15,
        )
        for i, axis in enumerate(axes.flat):
            axis.stairs(
                values[i], edges,
                color="#009E73" if i == 0 else "#0072B2",
                linewidth=1.5, alpha=.85,
            )
            axis.set_yscale("symlog", linthresh=10)
            axis.set(
                xlim=(-.01, 1.01), ylim=(0, maximum),
                title=f"{samples[i]} | n={values[i].sum():,}",
                xlabel="ALT-read fraction", ylabel="Sites",
            )
            axis.spines[["top", "right"]].set_visible(False)

        fig.suptitle(description, fontsize=20, weight="bold", y=.95)
        fig.text(
            .07, .89,
            "Chromosome 1 | 100 bins | raw counts | identical axes across samples and pages",
            fontsize=13,
        )
        fig.text(
            .07, .045,
            "Biallelic SNPs; REF+ALT depth >=10; no GT or VCF FILTER restriction. "
            "Y-scale: linear to 10 sites, logarithmic above.\n"
            "ERR1307111 plus five alphabetically selected controls; "
            "selection does not use AF shape. No smoothing.",
            fontsize=10,
        )
        pdf.savefig(fig)
        fig.savefig(OUT / f"six_controls_{label}.png", dpi=160)
        plt.close(fig)

print(f"PDF\t{pdf_path.resolve()}")
print((OUT / "control_AF_summary.tsv").read_text())
