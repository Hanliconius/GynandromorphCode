#!/usr/bin/env python3
"""Reuse the existing Band scorer and plotter with an explicit new cohort."""
import argparse
import csv
import hashlib
import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy

EXCLUDE = ("amaryllis", "burchelli", "nanna", "xenoclea", "plesseni", "besckei")
OLD = "B_heurippa_besckei_no_DR_positive_reanalysis"
NEW = "B_six_taxon_exclusions_reanalysis"
COHORTS = ("WGS_only", "WGS_plus_CS")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def header(path):
    with path.open() as handle:
        for line in handle:
            if line.startswith("#CHROM\t"):
                samples = line.rstrip("\r\n").split("\t")[9:]
                if len(samples) != len(set(samples)):
                    raise RuntimeError("Duplicate VCF sample IDs")
                return samples
    raise RuntimeError(f"Missing VCF sample header: {path}")


def save(frame, path):
    frame.to_csv(path, sep="\t", index=False, na_rep="NA")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--project", type=Path,
                        default=Path("/scratch/martinlab/joe/GynandroCalling"))
    args = parser.parse_args()
    project = args.project.resolve()
    scripts = project / "scripts"
    out = project / "analysis" / NEW / "checkpoints"
    original_path = scripts / "41_score_original_502_optix_callsets.py"
    helper_path = scripts / "110b_base_scorer_B_no_DR_positive.py"
    scan_path = scripts / "110c_scan_B_no_DR_positive.py"
    plot_source = scripts / "110d_plot_B_no_DR_positive_BA_and_focal.R"
    plot_target = scripts / "139_plot_B_six_taxon_exclusions.R"
    meta_path = project / "data/derived/region_panel_metadata.tsv"
    vcf_dir = project / "external/krista_genotype_plots_original/Genotype_Plots"
    broad = vcf_dir / "Optix_Hmr18-938195-1239421_biallelic.vcf"
    band = vcf_dir / "optix_band-1112651-1114314.vcf"
    inputs = [original_path, helper_path, scan_path, plot_source, meta_path, broad, band]
    for path in inputs:
        if not path.is_file():
            raise RuntimeError(f"Missing existing input: {path}")
    out.mkdir(parents=True, exist_ok=True)

    original = load(original_path, "band139_original")
    helper = load(helper_path, "band139_helper")
    scanner = load(scan_path, "band139_scanner")
    original.METADATA_FILE = meta_path
    helper.META = meta_path
    metadata_raw = original.load_metadata()
    metadata = helper.load_metadata(original)
    if "ERR3652402" not in metadata:
        raise RuntimeError("Expected ERR3652402 metadata was not found")
    metadata["ERR3652402"]["B"] = 0
    for entry in metadata.values():
        entry["taxon"] = re.sub(r"\s+", " ", entry["taxon"]).strip()
        entry["data_type"] = entry["data_type"].strip().upper()

    corrections = []
    for sample, entry in metadata.items():
        before = metadata_raw[sample]["B"]
        if before != entry["B"]:
            corrections.append(dict(sample=sample, taxon=entry["taxon"],
                                    B_in_source=before, B_for_analysis=entry["B"]))
    save(pd.DataFrame(corrections, columns=["sample", "taxon", "B_in_source", "B_for_analysis"]),
         out / "B_phenotype_corrections.tsv")

    samples = header(broad)
    if len(samples) != 502 or set(samples) != set(header(band)):
        raise RuntimeError("The broad and regional VCFs must contain the same 502 samples")
    if set(samples) - set(metadata):
        raise RuntimeError("Some VCF samples have no metadata")

    manifest_rows = []
    for selection in ("baseline", "six_taxon_exclusions"):
        for cohort in COHORTS:
            for sample in samples:
                entry = metadata[sample]
                terms = [t for t in EXCLUDE if t in entry["taxon"].lower()]
                reasons = []
                if sample in original.FOCAL_SAMPLES:
                    reasons.append("focal")
                if re.search(r"piano\s*keys?", entry["taxon"], re.I):
                    reasons.append("piano_key")
                if entry["B"] not in (0, 1):
                    reasons.append("missing_B")
                permitted = {"WGS"} if cohort == "WGS_only" else {"WGS", "CS"}
                if entry["data_type"] not in permitted:
                    reasons.append("sequencing_type")
                if selection == "six_taxon_exclusions" and terms:
                    reasons.append("taxon:" + ",".join(terms))
                manifest_rows.append(dict(
                    selection=selection, cohort=cohort, module="Band", sample=sample,
                    taxon_name=entry["taxon"], data_type=entry["data_type"], B=entry["B"],
                    named_exclusion=",".join(terms), keep=not reasons, keep_bool=not reasons,
                    reason=";".join(reasons) if reasons else "included"))

    manifest = pd.DataFrame(manifest_rows)
    save(manifest, out / "B_comparative_cohort_manifest.tsv")
    chosen = manifest.loc[manifest.selection.eq("six_taxon_exclusions")].copy()
    kept = manifest.loc[manifest.keep]
    cohort_counts = kept.groupby(["selection", "cohort", "B"]).size().rename("samples").reset_index()
    save(cohort_counts, out / "B_comparative_cohort_counts.tsv")
    for selection in ("baseline", "six_taxon_exclusions"):
        for cohort in COHORTS:
            group = kept.loc[kept.selection.eq(selection) & kept.cohort.eq(cohort)]
            if set(group.B) != {0, 1}:
                raise RuntimeError(f"Both B classes required: {selection}, {cohort}")

    print("\nEXPLICIT NEW BAND COHORT", flush=True)
    print(cohort_counts.to_string(index=False), flush=True)
    print("\nSix excluded taxon terms: " + ", ".join(EXCLUDE), flush=True)
    print("No D/R-positive filter and no old keep-list or individual exclusions are inherited.", flush=True)
    print("Existing malleti/ecuadoriensis B=0 corrections are preserved and audited.", flush=True)
    print("Named-taxon exclusion is a sensitivity analysis, not a verified geographic assignment.", flush=True)

    # Change only the routing and cohort description in the existing R figure.
    plot_text = plot_source.read_text()
    changes = {
        OLD: NEW,
        "broad_B_heurippa_besckei_no_DR_positive_BA_manhattan":
            "broad_B_six_taxon_exclusions_BA_manhattan",
        "Band association after excluding D-positive and R-positive samples":
            "Band association after six named taxon exclusions",
        "SNP selection and balanced accuracy use comparative samples only.":
            "BA uses comparative samples only; this strip shows the expected focal pairs.",
    }
    for before, after in changes.items():
        if before not in plot_text:
            raise RuntimeError(f"Expected plotter text not found: {before}")
        plot_text = plot_text.replace(before, after)
    plot_target.write_text(plot_text)
    print(f"\nPrepared existing-style plotter: {plot_target}", flush=True)
    print(f"Cohort manifest: {out / 'B_comparative_cohort_manifest.tsv'}", flush=True)
    if args.prepare_only:
        return

    print("\nHashing inputs for reproducibility...", flush=True)
    hashes = {str(p): digest(p) for p in inputs + [Path(__file__).resolve()]}
    config = dict(
        version="139.1", excluded_taxa=EXCLUDE, cohorts=COHORTS,
        old_sample_filters_inherited=False, D_R_filter=False,
        original_eligibility=dict(
            overall_call_rate=original.MINIMUM_OVERALL_CALL_RATE,
            state_call_rate=original.MINIMUM_STATE_CALL_RATE,
            called_per_phenotype=original.MINIMUM_CALLED_PER_PHENOTYPE,
            predicted_class=original.MINIMUM_PREDICTED_CLASS),
        schemes=original.SCHEMES,
        orientation="Original best dominant REF/ALT orientation; REF expectation recorded",
        focal_calls="Original VCF TN_Left/TN_Right; no reciprocal-call substitution",
        python=sys.version, numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__,
        input_sha256=hashes)
    signature = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()
    (out / "B_reanalysis_provenance.json").write_text(json.dumps(config, indent=2))
    done_file = out / "scoring_complete.json"
    if done_file.exists() and not args.force:
        done = json.loads(done_file.read_text())
        if done.get("signature") == signature and all(
            Path(p).is_file() and digest(p) == sha for p, sha in done.get("outputs", {}).items()
        ) and done.get("outputs"):
            print("Verified matching completed scores: reusing them. Plotting can proceed.", flush=True)
            return

    # Compare all 153 original-region SNPs before and after the six exclusions.
    # Both comparisons use the same phenotype corrections and scoring rules.
    specification = dict(original.PANELS["band"], vcf=band)
    all_regional_rows = []
    old_make_cohort = original.make_cohort
    try:
        for selection in ("baseline", "six_taxon_exclusions"):
            for cohort in COHORTS:
                allowed = set(kept.loc[
                    kept.selection.eq(selection) & kept.cohort.eq(cohort), "sample"])
                original.make_cohort = helper.make_filtered_cohort_function(allowed)
                result = original.score_panel("band", specification, metadata)
                rows, scored_manifest, summary = helper.unpack_score_panel_result(result)
                recovered = {r["sample"] for r in scored_manifest if r["reason"] == "included"}
                if recovered != allowed:
                    raise RuntimeError("Regional scoring cohort does not match the new manifest")
                for row in rows:
                    all_regional_rows.append(dict(selection=selection, analysis_cohort=cohort, **row))
    finally:
        original.make_cohort = old_make_cohort

    regional = pd.DataFrame(all_regional_rows)
    keys = ["analysis_cohort", "scheme", "chrom", "pos", "ref", "alt"]
    columns = keys + [
        "selected_orientation", "balanced_accuracy", "sensitivity", "specificity",
        "misclassified_n", "call_rate", "phenotype_0_call_rate", "phenotype_1_call_rate",
        "score_eligible", "matches_reference_expectation"
    ]
    before = regional.loc[regional.selection.eq("baseline"), columns]
    after = regional.loc[regional.selection.eq("six_taxon_exclusions"), columns]
    comparison = before.merge(
        after, on=keys, suffixes=("_before", "_after"), validate="one_to_one")
    comparison["BA_change"] = (
        comparison.balanced_accuracy_after - comparison.balanced_accuracy_before)
    comparison = comparison.sort_values(
        ["analysis_cohort", "scheme", "score_eligible_after", "balanced_accuracy_after", "pos"],
        ascending=[True, True, False, False, True])
    regional_file = out / "Krista_band_all_153_SNP_scores.tsv"
    comparison_file = out / "Krista_band_153_SNPs_before_vs_after.tsv"
    save(regional, regional_file)
    save(comparison, comparison_file)
    focus = comparison.loc[
        comparison.analysis_cohort.eq("WGS_plus_CS") & comparison.scheme.eq("DP2_GQ20")]
    if len(focus) != 153:
        raise RuntimeError(f"Expected 153 regional SNPs; recovered {len(focus)}")
    print("\nORIGINAL BAND REGION: TOP SITES AFTER EXCLUSIONS", flush=True)
    print(focus[[
        "pos", "selected_orientation_after", "balanced_accuracy_before",
        "balanced_accuracy_after", "misclassified_n_after", "score_eligible_after"
    ]].head(12).to_string(index=False), flush=True)

    # Route the unchanged broad scanner to the explicitly constructed cohort.
    # It does not read the old filtered manifest or inherit its restrictions.
    helper.load_original_module = lambda: original
    helper.load_metadata = lambda scorer: metadata
    helper.load_final_filtered_manifest = lambda: chosen.copy()
    scanner.load_helper = lambda: helper
    scanner.PROJECT = project
    scanner.OUTDIR = out
    for name, value in list(vars(scanner).items()):
        if name.startswith("OUT_") and isinstance(value, Path):
            setattr(scanner, name, out / value.name)
    print("\nRESCORING EXISTING BROAD VCF WITH THE NEW MANIFEST", flush=True)
    scanner.main()

    # Validate the full interval, both schemes, and both intended cohorts.
    raw = pd.read_csv(scanner.OUT_RAW, sep="\t", low_memory=False)
    for (cohort, scheme), group in raw.groupby(["analysis_cohort", "scheme"]):
        expected_n = int((chosen.cohort.eq(cohort) & chosen.keep_bool).sum())
        if set(group.comparative_cohort_n) != {expected_n}:
            raise RuntimeError(f"Broad cohort mismatch: {cohort}/{scheme}")
        if len(group) != 57797 or group.duplicated(["chrom", "pos", "ref", "alt"]).any():
            raise RuntimeError(f"Expected 57797 unique broad SNPs: {cohort}/{scheme}")
    if raw.groupby(["analysis_cohort", "scheme"]).ngroups != 4:
        raise RuntimeError("Expected two cohorts and two QC schemes")

    # The original-region and broad VCFs must reproduce the same scores.
    regional_after = regional.loc[regional.selection.eq("six_taxon_exclusions")]
    check = regional_after.merge(
        raw, on=keys, suffixes=("_region", "_broad"), validate="one_to_one")
    numeric_fields = ["balanced_accuracy", "tp", "tn", "fp", "fn", "comparative_cohort_n"]
    concordant = np.ones(len(check), dtype=bool)
    for field in numeric_fields:
        concordant &= np.isclose(
            check[field + "_region"], check[field + "_broad"],
            equal_nan=True, atol=1e-10, rtol=0)
    check["regional_broad_scores_agree"] = concordant
    crosscheck = out / "Krista_band_regional_vs_broad_score_check.tsv"
    save(check, crosscheck)
    print(f"\nRegional/broad matching score rows: {int(concordant.sum())}/{len(check)}", flush=True)
    if len(check) != 153 * 4 or not concordant.all():
        raise RuntimeError("Regional and broad scores differ; audit saved. Plotting stopped.")

    products = [regional_file, comparison_file, crosscheck]
    products += [
        v for k, v in vars(scanner).items()
        if k.startswith("OUT_") and isinstance(v, Path)
    ]
    done_file.write_text(json.dumps(
        dict(signature=signature, outputs={str(p): digest(p) for p in products}),
        indent=2))
    print(f"\nCompleted, validated scores: {out}", flush=True)
    print(f"Run the prepared plotter: {plot_target}", flush=True)


if __name__ == "__main__":
    main()
