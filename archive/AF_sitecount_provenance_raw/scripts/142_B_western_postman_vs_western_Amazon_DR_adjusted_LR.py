#!/usr/bin/env python3
"""Analysis 142: western postman B=1 versus western-Amazon D/R-positive B=0."""
import argparse
import hashlib
import importlib.util
import json
import re
import sys
import unicodedata
from pathlib import Path

import numpy as np
import pandas as pd
import scipy

NEW = "B_western_postman_vs_western_Amazon_DR_reanalysis"
OLD = "B_heurippa_besckei_no_DR_positive_reanalysis"
KEY = ["chrom", "pos", "ref", "alt"]
COHORTS = ("WGS_only", "WGS_plus_CS")
LOCKED_EXCLUSION = "ERR3653308"

SOURCES = {
    "Hausdorf2025": "https://onlinelibrary.wiley.com/doi/full/10.1111%2Fmec.17717",
    "Salazar2010": "https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000930",
    "Wallbank2016": "https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1002353",
}

# Exact literature-defined comparison. Missing country remains permissible for a
# taxon with a restricted rule, matching the validated analysis-141 handling.
RULES = {
    "melpomene cythera": (1, "Pacific_west_postman", ("ecuador",), "Hausdorf2025"),
    "melpomene vulcanus": (1, "Pacific_west_or_Panama_postman", ("colombia", "panama"), "Hausdorf2025"),
    "melpomene rosina": (1, "Central_America_postman", ("panama", "costa rica"), "Hausdorf2025"),
    "melpomene aglaope": (0, "western_Amazon_Dennis_Ray", ("peru",), "Wallbank2016"),
    "melpomene ecuadoriensis": (0, "western_Amazon_Dennis_Ray", ("ecuador",), "Wallbank2016"),
    "melpomene malleti": (0, "western_Amazon_Dennis_Ray", ("colombia", "ecuador"), "Wallbank2016"),
    "timareta florencia": (0, "western_Amazon_Dennis_Ray", ("colombia",), "Wallbank2016"),
    "timareta timareta f. contigua": (0, "western_Amazon_Ray_only", ("ecuador",), "Wallbank2016"),
}

REQUIRED_RULE_KEYS = set(RULES) | {"heurippa"}


def clean(value):
    return re.sub(r"\s+", " ", str(value)).strip()


def taxon_key(value):
    value = clean(value).lower().replace("'", "").replace('"', "")
    return re.sub(r"\b(?:heliconius|ssp|spp|subsp)\.?\s*", "", value).strip()


def country_key(value):
    value = unicodedata.normalize("NFKD", clean(value)).encode("ascii", "ignore").decode().lower()
    if value.startswith("panam"):
        return "panama"
    return "" if value in ("", "na", "nan", "none", "unknown") else value


def binary(value):
    value = clean(value).lower()
    if value in {"0", "0.0", "false"}:
        return 0
    if value in {"1", "1.0", "true"}:
        return 1
    return None


def comparative_rule(taxon, country):
    key, c = taxon_key(taxon), country_key(country)
    rule_key = key
    if key.startswith("heurippa"):
        rule_key = "heurippa"
        expected_b, region, countries, source = (
            1, "eastern_Andes_sensitivity_group", ("colombia",), "Salazar2010"
        )
    elif key == "melpomene melpomene" and c == "panama":
        rule_key = "melpomene melpomene"
        expected_b, region, countries, source = (
            1, "Central_America_postman", ("panama",), "Hausdorf2025"
        )
    elif key in RULES:
        expected_b, region, countries, source = RULES[key]
    else:
        return False, "not_in_predefined_comparison", "unassigned", "", None, rule_key

    if c and c not in countries:
        return False, "taxon_country_conflict", region, source, expected_b, rule_key
    decision = "predefined_taxon_country" if c else "predefined_taxon_country_missing"
    return True, decision, region, source, expected_b, rule_key


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8388608), b""):
            h.update(block)
    return h.hexdigest()


def save(frame, path):
    frame.to_csv(path, sep="\t", index=False, na_rep="NA")


def apply_locked_metadata_corrections(metadata):
    affected = []
    changed = []
    for accession, record in metadata.items():
        key = taxon_key(record.get("taxon", record.get("taxon_name", "")))
        if key in {"melpomene malleti", "melpomene ecuadoriensis"}:
            affected.append(accession)
            old = binary(record.get("B"))
            record["B"] = 0
            if old != 0:
                changed.append((accession, "B", old, 0))

    if len(affected) != 63:
        raise RuntimeError(f"Expected 63 malleti/ecuadoriensis metadata rows; found {len(affected)}")
    # The validated analysis-141 metadata loader already applies this
    # correction. Therefore zero changes here is expected; one change is also
    # accepted if the upstream loader is later run without that correction.
    if len(changed) not in (0, 1):
        raise RuntimeError(
            "Expected zero or one remaining malleti/ecuadoriensis B "
            f"value to change; found {len(changed)}"
        )
    if "ERR3652402" not in metadata:
        raise RuntimeError("ERR3652402 is absent from metadata")

    locked = []
    record = metadata["ERR3652402"]
    for trait, value in (("B", 0), ("D", 1), ("R", 1)):
        old = binary(record.get(trait))
        record[trait] = value
        locked.append(("ERR3652402", trait, old, value, old != value))

    print("LOCKED METADATA CORRECTIONS", flush=True)
    print(f"malleti/ecuadoriensis affected rows: {len(affected)}", flush=True)
    print(f"malleti/ecuadoriensis stored B values changed: {len(changed)}", flush=True)
    for row in locked:
        print("\t".join(map(str, row)), flush=True)


def prepare_plot(source, target):
    text = source.read_text()
    changes = {
        OLD: NEW,
        "broad_B_heurippa_besckei_no_DR_positive_BA_manhattan": "broad_B_western_postman_vs_western_Amazon_DR_BA_manhattan",
        "Band association after excluding D-positive and R-positive samples": "Band association: western postman versus western-Amazon Dennis/Ray",
        "Pooled genotypes at exact expected left-right sites": "Mixture-adjusted genotypes at exact expected left-right sites",
        "SNP selection and balanced accuracy use comparative samples only.": "Comparative SNP scoring excludes focal samples; focal calls are reciprocal-mixture posterior >=0.95 only.",
        "All other genotype combinations and QC failures are omitted.": "Other pairs and uncertain adjusted calls are omitted; no original-call fallback.",
    }
    for before, after in changes.items():
        if before not in text:
            raise RuntimeError(f"Plot source has changed; missing: {before}")
        text = text.replace(before, after)
    start = text.index("score_file <- file.path(")
    end = text.index("if (!file.exists(score_file))", start)
    text = text[:start] + 'score_file <- file.path(checkpoint_dir, "B_western_BA_and_adjusted_LR.tsv")\n\n' + text[end:]
    start = text.index("d$left_QC <- (")
    end = text.index("d$same_GT <- (", start)
    text = text[:start] + '''if (!("mixture_pair_posterior" %in% names(d))) stop("Missing adjusted-call posterior")
d$adjusted_confident <- is.finite(d$mixture_pair_posterior) & d$mixture_pair_posterior >= 0.95
d$left_QC <- d$adjusted_confident & !is.na(d$left_GT)
d$right_QC <- d$adjusted_confident & !is.na(d$right_GT)
d$both_QC <- d$left_QC & d$right_QC

''' + text[end:]
    text = text.replace('side = "TN_Left",\n        y = 1,', 'side = rep("TN_Left", nrow(difference_sites)),\n        y = rep(1, nrow(difference_sites)),')
    text = text.replace('side = "TN_Right",\n        y = 0,', 'side = rep("TN_Right", nrow(difference_sites)),\n        y = rep(0, nrow(difference_sites)),')
    compact = text.replace(",", "").replace("_", "")
    required_coordinates = {
        "zoom_start": "1060000",
        "zoom_end": "1240000",
        "optix_start": "1088602",
        "optix_end": "1097846",
        "Krista_start": "1112651",
        "Krista_end": "1114314",
    }
    missing = [label for label, value in required_coordinates.items() if value not in compact]
    if missing:
        raise RuntimeError("Generated plotter lost required coordinates: " + ", ".join(missing))
    target.write_text(text)


def validate_retained(retained, original):
    if retained.empty:
        raise RuntimeError("No comparative samples retained")
    taxa = retained["taxon_name"].astype(str).str.lower()
    if taxa.str.contains(r"cydno|pachinus", regex=True, na=False).any():
        raise RuntimeError("Cydno or pachinus survived the cohort definition")
    if retained["sample"].isin(original.FOCAL_SAMPLES).any():
        raise RuntimeError("A focal sample survived comparative cohort selection")
    if LOCKED_EXCLUSION in set(retained["sample"]):
        raise RuntimeError(f"Locked exclusion {LOCKED_EXCLUSION} survived")

    present = set(retained.loc[retained["cohort"].eq("WGS_plus_CS"), "rule_key"])
    missing = REQUIRED_RULE_KEYS - present
    if missing:
        raise RuntimeError("Required comparison groups absent from WGS_plus_CS: " + ", ".join(sorted(missing)))

    b0 = retained.loc[retained["B"].eq(0)]
    b1 = retained.loc[retained["B"].eq(1)]
    if b0.empty or b1.empty:
        raise RuntimeError("Both B classes must be represented")
    if (~(b0["D"].eq(1) | b0["R"].eq(1))).any():
        bad = b0.loc[~(b0["D"].eq(1) | b0["R"].eq(1))]
        raise RuntimeError("Retained B=0 sample lacks D=1/R=1:\n" + bad.to_string(index=False))
    if (b1["D"].eq(1) | b1["R"].eq(1)).any():
        bad = b1.loc[b1["D"].eq(1) | b1["R"].eq(1)]
        raise RuntimeError("Retained B=1 sample has D=1/R=1:\n" + bad.to_string(index=False))

    check = retained.loc[retained["sample"].eq("ERR3652402")]
    if check.empty or not (check["B"].eq(0) & check["D"].eq(1) & check["R"].eq(1)).all():
        raise RuntimeError("ERR3652402 is not retained as B=0,D=1,R=1")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--project", type=Path, default=Path("/scratch/martinlab/joe/GynandroCalling"))
    args = parser.parse_args()
    project = args.project.resolve()
    scripts = project / "scripts"
    out = project / "analysis" / NEW / "checkpoints"
    audit = out / "scorer_audit_original_fields"
    meta_path = project / "data/derived/region_panel_metadata.tsv"
    vcf = project / "external/krista_genotype_plots_original/Genotype_Plots/Optix_Hmr18-938195-1239421_biallelic.vcf"
    mix_path = project / "analysis/original_502/reciprocal_L_R_mixture/optix_reciprocal_L_R_genotype_pairs.tsv.gz"
    source_paths = [scripts / name for name in (
        "41_score_original_502_optix_callsets.py", "110b_base_scorer_B_no_DR_positive.py",
        "110c_scan_B_no_DR_positive.py", "110d_plot_B_no_DR_positive_BA_and_focal.R")]
    for path in source_paths + [meta_path, vcf, mix_path]:
        if not path.is_file():
            raise RuntimeError(f"Missing existing input: {path}")

    audit.mkdir(parents=True, exist_ok=True)
    original = load(source_paths[0], "west142_original")
    helper = load(source_paths[1], "west142_helper")
    scanner = load(source_paths[2], "west142_scanner")
    original.METADATA_FILE = meta_path
    helper.META = meta_path
    metadata = helper.load_metadata(original)
    apply_locked_metadata_corrections(metadata)

    source_meta = pd.read_csv(meta_path, sep="\t", dtype=str, keep_default_na=False)
    source_meta.columns = source_meta.columns.str.strip()
    required_meta = {"Accession", "Country", "Taxon Name"}
    if not required_meta.issubset(source_meta.columns):
        raise RuntimeError(f"Missing geographic metadata columns: {required_meta - set(source_meta.columns)}")
    source_meta["Accession"] = source_meta["Accession"].map(clean)
    if source_meta["Accession"].duplicated().any():
        raise RuntimeError("Duplicate accessions in metadata; resolve before geographic selection")
    locations = source_meta.set_index("Accession").to_dict("index")

    samples = None
    with vcf.open() as handle:
        for line in handle:
            if line.startswith("#CHROM\t"):
                samples = line.rstrip("\r\n").split("\t")[9:]
                break
    if samples is None or len(samples) != 502 or len(set(samples)) != 502:
        raise RuntimeError("Expected 502 unique VCF samples")
    if set(samples) - set(metadata) or set(samples) - set(locations):
        raise RuntimeError("VCF samples missing metadata")

    rows = []
    for cohort in COHORTS:
        for sample in samples:
            m, loc = metadata[sample], locations[sample]
            taxon, country = clean(m["taxon"]), clean(loc["Country"])
            include, decision, region, source, expected_b, rule_key = comparative_rule(taxon, country)
            b, d, r = binary(m.get("B")), binary(m.get("D")), binary(m.get("R"))
            seqtype = clean(m["data_type"]).upper()
            reasons = [] if include else [decision]
            if sample in original.FOCAL_SAMPLES:
                reasons.append("focal")
            if sample == LOCKED_EXCLUSION:
                reasons.append("locked_exclusion_contigua_not_DR_positive")
            if m.get("piano_key") or re.search(r"piano\s*keys?", taxon, re.I):
                reasons.append("piano_key")
            if b not in (0, 1):
                reasons.append("missing_B")
            elif include and b != expected_b:
                reasons.append(f"B_not_expected_{expected_b}")
            if include and expected_b == 0 and not (d == 1 or r == 1):
                reasons.append("B0_not_Dennis_or_Ray_positive")
            if rule_key == "timareta timareta f. contigua" and r != 1:
                reasons.append("contigua_not_R_positive")
            if include and expected_b == 1 and (d == 1 or r == 1):
                reasons.append("B1_Dennis_or_Ray_positive")
            if seqtype not in ({"WGS"} if cohort == "WGS_only" else {"WGS", "CS"}):
                reasons.append("sequencing_type")
            rows.append(dict(
                cohort=cohort, module="Band", sample=sample, taxon_name=taxon,
                country=country, latitude=loc.get("Latitude", ""), longitude=loc.get("Longitude", ""),
                data_type=seqtype, B=b, D=d, R=r, expected_B=expected_b, rule_key=rule_key,
                geographic_region=region, geographic_decision=decision,
                source=SOURCES.get(source, ""), keep=not reasons, keep_bool=not reasons,
                reason=";".join(reasons) if reasons else "included"))

    manifest = pd.DataFrame(rows)
    manifest_file = out / "B_western_cohort_manifest.tsv"
    save(manifest, manifest_file)
    retained = manifest.loc[manifest["keep"]].copy()
    validate_retained(retained, original)

    counts = retained.groupby(["cohort", "B"]).size().rename("samples").reset_index()
    taxon_counts = retained.groupby(
        ["cohort", "rule_key", "taxon_name", "country", "B", "D", "R"], dropna=False
    ).size().rename("samples").reset_index()
    trait_counts = retained.groupby(
        ["cohort", "B", "D", "R"], dropna=False
    ).size().rename("samples").reset_index()
    save(counts, out / "B_western_cohort_counts.tsv")
    save(taxon_counts, out / "B_western_retained_taxa.tsv")
    save(trait_counts, out / "B_western_retained_B_D_R_counts.tsv")
    save(retained, out / "B_western_retained_sample_manifest.tsv")
    print("WESTERN POSTMAN VERSUS WESTERN-AMAZON D/R COHORT\n", counts.to_string(index=False), flush=True)
    print(taxon_counts.to_string(index=False), flush=True)
    print("\nRetained B/D/R combinations\n", trait_counts.to_string(index=False), flush=True)
    print("\nLocked ERR3653308 manifest rows\n", manifest.loc[
        manifest["sample"].eq(LOCKED_EXCLUSION),
        ["cohort", "sample", "taxon_name", "B", "D", "R", "keep", "reason"],
    ].to_string(index=False), flush=True)
    print(
        "\nCydno and pachinus excluded; ERR3653308 locked out by group definition.\n"
        "Comparative sample/SNP selection excludes the focal samples.\n"
        "Focal calls are joined only after comparative scoring and use reciprocal-mixture posterior >=0.95 only.\n"
        "Original comparative callability and minimum class-count thresholds are unchanged.", flush=True)
    for cohort in COHORTS:
        n = retained.loc[retained["cohort"].eq(cohort), "B"].value_counts()
        if any(n.get(b, 0) < original.MINIMUM_CALLED_PER_PHENOTYPE for b in (0, 1)):
            raise RuntimeError(f"{cohort}: too few individuals for existing per-class eligibility rule")

    plotter = scripts / "142_plot_B_western_postman_vs_western_Amazon_DR_adjusted_LR.R"
    prepare_plot(source_paths[3], plotter)
    if args.prepare_only:
        print(f"\nPrepared: {plotter}\nManifest: {manifest_file}", flush=True)
        return

    provenance = dict(
        version="142.1", definition="western postman plus heurippa B=1 versus western-Amazon D/R-positive B=0",
        comparative_taxon_rules=RULES, sources=SOURCES,
        locked_metadata={"malleti_and_ecuadoriensis_B": 0, "ERR3652402": {"B": 0, "D": 1, "R": 1}},
        locked_exclusion={LOCKED_EXCLUSION: "contigua but not Dennis/Ray positive"},
        comparative_QC="Analysis-141 called_GT and DP2_GQ20 schemes; unchanged eligibility thresholds",
        focal_calls="Reciprocal mixture only; pair posterior >=0.95; no raw-call fallback",
        python=sys.version, numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__,
        input_SHA256={str(p): digest(p) for p in source_paths + [meta_path, vcf, mix_path, Path(__file__).resolve()]})
    signature = hashlib.sha256(json.dumps(provenance, sort_keys=True).encode()).hexdigest()
    (out / "B_western_provenance.json").write_text(json.dumps(provenance, indent=2))
    done = out / "scoring_complete.json"

    helper.load_original_module = lambda: original
    helper.load_metadata = lambda scorer: metadata
    helper.load_final_filtered_manifest = lambda: manifest.copy()
    scanner.load_helper = lambda: helper
    scanner.PROJECT, scanner.OUTDIR = project, audit
    for name, value in list(vars(scanner).items()):
        if name.startswith("OUT_") and isinstance(value, Path):
            setattr(scanner, name, audit / value.name)
    cached = json.loads(done.read_text()) if done.exists() else {}
    reuse = (cached.get("signature") == signature and scanner.OUT_RAW.is_file()
             and cached.get("raw_SHA256") == digest(scanner.OUT_RAW))
    if not reuse:
        scanner.main()
    else:
        print("Verified existing analysis-142 scores; rebuilding adjusted-call display input.", flush=True)

    raw = pd.read_csv(scanner.OUT_RAW, sep="\t", low_memory=False)
    if raw.groupby(["analysis_cohort", "scheme"]).ngroups != 4:
        raise RuntimeError("Expected two cohorts and two schemes")
    for (cohort, scheme), group in raw.groupby(["analysis_cohort", "scheme"]):
        expected = int(retained["cohort"].eq(cohort).sum())
        if len(group) != 57797 or group.duplicated(KEY).any() or set(group["comparative_cohort_n"]) != {expected}:
            raise RuntimeError(f"Scoring validation failed: {cohort}/{scheme}")

    # Focal data are first read here, after comparative scoring and its validation.
    mixture = pd.read_csv(mix_path, sep="\t", low_memory=False)
    columns = KEY + ["mixture_left_GT", "mixture_right_GT", "pair_posterior"]
    if not set(columns).issubset(mixture.columns) or mixture.duplicated(KEY).any():
        raise RuntimeError("Mixture table schema or site uniqueness failure")
    join = raw.merge(mixture[columns], on=KEY, how="left", validate="many_to_one", sort=False)
    if len(join) != len(raw) or not np.allclose(join["balanced_accuracy"], raw["balanced_accuracy"], equal_nan=True):
        raise RuntimeError("Focal join changed comparative BA or row count")
    confidence = pd.to_numeric(join["pair_posterior"], errors="coerce").ge(0.95)
    valid = {"0/0", "0/1", "1/1"}
    if confidence.any() and not (
        join.loc[confidence, "mixture_left_GT"].isin(valid).all()
        and join.loc[confidence, "mixture_right_GT"].isin(valid).all()
    ):
        raise RuntimeError("Invalid confident mixture genotype")
    if not confidence.any():
        raise RuntimeError("No confident adjusted calls joined; no original-call fallback")
    for side in ("Left", "Right"):
        join[f"TN_{side}_GT"] = join[f"mixture_{side.lower()}_GT"].where(confidence, "./.")
    join["mixture_pair_posterior"] = join["pair_posterior"]
    join = join.drop(columns=[
        "pooled_both_DP2_GQ20", "pooled_GT_different", "left_11_right_01",
        "dominant_focal_model_match", "mixture_left_GT", "mixture_right_GT", "pair_posterior"])
    join["focal_call_source"] = "reciprocal_mixture_posterior95_only"
    join["focal_pair_confident"] = confidence
    join["focal_exact_B_pair"] = confidence & join["TN_Left_GT"].eq("1/1") & join["TN_Right_GT"].eq("0/1")

    plot_input = out / "B_western_BA_and_adjusted_LR.tsv"
    save(join, plot_input)
    primary = join.loc[
        join["analysis_cohort"].eq("WGS_plus_CS") & join["scheme"].eq("DP2_GQ20")
    ].copy()
    primary["eligible"] = primary["score_eligible"].astype(str).str.lower().eq("true")
    regional = primary.loc[primary["pos"].between(1112651, 1114314)].sort_values(
        ["eligible", "balanced_accuracy", "pos"], ascending=[False, False, True])
    save(regional, out / "Krista_band_western_scores_and_adjusted_LR.tsv")
    show = ["pos", "selected_orientation", "balanced_accuracy", "misclassified_n",
            "phenotype_0_called", "phenotype_1_called", "eligible",
            "TN_Left_GT", "TN_Right_GT", "focal_pair_confident"]
    print("\nORIGINAL BAND REGION: RESTRICTED COHORT, DP2/GQ20", flush=True)
    print(regional[show].head(15).to_string(index=False), flush=True)

    summary_rows = []
    for region, frame in (
        ("full_VCF", primary),
        ("display_1060k_1240k", primary.loc[primary["pos"].between(1060000, 1240000)]),
        ("original_Band", regional),
    ):
        summary_rows.append(dict(
            region=region, SNPs=len(frame), eligible=int(frame["eligible"].sum()),
            eligible_BA90=int((frame["eligible"] & frame["balanced_accuracy"].ge(0.9)).sum()),
            eligible_BA1=int((frame["eligible"] & np.isclose(frame["balanced_accuracy"], 1)).sum()),
            adjusted_confident=int(frame["focal_pair_confident"].sum()),
            adjusted_exact_B=int(frame["focal_exact_B_pair"].sum()),
            adjusted_exact_and_BA90=int((frame["focal_exact_B_pair"] & frame["eligible"] & frame["balanced_accuracy"].ge(0.9)).sum())))
    summary = pd.DataFrame(summary_rows)
    save(summary, out / "B_western_summary.tsv")
    print("\n", summary.to_string(index=False), flush=True)
    done.write_text(json.dumps(dict(signature=signature, raw_SHA256=digest(scanner.OUT_RAW)), indent=2))
    save(pd.DataFrame([dict(file=str(p), SHA256=digest(p)) for p in (
        plot_input, plotter, manifest_file)]), out / "B_western_output_checksums.tsv")
    print(f"\nNew plot input: {plot_input}\nOriginal focal fields kept only in: {audit}\nPlotter: {plotter}", flush=True)


if __name__ == "__main__":
    main()
