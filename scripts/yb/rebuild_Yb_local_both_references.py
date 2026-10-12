#!/usr/bin/env python3
"""Submit or execute a matched, local-only two-reference Yb rebuild."""
import argparse
import csv
import uuid
import importlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone

ROOT = Path("/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling")
SAMPLES = {"TN_Left_T2": "TN2L_S1", "TN_Right_T2": "TN2R_S2",
           "TN_Left_T3": "TN3L_S3", "TN_Right_T3": "TN3R_S4"}
HERE = Path(__file__).resolve().parent

def require(path):
    p = Path(path)
    if not p.is_absolute() or not p.is_file() or p.stat().st_size == 0:
        raise RuntimeError(f"Required nonempty absolute input missing: {p}")
    return str(p.resolve())

def tool(name):
    found = shutil.which(name)
    if found:
        return str(Path(found).absolute())
    envroot = Path.home() / ".local/share/mamba/envs"
    preferred = ("gynandrocalling", "colias_w_popgen")
    candidates = [envroot / env / "bin" / name for env in preferred]
    candidates += sorted(envroot.glob(f"*/bin/{name}"))
    for p in candidates:
        if p.is_file() and os.access(p, os.X_OK):
            return str(p.absolute())
    raise RuntimeError(f"Required tool unavailable: {name}; see tool_setup.log")

def run(command):
    subprocess.run([str(x) for x in command], check=True)

def check_dependencies():
    """Collect every missing/broken dependency before reporting failure."""
    errors, tools, checks = [], {}, []
    directory = ROOT / "analysis/yb_local_both_references"
    directory.mkdir(parents=True, exist_ok=True)
    logpath = directory / "dependency_check.log"
    with logpath.open("w") as log:
        for package in ("matplotlib", "numpy", "pandas", "pyBigWig"):
            try:
                imported = importlib.import_module(package)
                checks.append((package, "OK", str(getattr(imported, "__version__", "importable"))))
            except Exception as error:
                errors.append(package)
                checks.append((package, "FAIL", str(error)))
        for name in ("bwa", "samtools", "bamCoverage", "bash", "sbatch", "scancel"):
            try:
                tools[name] = tool(name)
                command = [tools[name]] if name == "bwa" else [tools[name], "--version"]
                result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                        text=True, timeout=20)
                log.write(shlex.join(command) + "\n" + result.stdout + "\n")
                expected = (0, 1) if name == "bwa" else (0,)
                if result.returncode not in expected:
                    raise RuntimeError(f"exit {result.returncode}; see {logpath}")
                if name == "bwa" and "bwa" not in result.stdout.lower():
                    raise RuntimeError(f"BWA help not detected; see {logpath}")
                checks.append((name, "OK", tools[name]))
            except Exception as error:
                errors.append(name)
                checks.append((name, "FAIL", str(error)))
        try:
            picard = shutil.which("picard")
            if picard:
                tools["picard"] = str(Path(picard).absolute())
                command = [tools["picard"], "MarkDuplicates", "--help"]
                detail = tools["picard"]
            else:
                candidates = [os.environ.get("YB_PICARD_JAR", ""), os.environ.get("PICARD_JAR", ""),
                              str(ROOT / "picard.jar"), str(HERE.parent.parent / "picard.jar")]
                candidates += [str(p) for p in (Path.home() / ".local/share/mamba/envs")
                               .glob("*/share/picard*/picard.jar")]
                jar = next((p for p in candidates if p and Path(p).is_file()), None)
                if not jar:
                    raise RuntimeError("No picard command/jar; set YB_PICARD_JAR to an existing jar")
                tools["java"], tools["picard_jar"] = tool("java"), require(jar)
                command = [tools["java"], "-jar", tools["picard_jar"], "MarkDuplicates", "--help"]
                detail = tools["picard_jar"]
            result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                    text=True, timeout=30)
            log.write(shlex.join(command) + "\n" + result.stdout + "\n")
            if result.returncode != 0 or "MarkDuplicates" not in result.stdout:
                raise RuntimeError(f"MarkDuplicates help failed (exit {result.returncode}); see {logpath}")
            checks.append(("Picard/Java", "OK", detail))
        except Exception as error:
            errors.append("Picard/Java")
            checks.append(("Picard/Java", "FAIL", str(error)))
    print(f"PYTHON\t{sys.executable}", flush=True)
    for name, status, detail in checks:
        print(f"{status}\t{name}\t{detail}", flush=True)
    print(f"DETAIL_LOG\t{logpath}", flush=True)
    if errors:
        raise RuntimeError("Missing or broken dependencies: " + ", ".join(errors))
    print("ALL_DEPENDENCIES_OK", flush=True)
    return tools

def submit():
    source = ROOT / "analysis/genomewide_mixing_input/metadata/resolved_download_manifest.tsv"
    inputs = {}
    with source.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("role") == "focal":
                name = row["library_id"]
                if name in inputs or name not in SAMPLES:
                    raise RuntimeError(f"Unexpected/duplicate focal library: {name}")
                inputs[name] = {"sample": SAMPLES[name], "library_id": name,
                                "r1": require(row["fastq1"]), "r2": require(row["fastq2"]),
                                "side": row["side"], "leg": row["leg"]}
    if set(inputs) != set(SAMPLES):
        raise RuntimeError("Manifest must contain exactly the four focal libraries")
    hmros = require((ROOT / "analysis/genomewide_mixing_input/metadata/reference_fasta.path")
                    .read_text().strip().splitlines()[0])
    hm25 = require(ROOT / "analysis/yb_exact_coordinate_transfer/reference/hmel2_5.chromosomes.fa")
    # Confirm native Yb contigs using existing FASTA indexes before submitting.
    for fasta, contig in ((hmros, "NC_140403.1"), (hm25, "Hmel215003o")):
        fai = require(fasta + ".fai")
        if not any(line.split("\t")[0] == contig for line in Path(fai).read_text().splitlines()):
            raise RuntimeError(f"{contig} absent from {fai}")
    tools = check_dependencies()
    runid = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:8]
    base = ROOT / "analysis/yb_local_both_references" / runid
    base.mkdir(parents=True)
    (base / "logs").mkdir()
    refs = [{"name": "Hmel2.5", "source": hm25, "contig": "Hmel215003o"},
            {"name": "HmRos", "source": hmros, "contig": "NC_140403.1"}]
    tasks = [{"reference": ref, **inputs[name]} for ref in refs for name in SAMPLES]
    config = {"run_dir": str(base), "references": refs, "tasks": tasks, "tools": tools,
              "python": sys.executable, "code_dir": str(HERE),
              "manifest_source": str(source),
              "recipe": {"alignment": "BWA-MEM -M; Picard MarkDuplicates; retain marked duplicates",
                         "coverage": "bamCoverage 50-bp bins; no normalization or MAPQ/duplicate filter"},
              "upstream_commit": "65f4313299097fafed0df9e44654ddb98b6dcd61"}
    configpath = base / "run.json"
    configpath.write_text(json.dumps(config, indent=2) + "\n")
    # Freeze executable code to keep queued jobs independent of later git pulls.
    shutil.copy2(__file__, base / "rebuild_Yb_local_both_references.py")
    shutil.copy2(HERE / "plot_Yb_local_both_references.py", base / "plot_Yb_local_both_references.py")
    config["code_dir"] = str(base)
    configpath.write_text(json.dumps(config, indent=2) + "\n")
    worker = base / "worker.sh"
    worker.write_text("#!/usr/bin/env bash\nset -euo pipefail\n" +
                     "exec " + shlex.join([sys.executable, str(base / Path(__file__).name),
                     "--config", str(configpath)]) + ' --stage "$1" --task "$SLURM_ARRAY_TASK_ID"\n')
    log = base / "logs"
    jobs = []
    try:
        for stage, array, cpus, mem, limit, dependency in (
            ("prepare", "0-1%2", 8, "32G", "06:00:00", None),
            ("align", "0-7%4", 16, "64G", "2-00:00:00", "prepare"),
            ("plot", "0-1%2", 1, "8G", "00:30:00", "align")):
            command = ["sbatch", "--parsable", "--partition=cpu", f"--job-name=YbLocal_{stage}",
                       f"--array={array}", f"--cpus-per-task={cpus}", f"--mem={mem}",
                       f"--time={limit}", "--export=ALL", f"--chdir={base}",
                       f"--output={log}/{stage}.%A_%a.out", f"--error={log}/{stage}.%A_%a.err"]
            if dependency:
                command.append("--dependency=afterok:" + jobs[-1][1])
            command += [str(worker), stage]
            job = subprocess.check_output(command, text=True).strip().split(";")[0]
            jobs.append((stage, job))
            (base / "jobs.json").write_text(json.dumps(dict(jobs), indent=2) + "\n")
            print(f"{stage.upper()}_ARRAY\t{job}", flush=True)
    except Exception:
        for _, job in jobs:
            subprocess.run(["scancel", job], check=False)
        raise
    (ROOT / "analysis/yb_local_both_references/latest_run.path").write_text(str(base) + "\n")
    print(f"RUN_DIR\t{base}")
    print("ALIGNMENTS\t8 tasks; at most 4 concurrently")
    print("INPUTS\texisting local FASTQs and FASTAs only")

def execute(configpath, stage, task):
    config = json.loads(configpath.read_text())
    base = Path(config["run_dir"])
    tools = config["tools"]
    if stage == "prepare":
        ref = config["references"][task]
        directory = base / ref["name"] / "reference"
        directory.mkdir(parents=True, exist_ok=True)
        fasta = directory / "genome.fa"
        if not fasta.exists():
            fasta.symlink_to(ref["source"])
        shutil.copy2(ref["source"] + ".fai", str(fasta) + ".fai")
        run([tools["bwa"], "index", str(fasta)])
        print(f"REFERENCE_READY\t{ref['name']}")
    elif stage == "align":
        entry = config["tasks"][task]
        ref, sample = entry["reference"], entry["sample"]
        directory = base / ref["name"] / "coverage"
        directory.mkdir(parents=True, exist_ok=True)
        fasta = base / ref["name"] / "reference/genome.fa"
        cpus = os.environ.get("SLURM_CPUS_PER_TASK", "16")
        with tempfile.TemporaryDirectory(prefix=sample + ".", dir=directory) as tmp:
            temp = Path(tmp)
            raw, dedup, bw = temp / "raw.bam", temp / "dedup.bam", temp / "coverage.bw"
            # Literal backslash-t is the BWA read-group syntax.
            rg = "@RG\\tID:" + sample + "\\tSM:" + sample
            align = [tools["bwa"], "mem", "-M", "-t", cpus, "-R", rg,
                     str(fasta), entry["r1"], entry["r2"]]
            sort = [tools["samtools"], "sort", "-@", "4", "-m", "2G", "-o", str(raw), "-"]
            subprocess.run(["bash", "-o", "pipefail", "-c",
                            shlex.join(align) + " | " + shlex.join(sort)], check=True)
            picard = ([tools["picard"]] if "picard" in tools else
                      [tools["java"], "-Xmx48g", "-jar", tools["picard_jar"]])
            run(picard + ["MarkDuplicates", f"I={raw}", f"O={dedup}",
                          f"M={temp / 'metrics.txt'}", "VALIDATION_STRINGENCY=LENIENT"])
            run([tools["samtools"], "index", "-@", "4", str(dedup)])
            run([tools["bamCoverage"], "--bam", str(dedup), "--outFileName", str(bw),
                 "--outFileFormat", "bigwig", "--binSize", "50", "--numberOfProcessors", cpus])
            with (temp / "flagstat.txt").open("w") as handle:
                subprocess.run([tools["samtools"], "flagstat", str(dedup)], stdout=handle, check=True)
            for source, suffix in ((dedup, ".dedup.bam"), (Path(str(dedup) + ".bai"), ".dedup.bam.bai"),
                                   (bw, ".dedup.50bp.bw"), (temp / "metrics.txt", ".picard_metrics.txt"),
                                   (temp / "flagstat.txt", ".flagstat.txt")):
                source.replace(directory / (sample + suffix))
            (directory / (sample + ".complete.json")).write_text(json.dumps(entry, indent=2) + "\n")
        print(f"ALIGNMENT_COMPLETE\t{ref['name']}\t{sample}")
    else:
        ref = config["references"][task]
        directory = base / ref["name"]
        run([config["python"], str(Path(config["code_dir"]) / "plot_Yb_local_both_references.py"),
             "--reference", ref["name"], "--input-dir", str(directory / "coverage"),
             "--output-dir", str(directory / "plots")])

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-tools", action="store_true")
    parser.add_argument("--submit", action="store_true")
    parser.add_argument("--config", type=Path)
    parser.add_argument("--stage", choices=("prepare", "align", "plot"))
    parser.add_argument("--task", type=int)
    args = parser.parse_args()
    if args.check_tools:
        check_dependencies()
    elif args.submit:
        submit()
    elif args.config and args.stage and args.task is not None:
        execute(args.config, args.stage, args.task)
    else:
        parser.error("Use --submit or --config/--stage/--task")

if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(2)
