#!/bin/bash
#SBATCH --job-name=smudgeLowL
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=1-00:00:00
#SBATCH --array=0-5%3
#SBATCH --output=logs/smudgeLowL_%A_%a.out
#SBATCH --error=logs/smudgeLowL_%A_%a.err

set -euo pipefail
export PATH="$HOME/.local/share/mamba/envs/smudgeplot/bin:$PATH"
export MPLBACKEND=Agg
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

python -u - <<'PY'
import os
import pathlib
import shutil
import struct
import subprocess
import tempfile

P = pathlib.Path(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
)
datasets = [
    "pseudoTetraploid_ERR1307111_ERR1307112",
    "TN_T2_merge",
    "TN_T3_merge",
]
task = int(os.environ["SLURM_ARRAY_TASK_ID"])
name = datasets[task // 2]
lower = [4, 6][task % 2]

db = P / "analysis/smudgeplot_arrangements" / name / "FastK_Table"
out = P / "analysis/smudgeplot_lowerL_20261002" / name / f"L{lower}"
out.mkdir(parents=True, exist_ok=True)
os.chdir(out)
os.environ["MPLCONFIGDIR"] = str(out / "mplconfig")

for tool in ["smudgeplot", "Histex"]:
    if not shutil.which(tool):
        raise RuntimeError(f"{tool} missing from smudgeplot environment")

# Check the database's retained count threshold and all table parts.
with db.with_suffix(".ktab").open("rb") as f:
    k, parts, minimum, prefix = struct.unpack("=4i", f.read(16))

if not (5 <= k <= 1000 and 1 <= parts <= 10000
        and 1 <= minimum <= 32767):
    raise RuntimeError(f"Unexpected FastK header: {k, parts, minimum, prefix}")

if minimum > lower:
    raise RuntimeError(
        f"{name}: database only retains counts >= {minimum}; "
        f"L={lower} requires recounting."
    )

for i in range(1, parts + 1):
    shard = db.parent / f".{db.name}.ktab.{i}"
    if not shard.is_file():
        raise FileNotFoundError(shard)

print(f"{name}: k={k}; database minimum={minimum}; L={lower}", flush=True)
pathlib.Path("database_info.txt").write_text(
    f"database={db}\nk={k}\nparts={parts}\nminimum={minimum}\nL={lower}\n"
)

def run(args, log, stdout=None):
    args = list(map(str, args))
    print("RUN", " ".join(args), flush=True)
    with open(log, "w") as err:
        if stdout:
            with open(stdout, "w") as dest:
                subprocess.run(args, stdout=dest, stderr=err, check=True)
        else:
            subprocess.run(
                args, stdout=err, stderr=subprocess.STDOUT, check=True
            )

if pathlib.Path("COMPLETE").is_file():
    print("Already complete:", out, flush=True)
else:
    run(["smudgeplot", "-v"], "version.txt")
    run(
        ["Histex", "-A", "-h1:300", db],
        "histex.log", "FastK_histogram.tsv"
    )

    # Only publish the new SMU after successful extraction.
    if not pathlib.Path("kmerpairs.smu").is_file():
        with tempfile.TemporaryDirectory(prefix="extract_", dir=out) as tmp:
            root = pathlib.Path(tmp) / "kmerpairs"
            run(
                [
                    "smudgeplot", "hetmers",
                    "-L", lower,
                    "-t", os.environ["SLURM_CPUS_PER_TASK"],
                    "-o", root,
                    "-tmp", tmp,
                    "--verbose", db,
                ],
                "hetmers.log"
            )
            smu = root.with_suffix(".smu")
            if not smu.is_file() or smu.stat().st_size == 0:
                raise RuntimeError("No usable SMU; inspect hetmers.log")
            smu.replace("kmerpairs.smu")

    # Keep diagnostics showing exactly what the error mask excludes.
    for label, extra in [
        ("unmasked", []),
        ("masked", ["--mask_errors"]),
    ]:
        run(
            [
                "smudgeplot", "peak_aggregation",
                "-d", "2", "-nf", "1000",
                *extra, "kmerpairs.smu",
            ],
            f"{label}_peaks.log",
            f"{label}_peaks.tsv"
        )

    # Infer coverage automatically; do not force a ploidy interpretation.
    run(
        [
            "smudgeplot", "all",
            "--format", "pdf",
            "--json_report",
            "-ylim", "120",
            "-t", f"{name}; L={lower}; automatic coverage",
            "-o", "output",
            "kmerpairs.smu",
        ],
        "smudgeplot_all.log"
    )

    expected = [
        "output_smudgeplot.pdf",
        "output_smudgeplot_log10.pdf",
        "output.smudge_report.tsv",
    ]
    if not all(
        pathlib.Path(f).is_file() and pathlib.Path(f).stat().st_size
        for f in expected
    ):
        raise RuntimeError("Missing outputs; inspect smudgeplot_all.log")

    pathlib.Path("COMPLETE").write_text(
        "Extraction and automatic Smudgeplot completed; review required.\n"
    )
    print("COMPLETE:", out, flush=True)
PY
