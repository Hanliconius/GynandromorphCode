#!/usr/bin/env python3
"""Figure 2C — stage the APPROVED native Smudgeplot PDFs without redrawing them.

The final approved Smudgeplot aesthetic is the stock output from
`smudgeplot all --format pdf -ylim 120`.  This script deliberately does
not re-render .smu files with matplotlib.
"""

from __future__ import annotations
import argparse
from pathlib import Path
import shutil

DEFAULT_ROOT = Path("/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling")
SOURCE_DIR = Path("analysis/smudgeplot_lowerL_20261002/L6_download")
OUT_DIR = Path("analysis/figure2/panels")

SOURCES = [
    ("Fig2C1_diploid_control.pdf", "ERR1051012_diploid_control_L6_smudgeplot.pdf"),
    ("Fig2C2_pseudo_tetraploid.pdf", "pseudoTetraploid_ERR1307111_ERR1307112_L6_smudgeplot.pdf"),
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", default=None)
    parser.add_argument("--output", default=None, help="Ignored; native component PDFs are staged individually.")
    args = parser.parse_args()

    root = Path(args.project_root).expanduser().resolve() if args.project_root else DEFAULT_ROOT
    source_dir = root / SOURCE_DIR
    out_dir = root / OUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)

    # Remove the previously generated custom composite so it cannot be
    # mistaken for the approved native Smudgeplot output.
    bad_composite = out_dir / "Fig2C.pdf"
    if bad_composite.exists():
        bad_composite.unlink()

    staged = []
    for output_name, source_name in SOURCES:
        source = source_dir / source_name
        if not source.is_file() or source.stat().st_size == 0:
            raise FileNotFoundError(
                f"Approved native Smudgeplot PDF not found: {source}"
            )
        destination = out_dir / output_name
        shutil.copy2(source, destination)
        staged.append(destination)
        print(f"NATIVE_SOURCE\t{source}")
        print(f"WROTE\t{destination}")

    print("AESTHETIC\tstock smudgeplot all PDF; no custom matplotlib redraw")
    print(f"COMPONENTS\t{len(staged)}")


if __name__ == "__main__":
    main()
