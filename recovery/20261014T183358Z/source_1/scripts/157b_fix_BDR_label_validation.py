#!/usr/bin/env python3

from pathlib import Path

project = Path(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
)

target = project / (
    "scripts/"
    "157_plot_population_frequency_plus_LR_BWPK_exact_BA_"
    "nonfixed_ref_taxon_canonical_BDR_labels.R"
)

if not target.is_file():
    raise RuntimeError(
        f"Missing Analysis-157 plotter: {target}"
    )

text = target.read_text()


old = """    observed_populations <- sort(
        unique(
            comparative$population_label
        )
    )"""


new = """    # Validate the canonical taxon identity, not the decorated
    # B/D/R display label. population_key remains exactly equal
    # to the canonical taxon designation from the BA manifest.
    observed_populations <- sort(
        unique(
            comparative$population_key
        )
    )"""


count = text.count(old)

if count != 1:
    raise RuntimeError(
        "Expected exactly one population-set validation block "
        f"using population_label; found {count}."
    )

text = text.replace(
    old,
    new,
    1,
)

target.write_text(text)

print(f"Patched: {target}")
print(
    "Population-set validation now uses canonical population_key; "
    "population_label remains free to contain B/D/R display tags."
)
