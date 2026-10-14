#!/usr/bin/env python3

from pathlib import Path

script = Path(
    "/scratch/martinlab/joe/GynandroCalling/"
    "scripts/116b_plot_Wallbank_Dennis_Ray_BA.R"
)

backup = script.with_suffix(
    script.suffix + ".before_hline_fix"
)

old = '''    geom_hline(
        yintercept = c(
            0.50,
            0.80,
            0.90,
            1.00
        ),
        colour = c(
            "#888888",
            "#999999",
            "#999999",
            "#555555"
        ),
        linewidth = c(
            0.35,
            0.25,
            0.25,
            0.35
        ),
        linetype = c(
            "solid",
            "dotted",
            "dotted",
            "solid"
        )
    ) +
'''

new = '''    geom_hline(
        yintercept = 0.50,
        colour = "#888888",
        linewidth = 0.35,
        linetype = "solid"
    ) +

    geom_hline(
        yintercept = 0.80,
        colour = "#999999",
        linewidth = 0.25,
        linetype = "dotted"
    ) +

    geom_hline(
        yintercept = 0.90,
        colour = "#999999",
        linewidth = 0.25,
        linetype = "dotted"
    ) +

    geom_hline(
        yintercept = 1.00,
        colour = "#555555",
        linewidth = 0.35,
        linetype = "solid"
    ) +
'''

if not script.is_file():
    raise RuntimeError(f"Missing plotting script: {script}")

text = script.read_text()

match_count = text.count(old)

if match_count != 1:
    raise RuntimeError(
        "Expected exactly one vectorized geom_hline block; "
        f"found {match_count}. No file was changed."
    )

if not backup.exists():
    backup.write_text(text)
    print(f"Backup created: {backup}")

script.write_text(
    text.replace(old, new, 1)
)

print(f"Patched: {script}")
