#!/usr/bin/env python3
from pathlib import Path
import shutil


project = Path(__file__).resolve().parents[1]
target = project / "scripts/162_extend_founder_mixing_visuals.py"
backup = project / "scripts/162_extend_founder_mixing_visuals.pre_contrast_trends.py"

if not target.is_file():
    raise RuntimeError(f"Missing target script: {target}")

text = target.read_text()
marker = 'primary["contrast_bin"] = np.digitize('

if marker in text:
    print("Genotype-contrast trend update already present; no edit made.")
else:
    old_calibration = '''    grouped = primary.groupby(x_column)[y_column]
    x = np.array(sorted(primary[x_column].unique()))
    median = grouped.median().reindex(x).to_numpy()
    low = grouped.quantile(0.025).reindex(x).to_numpy()
    high = grouped.quantile(0.975).reindex(x).to_numpy()
    axis.fill_between(x, low, high, color=color, alpha=0.13, linewidth=0)
    axis.plot(x, median, color=color, linewidth=2, label="high-call median")
'''
    new_calibration = '''    grouped = primary.groupby(x_column)[y_column]
    x = np.array(sorted(primary[x_column].unique()))
    low = grouped.quantile(0.025).reindex(x).to_numpy()
    high = grouped.quantile(0.975).reindex(x).to_numpy()
    axis.fill_between(x, low, high, color=color, alpha=0.13, linewidth=0)

    contrast_edges = np.quantile(
        primary["mean_genotype_contrast"], [0.0, 1 / 3, 2 / 3, 1.0]
    )
    if np.any(np.diff(contrast_edges) <= 0):
        raise RuntimeError(
            f"Could not form three genotype-contrast bins for {founder}."
        )
    primary["contrast_bin"] = np.digitize(
        primary["mean_genotype_contrast"], contrast_edges[1:-1], right=True
    )
    bin_names = ["Low", "Middle", "High"]
    colour_map = plt.get_cmap(CMAP)
    for bin_index, bin_name in enumerate(bin_names):
        bin_data = primary.loc[primary["contrast_bin"].eq(bin_index)].copy()
        bin_grouped = bin_data.groupby(x_column)[y_column]
        bin_x = np.array(sorted(bin_data[x_column].unique()))
        bin_median = bin_grouped.median().reindex(bin_x).to_numpy()
        bin_colour = colour_map(
            contrast_norm(
                bin_data["mean_genotype_contrast"].median()
            )
        )
        axis.plot(
            bin_x,
            bin_median,
            color=bin_colour,
            linewidth=2.2,
            marker="o",
            markersize=3.2,
            label=(
                f"{bin_name} contrast "
                f"({bin_data['mean_genotype_contrast'].min():.3f}-"
                f"{bin_data['mean_genotype_contrast'].max():.3f})"
            ),
        )
'''
    old_subtitle = '''        "High-call simulations are coloured by mean genotype contrast; lower-call pairs remain as grey context.",
'''
    new_subtitle = '''        "High-call simulations are split into low, middle and high genotype-contrast thirds; "
        "each coloured line is the bin-specific median.",
'''
    old_legend = '''    axes[0, 0].legend(frameon=False, fontsize=6.5, loc="upper left")
'''
    new_legend = '''    axes[0, 0].legend(
        frameon=False,
        fontsize=6.2,
        loc="upper left",
        ncol=2,
        columnspacing=0.8,
        handlelength=1.8,
    )
'''

    replacements = [
        (
            "calibration block",
            old_calibration,
            new_calibration,
        ),
        (
            "subtitle",
            old_subtitle,
            new_subtitle,
        ),
        (
            "legend",
            old_legend,
            new_legend,
        ),
    ]

    for label, old, new in replacements:
        count = text.count(old)
        if count != 1:
            raise RuntimeError(
                f"Expected exactly one {label} anchor; "
                f"found {count}. No edit made."
            )

    if not backup.exists():
        shutil.copy2(
            target,
            backup,
        )

    for _, old, new in replacements:
        text = text.replace(
            old,
            new,
            1,
        )

    temporary = target.with_suffix(
        target.suffix + ".tmp"
    )
    temporary.write_text(text)
    temporary.replace(target)

    print(f"Updated: {target}")
    print(
        "Recoverable pre-update copy: "
        f"{backup}"
    )

final_text = target.read_text()

if final_text.count(marker) != 1:
    raise RuntimeError(
        "Final genotype-contrast trend "
        "marker count is not one."
    )

if final_text.count("bin-specific median") != 1:
    raise RuntimeError(
        "Final bin-specific subtitle "
        "marker count is not one."
    )

print("Trend update validation passed.")
