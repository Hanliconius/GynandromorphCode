#!/usr/bin/env python3

from pathlib import Path
import shutil

script = Path(
    "scripts/118_plot_simple_optix_comparative_and_focal_pattern.R"
)
backup = Path(
    "scripts/118_plot_simple_optix_comparative_and_focal_pattern.R"
    ".before_krista_windows"
)

if not script.is_file():
    raise RuntimeError(f"Missing plotting script: {script}")

text = script.read_text()

definition_tag = "# KRISTA_PLOT_WINDOWS_DEFINITION"
layer_tag = "# KRISTA_PLOT_WINDOWS_LAYER"

if definition_tag in text and layer_tag in text:
    print("Krista-window overlay is already present; no patch required.")

else:
    if not backup.exists():
        shutil.copy2(script, backup)
        print(f"Backup created: {backup}")

    definition = r'''
# KRISTA_PLOT_WINDOWS_DEFINITION
#
# Original regional windows used for Krista's genotype plots.
# Coordinates are converted to Mb because the aligned figure uses Mb.
#
krista_regions <- data.frame(
    module = c(
        "Band",
        "Dennis",
        "Dennis",
        "Ray",
        "Ray"
    ),
    region_label = c(
        "B",
        "D1",
        "D2",
        "R1",
        "R2"
    ),
    xmin_mb = c(
        1112651,
        1215210,
        1218557,
        1156746,
        1174547
    ) / 1e6,
    xmax_mb = c(
        1114314,
        1216124,
        1221494,
        1160220,
        1187833
    ) / 1e6,
    stringsAsFactors = FALSE
)

krista_regions$mid_mb <- (
    krista_regions$xmin_mb +
    krista_regions$xmax_mb
) / 2

'''

    ba_marker = "ba_plot <- ggplot() +"

    if ba_marker not in text:
        raise RuntimeError(
            "Could not find the ba_plot construction; no change was made."
        )

    text = text.replace(
        ba_marker,
        definition
        + ba_marker
        + r'''

    # KRISTA_PLOT_WINDOWS_LAYER
    geom_rect(
        data = krista_regions,
        aes(
            xmin = xmin_mb,
            xmax = xmax_mb,
            ymin = -Inf,
            ymax = Inf
        ),
        inherit.aes = FALSE,
        fill = "#6F6F6F",
        colour = "#505050",
        alpha = 0.13,
        linewidth = 0.28
    ) +

    geom_text(
        data = krista_regions,
        aes(
            x = mid_mb,
            y = 1.008,
            label = region_label
        ),
        inherit.aes = FALSE,
        colour = "#3A3A3A",
        size = 2.7,
        fontface = "bold",
        vjust = 0.5
    ) +''',
        1,
    )

    focal_marker = "focal_plot <- ggplot() +"

    if focal_marker not in text:
        raise RuntimeError(
            "Could not find the focal_plot construction; no change was made."
        )

    text = text.replace(
        focal_marker,
        focal_marker
        + r'''

    # KRISTA_PLOT_WINDOWS_LAYER
    geom_rect(
        data = krista_regions,
        aes(
            xmin = xmin_mb,
            xmax = xmax_mb,
            ymin = -Inf,
            ymax = Inf
        ),
        inherit.aes = FALSE,
        fill = "#6F6F6F",
        colour = "#505050",
        alpha = 0.13,
        linewidth = 0.28
    ) +''',
        1,
    )

    script.write_text(text)
    print(f"Patched: {script}")

print()
print("Krista windows:")
print("Band     Hmr18:1112651-1114314")
print("Dennis1  Hmr18:1215210-1216124")
print("Dennis2  Hmr18:1218557-1221494")
print("Ray1     Hmr18:1156746-1160220")
print("Ray2     Hmr18:1174547-1187833")
