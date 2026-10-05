#!/usr/bin/env python3

from pathlib import Path


script = Path(
    "scripts/110d_plot_B_no_DR_positive_BA_and_focal.R"
)

text = script.read_text()


def replace_required(
    text,
    old,
    new,
    label,
    count=1,
):
    if old not in text:
        raise RuntimeError(
            f"Required block not found: {label}"
        )

    return text.replace(
        old,
        new,
        count,
    )


# ------------------------------------------------------------
# Restrict the physical interval.
# ------------------------------------------------------------

text = replace_required(
    text,
    '''x_limits <- range(
    d$position_Mb,
    finite = TRUE
)
''',
    '''# Requested display interval:
# Hmr18:1,060,000-1,240,000
x_limits <- c(
    1060000 / 1000000,
    1240000 / 1000000
)

# optix annotation in the same Hmr18 coordinate system.
optix_start_Mb <- 1088602 / 1000000
optix_end_Mb <- 1097846 / 1000000
optix_midpoint_Mb <- (
    optix_start_Mb
    + optix_end_Mb
) / 2
''',
    "x-axis limits",
)


# ------------------------------------------------------------
# Add optix shading and gene bar to the BA panel.
# ------------------------------------------------------------

text = replace_required(
    text,
    '''ba_plot <- ggplot() +
''',
    '''ba_plot <- ggplot() +

    annotate(
        "rect",
        xmin = optix_start_Mb,
        xmax = optix_end_Mb,
        ymin = -Inf,
        ymax = Inf,
        fill = "#56B4E9",
        alpha = 0.12
    ) +

    annotate(
        "segment",
        x = optix_start_Mb,
        xend = optix_end_Mb,
        y = 1.017,
        yend = 1.017,
        colour = "#276B8A",
        linewidth = 1.25,
        lineend = "butt"
    ) +

    annotate(
        "segment",
        x = optix_start_Mb,
        xend = optix_start_Mb,
        y = 1.012,
        yend = 1.022,
        colour = "#276B8A",
        linewidth = 0.55
    ) +

    annotate(
        "segment",
        x = optix_end_Mb,
        xend = optix_end_Mb,
        y = 1.012,
        yend = 1.022,
        colour = "#276B8A",
        linewidth = 0.55
    ) +

    annotate(
        "text",
        x = optix_midpoint_Mb,
        y = 1.027,
        label = "optix",
        colour = "#194B63",
        fontface = "italic",
        size = 3.4,
        vjust = 0
    ) +
''',
    "BA-panel optix annotation",
)


# Make vertical room for the annotation.
text = replace_required(
    text,
    '''ylim = c(
            0.48,
            1.02
        )''',
    '''ylim = c(
            0.48,
            1.045
        )''',
    "BA-panel upper y limit",
)


# ------------------------------------------------------------
# Add matching optix shading to the focal genotype panel.
# ------------------------------------------------------------

text = replace_required(
    text,
    '''focal_plot <- ggplot() +
''',
    '''focal_plot <- ggplot() +

    annotate(
        "rect",
        xmin = optix_start_Mb,
        xmax = optix_end_Mb,
        ymin = -Inf,
        ymax = Inf,
        fill = "#56B4E9",
        alpha = 0.12
    ) +
''',
    "focal-panel optix shading",
)


# Add coordinate information to the figure subtitle.
text = replace_required(
    text,
    '''subtitle = paste(
            "Comparative balanced accuracy was calculated without",
            "TN_Left or TN_Right; focal genotypes are displayed below."
        )''',
    '''subtitle = paste(
            "Hmr18:1.060-1.240 Mb.",
            "Comparative balanced accuracy was calculated without",
            "TN_Left or TN_Right; focal genotypes are displayed below."
        )''',
    "figure subtitle",
)


script.write_text(
    text
)

print(
    f"Patched: {script}"
)
print(
    "Display interval: Hmr18:1,060,000-1,240,000"
)
print(
    "optix annotation: Hmr18:1,088,602-1,097,846"
)
print(
    "The timareta-restored 111 analysis was not modified."
)
