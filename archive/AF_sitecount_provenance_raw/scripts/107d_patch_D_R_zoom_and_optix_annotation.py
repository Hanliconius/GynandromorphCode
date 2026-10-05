#!/usr/bin/env python3

from pathlib import Path
import shutil
import subprocess


project = Path("/scratch/martinlab/joe/GynandroCalling")

target = (
    project
    / "scripts"
    / "107c_plot_BA_with_separate_focal_genotype_tracks.R"
)

backup = target.with_name(
    target.name + ".before_zoom_optix"
)

marker = "# OPTIX_ZOOM_1060_1240_AND_GENE_ANNOTATION"


if not target.is_file():
    raise RuntimeError(
        f"Missing D/R plotting script: {target}"
    )


text = target.read_text()


if marker in text:
    print(
        "D/R zoom and optix annotation are already present."
    )

else:
    required_anchors = [
        "ba_plot <- ggplot() +",
        "focal_plot <- ggplot() +",
        """    coord_cartesian(
        ylim = c(
            0.48,
            1.02
        ),
        clip = "off"
    ) +""",
        """        subtitle = paste(
            "Comparative balanced accuracy was calculated without",
            "TN_Left or TN_Right; focal genotypes are displayed below."
        )""",
    ]

    missing = [
        anchor
        for anchor in required_anchors
        if anchor not in text
    ]

    if missing:
        raise RuntimeError(
            "Expected plotting anchors were not found. "
            "No file was changed.\n\nMissing anchors:\n"
            + "\n\n".join(missing)
        )

    if not backup.exists():
        shutil.copy2(
            target,
            backup,
        )
        print(f"Backup created: {backup}")

    # --------------------------------------------------------
    # Override the broad plotting limits immediately before
    # constructing the plots. Positions in this script are Mb.
    # --------------------------------------------------------

    ba_anchor = "ba_plot <- ggplot() +"

    ba_replacement = f"""{marker}

# Presentation-only zoom and gene annotation.
# These values are in Mb because the plotting tables use pos_mb.
x_limits <- c(
    1.060,
    1.240
)

optix_start_mb <- 1088602 / 1e6
optix_end_mb <- 1097846 / 1e6
optix_midpoint_mb <- (
    optix_start_mb + optix_end_mb
) / 2

ba_plot <- ggplot() +

    annotate(
        "rect",
        xmin = optix_start_mb,
        xmax = optix_end_mb,
        ymin = -Inf,
        ymax = Inf,
        fill = "#56B4E9",
        alpha = 0.11
    ) +

    annotate(
        "segment",
        x = optix_start_mb,
        xend = optix_end_mb,
        y = 1.027,
        yend = 1.027,
        colour = "#2166AC",
        linewidth = 2.2,
        lineend = "butt"
    ) +

    annotate(
        "text",
        x = optix_midpoint_mb,
        y = 1.038,
        label = "optix",
        colour = "#2166AC",
        fontface = "italic",
        size = 3.5
    ) +"""

    text = text.replace(
        ba_anchor,
        ba_replacement,
        1,
    )

    # Add the same lightly shaded interval behind the focal tracks,
    # but omit a second label/bar to keep that panel clean.
    focal_anchor = "focal_plot <- ggplot() +"

    focal_replacement = """focal_plot <- ggplot() +

    annotate(
        "rect",
        xmin = optix_start_mb,
        xmax = optix_end_mb,
        ymin = -Inf,
        ymax = Inf,
        fill = "#56B4E9",
        alpha = 0.11
    ) +"""

    text = text.replace(
        focal_anchor,
        focal_replacement,
        1,
    )

    # Make enough vertical room for the optix bar and label.
    old_coord = """    coord_cartesian(
        ylim = c(
            0.48,
            1.02
        ),
        clip = "off"
    ) +"""

    new_coord = """    coord_cartesian(
        ylim = c(
            0.48,
            1.045
        ),
        clip = "off"
    ) +"""

    text = text.replace(
        old_coord,
        new_coord,
        1,
    )

    # Make the displayed range explicit.
    old_subtitle = """        subtitle = paste(
            "Comparative balanced accuracy was calculated without",
            "TN_Left or TN_Right; focal genotypes are displayed below."
        )"""

    new_subtitle = """        subtitle = paste(
            "Hmr18:1,060,000-1,240,000.",
            "Comparative balanced accuracy was calculated without",
            "TN_Left or TN_Right; focal genotypes are displayed below."
        )"""

    text = text.replace(
        old_subtitle,
        new_subtitle,
        1,
    )

    target.write_text(text)

    print(f"Patched: {target}")
    print(
        "Zoom: Hmr18:1,060,000-1,240,000"
    )
    print(
        "optix annotation: Hmr18:1,088,602-1,097,846"
    )


rscript = Path(
    "/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/"
    ".local/share/mamba/envs/gynandrocalling/bin/Rscript"
)

syntax = subprocess.run(
    [
        str(rscript),
        "-e",
        (
            f'parse(file="{target}"); '
            'cat("R syntax check passed.\\n")'
        ),
    ],
    text=True,
    capture_output=True,
)

print(syntax.stdout, end="")

if syntax.returncode != 0:
    print(syntax.stderr, end="")
    raise RuntimeError(
        "The patched R script did not pass its syntax check."
    )
