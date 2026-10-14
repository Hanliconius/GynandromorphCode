#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)

project_dir <- "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"

dr_script <- file.path(
    project_dir,
    "scripts",
    "115b_plot_D_R_Wallbank_fixed_dominant_adjusted_LR.R"
)

band_script <- file.path(
    project_dir,
    "scripts",
    "142_plot_B_western_postman_vs_western_Amazon_DR_adjusted_LR.R"
)

output_dir <- file.path(
    project_dir,
    "analysis",
    "figures",
    "final_B_D_R_combined"
)

output_pdf <- file.path(
    output_dir,
    "Wallbank_Dennis_Ray_adjusted_LR_plus_restricted_Band_adjusted_LR.pdf"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

required_scripts <- c(dr_script, band_script)
missing_scripts <- required_scripts[!file.exists(required_scripts)]

if (length(missing_scripts) > 0) {
    stop(
        paste(
            "Missing required plotter scripts:",
            paste(missing_scripts, collapse = "\n"),
            sep = "\n"
        )
    )
}

suppressPackageStartupMessages(
    library(grid)
)

load_plot_environment <- function(path) {
    lines <- readLines(path, warn = FALSE)

    device_start <- which(trimws(lines) == "cairo_pdf(")

    if (length(device_start) != 1) {
        stop(
            paste0(
                basename(path),
                ": expected exactly one cairo_pdf() call; found ",
                length(device_start)
            )
        )
    }

    source_text <- paste(
        lines[seq_len(device_start[[1]] - 1)],
        collapse = "\n"
    )

    env <- new.env(parent = globalenv())

    eval(
        parse(
            text = source_text,
            keep.source = TRUE
        ),
        envir = env
    )

    required_objects <- c(
        "ba_grob",
        "focal_grob",
        "draw_combined"
    )

    missing_objects <- required_objects[
        !vapply(
            required_objects,
            exists,
            logical(1),
            envir = env,
            inherits = FALSE
        )
    ]

    if (length(missing_objects) > 0) {
        stop(
            paste(
                basename(path),
                "is missing:",
                paste(missing_objects, collapse = ", ")
            )
        )
    }

    env
}

dr_env <- load_plot_environment(dr_script)
band_env <- load_plot_environment(band_script)

dr_panel <- grid.grabExpr(
    dr_env$draw_combined(),
    wrap = TRUE,
    wrap.grobs = TRUE,
    width = 11.5,
    height = 10.2
)

band_panel <- grid.grabExpr(
    band_env$draw_combined(),
    wrap = TRUE,
    wrap.grobs = TRUE,
    width = 11.5,
    height = 8.2
)

grDevices::cairo_pdf(
    filename = output_pdf,
    width = 11.7,
    height = 18.6,
    onefile = TRUE
)

grid.newpage()

pushViewport(
    viewport(
        layout = grid.layout(
            nrow = 2,
            ncol = 1,
            heights = unit(c(10.2, 8.2), "null")
        )
    )
)

pushViewport(
    viewport(
        layout.pos.row = 1,
        layout.pos.col = 1
    )
)

grid.draw(dr_panel)

grid.text(
    "A",
    x = unit(0.008, "npc"),
    y = unit(0.992, "npc"),
    just = c("left", "top"),
    gp = gpar(
        fontsize = 15,
        fontface = "bold"
    )
)

popViewport()

pushViewport(
    viewport(
        layout.pos.row = 2,
        layout.pos.col = 1
    )
)

grid.draw(band_panel)

grid.text(
    "B",
    x = unit(0.008, "npc"),
    y = unit(0.992, "npc"),
    just = c("left", "top"),
    gp = gpar(
        fontsize = 15,
        fontface = "bold"
    )
)

popViewport(2)

invisible(dev.off())

if (!file.exists(output_pdf)) {
    stop(
        paste(
            "Combined PDF was not created:",
            output_pdf
        )
    )
}

output_size <- file.info(output_pdf)$size

if (!is.finite(output_size) || output_size <= 0) {
    stop(
        paste(
            "Combined PDF is empty:",
            output_pdf
        )
    )
}

cat(
    "COMBINED FINAL BAND / DENNIS / RAY FIGURE\n\n",
    "Panel A: Wallbank Dennis/Ray analysis with reciprocal-mixture focal calls.\n",
    "Panel B: restricted Band analysis 142.\n\n",
    "PDF:\n",
    output_pdf,
    "\n\nBytes: ",
    output_size,
    "\n",
    sep = ""
)
