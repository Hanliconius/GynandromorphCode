#!/usr/bin/env Rscript

options(
    stringsAsFactors = FALSE,
    scipen = 999
)

project_dir <- "/scratch/martinlab/joe/GynandroCalling"

checkpoint_dir <- file.path(
    project_dir,
    "analysis",
    "drop_problematic_reanalysis",
    "checkpoints"
)

figure_dir <- file.path(
    project_dir,
    "analysis",
    "figures",
    "drop_problematic_reanalysis",
    "perfect_D_R_all_comparative_genoplots"
)

comparative_file <- file.path(
    checkpoint_dir,
    "perfect_D_R_all_comparative_genotypes.long.tsv"
)

focal_file <- file.path(
    checkpoint_dir,
    "perfect_D_R_focal_pooled_genotypes.long.tsv"
)

dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

if (!file.exists(comparative_file)) {
    stop(
        paste(
            "Missing comparative genotype table:",
            comparative_file
        )
    )
}

if (!file.exists(focal_file)) {
    stop(
        paste(
            "Missing focal genotype table:",
            focal_file
        )
    )
}

suppressPackageStartupMessages({
    library(ggplot2)
})


comparative <- read.delim(
    comparative_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE
)

focal <- read.delim(
    focal_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE
)

missing_from_focal <- setdiff(
    names(comparative),
    names(focal)
)

missing_from_comparative <- setdiff(
    names(focal),
    names(comparative)
)

if (length(missing_from_focal) > 0) {
    stop(
        paste(
            "Focal table is missing columns:",
            paste(
                missing_from_focal,
                collapse = ", "
            )
        )
    )
}

if (length(missing_from_comparative) > 0) {
    stop(
        paste(
            "Comparative table is missing columns:",
            paste(
                missing_from_comparative,
                collapse = ", "
            )
        )
    )
}

focal <- focal[
    ,
    names(comparative),
    drop = FALSE
]

is_true <- function(x) {
    tolower(
        trimws(
            as.character(x)
        )
    ) %in% c(
        "true",
        "t",
        "1",
        "yes",
        "y"
    )
}

genotype_colours <- c(
    "A/A" = "#3B6FB6",
    "P/A" = "#E69F00",
    "P/P" = "#C51B7D",
    "No call" = "#D9D9D9"
)


write_overview <- function(module_name) {

    comparative_module <- comparative[
        comparative$module == module_name,
        ,
        drop = FALSE
    ]

    focal_module <- focal[
        focal$module == module_name,
        ,
        drop = FALSE
    ]

    if (nrow(comparative_module) == 0) {
        stop(
            paste(
                "No comparative rows for module:",
                module_name
            )
        )
    }

    expected_focal_cells <- (
        length(
            unique(
                comparative_module$site_key
            )
        )
        * 2
    )

    if (nrow(focal_module) != expected_focal_cells) {
        stop(
            paste(
                module_name,
                "has",
                nrow(focal_module),
                "focal cells but expected",
                expected_focal_cells
            )
        )
    }

    trait <- unique(
        comparative_module$trait
    )

    comparative_n <- length(
        unique(
            comparative_module$sample
        )
    )

    site_n <- length(
        unique(
            comparative_module$site_key
        )
    )

    comparative_module$panel_group <- paste0(
        trait,
        " = ",
        comparative_module$phenotype
    )

    focal_module$panel_group <- "Focal gynandromorph"

    comparative_module$display_sample <- (
        comparative_module$sample_label
    )

    focal_module$display_sample <- ifelse(
        focal_module$sample == "TN_Left",
        "TN_Left | pooled left",
        "TN_Right | pooled right"
    )

    x <- rbind(
        focal_module,
        comparative_module
    )

    x$pos <- as.numeric(
        x$pos
    )

    x$display_GT <- factor(
        x$biological_GT,
        levels = c(
            "A/A",
            "P/A",
            "P/P",
            "No call"
        )
    )

    site_manifest <- unique(
        x[
            order(
                x$pos
            ),
            c(
                "site_key",
                "pos",
                "perfect_both_schemes"
            )
        ]
    )

    site_manifest$site_label <- paste0(
        site_manifest$pos,
        ifelse(
            is_true(
                site_manifest$perfect_both_schemes
            ),
            "*",
            ""
        )
    )

    site_label_map <- setNames(
        site_manifest$site_label,
        site_manifest$site_key
    )

    x$site_label <- unname(
        site_label_map[
            x$site_key
        ]
    )

    x$site_label <- factor(
        x$site_label,
        levels = site_manifest$site_label
    )

    focal_samples <- c(
        "TN_Left | pooled left",
        "TN_Right | pooled right"
    )

    comparative_manifest <- unique(
        comparative_module[
            order(
                comparative_module$phenotype,
                comparative_module$taxon,
                comparative_module$data_type,
                comparative_module$sample
            ),
            c(
                "sample",
                "display_sample",
                "phenotype"
            )
        ]
    )

    intended_top_to_bottom <- c(
        focal_samples,
        comparative_manifest$display_sample
    )

    x$display_sample <- factor(
        x$display_sample,
        levels = rev(
            intended_top_to_bottom
        )
    )

    panel_levels <- c(
        "Focal gynandromorph",
        paste0(
            trait,
            " = 0"
        ),
        paste0(
            trait,
            " = 1"
        )
    )

    x$panel_group <- factor(
        x$panel_group,
        levels = panel_levels
    )

    width <- max(
        10,
        min(
            28,
            5 + 0.38 * site_n
        )
    )

    height <- max(
        24,
        4.5 + 0.105 * (
            comparative_n + 2
        )
    )

    plot <- ggplot(
        x,
        aes(
            x = site_label,
            y = display_sample
        )
    ) +

        geom_tile(
            aes(
                fill = display_GT
            ),
            colour = "white",
            linewidth = 0.08,
            width = 0.97,
            height = 0.96
        ) +

        facet_grid(
            rows = vars(
                panel_group
            ),
            scales = "free_y",
            space = "free_y",
            switch = "y",
            drop = TRUE
        ) +

        scale_fill_manual(
            values = genotype_colours,
            drop = FALSE,
            name = "Genotype"
        ) +

        scale_x_discrete(
            name = paste(
                "Hmr18 position",
                "(* perfect under both genotype schemes)"
            ),
            drop = FALSE
        ) +

        scale_y_discrete(
            name = NULL,
            drop = TRUE
        ) +

        labs(
            title = paste0(
                module_name,
                " (",
                trait,
                "): pooled left/right and ",
                comparative_n,
                " comparative samples"
            ),
            subtitle = paste(
                "Perfect SNPs were selected using comparative",
                "butterflies only; focal genotypes are displayed",
                "afterward as an independent comparison."
            ),
            caption = paste(
                "A/A: non-pattern homozygote;",
                "P/A: heterozygote;",
                "P/P: pattern-associated homozygote.",
                "Gray cells failed DP >= 2 or GQ >= 20.",
                "TN_Left and TN_Right were not used to select SNPs",
                "or calculate balanced accuracy."
            )
        ) +

        theme_minimal(
            base_size = 9
        ) +

        theme(
            plot.title = element_text(
                face = "bold",
                size = 12
            ),
            plot.subtitle = element_text(
                size = 9,
                colour = "#444444",
                margin = margin(
                    b = 8
                )
            ),
            plot.caption = element_text(
                size = 8,
                colour = "#555555",
                hjust = 0
            ),
            axis.text.x = element_text(
                angle = 90,
                hjust = 1,
                vjust = 0.5,
                size = 6
            ),
            axis.text.y = element_text(
                size = 4.2,
                colour = "#222222"
            ),
            panel.grid = element_blank(),
            panel.spacing.y = grid::unit(
                0.12,
                "lines"
            ),
            strip.placement = "outside",
            strip.background = element_rect(
                fill = "#F2F2F2",
                colour = NA
            ),
            strip.text.y.left = element_text(
                angle = 0,
                face = "bold",
                size = 8
            ),
            legend.position = "top",
            legend.justification = "left",
            legend.title = element_text(
                face = "bold"
            ),
            plot.margin = margin(
                10,
                12,
                10,
                12
            )
        )

    overview_pdf <- file.path(
        figure_dir,
        paste0(
            "perfect_",
            module_name,
            "_all_comparative_genotypes.overview.pdf"
        )
    )

    overview_svg <- file.path(
        figure_dir,
        paste0(
            "perfect_",
            module_name,
            "_all_comparative_genotypes.overview.svg"
        )
    )

    ggsave(
        filename = overview_pdf,
        plot = plot,
        width = width,
        height = height,
        units = "in",
        device = cairo_pdf,
        limitsize = FALSE
    )

    ggsave(
        filename = overview_svg,
        plot = plot,
        width = width,
        height = height,
        units = "in",
        device = "svg",
        limitsize = FALSE
    )

    data.frame(
        module = module_name,
        comparative_samples = comparative_n,
        focal_samples = length(
            unique(
                focal_module$sample
            )
        ),
        perfect_SNPs = site_n,
        overview_pdf = overview_pdf
    )
}


summary <- rbind(
    write_overview(
        "Dennis"
    ),
    write_overview(
        "Ray"
    )
)

cat(
    "OVERWROTE D/R OVERVIEWS WITH FOCAL LEFT/RIGHT\n\n"
)

print(
    summary,
    row.names = FALSE
)

cat(
    "\nPaginated PDFs were not modified.\n"
)
