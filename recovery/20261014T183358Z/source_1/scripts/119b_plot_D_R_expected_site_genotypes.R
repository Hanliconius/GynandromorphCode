#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
})

project <- "/scratch/martinlab/joe/GynandroCalling"

input_file <- file.path(
    project,
    "analysis/original_502/checkpoints",
    "optix_D_R_expected_focal_site_genotypes.long.tsv"
)

figure_dir <- file.path(
    project,
    "analysis/figures",
    "optix_D_R_expected_site_genoplots"
)

dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

if (!file.exists(input_file)) {
    stop(
        paste(
            "Missing genotype table:",
            input_file
        )
    )
}

d <- read.delim(
    input_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

genotype_colours <- c(
    "A/A" = "#B2188B",
    "P/A" = "#F39C34",
    "P/P" = "#2B9CC4",
    "N/QC" = "#E2E2E2"
)

make_plot <- function(module_name) {

    x <- d[
        d$module == module_name,
        ,
        drop = FALSE
    ]

    x$pos <- as.integer(x$pos)

    site_order <- unique(
        x[
            order(x$pos),
            c("pos", "site_label")
        ]
    )

    x$site_label <- factor(
        x$site_label,
        levels = site_order$site_label
    )

    focal_labels <- unique(
        x$sample_label[
            x$sample_type == "focal"
        ]
    )

    comparative_rows <- unique(
        x[
            x$sample_type == "comparative",
            c(
                "sample_label",
                "row_group",
                "phenotype",
                "taxon_name",
                "sample"
            )
        ]
    )

    comparative_rows$phenotype_order <- ifelse(
        comparative_rows$phenotype == "1",
        1,
        ifelse(
            comparative_rows$phenotype == "0",
            2,
            3
        )
    )

    comparative_rows <- comparative_rows[
        order(
            comparative_rows$phenotype_order,
            comparative_rows$taxon_name,
            comparative_rows$sample
        ),
        ,
        drop = FALSE
    ]

    sample_order <- c(
        focal_labels,
        comparative_rows$sample_label
    )

    x$sample_label <- factor(
        x$sample_label,
        levels = rev(sample_order)
    )

    focal_group <- "Focal gynandromorph"

    comparative_groups <- unique(
        comparative_rows$row_group
    )

    x$row_group <- factor(
        x$row_group,
        levels = c(
            focal_group,
            comparative_groups
        )
    )

    site_count <- length(
        unique(x$site_label)
    )

    comparative_count <- length(
        unique(
            x$sample[
                x$sample_type == "comparative"
            ]
        )
    )

    ggplot(
        x,
        aes(
            x = site_label,
            y = sample_label,
            fill = display_GT
        )
    ) +

        geom_tile(
            colour = NA,
            width = 0.96,
            height = 0.96
        ) +

        facet_grid(
            rows = vars(row_group),
            scales = "free_y",
            space = "free_y",
            switch = "y"
        ) +

        scale_fill_manual(
            name = "Biological genotype",
            values = genotype_colours,
            breaks = c(
                "A/A",
                "P/A",
                "P/P",
                "N/QC"
            ),
            drop = FALSE
        ) +

        scale_x_discrete(
            position = "top",
            expand = expansion(
                mult = c(0.002, 0.002)
            )
        ) +

        labs(
            title = paste0(
                module_name,
                " comparative-associated SNP genotypes"
            ),
            subtitle = paste0(
                site_count,
                " comparative-selected SNPs that independently ",
                "show TN_Left=P/P and TN_Right=P/A; ",
                comparative_count,
                " comparative individuals"
            ),
            x = "Hmr18 position",
            y = NULL,
            caption = paste(
                "P denotes the allele associated with pattern presence",
                "in the comparative analysis.",
                "Calls require DP >= 2 and GQ >= 20.",
                "Focal genotypes were not used for comparative SNP selection."
            )
        ) +

        theme_minimal(
            base_size = 8
        ) +

        theme(
            panel.grid = element_blank(),

            axis.text.x = element_text(
                angle = 90,
                hjust = 0,
                vjust = 0.5,
                size = 6
            ),

            axis.text.y = element_text(
                size = 2.8,
                colour = "#222222"
            ),

            axis.title.x = element_text(
                margin = margin(
                    t = 6
                )
            ),

            strip.placement = "outside",

            strip.text.y.left = element_text(
                angle = 0,
                face = "bold",
                size = 8
            ),

            strip.background = element_rect(
                fill = "#F2F2F2",
                colour = "#B0B0B0",
                linewidth = 0.35
            ),

            legend.position = "top",
            legend.justification = "left",

            legend.title = element_text(
                face = "bold"
            ),

            plot.title = element_text(
                face = "bold",
                size = 13
            ),

            plot.subtitle = element_text(
                size = 9,
                colour = "#444444"
            ),

            plot.caption = element_text(
                size = 7,
                colour = "#4A4A4A",
                hjust = 0
            ),

            plot.margin = margin(
                10,
                12,
                10,
                12
            )
        )
}


dennis_plot <- make_plot("Dennis")
ray_plot <- make_plot("Ray")

dennis_pdf <- file.path(
    figure_dir,
    "optix_Dennis_expected_focal_site_genotypes.pdf"
)

ray_pdf <- file.path(
    figure_dir,
    "optix_Ray_expected_focal_site_genotypes.pdf"
)

combined_pdf <- file.path(
    figure_dir,
    "optix_D_R_expected_focal_site_genotypes.pdf"
)

ggsave(
    filename = dennis_pdf,
    plot = dennis_plot,
    device = cairo_pdf,
    width = 15,
    height = 28,
    units = "in",
    limitsize = FALSE
)

ggsave(
    filename = ray_pdf,
    plot = ray_plot,
    device = cairo_pdf,
    width = 15,
    height = 28,
    units = "in",
    limitsize = FALSE
)

cairo_pdf(
    filename = combined_pdf,
    width = 15,
    height = 28,
    onefile = TRUE
)

print(dennis_plot)
print(ray_plot)

invisible(dev.off())

cat("D/R EXPECTED-SITE SUPPLEMENTARY GENOTYPE PLOTS\n\n")
cat("Dennis PDF:\n", dennis_pdf, "\n\n")
cat("Ray PDF:\n", ray_pdf, "\n\n")
cat("Combined two-page PDF:\n", combined_pdf, "\n")
