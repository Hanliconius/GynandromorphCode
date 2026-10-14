#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
})

project <- "/scratch/martinlab/joe/GynandroCalling"

input_file <- file.path(
    project,
    "analysis/original_502/checkpoints",
    "optix_R_all_BA90_comparative_selected_genotypes.long.tsv"
)

figure_dir <- file.path(
    project,
    "analysis/figures",
    "optix_R_all_BA90_comparative_selected_genoplot"
)

dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

d <- read.delim(
    input_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

site_information <- unique(
    d[
        ,
        c(
            "pos",
            "site_label",
            "focal_outcome"
        )
    ]
)

site_information <- site_information[
    order(site_information$pos),
    ,
    drop = FALSE
]

d$site_label <- factor(
    d$site_label,
    levels = site_information$site_label
)

focal_rows <- unique(
    d[
        d$sample_type == "focal",
        c("sample", "sample_label")
    ]
)

focal_rows$order <- match(
    focal_rows$sample,
    c("TN_Left", "TN_Right")
)

focal_rows <- focal_rows[
    order(focal_rows$order),
    ,
    drop = FALSE
]

comparative_rows <- unique(
    d[
        d$sample_type == "comparative",
        c(
            "sample",
            "sample_label",
            "phenotype",
            "taxon_name"
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
    focal_rows$sample_label,
    comparative_rows$sample_label
)

d$sample_label <- factor(
    d$sample_label,
    levels = rev(sample_order)
)

d$row_group <- factor(
    d$row_group,
    levels = c(
        "Focal gynandromorph",
        "Ray present",
        "Ray absent",
        "Phenotype unavailable"
    )
)

genotype_colours <- c(
    "A/A" = "#B2188B",
    "P/A" = "#F39C34",
    "P/P" = "#2B9CC4",
    "N/QC" = "#E2E2E2"
)

site_count <- length(
    unique(d$site_label)
)

sample_count <- length(
    unique(
        d$sample[
            d$sample_type == "comparative"
        ]
    )
)

plot_width <- max(
    15,
    min(
        30,
        5 + 0.16 * site_count
    )
)

p <- ggplot(
    d,
    aes(
        x = site_label,
        y = sample_label,
        fill = display_GT
    )
) +

    geom_tile(
        width = 0.96,
        height = 0.96,
        colour = NA
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
            mult = c(0.001, 0.001)
        )
    ) +

    labs(
        title = "Broad optix Ray-associated BA90 genotype plot",
        subtitle = paste0(
            site_count,
            " comparative-selected Ray SNPs with minimum BA >= 0.90; ",
            sample_count,
            " comparative individuals. ",
            "Focal L/R agreement was not required."
        ),
        x = "Hmr18 position",
        y = NULL,
        caption = paste(
            "P is the Ray-associated allele and A is the",
            "non-Ray allele.",
            "For the dominant Ray model, P/A and P/P predict",
            "Ray presence.",
            "Calls require DP >= 2 and GQ >= 20.",
            "Comparative SNP selection did not use TN_Left",
            "or TN_Right."
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
            size = 5.5
        ),

        axis.text.y = element_text(
            size = 2.8,
            colour = "#222222"
        ),

        strip.placement = "outside",

        strip.text.y.left = element_text(
            angle = 0,
            face = "bold",
            size = 8
        ),

        strip.background = element_rect(
            fill = "#F2F2F2",
            colour = "#A0A0A0",
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
            hjust = 0,
            colour = "#444444"
        ),

        plot.margin = margin(
            10,
            12,
            10,
            12
        )
    )

pdf_file <- file.path(
    figure_dir,
    "optix_R_all_BA90_comparative_selected_genotypes.pdf"
)

ggsave(
    filename = pdf_file,
    plot = p,
    device = cairo_pdf,
    width = plot_width,
    height = 28,
    units = "in",
    limitsize = FALSE
)

cat("BROAD RAY BA90 COMPARATIVE-SELECTED GENOTYPE PLOT\n\n")
cat("Sites:", site_count, "\n")
cat("Comparative samples:", sample_count, "\n")
cat("PDF width:", plot_width, "inches\n")
cat("PDF:\n", pdf_file, "\n")
