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

dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

input_file <- file.path(
    checkpoint_dir,
    "perfect_D_R_all_comparative_genotypes.long.tsv"
)

if (!file.exists(input_file)) {
    stop(
        paste(
            "Missing genotype table:",
            input_file
        )
    )
}

suppressPackageStartupMessages({
    library(ggplot2)
})

d <- read.delim(
    input_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE
)

d$pos <- as.numeric(
    d$pos
)

d$phenotype_label <- paste0(
    d$trait,
    " = ",
    d$phenotype
)

d$display_GT <- factor(
    d$biological_GT,
    levels = c(
        "A/A",
        "P/A",
        "P/P",
        "No call"
    )
)

d$cell_label <- c(
    "A/A" = "AA",
    "P/A" = "PA",
    "P/P" = "PP",
    "No call" = "-"
)[
    as.character(
        d$display_GT
    )
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

d$site_label <- paste0(
    d$pos,
    ifelse(
        is_true(
            d$perfect_both_schemes
        ),
        "*",
        ""
    )
)

genotype_colours <- c(
    "A/A" = "#3B6FB6",
    "P/A" = "#E69F00",
    "P/P" = "#C51B7D",
    "No call" = "#D9D9D9"
)

text_colours <- c(
    "A/A" = "white",
    "P/A" = "#202020",
    "P/P" = "white",
    "No call" = "#777777"
)


prepare_module <- function(module_name) {

    x <- d[
        d$module == module_name,
        ,
        drop = FALSE
    ]

    site_manifest <- unique(
        x[
            order(
                x$pos
            ),
            c(
                "site_key",
                "site_label",
                "pos"
            )
        ]
    )

    sample_manifest <- unique(
        x[
            order(
                x$phenotype,
                x$taxon,
                x$data_type,
                x$sample
            ),
            c(
                "sample",
                "sample_label",
                "phenotype",
                "phenotype_label"
            )
        ]
    )

    x$site_label <- factor(
        x$site_label,
        levels = site_manifest$site_label
    )

    x$sample_label <- factor(
        x$sample_label,
        levels = rev(
            sample_manifest$sample_label
        )
    )

    x
}


make_plot <- function(
    x,
    title,
    detailed = FALSE
) {

    plot <- ggplot(
        x,
        aes(
            x = site_label,
            y = sample_label
        )
    ) +

        geom_tile(
            aes(
                fill = display_GT
            ),
            colour = "white",
            linewidth = ifelse(
                detailed,
                0.25,
                0.08
            ),
            width = 0.97,
            height = 0.96
        )

    if (detailed) {

        plot <- plot +
            geom_text(
                aes(
                    label = cell_label,
                    colour = display_GT
                ),
                size = 2.1,
                show.legend = FALSE
            )
    }

    plot +

        facet_grid(
            rows = vars(
                phenotype_label
            ),
            scales = "free_y",
            space = "free_y",
            switch = "y"
        ) +

        scale_fill_manual(
            values = genotype_colours,
            drop = FALSE,
            name = "Genotype"
        ) +

        scale_colour_manual(
            values = text_colours,
            guide = "none"
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
            title = title,
            subtitle = paste(
                "All comparative samples used in the",
                "WGS+CS balanced-accuracy calculation"
            ),
            caption = paste(
                "A/A: non-pattern homozygote;",
                "P/A: heterozygote;",
                "P/P: pattern-associated homozygote.",
                "Gray cells failed DP >= 2 or GQ >= 20."
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
                size = ifelse(
                    detailed,
                    7,
                    6
                )
            ),
            axis.text.y = element_text(
                size = ifelse(
                    detailed,
                    7,
                    4.2
                ),
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
}


write_module_plots <- function(module_name) {

    x <- prepare_module(
        module_name
    )

    trait <- unique(
        x$trait
    )

    sample_n <- length(
        unique(
            x$sample
        )
    )

    site_n <- length(
        unique(
            x$site_key
        )
    )

    width <- max(
        10,
        min(
            28,
            5 + 0.38 * site_n
        )
    )

    overview_height <- max(
        24,
        4 + 0.105 * sample_n
    )

    overview <- make_plot(
        x,
        paste0(
            module_name,
            " (",
            trait,
            "): all ",
            sample_n,
            " comparative samples"
        ),
        detailed = FALSE
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
        plot = overview,
        width = width,
        height = overview_height,
        units = "in",
        device = cairo_pdf,
        limitsize = FALSE
    )

    ggsave(
        filename = overview_svg,
        plot = overview,
        width = width,
        height = overview_height,
        units = "in",
        device = "svg",
        limitsize = FALSE
    )


    # --------------------------------------------------------
    # Paginated version: phenotype classes are kept separate
    # and every page contains at most 45 labelled samples.
    # --------------------------------------------------------

    page_size <- 45
    page_sets <- list()

    for (
        phenotype_value
        in sort(
            unique(
                x$phenotype
            )
        )
    ) {

        phenotype_samples <- unique(
            as.character(
                x$sample[
                    x$phenotype
                    == phenotype_value
                ]
            )
        )

        phenotype_samples <- phenotype_samples[
            order(
                match(
                    phenotype_samples,
                    unique(
                        x$sample
                    )
                )
            )
        ]

        chunks <- split(
            phenotype_samples,
            ceiling(
                seq_along(
                    phenotype_samples
                )
                / page_size
            )
        )

        for (chunk in chunks) {

            page_sets[[
                    length(page_sets) + 1
                ]
            ] <- list(
                phenotype = phenotype_value,
                samples = chunk
            )
        }
    }

    paginated_pdf <- file.path(
        figure_dir,
        paste0(
            "perfect_",
            module_name,
            "_all_comparative_genotypes.paginated.pdf"
        )
    )

    grDevices::cairo_pdf(
        filename = paginated_pdf,
        width = width,
        height = 13,
        onefile = TRUE
    )

    for (
        page_number
        in seq_along(
            page_sets
        )
    ) {

        page_definition <- page_sets[[page_number]]

        page_data <- x[
            x$sample
            %in% page_definition$samples,
            ,
            drop = FALSE
        ]

        sample_order <- unique(
            as.character(
                x$sample_label[
                    x$sample
                    %in% page_definition$samples
                ]
            )
        )

        page_data$sample_label <- factor(
            page_data$sample_label,
            levels = rev(
                sample_order
            )
        )

        page_plot <- make_plot(
            page_data,
            paste0(
                module_name,
                " (",
                trait,
                ") — phenotype ",
                trait,
                "=",
                page_definition$phenotype,
                " — page ",
                page_number,
                "/",
                length(page_sets)
            ),
            detailed = TRUE
        )

        print(
            page_plot
        )
    }

    dev.off()

    data.frame(
        module = module_name,
        trait = trait,
        samples = sample_n,
        perfect_SNPs = site_n,
        pages = length(page_sets),
        overview_pdf = overview_pdf,
        paginated_pdf = paginated_pdf
    )
}


summaries <- rbind(
    write_module_plots(
        "Dennis"
    ),
    write_module_plots(
        "Ray"
    )
)

summary_file <- file.path(
    checkpoint_dir,
    "perfect_D_R_all_comparative_genotype_plot_summary.tsv"
)

write.table(
    summaries,
    file = summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat(
    "PERFECT D/R ALL-COMPARATIVE GENOTYPE PLOTS\n\n"
)

print(
    summaries,
    row.names = FALSE
)

cat(
    "\nSummary:\n",
    summary_file,
    "\n",
    sep = ""
)
