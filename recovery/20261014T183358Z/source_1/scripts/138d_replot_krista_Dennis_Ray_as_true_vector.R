#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(grid)
})

project <- "/scratch/martinlab/joe/GynandroCalling"

input_dir <- file.path(
    project,
    "analysis/original_502/checkpoints",
    "krista_DR_regions_reciprocal_genoplots"
)

figure_dir <- file.path(
    project,
    "analysis/figures",
    "krista_DR_regions_reciprocal_genoplots"
)

summary_file <- file.path(
    input_dir,
    "krista_DR_regions_reciprocal_genoplots.vector_plot_summary.tsv"
)

dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

panels <- data.frame(
    panel = c(
        "Dennis1",
        "Dennis2",
        "Ray1",
        "Ray2"
    ),
    trait = c(
        "D",
        "D",
        "R",
        "R"
    ),
    expected_snps = c(
        204L,
        638L,
        850L,
        2459L
    ),
    stringsAsFactors = FALSE
)

genotype_colours <- c(
    "A/A" = "#2C9BCB",
    "P/A" = "#F28E2B",
    "P/P" = "#B61F72",
    "N" = "#F2F2F2",
    "Uncertain" = "#777777"
)

genotype_levels <- names(
    genotype_colours
)


make_vector_runs <- function(data) {

    data <- data[
        order(
            data$display_order,
            data$site_index
        ),
        ,
        drop = FALSE
    ]

    split_key <- interaction(
        as.character(
            data$plot_group
        ),
        as.character(
            data$sample_label
        ),
        drop = TRUE,
        lex.order = TRUE
    )

    sample_pieces <- split(
        data,
        split_key
    )

    run_tables <- lapply(
        sample_pieces,
        function(sample_piece) {

            sample_piece <- sample_piece[
                order(
                    sample_piece$site_index
                ),
                ,
                drop = FALSE
            ]

            run_encoding <- rle(
                as.character(
                    sample_piece$display_GT
                )
            )

            run_end <- cumsum(
                run_encoding$lengths
            )

            run_start <- (
                run_end
                -
                run_encoding$lengths
                +
                1L
            )

            data.frame(
                plot_group = as.character(
                    sample_piece$plot_group[1]
                ),
                sample_label = as.character(
                    sample_piece$sample_label[1]
                ),
                display_GT = run_encoding$values,
                run_start = run_start,
                run_end = run_end,
                run_midpoint = (
                    run_start
                    +
                    run_end
                ) / 2,
                run_width = (
                    run_end
                    -
                    run_start
                    +
                    1
                ),
                stringsAsFactors = FALSE
            )
        }
    )

    do.call(
        rbind,
        run_tables
    )
}


summary_rows <- list()


for (
    panel_index
    in seq_len(
        nrow(
            panels
        )
    )
) {

    panel_name <- panels$panel[
        panel_index
    ]

    trait <- panels$trait[
        panel_index
    ]

    expected_snps <- panels$expected_snps[
        panel_index
    ]

    input_file <- file.path(
        input_dir,
        paste0(
            "optix_",
            panel_name,
            "_original_region_reciprocal_genoplot.long.tsv.gz"
        )
    )

    pdf_file <- file.path(
        figure_dir,
        paste0(
            "optix_",
            panel_name,
            "_original_region_genotypes_with_",
            "reciprocal_and_piano_keys.pdf"
        )
    )

    if (
        !file.exists(
            input_file
        )
    ) {

        stop(
            "Missing plot input: ",
            input_file
        )
    }


    genotype_data <- read.delim(
        gzfile(
            input_file
        ),
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE,
        quote = "",
        comment.char = ""
    )

    genotype_data$pos <- as.integer(
        genotype_data$pos
    )

    genotype_data$display_order <- as.integer(
        genotype_data$display_order
    )

    genotype_data$changed_from_original <- (
        tolower(
            as.character(
                genotype_data$changed_from_original
            )
        )
        %in%
        c(
            "true",
            "t",
            "1",
            "yes"
        )
    )


    sites <- unique(
        genotype_data[
            order(
                genotype_data$pos,
                genotype_data$ref,
                genotype_data$alt
            ),
            c(
                "chrom",
                "pos",
                "ref",
                "alt",
                "site_label"
            )
        ]
    )

    if (
        nrow(
            sites
        ) != expected_snps
    ) {

        stop(
            panel_name,
            ": expected ",
            expected_snps,
            " SNPs; recovered ",
            nrow(
                sites
            )
        )
    }


    sites$site_key <- paste(
        sites$chrom,
        sites$pos,
        sites$ref,
        sites$alt,
        sep = ":"
    )

    genotype_data$site_key <- paste(
        genotype_data$chrom,
        genotype_data$pos,
        genotype_data$ref,
        genotype_data$alt,
        sep = ":"
    )

    genotype_data$site_index <- match(
        genotype_data$site_key,
        sites$site_key
    )


    sample_order <- unique(
        genotype_data[
            order(
                genotype_data$display_order
            ),
            c(
                "display_order",
                "sample_label",
                "plot_group"
            )
        ]
    )

    group_order <- unique(
        sample_order$plot_group
    )

    sample_levels <- rev(
        sample_order$sample_label
    )


    genotype_data$display_GT[
        !genotype_data$display_GT
        %in%
        genotype_levels
    ] <- "N"

    genotype_data$plot_group <- factor(
        genotype_data$plot_group,
        levels = group_order,
        ordered = TRUE
    )

    genotype_data$sample_label <- factor(
        genotype_data$sample_label,
        levels = sample_levels,
        ordered = TRUE
    )

    genotype_data$display_GT <- factor(
        genotype_data$display_GT,
        levels = genotype_levels,
        ordered = TRUE
    )


    # Collapse only directly adjacent cells having the same genotype.
    # Every SNP column remains represented by its exact original width.
    vector_runs <- make_vector_runs(
        genotype_data
    )

    vector_runs$plot_group <- factor(
        vector_runs$plot_group,
        levels = group_order,
        ordered = TRUE
    )

    vector_runs$sample_label <- factor(
        vector_runs$sample_label,
        levels = sample_levels,
        ordered = TRUE
    )

    vector_runs$display_GT <- factor(
        vector_runs$display_GT,
        levels = genotype_levels,
        ordered = TRUE
    )


    changed_calls <- genotype_data[
        (
            genotype_data$call_source
            ==
            "reciprocal mixture (posterior >= 0.95)"
        )
        &
        genotype_data$changed_from_original,
        ,
        drop = FALSE
    ]


    tick_count <- min(
        26L,
        expected_snps
    )

    tick_index <- unique(
        round(
            seq(
                1,
                expected_snps,
                length.out = tick_count
            )
        )
    )

    tick_labels <- format(
        sites$pos[
            tick_index
        ],
        big.mark = ",",
        scientific = FALSE,
        trim = TRUE
    )

    plot_height <- max(
        22,
        4
        +
        0.105
        *
        nrow(
            sample_order
        )
    )


    genotype_plot <- ggplot() +

        geom_tile(
            data = vector_runs,
            aes(
                x = run_midpoint,
                y = sample_label,
                width = run_width,
                fill = display_GT
            ),
            height = 0.98,
            colour = NA
        ) +

        geom_tile(
            data = changed_calls,
            aes(
                x = site_index,
                y = sample_label
            ),
            width = 1,
            height = 0.98,
            fill = NA,
            colour = "black",
            linewidth = 0.24
        ) +

        facet_grid(
            rows = vars(
                plot_group
            ),
            scales = "free_y",
            space = "free_y",
            switch = "y"
        ) +

        scale_fill_manual(
            values = genotype_colours,
            limits = genotype_levels,
            drop = FALSE,
            name = "Pattern-associated genotype"
        ) +

        scale_x_continuous(
            breaks = tick_index,
            labels = tick_labels,
            limits = c(
                0.5,
                expected_snps + 0.5
            ),
            expand = expansion(
                add = c(
                    0,
                    0
                )
            )
        ) +

        scale_y_discrete(
            expand = expansion(
                add = c(
                    0.05,
                    0.05
                )
            )
        ) +

        labs(
            title = paste(
                "Genotypes across Krista's original",
                panel_name,
                "interval"
            ),
            subtitle = paste(
                "All",
                expected_snps,
                "biallelic SNPs are shown; no BA threshold was used.",
                "The original 481-sample",
                trait,
                "cohort is retained."
            ),
            x = "Genomic position on Hmr18 (bp)",
            y = NULL,
            caption = paste(
                "Piano-key and pooled focal samples were not used",
                "to select SNPs.",
                "TN_Left and TN_Right are shown as original and",
                "reciprocal-mixture calls.",
                "Black outlines mark posterior >= 0.95 calls that",
                "differ from the original genotype."
            )
        ) +

        guides(
            fill = guide_legend(
                nrow = 1,
                byrow = TRUE,
                override.aes = list(
                    colour = "#555555",
                    linewidth = 0.25
                )
            )
        ) +

        theme_minimal(
            base_size = 9
        ) +

        theme(
            plot.title = element_text(
                face = "bold",
                size = 14
            ),
            plot.subtitle = element_text(
                size = 9.5,
                colour = "#444444",
                margin = margin(
                    b = 8
                )
            ),
            plot.caption = element_text(
                size = 8,
                colour = "#444444",
                hjust = 0,
                margin = margin(
                    t = 8
                )
            ),
            axis.text.x = element_text(
                angle = 90,
                hjust = 1,
                vjust = 0.5,
                size = 5
            ),
            axis.text.y = element_text(
                size = 4.5,
                colour = "#222222"
            ),
            axis.title.x = element_text(
                size = 10,
                margin = margin(
                    t = 8
                )
            ),
            panel.grid = element_blank(),
            panel.spacing.y = unit(
                0.10,
                "lines"
            ),
            strip.placement = "outside",
            strip.background = element_rect(
                fill = "#E6E6E6",
                colour = NA
            ),
            strip.text.y.left = element_text(
                angle = 0,
                face = "bold",
                size = 7
            ),
            legend.position = "top",
            legend.justification = "left",
            legend.title = element_text(
                face = "bold"
            ),
            plot.margin = margin(
                12,
                16,
                12,
                12
            )
        )


    ggsave(
        filename = pdf_file,
        plot = genotype_plot,
        device = cairo_pdf,
        width = 20,
        height = plot_height,
        units = "in",
        limitsize = FALSE
    )


    summary_rows[[
        panel_index
    ]] <- data.frame(
        panel = panel_name,
        trait = trait,
        SNPs = nrow(
            sites
        ),
        display_rows = nrow(
            sample_order
        ),
        original_cells = nrow(
            genotype_data
        ),
        vector_rectangles = nrow(
            vector_runs
        ),
        changed_reciprocal_cells = nrow(
            changed_calls
        ),
        PDF_height_inches = plot_height,
        PDF = pdf_file,
        stringsAsFactors = FALSE
    )


    cat(
        panel_name,
        ": ",
        format(
            nrow(
                genotype_data
            ),
            big.mark = ","
        ),
        " cells compressed to ",
        format(
            nrow(
                vector_runs
            ),
            big.mark = ","
        ),
        " vector rectangles\n",
        sep = ""
    )


    rm(
        genotype_data,
        vector_runs,
        genotype_plot
    )

    invisible(
        gc()
    )
}


summary <- do.call(
    rbind,
    summary_rows
)

write.table(
    summary,
    summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat(
    "\nTRUE-VECTOR KRISTA DENNIS/RAY GENOPLOTS\n\n"
)

print(
    summary[
        ,
        setdiff(
            names(
                summary
            ),
            "PDF"
        )
    ],
    row.names = FALSE
)

cat(
    "\nOverwritten PDFs:\n",
    paste(
        summary$PDF,
        collapse = "\n"
    ),
    "\n",
    sep = ""
)

cat(
    "\nSummary:\n",
    summary_file,
    "\n",
    sep = ""
)
