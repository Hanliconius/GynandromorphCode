#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
})

project_dir <- "/scratch/martinlab/joe/GynandroCalling"

input_file <- file.path(
    project_dir,
    "analysis/checkpoints/",
    "optix_focal_pseudopool_evidence.wide.tsv"
)

plotted_sites_file <- file.path(
    project_dir,
    "analysis/checkpoints/",
    "optix_minimal_main_figure_sites.tsv"
)

figure_dir <- file.path(
    project_dir,
    "analysis/figures"
)

pdf_file <- file.path(
    figure_dir,
    "optix_focal_pseudopool_minimal_v3.pdf"
)

svg_file <- file.path(
    figure_dir,
    "optix_focal_pseudopool_minimal_v3.svg"
)

dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

if (!file.exists(input_file)) {
    message("ERROR: input file is missing:")
    message(input_file)
} else {
    focal <- read.delim(
        input_file,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE,
        quote = ""
    )

    panel_order <- c(
        "band",
        "ray1",
        "ray2",
        "dennis2"
    )

    panel_labels <- c(
        band = "Band",
        ray1 = "Ray 1",
        ray2 = "Ray 2",
        dennis2 = "Dennis 2"
    )

    focal <- focal[
        focal$panel %in% panel_order,
    ]

    focal$is_representative <- (
        tolower(as.character(
            focal$is_representative
        )) == "true"
    )

    # Use one nonredundant SNP per comparative
    # genotype-pattern group.
    plot_sites <- focal[
        focal$is_representative,
    ]

    plot_sites$panel <- factor(
        plot_sites$panel,
        levels = panel_order
    )

    plot_sites$position_Mb <- (
        as.numeric(plot_sites$pos) / 1e6
    )

    plot_sites$left_alt_fraction <- as.numeric(
        plot_sites$left_alt_fraction
    )

    plot_sites$right_alt_fraction <- as.numeric(
        plot_sites$right_alt_fraction
    )

    # Polarize toward the marker allele associated with
    # presence of the relevant wing-pattern element.
    #
    # Band presence is REF-associated:
    # pattern fraction = 1 - ALT fraction.
    #
    # Dennis and Ray presence are ALT-associated:
    # pattern fraction = ALT fraction.
    plot_sites$left_pattern_fraction <- ifelse(
        plot_sites$trait == "B",
        1 - plot_sites$left_alt_fraction,
        plot_sites$left_alt_fraction
    )

    plot_sites$right_pattern_fraction <- ifelse(
        plot_sites$trait == "B",
        1 - plot_sites$right_alt_fraction,
        plot_sites$right_alt_fraction
    )

    plot_sites <- plot_sites[
        order(
            plot_sites$panel,
            plot_sites$pos
        ),
    ]

    write.table(
        plot_sites,
        file = plotted_sites_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE,
        na = ""
    )

    left_points <- data.frame(
        panel = plot_sites$panel,
        pos = plot_sites$pos,
        position_Mb = plot_sites$position_Mb,
        side = "Left",
        pattern_fraction = (
            plot_sites$left_pattern_fraction
        ),
        stringsAsFactors = FALSE
    )

    right_points <- data.frame(
        panel = plot_sites$panel,
        pos = plot_sites$pos,
        position_Mb = plot_sites$position_Mb,
        side = "Right",
        pattern_fraction = (
            plot_sites$right_pattern_fraction
        ),
        stringsAsFactors = FALSE
    )

    point_data <- rbind(
        left_points,
        right_points
    )

    point_data$panel <- factor(
        point_data$panel,
        levels = panel_order
    )

    point_data$side <- factor(
        point_data$side,
        levels = c("Left", "Right")
    )

    theme_minimal_gynandro <- theme_classic(
        base_size = 10,
        base_family = "sans"
    ) +
        theme(
            plot.title = element_text(
                face = "bold",
                size = 13,
                hjust = 0
            ),
            strip.background = element_rect(
                fill = "grey94",
                colour = NA
            ),
            strip.text = element_text(
                face = "bold",
                size = 10
            ),
            axis.title = element_text(size = 10),
            axis.text = element_text(
                size = 8,
                colour = "black"
            ),
            legend.title = element_blank(),
            legend.text = element_text(size = 9),
            legend.position = "bottom",
            panel.spacing = grid::unit(
                0.8,
                "lines"
            ),
            plot.caption = element_text(
                size = 8,
                colour = "grey30",
                hjust = 0
            ),
            plot.margin = margin(
                8, 10, 8, 8
            )
        )

    p <- ggplot() +
        # A single subtle balanced/heterozygote-like zone.
        geom_rect(
            data = plot_sites,
            aes(
                xmin = -Inf,
                xmax = Inf,
                ymin = 0.40,
                ymax = 0.60
            ),
            inherit.aes = FALSE,
            fill = "#DCEBF7",
            alpha = 0.35
        ) +
        geom_hline(
            yintercept = c(0, 0.5, 1),
            linewidth = 0.28,
            colour = "grey72"
        ) +
        # Each line and both endpoints are derived only
        # from the four focal leg samples.
        geom_segment(
            data = plot_sites,
            aes(
                x = position_Mb,
                xend = position_Mb,
                y = left_pattern_fraction,
                yend = right_pattern_fraction
            ),
            inherit.aes = FALSE,
            colour = "grey60",
            linewidth = 0.45
        ) +
        geom_point(
            data = point_data,
            aes(
                x = position_Mb,
                y = pattern_fraction,
                colour = side
            ),
            size = 3.0,
            shape = 16
        ) +
        # Comparative information is represented only
        # by a downward arrowhead at the SNP position.
        geom_point(
            data = plot_sites,
            aes(
                x = position_Mb,
                y = 1.065
            ),
            inherit.aes = FALSE,
            shape = 25,
            size = 2.7,
            stroke = 0.25,
            colour = "black",
            fill = "black"
        ) +
        facet_wrap(
            ~panel,
            scales = "free_x",
            ncol = 1,
            labeller = as_labeller(panel_labels)
        ) +
        scale_colour_manual(
            values = c(
                Left = "#D55E00",
                Right = "#0072B2"
            )
        ) +
        scale_y_continuous(
            limits = c(-0.03, 1.10),
            breaks = c(0, 0.5, 1),
            labels = c(
                "0",
                "0.5",
                "1"
            ),
            expand = expansion(
                mult = c(0, 0)
            )
        ) +
        scale_x_continuous(
            expand = expansion(
                mult = c(0.03, 0.03)
            )
        ) +
        labs(
            title = "Left-right allele balance at informative Optix markers",
            x = "Hmr18 position (Mb)",
            y = "Pattern-associated marker-allele fraction",
            colour = NULL,
            caption = paste(
                "Downward arrowheads mark nonredundant SNPs",
                "identified by the comparative genotype-phenotype",
                "analysis; points show only the temporary",
                "left and right focal pseudo-pools."
            )
        ) +
        theme_minimal_gynandro

    ggsave(
        pdf_file,
        p,
        width = 10,
        height = 9,
        units = "in",
        device = cairo_pdf
    )

    ggsave(
        svg_file,
        p,
        width = 10,
        height = 9,
        units = "in",
        device = "svg"
    )

    message("MINIMAL POLARIZED OPTIX FIGURE COMPLETED")
    message("PDF: ", pdf_file)
    message("SVG: ", svg_file)
    message("Plotted-site table: ", plotted_sites_file)
    message("Plotted sites: ", nrow(plot_sites))
}
