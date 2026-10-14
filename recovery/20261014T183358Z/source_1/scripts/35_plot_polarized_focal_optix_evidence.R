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

label_table_file <- file.path(
    project_dir,
    "analysis/checkpoints/",
    "optix_focal_pseudopool_labeled_sites.tsv"
)

figure_dir <- file.path(
    project_dir,
    "analysis/figures"
)

pdf_file <- file.path(
    figure_dir,
    "optix_focal_pseudopool_polarized_v2.pdf"
)

svg_file <- file.path(
    figure_dir,
    "optix_focal_pseudopool_polarized_v2.svg"
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

    focal$panel <- factor(
        focal$panel,
        levels = panel_order
    )

    focal$position_Mb <- focal$pos / 1e6

    focal$minimum_BA <- as.numeric(
        focal$minimum_BA
    )

    focal$minimum_MCC <- as.numeric(
        focal$minimum_MCC
    )

    focal$left_alt_fraction <- as.numeric(
        focal$left_alt_fraction
    )

    focal$right_alt_fraction <- as.numeric(
        focal$right_alt_fraction
    )

    focal$is_representative <- as.logical(
        focal$is_representative
    )

    focal$exact_left_11_right_01_match <- as.logical(
        focal$exact_left_11_right_01_match
    )

    focal$AB_match_threshold_0.80 <- as.logical(
        focal$AB_match_threshold_0.80
    )

    focal$tier_label <- ifelse(
        focal$tier == "tier_1",
        "Tier 1",
        "Tier 2"
    )

    focal$tier_label <- factor(
        focal$tier_label,
        levels = c("Tier 1", "Tier 2")
    )

    # Polarize every site toward the marker allele associated
    # with presence of its relevant wing-pattern element.
    focal$left_pattern_fraction <- ifelse(
        focal$trait == "B",
        1 - focal$left_alt_fraction,
        focal$left_alt_fraction
    )

    focal$right_pattern_fraction <- ifelse(
        focal$trait == "B",
        1 - focal$right_alt_fraction,
        focal$right_alt_fraction
    )

    focal$focal_support_class <- ifelse(
        focal$exact_left_11_right_01_match,
        "Formal pseudo-GT support",
        ifelse(
            focal$AB_match_threshold_0.80,
            "Allele-balance support",
            "No focal-asymmetry support"
        )
    )

    focal$focal_support_class <- factor(
        focal$focal_support_class,
        levels = c(
            "Formal pseudo-GT support",
            "Allele-balance support",
            "No focal-asymmetry support"
        )
    )

    # Select no more than one labelled SNP per comparative
    # genotype-pattern group. Focal evidence is prioritized,
    # followed by comparative tier and association strength.
    focal$label_eligible <- (
        focal$exact_left_11_right_01_match |
        focal$AB_match_threshold_0.80
    )

    focal$formal_rank <- ifelse(
        focal$exact_left_11_right_01_match,
        1,
        2
    )

    focal$tier_rank <- ifelse(
        focal$tier == "tier_1",
        1,
        2
    )

    focal_ordered <- focal[
        order(
            focal$group_id,
            !focal$label_eligible,
            focal$formal_rank,
            focal$tier_rank,
            -focal$minimum_BA,
            -focal$minimum_MCC,
            focal$pos
        ),
    ]

    label_sites <- focal_ordered[
        !duplicated(focal_ordered$group_id) &
        focal_ordered$label_eligible,
    ]

    label_sites <- label_sites[
        order(
            label_sites$panel,
            label_sites$pos
        ),
    ]

    # Alternate label heights within each panel.
    label_sites$label_index <- ave(
        label_sites$pos,
        label_sites$panel,
        FUN = seq_along
    )

    label_sites$label_y <- ifelse(
        label_sites$label_index %% 2 == 1,
        1.045,
        1.115
    )

    label_sites$position_label <- format(
        label_sites$pos,
        scientific = FALSE,
        trim = TRUE
    )

    write.table(
        label_sites,
        file = label_table_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE,
        na = ""
    )

    # Long table for plotting the two sides.
    left_points <- data.frame(
        panel = focal$panel,
        panel_text = as.character(focal$panel),
        trait = focal$trait,
        pos = focal$pos,
        position_Mb = focal$position_Mb,
        side = "Left",
        pattern_fraction = focal$left_pattern_fraction,
        tier_label = focal$tier_label,
        is_representative = focal$is_representative,
        focal_support_class = focal$focal_support_class,
        stringsAsFactors = FALSE
    )

    right_points <- data.frame(
        panel = focal$panel,
        panel_text = as.character(focal$panel),
        trait = focal$trait,
        pos = focal$pos,
        position_Mb = focal$position_Mb,
        side = "Right",
        pattern_fraction = focal$right_pattern_fraction,
        tier_label = focal$tier_label,
        is_representative = focal$is_representative,
        focal_support_class = focal$focal_support_class,
        stringsAsFactors = FALSE
    )

    focal_long <- rbind(
        left_points,
        right_points
    )

    focal_long$panel <- factor(
        focal_long$panel,
        levels = panel_order
    )

    focal_long$side <- factor(
        focal_long$side,
        levels = c("Left", "Right")
    )

    # Calculate a right-side annotation position separately
    # for each genomic interval.
    panel_ranges <- aggregate(
        position_Mb ~ panel,
        data = focal,
        FUN = function(values) {
            c(
                minimum = min(values),
                maximum = max(values)
            )
        }
    )

    range_rows <- do.call(
        rbind,
        lapply(
            seq_len(nrow(panel_ranges)),
            function(index) {
                values <- panel_ranges$position_Mb[
                    index,
                ][[1]]

                interval_width <- (
                    values["maximum"] -
                    values["minimum"]
                )

                if (
                    !is.finite(interval_width) ||
                    interval_width <= 0
                ) {
                    interval_width <- 0.0001
                }

                data.frame(
                    panel = panel_ranges$panel[index],
                    minimum = values["minimum"],
                    maximum = values["maximum"],
                    width = interval_width,
                    annotation_x = (
                        values["maximum"] +
                        interval_width * 0.12
                    ),
                    text_x = (
                        values["maximum"] +
                        interval_width * 0.17
                    )
                )
            }
        )
    )

    range_rows$panel <- factor(
        range_rows$panel,
        levels = panel_order
    )

    expected_rows <- do.call(
        rbind,
        lapply(
            seq_len(nrow(range_rows)),
            function(index) {
                panel_name <- as.character(
                    range_rows$panel[index]
                )

                expected_left <- ifelse(
                    panel_name == "band",
                    0,
                    1
                )

                data.frame(
                    panel = range_rows$panel[index],
                    annotation_x = range_rows$annotation_x[index],
                    side = c("Left", "Right"),
                    expected_fraction = c(
                        expected_left,
                        0.5
                    ),
                    stringsAsFactors = FALSE
                )
            }
        )
    )

    expected_rows$panel <- factor(
        expected_rows$panel,
        levels = panel_order
    )

    expected_rows$side <- factor(
        expected_rows$side,
        levels = c("Left", "Right")
    )

    # Panel-specific phenotype interpretation of 0, 0.5 and 1.
    phenotype_rows <- rbind(
        data.frame(
            panel = "band",
            state_y = c(0, 0.5, 1),
            state_label = c(
                "no Band\nnon-B/non-B",
                "Band\nB/non-B",
                "Band\nB/B"
            )
        ),
        data.frame(
            panel = "ray1",
            state_y = c(0, 0.5, 1),
            state_label = c(
                "no Rays\nnon-R/non-R",
                "Rays\nR/non-R",
                "Rays\nR/R"
            )
        ),
        data.frame(
            panel = "ray2",
            state_y = c(0, 0.5, 1),
            state_label = c(
                "no Rays\nnon-R/non-R",
                "Rays\nR/non-R",
                "Rays\nR/R"
            )
        ),
        data.frame(
            panel = "dennis2",
            state_y = c(0, 0.5, 1),
            state_label = c(
                "no Dennis\nnon-D/non-D",
                "Dennis\nD/non-D",
                "Dennis\nD/D"
            )
        )
    )

    phenotype_rows$panel <- factor(
        phenotype_rows$panel,
        levels = panel_order
    )

    phenotype_rows <- merge(
        phenotype_rows,
        range_rows[
            ,
            c("panel", "text_x")
        ],
        by = "panel",
        all.x = TRUE,
        sort = FALSE
    )

    phenotype_rows$panel <- factor(
        phenotype_rows$panel,
        levels = panel_order
    )

    # Plot-size and opacity vectors.
    focal_long$point_size <- ifelse(
        focal_long$is_representative,
        3.0,
        1.6
    )

    focal_long$point_alpha <- ifelse(
        focal_long$is_representative,
        1.0,
        0.50
    )

    theme_gynandro <- theme_classic(
        base_size = 10,
        base_family = "sans"
    ) +
        theme(
            plot.title = element_text(
                face = "bold",
                size = 13,
                hjust = 0
            ),
            plot.subtitle = element_text(
                size = 9,
                colour = "grey25",
                hjust = 0,
                margin = margin(b = 6)
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
            legend.title = element_text(
                face = "bold",
                size = 9
            ),
            legend.text = element_text(size = 8),
            legend.position = "bottom",
            panel.spacing = grid::unit(
                0.8,
                "lines"
            ),
            plot.margin = margin(
                10, 16, 8, 8
            )
        )

    p <- ggplot() +
        # Heterozygote-like/balanced zone.
        geom_rect(
            data = range_rows,
            aes(
                xmin = -Inf,
                xmax = Inf,
                ymin = 0.30,
                ymax = 0.70
            ),
            inherit.aes = FALSE,
            fill = "#DCEBF7",
            alpha = 0.32
        ) +
        # Pattern-associated homozygote-like zone.
        geom_rect(
            data = range_rows,
            aes(
                xmin = -Inf,
                xmax = Inf,
                ymin = 0.80,
                ymax = 1.00
            ),
            inherit.aes = FALSE,
            fill = "#F6E8A6",
            alpha = 0.32
        ) +
        geom_hline(
            yintercept = c(0, 0.5, 1),
            linewidth = 0.25,
            colour = "grey72"
        ) +
        # Join left and right observations at each SNP.
        geom_segment(
            data = focal,
            aes(
                x = position_Mb,
                xend = position_Mb,
                y = left_pattern_fraction,
                yend = right_pattern_fraction
            ),
            colour = "grey62",
            linewidth = 0.38,
            alpha = 0.60
        ) +
        geom_point(
            data = focal_long,
            aes(
                x = position_Mb,
                y = pattern_fraction,
                colour = side,
                shape = tier_label
            ),
            size = focal_long$point_size,
            alpha = focal_long$point_alpha,
            stroke = 0.6
        ) +
        # Ring sites supporting the focal haplotype model.
        geom_point(
            data = focal_long[
                focal_long$focal_support_class !=
                "No focal-asymmetry support",
            ],
            aes(
                x = position_Mb,
                y = pattern_fraction
            ),
            shape = 21,
            fill = NA,
            colour = "black",
            size = 4.1,
            stroke = 0.7
        ) +
        # Position labels: one focal-supporting SNP per
        # comparative genotype-pattern group.
        geom_segment(
            data = label_sites,
            aes(
                x = position_Mb,
                xend = position_Mb,
                y = pmax(
                    left_pattern_fraction,
                    right_pattern_fraction
                ) + 0.02,
                yend = label_y - 0.01
            ),
            inherit.aes = FALSE,
            linewidth = 0.3,
            colour = "grey35"
        ) +
        geom_text(
            data = label_sites,
            aes(
                x = position_Mb,
                y = label_y,
                label = position_label
            ),
            inherit.aes = FALSE,
            angle = 90,
            hjust = 0,
            vjust = 0.5,
            size = 2.35,
            family = "sans"
        ) +
        # Expected left/right positions under the linked
        # B-versus-DR haplotype model.
        geom_segment(
            data = range_rows,
            aes(
                x = annotation_x,
                xend = annotation_x,
                y = ifelse(
                    as.character(panel) == "band",
                    0,
                    1
                ),
                yend = 0.5
            ),
            inherit.aes = FALSE,
            linewidth = 0.55,
            colour = "grey30"
        ) +
        geom_point(
            data = expected_rows,
            aes(
                x = annotation_x,
                y = expected_fraction,
                colour = side
            ),
            inherit.aes = FALSE,
            size = 3.1,
            shape = 18
        ) +
        geom_text(
            data = expected_rows,
            aes(
                x = annotation_x,
                y = expected_fraction,
                label = substr(side, 1, 1),
                colour = side
            ),
            inherit.aes = FALSE,
            nudge_x = 0,
            nudge_y = 0.055,
            size = 2.6,
            fontface = "bold"
        ) +
        # Text placeholders for eventual phenotype cartoons.
        geom_label(
            data = phenotype_rows,
            aes(
                x = text_x,
                y = state_y,
                label = state_label
            ),
            inherit.aes = FALSE,
            hjust = 0,
            size = 2.15,
            lineheight = 0.9,
            label.size = 0,
            fill = scales::alpha("white", 0.70),
            colour = "grey20"
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
            ),
            name = "Pseudo-pooled side"
        ) +
        scale_shape_manual(
            values = c(
                "Tier 1" = 16,
                "Tier 2" = 17
            ),
            name = "Comparative association"
        ) +
        scale_x_continuous(
            expand = expansion(
                mult = c(0.02, 0.30)
            )
        ) +
        scale_y_continuous(
            breaks = c(0, 0.5, 1),
            labels = c("0", "0.5", "1"),
            limits = c(-0.04, 1.18),
            expand = expansion(mult = c(0, 0))
        ) +
        labs(
            title = "Comparative Optix markers reveal focal left-right haplotype mixing",
            subtitle = paste(
                "Values are polarized toward the marker allele",
                "associated with presence of Band, Rays or Dennis.",
                "Black rings mark sites supporting an ALT-enriched",
                "left and balanced right pattern before polarization."
            ),
            x = "Hmr18 position (Mb)",
            y = "Pattern-associated marker-allele fraction"
        ) +
        theme_gynandro

    ggsave(
        pdf_file,
        p,
        width = 12,
        height = 10.5,
        units = "in",
        device = cairo_pdf
    )

    ggsave(
        svg_file,
        p,
        width = 12,
        height = 10.5,
        units = "in",
        device = "svg"
    )

    message("POLARIZED FOCAL OPTIX FIGURE COMPLETED")
    message("PDF: ", pdf_file)
    message("SVG: ", svg_file)
    message("Labelled-site table: ", label_table_file)
    message("Labelled SNPs: ", nrow(label_sites))
}
