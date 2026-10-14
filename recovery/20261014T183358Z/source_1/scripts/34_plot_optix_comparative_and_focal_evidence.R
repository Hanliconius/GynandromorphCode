#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
})

project_dir <- "/scratch/martinlab/joe/GynandroCalling"

score_file <- file.path(
    project_dir,
    "analysis/checkpoints/optix_51sample_comparative_scores.tsv"
)

pseudopool_file <- file.path(
    project_dir,
    "analysis/checkpoints/optix_focal_pseudopool_evidence.wide.tsv"
)

selection_file <- file.path(
    project_dir,
    "analysis/checkpoints/optix_joint_evidence_figure_sites.tsv"
)

figure_dir <- file.path(
    project_dir,
    "analysis/figures"
)

combined_pdf <- file.path(
    figure_dir,
    "optix_comparative_and_focal_evidence_preliminary.pdf"
)

association_svg <- file.path(
    figure_dir,
    "optix_comparative_association_landscape.svg"
)

balance_svg <- file.path(
    figure_dir,
    "optix_focal_pseudopool_allele_balance.svg"
)

selected_svg <- file.path(
    figure_dir,
    "optix_joint_evidence_selected_sites.svg"
)

dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

required_files <- c(
    score_file,
    pseudopool_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
    message("ERROR: required files are missing:")
    message(paste(missing_files, collapse = "\n"))
} else {
    # Repair the one malformed delimiter in the existing score-table
    # header without modifying the source file.
    score_lines <- readLines(
        score_file,
        warn = FALSE
    )

    score_lines[1] <- sub(
        "phenotype_0_totaphenotype_1_total",
        "phenotype_0_total\tphenotype_1_total",
        score_lines[1],
        fixed = TRUE
    )

    score_connection <- textConnection(score_lines)

    scores <- read.delim(
        score_connection,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE,
        quote = ""
    )

    close(score_connection)

    focal <- read.delim(
        pseudopool_file,
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
        "dennis1",
        "dennis2"
    )

    panel_labels <- c(
        band = "Band",
        ray1 = "Ray 1",
        ray2 = "Ray 2",
        dennis1 = "Dennis 1",
        dennis2 = "Dennis 2"
    )

    scores$panel <- factor(
        scores$panel,
        levels = panel_order
    )

    focal$panel <- factor(
        focal$panel,
        levels = panel_order
    )

    scores$scheme_label <- ifelse(
        scores$scheme == "DP2_GQ20",
        "DP >= 2; GQ >= 20",
        "Called GT"
    )

    scores$scheme_label <- factor(
        scores$scheme_label,
        levels = c(
            "Called GT",
            "DP >= 2; GQ >= 20"
        )
    )

    scores$position_Mb <- scores$pos / 1e6
    focal$position_Mb <- focal$pos / 1e6

    scores$score_status <- ifelse(
        scores$score_eligible,
        "Score eligible",
        ifelse(
            scores$callability_eligible,
            "Callable but genotype-class limited",
            "Callability limited"
        )
    )

    scores$score_status <- factor(
        scores$score_status,
        levels = c(
            "Score eligible",
            "Callable but genotype-class limited",
            "Callability limited"
        )
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

    focal$minimum_BA <- as.numeric(focal$minimum_BA)
    focal$minimum_MCC <- as.numeric(focal$minimum_MCC)
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

    focal$joint_evidence_class <- ifelse(
        focal$exact_left_11_right_01_match,
        "Formal 1/1 left; 0/1 right",
        ifelse(
            focal$AB_match_threshold_0.80,
            "ALT-enriched left; balanced right",
            "No focal asymmetry match"
        )
    )

    focal$joint_evidence_class <- factor(
        focal$joint_evidence_class,
        levels = c(
            "Formal 1/1 left; 0/1 right",
            "ALT-enriched left; balanced right",
            "No focal asymmetry match"
        )
    )

    theme_gynandro <- theme_classic(base_size = 10) +
        theme(
            plot.title = element_text(
                face = "bold",
                size = 12,
                hjust = 0
            ),
            plot.subtitle = element_text(
                size = 9,
                colour = "grey25",
                hjust = 0
            ),
            strip.background = element_rect(
                fill = "grey94",
                colour = NA
            ),
            strip.text = element_text(
                face = "bold",
                size = 9
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
                0.7,
                "lines"
            ),
            plot.margin = margin(
                8, 10, 8, 8
            )
        )

    # Page 1: complete comparative association landscape.
    p_association <- ggplot(
        scores[
            !is.na(scores$balanced_accuracy),
        ],
        aes(
            x = position_Mb,
            y = balanced_accuracy,
            colour = scheme_label
        )
    ) +
        geom_hline(
            yintercept = 0.5,
            linewidth = 0.35,
            linetype = "dashed",
            colour = "grey55"
        ) +
        geom_point(
            aes(
                alpha = score_status,
                size = score_status
            )
        ) +
        geom_vline(
            data = focal[focal$is_representative, ],
            aes(
                xintercept = position_Mb,
                linetype = tier_label
            ),
            inherit.aes = FALSE,
            linewidth = 0.35,
            colour = "black",
            alpha = 0.45
        ) +
        facet_wrap(
            ~panel,
            scales = "free_x",
            ncol = 1,
            labeller = as_labeller(panel_labels)
        ) +
        scale_colour_manual(
            values = c(
                "Called GT" = "#7A7A7A",
                "DP >= 2; GQ >= 20" = "#0072B2"
            ),
            name = "Genotype treatment"
        ) +
        scale_alpha_manual(
            values = c(
                "Score eligible" = 0.85,
                "Callable but genotype-class limited" = 0.35,
                "Callability limited" = 0.12
            ),
            name = "Site status"
        ) +
        scale_size_manual(
            values = c(
                "Score eligible" = 1.8,
                "Callable but genotype-class limited" = 1.2,
                "Callability limited" = 0.8
            ),
            name = "Site status"
        ) +
        scale_linetype_manual(
            values = c(
                "Tier 1" = "solid",
                "Tier 2" = "dotted"
            ),
            name = "Stable candidate"
        ) +
        coord_cartesian(
            ylim = c(0.45, 1.01)
        ) +
        labs(
            title = "Comparative Optix genotype-phenotype correspondence",
            subtitle = paste(
                "Balanced accuracy among 38 phenotype-scored",
                "non-focal butterflies; focal samples excluded"
            ),
            x = "Hmr18 position (Mb)",
            y = "Balanced accuracy"
        ) +
        theme_gynandro

    # Page 2: all stable candidates and focal allele balance.
    p_balance <- ggplot(
        focal,
        aes(x = position_Mb)
    ) +
        annotate(
            "rect",
            xmin = -Inf,
            xmax = Inf,
            ymin = BALANCED_LOWER <- 0.30,
            ymax = BALANCED_UPPER <- 0.70,
            fill = "#DDEBF7",
            alpha = 0.32
        ) +
        annotate(
            "rect",
            xmin = -Inf,
            xmax = Inf,
            ymin = 0.80,
            ymax = Inf,
            fill = "#F6E8A6",
            alpha = 0.30
        ) +
        geom_segment(
            aes(
                xend = position_Mb,
                y = left_alt_fraction,
                yend = right_alt_fraction
            ),
            colour = "grey65",
            linewidth = 0.35,
            alpha = 0.55
        ) +
        geom_point(
            aes(
                y = left_alt_fraction,
                shape = tier_label
            ),
            colour = "#D55E00",
            fill = "#D55E00",
            size = ifelse(
                focal$is_representative,
                3.0,
                1.6
            ),
            alpha = ifelse(
                focal$is_representative,
                1.0,
                0.55
            ),
            stroke = 0.5
        ) +
        geom_point(
            aes(
                y = right_alt_fraction,
                shape = tier_label
            ),
            colour = "#0072B2",
            fill = "#0072B2",
            size = ifelse(
                focal$is_representative,
                3.0,
                1.6
            ),
            alpha = ifelse(
                focal$is_representative,
                1.0,
                0.55
            ),
            stroke = 0.5
        ) +
        geom_point(
            data = focal[
                focal$exact_left_11_right_01_match |
                focal$AB_match_threshold_0.80,
            ],
            aes(
                y = left_alt_fraction
            ),
            shape = 21,
            fill = NA,
            colour = "black",
            size = 4.2,
            stroke = 0.75
        ) +
        facet_wrap(
            ~panel,
            scales = "free_x",
            ncol = 1,
            labeller = as_labeller(panel_labels)
        ) +
        scale_shape_manual(
            values = c(
                "Tier 1" = 16,
                "Tier 2" = 17
            ),
            name = "Comparative tier"
        ) +
        coord_cartesian(
            ylim = c(0, 1)
        ) +
        labs(
            title = "Temporary focal pseudo-pool evidence",
            subtitle = paste(
                "Orange: left T2+T3; blue: right T2+T3;",
                "black rings mark formal or allele-balance",
                "support for ALT-enriched left and balanced right"
            ),
            x = "Hmr18 position (Mb)",
            y = "Combined ALT-read fraction"
        ) +
        theme_gynandro

    # Select one display site per comparative-pattern group.
    focal$joint_rank <- ifelse(
        focal$exact_left_11_right_01_match,
        3,
        ifelse(
            focal$AB_match_threshold_0.80,
            2,
            1
        )
    )

    focal$tier_rank <- ifelse(
        focal$tier == "tier_1",
        1,
        2
    )

    ordered_focal <- focal[
        order(
            focal$group_id,
            -focal$joint_rank,
            focal$tier_rank,
            -focal$minimum_BA,
            -focal$minimum_MCC,
            focal$pos
        ),
    ]

    selected <- ordered_focal[
        !duplicated(ordered_focal$group_id),
    ]

    selected <- selected[
        order(
            selected$panel,
            selected$pos
        ),
    ]

    write.table(
        selected,
        file = selection_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE,
        na = ""
    )

    selected$site_label <- paste0(
        selected$group_id,
        "  ",
        selected$pos,
        "  ",
        ifelse(
            selected$tier == "tier_1",
            "T1",
            "T2"
        ),
        "  BA=",
        sprintf("%.3f", selected$minimum_BA)
    )

    selected$site_label <- factor(
        selected$site_label,
        levels = rev(selected$site_label)
    )

    selected_long <- rbind(
        data.frame(
            panel = selected$panel,
            site_label = selected$site_label,
            side = "Left",
            ALT_fraction = selected$left_alt_fraction,
            pseudo_GT = selected$left_pseudo_GT,
            pseudo_GQ = selected$left_pseudo_GQ,
            joint_class = selected$joint_evidence_class,
            stringsAsFactors = FALSE
        ),
        data.frame(
            panel = selected$panel,
            site_label = selected$site_label,
            side = "Right",
            ALT_fraction = selected$right_alt_fraction,
            pseudo_GT = selected$right_pseudo_GT,
            pseudo_GQ = selected$right_pseudo_GQ,
            joint_class = selected$joint_evidence_class,
            stringsAsFactors = FALSE
        )
    )

    selected_long$side <- factor(
        selected_long$side,
        levels = c("Left", "Right")
    )

    selected_long$cell_text <- paste0(
        selected_long$pseudo_GT,
        "\n",
        sprintf(
            "%.2f",
            selected_long$ALT_fraction
        )
    )

    p_selected <- ggplot(
        selected_long,
        aes(
            x = side,
            y = site_label,
            fill = ALT_fraction
        )
    ) +
        geom_tile(
            colour = "white",
            linewidth = 0.7
        ) +
        geom_text(
            aes(label = cell_text),
            size = 2.7,
            lineheight = 0.9
        ) +
        facet_grid(
            panel ~ .,
            scales = "free_y",
            space = "free_y",
            labeller = as_labeller(panel_labels)
        ) +
        scale_fill_gradientn(
            colours = c(
                "#440154",
                "#3B528B",
                "#21918C",
                "#5EC962",
                "#FDE725"
            ),
            limits = c(0, 1),
            name = "ALT-read\nfraction"
        ) +
        labs(
            title = "Joint comparative and focal evidence",
            subtitle = paste(
                "One display site per comparative genotype-pattern",
                "group; values show pseudo-GT and combined",
                "ALT-read fraction"
            ),
            x = NULL,
            y = NULL
        ) +
        theme_gynandro +
        theme(
            axis.text.x = element_text(
                face = "bold",
                size = 9
            ),
            axis.text.y = element_text(
                family = "mono",
                size = 7
            ),
            axis.ticks = element_blank(),
            panel.border = element_blank()
        )

    pdf(
        combined_pdf,
        width = 10,
        height = 11,
        onefile = TRUE,
        useDingbats = FALSE
    )

    print(p_association)
    print(p_balance)
    print(p_selected)

    dev.off()

    ggsave(
        association_svg,
        p_association,
        width = 10,
        height = 11,
        units = "in",
        device = "svg"
    )

    ggsave(
        balance_svg,
        p_balance,
        width = 10,
        height = 11,
        units = "in",
        device = "svg"
    )

    ggsave(
        selected_svg,
        p_selected,
        width = 9,
        height = 10,
        units = "in",
        device = "svg"
    )

    message("PRELIMINARY OPTIX EVIDENCE FIGURES COMPLETED")
    message("Combined PDF: ", combined_pdf)
    message("Association SVG: ", association_svg)
    message("Allele-balance SVG: ", balance_svg)
    message("Selected-site SVG: ", selected_svg)
    message("Figure-site table: ", selection_file)
}
