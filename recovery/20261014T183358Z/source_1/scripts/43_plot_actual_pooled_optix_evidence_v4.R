#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
})

project_dir <- "/scratch/martinlab/joe/GynandroCalling"

paired_file <- file.path(
    project_dir,
    "analysis/original_502/checkpoints/",
    "optix_original_502_paired_site_evidence.tsv"
)

joint_file <- file.path(
    project_dir,
    "analysis/original_502/checkpoints/",
    "optix_original_502_joint_focal_candidates.tsv"
)

output_dir <- file.path(
    project_dir,
    "analysis/original_502/figures"
)

strict_table <- file.path(
    project_dir,
    "analysis/original_502/checkpoints/",
    "optix_main_figure_strict_sites.tsv"
)

inclusive_table <- file.path(
    project_dir,
    "analysis/original_502/checkpoints/",
    "optix_main_figure_inclusive_sites.tsv"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

required <- c(paired_file, joint_file)
missing <- required[!file.exists(required)]

if (length(missing) > 0) {
    message("ERROR: required inputs are missing:")
    message(paste(missing, collapse = "\n"))
} else {
    paired <- read.delim(
        paired_file,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    joint <- read.delim(
        joint_file,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    logical_columns <- c(
        "comparative_stable",
        "both_missingness_ok",
        "joint_exact",
        "joint_exact_missingness_ok",
        "pooled_both_DP2_GQ20"
    )

    for (column in logical_columns) {
        if (column %in% names(joint)) {
            joint[[column]] <- (
                tolower(as.character(
                    joint[[column]]
                )) == "true"
            )
        }

        if (column %in% names(paired)) {
            paired[[column]] <- (
                tolower(as.character(
                    paired[[column]]
                )) == "true"
            )
        }
    }

    numeric_columns <- c(
        "pos",
        "raw_BA",
        "QC_BA",
        "minimum_BA",
        "raw_MCC",
        "QC_MCC",
        "minimum_MCC",
        "TN_Left_ALT_fraction",
        "TN_Right_ALT_fraction"
    )

    for (column in numeric_columns) {
        if (column %in% names(joint)) {
            joint[[column]] <- as.numeric(
                joint[[column]]
            )
        }

        if (column %in% names(paired)) {
            paired[[column]] <- as.numeric(
                paired[[column]]
            )
        }
    }

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

    # Inclusive: every site with stable comparative evidence
    # and exact high-quality pooled L 1/1 -> R 0/1 calls.
    inclusive <- joint[
        joint$joint_exact &
        joint$pooled_both_DP2_GQ20,
    ]

    # Strict: well-performing, missingness-unbiased sites.
    strict <- inclusive[
        inclusive$both_missingness_ok &
        inclusive$minimum_BA >= 0.80 &
        inclusive$minimum_MCC >= 0.55,
    ]

    # Band is known to show weaker correspondence across the
    # broad multispecies panel. Retain its strongest exact focal
    # site as a visibly documented exception.
    band_sites <- inclusive[
        inclusive$panel == "band",
    ]

    if (nrow(band_sites) > 0) {
        band_sites <- band_sites[
            order(
                -band_sites$minimum_BA,
                -band_sites$minimum_MCC,
                band_sites$pos
            ),
        ]

        band_exception <- band_sites[1, ]

        strict <- rbind(
            strict[strict$panel != "band", ],
            band_exception
        )
    }

    strict <- strict[
        order(
            match(strict$panel, panel_order),
            strict$pos
        ),
    ]

    inclusive <- inclusive[
        order(
            match(inclusive$panel, panel_order),
            inclusive$pos
        ),
    ]

    strict$selection_class <- ifelse(
        strict$panel == "band" &
        !strict$both_missingness_ok,
        "Band exception",
        "Strict joint evidence"
    )

    inclusive$selection_class <- ifelse(
        inclusive$both_missingness_ok,
        "Missingness-unbiased",
        "Missingness-flagged"
    )

    write.table(
        strict,
        strict_table,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE,
        na = ""
    )

    write.table(
        inclusive,
        inclusive_table,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = TRUE,
        na = ""
    )

    prepare_plot_data <- function(data) {
        data$panel <- factor(
            data$panel,
            levels = panel_order
        )

        data$position_Mb <- data$pos / 1e6

        # Polarize toward the allele associated with presence
        # of the relevant phenotype.
        data$left_pattern_fraction <- ifelse(
            data$selected_orientation ==
            "REF_carrier_predicts_presence",
            1 - data$TN_Left_ALT_fraction,
            data$TN_Left_ALT_fraction
        )

        data$right_pattern_fraction <- ifelse(
            data$selected_orientation ==
            "REF_carrier_predicts_presence",
            1 - data$TN_Right_ALT_fraction,
            data$TN_Right_ALT_fraction
        )

        left <- data.frame(
            panel = data$panel,
            pos = data$pos,
            position_Mb = data$position_Mb,
            side = "Left",
            pattern_fraction = (
                data$left_pattern_fraction
            ),
            stringsAsFactors = FALSE
        )

        right <- data.frame(
            panel = data$panel,
            pos = data$pos,
            position_Mb = data$position_Mb,
            side = "Right",
            pattern_fraction = (
                data$right_pattern_fraction
            ),
            stringsAsFactors = FALSE
        )

        points <- rbind(left, right)

        points$panel <- factor(
            points$panel,
            levels = panel_order
        )

        points$side <- factor(
            points$side,
            levels = c("Left", "Right")
        )

        list(
            sites = data,
            points = points
        )
    }

    figure_theme <- theme_classic(
        base_size = 10,
        base_family = "sans"
    ) +
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
                size = 10
            ),
            axis.title = element_text(size = 10),
            axis.text = element_text(
                size = 8,
                colour = "black"
            ),
            legend.position = "bottom",
            legend.title = element_blank(),
            legend.text = element_text(size = 9),
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

    make_main_plot <- function(
        data,
        title,
        subtitle
    ) {
        prepared <- prepare_plot_data(data)
        sites <- prepared$sites
        points <- prepared$points

        panel_background <- unique(
            data.frame(panel = sites$panel)
        )

        ggplot() +
            # Subtle balanced/heterozygote-like region.
            geom_rect(
                data = panel_background,
                aes(
                    xmin = -Inf,
                    xmax = Inf,
                    ymin = 0.40,
                    ymax = 0.60
                ),
                inherit.aes = FALSE,
                fill = "#DCEBF7",
                alpha = 0.32
            ) +
            geom_hline(
                yintercept = c(0, 0.5, 1),
                linewidth = 0.28,
                colour = "grey72"
            ) +
            # Actual pooled left/right observations only.
            geom_segment(
                data = sites,
                aes(
                    x = position_Mb,
                    xend = position_Mb,
                    y = left_pattern_fraction,
                    yend = right_pattern_fraction
                ),
                inherit.aes = FALSE,
                colour = "grey58",
                linewidth = 0.45
            ) +
            geom_point(
                data = points,
                aes(
                    x = position_Mb,
                    y = pattern_fraction,
                    colour = side
                ),
                size = 3.0
            ) +
            # Arrowheads indicate comparative support.
            geom_point(
                data = sites,
                aes(
                    x = position_Mb,
                    y = 1.065
                ),
                inherit.aes = FALSE,
                shape = 25,
                size = 2.7,
                stroke = 0.2,
                colour = "black",
                fill = "black"
            ) +
            facet_wrap(
                ~panel,
                scales = "free_x",
                ncol = 1,
                labeller = as_labeller(
                    panel_labels
                )
            ) +
            scale_colour_manual(
                values = c(
                    Left = "#D55E00",
                    Right = "#0072B2"
                )
            ) +
            scale_y_continuous(
                limits = c(-0.04, 1.10),
                breaks = c(0, 0.5, 1),
                labels = c(
                    "0\nnon-pattern\nhom.-like",
                    "0.5\nheterozygous\nor mixed",
                    "1\npattern\nhom.-like"
                ),
                expand = expansion(
                    mult = c(0, 0)
                )
            ) +
            scale_x_continuous(
                labels = scales::label_number(
                    accuracy = 0.0001
                ),
                expand = expansion(
                    mult = c(0.04, 0.04)
                )
            ) +
            labs(
                title = title,
                subtitle = subtitle,
                x = "Hmr18 position (Mb)",
                y = (
                    "Pattern-associated marker-allele fraction"
                ),
                caption = paste(
                    "Arrowheads mark SNPs supported by the",
                    "comparative analysis; points are actual",
                    "TN_Left and TN_Right pooled calls."
                )
            ) +
            figure_theme
    }

    strict_plot <- make_main_plot(
        strict,
        "Left-right allele balance at informative Optix markers",
        paste(
            "Conservative site set; the sole Band marker is",
            "retained as a missingness-flagged regional exception."
        )
    )

    inclusive_plot <- make_main_plot(
        inclusive,
        "Left-right allele balance at all joint Optix markers",
        paste(
            "Inclusive site set: stable comparative correspondence",
            "and exact high-quality left 1/1 to right 0/1 calls."
        )
    )

    strict_pdf <- file.path(
        output_dir,
        "optix_actual_pooled_strict_v4.pdf"
    )

    strict_svg <- file.path(
        output_dir,
        "optix_actual_pooled_strict_v4.svg"
    )

    inclusive_pdf <- file.path(
        output_dir,
        "optix_actual_pooled_inclusive_v4.pdf"
    )

    inclusive_svg <- file.path(
        output_dir,
        "optix_actual_pooled_inclusive_v4.svg"
    )

    ggsave(
        strict_pdf,
        strict_plot,
        width = 10,
        height = 9.5,
        units = "in",
        device = cairo_pdf
    )

    ggsave(
        strict_svg,
        strict_plot,
        width = 10,
        height = 9.5,
        units = "in",
        device = "svg"
    )

    ggsave(
        inclusive_pdf,
        inclusive_plot,
        width = 10,
        height = 9.5,
        units = "in",
        device = cairo_pdf
    )

    ggsave(
        inclusive_svg,
        inclusive_plot,
        width = 10,
        height = 9.5,
        units = "in",
        device = "svg"
    )

    # Narrow comparative trace for eventual placement beneath
    # each supplementary genotype plot.
    paired$panel <- factor(
        paired$panel,
        levels = panel_order
    )

    paired$position_Mb <- paired$pos / 1e6

    trace_raw <- data.frame(
        panel = paired$panel,
        position_Mb = paired$position_Mb,
        scheme = "Called GT",
        BA = paired$raw_BA,
        eligible = paired$both_score_eligible
    )

    trace_qc <- data.frame(
        panel = paired$panel,
        position_Mb = paired$position_Mb,
        scheme = "DP >= 2; GQ >= 20",
        BA = paired$QC_BA,
        eligible = paired$both_score_eligible
    )

    trace <- rbind(trace_raw, trace_qc)

    trace$scheme <- factor(
        trace$scheme,
        levels = c(
            "Called GT",
            "DP >= 2; GQ >= 20"
        )
    )

    trace_plot <- ggplot(
        trace,
        aes(
            x = position_Mb,
            y = BA,
            colour = scheme
        )
    ) +
        geom_hline(
            yintercept = 0.5,
            linewidth = 0.25,
            linetype = "dashed",
            colour = "grey60"
        ) +
        geom_point(
            aes(
                alpha = eligible
            ),
            size = 0.9
        ) +
        facet_wrap(
            ~panel,
            scales = "free_x",
            ncol = 1,
            labeller = as_labeller(
                panel_labels
            )
        ) +
        scale_colour_manual(
            values = c(
                "Called GT" = "#777777",
                "DP >= 2; GQ >= 20" = "#0072B2"
            )
        ) +
        scale_alpha_manual(
            values = c(
                "FALSE" = 0.15,
                "TRUE" = 0.80
            ),
            guide = "none"
        ) +
        coord_cartesian(
            ylim = c(0.45, 1.01)
        ) +
        labs(
            x = "Hmr18 position (Mb)",
            y = "Balanced accuracy",
            colour = NULL
        ) +
        figure_theme +
        theme(
            plot.title = element_blank(),
            plot.subtitle = element_blank(),
            plot.caption = element_blank()
        )

    trace_pdf <- file.path(
        output_dir,
        "optix_comparative_score_trace_v4.pdf"
    )

    trace_svg <- file.path(
        output_dir,
        "optix_comparative_score_trace_v4.svg"
    )

    ggsave(
        trace_pdf,
        trace_plot,
        width = 9,
        height = 6.5,
        units = "in",
        device = cairo_pdf
    )

    ggsave(
        trace_svg,
        trace_plot,
        width = 9,
        height = 6.5,
        units = "in",
        device = "svg"
    )

    message("ACTUAL POOLED OPTIX FIGURES COMPLETED")
    message("Strict sites: ", nrow(strict))
    message("Inclusive sites: ", nrow(inclusive))
    message("Strict PDF: ", strict_pdf)
    message("Inclusive PDF: ", inclusive_pdf)
    message("Comparative trace: ", trace_pdf)
    message("Strict table: ", strict_table)
    message("Inclusive table: ", inclusive_table)
}
