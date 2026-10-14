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

output_pdf <- file.path(
    figure_dir,
    "optix_R_BA90_focal_and_high_depth_pattern_AF_rows.pdf"
)

control_file <- file.path(
    project,
    "analysis/original_502/checkpoints",
    "optix_R_BA90_high_depth_control_samples.tsv"
)

summary_file <- file.path(
    project,
    "analysis/original_502/checkpoints",
    "optix_R_BA90_focal_and_high_depth_AF_summary.tsv"
)

if (!file.exists(input_file)) {
    stop(
        paste(
            "Missing genotype table:",
            input_file
        )
    )
}

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

d$pos <- as.integer(d$pos)
d$DP <- suppressWarnings(as.numeric(d$DP))
d$GQ <- suppressWarnings(as.numeric(d$GQ))


parse_ad_component <- function(value, component) {

    if (
        is.na(value) ||
        value == "" ||
        value == "." ||
        value == ".,."
    ) {
        return(NA_real_)
    }

    pieces <- strsplit(
        value,
        ",",
        fixed = TRUE
    )[[1]]

    if (length(pieces) < component) {
        return(NA_real_)
    }

    suppressWarnings(
        as.numeric(pieces[component])
    )
}


d$ref_depth <- vapply(
    d$AD,
    parse_ad_component,
    numeric(1),
    component = 1
)

d$alt_depth <- vapply(
    d$AD,
    parse_ad_component,
    numeric(1),
    component = 2
)

d$AD_total <- d$ref_depth + d$alt_depth

# Ray-pattern allele is ALT.
d$pattern_AF <- ifelse(
    !is.na(d$AD_total) & d$AD_total > 0,
    d$alt_depth / d$AD_total,
    NA_real_
)

d$depth_QC <- (
    !is.na(d$DP) &
    !is.na(d$GQ) &
    d$DP >= 2 &
    d$GQ >= 20
)


# ------------------------------------------------------------
# Choose three comparative controls using depth only
# ------------------------------------------------------------

comparative <- d[
    d$sample_type == "comparative",
    ,
    drop = FALSE
]

comparative_samples <- unique(comparative$sample)

control_statistics <- do.call(
    rbind,
    lapply(
        comparative_samples,
        function(sample_id) {

            x <- comparative[
                comparative$sample == sample_id,
                ,
                drop = FALSE
            ]

            usable_depth <- x$DP[
                !is.na(x$DP)
            ]

            median_depth <- if (
                length(usable_depth) > 0
            ) {
                median(usable_depth)
            } else {
                NA_real_
            }

            data.frame(
                sample = sample_id,
                phenotype = unique(x$phenotype)[1],
                taxon_name = unique(x$taxon_name)[1],
                data_type = unique(x$data_type)[1],
                sites = nrow(x),
                depth_recorded = sum(!is.na(x$DP)),
                depth_QC_sites = sum(x$depth_QC),
                median_DP = median_depth,
                mean_DP = mean(
                    x$DP,
                    na.rm = TRUE
                ),
                maximum_DP = max(
                    x$DP,
                    na.rm = TRUE
                ),
                stringsAsFactors = FALSE
            )
        }
    )
)

control_statistics <- control_statistics[
    is.finite(control_statistics$median_DP),
    ,
    drop = FALSE
]

control_statistics <- control_statistics[
    order(
        -control_statistics$median_DP,
        -control_statistics$depth_QC_sites,
        control_statistics$sample
    ),
    ,
    drop = FALSE
]

present_candidates <- control_statistics[
    control_statistics$phenotype == "1",
    ,
    drop = FALSE
]

absent_candidates <- control_statistics[
    control_statistics$phenotype == "0",
    ,
    drop = FALSE
]

if (
    nrow(present_candidates) == 0 ||
    nrow(absent_candidates) == 0
) {
    stop(
        "Could not identify both R-present and R-absent controls."
    )
}

selected_present <- present_candidates[1, , drop = FALSE]
selected_absent <- absent_candidates[1, , drop = FALSE]

already_selected <- c(
    selected_present$sample,
    selected_absent$sample
)

remaining_candidates <- control_statistics[
    !control_statistics$sample %in% already_selected,
    ,
    drop = FALSE
]

if (nrow(remaining_candidates) == 0) {
    stop(
        "No third high-depth comparative control was available."
    )
}

selected_extra <- remaining_candidates[1, , drop = FALSE]

selected_controls <- rbind(
    transform(
        selected_present,
        selection_role = "highest_median_DP_R_present"
    ),
    transform(
        selected_absent,
        selection_role = "highest_median_DP_R_absent"
    ),
    transform(
        selected_extra,
        selection_role = "highest_median_DP_remaining"
    )
)

write.table(
    selected_controls,
    file = control_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ------------------------------------------------------------
# Assemble focal and control rows
# ------------------------------------------------------------

selected_samples <- c(
    "TN_Left",
    "TN_Right",
    selected_controls$sample
)

plot_data <- d[
    d$sample %in% selected_samples,
    ,
    drop = FALSE
]

if (
    length(unique(plot_data$sample)) != 5
) {
    stop(
        paste(
            "Expected five plot samples but recovered",
            length(unique(plot_data$sample))
        )
    )
}

site_information <- unique(
    plot_data[
        ,
        c(
            "pos",
            "site_label"
        )
    ]
)

site_information <- site_information[
    order(site_information$pos),
    ,
    drop = FALSE
]

site_levels <- site_information$site_label

plot_data$site_label <- factor(
    plot_data$site_label,
    levels = site_levels
)

control_labels <- setNames(
    paste0(
        selected_controls$sample,
        " | R=",
        selected_controls$phenotype,
        " | ",
        selected_controls$taxon_name,
        " | median DP=",
        format(
            selected_controls$median_DP,
            trim = TRUE,
            digits = 4
        )
    ),
    selected_controls$sample
)

panel_labels <- c(
    "TN_Left" = "TN_Left | pooled left",
    "TN_Right" = "TN_Right | pooled right",
    control_labels
)

plot_data$sample_panel <- unname(
    panel_labels[plot_data$sample]
)

panel_order <- unname(
    panel_labels[selected_samples]
)

plot_data$sample_panel <- factor(
    plot_data$sample_panel,
    levels = panel_order
)

plot_data$sample_key <- factor(
    plot_data$sample,
    levels = selected_samples
)

sample_colours <- setNames(
    c(
        "#2B8CBE",
        "#D95F0E",
        "#7B3294",
        "#008B8B",
        "#C28E00"
    ),
    selected_samples
)


# ------------------------------------------------------------
# Focal mixing expectation bands
# ------------------------------------------------------------

mixing_bands <- data.frame(
    sample_panel = factor(
        c(
            panel_labels["TN_Left"],
            panel_labels["TN_Right"]
        ),
        levels = panel_order
    ),
    ymin = c(
        0.85,
        0.50
    ),
    ymax = c(
        1.00,
        0.65
    ),
    band_colour = c(
        "#2B8CBE",
        "#D95F0E"
    ),
    stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# Quantitative summary for the five rows
# ------------------------------------------------------------

row_summary <- do.call(
    rbind,
    lapply(
        selected_samples,
        function(sample_id) {

            x <- plot_data[
                plot_data$sample == sample_id,
                ,
                drop = FALSE
            ]

            canonical_distance <- vapply(
                x$pattern_AF,
                function(value) {

                    if (is.na(value)) {
                        return(NA_real_)
                    }

                    min(
                        abs(
                            value - c(
                                0,
                                0.5,
                                1
                            )
                        )
                    )
                },
                numeric(1)
            )

            data.frame(
                sample = sample_id,
                sample_panel = as.character(
                    unique(x$sample_panel)[1]
                ),
                phenotype = unique(x$phenotype)[1],
                sites = nrow(x),
                sites_with_AD = sum(
                    !is.na(x$pattern_AF)
                ),
                DP2_GQ20_sites = sum(
                    x$depth_QC
                ),
                median_DP = median(
                    x$DP,
                    na.rm = TRUE
                ),
                mean_DP = mean(
                    x$DP,
                    na.rm = TRUE
                ),
                median_pattern_AF = median(
                    x$pattern_AF,
                    na.rm = TRUE
                ),
                median_distance_from_0_0.5_1 = median(
                    canonical_distance,
                    na.rm = TRUE
                ),
                fraction_within_0.10_of_0_0.5_1 = mean(
                    canonical_distance <= 0.10,
                    na.rm = TRUE
                ),
                stringsAsFactors = FALSE
            )
        }
    )
)

write.table(
    row_summary,
    file = summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

p <- ggplot(
    plot_data,
    aes(
        x = site_label,
        y = pattern_AF,
        colour = sample_key,
        group = sample_key
    )
) +

    geom_rect(
        data = mixing_bands,
        aes(
            xmin = -Inf,
            xmax = Inf,
            ymin = ymin,
            ymax = ymax,
            fill = band_colour
        ),
        inherit.aes = FALSE,
        alpha = 0.12,
        colour = NA
    ) +

    geom_hline(
        yintercept = c(
            0,
            0.5,
            1
        ),
        colour = "#777777",
        linewidth = 0.28,
        linetype = "dotted"
    ) +

    geom_line(
        linewidth = 0.42,
        alpha = 0.68,
        na.rm = TRUE
    ) +

    geom_point(
        aes(
            size = pmax(
                AD_total,
                1
            )
        ),
        alpha = 0.90,
        na.rm = TRUE
    ) +

    facet_grid(
        rows = vars(sample_panel)
    ) +

    scale_colour_manual(
        values = sample_colours,
        guide = "none"
    ) +

    scale_fill_identity() +

    scale_size_continuous(
        name = "Allelic depth",
        range = c(
            0.65,
            2.9
        )
    ) +

    scale_y_continuous(
        name = "Pattern-associated allele fraction",
        limits = c(
            0,
            1
        ),
        breaks = c(
            0,
            0.25,
            0.5,
            0.75,
            1
        ),
        expand = expansion(
            mult = c(
                0.01,
                0.02
            )
        )
    ) +

    scale_x_discrete(
        name = "Hmr18 position",
        expand = expansion(
            mult = c(
                0.001,
                0.001
            )
        )
    ) +

    labs(
        title = paste(
            "Pattern-allele balance in pooled focal and",
            "high-depth comparative samples"
        ),
        subtitle = paste(
            "Ray-associated SNPs with minimum comparative",
            "balanced accuracy >= 0.90"
        ),
        caption = paste(
            "Point size represents allelic depth.",
            "Comparative controls were selected using median DP only.",
            "The blue and orange shaded bands apply only to TN_Left",
            "and TN_Right and show the expected observations for a",
            "true P/P to P/A transition under up to 30% reciprocal mixing."
        )
    ) +

    theme_minimal(
        base_size = 9
    ) +

    theme(
        panel.grid.minor = element_blank(),

        panel.grid.major.x = element_blank(),

        axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5,
            size = 5.2
        ),

        strip.text.y = element_text(
            angle = 0,
            hjust = 0,
            face = "bold",
            size = 8
        ),

        strip.background = element_rect(
            fill = "#F3F3F3",
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
            size = 7.2,
            colour = "#444444",
            hjust = 0
        ),

        plot.margin = margin(
            10,
            14,
            10,
            12
        )
    )

site_count <- length(site_levels)

plot_width <- max(
    15,
    min(
        30,
        5 + 0.16 * site_count
    )
)

ggsave(
    filename = output_pdf,
    plot = p,
    device = cairo_pdf,
    width = plot_width,
    height = 14,
    units = "in",
    limitsize = FALSE
)

cat("FOCAL AND HIGH-DEPTH COMPARATIVE ALLELE-FRACTION ROWS\n\n")

cat("CONTROL SELECTION\n")
print(
    selected_controls[
        ,
        c(
            "selection_role",
            "sample",
            "phenotype",
            "taxon_name",
            "median_DP",
            "depth_QC_sites"
        )
    ],
    row.names = FALSE
)

cat("\nROW SUMMARY\n")
print(row_summary, row.names = FALSE)

cat("\nPDF:\n", output_pdf, "\n")
cat("\nControl table:\n", control_file, "\n")
cat("\nSummary:\n", summary_file, "\n")
