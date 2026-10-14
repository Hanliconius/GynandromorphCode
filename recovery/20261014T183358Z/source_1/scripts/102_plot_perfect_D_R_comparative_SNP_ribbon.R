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
    "drop_problematic_reanalysis"
)

dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

suppressPackageStartupMessages({
    library(ggplot2)
})


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


site_key <- function(d) {

    paste(
        d$chrom,
        d$pos,
        d$ref,
        d$alt,
        sep = ":"
    )
}


load_perfect_sites <- function(
    module,
    trait,
    prefix
) {

    score_file <- file.path(
        checkpoint_dir,
        paste0(
            prefix,
            "_original502_expected_pattern_sites.",
            "WGS_vs_WGSplusCS.tsv"
        )
    )

    paired_file <- file.path(
        checkpoint_dir,
        paste0(
            prefix,
            "_original502_allSNP_scores.filtered.",
            "WGS_vs_WGSplusCS.minimum_across_schemes.tsv"
        )
    )

    if (!file.exists(score_file)) {
        stop(
            paste(
                "Missing score file:",
                score_file
            )
        )
    }

    if (!file.exists(paired_file)) {
        stop(
            paste(
                "Missing paired-scheme file:",
                paired_file
            )
        )
    }

    scores <- read.delim(
        score_file,
        sep = "\t",
        header = TRUE,
        quote = "",
        comment.char = "",
        check.names = FALSE
    )

    scores$balanced_accuracy <- suppressWarnings(
        as.numeric(
            scores$balanced_accuracy
        )
    )

    keep <- (
        scores$analysis_cohort == "WGS_plus_CS"
        & scores$scheme == "DP2_GQ20"
        & !is.na(scores$balanced_accuracy)
        & scores$balanced_accuracy >= (
            1 - 1e-12
        )
    )

    if (
        "expected_pattern_eligible"
        %in% names(scores)
    ) {
        keep <- (
            keep
            & is_true(
                scores$expected_pattern_eligible
            )
        )
    }

    perfect <- scores[
        keep,
        ,
        drop = FALSE
    ]

    paired <- read.delim(
        paired_file,
        sep = "\t",
        header = TRUE,
        quote = "",
        comment.char = "",
        check.names = FALSE
    )

    paired$minimum_BA_across_schemes <- suppressWarnings(
        as.numeric(
            paired$minimum_BA_across_schemes
        )
    )

    paired_keep <- (
        paired$analysis_cohort == "WGS_plus_CS"
        & !is.na(
            paired$minimum_BA_across_schemes
        )
        & paired$minimum_BA_across_schemes >= (
            1 - 1e-12
        )
    )

    if (
        "all_schemes_expected_pattern_eligible"
        %in% names(paired)
    ) {
        paired_keep <- (
            paired_keep
            & is_true(
                paired[[
                    "all_schemes_expected_pattern_eligible"
                ]]
            )
        )
    }

    paired_perfect <- paired[
        paired_keep,
        ,
        drop = FALSE
    ]

    both_scheme_keys <- site_key(
        paired_perfect
    )

    perfect$module <- module
    perfect$trait <- trait
    perfect$site <- site_key(
        perfect
    )
    perfect$perfect_both_schemes <- (
        perfect$site
        %in% both_scheme_keys
    )

    perfect <- perfect[
        order(
            perfect$pos
        ),
        ,
        drop = FALSE
    ]

    perfect
}


dennis <- load_perfect_sites(
    module = "Dennis",
    trait = "D",
    prefix = "longDennis"
)

ray <- load_perfect_sites(
    module = "Ray",
    trait = "R",
    prefix = "longRay"
)

perfect <- rbind(
    dennis,
    ray
)

if (nrow(perfect) == 0) {
    stop(
        "No perfect Dennis or Ray SNPs were found."
    )
}

perfect$pos <- as.numeric(
    perfect$pos
)

perfect$position_Mb <- (
    perfect$pos / 1e6
)

perfect$row_y <- ifelse(
    perfect$module == "Dennis",
    2,
    1
)

perfect$scheme_class <- ifelse(
    perfect$perfect_both_schemes,
    "Perfect in both schemes",
    "Perfect after DP2/GQ20"
)


broad_start <- 938195
broad_end <- 1239421

ribbon_rows <- data.frame(
    module = c(
        "Dennis",
        "Ray"
    ),
    row_y = c(
        2,
        1
    ),
    xmin = broad_start / 1e6,
    xmax = broad_end / 1e6
)

prior_intervals <- data.frame(
    module = c(
        "Dennis",
        "Dennis",
        "Ray",
        "Ray"
    ),
    interval = c(
        "Dennis1",
        "Dennis2",
        "Ray1",
        "Ray2"
    ),
    start = c(
        1215210,
        1218557,
        1156746,
        1174547
    ),
    end = c(
        1216124,
        1221494,
        1160220,
        1187833
    ),
    row_y = c(
        2,
        2,
        1,
        1
    )
)

prior_intervals$xmin <- (
    prior_intervals$start / 1e6
)

prior_intervals$xmax <- (
    prior_intervals$end / 1e6
)

prior_intervals$xmid <- (
    (
        prior_intervals$start
        + prior_intervals$end
    )
    / 2e6
)

summary_table <- data.frame(
    module = c(
        "Dennis",
        "Ray"
    ),
    trait = c(
        "D",
        "R"
    ),
    perfect_DP2_GQ20 = c(
        nrow(dennis),
        nrow(ray)
    ),
    perfect_both_schemes = c(
        sum(
            dennis$perfect_both_schemes
        ),
        sum(
            ray$perfect_both_schemes
        )
    ),
    minimum_position = c(
        min(dennis$pos),
        min(ray$pos)
    ),
    maximum_position = c(
        max(dennis$pos),
        max(ray$pos)
    )
)

count_labels <- data.frame(
    module = summary_table$module,
    row_y = c(
        2,
        1
    ),
    label = paste0(
        "n = ",
        summary_table$perfect_DP2_GQ20,
        " perfect"
    )
)


module_colours <- c(
    "Dennis" = "#E69F00",
    "Ray" = "#C51B7D"
)


plot_object <- ggplot() +

    geom_rect(
        data = ribbon_rows,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = row_y - 0.12,
            ymax = row_y + 0.12
        ),
        fill = "#E6E6E6",
        colour = "#555555",
        linewidth = 0.35
    ) +

    geom_rect(
        data = prior_intervals,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = row_y - 0.12,
            ymax = row_y + 0.12,
            fill = module
        ),
        alpha = 0.38,
        colour = NA
    ) +

    geom_segment(
        data = perfect,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = row_y - 0.18,
            yend = row_y + 0.18,
            colour = module
        ),
        linewidth = 0.55,
        lineend = "butt"
    ) +

    geom_point(
        data = perfect[
            perfect$perfect_both_schemes,
            ,
            drop = FALSE
        ],
        aes(
            x = position_Mb,
            y = row_y + 0.25
        ),
        shape = 21,
        size = 2.0,
        stroke = 0.35,
        fill = "#1A1A1A",
        colour = "#1A1A1A"
    ) +

    geom_text(
        data = prior_intervals,
        aes(
            x = xmid,
            y = row_y - 0.28,
            label = interval
        ),
        size = 2.8,
        colour = "#333333"
    ) +

    geom_text(
        data = count_labels,
        aes(
            x = broad_end / 1e6,
            y = row_y + 0.40,
            label = label,
            colour = module
        ),
        hjust = 1,
        size = 3.1,
        fontface = "bold"
    ) +

    scale_colour_manual(
        values = module_colours,
        breaks = c(
            "Dennis",
            "Ray"
        ),
        name = "Comparative phenotype"
    ) +

    scale_fill_manual(
        values = module_colours,
        guide = "none"
    ) +

    scale_x_continuous(
        name = "Position on Hmr18 (Mb)",
        limits = c(
            broad_start / 1e6,
            broad_end / 1e6
        ),
        breaks = seq(
            0.95,
            1.25,
            by = 0.05
        ),
        expand = expansion(
            mult = c(
                0.005,
                0.005
            )
        )
    ) +

    scale_y_continuous(
        name = NULL,
        breaks = c(
            1,
            2
        ),
        labels = c(
            "Ray (R)",
            "Dennis (D)"
        ),
        limits = c(
            0.55,
            2.55
        ),
        expand = expansion(
            mult = c(
                0,
                0
            )
        )
    ) +

    labs(
        title = paste(
            "Perfect comparative Dennis and Ray",
            "genotype–phenotype markers"
        ),
        subtitle = paste(
            "Ticks: WGS+CS balanced accuracy = 1.0",
            "after DP >= 2 and GQ >= 20;",
            "black dots: also perfect using all called genotypes"
        ),
        caption = paste(
            "Shaded blocks mark the original Dennis1,",
            "Dennis2, Ray1 and Ray2 plotting intervals.",
            "The gynandromorph was not used to select SNPs."
        )
    ) +

    coord_cartesian(
        clip = "off"
    ) +

    theme_minimal(
        base_size = 10
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
                b = 10
            )
        ),
        plot.caption = element_text(
            size = 8,
            colour = "#555555",
            hjust = 0
        ),
        axis.text.y = element_text(
            face = "bold",
            colour = "#222222",
            size = 10
        ),
        axis.text.x = element_text(
            colour = "#333333"
        ),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(
            colour = "#E1E1E1",
            linewidth = 0.3
        ),
        legend.position = "top",
        legend.justification = "left",
        legend.title = element_text(
            face = "bold"
        ),
        plot.margin = margin(
            12,
            15,
            8,
            12
        )
    )


perfect_output <- file.path(
    checkpoint_dir,
    "perfect_D_R_comparative_SNPs.DP2_GQ20.tsv"
)

summary_output <- file.path(
    checkpoint_dir,
    "perfect_D_R_comparative_SNPs.summary.tsv"
)

pdf_output <- file.path(
    figure_dir,
    "perfect_D_R_comparative_SNP_ribbon.pdf"
)

svg_output <- file.path(
    figure_dir,
    "perfect_D_R_comparative_SNP_ribbon.svg"
)

write.table(
    perfect,
    file = perfect_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    summary_table,
    file = summary_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

ggsave(
    filename = pdf_output,
    plot = plot_object,
    width = 11,
    height = 4.2,
    units = "in",
    device = cairo_pdf
)

ggsave(
    filename = svg_output,
    plot = plot_object,
    width = 11,
    height = 4.2,
    units = "in",
    device = "svg"
)

cat(
    "PERFECT D/R COMPARATIVE SNP RIBBON\n\n"
)

print(
    summary_table,
    row.names = FALSE
)

cat(
    "\nPerfect-site table:\n",
    perfect_output,
    "\n\nPDF:\n",
    pdf_output,
    "\n\nSVG:\n",
    svg_output,
    "\n",
    sep = ""
)
