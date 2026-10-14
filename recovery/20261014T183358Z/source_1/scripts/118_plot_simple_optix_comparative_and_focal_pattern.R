#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(grid)
})


project_dir <- "/scratch/martinlab/joe/GynandroCalling"

checkpoint_dir <- file.path(
    project_dir,
    "analysis",
    "original_502",
    "checkpoints"
)

selected_file <- file.path(
    checkpoint_dir,
    "broad_optix_strict_comparative_SNPs.tsv"
)

focal_file <- file.path(
    checkpoint_dir,
    "core_explanatory_optix_genotypes.long.tsv"
)

output_dir <- file.path(
    project_dir,
    "analysis",
    "figures",
    "optix_simple_comparative_focal_pattern"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

pdf_file <- file.path(
    output_dir,
    "optix_comparative_BA_and_focal_DR_B_pattern.pdf"
)

png_file <- file.path(
    output_dir,
    "optix_comparative_BA_and_focal_DR_B_pattern.png"
)

summary_file <- file.path(
    checkpoint_dir,
    "optix_simple_comparative_focal_pattern.summary.tsv"
)

site_file <- file.path(
    checkpoint_dir,
    "optix_simple_comparative_focal_pattern.sites.tsv"
)

background_manifest_file <- file.path(
    checkpoint_dir,
    "optix_simple_comparative_background_files.tsv"
)


if (!file.exists(selected_file)) {
    stop(
        paste(
            "Missing comparative-selected SNP table:",
            selected_file
        )
    )
}

if (!file.exists(focal_file)) {
    stop(
        paste(
            "Missing focal genotype table:",
            focal_file
        )
    )
}


as_bool <- function(x) {
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


canonical_gt <- function(x) {
    x <- gsub(
        "\\|",
        "/",
        as.character(x)
    )

    x[
        x %in% c(
            "",
            ".",
            "./.",
            ".|.",
            "NA"
        )
    ] <- NA_character_

    x[
        x == "1/0"
    ] <- "0/1"

    x
}


normalise_module <- function(x) {
    x <- as.character(x)

    result <- rep(
        NA_character_,
        length(x)
    )

    result[
        grepl(
            "band",
            x,
            ignore.case = TRUE
        )
    ] <- "Band"

    result[
        grepl(
            "dennis",
            x,
            ignore.case = TRUE
        )
    ] <- "Dennis"

    result[
        grepl(
            "ray",
            x,
            ignore.case = TRUE
        )
    ] <- "Ray"

    result
}


safe_min <- function(x) {
    x <- suppressWarnings(
        as.numeric(x)
    )

    x <- x[
        is.finite(x)
    ]

    if (length(x) == 0) {
        return(
            NA_real_
        )
    }

    min(x)
}


selected <- read.delim(
    selected_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

focal <- read.delim(
    focal_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
)


required_selected <- c(
    "module",
    "chrom",
    "pos",
    "ref",
    "alt",
    "minimum_BA_all_analyses"
)

missing_selected <- setdiff(
    required_selected,
    names(selected)
)

if (length(missing_selected) > 0) {
    stop(
        paste(
            "Comparative-selected table lacks columns:",
            paste(
                missing_selected,
                collapse = ", "
            )
        )
    )
}


required_focal <- c(
    "module",
    "chrom",
    "pos",
    "ref",
    "alt",
    "sample",
    "canonical_GT",
    "pattern_dosage",
    "DP",
    "GQ",
    "DP2_GQ20_pass"
)

missing_focal <- setdiff(
    required_focal,
    names(focal)
)

if (length(missing_focal) > 0) {
    stop(
        paste(
            "Focal genotype table lacks columns:",
            paste(
                missing_focal,
                collapse = ", "
            )
        )
    )
}


selected$module <- normalise_module(
    selected$module
)

focal$module <- normalise_module(
    focal$module
)

selected <- selected[
    !is.na(selected$module),
    ,
    drop = FALSE
]

focal <- focal[
    !is.na(focal$module),
    ,
    drop = FALSE
]

selected$pos <- as.numeric(
    selected$pos
)

selected$minimum_BA_all_analyses <- as.numeric(
    selected$minimum_BA_all_analyses
)

selected$position_Mb <- (
    selected$pos / 1000000
)

selected$site_key <- paste(
    selected$module,
    selected$chrom,
    selected$pos,
    selected$ref,
    selected$alt,
    sep = ":"
)

if (anyDuplicated(selected$site_key)) {
    stop(
        "Comparative-selected SNP table contains duplicate module-site rows."
    )
}


focal$pos <- as.numeric(
    focal$pos
)

focal$canonical_GT <- canonical_gt(
    focal$canonical_GT
)

focal$pattern_dosage <- as.numeric(
    focal$pattern_dosage
)

focal$DP <- as.numeric(
    focal$DP
)

focal$GQ <- as.numeric(
    focal$GQ
)

focal$QC <- as_bool(
    focal$DP2_GQ20_pass
)

focal$site_key <- paste(
    focal$module,
    focal$chrom,
    focal$pos,
    focal$ref,
    focal$alt,
    sep = ":"
)

focal <- focal[
    focal$sample %in% c(
        "TN_Left",
        "TN_Right"
    ),
    ,
    drop = FALSE
]


left <- focal[
    focal$sample == "TN_Left",
    c(
        "site_key",
        "canonical_GT",
        "pattern_dosage",
        "DP",
        "GQ",
        "QC"
    ),
    drop = FALSE
]

right <- focal[
    focal$sample == "TN_Right",
    c(
        "site_key",
        "canonical_GT",
        "pattern_dosage",
        "DP",
        "GQ",
        "QC"
    ),
    drop = FALSE
]

names(left) <- c(
    "site_key",
    "left_GT",
    "left_pattern_dosage",
    "left_DP",
    "left_GQ",
    "left_QC"
)

names(right) <- c(
    "site_key",
    "right_GT",
    "right_pattern_dosage",
    "right_DP",
    "right_GQ",
    "right_QC"
)


if (anyDuplicated(left$site_key)) {
    stop(
        "Duplicate TN_Left module-site rows were found."
    )
}

if (anyDuplicated(right$site_key)) {
    stop(
        "Duplicate TN_Right module-site rows were found."
    )
}


site_data <- merge(
    selected,
    left,
    by = "site_key",
    all.x = TRUE,
    sort = FALSE
)

site_data <- merge(
    site_data,
    right,
    by = "site_key",
    all.x = TRUE,
    sort = FALSE
)


missing_left <- sum(
    is.na(site_data$left_GT)
    & is.na(site_data$left_DP)
)

missing_right <- sum(
    is.na(site_data$right_GT)
    & is.na(site_data$right_DP)
)

if (
    missing_left > 0
    || missing_right > 0
) {
    stop(
        paste(
            "Missing focal rows after joining:",
            missing_left,
            "left and",
            missing_right,
            "right."
        )
    )
}


# For Dennis and Ray, pattern dosage already measures the
# D- or R-associated allele:
#     expected left = 1, right = 0.5.
#
# For Band, pattern dosage measures the B-associated allele:
#     expected left = 0, right = 0.5.
#
# Invert Band so every module shares one DR-associated axis:
#     expected left = 1, right = 0.5.

site_data$left_DR_dosage <- ifelse(
    site_data$module == "Band",
    1 - site_data$left_pattern_dosage,
    site_data$left_pattern_dosage
)

site_data$right_DR_dosage <- ifelse(
    site_data$module == "Band",
    1 - site_data$right_pattern_dosage,
    site_data$right_pattern_dosage
)

site_data$both_QC <- (
    site_data$left_QC
    & site_data$right_QC
)

site_data$identical_GT <- (
    site_data$both_QC
    & !is.na(site_data$left_GT)
    & !is.na(site_data$right_GT)
    & site_data$left_GT
        == site_data$right_GT
)

site_data$different_GT <- (
    site_data$both_QC
    & !is.na(site_data$left_GT)
    & !is.na(site_data$right_GT)
    & site_data$left_GT
        != site_data$right_GT
)

site_data$expected_focal_pattern <- (
    site_data$both_QC
    & is.finite(
        site_data$left_DR_dosage
    )
    & is.finite(
        site_data$right_DR_dosage
    )
    & abs(
        site_data$left_DR_dosage - 1
    ) < 1e-8
    & abs(
        site_data$right_DR_dosage - 0.5
    ) < 1e-8
)

site_data$other_difference <- (
    site_data$different_GT
    & !site_data$expected_focal_pattern
)

site_data$focal_state <- "One or both calls fail QC"

site_data$focal_state[
    site_data$identical_GT
] <- "Identical left/right genotype"

site_data$focal_state[
    site_data$other_difference
] <- "Other left/right difference"

site_data$focal_state[
    site_data$expected_focal_pattern
] <- "Expected left DR/DR; right DR/B"


module_levels <- c(
    "Band",
    "Dennis",
    "Ray"
)

site_data$module <- factor(
    site_data$module,
    levels = module_levels
)

site_data <- site_data[
    order(
        site_data$module,
        site_data$pos
    ),
    ,
    drop = FALSE
]


write.table(
    site_data,
    file = site_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


summary_rows <- do.call(
    rbind,
    lapply(
        module_levels,
        function(module_name) {
            d <- site_data[
                site_data$module
                    == module_name,
                ,
                drop = FALSE
            ]

            data.frame(
                module = module_name,
                comparative_selected_SNPs = nrow(d),
                both_focal_QC = sum(
                    d$both_QC
                ),
                identical_left_right = sum(
                    d$identical_GT
                ),
                different_left_right = sum(
                    d$different_GT
                ),
                expected_left1_right05 = sum(
                    d$expected_focal_pattern
                ),
                other_left_right_difference = sum(
                    d$other_difference
                ),
                one_or_both_fail_QC = sum(
                    !d$both_QC
                ),
                expected_fraction_of_all_selected = if (
                    nrow(d) > 0
                ) {
                    sum(
                        d$expected_focal_pattern
                    ) / nrow(d)
                } else {
                    NA_real_
                },
                expected_fraction_of_differences = if (
                    sum(
                        d$different_GT
                    ) > 0
                ) {
                    sum(
                        d$expected_focal_pattern
                    ) / sum(
                        d$different_GT
                    )
                } else {
                    NA_real_
                },
                stringsAsFactors = FALSE
            )
        }
    )
)


unique_expected <- unique(
    site_data[
        site_data$expected_focal_pattern,
        c(
            "chrom",
            "pos",
            "ref",
            "alt"
        ),
        drop = FALSE
    ]
)

unique_differences <- unique(
    site_data[
        site_data$different_GT,
        c(
            "chrom",
            "pos",
            "ref",
            "alt"
        ),
        drop = FALSE
    ]
)


write.table(
    summary_rows,
    file = summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ------------------------------------------------------------
# Recover the frozen all-SNP BA background where available.
# This is presentation-only; selected SNPs always come from
# broad_optix_strict_comparative_SNPs.tsv.
# ------------------------------------------------------------

load_background <- function(
    module_name
) {
    pattern <- paste0(
        "^long",
        module_name,
        ".*minimum_across_schemes\\.tsv$"
    )

    files <- list.files(
        checkpoint_dir,
        pattern = pattern,
        full.names = TRUE,
        ignore.case = TRUE
    )

    files <- files[
        !grepl(
            "summary",
            basename(files),
            ignore.case = TRUE
        )
    ]

    if (length(files) == 0) {
        return(
            list(
                data = data.frame(),
                files = character()
            )
        )
    }

    input_rows <- list()

    for (index in seq_along(files)) {
        d <- read.delim(
            files[index],
            sep = "\t",
            header = TRUE,
            quote = "",
            comment.char = "",
            stringsAsFactors = FALSE,
            check.names = FALSE
        )

        required <- c(
            "chrom",
            "pos",
            "ref",
            "alt",
            "minimum_BA_across_schemes"
        )

        if (
            !all(
                required %in% names(d)
            )
        ) {
            next
        }

        d$source_file <- basename(
            files[index]
        )

        input_rows[[length(input_rows) + 1]
        ] <- d
    }

    if (length(input_rows) == 0) {
        return(
            list(
                data = data.frame(),
                files = files
            )
        )
    }

    combined <- do.call(
        rbind,
        input_rows
    )

    combined$pos <- as.numeric(
        combined$pos
    )

    combined$minimum_BA_across_schemes <- as.numeric(
        combined$minimum_BA_across_schemes
    )

    collapsed <- aggregate(
        minimum_BA_across_schemes
            ~ chrom + pos + ref + alt,
        data = combined,
        FUN = safe_min
    )

    names(collapsed)[
        names(collapsed)
            == "minimum_BA_across_schemes"
    ] <- "minimum_BA_all_available_cohorts"

    collapsed$module <- module_name

    list(
        data = collapsed,
        files = files
    )
}


background_results <- lapply(
    module_levels,
    load_background
)

background_tables <- lapply(
    background_results,
    function(x) {
        x$data
    }
)

background_tables <- background_tables[
    vapply(
        background_tables,
        nrow,
        integer(1)
    ) > 0
]

if (length(background_tables) > 0) {
    ba_background <- do.call(
        rbind,
        background_tables
    )

    ba_background$position_Mb <- (
        ba_background$pos / 1000000
    )

    ba_background$module <- factor(
        ba_background$module,
        levels = module_levels
    )
} else {
    ba_background <- data.frame(
        module = factor(
            levels = module_levels
        ),
        position_Mb = numeric(),
        minimum_BA_all_available_cohorts = numeric()
    )
}


background_manifest <- do.call(
    rbind,
    lapply(
        seq_along(module_levels),
        function(index) {
            files <- background_results[[index]
            ]$files

            if (length(files) == 0) {
                return(
                    data.frame(
                        module = module_levels[index],
                        file = "none recovered",
                        stringsAsFactors = FALSE
                    )
                )
            }

            data.frame(
                module = module_levels[index],
                file = files,
                stringsAsFactors = FALSE
            )
        }
    )
)

write.table(
    background_manifest,
    file = background_manifest_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ------------------------------------------------------------
# Plotting data: show only left/right differences below.
# Identical and QC-failed calls remain in the summary table.
# ------------------------------------------------------------

pair_data <- site_data[
    site_data$different_GT,
    ,
    drop = FALSE
]

expected_pairs <- pair_data[
    pair_data$expected_focal_pattern,
    ,
    drop = FALSE
]

other_pairs <- pair_data[
    pair_data$other_difference,
    ,
    drop = FALSE
]


focal_long <- rbind(
    data.frame(
        module = pair_data$module,
        position_Mb = pair_data$position_Mb,
        side = "TN_Left",
        DR_dosage = pair_data$left_DR_dosage,
        expected = pair_data$expected_focal_pattern,
        stringsAsFactors = FALSE
    ),
    data.frame(
        module = pair_data$module,
        position_Mb = pair_data$position_Mb,
        side = "TN_Right",
        DR_dosage = pair_data$right_DR_dosage,
        expected = pair_data$expected_focal_pattern,
        stringsAsFactors = FALSE
    )
)

focal_long$module <- factor(
    focal_long$module,
    levels = module_levels
)

focal_long$side <- factor(
    focal_long$side,
    levels = c(
        "TN_Left",
        "TN_Right"
    )
)

expected_points <- focal_long[
    focal_long$expected,
    ,
    drop = FALSE
]

other_points <- focal_long[
    !focal_long$expected,
    ,
    drop = FALSE
]


module_colours <- c(
    "Band" = "#3C9DC7",
    "Dennis" = "#B51A7D",
    "Ray" = "#E88924"
)

side_colours <- c(
    "TN_Left" = "#0072B2",
    "TN_Right" = "#D55E00"
)

x_limits <- c(
    1.060,
    1.240
)

gene_annotation <- data.frame(
    module = factor(
        module_levels,
        levels = module_levels
    ),
    xmin = 1088602 / 1000000,
    xmax = 1097846 / 1000000
)

intervals <- data.frame(
    module = factor(
        c(
            "Band",
            "Dennis",
            "Dennis",
            "Ray",
            "Ray"
        ),
        levels = module_levels
    ),
    interval = c(
        "B",
        "D1",
        "D2",
        "R1",
        "R2"
    ),
    xmin = c(
        1112651,
        1215210,
        1218557,
        1156746,
        1174547
    ) / 1000000,
    xmax = c(
        1114314,
        1216124,
        1221494,
        1160220,
        1187833
    ) / 1000000
)



# KRISTA_PLOT_WINDOWS_DEFINITION
#
# Original regional windows used for Krista's genotype plots.
# Coordinates are converted to Mb because the aligned figure uses Mb.
#
krista_regions <- data.frame(
    module = c(
        "Band",
        "Dennis",
        "Dennis",
        "Ray",
        "Ray"
    ),
    region_label = c(
        "B",
        "D1",
        "D2",
        "R1",
        "R2"
    ),
    xmin_mb = c(
        1112651,
        1215210,
        1218557,
        1156746,
        1174547
    ) / 1e6,
    xmax_mb = c(
        1114314,
        1216124,
        1221494,
        1160220,
        1187833
    ) / 1e6,
    stringsAsFactors = FALSE
)

krista_regions$mid_mb <- (
    krista_regions$xmin_mb +
    krista_regions$xmax_mb
) / 2

ba_plot <- ggplot() +

    # KRISTA_PLOT_WINDOWS_LAYER
    geom_rect(
        data = krista_regions,
        aes(
            xmin = xmin_mb,
            xmax = xmax_mb,
            ymin = -Inf,
            ymax = Inf
        ),
        inherit.aes = FALSE,
        fill = "#6F6F6F",
        colour = "#505050",
        alpha = 0.13,
        linewidth = 0.28
    ) +

    geom_text(
        data = krista_regions,
        aes(
            x = mid_mb,
            y = 1.008,
            label = region_label
        ),
        inherit.aes = FALSE,
        colour = "#3A3A3A",
        size = 2.7,
        fontface = "bold",
        vjust = 0.5
    ) +

    geom_rect(
        data = gene_annotation,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = -Inf,
            ymax = Inf
        ),
        inherit.aes = FALSE,
        fill = "#56B4E9",
        alpha = 0.10
    ) +

    geom_rect(
        data = intervals,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = 0.48,
            ymax = 0.505
        ),
        inherit.aes = FALSE,
        fill = "#777777",
        alpha = 0.24
    ) +

    geom_text(
        data = intervals,
        aes(
            x = (
                xmin + xmax
            ) / 2,
            y = 0.518,
            label = interval
        ),
        inherit.aes = FALSE,
        size = 2.45,
        colour = "#555555"
    ) +

    geom_hline(
        yintercept = 0.50,
        colour = "#777777",
        linewidth = 0.32
    ) +

    geom_hline(
        yintercept = 0.80,
        colour = "#888888",
        linewidth = 0.30,
        linetype = "dotted"
    ) +

    geom_hline(
        yintercept = 0.90,
        colour = "#888888",
        linewidth = 0.30,
        linetype = "dotted"
    ) +

    geom_hline(
        yintercept = 1.00,
        colour = "#555555",
        linewidth = 0.32
    ) +

    geom_point(
        data = ba_background,
        aes(
            x = position_Mb,
            y = minimum_BA_all_available_cohorts
        ),
        colour = "#BDBDBD",
        alpha = 0.16,
        size = 0.34,
        stroke = 0
    ) +

    geom_point(
        data = site_data,
        aes(
            x = position_Mb,
            y = minimum_BA_all_analyses,
            colour = module
        ),
        shape = 25,
        size = 1.55,
        stroke = 0.20,
        alpha = 0.80
    ) +

    geom_text(
        data = gene_annotation,
        aes(
            x = (
                xmin + xmax
            ) / 2,
            y = 1.012,
            label = "optix"
        ),
        inherit.aes = FALSE,
        fontface = "italic",
        colour = "#277DA1",
        size = 2.8
    ) +

    facet_grid(
        module ~ .
    ) +

    scale_colour_manual(
        values = module_colours,
        guide = "none"
    ) +

    scale_x_continuous(
        name = NULL,
        limits = x_limits,
        breaks = seq(
            1.06,
            1.24,
            by = 0.02
        ),
        labels = function(x) {
            sprintf(
                "%.2f",
                x
            )
        },
        expand = expansion(
            mult = c(
                0,
                0
            )
        )
    ) +

    scale_y_continuous(
        name = "Balanced accuracy",
        limits = c(
            0.47,
            1.025
        ),
        breaks = c(
            0.5,
            0.6,
            0.7,
            0.8,
            0.9,
            1.0
        )
    ) +

    labs(
        title = paste(
            "Comparative pattern association across optix"
        ),
        subtitle = paste(
            "Triangles are independently selected comparative SNPs",
            "with minimum BA >= 0.80.",
            "TN_Left and TN_Right were excluded from selection."
        )
    ) +

    theme_minimal(
        base_size = 10
    ) +

    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(
            colour = "#ECECEC",
            linewidth = 0.25
        ),
        strip.text.y = element_text(
            face = "bold",
            size = 10
        ),
        strip.background = element_rect(
            fill = "#F2F2F2",
            colour = NA
        ),
        plot.title = element_text(
            face = "bold",
            size = 13
        ),
        plot.subtitle = element_text(
            size = 9,
            colour = "#444444"
        ),
        axis.title.y = element_text(
            face = "bold"
        ),
        plot.margin = margin(
            8,
            12,
            3,
            8
        )
    )


focal_plot <- ggplot() +

    # KRISTA_PLOT_WINDOWS_LAYER
    geom_rect(
        data = krista_regions,
        aes(
            xmin = xmin_mb,
            xmax = xmax_mb,
            ymin = -Inf,
            ymax = Inf
        ),
        inherit.aes = FALSE,
        fill = "#6F6F6F",
        colour = "#505050",
        alpha = 0.13,
        linewidth = 0.28
    ) +

    geom_rect(
        data = gene_annotation,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = -Inf,
            ymax = Inf
        ),
        inherit.aes = FALSE,
        fill = "#56B4E9",
        alpha = 0.10
    ) +

    geom_hline(
        yintercept = 0.50,
        colour = side_colours[
            "TN_Right"
        ],
        linewidth = 0.30,
        linetype = "dashed",
        alpha = 0.52
    ) +

    geom_hline(
        yintercept = 1.00,
        colour = side_colours[
            "TN_Left"
        ],
        linewidth = 0.30,
        linetype = "dashed",
        alpha = 0.52
    ) +

    geom_segment(
        data = other_pairs,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = right_DR_dosage,
            yend = left_DR_dosage
        ),
        colour = "#A8A8A8",
        linewidth = 0.27,
        alpha = 0.40
    ) +

    geom_segment(
        data = expected_pairs,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = right_DR_dosage,
            yend = left_DR_dosage
        ),
        colour = "#009E73",
        linewidth = 0.40,
        alpha = 0.64
    ) +

    geom_point(
        data = other_points,
        aes(
            x = position_Mb,
            y = DR_dosage
        ),
        colour = "#A8A8A8",
        shape = 124,
        size = 1.9,
        alpha = 0.50
    ) +

    geom_point(
        data = expected_points,
        aes(
            x = position_Mb,
            y = DR_dosage,
            colour = side
        ),
        shape = 124,
        size = 2.35,
        alpha = 0.92
    ) +

    facet_grid(
        module ~ .
    ) +

    scale_colour_manual(
        name = "Pooled side",
        values = side_colours
    ) +

    scale_x_continuous(
        name = "Position on Hmr18 (Mb)",
        limits = x_limits,
        breaks = seq(
            1.06,
            1.24,
            by = 0.02
        ),
        labels = function(x) {
            sprintf(
                "%.2f",
                x
            )
        },
        expand = expansion(
            mult = c(
                0,
                0
            )
        )
    ) +

    scale_y_continuous(
        name = "DR-associated allele dosage",
        limits = c(
            -0.08,
            1.08
        ),
        breaks = c(
            0,
            0.5,
            1
        ),
        labels = c(
            "0\nB/B",
            "0.5\nDR/B",
            "1\nDR/DR"
        )
    ) +

    labs(
        title = paste(
            "Pooled focal genotypes at comparative-selected SNPs"
        ),
        subtitle = paste(
            "Only left-right differences are drawn.",
            "Green marks the expected left DR/DR and right DR/B state;",
            "other differences are gray."
        ),
        caption = paste(
            "Identical genotypes and QC failures are omitted visually",
            "but retained in the accompanying summary table."
        )
    ) +

    theme_minimal(
        base_size = 10
    ) +

    theme(
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_line(
            colour = "#ECECEC",
            linewidth = 0.25
        ),
        panel.grid.major.y = element_line(
            colour = "#E2E2E2",
            linewidth = 0.30
        ),
        strip.text.y = element_text(
            face = "bold",
            size = 10
        ),
        strip.background = element_rect(
            fill = "#F2F2F2",
            colour = NA
        ),
        plot.title = element_text(
            face = "bold",
            size = 11
        ),
        plot.subtitle = element_text(
            size = 8.5,
            colour = "#444444"
        ),
        plot.caption = element_text(
            size = 7.7,
            colour = "#555555",
            hjust = 0
        ),
        legend.position = "top",
        legend.justification = "left",
        legend.title = element_text(
            face = "bold"
        ),
        axis.title = element_text(
            face = "bold"
        ),
        axis.text.y = element_text(
            size = 8
        ),
        plot.margin = margin(
            3,
            12,
            8,
            8
        )
    )


ba_grob <- ggplotGrob(
    ba_plot
)

focal_grob <- ggplotGrob(
    focal_plot
)

common_widths <- unit.pmax(
    ba_grob$widths,
    focal_grob$widths
)

ba_grob$widths <- common_widths
focal_grob$widths <- common_widths


draw_combined <- function() {
    grid.newpage()

    pushViewport(
        viewport(
            layout = grid.layout(
                nrow = 2,
                ncol = 1,
                heights = unit(
                    c(
                        3,
                        2.15
                    ),
                    "null"
                )
            )
        )
    )

    pushViewport(
        viewport(
            layout.pos.row = 1,
            layout.pos.col = 1
        )
    )

    grid.draw(
        ba_grob
    )

    popViewport()

    pushViewport(
        viewport(
            layout.pos.row = 2,
            layout.pos.col = 1
        )
    )

    grid.draw(
        focal_grob
    )

    popViewport(2)
}


cairo_pdf(
    filename = pdf_file,
    width = 11.5,
    height = 11.0
)

draw_combined()

dev.off()


png(
    filename = png_file,
    width = 3450,
    height = 3300,
    res = 300,
    type = "cairo"
)

draw_combined()

dev.off()


cat(
    "SIMPLE OPTIX COMPARATIVE AND FOCAL PATTERN FIGURE\n\n"
)

print(
    summary_rows,
    row.names = FALSE
)

cat(
    "\nUnique physical SNPs with expected focal pattern: ",
    nrow(unique_expected),
    "\n",
    sep = ""
)

cat(
    "Unique physical SNPs with any left-right difference: ",
    nrow(unique_differences),
    "\n",
    sep = ""
)

cat(
    "\nPDF:\n",
    pdf_file,
    "\n\nPNG:\n",
    png_file,
    "\n\nSite table:\n",
    site_file,
    "\n\nSummary:\n",
    summary_file,
    "\n\nBA background manifest:\n",
    background_manifest_file,
    "\n",
    sep = ""
)
