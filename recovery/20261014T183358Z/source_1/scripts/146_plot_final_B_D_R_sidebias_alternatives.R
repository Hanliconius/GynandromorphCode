#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)

project <- "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"

output_dir <- file.path(
    project,
    "analysis",
    "figures",
    "final_B_D_R_sidebias_alternatives"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

suppressPackageStartupMessages({
    library(ggplot2)
    library(grid)
})

as_bool <- function(x) {
    tolower(trimws(as.character(x))) %in%
        c("true", "t", "1", "yes", "y")
}

load_before_device <- function(path) {
    lines <- readLines(path, warn = FALSE)

    cut_point <- which(
        trimws(lines) == "cairo_pdf("
    )

    if (length(cut_point) != 1) {
        stop(
            paste(
                basename(path),
                "has",
                length(cut_point),
                "cairo_pdf() cut points"
            )
        )
    }

    environment <- new.env(
        parent = globalenv()
    )

    eval(
        parse(
            text = paste(
                lines[seq_len(cut_point[[1]] - 1)],
                collapse = "\n"
            )
        ),
        envir = environment
    )

    environment
}

dr_environment <- load_before_device(
    file.path(
        project,
        "scripts",
        "115b_plot_D_R_Wallbank_fixed_dominant_adjusted_LR.R"
    )
)

band_environment <- load_before_device(
    file.path(
        project,
        "scripts",
        "142_plot_B_western_postman_vs_western_Amazon_DR_adjusted_LR.R"
    )
)

mixture_file <- file.path(
    project,
    "analysis",
    "original_502",
    "reciprocal_L_R_mixture",
    "optix_reciprocal_L_R_genotype_pairs.posterior95.tsv"
)

if (!file.exists(mixture_file)) {
    stop(
        paste(
            "Missing mixture table:",
            mixture_file
        )
    )
}

mixture <- read.delim(
    mixture_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

mixture_required <- c(
    "chrom",
    "pos",
    "ref",
    "alt",
    "mixture_left_GT",
    "mixture_right_GT",
    "pair_posterior",
    "posterior95"
)

missing_mixture <- setdiff(
    mixture_required,
    names(mixture)
)

if (length(missing_mixture) > 0) {
    stop(
        paste(
            "Mixture table is missing:",
            paste(
                missing_mixture,
                collapse = ", "
            )
        )
    )
}

mixture$pos <- as.numeric(
    mixture$pos
)

mixture$site_key <- paste(
    mixture$chrom,
    mixture$pos,
    mixture$ref,
    mixture$alt,
    sep = ":"
)

if (anyDuplicated(mixture$site_key)) {
    stop(
        "Mixture table contains duplicate physical SNPs"
    )
}

mixture$pair_posterior <- as.numeric(
    mixture$pair_posterior
)

mixture$confident <- (
    as_bool(mixture$posterior95)
    & is.finite(mixture$pair_posterior)
    & mixture$pair_posterior >= 0.95
)

# This is the raw molecular state used by the Band exact-site filter.
# Raw labels avoid reversing the meaning between Band and Dennis/Ray.
mixture$exact_raw_11_01 <- (
    mixture$confident
    & mixture$mixture_left_GT == "1/1"
    & mixture$mixture_right_GT == "0/1"
)

make_association <- function(
    data,
    trait_name
) {
    required <- c(
        "chrom",
        "pos",
        "ref",
        "alt",
        "balanced_accuracy"
    )

    missing <- setdiff(
        required,
        names(data)
    )

    if (length(missing) > 0) {
        stop(
            paste(
                trait_name,
                "table is missing:",
                paste(missing, collapse = ", ")
            )
        )
    }

    if (
        all(
            c(
                "analysis_cohort",
                "scheme"
            ) %in% names(data)
        )
    ) {
        primary <- (
            as.character(data$analysis_cohort)
                == "WGS_plus_CS"
            & as.character(data$scheme)
                == "DP2_GQ20"
        )

        if (any(primary)) {
            data <- data[
                primary,
                ,
                drop = FALSE
            ]
        }
    }

    eligible <- if (
        "score_eligible_bool" %in% names(data)
    ) {
        as_bool(data$score_eligible_bool)
    } else if (
        "score_eligible" %in% names(data)
    ) {
        as_bool(data$score_eligible)
    } else {
        stop(
            paste(
                trait_name,
                "table lacks a score-eligibility column"
            )
        )
    }

    output <- data.frame(
        trait = trait_name,
        chrom = as.character(data$chrom),
        pos = as.numeric(data$pos),
        ref = as.character(data$ref),
        alt = as.character(data$alt),
        balanced_accuracy = as.numeric(
            data$balanced_accuracy
        ),
        eligible = eligible,
        stringsAsFactors = FALSE
    )

    output$site_key <- paste(
        output$chrom,
        output$pos,
        output$ref,
        output$alt,
        sep = ":"
    )

    if (anyDuplicated(output$site_key)) {
        stop(
            paste(
                trait_name,
                "has duplicate physical SNPs after primary filtering"
            )
        )
    }

    mixture_index <- match(
        output$site_key,
        mixture$site_key
    )

    output$focal_confident <- (
        !is.na(mixture_index)
        & mixture$confident[mixture_index]
    )

    output$focal_confident[
        is.na(output$focal_confident)
    ] <- FALSE

    output$focal_exact <- (
        !is.na(mixture_index)
        & mixture$exact_raw_11_01[mixture_index]
    )

    output$focal_exact[
        is.na(output$focal_exact)
    ] <- FALSE

    output$position_Mb <- (
        output$pos / 1e6
    )

    output$high_BA <- (
        output$eligible
        & is.finite(output$balanced_accuracy)
        & output$balanced_accuracy >= 0.90
    )

    output$perfect_BA <- (
        output$eligible
        & is.finite(output$balanced_accuracy)
        & abs(
            output$balanced_accuracy - 1
        ) < 1e-12
    )

    output
}

dennis <- make_association(
    dr_environment$dennis,
    "Dennis"
)

ray <- make_association(
    dr_environment$ray,
    "Ray"
)

band <- make_association(
    band_environment$d,
    "Band"
)

association <- rbind(
    dennis,
    ray,
    band
)

association$trait <- factor(
    association$trait,
    levels = c(
        "Dennis",
        "Ray",
        "Band"
    )
)

x_start <- 1060000
x_end <- 1240000

x_limits <- c(
    x_start,
    x_end
) / 1e6

optix <- data.frame(
    xmin = 1088602 / 1e6,
    xmax = 1097846 / 1e6
)

region_data <- data.frame(
    xmin = numeric(),
    xmax = numeric()
)

if (
    exists(
        "plot_intervals",
        envir = dr_environment,
        inherits = FALSE
    )
) {
    intervals <- dr_environment$plot_intervals

    if (
        all(
            c("xmin", "xmax")
            %in% names(intervals)
        )
    ) {
        region_data <- unique(
            intervals[
                ,
                c("xmin", "xmax"),
                drop = FALSE
            ]
        )
    }
}

if (
    exists(
        "band_interval",
        envir = band_environment,
        inherits = FALSE
    )
) {
    interval <- band_environment$band_interval

    if (
        all(
            c("xmin", "xmax")
            %in% names(interval)
        )
    ) {
        region_data <- unique(
            rbind(
                region_data,
                interval[
                    ,
                    c("xmin", "xmax"),
                    drop = FALSE
                ]
            )
        )
    }
}

trait_colours <- c(
    Dennis = "#E69F00",
    Ray = "#CC79A7",
    Band = "#0072B2"
)

base_ba_plot <- function(data) {
    ggplot() +

        geom_rect(
            data = optix,
            aes(
                xmin = xmin,
                xmax = xmax,
                ymin = -Inf,
                ymax = Inf
            ),
            inherit.aes = FALSE,
            fill = "#56B4E9",
            alpha = 0.12
        ) +

        geom_rect(
            data = region_data,
            aes(
                xmin = xmin,
                xmax = xmax,
                ymin = -Inf,
                ymax = Inf
            ),
            inherit.aes = FALSE,
            fill = "#8C8C8C",
            alpha = 0.10
        ) +

        geom_point(
            data = data[
                is.finite(data$balanced_accuracy)
                & !data$eligible,
                ,
                drop = FALSE
            ],
            aes(
                x = position_Mb,
                y = balanced_accuracy
            ),
            colour = "#BDBDBD",
            alpha = 0.24,
            size = 0.34
        ) +

        geom_point(
            data = data[
                is.finite(data$balanced_accuracy)
                & data$eligible,
                ,
                drop = FALSE
            ],
            aes(
                x = position_Mb,
                y = balanced_accuracy,
                colour = trait
            ),
            alpha = 0.72,
            size = 0.42
        ) +

        geom_hline(
            yintercept = 0.90,
            linewidth = 0.30,
            linetype = "dotted",
            colour = "#555555"
        ) +

        geom_hline(
            yintercept = 0.80,
            linewidth = 0.35,
            linetype = "dashed",
            colour = "#222222"
        ) +

        annotate(
            "segment",
            x = optix$xmin,
            xend = optix$xmax,
            y = 1.018,
            yend = 1.018,
            colour = "#1764A3",
            linewidth = 2.2
        ) +

        annotate(
            "text",
            x = mean(
                c(
                    optix$xmin,
                    optix$xmax
                )
            ),
            y = 1.025,
            label = "optix",
            colour = "#1764A3",
            fontface = "italic",
            size = 3
        ) +

        facet_grid(
            trait ~ .,
            switch = "y"
        ) +

        scale_colour_manual(
            values = trait_colours,
            guide = "none"
        ) +

        scale_x_continuous(
            limits = x_limits,
            expand = expansion(
                mult = c(
                    0.004,
                    0.004
                )
            )
        ) +

        theme_minimal(
            base_size = 10
        ) +

        theme(
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_line(
                colour = "#E8E8E8",
                linewidth = 0.3
            ),
            panel.grid.major.y = element_line(
                colour = "#E2E2E2",
                linewidth = 0.3
            ),
            strip.placement = "outside",
            strip.text.y.left = element_text(
                face = "bold"
            ),
            axis.title.y = element_text(
                face = "bold"
            )
        )
}


# Alternative 1: BA panels with nested intersection rugs.

rug_all <- association[
    association$focal_exact,
    ,
    drop = FALSE
]

rug_all$rug_y <- 0.476

rug_high <- association[
    association$focal_exact
    & association$high_BA,
    ,
    drop = FALSE
]

rug_high$rug_y <- 0.458

rug_perfect <- association[
    association$focal_exact
    & association$perfect_BA,
    ,
    drop = FALSE
]

rug_perfect$rug_y <- 0.440

plot_1 <- base_ba_plot(
    association
) +

    geom_point(
        data = rug_all,
        aes(
            x = position_Mb,
            y = rug_y
        ),
        shape = 124,
        colour = "#777777",
        size = 2.0,
        alpha = 0.45
    ) +

    geom_point(
        data = rug_high,
        aes(
            x = position_Mb,
            y = rug_y,
            colour = trait
        ),
        shape = 124,
        size = 2.4,
        alpha = 0.90
    ) +

    geom_point(
        data = rug_perfect,
        aes(
            x = position_Mb,
            y = rug_y
        ),
        shape = 124,
        colour = "#111111",
        size = 2.7,
        alpha = 0.95
    ) +

    coord_cartesian(
        ylim = c(
            0.435,
            1.035
        ),
        clip = "off"
    ) +

    labs(
        title = paste(
            "Comparative BA with reciprocal-mixture",
            "focal intersections"
        ),
        subtitle = paste(
            paste(
                "Rugs use the same raw focal state as Band:",
                "TN_Left=1/1 and TN_Right=0/1, posterior >=0.95."
            ),
            paste(
                "Grey = all exact focal sites;",
                "colour = exact plus BA>=0.90;",
                "black = exact plus BA=1."
            )
        ),
        x = "Position on Hmr18 (Mb)",
        y = "Balanced accuracy",
        caption = paste(
            "Comparative scoring excludes the focal samples.",
            "Raw genotype labels avoid the reversed biological",
            "allele orientation between Band and Dennis/Ray."
        )
    ) +

    theme(
        plot.title = element_text(
            face = "bold",
            size = 13
        ),
        plot.subtitle = element_text(
            size = 9,
            colour = "#444444"
        ),
        plot.caption = element_text(
            size = 8,
            hjust = 0,
            colour = "#555555"
        )
    )


# Alternative 2: 10-kb binned prevalence and overlap.

bin_width <- 10000

bin_starts <- seq(
    x_start,
    x_end - bin_width,
    by = bin_width
)

bin_template <- data.frame(
    bin_start = bin_starts,
    bin_end = bin_starts + bin_width,
    bin_mid_Mb = (
        bin_starts
        + bin_width / 2
    ) / 1e6
)

mixture_region <- mixture[
    mixture$chrom == "Hmr18"
    & mixture$pos >= x_start
    & mixture$pos < x_end
    & mixture$confident,
    ,
    drop = FALSE
]

mixture_region$bin_start <- (
    floor(
        (
            mixture_region$pos
            - x_start
        ) / bin_width
    ) * bin_width
    + x_start
)

focal_bins <- bin_template

focal_bins$n_confident <- vapply(
    bin_starts,
    function(bin) {
        sum(
            mixture_region$bin_start
            == bin
        )
    },
    integer(1)
)

focal_bins$n_exact <- vapply(
    bin_starts,
    function(bin) {
        sum(
            mixture_region$bin_start
                == bin
            & mixture_region$exact_raw_11_01
        )
    },
    integer(1)
)

focal_bins$exact_fraction <- ifelse(
    focal_bins$n_confident > 0,
    focal_bins$n_exact
        / focal_bins$n_confident,
    NA_real_
)

make_binned_overlap <- function(
    data,
    trait_name
) {
    selected <- data[
        data$eligible
        & data$focal_confident
        & data$pos >= x_start
        & data$pos < x_end,
        ,
        drop = FALSE
    ]

    selected$bin_start <- (
        floor(
            (
                selected$pos
                - x_start
            ) / bin_width
        ) * bin_width
        + x_start
    )

    result <- bin_template
    result$trait <- trait_name

    result$n_confident_eligible <- vapply(
        bin_starts,
        function(bin) {
            sum(
                selected$bin_start
                == bin
            )
        },
        integer(1)
    )

    result$n_high <- vapply(
        bin_starts,
        function(bin) {
            sum(
                selected$bin_start
                    == bin
                & selected$high_BA
            )
        },
        integer(1)
    )

    result$n_exact_high <- vapply(
        bin_starts,
        function(bin) {
            sum(
                selected$bin_start
                    == bin
                & selected$high_BA
                & selected$focal_exact
            )
        },
        integer(1)
    )

    result$high_exact_fraction <- ifelse(
        result$n_high > 0,
        result$n_exact_high
            / result$n_high,
        NA_real_
    )

    denominator <- sum(
        selected$high_BA
    )

    result$overall_high_exact_fraction <- if (
        denominator > 0
    ) {
        sum(
            selected$high_BA
            & selected$focal_exact
        ) / denominator
    } else {
        NA_real_
    }

    result
}

binned_overlap <- rbind(
    make_binned_overlap(
        dennis,
        "Dennis"
    ),
    make_binned_overlap(
        ray,
        "Ray"
    ),
    make_binned_overlap(
        band,
        "Band"
    )
)

binned_overlap$trait <- factor(
    binned_overlap$trait,
    levels = c(
        "Dennis",
        "Ray",
        "Band"
    )
)

baseline <- unique(
    binned_overlap[
        ,
        c(
            "trait",
            "overall_high_exact_fraction"
        ),
        drop = FALSE
    ]
)

plot_2a <- ggplot(
    focal_bins,
    aes(
        x = bin_mid_Mb,
        y = exact_fraction
    )
) +

    geom_rect(
        data = optix,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = -Inf,
            ymax = Inf
        ),
        inherit.aes = FALSE,
        fill = "#56B4E9",
        alpha = 0.12
    ) +

    geom_line(
        colour = "#009E73",
        linewidth = 0.65,
        na.rm = TRUE
    ) +

    geom_point(
        aes(
            size = n_confident
        ),
        colour = "#006B50",
        alpha = 0.82,
        na.rm = TRUE
    ) +

    scale_size_continuous(
        range = c(
            1.2,
            4.5
        ),
        name = "Posterior95 sites"
    ) +

    scale_x_continuous(
        limits = x_limits,
        expand = expansion(
            mult = c(
                0.004,
                0.004
            )
        )
    ) +

    scale_y_continuous(
        limits = c(
            0,
            NA
        )
    ) +

    labs(
        title = paste(
            "10-kb prevalence of the exact focal",
            "1/1 to 0/1 state"
        ),
        subtitle = paste(
            "Numerator: exact raw focal state.",
            "Denominator: every reciprocal-mixture",
            "posterior >=0.95 site."
        ),
        x = NULL,
        y = "Exact-state fraction"
    ) +

    theme_minimal(
        base_size = 10
    ) +

    theme(
        plot.title = element_text(
            face = "bold",
            size = 12
        ),
        panel.grid.minor = element_blank()
    )

plot_2b <- ggplot(
    binned_overlap,
    aes(
        x = bin_mid_Mb,
        y = high_exact_fraction
    )
) +

    geom_rect(
        data = optix,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = -Inf,
            ymax = Inf
        ),
        inherit.aes = FALSE,
        fill = "#56B4E9",
        alpha = 0.12
    ) +

    geom_hline(
        data = baseline,
        aes(
            yintercept =
                overall_high_exact_fraction,
            colour = trait
        ),
        linewidth = 0.45,
        linetype = "dashed",
        show.legend = FALSE
    ) +

    geom_line(
        aes(
            colour = trait
        ),
        linewidth = 0.55,
        na.rm = TRUE
    ) +

    geom_point(
        aes(
            colour = trait,
            size = n_high
        ),
        alpha = 0.85,
        na.rm = TRUE
    ) +

    facet_grid(
        trait ~ .,
        switch = "y"
    ) +

    scale_colour_manual(
        values = trait_colours,
        guide = "none"
    ) +

    scale_size_continuous(
        range = c(
            1.0,
            4.5
        ),
        name = "BA>=0.90 sites"
    ) +

    scale_x_continuous(
        limits = x_limits,
        expand = expansion(
            mult = c(
                0.004,
                0.004
            )
        )
    ) +

    scale_y_continuous(
        limits = c(
            0,
            1
        )
    ) +

    labs(
        title = "Exact focal state among high-BA sites",
        subtitle = paste(
            "Points are 10-kb bins;",
            "dashed lines are interval-wide fractions.",
            "Only high-BA sites with a posterior95 focal call",
            "enter the denominator."
        ),
        x = "Position on Hmr18 (Mb)",
        y = "Exact fraction of BA>=0.90 sites"
    ) +

    theme_minimal(
        base_size = 10
    ) +

    theme(
        plot.title = element_text(
            face = "bold",
            size = 12
        ),
        panel.grid.minor = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_text(
            face = "bold"
        )
    )


# Alternative 4: aligned BA panels and one shared focal track.

plot_4a <- base_ba_plot(
    association
) +

    coord_cartesian(
        ylim = c(
            0.49,
            1.035
        ),
        clip = "off"
    ) +

    labs(
        title = paste(
            "Dennis, Ray, and Band",
            "comparative association"
        ),
        subtitle = paste(
            "All BA calculations use comparative butterflies only;",
            "the shared focal track is joined afterward."
        ),
        x = NULL,
        y = "Balanced accuracy"
    ) +

    theme(
        plot.title = element_text(
            face = "bold",
            size = 13
        ),
        plot.subtitle = element_text(
            size = 9,
            colour = "#444444"
        )
    )

exact_mixture <- mixture_region[
    mixture_region$exact_raw_11_01,
    ,
    drop = FALSE
]

exact_mixture$position_Mb <- (
    exact_mixture$pos / 1e6
)

focal_long <- rbind(
    data.frame(
        position_Mb =
            exact_mixture$position_Mb,
        side = "TN_Left",
        y = 1,
        genotype = "1/1",
        stringsAsFactors = FALSE
    ),
    data.frame(
        position_Mb =
            exact_mixture$position_Mb,
        side = "TN_Right",
        y = 0,
        genotype = "0/1",
        stringsAsFactors = FALSE
    )
)

overlap <- association[
    association$focal_exact
    & association$high_BA,
    ,
    drop = FALSE
]

overlap$rug_y <- c(
    Dennis = -0.24,
    Ray = -0.40,
    Band = -0.56
)[
    as.character(overlap$trait)
]

plot_4b <- ggplot() +

    geom_rect(
        data = optix,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = -Inf,
            ymax = Inf
        ),
        inherit.aes = FALSE,
        fill = "#56B4E9",
        alpha = 0.12
    ) +

    geom_rect(
        data = region_data,
        aes(
            xmin = xmin,
            xmax = xmax,
            ymin = -Inf,
            ymax = Inf
        ),
        inherit.aes = FALSE,
        fill = "#8C8C8C",
        alpha = 0.10
    ) +

    geom_segment(
        data = exact_mixture,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = 0,
            yend = 1
        ),
        colour = "#009E73",
        alpha = 0.30,
        linewidth = 0.28
    ) +

    geom_point(
        data = focal_long,
        aes(
            x = position_Mb,
            y = y,
            colour = genotype
        ),
        shape = 124,
        size = 2.1,
        alpha = 0.80
    ) +

    geom_point(
        data = overlap,
        aes(
            x = position_Mb,
            y = rug_y,
            colour = trait
        ),
        shape = 124,
        size = 2.4,
        alpha = 0.90
    ) +

    scale_colour_manual(
        values = c(
            "1/1" = "#0072B2",
            "0/1" = "#E69F00",
            trait_colours
        ),
        breaks = c(
            "1/1",
            "0/1",
            "Dennis",
            "Ray",
            "Band"
        ),
        labels = c(
            "Left 1/1",
            "Right 0/1",
            "Dennis BA>=0.90",
            "Ray BA>=0.90",
            "Band BA>=0.90"
        ),
        name = NULL
    ) +

    scale_x_continuous(
        limits = x_limits,
        expand = expansion(
            mult = c(
                0.004,
                0.004
            )
        )
    ) +

    scale_y_continuous(
        breaks = c(
            -0.56,
            -0.40,
            -0.24,
            0,
            1
        ),
        labels = c(
            "Band overlap",
            "Ray overlap",
            "Dennis overlap",
            "TN_Right 0/1",
            "TN_Left 1/1"
        ),
        limits = c(
            -0.68,
            1.18
        )
    ) +

    labs(
        title = "Shared reciprocal-mixture focal track",
        subtitle = paste(
            paste(
                "Every connector is TN_Left=1/1 and",
                "TN_Right=0/1 with posterior >=0.95."
            ),
            paste(
                "Lower rugs mark exact sites that also have",
                "BA>=0.90 in each comparative analysis."
            )
        ),
        x = "Position on Hmr18 (Mb)",
        y = NULL,
        caption = paste(
            "Raw genotype labels are used so that",
            "allele-orientation differences between Band",
            "and Dennis/Ray do not alter the definition",
            "of focal side bias."
        )
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
            colour = "#444444"
        ),
        plot.caption = element_text(
            size = 8,
            hjust = 0,
            colour = "#555555"
        ),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(
            face = "bold"
        ),
        legend.position = "top",
        legend.justification = "left"
    )

draw_stack <- function(
    top_plot,
    bottom_plot,
    heights = c(1, 1)
) {
    top_grob <- ggplotGrob(
        top_plot
    )

    bottom_grob <- ggplotGrob(
        bottom_plot
    )

    common_widths <- unit.pmax(
        top_grob$widths,
        bottom_grob$widths
    )

    top_grob$widths <- common_widths
    bottom_grob$widths <- common_widths

    grid.newpage()

    pushViewport(
        viewport(
            layout = grid.layout(
                2,
                1,
                heights = unit(
                    heights,
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
        top_grob
    )

    popViewport()

    pushViewport(
        viewport(
            layout.pos.row = 2,
            layout.pos.col = 1
        )
    )

    grid.draw(
        bottom_grob
    )

    popViewport(2)
}

file_1 <- file.path(
    output_dir,
    "146a_B_D_R_intersection_rugs.pdf"
)

file_2 <- file.path(
    output_dir,
    "146b_B_D_R_binned_sidebias_enrichment.pdf"
)

file_4 <- file.path(
    output_dir,
    "146c_B_D_R_aligned_BA_shared_focal_track.pdf"
)

file_all <- file.path(
    output_dir,
    "146_B_D_R_sidebias_alternatives_three_pages.pdf"
)

cairo_pdf(
    file_1,
    width = 11.5,
    height = 9.2
)

print(plot_1)

invisible(
    dev.off()
)

cairo_pdf(
    file_2,
    width = 11.5,
    height = 10
)

draw_stack(
    plot_2a,
    plot_2b,
    c(1, 1.8)
)

invisible(
    dev.off()
)

cairo_pdf(
    file_4,
    width = 11.5,
    height = 10.5
)

draw_stack(
    plot_4a,
    plot_4b,
    c(2.2, 1.25)
)

invisible(
    dev.off()
)

cairo_pdf(
    file_all,
    width = 11.5,
    height = 10.5,
    onefile = TRUE
)

print(plot_1)

draw_stack(
    plot_2a,
    plot_2b,
    c(1, 1.8)
)

draw_stack(
    plot_4a,
    plot_4b,
    c(2.2, 1.25)
)

invisible(
    dev.off()
)

summary_table <- do.call(
    rbind,
    lapply(
        split(
            association,
            association$trait
        ),
        function(data) {
            data.frame(
                trait = as.character(
                    data$trait[[1]]
                ),
                scored_sites = nrow(data),
                eligible_sites = sum(
                    data$eligible
                ),
                focal_posterior95_sites = sum(
                    data$focal_confident
                ),
                exact_raw_11_01_sites = sum(
                    data$focal_exact
                ),
                BA_ge_0_90_sites = sum(
                    data$high_BA
                ),
                exact_and_BA_ge_0_90_sites = sum(
                    data$focal_exact
                    & data$high_BA
                ),
                exact_and_BA_1_sites = sum(
                    data$focal_exact
                    & data$perfect_BA
                ),
                stringsAsFactors = FALSE
            )
        }
    )
)

write.table(
    summary_table,
    file.path(
        output_dir,
        "146_B_D_R_sidebias_overlap_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    focal_bins,
    file.path(
        output_dir,
        "146_focal_exact_fraction_by_10kb.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    binned_overlap,
    file.path(
        output_dir,
        "146_trait_high_BA_exact_overlap_by_10kb.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat(
    "\nSIDE-BIAS ALTERNATIVES COMPLETE\n\n"
)

print(
    summary_table,
    row.names = FALSE
)

cat(
    "\nPDFs:\n",
    paste(
        c(
            file_1,
            file_2,
            file_4,
            file_all
        ),
        collapse = "\n"
    ),
    "\n",
    sep = ""
)
