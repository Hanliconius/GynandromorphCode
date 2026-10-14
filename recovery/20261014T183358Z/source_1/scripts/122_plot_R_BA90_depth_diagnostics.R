#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(grid)
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

audit_file <- file.path(
    project,
    "analysis/original_502/checkpoints",
    "optix_R_BA90_focal_depth_audit.tsv"
)

summary_file <- file.path(
    project,
    "analysis/original_502/checkpoints",
    "optix_R_BA90_depth_diagnostic.summary.tsv"
)

focal_pdf <- file.path(
    figure_dir,
    "optix_R_BA90_focal_depth_and_pattern_AF.pdf"
)

heatmap_pdf <- file.path(
    figure_dir,
    "optix_R_BA90_all_sample_depth_heatmap.pdf"
)

combined_pdf <- file.path(
    figure_dir,
    "optix_R_BA90_depth_diagnostics.pdf"
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

# For Ray, ALT is the pattern-associated allele.
d$pattern_AF <- ifelse(
    !is.na(d$AD_total) & d$AD_total > 0,
    d$alt_depth / d$AD_total,
    NA_real_
)

site_information <- unique(
    d[
        ,
        c(
            "pos",
            "site_label",
            "minimum_BA_all_analyses",
            "focal_outcome"
        )
    ]
)

site_information <- site_information[
    order(site_information$pos),
    ,
    drop = FALSE
]

site_levels <- site_information$site_label

d$site_label <- factor(
    d$site_label,
    levels = site_levels
)

# ------------------------------------------------------------
# Focal depth and allele-fraction table
# ------------------------------------------------------------

focal <- d[
    d$sample %in% c(
        "TN_Left",
        "TN_Right"
    ),
    ,
    drop = FALSE
]

focal$side <- ifelse(
    focal$sample == "TN_Left",
    "Left",
    "Right"
)

focal$side <- factor(
    focal$side,
    levels = c(
        "Left",
        "Right"
    )
)

focal <- focal[
    order(
        focal$pos,
        focal$side
    ),
    ,
    drop = FALSE
]

write.table(
    focal[
        ,
        c(
            "chrom",
            "pos",
            "ref",
            "alt",
            "minimum_BA_all_analyses",
            "focal_outcome",
            "side",
            "raw_GT",
            "AD",
            "ref_depth",
            "alt_depth",
            "AD_total",
            "DP",
            "GQ",
            "QC_pass",
            "pattern_AF",
            "biological_GT"
        )
    ],
    file = audit_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

side_colours <- c(
    "Left" = "#2B8CBE",
    "Right" = "#D95F0E"
)

depth_plot <- ggplot(
    focal,
    aes(
        x = site_label,
        y = DP,
        colour = side,
        group = side
    )
) +

    geom_line(
        linewidth = 0.45,
        alpha = 0.72
    ) +

    geom_point(
        size = 1.35,
        alpha = 0.90
    ) +

    scale_colour_manual(
        name = "Pooled side",
        values = side_colours
    ) +

    scale_y_continuous(
        name = "Total read depth (DP)",
        trans = "sqrt",
        expand = expansion(
            mult = c(0.03, 0.08)
        )
    ) +

    scale_x_discrete(
        name = NULL,
        labels = NULL,
        breaks = NULL,
        expand = expansion(
            mult = c(0.002, 0.002)
        )
    ) +

    labs(
        title = "Pooled focal depth across Ray-associated BA90 SNPs",
        subtitle = paste(
            "These are the same SNP columns as the broad",
            "Ray BA >= 0.90 genotype plot."
        )
    ) +

    theme_minimal(
        base_size = 9
    ) +

    theme(
        panel.grid.minor = element_blank(),

        panel.grid.major.x = element_blank(),

        plot.title = element_text(
            face = "bold",
            size = 12
        ),

        plot.subtitle = element_text(
            size = 8.5,
            colour = "#444444"
        ),

        legend.position = "top",
        legend.justification = "left",

        plot.margin = margin(
            8,
            10,
            2,
            10
        )
    )

af_plot <- ggplot(
    focal,
    aes(
        x = site_label,
        y = pattern_AF,
        colour = side,
        group = side
    )
) +

    annotate(
        "rect",
        xmin = -Inf,
        xmax = Inf,
        ymin = 0.85,
        ymax = 1.00,
        fill = side_colours["Left"],
        alpha = 0.10
    ) +

    annotate(
        "rect",
        xmin = -Inf,
        xmax = Inf,
        ymin = 0.50,
        ymax = 0.65,
        fill = side_colours["Right"],
        alpha = 0.10
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
        linewidth = 0.45,
        alpha = 0.72
    ) +

    geom_point(
        aes(
            size = pmax(
                AD_total,
                1
            )
        ),
        alpha = 0.88
    ) +

    scale_colour_manual(
        name = "Pooled side",
        values = side_colours
    ) +

    scale_size_continuous(
        name = "AD depth",
        range = c(
            0.8,
            3.2
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
            mult = c(0.01, 0.02)
        )
    ) +

    scale_x_discrete(
        name = "Hmr18 position",
        expand = expansion(
            mult = c(0.002, 0.002)
        )
    ) +

    labs(
        caption = paste(
            "Shaded bands show the expected observed fractions",
            "for a true Left=P/P and Right=P/A state under up to",
            "30% reciprocal mixing: Left 0.85-1.00;",
            "Right 0.50-0.65.",
            "Point size represents total allelic depth."
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
            size = 5.5
        ),

        legend.position = "top",
        legend.justification = "left",

        plot.caption = element_text(
            size = 7.2,
            colour = "#444444",
            hjust = 0
        ),

        plot.margin = margin(
            2,
            10,
            8,
            10
        )
    )

# ------------------------------------------------------------
# All-sample depth heatmap
# ------------------------------------------------------------

d$depth_bin <- ifelse(
    is.na(d$DP),
    "Missing",
    as.character(
        cut(
            d$DP,
            breaks = c(
                -Inf,
                0,
                1,
                4,
                9,
                19,
                49,
                99,
                Inf
            ),
            labels = c(
                "0",
                "1",
                "2-4",
                "5-9",
                "10-19",
                "20-49",
                "50-99",
                "100+"
            ),
            right = TRUE
        )
    )
)

depth_levels <- c(
    "Missing",
    "0",
    "1",
    "2-4",
    "5-9",
    "10-19",
    "20-49",
    "50-99",
    "100+"
)

d$depth_bin <- factor(
    d$depth_bin,
    levels = depth_levels
)

depth_colours <- c(
    "Missing" = "#D9D9D9",
    "0" = "#F2F2F2",
    "1" = "#2D114B",
    "2-4" = "#57106E",
    "5-9" = "#8C2981",
    "10-19" = "#C43C78",
    "20-49" = "#ED6925",
    "50-99" = "#FCA636",
    "100+" = "#F6E65C"
)

focal_rows <- unique(
    d[
        d$sample_type == "focal",
        c(
            "sample",
            "sample_label"
        )
    ]
)

focal_rows$order <- match(
    focal_rows$sample,
    c(
        "TN_Left",
        "TN_Right"
    )
)

focal_rows <- focal_rows[
    order(focal_rows$order),
    ,
    drop = FALSE
]

comparative_rows <- unique(
    d[
        d$sample_type == "comparative",
        c(
            "sample",
            "sample_label",
            "phenotype",
            "taxon_name"
        )
    ]
)

comparative_rows$phenotype_order <- ifelse(
    comparative_rows$phenotype == "1",
    1,
    ifelse(
        comparative_rows$phenotype == "0",
        2,
        3
    )
)

comparative_rows <- comparative_rows[
    order(
        comparative_rows$phenotype_order,
        comparative_rows$taxon_name,
        comparative_rows$sample
    ),
    ,
    drop = FALSE
]

sample_order <- c(
    focal_rows$sample_label,
    comparative_rows$sample_label
)

d$sample_label <- factor(
    d$sample_label,
    levels = rev(sample_order)
)

d$row_group <- factor(
    d$row_group,
    levels = c(
        "Focal gynandromorph",
        "Ray present",
        "Ray absent",
        "Phenotype unavailable"
    )
)

heatmap_plot <- ggplot(
    d,
    aes(
        x = site_label,
        y = sample_label,
        fill = depth_bin
    )
) +

    geom_tile(
        width = 0.96,
        height = 0.96,
        colour = NA
    ) +

    facet_grid(
        rows = vars(row_group),
        scales = "free_y",
        space = "free_y",
        switch = "y"
    ) +

    scale_fill_manual(
        name = "Read depth (DP)",
        values = depth_colours,
        breaks = depth_levels,
        drop = FALSE
    ) +

    scale_x_discrete(
        position = "top",
        expand = expansion(
            mult = c(0.001, 0.001)
        )
    ) +

    labs(
        title = "Depth across Ray-associated BA90 SNPs",
        subtitle = paste(
            "Rows and SNP columns match the broad Ray",
            "BA >= 0.90 genotype plot."
        ),
        x = "Hmr18 position",
        y = NULL
    ) +

    theme_minimal(
        base_size = 8
    ) +

    theme(
        panel.grid = element_blank(),

        axis.text.x = element_text(
            angle = 90,
            hjust = 0,
            vjust = 0.5,
            size = 5.5
        ),

        axis.text.y = element_text(
            size = 2.8,
            colour = "#222222"
        ),

        strip.placement = "outside",

        strip.text.y.left = element_text(
            angle = 0,
            face = "bold",
            size = 8
        ),

        strip.background = element_rect(
            fill = "#F2F2F2",
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

        plot.margin = margin(
            10,
            12,
            10,
            12
        )
    )

# ------------------------------------------------------------
# Summaries
# ------------------------------------------------------------

left_depth <- focal[
    focal$side == "Left",
    c(
        "pos",
        "DP",
        "pattern_AF"
    )
]

names(left_depth)[2:3] <- c(
    "left_DP",
    "left_pattern_AF"
)

right_depth <- focal[
    focal$side == "Right",
    c(
        "pos",
        "DP",
        "pattern_AF"
    )
]

names(right_depth)[2:3] <- c(
    "right_DP",
    "right_pattern_AF"
)

paired <- merge(
    left_depth,
    right_depth,
    by = "pos",
    all = TRUE
)

depth_correlation <- suppressWarnings(
    cor(
        paired$left_DP,
        paired$right_DP,
        method = "spearman",
        use = "complete.obs"
    )
)

summary <- data.frame(
    statistic = c(
        "SNPs",
        "comparative_samples",
        "left_median_DP",
        "right_median_DP",
        "left_minimum_DP",
        "right_minimum_DP",
        "left_maximum_DP",
        "right_maximum_DP",
        "left_right_DP_Spearman",
        "left_median_pattern_AF",
        "right_median_pattern_AF"
    ),
    value = c(
        length(unique(d$pos)),
        length(unique(d$sample[d$sample_type == "comparative"])),
        median(paired$left_DP, na.rm = TRUE),
        median(paired$right_DP, na.rm = TRUE),
        min(paired$left_DP, na.rm = TRUE),
        min(paired$right_DP, na.rm = TRUE),
        max(paired$left_DP, na.rm = TRUE),
        max(paired$right_DP, na.rm = TRUE),
        depth_correlation,
        median(paired$left_pattern_AF, na.rm = TRUE),
        median(paired$right_pattern_AF, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
)

write.table(
    summary,
    file = summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# ------------------------------------------------------------
# Write figures
# ------------------------------------------------------------

site_count <- length(unique(d$pos))

plot_width <- max(
    15,
    min(
        30,
        5 + 0.16 * site_count
    )
)

cairo_pdf(
    filename = focal_pdf,
    width = plot_width,
    height = 10
)

grid.newpage()

pushViewport(
    viewport(
        layout = grid.layout(
            nrow = 2,
            ncol = 1,
            heights = unit(
                c(
                    0.42,
                    0.58
                ),
                "npc"
            )
        )
    )
)

print(
    depth_plot,
    vp = viewport(
        layout.pos.row = 1,
        layout.pos.col = 1
    ),
    newpage = FALSE
)

print(
    af_plot,
    vp = viewport(
        layout.pos.row = 2,
        layout.pos.col = 1
    ),
    newpage = FALSE
)

invisible(dev.off())

ggsave(
    filename = heatmap_pdf,
    plot = heatmap_plot,
    device = cairo_pdf,
    width = plot_width,
    height = 28,
    units = "in",
    limitsize = FALSE
)

cairo_pdf(
    filename = combined_pdf,
    width = plot_width,
    height = 28,
    onefile = TRUE
)

grid.newpage()

pushViewport(
    viewport(
        layout = grid.layout(
            nrow = 2,
            ncol = 1,
            heights = unit(
                c(
                    0.42,
                    0.58
                ),
                "npc"
            )
        )
    )
)

print(
    depth_plot,
    vp = viewport(
        layout.pos.row = 1,
        layout.pos.col = 1
    ),
    newpage = FALSE
)

print(
    af_plot,
    vp = viewport(
        layout.pos.row = 2,
        layout.pos.col = 1
    ),
    newpage = FALSE
)

grid.newpage()
print(heatmap_plot, newpage = FALSE)

invisible(dev.off())

cat("RAY BA90 DEPTH DIAGNOSTICS\n\n")
print(summary, row.names = FALSE)
cat("\nFocal depth/AF PDF:\n", focal_pdf, "\n")
cat("\nAll-sample depth heatmap:\n", heatmap_pdf, "\n")
cat("\nCombined two-page PDF:\n", combined_pdf, "\n")
cat("\nFocal audit:\n", audit_file, "\n")
cat("\nSummary:\n", summary_file, "\n")
