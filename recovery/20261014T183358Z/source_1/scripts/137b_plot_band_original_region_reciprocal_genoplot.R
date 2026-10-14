#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(grid)
})

project <- "/scratch/martinlab/joe/GynandroCalling"

input_file <- file.path(
    project,
    "analysis/original_502/checkpoints",
    "optix_band_original_region_reciprocal_genoplot.long.tsv.gz"
)

figure_dir <- file.path(
    project,
    "analysis/figures",
    "optix_band_original_region_reciprocal_genoplot"
)

pdf_file <- file.path(
    figure_dir,
    paste0(
        "optix_band_original_region_genotypes_",
        "with_reciprocal_and_piano_keys.pdf"
    )
)

summary_file <- file.path(
    project,
    "analysis/original_502/checkpoints",
    "optix_band_original_region_reciprocal_genoplot.plot_summary.tsv"
)

if (!file.exists(input_file)) {
    stop(
        "Missing plot input: ",
        input_file
    )
}

dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

d <- read.delim(
    gzfile(input_file),
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
)

d$pos <- as.integer(
    d$pos
)

d$display_order <- as.integer(
    d$display_order
)

d$changed_from_original <- tolower(
    as.character(
        d$changed_from_original
    )
) %in% c(
    "true",
    "t",
    "1",
    "yes"
)

sites <- unique(
    d[
        order(
            d$pos,
            d$ref,
            d$alt
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

if (nrow(sites) != 153) {
    stop(
        "Expected 153 Band SNPs; recovered ",
        nrow(sites)
    )
}

sites$site_key <- paste(
    sites$chrom,
    sites$pos,
    sites$ref,
    sites$alt,
    sep = ":"
)

d$site_key <- paste(
    d$chrom,
    d$pos,
    d$ref,
    d$alt,
    sep = ":"
)

d$site_key <- factor(
    d$site_key,
    levels = sites$site_key,
    ordered = TRUE
)

sample_order <- unique(
    d[
        order(
            d$display_order
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

d$plot_group <- factor(
    d$plot_group,
    levels = group_order,
    ordered = TRUE
)

d$sample_label <- factor(
    d$sample_label,
    levels = rev(
        sample_order$sample_label
    ),
    ordered = TRUE
)

genotype_levels <- c(
    "A/A",
    "P/A",
    "P/P",
    "N",
    "Uncertain"
)

d$display_GT[
    !d$display_GT
    %in% genotype_levels
] <- "N"

d$display_GT <- factor(
    d$display_GT,
    levels = genotype_levels,
    ordered = TRUE
)

tick_indices <- unique(
    c(
        1L,
        seq.int(
            5L,
            nrow(sites),
            by = 5L
        ),
        nrow(sites)
    )
)

tick_breaks <- sites$site_key[
    tick_indices
]

tick_labels <- format(
    sites$pos[
        tick_indices
    ],
    big.mark = ",",
    scientific = FALSE,
    trim = TRUE
)

changed <- d[
    d$call_source
        == "reciprocal mixture (posterior >= 0.95)"
    & d$changed_from_original,
]

genotype_colours <- c(
    "A/A" = "#2C9BCB",
    "P/A" = "#F28E2B",
    "P/P" = "#B61F72",
    "N" = "#F2F2F2",
    "Uncertain" = "#777777"
)

plot_height <- max(
    22,
    4
    + 0.105
    * nrow(sample_order)
)

plot_object <- ggplot(
    d,
    aes(
        x = site_key,
        y = sample_label,
        fill = display_GT
    )
) +

    geom_tile(
        linewidth = 0
    ) +

    geom_tile(
        data = changed,
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
        name = "Band reference-state genotype"
    ) +

    scale_x_discrete(
        breaks = tick_breaks,
        labels = tick_labels,
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
            "Genotypes across Krista's",
            "original Band interval"
        ),
        subtitle = paste(
            "All 153 biallelic SNPs in",
            "Hmr18:1,112,651-1,114,314 are shown;",
            "no BA threshold was used.",
            "The original 481-sample Band cohort is retained."
        ),
        x = "Genomic position on Hmr18 (bp)",
        y = NULL,
        caption = paste(
            "Piano-key and pooled focal samples were not used",
            "to select SNPs.",
            "TN_Left and TN_Right are shown as both original",
            "VCF and reciprocal-mixture calls.",
            "Black outlines mark posterior >= 0.95 calls",
            "that differ from the original genotype."
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
    plot = plot_object,
    device = cairo_pdf,
    width = 20,
    height = plot_height,
    units = "in",
    limitsize = FALSE
)

plot_summary <- data.frame(
    statistic = c(
        "SNPs",
        "display_rows",
        "changed_reciprocal_cells",
        "PDF_height_inches"
    ),
    value = c(
        nrow(sites),
        nrow(sample_order),
        nrow(changed),
        plot_height
    )
)

write.table(
    plot_summary,
    summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat(
    "KRISTA ORIGINAL BAND-REGION GENOPLOT\n\n"
)

print(
    plot_summary,
    row.names = FALSE
)

cat(
    "\nPDF:\n",
    pdf_file,
    "\n"
)
