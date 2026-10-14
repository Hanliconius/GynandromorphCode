#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(scales)
    library(grid)
})

project <- "/scratch/martinlab/joe/GynandroCalling"

input_file <- file.path(
    project,
    "analysis/original_502/checkpoints",
    "optix_fullVCF_focal_and_WGS_control_AF.long.tsv.gz"
)

figure_dir <- file.path(
    project,
    "analysis/figures",
    "optix_fullVCF_AF_histograms"
)

pdf_file <- file.path(
    figure_dir,
    "optix_fullVCF_focal_and_WGS_control_AF_histograms.pdf"
)

png_file <- file.path(
    figure_dir,
    "optix_fullVCF_focal_and_WGS_control_AF_histograms.png"
)

het_file <- file.path(
    project,
    "analysis/original_502/checkpoints",
    "optix_fullVCF_heterozygous_called_SNPs.tsv"
)

summary_file <- file.path(
    project,
    "analysis/original_502/checkpoints",
    "optix_fullVCF_AF_histogram_sample_summary.tsv"
)

dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

if (!file.exists(input_file)) {
    stop(
        paste(
            "Missing allele-fraction table:",
            input_file
        )
    )
}

d <- read.delim(
    gzfile(input_file),
    stringsAsFactors = FALSE,
    check.names = FALSE
)

d$alt_fraction <- suppressWarnings(
    as.numeric(d$alt_fraction)
)

d$DP <- suppressWarnings(
    as.numeric(d$DP)
)

d$GQ <- suppressWarnings(
    as.numeric(d$GQ)
)

d$AD_total <- suppressWarnings(
    as.numeric(d$AD_total)
)

d$sample_order <- as.integer(
    d$sample_order
)

d$QC_pass <- tolower(
    as.character(d$QC_pass)
) %in% c(
    "true",
    "t",
    "1",
    "yes"
)

d <- d[
    d$QC_pass &
    !is.na(d$alt_fraction),
    ,
    drop = FALSE
]

if (nrow(d) == 0) {
    stop("No QC-passing allele-fraction observations remain.")
}

sample_information <- unique(
    d[
        ,
        c(
            "sample",
            "sample_type",
            "sample_order",
            "phenotype",
            "taxon_name"
        )
    ]
)

sample_information <- sample_information[
    order(sample_information$sample_order),
    ,
    drop = FALSE
]

sample_information$panel_label <- ifelse(
    sample_information$sample == "TN_Left",
    "TN_Left | pooled left",
    ifelse(
        sample_information$sample == "TN_Right",
        "TN_Right | pooled right",
        paste0(
            sample_information$sample,
            " | WGS rank ",
            sample_information$sample_order - 2,
            " | R=",
            sample_information$phenotype
        )
    )
)

panel_lookup <- setNames(
    sample_information$panel_label,
    sample_information$sample
)

d$panel_label <- unname(
    panel_lookup[d$sample]
)

panel_levels <- sample_information$panel_label

d$panel_label <- factor(
    d$panel_label,
    levels = panel_levels
)

sample_colours <- setNames(
    c(
        "#2B8CBE",
        "#D95F0E",
        "#7B3294",
        "#008B8B",
        "#C28E00"
    ),
    panel_levels
)


# ------------------------------------------------------------
# Bin counts manually so all and heterozygous distributions
# use the same denominator within each sample.
# ------------------------------------------------------------

bin_width <- 0.02
bin_count <- as.integer(1 / bin_width)

d$bin_index <- floor(
    d$alt_fraction / bin_width
)

d$bin_index <- pmax(
    0,
    pmin(
        bin_count - 1,
        d$bin_index
    )
)

d$bin_midpoint <- (
    d$bin_index + 0.5
) * bin_width

sample_totals <- aggregate(
    x = list(
        total_QC_calls = d$sample
    ),
    by = list(
        panel_label = as.character(d$panel_label)
    ),
    FUN = length
)

all_counts <- aggregate(
    x = list(
        count = d$sample
    ),
    by = list(
        panel_label = as.character(d$panel_label),
        bin_index = d$bin_index,
        bin_midpoint = d$bin_midpoint
    ),
    FUN = length
)

all_counts <- merge(
    all_counts,
    sample_totals,
    by = "panel_label",
    all.x = TRUE
)

all_counts$fraction <- (
    all_counts$count /
    all_counts$total_QC_calls
)

all_counts$panel_label <- factor(
    all_counts$panel_label,
    levels = panel_levels
)

heterozygous <- d[
    d$GT_class == "0/1",
    ,
    drop = FALSE
]

het_counts <- aggregate(
    x = list(
        count = heterozygous$sample
    ),
    by = list(
        panel_label = as.character(
            heterozygous$panel_label
        ),
        bin_index = heterozygous$bin_index,
        bin_midpoint = heterozygous$bin_midpoint
    ),
    FUN = length
)

het_counts <- merge(
    het_counts,
    sample_totals,
    by = "panel_label",
    all.x = TRUE
)

het_counts$fraction <- (
    het_counts$count /
    het_counts$total_QC_calls
)

het_counts$panel_label <- factor(
    het_counts$panel_label,
    levels = panel_levels
)


# ------------------------------------------------------------
# Per-row heterozygote annotation
# ------------------------------------------------------------

sample_summary_rows <- lapply(
    panel_levels,
    function(panel_name) {

        x <- d[
            as.character(d$panel_label) == panel_name,
            ,
            drop = FALSE
        ]

        h <- x[
            x$GT_class == "0/1",
            ,
            drop = FALSE
        ]

        data.frame(
            panel_label = panel_name,
            sample = unique(x$sample)[1],
            QC_calls = nrow(x),
            heterozygous_calls = nrow(h),
            heterozygous_fraction = nrow(h) / nrow(x),
            median_heterozygous_AF = if (
                nrow(h) > 0
            ) {
                median(
                    h$alt_fraction,
                    na.rm = TRUE
                )
            } else {
                NA_real_
            },
            heterozygous_AF_Q25 = if (
                nrow(h) > 0
            ) {
                quantile(
                    h$alt_fraction,
                    0.25,
                    na.rm = TRUE
                )
            } else {
                NA_real_
            },
            heterozygous_AF_Q75 = if (
                nrow(h) > 0
            ) {
                quantile(
                    h$alt_fraction,
                    0.75,
                    na.rm = TRUE
                )
            } else {
                NA_real_
            },
            stringsAsFactors = FALSE
        )
    }
)

sample_summary <- do.call(
    rbind,
    sample_summary_rows
)

sample_summary$annotation <- paste0(
    "0/1 calls: ",
    format(
        sample_summary$heterozygous_calls,
        big.mark = ","
    ),
    " (",
    percent(
        sample_summary$heterozygous_fraction,
        accuracy = 0.1
    ),
    ")"
)

sample_summary$panel_label <- factor(
    sample_summary$panel_label,
    levels = panel_levels
)

write.table(
    sample_summary,
    file = summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

write.table(
    heterozygous[
        ,
        c(
            "sample",
            "panel_label",
            "chrom",
            "pos",
            "ref",
            "alt",
            "raw_GT",
            "AD",
            "DP",
            "GQ",
            "alt_fraction"
        )
    ],
    file = het_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

p <- ggplot() +

    # Gray background: every QC-passing SNP.
    geom_col(
        data = all_counts,
        aes(
            x = bin_midpoint,
            y = fraction
        ),
        width = bin_width * 0.94,
        fill = "#D4D4D4",
        colour = NA
    ) +

    # Coloured overlay: only SNPs called 0/1.
    geom_col(
        data = het_counts,
        aes(
            x = bin_midpoint,
            y = fraction,
            fill = panel_label
        ),
        width = bin_width * 0.70,
        colour = NA,
        alpha = 0.92
    ) +

    # Individual 0/1-called observations along the baseline.
    geom_rug(
        data = heterozygous,
        aes(
            x = alt_fraction,
            colour = panel_label
        ),
        sides = "b",
        length = unit(
            0.025,
            "npc"
        ),
        alpha = 0.12,
        linewidth = 0.20
    ) +

    geom_vline(
        xintercept = 0.5,
        colour = "#202020",
        linewidth = 0.55,
        linetype = "dashed"
    ) +

    geom_text(
        data = sample_summary,
        aes(
            x = 0.74,
            y = Inf,
            label = annotation
        ),
        inherit.aes = FALSE,
        hjust = 0,
        vjust = 1.35,
        size = 3,
        colour = "#303030"
    ) +

    facet_grid(
        rows = vars(panel_label),
        scales = "free_y",
        switch = "y"
    ) +

    scale_fill_manual(
        values = sample_colours,
        guide = "none"
    ) +

    scale_colour_manual(
        values = sample_colours,
        guide = "none"
    ) +

    scale_x_continuous(
        name = paste(
            "ALT-allele fraction",
            "(0 = REF/REF; 0.5 = REF/ALT; 1 = ALT/ALT)"
        ),
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
                0.005,
                0.005
            )
        )
    ) +

    scale_y_continuous(
        name = paste(
            "Fraction of QC-passing SNPs per bin",
            "(square-root scale)"
        ),
        trans = "identity",
        limits = c(0, 0.015),
        oob = scales::squish,
        labels = percent_format(
            accuracy = 0.1
        ),
        breaks = c(
            0,
            0.001,
            0.005,
            0.01,
            0.025,
            0.05,
            0.10,
            0.20
        ),
        expand = expansion(
            mult = c(
                0,
                0.08
            )
        )
    ) +

    labs(
        title = paste(
            "Allele-fraction distributions across the",
            "broad optix VCF"
        ),
        subtitle = paste(
            "Gray bars show all QC-passing SNPs;",
            "coloured overlays and rugs show SNPs called 0/1"
        ),
        caption = paste(
            "Calls require DP >= 2 and GQ >= 20.",
            "The square-root y-axis retains the large homozygous",
            "peaks while expanding the low-frequency interior bins.",
            "The dashed line marks AF = 0.5.",
            "The VCF contains sites variable somewhere in the",
            "502-sample panel rather than invariant genome positions."
        )
    ) +

    theme_minimal(
        base_size = 9
    ) +

    theme(
        panel.grid.minor = element_blank(),

        panel.grid.major.x = element_line(
            colour = "#E2E2E2",
            linewidth = 0.3
        ),

        strip.placement = "outside",

        strip.text.y.left = element_text(
            angle = 0,
            hjust = 1,
            face = "bold",
            size = 8
        ),

        strip.background = element_rect(
            fill = "#F3F3F3",
            colour = "#A0A0A0",
            linewidth = 0.35
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

ggsave(
    filename = pdf_file,
    plot = p,
    device = cairo_pdf,
    width = 12,
    height = 12,
    units = "in",
    limitsize = FALSE
)

ggsave(
    filename = png_file,
    plot = p,
    width = 12,
    height = 12,
    units = "in",
    dpi = 300,
    limitsize = FALSE
)

cat("FULL OPTIX AF HISTOGRAMS WITH HETEROZYGOTE OVERLAY\n\n")
print(sample_summary, row.names = FALSE)
cat("\nPDF:\n", pdf_file, "\n")
cat("\nPNG:\n", png_file, "\n")
cat("\nHeterozygous SNP list:\n", het_file, "\n")
cat("\nSample summary:\n", summary_file, "\n")
