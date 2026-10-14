#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
})

project <- "/scratch/martinlab/joe/GynandroCalling"

input_file <- file.path(
    project,
    "analysis/original_502/checkpoints",
    "optix_first5000bp_focal_and_WGS_controls.long.tsv"
)

selection_file <- file.path(
    project,
    "analysis/original_502/checkpoints",
    "optix_first5000bp_selection.tsv"
)

figure_dir <- file.path(
    project,
    "analysis/figures",
    "optix_first5000bp_focal_vs_high_depth_WGS_allele_balance"
)

pdf_file <- file.path(
    figure_dir,
    "optix_first5000bp_focal_vs_high_depth_WGS_allele_balance.pdf"
)

png_file <- file.path(
    figure_dir,
    "optix_first5000bp_focal_vs_high_depth_WGS_allele_balance.png"
)

dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

if (!file.exists(input_file)) {
    stop(
        paste(
            "Missing plot input:",
            input_file
        )
    )
}

d <- read.delim(
    input_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

selection <- read.delim(
    selection_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

d$pos <- as.integer(d$pos)
d$alt_fraction <- suppressWarnings(
    as.numeric(d$alt_fraction)
)
d$AD_total <- suppressWarnings(
    as.numeric(d$AD_total)
)

d$QC_pass <- tolower(
    as.character(d$QC_pass)
) %in% c(
    "true",
    "t",
    "1",
    "yes"
)

sample_information <- unique(
    d[
        ,
        c(
            "sample",
            "display_name",
            "sample_order"
        )
    ]
)

sample_information$sample_order <- as.integer(
    sample_information$sample_order
)

sample_information <- sample_information[
    order(sample_information$sample_order),
    ,
    drop = FALSE
]

d$display_name <- factor(
    d$display_name,
    levels = sample_information$display_name
)

sample_colours <- c(
    "#2B8CBE",
    "#D95F0E",
    "#7B3294",
    "#008B8B",
    "#C28E00"
)

names(sample_colours) <- sample_information$sample

d$sample <- factor(
    d$sample,
    levels = sample_information$sample
)

qc_data <- d[
    d$QC_pass & !is.na(d$alt_fraction),
    ,
    drop = FALSE
]

low_qc_data <- d[
    !d$QC_pass & !is.na(d$alt_fraction),
    ,
    drop = FALSE
]

genotype_bands <- data.frame(
    ymin = c(
        0.00,
        0.42,
        0.90
    ),
    ymax = c(
        0.10,
        0.58,
        1.00
    ),
    fill_colour = c(
        "#B2188B",
        "#F39C34",
        "#2B9CC4"
    ),
    stringsAsFactors = FALSE
)

selected_region <- as.character(
    selection$region[1]
)

snp_count <- length(unique(d$pos))

p <- ggplot() +

    geom_rect(
        data = genotype_bands,
        aes(
            xmin = -Inf,
            xmax = Inf,
            ymin = ymin,
            ymax = ymax,
            fill = fill_colour
        ),
        inherit.aes = FALSE,
        alpha = 0.09,
        colour = NA
    ) +

    geom_hline(
        yintercept = c(
            0,
            0.5,
            1
        ),
        colour = "#777777",
        linewidth = 0.3,
        linetype = "dotted"
    ) +

    geom_line(
        data = qc_data,
        aes(
            x = pos / 1e6,
            y = alt_fraction,
            colour = sample,
            group = sample
        ),
        linewidth = 0.30,
        alpha = 0.30
    ) +

    geom_point(
        data = qc_data,
        aes(
            x = pos / 1e6,
            y = alt_fraction,
            colour = sample,
            size = pmax(
                AD_total,
                1
            )
        ),
        alpha = 0.86
    ) +

    geom_point(
        data = low_qc_data,
        aes(
            x = pos / 1e6,
            y = alt_fraction
        ),
        shape = 4,
        colour = "#A8A8A8",
        size = 0.9,
        stroke = 0.3,
        alpha = 0.55
    ) +

    facet_grid(
        rows = vars(display_name)
    ) +

    scale_colour_manual(
        values = sample_colours,
        guide = "none"
    ) +

    scale_fill_identity() +

    scale_size_continuous(
        name = "Allelic depth",
        range = c(
            0.7,
            3.0
        )
    ) +

    scale_y_continuous(
        name = paste(
            "ALT-allele fraction\n",
            "0 = REF/REF; 0.5 = REF/ALT; 1 = ALT/ALT"
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
                0.01,
                0.02
            )
        )
    ) +

    scale_x_continuous(
        name = "Hmr18 position (Mb)",
        expand = expansion(
            mult = c(
                0.005,
                0.005
            )
        )
    ) +

    labs(
        title = paste(
            "Allele balance in pooled focal and",
            "high-depth WGS samples"
        ),
        subtitle = paste0(
            selected_region,
            ": all ",
            snp_count,
            " biallelic SNPs in the first 5000-bp span of the broad VCF"
        ),
        caption = paste(
            "The three controls are the same highest-median-depth",
            "non-CS WGS samples used in the regional comparison.",
            "Every biallelic SNP in the first 5000-bp span of the",
            "broad VCF is displayed; no BA or focal-genotype filter",
            "was applied.",
            "Colored points pass DP >= 2 and GQ >= 20;",
            "gray crosses retain allele-depth information but fail QC."
        )
    ) +

    theme_minimal(
        base_size = 9
    ) +

    theme(
        panel.grid.minor = element_blank(),

        panel.grid.major.x = element_blank(),

        strip.text.y = element_text(
            angle = 0,
            hjust = 0,
            face = "bold",
            size = 7.8
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
            size = 7,
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
    width = 14,
    height = 14,
    units = "in",
    limitsize = FALSE
)

ggsave(
    filename = png_file,
    plot = p,
    width = 14,
    height = 14,
    units = "in",
    dpi = 300,
    limitsize = FALSE
)

cat("FOCAL VERSUS HIGH-DEPTH WGS ALLELE BALANCE\n\n")
cat("Selected region:", selected_region, "\n")
cat("Biallelic SNPs:", snp_count, "\n")
cat("Selection basis:", selection$selection_basis[1], "\n")
cat("\nPDF:\n", pdf_file, "\n")
cat("\nPNG:\n", png_file, "\n")
