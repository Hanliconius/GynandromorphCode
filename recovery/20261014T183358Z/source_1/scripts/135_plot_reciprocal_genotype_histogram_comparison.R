#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(scales)
})

project <- "/scratch/martinlab/joe/GynandroCalling"

input_file <- file.path(
    project,
    "analysis/original_502/reciprocal_L_R_mixture",
    "optix_reciprocal_L_R_genotype_pairs.tsv.gz"
)

expected_input_sha256 <- paste0(
    "8dd2ebecdcd039ca643e9b68c122407bf",
    "a8f98745df2aac1e00a61901e0aa9a8"
)

output_directory <- file.path(
    project,
    "analysis/figures/optix_reciprocal_L_R_mixture"
)

checkpoint_directory <- file.path(
    project,
    "analysis/original_502/reciprocal_L_R_mixture"
)

pdf_file <- file.path(
    output_directory,
    "optix_original_vs_reciprocal_GT_AF_histograms.pdf"
)

png_file <- file.path(
    output_directory,
    "optix_original_vs_reciprocal_GT_AF_histograms.png"
)

summary_file <- file.path(
    checkpoint_directory,
    "optix_original_vs_reciprocal_GT_AF_histograms.summary.tsv"
)

manifest_file <- file.path(
    checkpoint_directory,
    "optix_original_vs_reciprocal_GT_AF_histograms.manifest.tsv"
)

analysis_version <- "1.0.0"
posterior_cutoff <- 0.95
bin_width <- 0.02
y_cap <- 0.015

dir.create(
    output_directory,
    recursive = TRUE,
    showWarnings = FALSE
)

sha256_file <- function(path) {

    answer <- system2(
        "sha256sum",
        args = path,
        stdout = TRUE,
        stderr = TRUE
    )

    if (length(answer) != 1L) {
        stop("Could not calculate SHA256 for: ", path)
    }

    strsplit(
        answer,
        "[[:space:]]+"
    )[[1]][1]
}

canonical_gt <- function(x) {

    x <- gsub(
        "|",
        "/",
        as.character(x),
        fixed = TRUE
    )

    x[x == "1/0"] <- "0/1"

    x
}

if (!file.exists(input_file)) {
    stop("Missing input: ", input_file)
}

observed_input_sha256 <- sha256_file(
    input_file
)

if (!identical(
    observed_input_sha256,
    expected_input_sha256
)) {
    stop(
        "Input SHA256 mismatch. Expected ",
        expected_input_sha256,
        "; observed ",
        observed_input_sha256
    )
}

d <- read.delim(
    gzfile(input_file),
    stringsAsFactors = FALSE,
    check.names = FALSE
)

required_columns <- c(
    "chrom",
    "pos",
    "left_original_GT",
    "right_original_GT",
    "left_observed_AF",
    "right_observed_AF",
    "mixture_left_GT",
    "mixture_right_GT",
    "pair_posterior",
    "posterior95"
)

missing_columns <- setdiff(
    required_columns,
    names(d)
)

if (length(missing_columns)) {
    stop(
        "Missing columns: ",
        paste(missing_columns, collapse = ", ")
    )
}

posterior95 <- (
    tolower(
        as.character(d$posterior95)
    ) %in% c(
        "true",
        "t",
        "1"
    )
)

d <- d[
    posterior95
    & is.finite(d$pair_posterior)
    & d$pair_posterior >= posterior_cutoff,
    ,
    drop = FALSE
]

if (!nrow(d)) {
    stop("No posterior-0.95 sites were retained.")
}

d$left_original_GT <- canonical_gt(
    d$left_original_GT
)

d$right_original_GT <- canonical_gt(
    d$right_original_GT
)

d$mixture_left_GT <- canonical_gt(
    d$mixture_left_GT
)

d$mixture_right_GT <- canonical_gt(
    d$mixture_right_GT
)

d$left_changed <- (
    d$left_original_GT
    != d$mixture_left_GT
)

d$right_changed <- (
    d$right_original_GT
    != d$mixture_right_GT
)

long <- rbind(
    data.frame(
        side = "TN_Left",
        calling_model = "Original VCF genotype",
        observed_AF = d$left_observed_AF,
        genotype = d$left_original_GT
    ),
    data.frame(
        side = "TN_Left",
        calling_model = "Reciprocal-mixture genotype",
        observed_AF = d$left_observed_AF,
        genotype = d$mixture_left_GT
    ),
    data.frame(
        side = "TN_Right",
        calling_model = "Original VCF genotype",
        observed_AF = d$right_observed_AF,
        genotype = d$right_original_GT
    ),
    data.frame(
        side = "TN_Right",
        calling_model = "Reciprocal-mixture genotype",
        observed_AF = d$right_observed_AF,
        genotype = d$mixture_right_GT
    )
)

long <- long[
    is.finite(long$observed_AF)
    & long$observed_AF >= 0
    & long$observed_AF <= 1
    & long$genotype %in% c(
        "0/0",
        "0/1",
        "1/1"
    ),
    ,
    drop = FALSE
]

side_levels <- c(
    "TN_Left",
    "TN_Right"
)

model_levels <- c(
    "Original VCF genotype",
    "Reciprocal-mixture genotype"
)

genotype_levels <- c(
    "0/0",
    "0/1",
    "1/1"
)

long$side <- factor(
    long$side,
    levels = side_levels
)

long$calling_model <- factor(
    long$calling_model,
    levels = model_levels
)

long$genotype <- factor(
    long$genotype,
    levels = genotype_levels
)

number_of_bins <- as.integer(
    1 / bin_width
)

long$bin_index <- pmin(
    floor(
        long$observed_AF
        / bin_width
    ),
    number_of_bins - 1L
)

long$bin_midpoint <- (
    long$bin_index + 0.5
) * bin_width

panel_totals <- aggregate(
    observed_AF ~ side + calling_model,
    data = long,
    FUN = length
)

names(panel_totals)[3] <- "panel_total"

binned <- aggregate(
    observed_AF ~
        side
        + calling_model
        + genotype
        + bin_midpoint,
    data = long,
    FUN = length
)

names(binned)[5] <- "SNPs"

binned <- merge(
    binned,
    panel_totals,
    by = c(
        "side",
        "calling_model"
    ),
    all.x = TRUE,
    sort = FALSE
)

binned$side <- factor(
    binned$side,
    levels = side_levels
)

binned$calling_model <- factor(
    binned$calling_model,
    levels = model_levels
)

binned$genotype <- factor(
    binned$genotype,
    levels = genotype_levels
)

binned$fraction <- (
    binned$SNPs
    / binned$panel_total
)

annotations <- data.frame(
    side = factor(
        c(
            "TN_Left",
            "TN_Left",
            "TN_Right",
            "TN_Right"
        ),
        levels = side_levels
    ),
    calling_model = factor(
        c(
            "Original VCF genotype",
            "Reciprocal-mixture genotype",
            "Original VCF genotype",
            "Reciprocal-mixture genotype"
        ),
        levels = model_levels
    ),
    label = c(
        sprintf(
            "n = %s SNPs",
            comma(nrow(d))
        ),
        sprintf(
            "%s left calls changed",
            comma(sum(d$left_changed))
        ),
        sprintf(
            "n = %s SNPs",
            comma(nrow(d))
        ),
        sprintf(
            "%s right calls changed",
            comma(sum(d$right_changed))
        )
    )
)

genotype_colours <- c(
    "0/0" = "#2C9BCB",
    "0/1" = "#F28E2B",
    "1/1" = "#B61F72"
)

plot_object <- ggplot(
    binned,
    aes(
        x = bin_midpoint,
        y = fraction,
        fill = genotype
    )
) +
    geom_col(
        width = bin_width * 0.92,
        colour = NA
    ) +
    geom_vline(
        xintercept = 0.5,
        linetype = "dashed",
        linewidth = 0.45,
        colour = "#222222"
    ) +
    geom_text(
        data = annotations,
        aes(
            x = Inf,
            y = Inf,
            label = label
        ),
        inherit.aes = FALSE,
        hjust = 1.08,
        vjust = 1.35,
        size = 3.2,
        colour = "#333333"
    ) +
    facet_grid(
        rows = vars(side),
        cols = vars(calling_model)
    ) +
    scale_fill_manual(
        name = "Assigned genotype",
        values = genotype_colours,
        drop = FALSE
    ) +
    scale_x_continuous(
        name = "Observed ALT-allele fraction",
        limits = c(0, 1),
        breaks = seq(0, 1, 0.25),
        expand = expansion(
            mult = c(0, 0)
        )
    ) +
    scale_y_continuous(
        name = paste0(
            "Fraction of SNPs per ",
            bin_width,
            " AF bin"
        ),
        limits = c(0, y_cap),
        oob = scales::squish,
        breaks = c(
            0,
            0.0025,
            0.005,
            0.01,
            0.015
        ),
        labels = percent_format(
            accuracy = 0.1
        ),
        expand = expansion(
            mult = c(0, 0.03)
        )
    ) +
    labs(
        title = paste(
            "Original versus reciprocal-mixture genotype assignments",
            "across the broad optix interval"
        ),
        subtitle = paste(
            "The same posterior >= 0.95 SNPs and observed allele fractions",
            "are shown in every panel; only the genotype assignment changes."
        ),
        caption = paste(
            "Bars are stacked by assigned genotype.",
            "The y-axis is capped at 1.5% per bin;",
            "the AF=0 and AF=1 endpoint bars exceed this cap.",
            "Focal GQ and B/D/R phenotype expectations were not used."
        )
    ) +
    theme_classic(
        base_size = 11
    ) +
    theme(
        plot.title = element_text(
            face = "bold",
            size = 14
        ),
        plot.subtitle = element_text(
            size = 10,
            colour = "#444444",
            margin = margin(b = 9)
        ),
        plot.caption = element_text(
            size = 8,
            colour = "#555555",
            hjust = 0
        ),
        strip.background = element_rect(
            fill = "#F1F1F1",
            colour = "#999999",
            linewidth = 0.35
        ),
        strip.text = element_text(
            face = "bold",
            size = 10
        ),
        panel.grid.major.y = element_line(
            colour = "#E3E3E3",
            linewidth = 0.3
        ),
        panel.grid.major.x = element_line(
            colour = "#ECECEC",
            linewidth = 0.25
        ),
        legend.position = "top",
        legend.justification = "left",
        legend.title = element_text(
            face = "bold"
        ),
        plot.margin = margin(
            10,
            14,
            10,
            10
        )
    )

ggsave(
    filename = pdf_file,
    plot = plot_object,
    width = 12,
    height = 7.8,
    device = cairo_pdf
)

ggsave(
    filename = png_file,
    plot = plot_object,
    width = 12,
    height = 7.8,
    units = "in",
    dpi = 320,
    bg = "white"
)

summary_table <- rbind(
    data.frame(
        side = "TN_Left",
        posterior95_SNPs = nrow(d),
        genotype_calls_changed =
            sum(d$left_changed),
        changed_fraction =
            mean(d$left_changed)
    ),
    data.frame(
        side = "TN_Right",
        posterior95_SNPs = nrow(d),
        genotype_calls_changed =
            sum(d$right_changed),
        changed_fraction =
            mean(d$right_changed)
    )
)

write.table(
    summary_table,
    file = summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

arguments <- commandArgs(
    trailingOnly = FALSE
)

script_file <- sub(
    "^--file=",
    "",
    arguments[
        grepl(
            "^--file=",
            arguments
        )
    ][1]
)

manifest <- data.frame(
    field = c(
        "analysis_version",
        "script",
        "script_SHA256",
        "input",
        "input_SHA256",
        "R_version",
        "ggplot2_version",
        "scales_version",
        "posterior_threshold",
        "AF_bin_width",
        "y_axis_cap",
        "observed_AF_changed",
        "phenotype_expectation_used"
    ),
    value = c(
        analysis_version,
        normalizePath(script_file),
        sha256_file(script_file),
        normalizePath(input_file),
        observed_input_sha256,
        as.character(getRversion()),
        as.character(
            packageVersion("ggplot2")
        ),
        as.character(
            packageVersion("scales")
        ),
        posterior_cutoff,
        bin_width,
        y_cap,
        FALSE,
        FALSE
    ),
    stringsAsFactors = FALSE
)

write.table(
    manifest,
    file = manifest_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat(
    "ORIGINAL VERSUS RECIPROCAL-GENOTYPE HISTOGRAMS\n\n"
)

print(
    summary_table,
    row.names = FALSE
)

cat(
    "\nPDF:\n",
    pdf_file,
    "\n\nPNG:\n",
    png_file,
    "\n\nSummary:\n",
    summary_file,
    "\n\nManifest:\n",
    manifest_file,
    "\n"
)
