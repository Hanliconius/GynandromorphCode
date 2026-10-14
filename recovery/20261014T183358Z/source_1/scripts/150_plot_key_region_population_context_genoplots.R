#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(ggplot2)
    library(grid)
})

project <- "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"

output_dir <- file.path(
    project,
    "analysis",
    "figures",
    "key_region_population_context_genoplots"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

log_prefix <- "KEY-REGION POPULATION-CONTEXT GENOPLOTS"

# ============================================================
# USER-LEVEL SETTINGS
# ============================================================

key_regions <- c(
    "Band",
    "Ray1",
    "Ray2",
    "Dennis1",
    "Dennis2"
)

dp_min <- 2
gq_min <- 20

focal_samples <- c(
    "TN_Left",
    "TN_Right"
)

bwpk_pattern <- "(?i)\\bBWPK\\b|\\bBWPK\\b"

# Keep the older genotype-plot colour language:
# 0 associated copies  -> blue
# 1 copy / balanced    -> orange
# 2 associated copies  -> magenta
# missing              -> grey
zero_colour <- "#3B6FB6"
half_colour <- "#E69F00"
one_colour <- "#C51B7D"
missing_colour <- "#D0D0D0"

# This is the broad Optix VCF used in the older plotting work.
vcf_file <- file.path(
    project,
    "external",
    "krista_genotype_plots_original",
    "Genotype_Plots",
    "Optix_Hmr18-938195-1239421_biallelic.vcf"
)

# Region metadata: expected to contain region names, coordinates,
# and allele-orientation information used to polarize dosages.
region_metadata_file <- file.path(
    project,
    "data",
    "derived",
    "region_panel_metadata.tsv"
)

# Try a few plausible sample-metadata locations.
sample_metadata_candidates <- c(
    file.path(project, "data", "derived", "sample_plot_metadata.tsv"),
    file.path(project, "data", "derived", "sample_metadata.tsv"),
    file.path(project, "metadata", "sample_metadata.tsv"),
    file.path(project, "metadata", "original_502_sample_metadata.tsv"),
    file.path(project, "analysis", "original_502", "checkpoints", "sample_metadata.tsv"),
    file.path(project, "analysis", "original_502", "checkpoints", "cohort_sample_metadata.tsv")
)


# ============================================================
# HELPERS
# ============================================================

choose_existing_file <- function(
    paths
) {
    hits <- paths[file.exists(paths)]

    if (length(hits) == 0) {
        return(NA_character_)
    }

    hits[[1]]
}


find_first_column <- function(
    d,
    candidates,
    required = TRUE
) {
    hit <- candidates[candidates %in% names(d)]

    if (length(hit) > 0) {
        return(hit[[1]])
    }

    if (required) {
        stop(
            paste(
                "Could not find any of the expected columns:",
                paste(candidates, collapse = ", ")
            )
        )
    }

    NA_character_
}


canonical_gt <- function(
    x
) {
    x <- trimws(as.character(x))
    x <- gsub("\\|", "/", x)

    x[
        !(x %in% c(
            "0/0",
            "0/1",
            "1/0",
            "1/1"
        ))
    ] <- NA_character_

    x[x == "1/0"] <- "0/1"

    x
}


dose_from_gt <- function(
    gt,
    associated_is_alt
) {
    gt <- canonical_gt(gt)

    out <- rep(
        NA_real_,
        length(gt)
    )

    if (associated_is_alt) {
        out[gt == "0/0"] <- 0
        out[gt == "0/1"] <- 1
        out[gt == "1/1"] <- 2
    } else {
        out[gt == "0/0"] <- 2
        out[gt == "0/1"] <- 1
        out[gt == "1/1"] <- 0
    }

    out
}


normalise_boolish <- function(
    x
) {
    if (is.logical(x)) {
        return(x)
    }

    y <- tolower(trimws(as.character(x)))

    out <- rep(NA, length(y))

    out[y %in% c("true", "t", "1", "yes", "y")] <- TRUE
    out[y %in% c("false", "f", "0", "no", "n")] <- FALSE

    out
}


interpret_orientation <- function(
    value
) {
    if (length(value) == 0 || is.na(value)) {
        return(NA)
    }

    x <- tolower(trimws(as.character(value[[1]])))

    if (x %in% c("alt", "1", "true", "t", "yes", "y")) {
        return(TRUE)
    }

    if (x %in% c("ref", "0", "false", "f", "no", "n")) {
        return(FALSE)
    }

    NA
}


safe_system2_stdout <- function(
    command,
    args
) {
    out <- tryCatch(
        system2(
            command,
            args = args,
            stdout = TRUE,
            stderr = TRUE
        ),
        warning = function(w) {
            conditionMessage(w)
        },
        error = function(e) {
            stop(
                paste(
                    "Command failed:",
                    command,
                    paste(args, collapse = " "),
                    "\n",
                    conditionMessage(e)
                )
            )
        }
    )

    status <- attr(out, "status")

    if (!is.null(status) && status != 0) {
        stop(
            paste(
                "Command returned non-zero status:",
                command,
                paste(args, collapse = " "),
                "\n",
                paste(out, collapse = "\n")
            )
        )
    }

    out
}


query_region_vcf <- function(
    vcf,
    chrom,
    start,
    end
) {
    sample_names <- safe_system2_stdout(
        "bcftools",
        c(
            "query",
            "-l",
            vcf
        )
    )

    sample_names <- sample_names[nzchar(sample_names)]

    if (length(sample_names) == 0) {
        stop("No sample names recovered from VCF.")
    }

    region_string <- paste0(
        chrom,
        ":",
        start,
        "-",
        end
    )

    fmt <- paste0(
        "%CHROM\t%POS\t%REF\t%ALT[\t%GT:%DP:%GQ]\n"
    )

    lines <- safe_system2_stdout(
        "bcftools",
        c(
            "query",
            "-r",
            region_string,
            "-f",
            fmt,
            vcf
        )
    )

    lines <- lines[nzchar(lines)]

    if (length(lines) == 0) {
        return(
            list(
                sites = data.frame(),
                calls = data.frame(),
                sample_names = sample_names
            )
        )
    }

    split_lines <- strsplit(
        lines,
        "\t",
        fixed = TRUE
    )

    n_expected <- 4 + length(sample_names)

    bad_line <- which(
        vapply(
            split_lines,
            length,
            integer(1)
        ) != n_expected
    )

    if (length(bad_line) > 0) {
        stop(
            paste(
                "Malformed bcftools-query output in region",
                region_string,
                "- line",
                bad_line[[1]],
                "had",
                length(split_lines[[bad_line[[1]]]]),
                "fields but expected",
                n_expected
            )
        )
    }

    site_records <- vector(
        "list",
        length(split_lines)
    )

    call_records <- vector(
        "list",
        length(split_lines) * length(sample_names)
    )

    call_i <- 0L

    for (i in seq_along(split_lines)) {

        fields <- split_lines[[i]]

        chrom_i <- fields[[1]]
        pos_i <- as.numeric(fields[[2]])
        ref_i <- fields[[3]]
        alt_i <- fields[[4]]

        if (grepl(",", alt_i, fixed = TRUE)) {
            next
        }

        site_id <- paste(
            chrom_i,
            pos_i,
            ref_i,
            alt_i,
            sep = ":"
        )

        site_records[[i]] <- data.frame(
            site_id = site_id,
            chrom = chrom_i,
            pos = pos_i,
            ref = ref_i,
            alt = alt_i,
            stringsAsFactors = FALSE
        )

        sample_fields <- fields[-seq_len(4)]

        for (j in seq_along(sample_names)) {

            tok <- strsplit(
                sample_fields[[j]],
                ":",
                fixed = TRUE
            )[[1]]

            gt <- if (length(tok) >= 1) tok[[1]] else NA_character_
            dp <- if (length(tok) >= 2) tok[[2]] else NA_character_
            gq <- if (length(tok) >= 3) tok[[3]] else NA_character_

            call_i <- call_i + 1L

            call_records[[call_i]] <- data.frame(
                site_id = site_id,
                chrom = chrom_i,
                pos = pos_i,
                ref = ref_i,
                alt = alt_i,
                sample_id = sample_names[[j]],
                GT = gt,
                DP = suppressWarnings(as.numeric(dp)),
                GQ = suppressWarnings(as.numeric(gq)),
                stringsAsFactors = FALSE
            )
        }
    }

    sites <- do.call(
        rbind,
        site_records[!vapply(site_records, is.null, logical(1))]
    )

    calls <- do.call(
        rbind,
        call_records[!vapply(call_records, is.null, logical(1))]
    )

    if (is.null(sites) || nrow(sites) == 0) {
        return(
            list(
                sites = data.frame(),
                calls = data.frame(),
                sample_names = sample_names
            )
        )
    }

    sites <- sites[
        order(sites$pos),
        ,
        drop = FALSE
    ]

    calls <- calls[
        match(calls$site_id, sites$site_id),
        ,
        drop = FALSE
    ]

    list(
        sites = sites,
        calls = calls,
        sample_names = sample_names
    )
}


add_section_divider_rows <- function(
    plot_df,
    population_rows,
    individual_rows
) {
    if (length(population_rows) == 0 || length(individual_rows) == 0) {
        return(plot_df)
    }

    divider_label <- " "

    extra <- unique(
        plot_df[
            ,
            c("site_id", "site_index", "pos"),
            drop = FALSE
        ]
    )

    extra$display_row <- divider_label
    extra$row_type <- "divider"
    extra$plot_value <- NA_real_

    out <- rbind(
        plot_df,
        extra
    )

    desired_levels <- c(
        population_rows,
        divider_label,
        individual_rows
    )

    out$display_row <- factor(
        out$display_row,
        levels = rev(desired_levels)
    )

    out
}


make_region_plot <- function(
    plot_df,
    region_name,
    chrom,
    start,
    end
) {
    plot_df <- plot_df[
        order(
            plot_df$site_index
        ),
        ,
        drop = FALSE
    ]

    n_sites <- length(unique(plot_df$site_id))
    x_break_idx <- unique(
        round(
            seq(
                1,
                n_sites,
                length.out = min(8, n_sites)
            )
        )
    )

    site_ref <- unique(
        plot_df[
            ,
            c("site_index", "pos"),
            drop = FALSE
        ]
    )

    site_ref <- site_ref[
        order(site_ref$site_index),
        ,
        drop = FALSE
    ]

    x_break_pos <- x_break_idx
    x_break_lab <- sprintf(
        "%.4f",
        site_ref$pos[match(x_break_idx, site_ref$site_index)] / 1e6
    )

    population_labels <- levels(droplevels(
        plot_df$display_row
    ))

    divider_idx <- which(population_labels == " ")

    separator_y <- NULL
    if (length(divider_idx) == 1) {
        separator_y <- divider_idx + 0.5
    }

    p <- ggplot(
        plot_df,
        aes(
            x = site_index,
            y = display_row,
            fill = plot_value
        )
    ) +

        geom_tile(
            width = 0.95,
            height = 0.92,
            colour = NA
        ) +

        scale_fill_gradientn(
            colours = c(
                zero_colour,
                half_colour,
                one_colour
            ),
            values = c(
                0,
                0.5,
                1
            ),
            limits = c(0, 1),
            oob = scales::squish,
            na.value = missing_colour,
            breaks = c(
                0,
                0.5,
                1
            ),
            labels = c(
                "0 associated-allele copies",
                "1 copy / balanced",
                "2 associated-allele copies"
            ),
            name = "Associated-allele state"
        ) +

        scale_x_continuous(
            breaks = x_break_pos,
            labels = x_break_lab,
            expand = c(0.01, 0.01),
            name = paste0(
                chrom,
                " position (Mb)"
            )
        ) +

        labs(
            title = paste0(
                region_name,
                ": population genotype context plus focal and BWPK individuals"
            ),
            subtitle = paste(
                chrom,
                paste0(start, "-", end, "."),
                "Top rows: per-population mean associated-allele fraction among comparative butterflies.",
                paste0(
                    "Bottom rows: TN_Left, TN_Right, and BWPK individuals as QC-filtered genotypes (DP≥",
                    dp_min,
                    ", GQ≥",
                    gq_min,
                    ")."
                )
            ),
            caption = paste(
                "Columns are SNPs in genomic order.",
                "Population rows use mean associated-allele fraction;",
                "individual rows use exact QC-passing genotype state;",
                "missing/failed calls are grey."
            )
        ) +

        theme_minimal(
            base_size = 9
        ) +

        theme(
            plot.title = element_text(
                face = "bold",
                size = 11
            ),
            plot.subtitle = element_text(
                size = 8.5,
                colour = "#444444",
                margin = margin(b = 8)
            ),
            plot.caption = element_text(
                size = 7.6,
                colour = "#555555",
                hjust = 0,
                margin = margin(t = 8)
            ),
            panel.grid = element_blank(),
            axis.text.x = element_text(
                size = 7.5,
                angle = 0,
                vjust = 1
            ),
            axis.text.y = element_text(
                size = 8,
                colour = "#222222"
            ),
            axis.title.x = element_text(
                size = 8.8,
                margin = margin(t = 6)
            ),
            axis.title.y = element_blank(),
            legend.position = "top",
            legend.justification = "left",
            legend.key.width = unit(1.25, "cm"),
            legend.title = element_text(
                face = "bold",
                size = 8.5
            ),
            legend.text = element_text(
                size = 8
            ),
            plot.margin = margin(
                10,
                14,
                10,
                10
            )
        )

    if (!is.null(separator_y)) {
        p <- p +
            geom_hline(
                yintercept = separator_y,
                colour = "#6E6E6E",
                linewidth = 0.35
            )
    }

    p
}


# ============================================================
# INPUT CHECKS
# ============================================================

if (!file.exists(vcf_file)) {
    stop(
        paste("Missing VCF:", vcf_file)
    )
}

if (!file.exists(region_metadata_file)) {
    stop(
        paste(
            "Missing region metadata:",
            region_metadata_file
        )
    )
}

sample_metadata_file <- choose_existing_file(
    sample_metadata_candidates
)

if (is.na(sample_metadata_file)) {
    stop(
        paste(
            "Could not find a sample-metadata file.",
            "Checked:",
            paste(sample_metadata_candidates, collapse = "\n")
        )
    )
}

cat(
    "\n========================================\n",
    log_prefix,
    "\n========================================\n\n",
    sep = ""
)

cat(
    "VCF: ", vcf_file, "\n",
    "Region metadata: ", region_metadata_file, "\n",
    "Sample metadata: ", sample_metadata_file, "\n\n",
    sep = ""
)


# ============================================================
# LOAD REGION METADATA
# ============================================================

region_meta <- read.delim(
    region_metadata_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

region_name_col <- find_first_column(
    region_meta,
    c(
        "region",
        "panel",
        "module",
        "interval",
        "region_name"
    )
)

chrom_col <- find_first_column(
    region_meta,
    c(
        "chrom",
        "chromosome",
        "chr"
    )
)

start_col <- find_first_column(
    region_meta,
    c(
        "start",
        "region_start",
        "window_start",
        "plot_start",
        "xmin_bp"
    )
)

end_col <- find_first_column(
    region_meta,
    c(
        "end",
        "region_end",
        "window_end",
        "plot_end",
        "xmax_bp"
    )
)

orientation_col <- find_first_column(
    region_meta,
    c(
        "selected_orientation",
        "associated_allele",
        "pattern_allele",
        "polarized_allele",
        "associated_state",
        "alt_is_associated"
    )
)

region_meta$region_name <- as.character(region_meta[[region_name_col]])
region_meta$chrom_use <- as.character(region_meta[[chrom_col]])
region_meta$start_use <- suppressWarnings(as.numeric(region_meta[[start_col]]))
region_meta$end_use <- suppressWarnings(as.numeric(region_meta[[end_col]]))
region_meta$associated_is_alt <- vapply(
    region_meta[[orientation_col]],
    interpret_orientation,
    logical(1)
)

region_meta <- region_meta[
    region_meta$region_name %in% key_regions,
    ,
    drop = FALSE
]

missing_regions <- setdiff(
    key_regions,
    region_meta$region_name
)

if (length(missing_regions) > 0) {
    stop(
        paste(
            "These key regions were not found in region metadata:",
            paste(missing_regions, collapse = ", ")
        )
    )
}

if (any(!is.finite(region_meta$start_use) | !is.finite(region_meta$end_use))) {
    stop("Region metadata contains non-numeric start/end coordinates.")
}

if (any(is.na(region_meta$associated_is_alt))) {
    bad <- region_meta$region_name[is.na(region_meta$associated_is_alt)]
    stop(
        paste(
            "Could not interpret allele orientation for:",
            paste(bad, collapse = ", ")
        )
    )
}

region_meta <- region_meta[
    match(key_regions, region_meta$region_name),
    ,
    drop = FALSE
]


# ============================================================
# LOAD SAMPLE METADATA
# ============================================================

sample_meta <- read.delim(
    sample_metadata_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

sample_id_col <- find_first_column(
    sample_meta,
    c(
        "sample_id",
        "sample",
        "Sample",
        "ID",
        "id"
    )
)

population_col <- find_first_column(
    sample_meta,
    c(
        "plot_group",
        "population_plot_group",
        "population_group",
        "population",
        "group",
        "locality_group",
        "cohort"
    )
)

sample_meta$sample_id_use <- as.character(sample_meta[[sample_id_col]])
sample_meta$population_use <- as.character(sample_meta[[population_col]])

sample_meta <- sample_meta[
    !is.na(sample_meta$sample_id_use) &
    nzchar(sample_meta$sample_id_use),
    ,
    drop = FALSE
]

if (anyDuplicated(sample_meta$sample_id_use)) {
    dup <- unique(sample_meta$sample_id_use[duplicated(sample_meta$sample_id_use)])
    stop(
        paste(
            "Duplicate sample IDs in sample metadata:",
            paste(head(dup, 20), collapse = ", ")
        )
    )
}

char_cols <- names(sample_meta)[vapply(sample_meta, is.character, logical(1))]

bwpk_hits <- rep(FALSE, nrow(sample_meta))
for (cc in char_cols) {
    bwpk_hits <- bwpk_hits | grepl(
        bwpk_pattern,
        sample_meta[[cc]],
        perl = TRUE
    )
}

bwpk_samples_all <- unique(sample_meta$sample_id_use[bwpk_hits])

cat(
    "Detected BWPK samples in metadata: ",
    length(bwpk_samples_all),
    "\n",
    sep = ""
)

if (length(bwpk_samples_all) > 0) {
    cat(
        paste(
            "BWPK IDs:",
            paste(bwpk_samples_all, collapse = ", ")
        ),
        "\n\n"
    )
} else {
    cat("No BWPK samples detected by regex; proceeding with TN_Left/TN_Right only.\n\n")
}


# ============================================================
# MAIN LOOP
# ============================================================

plot_list <- list()
summary_rows <- list()

for (region_name in key_regions) {

    cat(
        "----------------------------------------\n",
        "Region: ", region_name, "\n",
        "----------------------------------------\n",
        sep = ""
    )

    rr <- region_meta[
        region_meta$region_name == region_name,
        ,
        drop = FALSE
    ]

    chrom <- rr$chrom_use[[1]]
    start <- rr$start_use[[1]]
    end <- rr$end_use[[1]]
    associated_is_alt <- rr$associated_is_alt[[1]]

    cat(
        "Coordinates: ",
        chrom,
        ":",
        start,
        "-",
        end,
        "\n",
        "Associated allele is ",
        ifelse(associated_is_alt, "ALT", "REF"),
        "\n",
        sep = ""
    )

    q <- query_region_vcf(
        vcf = vcf_file,
        chrom = chrom,
        start = start,
        end = end
    )

    if (nrow(q$sites) == 0 || nrow(q$calls) == 0) {
        cat("No SNPs recovered; skipping.\n\n")
        next
    }

    vcf_samples <- unique(q$calls$sample_id)

    focal_present <- focal_samples[focal_samples %in% vcf_samples]
    missing_focal <- setdiff(focal_samples, focal_present)

    if (length(missing_focal) > 0) {
        stop(
            paste(
                "Focal sample(s) absent from VCF:",
                paste(missing_focal, collapse = ", ")
            )
        )
    }

    bwpk_present <- bwpk_samples_all[bwpk_samples_all %in% vcf_samples]

    sample_meta_region <- sample_meta[
        sample_meta$sample_id_use %in% vcf_samples,
        ,
        drop = FALSE
    ]

    q$calls$GT_use <- canonical_gt(q$calls$GT)
    q$calls$pass_qc <- (
        !is.na(q$calls$GT_use) &
        is.finite(q$calls$DP) &
        is.finite(q$calls$GQ) &
        q$calls$DP >= dp_min &
        q$calls$GQ >= gq_min
    )

    q$calls$dose <- dose_from_gt(
        q$calls$GT_use,
        associated_is_alt = associated_is_alt
    )

    q$calls$plot_value <- ifelse(
        q$calls$pass_qc,
        q$calls$dose / 2,
        NA_real_
    )

    comparative_meta <- sample_meta_region[
        !(sample_meta_region$sample_id_use %in% c(focal_present, bwpk_present)),
        ,
        drop = FALSE
    ]

    if (nrow(comparative_meta) == 0) {
        stop(
            paste(
                "No comparative samples remain for",
                region_name,
                "after excluding TN and BWPK."
            )
        )
    }

    comparative_calls <- q$calls[
        q$calls$sample_id %in% comparative_meta$sample_id_use,
        ,
        drop = FALSE
    ]

    pop_map <- comparative_meta[
        ,
        c("sample_id_use", "population_use"),
        drop = FALSE
    ]

    comparative_calls$population_use <- pop_map$population_use[
        match(
            comparative_calls$sample_id,
            pop_map$sample_id_use
        )
    ]

    comparative_calls <- comparative_calls[
        !is.na(comparative_calls$population_use) &
        nzchar(comparative_calls$population_use),
        ,
        drop = FALSE
    ]

    comparative_qc <- comparative_calls[
        comparative_calls$pass_qc,
        ,
        drop = FALSE
    ]

    if (nrow(comparative_qc) == 0) {
        stop(
            paste(
                "No QC-passing comparative calls for",
                region_name
            )
        )
    }

    pop_n <- aggregate(
        sample_id ~ population_use,
        data = unique(
            comparative_calls[
                ,
                c("population_use", "sample_id"),
                drop = FALSE
            ]
        ),
        FUN = function(x) length(unique(x))
    )

    names(pop_n)[names(pop_n) == "sample_id"] <- "population_n"

    pop_summary <- aggregate(
        plot_value ~ population_use + site_id + pos,
        data = comparative_qc,
        FUN = function(x) mean(x, na.rm = TRUE)
    )

    pop_summary <- merge(
        pop_summary,
        pop_n,
        by = "population_use",
        all.x = TRUE,
        sort = FALSE
    )

    pop_order <- pop_n$population_use[
        order(
            pop_n$population_n,
            pop_n$population_use,
            decreasing = TRUE
        )
    ]

    pop_labels <- paste0(
        pop_order,
        " (n=",
        pop_n$population_n[match(pop_order, pop_n$population_use)],
        ")"
    )

    pop_summary$display_row <- paste0(
        pop_summary$population_use,
        " (n=",
        pop_summary$population_n,
        ")"
    )

    pop_summary$row_type <- "population"

    individual_order <- c(
        focal_present,
        sort(bwpk_present)
    )

    individual_calls <- q$calls[
        q$calls$sample_id %in% individual_order,
        ,
        drop = FALSE
    ]

    individual_calls$display_row <- individual_calls$sample_id
    individual_calls$row_type <- ifelse(
        individual_calls$sample_id %in% focal_present,
        "focal",
        "bwpk"
    )

    site_order <- unique(
        q$sites[
            ,
            c("site_id", "pos"),
            drop = FALSE
        ]
    )
    site_order <- site_order[
        order(site_order$pos),
        ,
        drop = FALSE
    ]
    site_order$site_index <- seq_len(nrow(site_order))

    pop_summary$site_index <- site_order$site_index[
        match(pop_summary$site_id, site_order$site_id)
    ]
    individual_calls$site_index <- site_order$site_index[
        match(individual_calls$site_id, site_order$site_id)
    ]

    plot_df <- rbind(
        pop_summary[
            ,
            c(
                "site_id",
                "site_index",
                "pos",
                "display_row",
                "row_type",
                "plot_value"
            )
        ],
        individual_calls[
            ,
            c(
                "site_id",
                "site_index",
                "pos",
                "display_row",
                "row_type",
                "plot_value"
            )
        ]
    )

    plot_df <- add_section_divider_rows(
        plot_df = plot_df,
        population_rows = pop_labels,
        individual_rows = individual_order
    )

    desired_rows <- c(
        pop_labels,
        if (length(individual_order) > 0) " " else NULL,
        individual_order
    )

    plot_df$display_row <- factor(
        plot_df$display_row,
        levels = rev(desired_rows)
    )

    plot_df <- plot_df[
        order(
            plot_df$site_index,
            plot_df$display_row
        ),
        ,
        drop = FALSE
    ]

    p <- make_region_plot(
        plot_df = plot_df,
        region_name = region_name,
        chrom = chrom,
        start = start,
        end = end
    )

    pdf_file <- file.path(
        output_dir,
        paste0(
            "Optix_",
            region_name,
            "_population_context_genoplot.pdf"
        )
    )

    png_file <- file.path(
        output_dir,
        paste0(
            "Optix_",
            region_name,
            "_population_context_genoplot.png"
        )
    )

    site_audit_file <- file.path(
        output_dir,
        paste0(
            "Optix_",
            region_name,
            "_population_context.site_order.tsv"
        )
    )

    row_audit_file <- file.path(
        output_dir,
        paste0(
            "Optix_",
            region_name,
            "_population_context.row_order.tsv"
        )
    )

    row_height <- max(
        4.2,
        2.0 + 0.28 * length(desired_rows)
    )

    ggsave(
        filename = pdf_file,
        plot = p,
        width = 10.5,
        height = row_height,
        units = "in",
        useDingbats = FALSE
    )

    ggsave(
        filename = png_file,
        plot = p,
        width = 10.5,
        height = row_height,
        units = "in",
        dpi = 450
    )

    write.table(
        site_order,
        file = site_audit_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    write.table(
        data.frame(
            display_row = desired_rows,
            stringsAsFactors = FALSE
        ),
        file = row_audit_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    plot_list[[region_name]] <- list(
        plot = p,
        height = row_height
    )

    summary_rows[[region_name]] <- data.frame(
        region = region_name,
        chrom = chrom,
        start = start,
        end = end,
        n_sites = nrow(site_order),
        n_populations = length(pop_order),
        n_focal_rows = length(focal_present),
        n_bwpk_rows = length(bwpk_present),
        pdf = pdf_file,
        png = png_file,
        stringsAsFactors = FALSE
    )

    cat(
        "Sites plotted: ", nrow(site_order), "\n",
        "Population rows: ", length(pop_order), "\n",
        "Focal rows: ", length(focal_present), "\n",
        "BWPK rows: ", length(bwpk_present), "\n",
        "PDF: ", pdf_file, "\n\n",
        sep = ""
    )
}


# ============================================================
# WRITE COMBINED MULTIPAGE PDF + SUMMARY
# ============================================================

if (length(plot_list) == 0) {
    stop("No region plots were produced.")
}

combined_pdf <- file.path(
    output_dir,
    "Optix_key_regions_population_context_genoplots.combined.pdf"
)

grDevices::cairo_pdf(
    filename = combined_pdf,
    width = 10.5,
    height = 7.5,
    onefile = TRUE
)

for (region_name in names(plot_list)) {
    print(plot_list[[region_name]]$plot)
}

dev.off()

summary_table <- do.call(
    rbind,
    summary_rows
)

summary_file <- file.path(
    output_dir,
    "Optix_key_regions_population_context_genoplots.summary.tsv"
)

write.table(
    summary_table,
    file = summary_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat(
    "========================================\n",
    "DONE\n",
    "========================================\n\n",
    sep = ""
)

print(summary_table)

cat(
    "\nCombined PDF:\n",
    combined_pdf,
    "\n\nSummary table:\n",
    summary_file,
    "\n",
    sep = ""
)
