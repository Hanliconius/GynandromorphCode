#!/usr/bin/env Rscript

options(
    stringsAsFactors = FALSE,
    scipen = 999
)

project_dir <- "/scratch/martinlab/joe/GynandroCalling"

checkpoint_dir <- file.path(
    project_dir,
    "analysis",
    "original_502",
    "checkpoints"
)

figure_dir <- file.path(
    project_dir,
    "analysis",
    "figures",
    "core_explanatory_optix_genoplots"
)

dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

genotype_file <- file.path(
    checkpoint_dir,
    "core_explanatory_optix_genotypes.long.tsv"
)

sample_file <- file.path(
    checkpoint_dir,
    "core_explanatory_optix_samples.tsv"
)

callability_output <- file.path(
    checkpoint_dir,
    "core_explanatory_optix_genoplot_callability.tsv"
)

raw_pdf <- file.path(
    figure_dir,
    "core_explanatory_optix_genoplots_raw.pdf"
)

qc_pdf <- file.path(
    figure_dir,
    "core_explanatory_optix_genoplots_DP2_GQ20.pdf"
)


as_boolean <- function(x) {

    tolower(
        trimws(
            as.character(x)
        )
    ) %in% c(
        "true",
        "t",
        "1",
        "yes"
    )
}


genotypes <- read.delim(
    genotype_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

samples <- read.delim(
    sample_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

genotypes$pos <- as.numeric(
    genotypes$pos
)

genotypes$DP2_GQ20_pass_bool <- as_boolean(
    genotypes$DP2_GQ20_pass
)

samples <- samples[
    order(
        samples$plot_order
    ),
    ,
    drop = FALSE
]

module_order <- c(
    "Band",
    "Ray",
    "Dennis"
)

genotype_colours <- c(
    "A/A" = "#2B83BA",
    "P/A" = "#E68613",
    "P/P" = "#B21875",
    "missing" = "#E5E5E5"
)

group_colours <- c(
    "Focal gynandromorph" = "#222222",
    "H. besckei" = "#B21875",
    "H. t. contigua" = "#E68613",
    "H. m. meriana" = "#2B83BA"
)

original_intervals <- data.frame(
    label = c(
        "B",
        "R1",
        "R2",
        "D1",
        "D2"
    ),
    module = c(
        "Band",
        "Ray",
        "Ray",
        "Dennis",
        "Dennis"
    ),
    start = c(
        1112651,
        1156746,
        1174547,
        1215210,
        1218557
    ),
    end = c(
        1114314,
        1160220,
        1187833,
        1216124,
        1221494
    ),
    stringsAsFactors = FALSE
)


expected_focal_genotype <- function(module_name, sample) {

    expected <- list(
        Band = c(
            TN_Left = "A/A",
            TN_Right = "P/A"
        ),
        Ray = c(
            TN_Left = "P/P",
            TN_Right = "P/A"
        ),
        Dennis = c(
            TN_Left = "P/P",
            TN_Right = "P/A"
        )
    )

    if (
        !module_name %in% names(expected)
        ||
        !sample %in% names(expected[[module_name]])
    ) {
        return("")
    }

    expected[[module_name]][[sample]]
}


interval_for_position <- function(position) {

    matches <- which(
        original_intervals$start
        <=
        position
        &
        original_intervals$end
        >=
        position
    )

    if (length(matches) == 0) {
        return("outside")
    }

    original_intervals$label[
        matches[[1]]
    ]
}


make_matrix <- function(module_name, qc_mode) {

    d <- genotypes[
        genotypes$module
        ==
        module_name,
        ,
        drop = FALSE
    ]

    sites <- unique(
        d[
            ,
            c(
                "chrom",
                "pos",
                "ref",
                "alt",
                "minimum_BA_all_analyses"
            ),
            drop = FALSE
        ]
    )

    sites <- sites[
        order(
            sites$pos,
            sites$ref,
            sites$alt
        ),
        ,
        drop = FALSE
    ]

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

    matrix_values <- matrix(
        NA_character_,
        nrow = nrow(samples),
        ncol = nrow(sites),
        dimnames = list(
            samples$sample,
            sites$site_key
        )
    )

    row_index <- match(
        d$sample,
        samples$sample
    )

    column_index <- match(
        d$site_key,
        sites$site_key
    )

    values <- d$biological_GT

    values[
        !values
        %in%
        c(
            "A/A",
            "P/A",
            "P/P"
        )
    ] <- NA_character_

    if (qc_mode) {
        values[
            !d$DP2_GQ20_pass_bool
        ] <- NA_character_
    }

    valid <- (
        !is.na(row_index)
        &
        !is.na(column_index)
    )

    matrix_values[
        cbind(
            row_index[valid],
            column_index[valid]
        )
    ] <- values[valid]

    list(
        matrix = matrix_values,
        sites = sites
    )
}


make_sample_labels <- function(module_name) {

    vapply(
        seq_len(
            nrow(samples)
        ),
        function(index) {

            sample <- samples$sample[
                index
            ]

            label <- samples$sample_label[
                index
            ]

            expected <- expected_focal_genotype(
                module_name,
                sample
            )

            if (nzchar(expected)) {
                label <- paste0(
                    label,
                    " | expected ",
                    expected
                )
            }

            label
        },
        character(1)
    )
}


draw_genoplot <- function(module_name, qc_mode) {

    prepared <- make_matrix(
        module_name,
        qc_mode
    )

    genotype_matrix <- prepared$matrix
    sites <- prepared$sites

    sample_n <- nrow(
        genotype_matrix
    )

    site_n <- ncol(
        genotype_matrix
    )

    row_y <- sample_n:1
    column_x <- seq_len(
        site_n
    )

    sample_labels <- make_sample_labels(
        module_name
    )

    par(
        mar = c(
            5.1,
            32,
            4.4,
            1.4
        ),
        xpd = NA,
        las = 1
    )

    plot(
        NA,
        xlim = c(
            0.5,
            site_n + 0.5
        ),
        ylim = c(
            0.5,
            sample_n + 3.1
        ),
        axes = FALSE,
        xlab = "",
        ylab = "",
        bty = "n"
    )

    # Missing-data background.
    rect(
        0.5,
        0.5,
        site_n + 0.5,
        sample_n + 0.5,
        col = genotype_colours[
            "missing"
        ],
        border = NA
    )

    for (
        genotype_class
        in
        c(
            "A/A",
            "P/A",
            "P/P"
        )
    ) {

        matrix_positions <- which(
            genotype_matrix
            ==
            genotype_class,
            arr.ind = TRUE
        )

        if (nrow(matrix_positions) == 0) {
            next
        }

        x <- matrix_positions[
            ,
            2
        ]

        y <- sample_n - matrix_positions[
            ,
            1
        ] + 1

        rect(
            x - 0.5,
            y - 0.5,
            x + 0.5,
            y + 0.5,
            col = genotype_colours[
                genotype_class
            ],
            border = NA
        )
    }

    # Fine sample-row guides.
    abline(
        h = seq(
            0.5,
            sample_n + 0.5,
            by = 1
        ),
        col = adjustcolor(
            "#FFFFFF",
            alpha.f = 0.22
        ),
        lwd = 0.25
    )

    # Strong separators between biological groups.
    group_ends <- cumsum(
        as.numeric(
            table(
                factor(
                    samples$group,
                    levels = unique(
                        samples$group
                    )
                )
            )
        )
    )

    group_ends <- group_ends[
        group_ends < sample_n
    ]

    for (group_end in group_ends) {

        separator_y <- (
            sample_n
            -
            group_end
            +
            0.5
        )

        abline(
            h = separator_y,
            col = "#111111",
            lwd = 1.15
        )
    }

    # Original functional intervals as a strip above the matrix.
    interval_classes <- vapply(
        sites$pos,
        interval_for_position,
        character(1)
    )

    interval_colour_lookup <- c(
        "B" = "#B21875",
        "R1" = "#E68613",
        "R2" = "#E68613",
        "D1" = "#2B83BA",
        "D2" = "#2B83BA",
        "outside" = "#F1F1F1"
    )

    annotation_y_bottom <- sample_n + 1.05
    annotation_y_top <- sample_n + 1.45

    for (column_index in column_x) {

        interval_class <- interval_classes[
            column_index
        ]

        rect(
            column_index - 0.5,
            annotation_y_bottom,
            column_index + 0.5,
            annotation_y_top,
            col = interval_colour_lookup[
                interval_class
            ],
            border = NA
        )
    }

    for (
        interval_label
        in
        c(
            "B",
            "R1",
            "R2",
            "D1",
            "D2"
        )
    ) {

        columns <- which(
            interval_classes
            ==
            interval_label
        )

        if (length(columns) == 0) {
            next
        }

        text(
            mean(
                range(
                    columns
                )
            ),
            sample_n + 1.75,
            labels = interval_label,
            cex = 0.72,
            font = 2,
            col = interval_colour_lookup[
                interval_label
            ]
        )
    }

    axis(
        side = 2,
        at = row_y,
        labels = sample_labels,
        tick = FALSE,
        cex.axis = 0.58,
        line = -0.7
    )

    # Colour focal and diagnostic-taxon labels with group marks.
    for (index in seq_len(nrow(samples))) {

        text(
            x = -0.006 * site_n,
            y = row_y[
                index
            ],
            labels = "\u25A0",
            col = group_colours[
                samples$group[
                    index
                ]
            ],
            cex = 0.75,
            adj = 1
        )
    }

    tick_columns <- unique(
        round(
            seq(
                1,
                site_n,
                length.out = min(
                    8,
                    site_n
                )
            )
        )
    )

    axis(
        side = 1,
        at = tick_columns,
        labels = sprintf(
            "%.3f",
            sites$pos[
                tick_columns
            ] / 1000000
        ),
        cex.axis = 0.75,
        tck = -0.012,
        mgp = c(
            2.4,
            0.6,
            0
        )
    )

    mtext(
        paste0(
            "Hmr18 position (Mb); ",
            "SNP columns are equally spaced in physical order"
        ),
        side = 1,
        line = 2.8,
        cex = 0.82
    )

    mode_label <- if (qc_mode) {
        "DP >= 2 and GQ >= 20"
    } else {
        "all called genotypes"
    }

    title(
        main = paste0(
            module_name,
            ": core explanatory genotypes"
        ),
        sub = paste0(
            site_n,
            " independently selected SNPs; ",
            mode_label
        ),
        cex.main = 1.05,
        cex.sub = 0.78,
        line = 1.1
    )

    legend(
        x = 0.5,
        y = 0.02,
        legend = c(
            "A/A",
            "P/A",
            "P/P",
            "missing/failed QC"
        ),
        fill = genotype_colours[
            c(
                "A/A",
                "P/A",
                "P/P",
                "missing"
            )
        ],
        border = NA,
        horiz = TRUE,
        bty = "n",
        cex = 0.78,
        xpd = NA
    )
}


write_genoplot_pdf <- function(filename, qc_mode) {

    pdf(
        filename,
        width = 16,
        height = 10.5,
        useDingbats = FALSE,
        onefile = TRUE
    )

    for (module_name in module_order) {
        draw_genoplot(
            module_name,
            qc_mode
        )
    }

    dev.off()
}


write_genoplot_pdf(
    raw_pdf,
    qc_mode = FALSE
)

write_genoplot_pdf(
    qc_pdf,
    qc_mode = TRUE
)


callability_rows <- list()

row_counter <- 1

for (module_name in module_order) {

    d <- genotypes[
        genotypes$module
        ==
        module_name,
        ,
        drop = FALSE
    ]

    for (sample_name in samples$sample) {

        x <- d[
            d$sample
            ==
            sample_name,
            ,
            drop = FALSE
        ]

        callability_rows[[
            row_counter
        ]] <- data.frame(
            module = module_name,
            sample = sample_name,
            selected_sites = nrow(x),
            raw_called = sum(
                x$biological_GT
                %in%
                c(
                    "A/A",
                    "P/A",
                    "P/P"
                )
            ),
            DP2_GQ20_called = sum(
                x$DP2_GQ20_pass_bool
            ),
            DP2_GQ20_call_rate = mean(
                x$DP2_GQ20_pass_bool
            ),
            stringsAsFactors = FALSE
        )

        row_counter <- row_counter + 1
    }
}

callability <- do.call(
    rbind,
    callability_rows
)

write.table(
    callability,
    callability_output,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

cat(
    "CORE EXPLANATORY OPTIX GENOPLOTS\n\n"
)

cat(
    "Samples: ",
    nrow(samples),
    "\n",
    sep = ""
)

cat(
    "Band SNPs: ",
    sum(
        genotypes$module
        ==
        "Band"
    ) / nrow(samples),
    "\n",
    sep = ""
)

cat(
    "Ray SNPs: ",
    sum(
        genotypes$module
        ==
        "Ray"
    ) / nrow(samples),
    "\n",
    sep = ""
)

cat(
    "Dennis SNPs: ",
    sum(
        genotypes$module
        ==
        "Dennis"
    ) / nrow(samples),
    "\n\n",
    sep = ""
)

cat(
    "Raw PDF:\n",
    raw_pdf,
    "\n\nQC PDF:\n",
    qc_pdf,
    "\n\nCallability:\n",
    callability_output,
    "\n",
    sep = ""
)
