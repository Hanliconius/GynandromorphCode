options(stringsAsFactors = FALSE)

project_root <- "/scratch/martinlab/joe/GynandroCalling"

checkpoint_dir <- file.path(
    project_root,
    "analysis",
    "original_502",
    "checkpoints"
)

figure_dir <- file.path(
    project_root,
    "analysis",
    "figures",
    "main_biological_focal_evidence"
)

dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

wide_file <- file.path(
    checkpoint_dir,
    "focal_calls_at_frozen_comparative_sites_biological.wide.tsv"
)

span_file <- file.path(
    checkpoint_dir,
    "focal_module_locus_spans.tsv"
)

pdf_file <- file.path(
    figure_dir,
    "focal_locus_cartoon_biological_genotypes.pdf"
)

svg_file <- file.path(
    figure_dir,
    "focal_locus_cartoon_biological_genotypes.svg"
)

dat <- read.delim(
    wide_file,
    check.names = FALSE
)

spans <- read.delim(
    span_file,
    check.names = FALSE
)

module_order <- c(
    "Band",
    "Ray1",
    "Ray2",
    "Dennis1",
    "Dennis2",
    "Yb"
)

display_name <- c(
    Band = "Band",
    Ray1 = "Ray 1",
    Ray2 = "Ray 2",
    Dennis1 = "Dennis 1",
    Dennis2 = "Dennis 2",
    Yb = "Yb"
)

left_col <- "#0072B2"
right_col <- "#D55E00"

dark <- "#222222"
mid <- "#777777"
light <- "#D9D9D9"
lighter <- "#F4F4F4"

module_fill <- c(
    Band = "#E8E8E8",
    Ray1 = "#D9D9D9",
    Ray2 = "#C9C9C9",
    Dennis1 = "#B9B9B9",
    Dennis2 = "#A9A9A9",
    Yb = "#D9D9D9"
)

genotype_label <- function(x) {
    x <- as.character(x)

    switch(
        x,
        AA = "A/A",
        PA = "P/A",
        PP = "P/P",
        x
    )
}

scale_value <- function(
    x,
    in_min,
    in_max,
    out_min,
    out_max
) {
    if (!is.finite(in_min) ||
        !is.finite(in_max) ||
        in_max <= in_min) {
        return(
            rep(
                mean(c(out_min, out_max)),
                length(x)
            )
        )
    }

    out_min +
        (
            (x - in_min) /
                (in_max - in_min)
        ) *
        (out_max - out_min)
}

draw_overview <- function(
    modules,
    chromosome,
    y,
    x0,
    x1,
    label
) {
    s <- spans[
        spans$module %in% modules,
        ,
        drop = FALSE
    ]

    start <- min(s$start)
    end <- max(s$end)

    pad <- max(
        1,
        round(
            0.015 *
                (end - start)
        )
    )

    start <- start - pad
    end <- end + pad

    segments(
        x0,
        y,
        x1,
        y,
        col = dark,
        lwd = 1.2
    )

    text(
        x0,
        y + 0.035,
        label,
        adj = c(0, 0.5),
        cex = 0.90,
        font = 2
    )

    for (module in modules) {
        sm <- s[
            s$module == module,
            ,
            drop = FALSE
        ]

        if (nrow(sm) == 0) {
            next
        }

        mx0 <- scale_value(
            sm$start,
            start,
            end,
            x0,
            x1
        )

        mx1 <- scale_value(
            sm$end,
            start,
            end,
            x0,
            x1
        )

        rect(
            mx0,
            y - 0.012,
            mx1,
            y + 0.012,
            col = module_fill[[module]],
            border = dark,
            lwd = 0.8
        )

        text(
            mean(c(mx0, mx1)),
            y + 0.027,
            display_name[[module]],
            cex = 0.72,
            font = 2
        )

        d <- dat[
            dat$module == module,
            ,
            drop = FALSE
        ]

        if (nrow(d) > 0) {
            sx <- scale_value(
                d$pos,
                start,
                end,
                x0,
                x1
            )

            segments(
                sx,
                y - 0.021,
                sx,
                y - 0.013,
                col = dark,
                lwd = 0.55
            )
        }
    }

    text(
        x0,
        y - 0.035,
        paste0(
            chromosome,
            ": ",
            format(
                start,
                big.mark = ",",
                scientific = FALSE
            ),
            "–",
            format(
                end,
                big.mark = ",",
                scientific = FALSE
            )
        ),
        adj = c(0, 0.5),
        cex = 0.66,
        col = mid
    )

    invisible(
        list(
            start = start,
            end = end
        )
    )
}

allocate_zoom_boxes <- function(
    modules,
    x0 = 0.065,
    x1 = 0.965,
    gap = 0.012,
    minimum_span = 2000
) {
    s <- spans[
        match(
            modules,
            spans$module
        ),
        ,
        drop = FALSE
    ]

    weight <- pmax(
        s$span_bp,
        minimum_span
    )

    usable <- (
        x1 - x0
        - gap *
            (length(modules) - 1)
    )

    widths <- (
        usable *
            weight /
            sum(weight)
    )

    starts <- numeric(
        length(modules)
    )

    ends <- numeric(
        length(modules)
    )

    cursor <- x0

    for (i in seq_along(modules)) {
        starts[i] <- cursor
        ends[i] <- (
            cursor + widths[i]
        )

        cursor <- (
            ends[i] + gap
        )
    }

    data.frame(
        module = modules,
        x0 = starts,
        x1 = ends,
        stringsAsFactors = FALSE
    )
}

draw_zoom <- function(
    module,
    x0,
    x1,
    y0,
    y1
) {
    d <- dat[
        dat$module == module,
        ,
        drop = FALSE
    ]

    s <- spans[
        spans$module == module,
        ,
        drop = FALSE
    ]

    if (nrow(d) == 0 ||
        nrow(s) == 0) {
        return(
            invisible(NULL)
        )
    }

    rect(
        x0,
        y0,
        x1,
        y1,
        border = light,
        col = lighter,
        lwd = 0.8
    )

    plot_y0 <- y0 + 0.055
    plot_y1 <- y1 - 0.070

    y_of <- function(z) {
        scale_value(
            z,
            0,
            1,
            plot_y0,
            plot_y1
        )
    }

    for (z in c(0, 0.5, 1)) {
        segments(
            x0,
            y_of(z),
            x1,
            y_of(z),
            col = light,
            lwd = 0.55
        )
    }

    start <- s$start[1]
    end <- s$end[1]

    pad <- max(
        1,
        round(
            0.025 *
                max(
                    1,
                    end - start
                )
        )
    )

    sx <- scale_value(
        d$pos,
        start - pad,
        end + pad,
        x0 + 0.008,
        x1 - 0.008
    )

    left_expected <- unique(
        d$left_expected_presence_dosage
    )[1]

    right_expected <- unique(
        d$right_expected_presence_dosage
    )[1]

    segments(
        x0,
        y_of(left_expected),
        x1,
        y_of(left_expected),
        col = left_col,
        lty = 3,
        lwd = 1.0
    )

    segments(
        x0,
        y_of(right_expected),
        x1,
        y_of(right_expected),
        col = right_col,
        lty = 3,
        lwd = 1.0
    )

    left <- suppressWarnings(
        as.numeric(
            d$left_presence_dosage
        )
    )

    right <- suppressWarnings(
        as.numeric(
            d$right_presence_dosage
        )
    )

    both <- (
        is.finite(left)
        & is.finite(right)
    )

    if (any(both)) {
        segments(
            sx[both],
            y_of(left[both]),
            sx[both],
            y_of(right[both]),
            col = "#8A8A8A",
            lwd = 0.85
        )
    }

    left_ok <- is.finite(left)

    if (any(left_ok)) {
        points(
            sx[left_ok],
            y_of(left[left_ok]),
            pch = 16,
            cex = 0.72,
            col = left_col
        )
    }

    right_ok <- is.finite(right)

    if (any(right_ok)) {
        points(
            sx[right_ok],
            y_of(right[right_ok]),
            pch = 1,
            cex = 0.92,
            lwd = 1.2,
            col = right_col
        )
    }

    segments(
        sx,
        y0 + 0.006,
        sx,
        y0 + 0.018,
        col = dark,
        lwd = 0.45
    )

    left_gt <- genotype_label(
        unique(
            d$left_expected_genotype
        )[1]
    )

    right_gt <- genotype_label(
        unique(
            d$right_expected_genotype
        )[1]
    )

    text(
        mean(c(x0, x1)),
        y1 - 0.020,
        display_name[[module]],
        cex = 0.84,
        font = 2
    )

    text(
        mean(c(x0, x1)),
        y1 - 0.045,
        paste0(
            "expected  L ",
            left_gt,
            "   |   R ",
            right_gt
        ),
        cex = 0.62
    )

    text(
        x0 + 0.002,
        y_of(1),
        "1",
        adj = c(0, 0.5),
        cex = 0.55,
        col = mid
    )

    text(
        x0 + 0.002,
        y_of(0.5),
        "0.5",
        adj = c(0, 0.5),
        cex = 0.55,
        col = mid
    )

    text(
        x0 + 0.002,
        y_of(0),
        "0",
        adj = c(0, 0.5),
        cex = 0.55,
        col = mid
    )

    n_both <- sum(
        d$both_callable %in%
            c(TRUE, "TRUE", "True")
    )

    n_match <- sum(
        d$pair_matches_expectation %in%
            c(TRUE, "TRUE", "True"),
        na.rm = TRUE
    )

    text(
        mean(c(x0, x1)),
        y0 + 0.030,
        paste0(
            "n=",
            nrow(d),
            "  |  callable pairs=",
            n_both,
            "  |  exact=",
            n_match
        ),
        cex = 0.55,
        col = mid
    )
}

draw_page <- function() {
    par(
        mar = c(0, 0, 0, 0),
        xpd = NA,
        family = "sans"
    )

    plot.new()

    plot.window(
        xlim = c(0, 1),
        ylim = c(0, 1)
    )

    text(
        0.055,
        0.965,
        "Focal left–right genotypes at independently ascertained comparative sites",
        adj = c(0, 0.5),
        cex = 1.28,
        font = 2
    )

    text(
        0.055,
        0.932,
        "Genotypes are polarized to biological phenotype-associated presence: A/A = 0, P/A = 0.5, P/P = 1",
        adj = c(0, 0.5),
        cex = 0.76,
        col = mid
    )

    points(
        0.66,
        0.932,
        pch = 16,
        col = left_col,
        cex = 0.8
    )

    text(
        0.675,
        0.932,
        "Left",
        adj = c(0, 0.5),
        cex = 0.70
    )

    points(
        0.735,
        0.932,
        pch = 1,
        col = right_col,
        cex = 0.95,
        lwd = 1.2
    )

    text(
        0.750,
        0.932,
        "Right",
        adj = c(0, 0.5),
        cex = 0.70
    )

    segments(
        0.815,
        0.932,
        0.850,
        0.932,
        col = dark,
        lwd = 0.8
    )

    text(
        0.862,
        0.932,
        "same SNP",
        adj = c(0, 0.5),
        cex = 0.70
    )

    optix_modules <- c(
        "Band",
        "Ray1",
        "Ray2",
        "Dennis1",
        "Dennis2"
    )

    draw_overview(
        modules = optix_modules,
        chromosome = "Hmr18",
        y = 0.845,
        x0 = 0.065,
        x1 = 0.965,
        label = "Optix locus — genomic overview"
    )

    boxes <- allocate_zoom_boxes(
        optix_modules,
        x0 = 0.065,
        x1 = 0.965,
        gap = 0.012,
        minimum_span = 2000
    )

    for (i in seq_len(nrow(boxes))) {
        module <- boxes$module[i]

        sm <- spans[
            spans$module == module,
            ,
            drop = FALSE
        ]

        overview_all <- spans[
            spans$module %in%
                optix_modules,
            ,
            drop = FALSE
        ]

        overview_start <- min(
            overview_all$start
        )

        overview_end <- max(
            overview_all$end
        )

        locus_mid <- mean(
            c(
                sm$start[1],
                sm$end[1]
            )
        )

        leader_x <- scale_value(
            locus_mid,
            overview_start,
            overview_end,
            0.065,
            0.965
        )

        box_mid <- mean(
            c(
                boxes$x0[i],
                boxes$x1[i]
            )
        )

        # Orthogonal callout, deliberately not the
        # left/right genotype connector.
        segments(
            leader_x,
            0.823,
            leader_x,
            0.790,
            col = light,
            lwd = 0.7
        )

        segments(
            leader_x,
            0.790,
            box_mid,
            0.790,
            col = light,
            lwd = 0.7
        )

        segments(
            box_mid,
            0.790,
            box_mid,
            0.755,
            col = light,
            lwd = 0.7
        )

        draw_zoom(
            module,
            boxes$x0[i],
            boxes$x1[i],
            y0 = 0.455,
            y1 = 0.750
        )
    }

    text(
        0.065,
        0.425,
        "Expanded Optix intervals: widths scale with genomic span (2-kb minimum); inter-module gaps are compressed.",
        adj = c(0, 0.5),
        cex = 0.63,
        col = mid
    )

    draw_overview(
        modules = "Yb",
        chromosome = "Hmr15",
        y = 0.355,
        x0 = 0.065,
        x1 = 0.965,
        label = "Ivory / Yb locus — genomic overview"
    )

    draw_zoom(
        "Yb",
        x0 = 0.245,
        x1 = 0.785,
        y0 = 0.075,
        y1 = 0.295
    )

    text(
        0.055,
        0.025,
        "Dashed lines show the pre-specified left and right genotype expectations. QC-failing focal calls are gaps and never remove a comparative site.",
        adj = c(0, 0.5),
        cex = 0.65,
        col = mid
    )
}

pdf(
    pdf_file,
    width = 18,
    height = 10,
    useDingbats = FALSE
)

draw_page()

dev.off()

svg(
    svg_file,
    width = 18,
    height = 10,
    pointsize = 12
)

draw_page()

dev.off()

cat(
    "\nMain biological-space figure written:\n",
    pdf_file,
    "\n",
    svg_file,
    "\n",
    sep = ""
)
