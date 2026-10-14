#!/usr/bin/env Rscript

project <- "/scratch/martinlab/joe/GynandroCalling"

ck <- file.path(
    project,
    "analysis",
    "original_502",
    "checkpoints"
)

long_file <- file.path(
    ck,
    "optix_comparative_selected_sites_genotype_plot.long.tsv"
)

sites_file <- file.path(
    ck,
    "optix_comparative_selected_sites_for_plotting.tsv"
)

out_dir <- file.path(
    project,
    "analysis",
    "figures",
    "supplement_optix_predictive_genotype_plots"
)

dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

long <- read.delim(
    long_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

sites <- read.delim(
    sites_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
)

module_order <- c(
    "Band",
    "Ray1",
    "Ray2",
    "Dennis1",
    "Dennis2"
)

display_names <- c(
    Band = "Band",
    Ray1 = "Ray 1",
    Ray2 = "Ray 2",
    Dennis1 = "Dennis 1",
    Dennis2 = "Dennis 2"
)

# Requested palette: magenta / orange / blue
gt_cols <- c(
    "AA" = "#2F6BFF",   # blue
    "PA" = "#F39C12",   # orange
    "PP" = "#D81BCE",   # magenta
    "missing" = "#E5E7EB"
)

pheno_cols <- c(
    "absent" = "#FFFFFF",
    "present" = "#111111"
)

pdf_file <- file.path(
    out_dir,
    "optix_predictive_score_plus_genotype_matrix.pdf"
)

module_heights <- c()
modules_present <- c()

for (module in module_order) {
    d <- long[long$module == module, , drop = FALSE]
    if (nrow(d) == 0) {
        next
    }

    ns <- length(unique(d$sample))
    # BA track + genotype matrix
    module_heights <- c(module_heights, 1.0, max(3.0, ns * 0.22))
    modules_present <- c(modules_present, module)
}

if (length(modules_present) == 0) {
    stop("No modules with data found.")
}

pdf(
    pdf_file,
    width = 13.5,
    height = max(10, sum(module_heights) * 0.55),
    useDingbats = FALSE
)

layout(
    matrix(seq_len(length(module_heights)), ncol = 1),
    heights = module_heights
)

par(oma = c(1.2, 1.2, 2.2, 0.8))

draw_module <- function(module) {

    d <- long[long$module == module, , drop = FALSE]
    s <- sites[sites$module == module, , drop = FALSE]

    s <- s[order(s$site_order_within_module, s$pos), , drop = FALSE]

    samples <- unique(
        d[order(d$sample_order), "sample"]
    )

    sample_df <- unique(
        d[, c(
            "sample",
            "sample_order",
            "taxon_name",
            "phenotype_label",
            "phenotype_value"
        )]
    )

    sample_df <- sample_df[
        order(sample_df$sample_order),
        ,
        drop = FALSE
    ]

    x_pos <- seq_len(nrow(s))
    y_pos <- seq_len(nrow(sample_df))

    # --------------------------------------------------------
    # Top: BA track
    # --------------------------------------------------------

    par(mar = c(0.6, 8.2, 2.2, 4.8))

    plot(
        NA,
        xlim = c(0.5, nrow(s) + 0.5),
        ylim = c(0.48, 1.02),
        xaxt = "n",
        xlab = "",
        ylab = "BA",
        main = "",
        bty = "l"
    )

    abline(h = 0.80, lty = 2, lwd = 0.9, col = "grey35")
    abline(h = c(0.5, 1.0), lty = 3, lwd = 0.5, col = "grey88")

    segments(
        x0 = x_pos,
        y0 = 0.5,
        x1 = x_pos,
        y1 = s$plot_BA,
        col = "grey55",
        lwd = 0.8
    )

    points(
        x_pos,
        s$plot_BA,
        pch = 16,
        cex = 0.75
    )

    title(
        main = sprintf(
            "%s  (%d selected SNPs; %d comparative samples)",
            display_names[[module]],
            nrow(s),
            nrow(sample_df)
        ),
        line = 0.5,
        font.main = 2,
        cex.main = 0.95
    )

    mtext(
        "balanced accuracy",
        side = 2,
        line = 2.2,
        cex = 0.72
    )

    text(
        nrow(s) + 0.45,
        0.80,
        "0.80",
        pos = 2,
        cex = 0.60,
        col = "grey25"
    )

    # --------------------------------------------------------
    # Bottom: genotype matrix
    # --------------------------------------------------------

    par(mar = c(3.2, 8.2, 0.3, 4.8))

    plot(
        NA,
        xlim = c(-1.9, nrow(s) + 2.8),
        ylim = c(nrow(sample_df) + 0.5, 0.5),
        xaxt = "n",
        yaxt = "n",
        xlab = "comparative-selected SNPs in genomic order",
        ylab = "",
        bty = "n"
    )

    # phenotype strip title
    text(
        -1.25,
        0.25,
        "Phen.",
        xpd = NA,
        cex = 0.58,
        font = 2
    )

    # taxon / sample label title
    text(
        -0.05,
        0.25,
        "Sample / taxon",
        xpd = NA,
        adj = 0,
        cex = 0.58,
        font = 2
    )

    # genotype header key
    legend(
        "topright",
        inset = c(0.0, -0.02),
        horiz = TRUE,
        bty = "n",
        xpd = NA,
        legend = c("AA", "PA", "PP", "missing"),
        fill = gt_cols[c("AA", "PA", "PP", "missing")],
        border = "grey55",
        cex = 0.62
    )

    # row backgrounds / labels
    for (i in seq_len(nrow(sample_df))) {

        row <- sample_df[i, , drop = FALSE]
        yy <- i

        # phenotype strip
        rect(
            -1.6, yy - 0.5,
            -0.9, yy + 0.5,
            col = pheno_cols[[row$phenotype_label]],
            border = "grey55",
            lwd = 0.5
        )

        # row label
        label <- paste0(
            row$sample,
            "   ",
            row$taxon_name
        )

        text(
            -0.75,
            yy,
            labels = label,
            adj = 0,
            cex = 0.43
        )
    }

    # horizontal separators between absent and present groups
    phen <- sample_df$phenotype_label
    phen_breaks <- which(phen[-1] != phen[-length(phen)])
    if (length(phen_breaks) > 0) {
        for (b in phen_breaks) {
            segments(
                -1.8,
                b + 0.5,
                nrow(s) + 0.5,
                b + 0.5,
                lwd = 1.2,
                col = "grey35"
            )
        }
    }

    # vertical light guides
    for (x in x_pos) {
        segments(
            x + 0.5, 0.5,
            x + 0.5, nrow(sample_df) + 0.5,
            lwd = 0.35,
            col = "grey92"
        )
    }

    # draw genotype tiles
    for (j in seq_len(nrow(s))) {

        sj <- s[j, , drop = FALSE]

        dj <- d[
            d$site_key == sj$site_key,
            ,
            drop = FALSE
        ]

        dj <- dj[
            order(dj$sample_order),
            ,
            drop = FALSE
        ]

        for (i in seq_len(nrow(sample_df))) {
            yy <- i
            sample_i <- sample_df$sample[i]

            hit <- dj[dj$sample == sample_i, , drop = FALSE]

            gt <- if (nrow(hit) == 1) {
                as.character(hit$biological_genotype)
            } else {
                "missing"
            }

            if (!(gt %in% names(gt_cols))) {
                gt <- "missing"
            }

            rect(
                j - 0.5, yy - 0.5,
                j + 0.5, yy + 0.5,
                col = gt_cols[[gt]],
                border = "grey70",
                lwd = 0.35
            )
        }
    }

    # x-axis labels = genomic positions
    axis(
        side = 1,
        at = x_pos,
        labels = s$pos,
        las = 2,
        tick = FALSE,
        line = -0.5,
        cex.axis = 0.42
    )

    # phenotype group labels at left
    absent_idx <- which(sample_df$phenotype_label == "absent")
    present_idx <- which(sample_df$phenotype_label == "present")

    if (length(absent_idx) > 0) {
        text(
            -1.95,
            mean(absent_idx),
            "absent",
            srt = 90,
            cex = 0.52
        )
    }

    if (length(present_idx) > 0) {
        text(
            -1.95,
            mean(present_idx),
            "present",
            srt = 90,
            cex = 0.52,
            col = "grey15"
        )
    }
}

for (module in modules_present) {
    draw_module(module)
}

mtext(
    "Comparative genotype plots at frozen predictive sites",
    side = 3,
    outer = TRUE,
    line = 0.6,
    font = 2,
    cex = 1.05
)

mtext(
    "Rows are comparative butterflies only; columns are frozen selected SNPs. Genotypes are polarized into biological presence-space (AA blue, PA orange, PP magenta).",
    side = 3,
    outer = TRUE,
    line = -0.5,
    cex = 0.72
)

dev.off()

cat(
    "\nOutput:\n  ",
    pdf_file,
    "\n",
    sep = ""
)
