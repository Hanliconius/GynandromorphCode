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
    "supplement_full_focal_evidence"
)

dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

input_file <- file.path(
    checkpoint_dir,
    "focal_calls_at_frozen_comparative_sites_biological.wide.tsv"
)

dat <- read.delim(
    input_file,
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
connector_col <- "#8A8A8A"
grid_col <- "#D9D9D9"

draw_module <- function(
    module,
    value_type = c(
        "genotype",
        "read_fraction"
    )
) {
    value_type <- match.arg(
        value_type
    )

    d <- dat[
        dat$module == module,
        ,
        drop = FALSE
    ]

    d <- d[
        order(d$pos),
        ,
        drop = FALSE
    ]

    if (value_type == "genotype") {
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

        ylab <- (
            "Biological presence dosage"
        )
    } else {
        left <- suppressWarnings(
            as.numeric(
                d$left_presence_allele_fraction
            )
        )

        right <- suppressWarnings(
            as.numeric(
                d$right_presence_allele_fraction
            )
        )

        ylab <- (
            "Presence-associated read fraction"
        )
    }

    x <- d$pos / 1e6

    plot(
        x,
        rep(NA_real_, length(x)),
        type = "n",
        ylim = c(-0.05, 1.05),
        xlab = paste0(
            unique(d$chrom)[1],
            " position (Mb)"
        ),
        ylab = ylab,
        main = paste0(
            display_name[[module]],
            "   n=",
            nrow(d)
        ),
        bty = "l"
    )

    for (z in c(0, 0.5, 1)) {
        abline(
            h = z,
            col = grid_col,
            lwd = 0.6
        )
    }

    left_expected <- suppressWarnings(
        as.numeric(
            unique(
                d$left_expected_presence_dosage
            )[1]
        )
    )

    right_expected <- suppressWarnings(
        as.numeric(
            unique(
                d$right_expected_presence_dosage
            )[1]
        )
    )

    abline(
        h = left_expected,
        col = left_col,
        lty = 3,
        lwd = 1.1
    )

    abline(
        h = right_expected,
        col = right_col,
        lty = 3,
        lwd = 1.1
    )

    both <- (
        is.finite(left)
        & is.finite(right)
    )

    if (any(both)) {
        segments(
            x[both],
            left[both],
            x[both],
            right[both],
            col = connector_col,
            lwd = 0.8
        )
    }

    left_ok <- is.finite(left)

    if (any(left_ok)) {
        points(
            x[left_ok],
            left[left_ok],
            pch = 16,
            col = left_col,
            cex = 0.72
        )
    }

    right_ok <- is.finite(right)

    if (any(right_ok)) {
        points(
            x[right_ok],
            right[right_ok],
            pch = 1,
            col = right_col,
            cex = 0.92,
            lwd = 1.2
        )
    }

    legend(
        "bottomleft",
        legend = c(
            "Left",
            "Right"
        ),
        pch = c(16, 1),
        col = c(
            left_col,
            right_col
        ),
        pt.cex = c(
            0.8,
            1.0
        ),
        bty = "n",
        cex = 0.72
    )
}

draw_all <- function(
    value_type
) {
    old_par <- par(
        mfrow = c(3, 2),
        mar = c(4.0, 4.2, 2.6, 1.0),
        oma = c(0.5, 0.5, 2.7, 0.5),
        family = "sans"
    )

    on.exit(
        par(old_par),
        add = TRUE
    )

    for (module in module_order) {
        draw_module(
            module,
            value_type
        )
    }

    top_title <- if (
        value_type == "genotype"
    ) {
        "Full focal evidence: biological genotype dosage"
    } else {
        "Full focal evidence: presence-associated read balance"
    }

    mtext(
        top_title,
        outer = TRUE,
        side = 3,
        line = 1.0,
        cex = 1.18,
        font = 2
    )

    mtext(
        "Sites were frozen from comparative data before focal calls were inspected; QC failures remain as gaps.",
        outer = TRUE,
        side = 3,
        line = -0.1,
        cex = 0.70
    )
}

genotype_pdf <- file.path(
    figure_dir,
    "full_focal_biological_genotype_evidence.pdf"
)

genotype_svg <- file.path(
    figure_dir,
    "full_focal_biological_genotype_evidence.svg"
)

fraction_pdf <- file.path(
    figure_dir,
    "full_focal_presence_read_fraction_evidence.pdf"
)

fraction_svg <- file.path(
    figure_dir,
    "full_focal_presence_read_fraction_evidence.svg"
)

pdf(
    genotype_pdf,
    width = 12,
    height = 13,
    useDingbats = FALSE
)

draw_all(
    "genotype"
)

dev.off()

svg(
    genotype_svg,
    width = 12,
    height = 13,
    pointsize = 12
)

draw_all(
    "genotype"
)

dev.off()

pdf(
    fraction_pdf,
    width = 12,
    height = 13,
    useDingbats = FALSE
)

draw_all(
    "read_fraction"
)

dev.off()

svg(
    fraction_svg,
    width = 12,
    height = 13,
    pointsize = 12
)

draw_all(
    "read_fraction"
)

dev.off()

cat(
    "\nSupplementary focal-evidence figures written:\n",
    genotype_pdf,
    "\n",
    genotype_svg,
    "\n",
    fraction_pdf,
    "\n",
    fraction_svg,
    "\n",
    sep = ""
)
