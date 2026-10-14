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
    "supplement_comparative_site_ascertainment"
)

dir.create(
    figure_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

input_file <- file.path(
    checkpoint_dir,
    "comparative_site_selection_audit.tsv"
)

pdf_file <- file.path(
    figure_dir,
    "comparative_BA80_site_ascertainment.pdf"
)

svg_file <- file.path(
    figure_dir,
    "comparative_BA80_site_ascertainment.svg"
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

all_col <- "#BDBDBD"
selected_col <- "#0072B2"
threshold_col <- "#444444"

is_true <- function(x) {
    x %in% c(
        TRUE,
        "TRUE",
        "True",
        "true",
        1,
        "1"
    )
}

draw_figure <- function() {
    old_par <- par(
        mfrow = c(3, 2),
        mar = c(4.0, 4.5, 2.8, 1.0),
        oma = c(0.5, 0.5, 3.2, 0.5),
        family = "sans"
    )

    on.exit(
        par(old_par),
        add = TRUE
    )

    for (module in module_order) {
        d <- dat[
            dat$module == module,
            ,
            drop = FALSE
        ]

        d$selection_BA <- suppressWarnings(
            as.numeric(
                d$selection_BA
            )
        )

        d$pos <- suppressWarnings(
            as.numeric(
                d$pos
            )
        )

        d <- d[
            is.finite(d$selection_BA)
            & is.finite(d$pos),
            ,
            drop = FALSE
        ]

        selected <- is_true(
            d$selected
        )

        if (nrow(d) == 0) {
            plot.new()

            title(
                main = display_name[[module]]
            )

            text(
                0.5,
                0.5,
                "No comparative score rows"
            )

            next
        }

        x <- d$pos / 1e6

        ylim <- range(
            c(
                0.45,
                1.0,
                d$selection_BA
            ),
            finite = TRUE
        )

        plot(
            x,
            d$selection_BA,
            type = "n",
            xlab = paste0(
                unique(d$chrom)[1],
                " position (Mb)"
            ),
            ylab = "Comparative balanced accuracy",
            ylim = ylim,
            main = paste0(
                display_name[[module]],
                "   (selected n=",
                sum(selected),
                ")"
            ),
            bty = "l"
        )

        abline(
            h = 0.80,
            lty = 2,
            lwd = 1.2,
            col = threshold_col
        )

        points(
            x[!selected],
            d$selection_BA[!selected],
            pch = 16,
            cex = 0.42,
            col = all_col
        )

        if (any(selected)) {
            points(
                x[selected],
                d$selection_BA[selected],
                pch = 16,
                cex = 0.72,
                col = selected_col
            )
        }

        mtext(
            "BA >= 0.80",
            side = 4,
            at = 0.80,
            las = 1,
            line = 0.1,
            cex = 0.60,
            col = threshold_col
        )
    }

    mtext(
        "Comparative-only ascertainment of focal-site columns",
        outer = TRUE,
        side = 3,
        line = 1.3,
        cex = 1.20,
        font = 2
    )

    mtext(
        "Focal GT, AD, DP, GQ and expected left/right zygosity are not used in site selection",
        outer = TRUE,
        side = 3,
        line = 0.1,
        cex = 0.72
    )
}

pdf(
    pdf_file,
    width = 12,
    height = 13,
    useDingbats = FALSE
)

draw_figure()

dev.off()

svg(
    svg_file,
    width = 12,
    height = 13,
    pointsize = 12
)

draw_figure()

dev.off()

cat(
    "\nSupplementary ascertainment figure written:\n",
    pdf_file,
    "\n",
    svg_file,
    "\n",
    sep = ""
)
