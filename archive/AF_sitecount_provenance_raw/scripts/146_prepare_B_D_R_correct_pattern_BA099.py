#!/usr/bin/env python3

from pathlib import Path

project = Path(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
)

source = (
    project
    / "scripts"
    / "145_plot_B_D_R_correct_pattern_ribbons.R"
)

target = (
    project
    / "scripts"
    / "146_plot_B_D_R_correct_pattern_BA099_ribbons.R"
)


if not source.exists():
    raise FileNotFoundError(
        f"Missing source script: {source}"
    )


text = source.read_text()


def replace_section(
    text,
    start_heading,
    end_heading,
    replacement,
):
    start = text.index(start_heading)
    end = text.index(
        end_heading,
        start,
    )

    return (
        text[:start]
        + replacement.rstrip()
        + "\n\n"
        + text[end:]
    )


# ============================================================
# D/R:
# retain all correct-direction P/P -> P/A sites;
# highlight only those with Dennis OR Ray BA > 0.99.
# ============================================================

dr_block = r'''
# D/R RIBBON
#
# First retain the correct focal inheritance direction across
# the whole broad region:
#
#     TN_Left  = 1/1 = biological P/P
#     TN_Right = 0/1 = biological P/A
#
# Then annotate each site with BOTH Dennis and Ray comparative
# BA.  A site is highlighted when either BA is strictly >0.99.
# Correct-direction sites below that threshold remain grey.
# ============================================================

if (!exists(
    "mixture",
    envir = dr_env,
    inherits = FALSE
)) {
    stop(
        "D/R environment does not contain the full mixture table."
    )
}

dr_mix <- dr_env$mixture

required_dr_columns <- c(
    "chrom",
    "pos",
    "mixture_left_GT",
    "mixture_right_GT",
    "pair_posterior"
)

missing_dr_columns <- setdiff(
    required_dr_columns,
    names(dr_mix)
)

if (length(missing_dr_columns) > 0) {
    stop(
        paste(
            "D/R mixture table is missing:",
            paste(
                missing_dr_columns,
                collapse = ", "
            )
        )
    )
}

dr_mix$pos <- suppressWarnings(
    as.numeric(
        dr_mix$pos
    )
)

dr_mix$pair_posterior <- suppressWarnings(
    as.numeric(
        dr_mix$pair_posterior
    )
)

dr_correct <- dr_mix[
    dr_mix$chrom == chromosome &
    is.finite(
        dr_mix$pos
    ) &
    dr_mix$pos >= region_start &
    dr_mix$pos <= region_end &
    is.finite(
        dr_mix$pair_posterior
    ) &
    dr_mix$pair_posterior >= 0.95 &
    dr_mix$mixture_left_GT == "1/1" &
    dr_mix$mixture_right_GT == "0/1",
    ,
    drop = FALSE
]

dr_correct <- dr_correct[
    !duplicated(
        paste(
            dr_correct$chrom,
            dr_correct$pos,
            sep = ":"
        )
    ),
    ,
    drop = FALSE
]

dr_correct <- dr_correct[
    order(
        dr_correct$pos
    ),
    ,
    drop = FALSE
]

if (nrow(dr_correct) == 0) {
    stop(
        "No broad-region D/R P/P -> P/A sites survived."
    )
}

dr_correct$position_Mb <- (
    dr_correct$pos /
    1e6
)


# ------------------------------------------------------------
# Add Dennis and Ray BA independently.
# ------------------------------------------------------------

if (
    !exists(
        "dennis",
        envir = dr_env,
        inherits = FALSE
    ) ||
    !exists(
        "ray",
        envir = dr_env,
        inherits = FALSE
    )
) {
    stop(
        "D/R environment does not contain Dennis and Ray BA tables."
    )
}

dennis_ba <- dr_env$dennis
ray_ba <- dr_env$ray

required_ba_columns <- c(
    "chrom",
    "pos",
    "balanced_accuracy"
)

for (
    object_name
    in c(
        "dennis_ba",
        "ray_ba"
    )
) {

    object <- get(
        object_name
    )

    missing_columns <- setdiff(
        required_ba_columns,
        names(object)
    )

    if (length(missing_columns) > 0) {
        stop(
            paste(
                object_name,
                "is missing:",
                paste(
                    missing_columns,
                    collapse = ", "
                )
            )
        )
    }

    object$pos <- suppressWarnings(
        as.numeric(
            object$pos
        )
    )

    object$balanced_accuracy <- suppressWarnings(
        as.numeric(
            object$balanced_accuracy
        )
    )

    object$position_key <- paste(
        object$chrom,
        object$pos,
        sep = ":"
    )

    if (anyDuplicated(
        object$position_key
    )) {
        stop(
            paste(
                object_name,
                "contains duplicated chromosome-position keys;",
                "refusing to guess which BA belongs to a site."
            )
        )
    }

    assign(
        object_name,
        object
    )
}

dr_correct$position_key <- paste(
    dr_correct$chrom,
    dr_correct$pos,
    sep = ":"
)

dennis_index <- match(
    dr_correct$position_key,
    dennis_ba$position_key
)

ray_index <- match(
    dr_correct$position_key,
    ray_ba$position_key
)

dr_correct$Dennis_BA <- (
    dennis_ba$balanced_accuracy[
        dennis_index
    ]
)

dr_correct$Ray_BA <- (
    ray_ba$balanced_accuracy[
        ray_index
    ]
)


# Strict threshold: > 0.99, not >= 0.99.
dr_correct$BA_gt_099 <- (
    (
        is.finite(
            dr_correct$Dennis_BA
        ) &
        dr_correct$Dennis_BA > 0.99
    ) |
    (
        is.finite(
            dr_correct$Ray_BA
        ) &
        dr_correct$Ray_BA > 0.99
    )
)

dennis_for_max <- ifelse(
    is.finite(
        dr_correct$Dennis_BA
    ),
    dr_correct$Dennis_BA,
    -Inf
)

ray_for_max <- ifelse(
    is.finite(
        dr_correct$Ray_BA
    ),
    dr_correct$Ray_BA,
    -Inf
)

dr_correct$maximum_DR_BA <- pmax(
    dennis_for_max,
    ray_for_max
)

dr_correct$maximum_DR_BA[
    !is.finite(
        dr_correct$maximum_DR_BA
    )
] <- NA_real_


dr_long <- rbind(
    data.frame(
        position_Mb = dr_correct$position_Mb,
        side = "TN_Left",
        y = 1,
        biological_GT = "P/P",
        BA_gt_099 = dr_correct$BA_gt_099,
        stringsAsFactors = FALSE
    ),
    data.frame(
        position_Mb = dr_correct$position_Mb,
        side = "TN_Right",
        y = 0,
        biological_GT = "P/A",
        BA_gt_099 = dr_correct$BA_gt_099,
        stringsAsFactors = FALSE
    )
)
'''

text = replace_section(
    text,
    "# D/R RIBBON",
    "# BAND RIBBON",
    dr_block,
)


# ============================================================
# BAND:
# correct direction first; BA > 0.99 is highlight status.
# ============================================================

band_block = r'''
# BAND RIBBON
#
# Correct focal inheritance direction:
#
#     TN_Left  = biological A/A
#     TN_Right = biological P/A
#
# All such sites remain visible.  Those with Band BA >0.99
# are highlighted; the remainder are grey.
# ============================================================

if (!exists(
    "d",
    envir = band_env,
    inherits = FALSE
)) {
    stop(
        "Band environment does not contain d."
    )
}

band_all <- band_env$d

required_band_columns <- c(
    "chrom",
    "pos",
    "balanced_accuracy",
    "adjusted_confident",
    "exact_expected_LR",
    "left_biological_GT",
    "right_biological_GT"
)

missing_band_columns <- setdiff(
    required_band_columns,
    names(band_all)
)

if (length(missing_band_columns) > 0) {
    stop(
        paste(
            "Band table is missing:",
            paste(
                missing_band_columns,
                collapse = ", "
            )
        )
    )
}

band_all$pos <- suppressWarnings(
    as.numeric(
        band_all$pos
    )
)

band_all$balanced_accuracy <- suppressWarnings(
    as.numeric(
        band_all$balanced_accuracy
    )
)

band_correct <- band_all[
    band_all$chrom == chromosome &
    is.finite(
        band_all$pos
    ) &
    band_all$pos >= region_start &
    band_all$pos <= region_end &
    band_all$adjusted_confident %in% TRUE &
    band_all$exact_expected_LR %in% TRUE &
    band_all$left_biological_GT == "A/A" &
    band_all$right_biological_GT == "P/A",
    ,
    drop = FALSE
]

band_correct <- band_correct[
    !duplicated(
        paste(
            band_correct$chrom,
            band_correct$pos,
            sep = ":"
        )
    ),
    ,
    drop = FALSE
]

band_correct <- band_correct[
    order(
        band_correct$pos
    ),
    ,
    drop = FALSE
]

if (nrow(band_correct) == 0) {
    stop(
        "No broad-region Band A/A -> P/A sites survived."
    )
}

band_correct$position_Mb <- (
    band_correct$pos /
    1e6
)

# Strict threshold: >0.99.
band_correct$BA_gt_099 <- (
    is.finite(
        band_correct$balanced_accuracy
    ) &
    band_correct$balanced_accuracy > 0.99
)

band_long <- rbind(
    data.frame(
        position_Mb = band_correct$position_Mb,
        side = "TN_Left",
        y = 1,
        biological_GT = "A/A",
        BA_gt_099 = band_correct$BA_gt_099,
        stringsAsFactors = FALSE
    ),
    data.frame(
        position_Mb = band_correct$position_Mb,
        side = "TN_Right",
        y = 0,
        biological_GT = "P/A",
        BA_gt_099 = band_correct$BA_gt_099,
        stringsAsFactors = FALSE
    )
)
'''

text = replace_section(
    text,
    "# BAND RIBBON",
    "# GENOTYPE PALETTE",
    band_block,
)


# ============================================================
# Replace ribbon plotting:
# grey = correct direction but BA <=0.99 / unavailable
# colour = correct direction AND BA >0.99
# ============================================================

plot_block = r'''
# SIMPLE CORRECT-PATTERN RIBBON
#
# Every site shown already has the correct focal inheritance
# direction.
#
#   coloured = comparative BA > 0.99
#   grey     = correct direction but BA <= 0.99 / unavailable
#
# Thus BA is visual emphasis, NOT the first-stage site filter.
# ============================================================

make_correct_pattern_ribbon <- function(
    wide_data,
    long_data,
    title,
    subtitle,
    caption
) {

    grey_wide <- wide_data[
        !wide_data$BA_gt_099,
        ,
        drop = FALSE
    ]

    high_wide <- wide_data[
        wide_data$BA_gt_099,
        ,
        drop = FALSE
    ]

    grey_long <- long_data[
        !long_data$BA_gt_099,
        ,
        drop = FALSE
    ]

    high_long <- long_data[
        long_data$BA_gt_099,
        ,
        drop = FALSE
    ]

    ggplot() +

        # Correct-direction sites below BA threshold.
        geom_segment(
            data = grey_wide,
            aes(
                x = position_Mb,
                xend = position_Mb,
                y = 0,
                yend = 1
            ),
            colour = "#BDBDBD",
            alpha = 0.34,
            linewidth = 0.28
        ) +

        geom_point(
            data = grey_long,
            aes(
                x = position_Mb,
                y = y
            ),
            colour = "#BDBDBD",
            shape = 124,
            size = 2.05,
            alpha = 0.48
        ) +

        # Correct-direction sites with BA >0.99.
        geom_segment(
            data = high_wide,
            aes(
                x = position_Mb,
                xend = position_Mb,
                y = 0,
                yend = 1
            ),
            colour = "#009E73",
            alpha = 0.72,
            linewidth = 0.42
        ) +

        geom_point(
            data = high_long,
            aes(
                x = position_Mb,
                y = y,
                colour = biological_GT
            ),
            shape = 124,
            size = 2.55,
            alpha = 0.96
        ) +

        scale_colour_manual(
            values = genotype_colours,
            drop = FALSE,
            name = "Biological genotype"
        ) +

        scale_x_continuous(
            name = "Position across the broad Optix interval (Mb)",
            limits = x_limits,
            breaks = x_breaks,
            labels = sprintf(
                "%.2f",
                x_breaks
            ),
            expand = expansion(
                mult = c(
                    0.002,
                    0.002
                )
            )
        ) +

        scale_y_continuous(
            name = NULL,
            breaks = c(
                0,
                1
            ),
            labels = c(
                "TN_Right",
                "TN_Left"
            ),
            limits = c(
                -0.30,
                1.30
            )
        ) +

        labs(
            title = title,
            subtitle = subtitle,
            caption = caption
        ) +

        theme_minimal(
            base_size = 10
        ) +

        theme(
            plot.title = element_text(
                face = "bold",
                size = 11
            ),
            plot.subtitle = element_text(
                size = 8.5,
                colour = "#444444",
                margin = margin(
                    b = 5
                )
            ),
            plot.caption = element_text(
                size = 7.8,
                colour = "#555555",
                hjust = 0,
                margin = margin(
                    t = 5
                )
            ),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_line(
                colour = "#EEEEEE",
                linewidth = 0.3
            ),
            panel.grid.major.y = element_line(
                colour = "#E0E0E0",
                linewidth = 0.35
            ),
            axis.text.y = element_text(
                face = "bold"
            ),
            legend.position = "top",
            legend.justification = "left",
            legend.title = element_text(
                face = "bold"
            ),
            plot.margin = margin(
                2,
                16,
                9,
                10
            )
        )
}


dr_ribbon <- make_correct_pattern_ribbon(
    wide_data = dr_correct,
    long_data = dr_long,
    title = "Mixture-adjusted D/R-compatible left-right genotypes",
    subtitle = paste0(
        nrow(dr_correct),
        " correct-direction SNPs across Hmr18:1.060-1.240 Mb; ",
        sum(dr_correct$BA_gt_099),
        " have Dennis or Ray BA > 0.99."
    ),
    caption = paste(
        "Only TN_Left=P/P and TN_Right=P/A sites are shown.",
        "Sites with Dennis or Ray BA >0.99 are coloured;",
        "other correct-direction sites are retained in grey."
    )
)


band_ribbon <- make_correct_pattern_ribbon(
    wide_data = band_correct,
    long_data = band_long,
    title = "Mixture-adjusted Band-compatible left-right genotypes",
    subtitle = paste0(
        nrow(band_correct),
        " correct-direction SNPs across Hmr18:1.060-1.240 Mb; ",
        sum(band_correct$BA_gt_099),
        " have Band BA > 0.99."
    ),
    caption = paste(
        "Only TN_Left=A/A and TN_Right=P/A sites are shown.",
        "Sites with Band BA >0.99 are coloured;",
        "other correct-direction sites are retained in grey."
    )
)
'''

text = replace_section(
    text,
    "# SIMPLE CORRECT-PATTERN RIBBON",
    "# ALIGN GROBS",
    plot_block,
)


# ============================================================
# Rename version/output files.
# ============================================================

text = text.replace(
    "Wallbank_Dennis_Ray_plus_restricted_Band_correct_pattern_ribbons.pdf",
    "Wallbank_Dennis_Ray_plus_restricted_Band_correct_pattern_BAgt099_ribbons.pdf",
)

text = text.replace(
    "DR_broad_correct_pattern_ribbon_sites.tsv",
    "DR_broad_correct_pattern_BAgt099_ribbon_sites.tsv",
)

text = text.replace(
    "Band_broad_correct_pattern_ribbon_sites.tsv",
    "Band_broad_correct_pattern_BAgt099_ribbon_sites.tsv",
)

text = text.replace(
    "ANALYSIS 145 FIGURE COMPLETE",
    "ANALYSIS 146 FIGURE COMPLETE",
)


# Add high-BA counts to terminal summary.
old = '''cat(
    "D/R correct-pattern broad-region SNPs: ",
    nrow(dr_correct),
    "\\n",
    sep = ""
)

cat(
    "Band correct-pattern broad-region SNPs: ",
    nrow(band_correct),
    "\\n\\n",
    sep = ""
)
'''

new = '''cat(
    "D/R correct-direction SNPs: ",
    nrow(dr_correct),
    "\\n",
    "D/R correct-direction with Dennis or Ray BA >0.99: ",
    sum(dr_correct$BA_gt_099),
    "\\n\\n",
    sep = ""
)

cat(
    "Band correct-direction SNPs: ",
    nrow(band_correct),
    "\\n",
    "Band correct-direction with BA >0.99: ",
    sum(band_correct$BA_gt_099),
    "\\n\\n",
    sep = ""
)
'''

if old not in text:
    raise RuntimeError(
        "Could not find the Analysis-145 terminal summary block."
    )

text = text.replace(
    old,
    new,
)


target.write_text(
    text
)

print(
    f"Created: {target}"
)
