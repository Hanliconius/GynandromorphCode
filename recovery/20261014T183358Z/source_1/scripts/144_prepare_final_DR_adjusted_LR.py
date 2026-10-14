#!/usr/bin/env python3

from pathlib import Path

PROJECT = Path(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling"
)

SOURCE = (
    PROJECT
    / "scripts"
    / "115_plot_D_R_Wallbank_fixed_dominant_evidence.R"
)

OUTPUT = (
    PROJECT
    / "scripts"
    / "115b_plot_D_R_Wallbank_fixed_dominant_adjusted_LR.R"
)


def replace_once(text, old, new):
    count = text.count(old)

    if count != 1:
        raise RuntimeError(
            f"Expected one occurrence of {old!r}; found {count}"
        )

    return text.replace(old, new, 1)


def main():
    text = SOURCE.read_text()

    start_marker = (
        "# One physical copy of every SNP for the focal genotype tracks."
    )

    end_marker = (
        "# Add comparative BA values to the focal-site table."
    )

    if text.count(start_marker) != 1:
        raise RuntimeError(
            "Could not uniquely locate the start of the focal-call block"
        )

    if text.count(end_marker) != 1:
        raise RuntimeError(
            "Could not uniquely locate the end of the focal-call block"
        )

    start = text.index(start_marker)
    end = text.index(end_marker, start)

    adjusted_block = r'''# One physical copy of every SNP for the focal genotype tracks.
# Comparative SNP selection and BA values remain unchanged.
# Focal genotypes are replaced with reciprocal-mixture calls here.
focal <- dennis[
    !duplicated(dennis$site_key),
    ,
    drop = FALSE
]

mixture_file <- file.path(
    project_dir,
    "analysis",
    "original_502",
    "reciprocal_L_R_mixture",
    "optix_reciprocal_L_R_genotype_pairs.posterior95.tsv"
)

if (!file.exists(mixture_file)) {
    stop(
        paste(
            "Missing reciprocal-mixture table:",
            mixture_file
        )
    )
}

mixture <- read.delim(
    mixture_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

mixture_required <- c(
    "chrom",
    "pos",
    "ref",
    "alt",
    "mixture_left_GT",
    "mixture_right_GT",
    "pair_posterior",
    "posterior95"
)

mixture_missing <- setdiff(
    mixture_required,
    names(mixture)
)

if (length(mixture_missing) > 0) {
    stop(
        paste(
            "Reciprocal-mixture table is missing columns:",
            paste(
                mixture_missing,
                collapse = ", "
            )
        )
    )
}

mixture$pos <- as.integer(
    mixture$pos
)

mixture$site_key <- paste(
    mixture$chrom,
    mixture$pos,
    mixture$ref,
    mixture$alt,
    sep = ":"
)

if (anyDuplicated(mixture$site_key)) {
    stop(
        "Reciprocal-mixture table contains duplicate physical SNPs"
    )
}

mixture$pair_posterior <- as.numeric(
    mixture$pair_posterior
)

mixture$posterior95_bool <- as_bool(
    mixture$posterior95
)

mixture$adjusted_confident <- (
    mixture$posterior95_bool
    & is.finite(mixture$pair_posterior)
    & mixture$pair_posterior >= 0.95
)

mixture_index <- match(
    focal$site_key,
    mixture$site_key
)

focal$mixture_pair_posterior <- mixture$pair_posterior[
    mixture_index
]

focal$adjusted_confident <- mixture$adjusted_confident[
    mixture_index
]

focal$adjusted_confident[
    is.na(focal$adjusted_confident)
] <- FALSE

focal$left_GT <- canonical_gt(
    mixture$mixture_left_GT[
        mixture_index
    ]
)

focal$right_GT <- canonical_gt(
    mixture$mixture_right_GT[
        mixture_index
    ]
)

focal$left_QC <- (
    focal$adjusted_confident
    & !is.na(focal$left_GT)
)

focal$right_QC <- (
    focal$adjusted_confident
    & !is.na(focal$right_GT)
)

focal$both_QC <- (
    focal$left_QC
    & focal$right_QC
)

focal$same_GT <- (
    focal$both_QC
    & focal$left_GT == focal$right_GT
)

focal$different_GT <- (
    focal$both_QC
    & focal$left_GT != focal$right_GT
)

focal$exact_expected_LR <- (
    focal$both_QC
    & focal$left_GT == "1/1"
    & focal$right_GT == "0/1"
)

focal$pair_class <- (
    "One or both adjusted calls unavailable"
)

focal$pair_class[
    focal$same_GT
] <- "Same genotype"

focal$pair_class[
    focal$different_GT
] <- "Other left-right difference"

focal$pair_class[
    focal$exact_expected_LR
] <- "Exact P/P to P/A"

focal$focal_call_source <- ifelse(
    focal$adjusted_confident,
    "reciprocal_mixture_posterior95_only",
    "no_posterior95_call"
)

'''

    text = (
        text[:start]
        + adjusted_block
        + text[end:]
    )

    replacements = {
        (
            "Comparative association and pooled focal genotypes"
        ): (
            "Comparative association and mixture-adjusted focal genotypes"
        ),
        (
            "Pooled left and right genotypes at comparative-selected sites"
        ): (
            "Mixture-adjusted left and right genotypes "
            "at comparative-selected sites"
        ),
        (
            '"All 19 comparative-selected SNPs are shown.",'
        ): (
            '"All 19 comparative-selected SNPs are shown.",\n'
            '            "Reciprocal-mixture posterior >=0.95 only; '
            'no original-call fallback.",'
        ),
        (
            'name = "Pooled genotype"'
        ): (
            'name = "Mixture-adjusted genotype"'
        ),
        (
            "Genotype color is shown separately for each pooled side."
        ): (
            "Genotype color is shown separately for each focal side."
        ),
    }

    for old, new in replacements.items():
        text = replace_once(
            text,
            old,
            new
        )

    old_basename = (
        "broad_D_R_DP2_GQ20_BA_manhattan_"
        "Wallbank_fixed_dominant"
    )

    new_basename = (
        "broad_D_R_DP2_GQ20_BA_manhattan_"
        "Wallbank_fixed_dominant_adjusted_LR"
    )

    basename_count = text.count(
        old_basename
    )

    if basename_count != 2:
        raise RuntimeError(
            "Expected two D/R PDF/PNG basename occurrences; "
            f"found {basename_count}"
        )

    text = text.replace(
        old_basename,
        new_basename
    )

    OUTPUT.write_text(text)

    print(f"Prepared: {OUTPUT}")
    print(
        "Comparative Dennis/Ray BA and 19-SNP selection: unchanged"
    )
    print(
        "Focal calls: reciprocal mixture posterior >=0.95 only"
    )
    print(
        "Original focal GT fallback: disabled"
    )


if __name__ == "__main__":
    main()
