#!/usr/bin/env Rscript

source_file <- paste0(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/",
    "GynandroCalling/scripts/",
    "152_plot_population_frequency_plus_LR_BWPK_exact_BA.R"
)

target_file <- paste0(
    "/gpfs/automountdir/gpfs/scratch/martinlab/joe/",
    "GynandroCalling/scripts/",
    "153_plot_population_frequency_plus_LR_BWPK_exact_BA_nonfixed_ref.R"
)


count_fixed <- function(
    text,
    needle
) {

    hit <- gregexpr(
        needle,
        text,
        fixed = TRUE
    )[[1]]

    if (
        length(hit) == 1 &&
        hit[[1]] == -1
    ) {
        return(0L)
    }

    length(hit)
}


if (!file.exists(source_file)) {
    stop(
        paste(
            "Missing working Analysis-152 source:",
            source_file
        )
    )
}


text <- paste(
    readLines(
        source_file,
        warn = FALSE
    ),
    collapse = "\n"
)


# ============================================================
# NEW OUTPUT DIRECTORY
# ============================================================

old_output <- paste0(
    '"',
    "population_frequency_plus_LR_BWPK_exact_BA",
    '"'
)

new_output <- paste0(
    '"',
    "population_frequency_plus_LR_BWPK_exact_BA_nonfixed_ref",
    '"'
)

n_output <- count_fixed(
    text,
    old_output
)

if (n_output != 1L) {
    stop(
        paste(
            "Expected one Analysis-152 output-directory string;",
            "found",
            n_output
        )
    )
}

text <- sub(
    old_output,
    new_output,
    text,
    fixed = TRUE
)


# ============================================================
# INSERT FILTER BEFORE EXISTING COMBINATION STEP
# ============================================================

marker <- paste(
    "    # --------------------------------------------------------",
    "    # Combine population-frequency and individual rows.",
    "    # --------------------------------------------------------",
    sep = "\n"
)

n_marker <- count_fixed(
    text,
    marker
)

if (n_marker != 1L) {
    stop(
        paste(
            "Expected exactly one population/individual",
            "combination marker; found",
            n_marker
        )
    )
}


filter_code <- paste(
'    # ========================================================',
'    # REMOVE POSITIONS FIXED FOR RAW REFERENCE',
'    #',
'    # Scope:',
'    #   * exact comparative BA cohort for this pattern',
'    #   * TN_Left reciprocal-mixture row',
'    #   * TN_Right reciprocal-mixture row',
'    #   * all Butterfly World piano-key individuals',
'    #',
'    # Remove a SNP iff:',
'    #   * >=1 displayed genotype is callable',
'    #   * every callable displayed RAW genotype is 0/0',
'    #',
'    # N / missing / Uncertain calls are ignored.',
'    # An all-missing site is retained.',
'    # ========================================================',
'',
'    if (!("raw_GT" %in% names(d))) {',
'        stop(',
'            paste0(',
'                panel_name,',
'                ": raw_GT is absent from the established ",',
'                "genotype-plot input."',
'            )',
'        )',
'    }',
'',
'',
'    # --------------------------------------------------------',
'    # Select precisely the butterflies represented by this',
'    # figure.',
'    #',
'    # allowed_samples was already frozen by Analysis 152 to',
'    # the exact Band / Ray / Dennis BA cohort.',
'    # --------------------------------------------------------',
'',
'    in_display_scope <- (',
'        d$sample_id %in% allowed_samples |',
'        is_bwpk |',
'        is_left_recip |',
'        is_right_recip',
'    )',
'',
'    site_calls <- d[',
'        in_display_scope,',
'        c(',
'            "site_key",',
'            "display_GT",',
'            "raw_GT"',
'        ),',
'        drop = FALSE',
'    ]',
'',
'',
'    # --------------------------------------------------------',
'    # Callable means that the established genotype plot itself',
'    # displayed a resolved genotype.',
'    # --------------------------------------------------------',
'',
'    site_calls$callable <- (',
'        site_calls$display_GT %in%',
'            c(',
'                "A/A",',
'                "P/A",',
'                "P/P"',
'            )',
'    )',
'',
'',
'    # --------------------------------------------------------',
'    # Canonicalize RAW REF/ALT genotype.',
'    # --------------------------------------------------------',
'',
'    site_calls$raw_GT_canonical <- gsub(',
'        "\\\\|",',
'        "/",',
'        trimws(',
'            as.character(',
'                site_calls$raw_GT',
'            )',
'        )',
'    )',
'',
'    site_calls$raw_GT_canonical[',
'        site_calls$raw_GT_canonical == "1/0"',
'    ] <- "0/1"',
'',
'    site_calls$raw_GT_canonical[',
'        !(site_calls$raw_GT_canonical %in%',
'            c(',
'                "0/0",',
'                "0/1",',
'                "1/1"',
'            )',
'        )',
'    ] <- NA_character_',
'',
'',
'    site_calls$usable <- (',
'        site_calls$callable &',
'        !is.na(',
'            site_calls$raw_GT_canonical',
'        )',
'    )',
'',
'',
'    # --------------------------------------------------------',
'    # One fixed-reference classification per physical SNP.',
'    # --------------------------------------------------------',
'',
'    calls_by_site <- split(',
'        site_calls,',
'        site_calls$site_key',
'    )',
'',
'    fixed_ref_rows <- lapply(',
'        names(calls_by_site),',
'        function(site_key_value) {',
'',
'            z <- calls_by_site[[site_key_value]]',
'',
'            usable_gt <- z$raw_GT_canonical[',
'                z$usable',
'            ]',
'',
'            callable_n <- length(usable_gt)',
'',
'            ref_hom_n <- sum(',
'                usable_gt == "0/0"',
'            )',
'',
'            het_n <- sum(',
'                usable_gt == "0/1"',
'            )',
'',
'            alt_hom_n <- sum(',
'                usable_gt == "1/1"',
'            )',
'',
'            nonref_n <- (',
'                het_n +',
'                alt_hom_n',
'            )',
'',
'            data.frame(',
'                site_key = site_key_value,',
'                callable_n = callable_n,',
'                ref_hom_n = ref_hom_n,',
'                het_n = het_n,',
'                alt_hom_n = alt_hom_n,',
'                nonref_n = nonref_n,',
'',
'                fixed_reference = (',
'                    callable_n > 0 &&',
'                    nonref_n == 0',
'                ),',
'',
'                stringsAsFactors = FALSE',
'            )',
'        }',
'    )',
'',
'    fixed_ref_audit <- do.call(',
'        rbind,',
'        fixed_ref_rows',
'    )',
'',
'',
'    # --------------------------------------------------------',
'    # Join back to the complete original SNP list.',
'    # --------------------------------------------------------',
'',
'    fixed_ref_audit <- merge(',
'        sites[',
'            ,',
'            c(',
'                "site_key",',
'                "chrom",',
'                "pos",',
'                "ref",',
'                "alt"',
'            ),',
'            drop = FALSE',
'        ],',
'        fixed_ref_audit,',
'        by = "site_key",',
'        all.x = TRUE,',
'        sort = FALSE',
'    )',
'',
'    fixed_ref_audit$fixed_reference[',
'        is.na(',
'            fixed_ref_audit$fixed_reference',
'        )',
'    ] <- FALSE',
'',
'    count_columns <- c(',
'        "callable_n",',
'        "ref_hom_n",',
'        "het_n",',
'        "alt_hom_n",',
'        "nonref_n"',
'    )',
'',
'    for (count_column in count_columns) {',
'',
'        fixed_ref_audit[[count_column]][',
'            is.na(',
'                fixed_ref_audit[[count_column]]',
'            )',
'        ] <- 0L',
'    }',
'',
'    fixed_ref_audit <- fixed_ref_audit[',
'        order(',
'            fixed_ref_audit$pos,',
'            fixed_ref_audit$ref,',
'            fixed_ref_audit$alt',
'        ),',
'        ,',
'        drop = FALSE',
'    ]',
'',
'',
'    remove_site_keys <- fixed_ref_audit$site_key[',
'        fixed_ref_audit$fixed_reference',
'    ]',
'',
'    keep_site_keys <- fixed_ref_audit$site_key[',
'        !fixed_ref_audit$fixed_reference',
'    ]',
'',
'    original_site_n <- nrow(sites)',
'    removed_site_n <- length(remove_site_keys)',
'    retained_site_n <- length(keep_site_keys)',
'',
'    if (retained_site_n < 1L) {',
'        stop(',
'            paste0(',
'                panel_name,',
'                ": fixed-reference filtering removed every SNP."',
'            )',
'        )',
'    }',
'',
'',
'    # --------------------------------------------------------',
'    # Audit.',
'    # --------------------------------------------------------',
'',
'    fixed_ref_audit_file <- file.path(',
'        output_dir,',
'        paste0(',
'            "optix_",',
'            panel_name,',
'            "_fixed_reference_site_filter.tsv"',
'        )',
'    )',
'',
'    write.table(',
'        fixed_ref_audit,',
'        file = fixed_ref_audit_file,',
'        sep = "\\t",',
'        quote = FALSE,',
'        row.names = FALSE',
'    )',
'',
'',
'    # --------------------------------------------------------',
'    # Remove fixed-reference columns and re-index retained SNPs',
'    # consecutively so there are no blank visual gaps.',
'    # --------------------------------------------------------',
'',
'    sites <- sites[',
'        sites$site_key %in% keep_site_keys,',
'        ,',
'        drop = FALSE',
'    ]',
'',
'    sites <- sites[',
'        order(',
'            sites$pos,',
'            sites$ref,',
'            sites$alt',
'        ),',
'        ,',
'        drop = FALSE',
'    ]',
'',
'    sites$site_index <- seq_len(',
'        nrow(sites)',
'    )',
'',
'',
'    pop_grid <- pop_grid[',
'        pop_grid$site_key %in% keep_site_keys,',
'        ,',
'        drop = FALSE',
'    ]',
'',
'    pop_grid$site_index <- sites$site_index[',
'        match(',
'            pop_grid$site_key,',
'            sites$site_key',
'        )',
'    ]',
'',
'',
'    focal <- focal[',
'        focal$site_key %in% keep_site_keys,',
'        ,',
'        drop = FALSE',
'    ]',
'',
'    focal$site_index <- sites$site_index[',
'        match(',
'            focal$site_key,',
'            sites$site_key',
'        )',
'    ]',
'',
'',
'    bwpk <- bwpk[',
'        bwpk$site_key %in% keep_site_keys,',
'        ,',
'        drop = FALSE',
'    ]',
'',
'    bwpk$site_index <- sites$site_index[',
'        match(',
'            bwpk$site_key,',
'            sites$site_key',
'        )',
'    ]',
'',
'',
'    cat(',
'        sprintf(',
'            paste0(',
'                "%-8s | SNPs: %d original -> %d retained; ",',
'                "%d fixed-reference positions removed\\n"',
'            ),',
'            panel_name,',
'            original_site_n,',
'            retained_site_n,',
'            removed_site_n',
'        )',
'    )',
'',
    sep = "\n"
)


text <- sub(
    marker,
    paste0(
        filter_code,
        "\n\n",
        marker
    ),
    text,
    fixed = TRUE
)


# ============================================================
# WRITE ONLY AFTER BOTH PATCHES SUCCEEDED
# ============================================================

writeLines(
    text,
    con = target_file,
    useBytes = TRUE
)

cat(
    "Created:\n",
    target_file,
    "\n\n",
    sep = ""
)

cat(
    "Filter rule:\n",
    "  remove site iff >=1 callable displayed raw genotype exists\n",
    "  and every callable displayed raw genotype is 0/0.\n",
    sep = ""
)
