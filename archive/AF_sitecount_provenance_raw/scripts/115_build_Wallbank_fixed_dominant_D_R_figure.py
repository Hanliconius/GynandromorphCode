#!/usr/bin/env python3

from pathlib import Path
import re


PROJECT = Path("/scratch/martinlab/joe/GynandroCalling")

SOURCE = (
    PROJECT
    / "scripts"
    / "107c_plot_BA_with_separate_focal_genotype_tracks.R"
)

OUTPUT = (
    PROJECT
    / "scripts"
    / "115_plot_D_R_Wallbank_fixed_dominant_evidence.R"
)


def replace_once(text, old, new, label):
    count = text.count(old)

    if count != 1:
        raise RuntimeError(
            f"{label}: expected exactly one match, found {count}"
        )

    return text.replace(old, new, 1)


if not SOURCE.is_file():
    raise RuntimeError(f"Missing source script: {SOURCE}")

text = SOURCE.read_text()

if "OPTIX_ZOOM_1060_1240_AND_GENE_ANNOTATION" not in text:
    raise RuntimeError(
        "The source script does not contain the expected "
        "1.060-1.240 Mb zoom/optix annotation patch."
    )


# ------------------------------------------------------------
# 1. Load and define the comparative-only evidence sites.
# ------------------------------------------------------------

selection_anchor = (
    "# Long representation for separate left and right genotype tracks."
)

selection_block = r'''
# WALLBANK_FIXED_DOMINANT_COMPARATIVE_SELECTION
#
# These positions are selected entirely from comparative butterflies.
# TN_Left and TN_Right are inspected only after this selection.

fixed_replication_dir <- file.path(
    "/scratch/martinlab/joe/GynandroCalling",
    "analysis",
    "original_502",
    "fixed_SNP_replication"
)

dennis_fixed_file <- file.path(
    fixed_replication_dir,
    "Wallbank2016_fixed_SNP_replication_sites.tsv"
)

ray_dominant_file <- file.path(
    fixed_replication_dir,
    "Ray_contigua_dominant_pattern_sites.tsv"
)

if (!file.exists(dennis_fixed_file)) {
    stop(
        paste(
            "Missing Dennis fixed-SNP table:",
            dennis_fixed_file
        )
    )
}

if (!file.exists(ray_dominant_file)) {
    stop(
        paste(
            "Missing Ray dominant-pattern table:",
            ray_dominant_file
        )
    )
}

dennis_fixed_all <- read.delim(
    dennis_fixed_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

ray_dominant_all <- read.delim(
    ray_dominant_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
)

required_dennis_columns <- c(
    "analysis",
    "scheme",
    "chrom",
    "pos",
    "ref",
    "alt",
    "fixed_complete"
)

missing_dennis_columns <- setdiff(
    required_dennis_columns,
    names(dennis_fixed_all)
)

if (length(missing_dennis_columns) > 0) {
    stop(
        paste(
            "Dennis table lacks columns:",
            paste(missing_dennis_columns, collapse = ", ")
        )
    )
}

required_ray_columns <- c(
    "scheme",
    "chrom",
    "pos",
    "ref",
    "alt",
    "dominant_pattern_among_called_min2"
)

missing_ray_columns <- setdiff(
    required_ray_columns,
    names(ray_dominant_all)
)

if (length(missing_ray_columns) > 0) {
    stop(
        paste(
            "Ray table lacks columns:",
            paste(missing_ray_columns, collapse = ", ")
        )
    )
}

dennis_fixed_all$pos <- as.numeric(
    dennis_fixed_all$pos
)

ray_dominant_all$pos <- as.numeric(
    ray_dominant_all$pos
)

# The 16 paper-QC, complete fixed SNPs in the best Dennis
# 5-kb window reported by analysis 112.
dennis_selected <- dennis_fixed_all[
    dennis_fixed_all$analysis
        == "Dennis_exact_specimen_overlap"
    & dennis_fixed_all$scheme == "paper_QC"
    & as_bool(dennis_fixed_all$fixed_complete)
    & dennis_fixed_all$pos >= 1215000
    & dennis_fixed_all$pos <= 1219999,
    ,
    drop = FALSE
]

if (nrow(dennis_selected) != 16) {
    stop(
        paste(
            "Expected 16 complete Dennis fixed SNPs in",
            "Hmr18:1215000-1219999; recovered",
            nrow(dennis_selected)
        )
    )
}

# The three N-tolerant dominant-compatible sites in the
# broad Ray neighbourhood selected without focal data.
ray_target_positions <- c(
    1162302,
    1167014,
    1170389
)

ray_selected <- ray_dominant_all[
    ray_dominant_all$scheme == "paper_QC"
    & ray_dominant_all$pos %in% ray_target_positions
    & as_bool(
        ray_dominant_all$dominant_pattern_among_called_min2
    ),
    ,
    drop = FALSE
]

if (
    nrow(ray_selected) != 3
    || !setequal(
        ray_selected$pos,
        ray_target_positions
    )
) {
    stop(
        paste(
            "Expected the three Ray sites",
            paste(ray_target_positions, collapse = ", "),
            "but recovered",
            paste(sort(ray_selected$pos), collapse = ", ")
        )
    )
}

selected_evidence <- rbind(
    data.frame(
        module = "Dennis",
        chrom = dennis_selected$chrom,
        pos = dennis_selected$pos,
        ref = dennis_selected$ref,
        alt = dennis_selected$alt,
        evidence_class = (
            "Complete fixed SNP: Dennis 5-kb peak"
        ),
        stringsAsFactors = FALSE
    ),
    data.frame(
        module = "Ray",
        chrom = ray_selected$chrom,
        pos = ray_selected$pos,
        ref = ray_selected$ref,
        alt = ray_selected$alt,
        evidence_class = (
            "Dominant-compatible among callable Ray specimens"
        ),
        stringsAsFactors = FALSE
    )
)

selected_evidence$selection_key <- paste(
    selected_evidence$chrom,
    selected_evidence$pos,
    selected_evidence$ref,
    selected_evidence$alt,
    sep = ":"
)

if (anyDuplicated(selected_evidence$selection_key)) {
    stop("The selected evidence table contains duplicate physical SNPs.")
}

focal$selection_key <- paste(
    focal$chrom,
    focal$pos,
    focal$ref,
    focal$alt,
    sep = ":"
)

selection_index <- match(
    focal$selection_key,
    selected_evidence$selection_key
)

focal$evidence_module <- selected_evidence$module[
    selection_index
]

focal$evidence_class <- selected_evidence$evidence_class[
    selection_index
]

focal <- focal[
    !is.na(focal$evidence_module),
    ,
    drop = FALSE
]

if (nrow(focal) != 19) {
    stop(
        paste(
            "Expected focal calls at 19 selected SNPs;",
            "recovered",
            nrow(focal)
        )
    )
}

selected_evidence_file <- file.path(
    "/scratch/martinlab/joe/GynandroCalling",
    "analysis",
    "original_502",
    "fixed_SNP_replication",
    "D_R_Wallbank_fixed_dominant_figure_sites.tsv"
)

write.table(
    selected_evidence[
        ,
        c(
            "module",
            "chrom",
            "pos",
            "ref",
            "alt",
            "evidence_class"
        )
    ],
    file = selected_evidence_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


'''

text = replace_once(
    text,
    selection_anchor,
    selection_block + selection_anchor,
    "comparative-site selection insertion"
)


# ------------------------------------------------------------
# 2. Show every selected focal genotype pair, rather than only
#    pairs already matching the expected focal transition.
# ------------------------------------------------------------

old_different_rows = '''# Display only the exact expected Dennis/Ray focal state:
# TN_Left = P/P and TN_Right = P/A.
different_rows <- focal_long[
    focal_long$exact_expected_LR,
    ,
    drop = FALSE
]
'''

new_different_rows = '''# Display every non-identical, QC-passing genotype pair at the
# independently selected comparative sites. Exact focal matches
# are highlighted separately below.
different_rows <- focal_long[
    focal_long$pair_class %in% c(
        "Other left-right difference",
        "Exact P/P to P/A"
    ),
    ,
    drop = FALSE
]
'''

text = replace_once(
    text,
    old_different_rows,
    new_different_rows,
    "focal-row display rule"
)


old_exact_segments = '''exact_segments <- focal[
    focal$exact_expected_LR,
    ,
    drop = FALSE
]
'''

new_exact_segments = '''other_difference_segments <- focal[
    focal$both_QC
    & focal$different_GT
    & !focal$exact_expected_LR,
    ,
    drop = FALSE
]

exact_segments <- focal[
    focal$exact_expected_LR,
    ,
    drop = FALSE
]
'''

text = replace_once(
    text,
    old_exact_segments,
    new_exact_segments,
    "focal segment classification"
)


# ------------------------------------------------------------
# 3. Make the upper-panel markers represent the new
#    comparative-only evidence sets.
# ------------------------------------------------------------

ba_selected_block = r'''
# Match the comparative-only selected sites to their appropriate
# Dennis or Ray BA trace. Focal genotypes are not involved.
ba_selection_key <- paste(
    ba_data$module,
    ba_data$chrom,
    ba_data$pos,
    ba_data$ref,
    ba_data$alt,
    sep = ":"
)

selected_ba_key <- paste(
    selected_evidence$module,
    selected_evidence$chrom,
    selected_evidence$pos,
    selected_evidence$ref,
    selected_evidence$alt,
    sep = ":"
)

ba_selected <- ba_data[
    ba_selection_key %in% selected_ba_key,
    ,
    drop = FALSE
]

if (nrow(ba_selected) != 19) {
    stop(
        paste(
            "Expected 19 selected points in the BA traces;",
            "recovered",
            nrow(ba_selected)
        )
    )
}


'''

zoom_anchor = "# OPTIX_ZOOM_1060_1240_AND_GENE_ANNOTATION"

text = replace_once(
    text,
    zoom_anchor,
    ba_selected_block + zoom_anchor,
    "BA marker-data insertion"
)


marker_pattern = re.compile(
    r'''    geom_point\(
        data = ba_perfect,
        aes\(
            x = position_Mb,
            y = balanced_accuracy
        \),
        shape = 21,
        fill = "white",
        colour = "#111111",
        size = 1\.7,
(?P<tail>.*?)    \) \+''',
    flags=re.DOTALL
)

marker_match = marker_pattern.search(text)

if marker_match is None:
    raise RuntimeError(
        "Could not locate the existing ba_perfect marker layer."
    )

new_marker_layer = '''    geom_point(
        data = ba_selected,
        aes(
            x = position_Mb,
            y = balanced_accuracy
        ),
        shape = 25,
        fill = "#111111",
        colour = "#111111",
        size = 2.25,
        stroke = 0.25
    ) +'''

text = (
    text[:marker_match.start()]
    + new_marker_layer
    + text[marker_match.end():]
)


# Add neutral connectors for non-exact left/right differences.
exact_segment_anchor = '''    geom_segment(
        data = exact_segments,
'''

other_segment_layer = '''    geom_segment(
        data = other_difference_segments,
        aes(
            x = position_Mb,
            xend = position_Mb,
            y = 0,
            yend = 1
        ),
        colour = "#777777",
        alpha = 0.34,
        linewidth = 0.28
    ) +

'''

text = replace_once(
    text,
    exact_segment_anchor,
    other_segment_layer + exact_segment_anchor,
    "non-exact focal connector layer"
)


# ------------------------------------------------------------
# 4. Clarify the figure text.
# ------------------------------------------------------------

old_ba_subtitle = '''        subtitle = paste(
            "Hmr18:1,060,000-1,240,000.",
            "Comparative balanced accuracy was calculated without",
            "TN_Left or TN_Right; focal genotypes are displayed below."
        )
'''

new_ba_subtitle = '''        subtitle = paste(
            "Hmr18:1,060,000-1,240,000.",
            "Triangles mark 16 complete Dennis fixed SNPs and",
            "3 Ray-region dominant-compatible SNPs.",
            "TN_Left and TN_Right were excluded from SNP selection."
        )
'''

text = replace_once(
    text,
    old_ba_subtitle,
    new_ba_subtitle,
    "upper-panel subtitle"
)


old_focal_labs = '''        title = "Pooled genotypes at exact expected left-right sites",
        subtitle = paste(
            "Only positions with TN_Left=P/P and TN_Right=P/A are shown.",
            "SNP selection and balanced accuracy use comparative samples only."
        ),
        caption = paste(
            "All other genotype combinations, identical calls and",
            "QC-failed positions are omitted from this display.",
'''

new_focal_labs = '''        title = "Pooled left and right genotypes at comparative-selected sites",
        subtitle = paste(
            "All 19 comparative-selected SNPs are shown.",
            "Green connectors mark TN_Left=P/P and TN_Right=P/A;",
            "gray connectors mark other called left-right differences."
        ),
        caption = paste(
            "Dennis: 16 complete fixed SNPs in the best 5-kb window.",
            "Ray: 3 dominant-compatible sites in the Ray region;",
            "missing or low-quality comparative calls were tolerated.",
'''

text = replace_once(
    text,
    old_focal_labs,
    new_focal_labs,
    "lower-panel labels"
)


# ------------------------------------------------------------
# 5. Give every output after pdf_file a new name, protecting
#    the frozen figure and its summaries.
# ------------------------------------------------------------

output_anchor = "pdf_file <- file.path("

if text.count(output_anchor) != 1:
    raise RuntimeError(
        "Expected one pdf_file definition in the plotting script."
    )

prefix, suffix = text.split(output_anchor, 1)

suffix = re.sub(
    r'"(broad_D_R_[^"\n]+?)(\.(?:pdf|png|tsv))"',
    r'"\1_Wallbank_fixed_dominant\2"',
    suffix
)

text = prefix + output_anchor + suffix


# Add an explicit version tag to console output.
text = text.replace(
    'cat(\n    "BROAD D/R BA MANHATTAN WITH FOCAL L/R TRACK',
    'cat(\n    "BROAD D/R WALLBANK FIXED/DOMINANT EVIDENCE FIGURE',
    1
)


OUTPUT.write_text(text)

print(f"Source preserved: {SOURCE}")
print(f"New plotting script: {OUTPUT}")
print()
print("Comparative evidence encoded:")
print("Dennis: 16 complete fixed SNPs in Hmr18:1215000-1219999")
print("Ray: Hmr18:1162302, Hmr18:1167014, Hmr18:1170389")
print("Focal TN_Left/TN_Right genotypes are display-only.")
