from pathlib import Path

src = Path("scripts/169_plot_science_exact_L_DRDR_R_BDR.R")
dst = Path("scripts/170_plot_science_exact_L_DRDR_R_BDR_homhetcolours.R")

if not src.exists():
    raise RuntimeError(f"Missing source script: {src}")

text = src.read_text()


# ============================================================
# 1. Preserve 169; create separate 170 outputs
# ============================================================

text = text.replace(
    '"science_exact_L_DRDR_R_BDR"',
    '"science_exact_L_DRDR_R_BDR_homhetcolours"'
)

text = text.replace("Science_169", "Science_170")
text = text.replace("169A", "170A")
text = text.replace("169B", "170B")
text = text.replace("ANALYSIS 169", "ANALYSIS 170")


# ============================================================
# 2. Exact supplement genotype colours
#
# Established genoplot scheme:
#   A/A  #2C9BCB
#   P/A  #F28E2B
#   P/P  #B61F72
#
# For THIS focal model:
#   Left  DR/DR = P/P -> magenta
#   Right B/DR  = P/A -> orange
# ============================================================

old = '''B_colour <- "#0072B2"

DR_colour <- "#CC79A7"
'''

new = '''# Supplement genoplot genotype colours
HOM_AA_colour <- "#2C9BCB"
HET_colour    <- "#F28E2B"
HOM_PP_colour <- "#B61F72"

# Focal biological states used here
DR_colour <- HOM_PP_colour
'''

if old not in text:
    raise RuntimeError("Could not locate original lower-panel colour definitions.")

text = text.replace(old, new)


# ============================================================
# 3. Centre the Right tick at the actual SNP coordinate.
#
# 169 used +/- pair_offset to make two adjacent dashes.
# 170 uses ONE heterozygote tick, so offset = 0.
# ============================================================

if "pair_offset <- 0.00033" not in text:
    raise RuntimeError("Could not locate pair_offset in 169.")

text = text.replace(
    "pair_offset <- 0.00033",
    "pair_offset <- 0"
)


# ============================================================
# 4. Modify only the four RIGHT plotting layers.
#
# Keep *_right_B as the visible layer, but recolour it as
# heterozygote orange.
#
# Make the redundant *_right_DR layer invisible.
#
# This avoids brittle whole-block replacement.
# ============================================================

lines = text.splitlines()

visible_targets = {
    "all_right_B",
    "high_right_B",
}

hidden_targets = {
    "all_right_DR",
    "high_right_DR",
}

current_target = None
changed_colour = []
changed_alpha = []

for i, line in enumerate(lines):

    stripped = line.strip()

    # Identify which geom's data argument we just entered.
    for target in visible_targets | hidden_targets:
        if stripped == f"data = {target},":
            current_target = target
            break

    if current_target is None:
        continue

    # Visible RIGHT layer: orange heterozygote
    if current_target in visible_targets and stripped.startswith("colour = "):
        indent = line[:len(line) - len(line.lstrip())]
        lines[i] = f"{indent}colour = HET_colour,"
        changed_colour.append(current_target)

    # Hidden duplicate RIGHT layer
    if current_target in hidden_targets and stripped.startswith("alpha = "):
        indent = line[:len(line) - len(line.lstrip())]
        comma = "," if stripped.endswith(",") else ""
        lines[i] = f"{indent}alpha = 0{comma}"
        changed_alpha.append(current_target)
        current_target = None

    # Once visible colour was replaced, keep scanning until its alpha
    # merely so we know we've passed that geom.
    elif current_target in visible_targets and stripped.startswith("alpha = "):
        current_target = None


for target in visible_targets:
    if target not in changed_colour:
        raise RuntimeError(
            f"Did not patch visible colour for plotting layer {target}"
        )

for target in hidden_targets:
    if target not in changed_alpha:
        raise RuntimeError(
            f"Did not hide redundant plotting layer {target}"
        )

text = "\n".join(lines) + "\n"


# ============================================================
# 5. Clean explanatory wording
# ============================================================

replacements = {
    "paired blue + magenta R ticks = B/DR.":
        "single orange R tick = B/DR.",

    "paired blue + magenta ticks":
        "single heterozygote-coloured tick",

    "paired blue + magenta R ticks":
        "single orange R tick",

    "single magenta L tick = DR/DR;":
        "single magenta L tick = DR/DR;",

    "High-BA sites are taller and fully saturated.":
        "High-BA sites are taller and fully saturated."
}

for old_text, new_text in replacements.items():
    text = text.replace(old_text, new_text)


# ============================================================
# 6. Replace the ENTIRE sidecar function.
#
# This is robust because we replace between two function markers,
# rather than trying to identify individual grobs.
#
# The key now shows exactly ONE supportive pattern:
#
#       L   DR/DR   magenta
#       R   B/DR    orange
#
# No B/B.
# No double dash.
# ============================================================

start_marker = "make_supportive_key <- function() {"
end_marker = "add_sidecar <- function("

start = text.find(start_marker)
end = text.find(end_marker, start)

if start == -1 or end == -1:
    raise RuntimeError(
        "Could not locate make_supportive_key()/add_sidecar() boundary."
    )

new_sidecar = r'''make_supportive_key <- function() {

    grobTree(

        textGrob(
            "Supportive focal SNP",
            x = unit(0.06, "npc"),
            y = unit(0.94, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 8.2,
                fontface = "bold",
                col = "#303030"
            )
        ),

        textGrob(
            "Only this ordered L/R genotype pattern is highlighted",
            x = unit(0.06, "npc"),
            y = unit(0.875, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.7,
                col = "#666666"
            )
        ),

        # ====================================================
        # LEFT
        # ====================================================

        textGrob(
            "L",
            x = unit(0.07, "npc"),
            y = unit(0.69, "npc"),
            hjust = 1,
            gp = gpar(
                fontsize = 8,
                fontface = "bold",
                col = "#333333"
            )
        ),

        segmentsGrob(
            x0 = unit(0.12, "npc"),
            x1 = unit(0.12, "npc"),
            y0 = unit(0.64, "npc"),
            y1 = unit(0.74, "npc"),
            gp = gpar(
                col = DR_colour,
                lwd = 3.5
            )
        ),

        grobTree(
            DR_cartoon,
            vp = viewport(
                x = 0.26,
                y = 0.69,
                width = 0.20,
                height = 0.16
            )
        ),

        textGrob(
            "DR/DR",
            x = unit(0.39, "npc"),
            y = unit(0.71, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 7.0,
                fontface = "bold",
                col = DR_colour
            )
        ),

        textGrob(
            "DR only",
            x = unit(0.39, "npc"),
            y = unit(0.65, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.8,
                col = "#555555"
            )
        ),

        # ====================================================
        # Same SNP / transition
        # ====================================================

        segmentsGrob(
            x0 = unit(0.12, "npc"),
            x1 = unit(0.12, "npc"),
            y0 = unit(0.58, "npc"),
            y1 = unit(0.49, "npc"),
            arrow = arrow(
                length = unit(0.055, "in"),
                type = "closed"
            ),
            gp = gpar(
                col = "#808080",
                lwd = 0.9
            )
        ),

        textGrob(
            "same SNP",
            x = unit(0.18, "npc"),
            y = unit(0.535, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.7,
                col = "#777777"
            )
        ),

        # ====================================================
        # RIGHT
        # ====================================================

        textGrob(
            "R",
            x = unit(0.07, "npc"),
            y = unit(0.39, "npc"),
            hjust = 1,
            gp = gpar(
                fontsize = 8,
                fontface = "bold",
                col = "#333333"
            )
        ),

        segmentsGrob(
            x0 = unit(0.12, "npc"),
            x1 = unit(0.12, "npc"),
            y0 = unit(0.34, "npc"),
            y1 = unit(0.44, "npc"),
            gp = gpar(
                col = HET_colour,
                lwd = 3.5
            )
        ),

        grobTree(
            band_cartoon,
            vp = viewport(
                x = 0.22,
                y = 0.39,
                width = 0.14,
                height = 0.11
            )
        ),

        textGrob(
            "+",
            x = unit(0.285, "npc"),
            y = unit(0.39, "npc"),
            gp = gpar(
                fontsize = 7.5,
                fontface = "bold",
                col = "#555555"
            )
        ),

        grobTree(
            DR_cartoon,
            vp = viewport(
                x = 0.35,
                y = 0.39,
                width = 0.14,
                height = 0.11
            )
        ),

        textGrob(
            "B/DR",
            x = unit(0.45, "npc"),
            y = unit(0.41, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 7.0,
                fontface = "bold",
                col = HET_colour
            )
        ),

        textGrob(
            "Band + DR",
            x = unit(0.45, "npc"),
            y = unit(0.35, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.8,
                col = "#555555"
            )
        ),

        # ====================================================
        # Colour grammar
        # ====================================================

        textGrob(
            "Supplement genoplot colours",
            x = unit(0.06, "npc"),
            y = unit(0.20, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 6.7,
                fontface = "bold",
                col = "#444444"
            )
        ),

        segmentsGrob(
            x0 = unit(0.10, "npc"),
            x1 = unit(0.10, "npc"),
            y0 = unit(0.11, "npc"),
            y1 = unit(0.16, "npc"),
            gp = gpar(
                col = DR_colour,
                lwd = 3.0
            )
        ),

        textGrob(
            "DR/DR",
            x = unit(0.16, "npc"),
            y = unit(0.135, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.9,
                col = "#555555"
            )
        ),

        segmentsGrob(
            x0 = unit(0.37, "npc"),
            x1 = unit(0.37, "npc"),
            y0 = unit(0.11, "npc"),
            y1 = unit(0.16, "npc"),
            gp = gpar(
                col = HET_colour,
                lwd = 3.0
            )
        ),

        textGrob(
            "B/DR",
            x = unit(0.43, "npc"),
            y = unit(0.135, "npc"),
            hjust = 0,
            gp = gpar(
                fontsize = 5.9,
                col = "#555555"
            )
        )
    )
}


'''

text = text[:start] + new_sidecar + text[end:]


# ============================================================
# 7. Write derivative script
# ============================================================

dst.write_text(text)

print(f"Created: {dst}")
print()
print("170 lower-panel grammar:")
print("  Left  DR/DR = #B61F72")
print("  Right B/DR  = #F28E2B")
print("  Right is ONE tick, not a double dash.")
print()
print("Biological definition unchanged:")
print("  supportive = ordered Left DR/DR -> Right B/DR")
