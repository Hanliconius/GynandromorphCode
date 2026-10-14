#!/usr/bin/env python3

from pathlib import Path
import py_compile


project = Path(
    "/scratch/martinlab/joe/GynandroCalling"
)

source = (
    project
    / "scripts"
    / "112_replicate_wallbank2016_fixed_SNP_analysis.py"
)

target = (
    project
    / "scripts"
    / "113_replicate_Ray_with_current_contigua.py"
)


if not source.is_file():
    raise RuntimeError(
        f"Missing source scanner: {source}"
    )


text = source.read_text()


output_replacements = {
    "Wallbank2016_fixed_SNP_replication_samples.tsv":
        "Wallbank2016_Ray_contigua_replication_samples.tsv",

    "Wallbank2016_fixed_SNP_replication_sites.tsv":
        "Wallbank2016_Ray_contigua_replication_sites.tsv",

    "Wallbank2016_fixed_SNP_replication_5kb_windows.tsv":
        "Wallbank2016_Ray_contigua_replication_5kb_windows.tsv",

    "Wallbank2016_fixed_SNP_replication_summary.tsv":
        "Wallbank2016_Ray_contigua_replication_summary.tsv",
}


for old, new in output_replacements.items():

    if old not in text:
        raise RuntimeError(
            f"Expected output name was not found: {old}"
        )

    text = text.replace(
        old,
        new,
    )


analysis_start = text.find(
    "analyses = {"
)

analysis_end_anchor = (
    "\n\n# Read only the VCF header"
)

analysis_end = text.find(
    analysis_end_anchor,
    analysis_start,
)


if (
    analysis_start < 0
    or analysis_end < 0
):
    raise RuntimeError(
        "Could not locate the analysis dictionary "
        "in script 112."
    )


new_analysis = '''# ============================================================
# Hybrid Wallbank-ray replication.
#
# Presence group:
#   Current H. timareta f. contigua specimens with the
#   diagnostic D=0, R=1 phenotype.
#
# Absence group:
#   Original Wallbank et al. ray-absent biological specimens
#   that are represented in the current 502-sample VCF.
#
# The D=0, R=0 contigua record is excluded automatically.
# ============================================================

ray_present_current_contigua = sorted(
    set(
        metadata.loc[
            metadata[
                "taxon_clean"
            ].map(
                lambda value:
                taxon_matches(
                    value,
                    (
                        "timareta f. contigua",
                    ),
                )
            )
            & (
                metadata[
                    "D_value"
                ]
                == 0
            )
            & (
                metadata[
                    "R_value"
                ]
                == 1
            ),
            "Accession",
        ]
    )
    - {""}
)


ray_absent_original_overlap = [
    "ERR3653537", # CAM008533, H. t. timareta
    "ERR3653539", # CAM009178, H. t. timareta
    "ERR1143606", # CAM016293, H. m. plesseni
    "ERR1143605", # CAM009156, H. m. plesseni
    "ERR3653533", # CAM013819, H. m. meriana
]


analyses = {
    "Ray_current_contigua_vs_original_absent": {
        "trait": "Ray",
        "group_A_label": "original_no_R",
        "group_B_label": "current_contigua_D0R1",
        "group_A": ray_absent_original_overlap,
        "group_B": ray_present_current_contigua,
    },
}
'''


text = (
    text[
        :analysis_start
    ]
    + new_analysis
    + text[
        analysis_end:
    ]
)


target.write_text(
    text
)

py_compile.compile(
    str(target),
    doraise=True,
)

print(
    f"Created and syntax-checked: {target}"
)
