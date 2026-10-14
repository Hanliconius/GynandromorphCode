#!/usr/bin/env python3

import ast
import hashlib
from pathlib import Path


PROJECT = Path(
    "/scratch/martinlab/joe/GynandroCalling"
)

SOURCE = (
    PROJECT
    / "scripts"
    / "141_B_western_plus_heurippa_adjusted_LR.py"
)

TARGET = (
    PROJECT
    / "scripts"
    / "142_B_western_postman_vs_western_Amazon_DR_adjusted_LR.py"
)


if not SOURCE.is_file():
    raise RuntimeError(
        f"Missing working source: {SOURCE}"
    )


raw = SOURCE.read_text()

tree = ast.parse(
    raw,
    filename=str(
        SOURCE
    ),
)

lines = raw.splitlines(
    keepends=True
)


def target_name(
    node,
    name,
):

    return (
        isinstance(
            node,
            ast.Assign,
        )
        and any(
            isinstance(
                target,
                ast.Name,
            )
            and target.id == name
            for target in node.targets
        )
    )


def one(
    nodes,
    label,
):

    if len(
        nodes
    ) != 1:

        raise RuntimeError(
            f"Expected one {label}; found {len(nodes)}"
        )

    return nodes[
        0
    ]


def segment(
    node,
):

    return "".join(
        lines[
            node.lineno - 1:
            node.end_lineno
        ]
    ).rstrip(
        "\n"
    )


rules_node = one(
    [
        node
        for node in ast.walk(
            tree
        )
        if target_name(
            node,
            "RULES",
        )
    ],
    "RULES assignment",
)


metadata_node = one(
    [
        node
        for node in ast.walk(
            tree
        )
        if target_name(
            node,
            "metadata",
        )
        and "load_metadata" in ast.dump(
            node.value
        )
    ],
    "metadata load assignment",
)


kept_node = one(
    [
        node
        for node in ast.walk(
            tree
        )
        if target_name(
            node,
            "retained",
        )
        and "manifest" in ast.dump(
            node.value
        )
    ],
    "retained-cohort assignment",
)


message_node = one(
    [
        node
        for node in ast.walk(
            tree
        )
        if isinstance(
            node,
            ast.Expr,
        )
        and "No D/R filter" in (
            ast.get_source_segment(
                raw,
                node,
            )
            or ""
        )
    ],
    "old cohort-description print",
)


new_rules = '''SOURCES["Wallbank2016"] = "https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1002353"

RULES = {
    "melpomene cythera": (
        "Pacific_west_postman",
        ["ecuador"],
        "Hausdorf2025",
    ),

    "melpomene vulcanus": (
        "Pacific_west_or_Panama_postman",
        ["colombia", "panama"],
        "Hausdorf2025",
    ),

    "melpomene rosina": (
        "Central_America_postman",
        ["panama", "costa rica"],
        "Hausdorf2025",
    ),

    "melpomene aglaope": (
        "western_Amazon_Dennis_Ray",
        ["peru"],
        "Wallbank2016",
    ),

    "melpomene ecuadoriensis": (
        "western_Amazon_Dennis_Ray",
        ["ecuador"],
        "Wallbank2016",
    ),

    "melpomene malleti": (
        "western_Amazon_Dennis_Ray",
        ["colombia", "ecuador"],
        "Wallbank2016",
    ),

    "timareta florencia": (
        "western_Amazon_Dennis_Ray",
        ["colombia"],
        "Wallbank2016",
    ),

    "timareta timareta f. contigua": (
        "western_Amazon_Ray_only",
        ["ecuador"],
        "Wallbank2016",
    ),
}
'''


new_metadata = segment(
    metadata_node
) + '''

    corrections = []

    for accession, record in metadata.items():

        taxon = " ".join(
            str(
                record.get(
                    "taxon",
                    record.get(
                        "taxon_name",
                        "",
                    ),
                )
            )
            .lower()
            .replace(
                "\\n",
                " ",
            )
            .split()
        )

        if (
            "melpomene" in taxon
            and (
                "malleti" in taxon
                or "ecuadoriensis" in taxon
            )
        ):

            old_value = record.get(
                "B"
            )

            record[
                "B"
            ] = 0

            if old_value != 0:

                corrections.append(
                    (
                        accession,
                        "B",
                        old_value,
                        0,
                    )
                )

        if accession == "ERR3652402":

            for trait, value in (
                ("B", 0),
                ("D", 1),
                ("R", 1),
            ):

                old_value = record.get(
                    trait
                )

                record[
                    trait
                ] = value

                if old_value != value:

                    corrections.append(
                        (
                            accession,
                            trait,
                            old_value,
                            value,
                        )
                    )

    print(
        "LOCKED PHENOTYPE CORRECTIONS"
    )

    if corrections:

        print(
            "\\n".join(
                "\\t".join(
                    map(
                        str,
                        correction,
                    )
                )
                for correction in corrections
            )
        )

    else:

        print(
            "No stored values changed."
        )
'''


new_kept = '''    for trait in ("D", "R"):

        manifest[trait] = pd.to_numeric(
            manifest["sample"].map(
                lambda sample:
                    metadata.get(
                        sample,
                        {},
                    ).get(
                        trait
                    )
            ),
            errors="coerce",
        )

    # ERR3653308 is a contigua specimen coded
    # B=0, D=0, R=0 and therefore is not part of the
    # predefined Dennis/Ray-positive B=0 comparison.
    excluded_accessions = {
        "ERR3653308",
    }

    excluded_mask = manifest["sample"].isin(
        excluded_accessions
    )

    recovered_exclusions = sorted(
        manifest.loc[
            excluded_mask,
            "sample",
        ].unique()
    )

    if recovered_exclusions != sorted(
        excluded_accessions
    ):
        raise RuntimeError(
            "Failed to recover locked exclusion: "
            "ERR3653308"
        )

    manifest.loc[
        excluded_mask,
        "keep",
    ] = False

    manifest.loc[
        excluded_mask,
        "keep_bool",
    ] = False

    manifest.loc[
        excluded_mask,
        "reason",
    ] = "excluded_D0R0_contigua"

    print(
        "LOCKED COMPARISON-GROUP EXCLUSIONS"
    )

    print(
        ",".join(
            recovered_exclusions
        )
    )

''' + segment(
    kept_node
) + '''

    retained = retained.copy()

    for trait in ("B", "D", "R"):

        retained[trait] = pd.to_numeric(
            retained[trait],
            errors="coerce",
        )

    taxa = (
        retained["taxon_name"]
        .astype(str)
        .str.lower()
    )

    forbidden = taxa.str.contains(
        r"cydno|pachinus",
        regex=True,
        na=False,
    )

    if forbidden.any():

        raise RuntimeError(
            "Cydno or pachinus survived:\\n"
            + retained.loc[
                forbidden,
                [
                    "cohort",
                    "sample",
                    "taxon_name",
                    "B",
                    "D",
                    "R",
                ],
            ].to_string(
                index=False
            )
        )

    b0 = retained.loc[
        retained["B"].eq(0)
    ].copy()

    b1 = retained.loc[
        retained["B"].eq(1)
    ].copy()

    if b0.empty or b1.empty:

        raise RuntimeError(
            "Both B classes must be represented."
        )

    bad_b0 = ~(
        b0["D"].eq(1)
        |
        b0["R"].eq(1)
    )

    if bad_b0.any():

        raise RuntimeError(
            "B=0 but neither D nor R positive:\\n"
            + b0.loc[
                bad_b0,
                [
                    "cohort",
                    "sample",
                    "taxon_name",
                    "B",
                    "D",
                    "R",
                ],
            ].to_string(
                index=False
            )
        )

    bad_b1 = (
        b1["D"].eq(1)
        |
        b1["R"].eq(1)
    )

    if bad_b1.any():

        raise RuntimeError(
            "B=1 but D or R positive:\\n"
            + b1.loc[
                bad_b1,
                [
                    "cohort",
                    "sample",
                    "taxon_name",
                    "B",
                    "D",
                    "R",
                ],
            ].to_string(
                index=False
            )
        )

    heurippa = taxa.str.contains(
        r"heliconius\\s+heurippa",
        regex=True,
        na=False,
    )

    if not heurippa.any():

        raise RuntimeError(
            "Heurippa was not retained."
        )

    retained[
        "comparative_role"
    ] = "western_postman_or_heurippa_B1"

    retained.loc[
        retained["B"].eq(0),
        "comparative_role",
    ] = "western_Amazon_Dennis_or_Ray_B0"

    audit = (
        retained.groupby(
            [
                "cohort",
                "comparative_role",
                "B",
                "D",
                "R",
            ],
            dropna=False,
        )
        .size()
        .rename(
            "samples"
        )
        .reset_index()
    )

    audit_dir = (
        PROJECT
        / "analysis"
        / "B_western_postman_vs_western_Amazon_DR_reanalysis"
        / "checkpoints"
    )

    audit_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    audit.to_csv(
        audit_dir
        / "B_western_retained_B_D_R_counts.tsv",
        sep="\\t",
        index=False,
    )

    retained.to_csv(
        audit_dir
        / "B_western_retained_sample_manifest.tsv",
        sep="\\t",
        index=False,
    )

    print(
        "WESTERN POSTMAN VERSUS "
        "WESTERN-AMAZON D/R COHORT"
    )

    print(
        audit.to_string(
            index=False
        )
    )
'''


new_message = '''    print(
        "B=0: literature-defined western-Amazon "
        "Dennis/Ray taxa.\\n"
        "B=1: western postman taxa plus heurippa.\\n"
        "Cydno and pachinus excluded.\\n"
        "Comparative scoring uses no focal calls; "
        "the focal display uses reciprocal-mixture "
        "posterior >= 0.95 calls only.",
        flush=True,
    )
'''


replacements = [
    (
        rules_node.lineno,
        rules_node.end_lineno,
        new_rules,
    ),

    (
        metadata_node.lineno,
        metadata_node.end_lineno,
        new_metadata,
    ),

    (
        kept_node.lineno,
        kept_node.end_lineno,
        new_kept,
    ),

    (
        message_node.lineno,
        message_node.end_lineno,
        new_message,
    ),
]


for start, end, replacement in sorted(
    replacements,
    key=lambda item: item[0],
    reverse=True,
):

    lines[
        start - 1:
        end
    ] = [
        replacement.rstrip(
            "\n"
        )
        + "\n"
    ]


text = "".join(
    lines
)


literal_changes = {
    "B_western_plus_heurippa_reanalysis":
        "B_western_postman_vs_western_Amazon_DR_reanalysis",

    "141_plot_B_western_plus_heurippa_adjusted_LR.R":
        "142_plot_B_western_postman_vs_western_Amazon_DR_adjusted_LR.R",

    "Band association in western samples plus heurippa":
        "Band association: western postman versus western-Amazon Dennis/Ray",

    'version="141.1"':
        'version="142.1"',

    '"version": "141.1"':
        '"version": "142.1"',
}


for old, new in literal_changes.items():

    text = text.replace(
        old,
        new,
    )


compile(
    text,
    str(
        TARGET
    ),
    "exec",
)

TARGET.write_text(
    text
)


print(
    "ANALYSIS 142 CREATED STRUCTURALLY"
)

print(
    "source_SHA256\t"
    + hashlib.sha256(
        raw.encode()
    ).hexdigest()
)

print(
    "target_SHA256\t"
    + hashlib.sha256(
        text.encode()
    ).hexdigest()
)

print(
    f"RULES_lines\t{rules_node.lineno}-{rules_node.end_lineno}"
)

print(
    f"metadata_lines\t"
    f"{metadata_node.lineno}-{metadata_node.end_lineno}"
)

print(
    f"kept_lines\t{kept_node.lineno}-{kept_node.end_lineno}"
)
