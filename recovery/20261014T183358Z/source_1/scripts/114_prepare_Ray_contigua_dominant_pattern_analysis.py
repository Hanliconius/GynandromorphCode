#!/usr/bin/env python3

from pathlib import Path
import py_compile


project = Path(
    "/scratch/martinlab/joe/GynandroCalling"
)

source = (
    project
    / "scripts"
    / "113_replicate_Ray_with_current_contigua.py"
)

target = (
    project
    / "scripts"
    / "114_Ray_contigua_dominant_pattern_analysis.py"
)


if not source.is_file():
    raise RuntimeError(
        f"Missing source scanner: {source}"
    )


text = source.read_text()


filename_replacements = {
    "Wallbank2016_Ray_contigua_replication_samples.tsv":
        "Ray_contigua_dominant_pattern_samples.tsv",

    "Wallbank2016_Ray_contigua_replication_sites.tsv":
        "Ray_contigua_dominant_pattern_sites.tsv",

    "Wallbank2016_Ray_contigua_replication_5kb_windows.tsv":
        "Ray_contigua_dominant_pattern_5kb_windows.tsv",

    "Wallbank2016_Ray_contigua_replication_summary.tsv":
        "Ray_contigua_dominant_pattern_summary.tsv",
}


for old, new in filename_replacements.items():

    if old not in text:
        raise RuntimeError(
            f"Expected filename was not found: {old}"
        )

    text = text.replace(
        old,
        new,
    )


column_replacements = {
    "fixed_among_called_min2":
        "dominant_pattern_among_called_min2",

    "fixed_80pct_callable":
        "dominant_pattern_80pct_callable",

    "fixed_complete":
        "dominant_pattern_complete",

    "group_A_fixed_allele":
        "no_R_homozygous_allele",

    "group_B_fixed_allele":
        "R_associated_allele",

    "total_fixed_among_called_min2":
        "total_dominant_pattern_among_called_min2",

    "total_fixed_80pct_callable":
        "total_dominant_pattern_80pct_callable",

    "total_fixed_complete":
        "total_dominant_pattern_complete",

    "maximum_5kb_fixed_among_called_min2":
        "maximum_5kb_dominant_pattern_among_called_min2",

    "maximum_5kb_fixed_80pct_callable":
        "maximum_5kb_dominant_pattern_80pct_callable",

    "maximum_5kb_fixed_complete":
        "maximum_5kb_dominant_pattern_complete",
}


for old, new in column_replacements.items():
    text = text.replace(
        old,
        new,
    )


function_start = text.find(
    "def evaluate_fixed_difference("
)

function_end = text.find(
    "\n\nif not VCF.is_file():",
    function_start,
)


if (
    function_start < 0
    or function_end < 0
):
    raise RuntimeError(
        "Could not locate the fixed-difference "
        "evaluation function."
    )


new_function = '''def evaluate_dominant_ray_pattern(
    calls_no_R,
    calls_R,
):
    """
    Dominant Ray model.

    The no-Ray group must be homozygous for one allele.

    Every callable Ray-present contigua may be either:
      - heterozygous, or
      - homozygous for the opposite, Ray-associated allele.

    Missing/failed-QC calls do not create discordance, but reduce
    callability. Both possible allele orientations are evaluated.
    """

    total_no_R = len(
        calls_no_R
    )

    total_R = len(
        calls_R
    )

    callable_no_R = [
        call
        for call in calls_no_R
        if call is not None
    ]

    callable_R = [
        call
        for call in calls_R
        if call is not None
    ]

    call_rate_no_R = (
        len(callable_no_R)
        / total_no_R
        if total_no_R
        else 0.0
    )

    call_rate_R = (
        len(callable_R)
        / total_R
        if total_R
        else 0.0
    )

    best = None

    for no_R_allele in (
        0,
        1,
    ):

        associated_allele = (
            1 - no_R_allele
        )

        no_R_correct = [
            call == no_R_allele
            for call in callable_no_R
        ]

        R_correct = [
            (
                call == "heterozygous"
                or call == associated_allele
            )
            for call in callable_R
        ]

        no_R_discordances = sum(
            not value
            for value in no_R_correct
        )

        R_discordances = sum(
            not value
            for value in R_correct
        )

        total_discordances = (
            no_R_discordances
            + R_discordances
        )

        candidate = {
            "no_R_homozygous_allele":
                no_R_allele,
            "R_associated_allele":
                associated_allele,
            "no_R_discordances":
                no_R_discordances,
            "R_discordances":
                R_discordances,
            "total_callable_discordances":
                total_discordances,
        }

        if (
            best is None
            or candidate[
                "total_callable_discordances"
            ]
            < best[
                "total_callable_discordances"
            ]
        ):
            best = candidate

    enough_calls = (
        len(callable_no_R) >= 2
        and len(callable_R) >= 2
    )

    perfect_among_called = (
        enough_calls
        and best[
            "total_callable_discordances"
        ] == 0
    )

    perfect_80pct_callable = (
        perfect_among_called
        and call_rate_no_R >= 0.80
        and call_rate_R >= 0.80
    )

    perfect_complete = (
        perfect_among_called
        and len(callable_no_R)
            == total_no_R
        and len(callable_R)
            == total_R
    )

    return {
        "group_A_n":
            total_no_R,
        "group_B_n":
            total_R,
        "group_A_called":
            len(callable_no_R),
        "group_B_called":
            len(callable_R),
        "group_A_call_rate":
            call_rate_no_R,
        "group_B_call_rate":
            call_rate_R,
        "no_R_homozygous_allele":
            best[
                "no_R_homozygous_allele"
            ],
        "R_associated_allele":
            best[
                "R_associated_allele"
            ],
        "no_R_discordances":
            best[
                "no_R_discordances"
            ],
        "R_discordances":
            best[
                "R_discordances"
            ],
        "total_callable_discordances":
            best[
                "total_callable_discordances"
            ],
        "dominant_pattern_among_called_min2":
            perfect_among_called,
        "dominant_pattern_80pct_callable":
            perfect_80pct_callable,
        "dominant_pattern_complete":
            perfect_complete,
    }
'''


text = (
    text[
        :function_start
    ]
    + new_function
    + text[
        function_end:
    ]
)


text = text.replace(
    "evaluate_fixed_difference(",
    "evaluate_dominant_ray_pattern(",
)


text = text.replace(
    "WALLBANK ET AL. 2016 FIXED-SNP REPLICATION",
    "RAY CONTIGUA DOMINANT-PATTERN ANALYSIS",
)

text = text.replace(
    '"Original analysis: 5-kb windows, 1-kb step."',
    (
        '"Dominant Ray model: no-Ray homozygote versus "\n'
        '    "heterozygous or associated-homozygous contigua."'
    ),
)

text = text.replace(
    '"Paper-QC calls require DP 4-300 and GQ >=30."',
    (
        '"Windows are 5 kb with a 1-kb step; paper-QC "\n'
        '    "calls require DP 4-300 and GQ >=30."'
    ),
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
print()
print(
    "Ray-present contigua may now be heterozygous "
    "or associated homozygous."
)
print(
    "No-Ray samples must remain homozygous for "
    "the opposite allele."
)
