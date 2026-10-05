#!/usr/bin/env bash

set -u
set -o pipefail

main() {
    local project_root
    local checkpoint_dir
    local expectation_table
    local source_audit
    local current_focal
    local current_summary
    local candidate_list
    local tmp_header
    local f
    local header
    local rows
    local band_hits
    local module_hint

    project_root="/scratch/martinlab/joe/GynandroCalling"
    checkpoint_dir="${project_root}/analysis/original_502/checkpoints"

    expectation_table="${checkpoint_dir}/focal_module_expectations.tsv"
    source_audit="${checkpoint_dir}/comparative_site_score_source_candidates.tsv"

    current_focal="${checkpoint_dir}/focal_calls_at_independent_comparative_sites.tsv"
    current_summary="${checkpoint_dir}/focal_calls_at_independent_comparative_sites_summary.tsv"

    candidate_list="${checkpoint_dir}/.90a_candidate_tables.txt"
    tmp_header="${checkpoint_dir}/.90a_header.txt"

    mkdir -p \
        "${project_root}/scripts" \
        "${project_root}/logs" \
        "${checkpoint_dir}"

    # ------------------------------------------------------------------
    # CANONICAL, A-PRIORI BIOLOGICAL EXPECTATION TABLE
    #
    # P = biological presence allele / presence-associated haplotype
    # A = biological absence allele / absence-associated haplotype
    #
    # Biological presence dosage:
    #   AA = 0
    #   PA = 0.5
    #   PP = 1
    #
    # These expectations are fixed independently of focal genotype data.
    # ------------------------------------------------------------------

    cat > "${expectation_table}" <<'TABLE'
moduleblockleft_phenotyperight_phenotypeleft_expected_genotyperight_expected_genotypeleft_expected_presence_dosageright_expected_presence_dosagereference_biological_genotypereference_presence_dosagereference_trait_state
BandBBand absentBand presentAAPA00.5PP1presence
Dennis1DDennis presentDennis presentPPPA10.5AA0absence
Dennis2DDennis presentDennis presentPPPA10.5AA0absence
Ray1RRay presentRay presentPPPA10.5AA0absence
Ray2RRay presentRay presentPPPA10.5AA0absence
YbYbYb shadow-bandYb presentPAPP0.51PP1presence
TABLE

    printf '\nCANONICAL FOCAL EXPECTATION TABLE\n'
    column -t -s $'\t' "${expectation_table}"

    # ------------------------------------------------------------------
    # FIND COMPARATIVE SITE-SCORE TABLES.
    #
    # We deliberately inspect existing tables rather than regenerating
    # VCF-derived data.
    #
    # Candidate definition:
    #   TSV under analysis/
    #   header mentions a position field
    #   AND header mentions balanced accuracy / BA
    #
    # No file is modified.
    # ------------------------------------------------------------------

    find "${project_root}/analysis" \
        -type f \
        -name '*.tsv' \
        ! -path '*/focal_calls_at_independent_comparative_sites.tsv' \
        ! -path '*/focal_calls_at_independent_comparative_sites_summary.tsv' \
        ! -path '*/focal_module_expectations.tsv' \
        -print \
        | sort \
        > "${candidate_list}"

    printf '%s\t%s\t%s\t%s\t%s\n' \
        "file" \
        "rows" \
        "module_hint" \
        "known_band_marker_hits" \
        "header" \
        > "${source_audit}"

    while IFS= read -r f
    do
        if [[ ! -s "${f}" ]]
        then
            continue
        fi

        head -n 1 "${f}" > "${tmp_header}"
        header="$(cat "${tmp_header}")"

        if ! printf '%s\n' "${header}" |
            grep -Eqi '(^|)(pos|position|POS|site_pos|variant_position)(|$)|position'
        then
            continue
        fi

        if ! printf '%s\n' "${header}" |
            grep -Eqi 'balanced.?accuracy|(^|)BA(|$)|(^|)ba(|$)|min_BA|max_BA'
        then
            continue
        fi

        rows="$(
            awk 'END { print (NR > 0 ? NR - 1 : 0) }' "${f}"
        )"

        module_hint="unknown"

        case "${f}" in
            *band*|*Band*)
                module_hint="Band"
                ;;
            *dennis1*|*Dennis1*|*dennis_1*)
                module_hint="Dennis1"
                ;;
            *dennis2*|*Dennis2*|*dennis_2*)
                module_hint="Dennis2"
                ;;
            *ray1*|*Ray1*|*ray_1*)
                module_hint="Ray1"
                ;;
            *ray2*|*Ray2*|*ray_2*)
                module_hint="Ray2"
                ;;
            *Yb*|*yb*)
                module_hint="Yb"
                ;;
        esac

        # Exact known Band markers from the earlier comparative analysis.
        # This is diagnostic only: it does not select sites.
        band_hits="$(
            grep -E -c \
                '(^|)(1112841|1112866|1113027)(|$)' \
                "${f}" \
                2>/dev/null ||
            true
        )"

        printf '%s\t%s\t%s\t%s\t%s\n' \
            "${f}" \
            "${rows}" \
            "${module_hint}" \
            "${band_hits}" \
            "${header}" \
            >> "${source_audit}"

    done < "${candidate_list}"

    rm -f "${candidate_list}" "${tmp_header}"

    printf '\nCOMPARATIVE SCORE TABLE CANDIDATES\n'
    if [[ -s "${source_audit}" ]]
    then
        column -t -s $'\t' "${source_audit}"
    else
        printf 'No candidate tables detected.\n'
    fi

    # ------------------------------------------------------------------
    # CURRENT FROZEN TABLE: record its exact schema so the production
    # recoder can preserve all useful focal GT/AD/DP/GQ fields.
    # ------------------------------------------------------------------

    printf '\nCURRENT FROZEN FOCAL TABLE\n'

    if [[ -s "${current_focal}" ]]
    then
        printf 'Path: %s\n\n' "${current_focal}"

        printf 'Header:\n'
        head -n 1 "${current_focal}" |
            tr '\t' '\n' |
            nl -ba

        printf '\nFirst five data rows:\n'
        {
            head -n 1 "${current_focal}"
            sed -n '2,6p' "${current_focal}"
        } |
            column -t -s $'\t'
    else
        printf 'MISSING: %s\n' "${current_focal}"
    fi

    printf '\nCURRENT FROZEN SUMMARY\n'

    if [[ -s "${current_summary}" ]]
    then
        column -t -s $'\t' "${current_summary}"
    else
        printf 'MISSING: %s\n' "${current_summary}"
    fi

    printf '\nBAND SOURCE CHECK\n'

    awk -F '\t' '
        NR == 1 {
            next
        }
        $4 + 0 > 0 {
            print
            found = 1
        }
        END {
            if (!found) {
                print "NO_TABLE_CONTAINING_KNOWN_BAND_MARKERS_FOUND"
            }
        }
    ' "${source_audit}" |
        column -t -s $'\t'

    printf '\nOutputs written:\n'
    printf '  %s\n' "${expectation_table}"
    printf '  %s\n' "${source_audit}"

    printf '\nNo focal calls or comparative site sets were modified.\n'
}

main "$@"
