#!/usr/bin/env bash

project_dir="/scratch/martinlab/joe/GynandroCalling"
output="${project_dir}/chat_test_success_handoff_audit.txt"
since_date="2026-08-31 00:00:00"

cd "${project_dir}" || {
    printf '%s\n' "ERROR: could not enter ${project_dir}"
}

{
    printf '%s\n\n' "CHAT TEST SUCCESS HANDOFF AUDIT"

    printf '%s\n' "PROJECT"
    printf 'Directory: %s\n' "${project_dir}"
    printf 'Audit time: %s\n\n' "$(date -Is)"

    printf '%s\n' "SCRIPTS 90-99"
    find scripts \
        -maxdepth 1 \
        -type f \
        -newermt "${since_date}" \
        \( -name '9[0-9]*' -o -name '90*' \) \
        -printf '%TY-%Tm-%Td %TH:%TM:%TS\t%s bytes\t%p\n' \
        | sort
    printf '\n'

    printf '%s\n' "LOGS 90-99"
    find logs \
        -maxdepth 1 \
        -type f \
        -newermt "${since_date}" \
        \( -name '9[0-9]*' -o -name '90*' \) \
        -printf '%TY-%Tm-%Td %TH:%TM:%TS\t%s bytes\t%p\n' \
        | sort
    printf '\n'

    printf '%s\n' "RECENT CHECKPOINT FILES"
    find analysis \
        -type f \
        -newermt "${since_date}" \
        \( -name '*.tsv' -o -name '*.txt' -o -name '*.csv' \) \
        -printf '%TY-%Tm-%Td %TH:%TM:%TS\t%s bytes\t%p\n' \
        | sort
    printf '\n'

    printf '%s\n' "RECENT FIGURES"
    find analysis \
        -type f \
        -newermt "${since_date}" \
        \( -name '*.pdf' -o -name '*.png' -o -name '*.svg' \) \
        -printf '%TY-%Tm-%Td %TH:%TM:%TS\t%s bytes\t%p\n' \
        | sort
    printf '\n'

    printf '%s\n' "CONTENTS OF RECENT 90-99 LOGS"

    while IFS= read -r logfile
    do
        line_count="$(wc -l < "${logfile}")"

        printf '\n%s\n' "============================================================"
        printf 'LOG: %s\n' "${logfile}"
        printf 'LINES: %s\n' "${line_count}"
        printf '%s\n' "============================================================"

        if [[ "${line_count}" -le 300 ]]
        then
            cat "${logfile}"
        else
            printf '%s\n' "[FIRST 40 LINES]"
            head -n 40 "${logfile}"
            printf '%s\n' "[LAST 160 LINES]"
            tail -n 160 "${logfile}"
        fi
    done < <(
        find logs \
            -maxdepth 1 \
            -type f \
            -newermt "${since_date}" \
            \( -name '9[0-9]*.log' -o -name '90*.log' \) \
            | sort
    )

    printf '\n%s\n' "HEADERS AND FIRST ROWS OF RECENT CHECKPOINT TABLES"

    while IFS= read -r table
    do
        printf '\n%s\n' "------------------------------------------------------------"
        printf 'TABLE: %s\n' "${table}"
        head -n 3 "${table}"
    done < <(
        find analysis \
            -type f \
            -newermt "${since_date}" \
            \( -name '*.tsv' -o -name '*.csv' \) \
            | sort
    )
} > "${output}"

printf 'Audit written to:\n%s\n' "${output}"
printf 'Size: '
du -h "${output}" | cut -f1
