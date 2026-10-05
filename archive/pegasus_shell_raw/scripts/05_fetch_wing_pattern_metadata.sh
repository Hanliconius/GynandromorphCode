#!/usr/bin/env bash
set -u
set -o pipefail

main() {
    local project_root
    local metadata_dir
    local spreadsheet_id
    local export_url
    local output
    local temporary

    project_root="/scratch/martinlab/joe/GynandroCalling"
    metadata_dir="${project_root}/metadata"

    spreadsheet_id="1wNaOR-9KqBlaSi5PgAEQuYKScU3AP8ZqDkabtoucUss"
    export_url="https://docs.google.com/spreadsheets/d/${spreadsheet_id}/export?format=xlsx"

    output="${metadata_dir}/Wing_Pattern_MetaData.xlsx"
    temporary="${output}.part"

    mkdir -p "${metadata_dir}"

    validate_workbook() {
        local workbook="$1"

        python - "${workbook}" <<'PY'
import sys
from openpyxl import load_workbook

path = sys.argv[1]
workbook = load_workbook(path, read_only=True, data_only=True)

if not workbook.sheetnames:
    raise ValueError("Workbook contains no worksheets")

print(
    "Validated workbook: "
    + ",".join(workbook.sheetnames)
)

workbook.close()
PY
    }

    if [[ -s "${output}" ]]
    then
        if validate_workbook "${output}" >/dev/null 2>&1
        then
            printf 'Reusing valid workbook: %s\n' "${output}"
            return 0
        else
            printf 'ERROR: existing workbook is invalid: %s\n' \
                "${output}" >&2
            return 1
        fi
    fi

    rm -f "${temporary}"

    if curl \
        --location \
        --fail \
        --silent \
        --show-error \
        --retry 3 \
        --output "${temporary}" \
        "${export_url}"
    then
        if [[ -s "${temporary}" ]] &&
            validate_workbook "${temporary}" >/dev/null 2>&1
        then
            mv "${temporary}" "${output}"

            printf 'Downloaded and validated native Google Sheet export:\n'
            printf '%s\n' "${output}"
            printf 'Bytes: %s\n' "$(stat -c '%s' "${output}")"
            return 0
        else
            printf 'ERROR: export was not a valid XLSX workbook.\n' >&2
            rm -f "${temporary}"
            return 1
        fi
    else
        printf 'ERROR: native Google Sheets export failed.\n' >&2
        rm -f "${temporary}"
        return 1
    fi
}

main "$@"
