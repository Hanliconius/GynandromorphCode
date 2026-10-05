#!/usr/bin/env bash

set -u
set -o pipefail

simulation_job=138265731

old_stage202=138265838
old_stage203=138265839
old_stage204=138265840
old_stage205=138265841

command_for_job() {
    scontrol show job "$1" -o 2>/dev/null \
        | sed -n 's/.* Command=\([^ ]*\).*/\1/p'
}

stage202_script=$(command_for_job "${old_stage202}")
stage203_script=$(command_for_job "${old_stage203}")
stage204_script=$(command_for_job "${old_stage204}")
stage205_script=$(command_for_job "${old_stage205}")

preflight_status=0

printf '%s\n' 'SCRIPT PREFLIGHT'

for specification in \
    "202:${stage202_script}" \
    "203:${stage203_script}" \
    "204:${stage204_script}" \
    "205:${stage205_script}"
do
    stage=${specification%%:*}
    script=${specification#*:}

    if [[ -z "${script}" ]]
    then
        printf 'FAIL\tstage=%s\tcould not recover script path\n' "${stage}"
        preflight_status=1
    elif [[ ! -s "${script}" ]]
    then
        printf 'FAIL\tstage=%s\tmissing_or_empty=%s\n' \
            "${stage}" "${script}"
        preflight_status=1
    elif bash -n "${script}"
    then
        printf 'PASS\tstage=%s\t%s\n' "${stage}" "${script}"
    else
        printf 'FAIL\tstage=%s\tsyntax=%s\n' "${stage}" "${script}"
        preflight_status=1
    fi
done

printf '\n%s\n' 'STAGE 201 FAILED-TASK CHECK'

failed_tasks=$(
    sacct -X -j "${simulation_job}" \
        --format=JobIDRaw,State \
        --parsable2 \
        | awk -F '|' '
            $1 ~ /_[0-9]+$/ &&
            $2 !~ /^(COMPLETED|RUNNING|PENDING)$/ {
                sub(/^.*_/, "", $1)
                print $1
            }
        '
)

if [[ -n "${failed_tasks}" ]]
then
    printf 'FAIL\tfailed_tasks=%s\n' \
        "$(printf '%s\n' "${failed_tasks}" | paste -sd, -)"
    preflight_status=1
else
    printf 'PASS\tno failed Stage 201 tasks detected\n'
fi

submission_status=127
new_stage202=""
new_stage203=""
new_stage204=""
new_stage205=""

if [[ "${preflight_status}" -eq 0 ]]
then
    submission=$(
        sbatch \
            --parsable \
            --dependency="afterok:${simulation_job}" \
            "${stage202_script}"
    )
    submission_status=$?

    if [[ "${submission_status}" -eq 0 ]]
    then
        new_stage202=${submission%%;*}

        submission=$(
            sbatch \
                --parsable \
                --dependency="afterok:${new_stage202}" \
                "${stage203_script}"
        )
        submission_status=$?
    fi

    if [[ "${submission_status}" -eq 0 ]]
    then
        new_stage203=${submission%%;*}

        submission=$(
            sbatch \
                --parsable \
                --dependency="afterok:${new_stage203}" \
                "${stage204_script}"
        )
        submission_status=$?
    fi

    if [[ "${submission_status}" -eq 0 ]]
    then
        new_stage204=${submission%%;*}

        submission=$(
            sbatch \
                --parsable \
                --dependency="afterok:${new_stage204}" \
                "${stage205_script}"
        )
        submission_status=$?
    fi

    if [[ "${submission_status}" -eq 0 ]]
    then
        new_stage205=${submission%%;*}
    fi
fi

printf '\n%s\n' 'FINAL STATUS'
printf 'preflight_status=%s\n' "${preflight_status}"
printf 'submission_status=%s\n' "${submission_status}"
printf 'new_stage202=%s\n' "${new_stage202}"
printf 'new_stage203=%s\n' "${new_stage203}"
printf 'new_stage204=%s\n' "${new_stage204}"
printf 'new_stage205=%s\n' "${new_stage205}"

if [[ \
    "${submission_status}" -eq 0 &&
    -n "${new_stage205}" \
]]
then
    scancel \
        "${old_stage203}" \
        "${old_stage204}" \
        "${old_stage205}" \
        2>/dev/null

    printf '\n%s\n' 'REPLACEMENT CHAIN'
    printf 'Stage 201 simulations\t%s\n' "${simulation_job}"
    printf 'Stage 202 collection\t%s\tafterok:%s\n' \
        "${new_stage202}" "${simulation_job}"
    printf 'Stage 203 evidence\t%s\tafterok:%s\n' \
        "${new_stage203}" "${new_stage202}"
    printf 'Stage 204 narrative\t%s\tafterok:%s\n' \
        "${new_stage204}" "${new_stage203}"
    printf 'Stage 205 validation\t%s\tafterok:%s\n' \
        "${new_stage205}" "${new_stage204}"

    printf '\n%s\n' 'CURRENT QUEUE'

    squeue -j \
        "${simulation_job},${new_stage202},${new_stage203},${new_stage204},${new_stage205}"
else
    printf '\nReplacement chain was not completed; old downstream jobs were preserved.\n'
fi
