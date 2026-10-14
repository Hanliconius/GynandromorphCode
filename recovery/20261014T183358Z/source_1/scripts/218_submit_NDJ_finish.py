#!/usr/bin/env python3
"""Submit stages after the NDJ alignment array; record each job immediately."""
import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--alignment-job', required=True)
    p.add_argument('--project', type=Path, default=Path.cwd())
    args = p.parse_args()
    if not args.alignment_job.isdigit():
        p.error('alignment-job must be a numeric parent array job ID')
    project = args.project.resolve()
    python = project / 'pairwise-af/.venv-af/bin/python'
    driver = project / 'scripts/218_NDJ_controls.py'
    for path in (python, driver):
        if not path.is_file():
            raise RuntimeError(f'Missing {path}')
    subprocess.run([str(python), str(driver), 'preflight', '--project', str(project)], check=True)
    logs = project / 'logs'
    logs.mkdir(exist_ok=True)
    root = project / 'analysis/NDJ_control_comparison_v1'
    root.mkdir(parents=True, exist_ok=True)
    state_file = root / 'finish_chain.submissions.jsonl'
    # Refuse duplicate chains. A partial chain is recorded for targeted recovery.
    with state_file.open('x') as state:
        dependency = args.alignment_job
        jobs = [('alignment', dependency)]
        stages = [('prepare', 'nano', '00:30:00', 1, '4G', 'NDJprep'),
                  ('call', 'cpu', '24:00:00', 4, '16G', 'NDJcall'),
                  ('finalize', 'cpu', '04:00:00', 4, '12G', 'NDJconcat'),
                  ('plot', 'cpu', '02:00:00', 1, '8G', 'NDJaf')]
        for stage, partition, duration, cpus, memory, name in stages:
            command = [str(python), str(driver), stage, '--project', str(project)]
            sbatch = ['sbatch', '--parsable', '--chdir', str(project), '--job-name', name,
                      '--partition', partition, '--time', duration, '--cpus-per-task', str(cpus),
                      '--mem', memory, '--dependency', f'afterok:{dependency}',
                      '--output', str(logs / f'218_{stage}.%A_%a.out'),
                      '--error', str(logs / f'218_{stage}.%A_%a.err')]
            if stage == 'call':
                sbatch += ['--array', '0-20%21']
            sbatch += ['--wrap', shlex.join(command)]
            reply = subprocess.run(sbatch, check=False, capture_output=True, text=True)
            raw = reply.stdout.strip()
            # Record Slurm's response BEFORE attempting to interpret it.
            state.write(json.dumps({'stage': stage, 'afterok': dependency, 'submission_response': raw,
                                    'returncode': reply.returncode, 'stderr': reply.stderr,
                                    'sbatch': sbatch}) + '\n')
            state.flush()
            os.fsync(state.fileno())
            if reply.returncode:
                raise RuntimeError(f'sbatch failed for {stage}: {reply.stderr.strip()}')
            job = raw.split(';', 1)[0]
            if not job.isdigit():
                raise RuntimeError(f'Unexpected sbatch output: {raw!r}; inspect {state_file} before retrying')
            print(f'SUBMITTED\t{stage}\t{job}\tafterok:{dependency}', flush=True)
            jobs.append((stage, job))
            dependency = job
    print('\nCHAIN')
    for stage, job in jobs:
        print(f'{stage}\t{job}')
    subprocess.run(['squeue', '-j', ','.join(job for _, job in jobs)], check=False)


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        print(f'ERROR: {error}\nAny submitted jobs are preserved; inspect finish_chain.submissions.jsonl before retrying.', file=sys.stderr)
        raise SystemExit(1)
