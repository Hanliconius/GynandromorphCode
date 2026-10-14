#!/usr/bin/env python3
"""Regenerate both AF figures using saved input/filter audits; no calling."""
import argparse
from datetime import datetime, timezone
import importlib.util
import json
import math
import os
from pathlib import Path
import shlex
import subprocess
import sys

AUDITS = [
    ('TN', 'analysis/figures/pairwise_AF_cli/TN_merged_control_sites.AF.audit.json'),
    ('NDJ', 'analysis/NDJ_control_comparison_v1/figures/NDJ_control_sites.AF.audit.json'),
]


def load_tool(project):
    path = project / 'pairwise-af/plot_pairwise_AF.py'
    spec = importlib.util.spec_from_file_location('pairwise_af_regenerate', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if module.__version__ != '0.2.0':
        raise RuntimeError('Install pairwise-af v0.2.0 before regenerating')
    module.colourmap('batlowW_r')
    return module


def make_plan(project, tool):
    entries = []
    for label, relative in AUDITS:
        path = project / relative
        if not path.is_file():
            raise RuntimeError(f'Original plot must finish first: {path}')
        audit = json.loads(path.read_text())
        params = audit['parameters']
        vcf = Path(audit['input']['path'])
        if tool.fingerprint(vcf) != audit['input']:
            raise RuntimeError(f'Input changed since original plot: {vcf}')
        # Reproduce every scientific selection from the original plot.
        prefix = Path(str(path)[:-len('.AF.audit.json')] + '_hex65_batlowW')
        for suffix in ('AF.pdf', 'AF.png', 'AF_bins.tsv.gz', 'AF.audit.json'):
            if Path(str(prefix) + '.' + suffix).exists():
                raise RuntimeError(f'New output already exists: {prefix}')
        command = ['--vcf', str(vcf), '--output-prefix', str(prefix),
                   '--bins', '65', '--bin-shape', 'hex', '--cmap', 'batlowW_r',
                   '--scale', 'count', '--vmin', '1']
        for key in ('left_sample', 'right_sample', 'min_depth', 'filter',
                    'max_depth', 'mixing_left', 'mixing_right', 'title'):
            value = params.get(key)
            if value is not None:
                command.extend(['--' + key.replace('_', '-'), str(value)])
        for contig in params.get('contig') or []:
            command.extend(['--contig', contig])
        whitelist = audit.get('whitelist')
        if params.get('sites'):
            if not whitelist:
                raise RuntimeError(f'Missing whitelist provenance: {path}')
            sites = Path(whitelist['path'])
            if tool.sha256(sites) != whitelist['sha256']:
                raise RuntimeError(f'Whitelist changed: {sites}')
            command.extend(['--sites', str(sites)])
        # vmax added after both retained totals are known.
        entries.append({'label': label, 'old_audit': str(path),
                        'old_audit_sha256': tool.sha256(path),
                        'output_prefix': str(prefix), 'argv': command,
                        'retained': audit['records_retained']})
    # Common unclipped count scale. Bound by the larger TOTAL, so every cell fits.
    vmax = float(10 ** max(1, math.ceil(math.log10(max(e['retained'] for e in entries)))))
    for entry in entries:
        entry['argv'].extend(['--vmax', str(vmax)])
        tool.check_args(tool.parser().parse_args(entry['argv']))
    return {'project': str(project), 'vmin': 1, 'vmax': vmax, 'entries': entries,
            'tool_sha256': tool.sha256(project / 'pairwise-af/plot_pairwise_AF.py'),
            'scale_note': 'Shared log count scale bounded by larger original retained total; not per-sample auto scaling.'}


def run_task(plan_path, task):
    plan = json.loads(plan_path.read_text())
    project = Path(plan['project'])
    tool = load_tool(project)
    if tool.sha256(project / 'pairwise-af/plot_pairwise_AF.py') != plan['tool_sha256']:
        raise RuntimeError('Plot tool changed after submission')
    entry = plan['entries'][task]
    oldpath = Path(entry['old_audit'])
    if tool.sha256(oldpath) != entry['old_audit_sha256']:
        raise RuntimeError('Original audit changed after submission')
    old = json.loads(oldpath.read_text())
    if tool.fingerprint(Path(old['input']['path'])) != old['input']:
        raise RuntimeError('VCF changed after submission')
    if old.get('whitelist'):
        white = old['whitelist']
        if tool.sha256(Path(white['path'])) != white['sha256']:
            raise RuntimeError('Whitelist changed after submission')
    code = tool.main(entry['argv'])
    if code:
        raise RuntimeError(f'Plot failed for {entry["label"]}')
    new = json.loads(Path(entry['output_prefix'] + '.AF.audit.json').read_text())
    for key in ('records_read', 'records_retained', 'excluded_first_reason', 'retained_by_contig'):
        if old[key] != new[key]:
            raise RuntimeError(f'Selection changed for {entry["label"]}: {key}; inspect audits')
    print(f'SELECTION_MATCH\t{entry["label"]}\tretained={new["records_retained"]}', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project', type=Path, default=Path.cwd())
    parser.add_argument('--plan', type=Path)
    parser.add_argument('--task', type=int, choices=[0, 1])
    args = parser.parse_args()
    if args.task is not None:
        if args.plan is None:
            parser.error('--task requires --plan')
        run_task(args.plan, args.task)
        return
    project = args.project.resolve()
    tool = load_tool(project)
    plan = make_plan(project, tool)
    logs = project / 'logs'
    logs.mkdir(exist_ok=True)
    tag = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    plan_path = logs / f'219_AF_hex65.{tag}.plan.json'
    plan_path.write_text(json.dumps(plan, indent=2) + '\n')
    python = project / 'pairwise-af/.venv-af/bin/python'
    if not python.is_file():
        raise RuntimeError(f'Missing plotting environment: {python}')
    runner = project / 'scripts/219_regenerate_TN_NDJ_AF.py'
    # One array, two independent tasks. Shell expansion only for Slurm task index.
    wrapped = shlex.join([str(python), str(runner), '--plan', str(plan_path)])
    wrapped += ' --task "$SLURM_ARRAY_TASK_ID"'
    cmd = ['sbatch', '--parsable', '--chdir', str(project), '--job-name', 'AFhex65',
           '--partition', 'cpu', '--time', '02:00:00', '--cpus-per-task', '1',
           '--mem', '8G', '--array', '0-1%2',
           '--output', str(logs / '219_AF_hex65.%A_%a.out'),
           '--error', str(logs / '219_AF_hex65.%A_%a.err'), '--wrap', wrapped]
    reply = subprocess.run(cmd, text=True, capture_output=True)
    plan_path.with_suffix('.submission.json').write_text(json.dumps({
        'command': cmd, 'returncode': reply.returncode,
        'stdout': reply.stdout, 'stderr': reply.stderr}, indent=2) + '\n')
    if reply.returncode:
        raise RuntimeError(reply.stderr.strip())
    job = reply.stdout.strip().split(';', 1)[0]
    print(f'SUBMITTED\t{job}\tarray=0-1%2\tshared_colour_max={plan["vmax"]:g}')
    for entry in plan['entries']:
        print(entry['output_prefix'] + '.AF.pdf')
    if job.isdigit():
        subprocess.run(['squeue', '-j', job], check=False)


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        print(f'ERROR: {error}', file=sys.stderr)
        raise SystemExit(1)
