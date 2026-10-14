#!/usr/bin/env python3
"""Chromosome coverage bars from indexed BAMs; no VCF ascertainment."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import pysam

VERSION = '1.1.0'
SLOTS = ['left_t2', 'left_t3', 'right_t2', 'right_t3']
SEX = {'Z', 'W', 'X', 'Y'}


def write_json(path, obj):
    path = Path(path)
    tmp = path.with_name(path.name + '.partial')
    tmp.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


def stamp(path):
    path = Path(path).resolve(strict=True)
    stat = path.stat()
    return dict(path=str(path), size=stat.st_size, mtime_ns=stat.st_mtime_ns)


def chromosome_map(report, lengths):
    rows = []
    for line in Path(report).read_text().splitlines():
        if not line or line.startswith('#'):
            continue
        f = line.split('\t')
        if len(f) < 9 or f[1] != 'assembled-molecule':
            continue
        name = f[2]
        if name.upper() in {'MT', 'M', 'MITOCHONDRION'}:
            continue
        if not (name.isdigit() or name.upper() in SEX):
            raise ValueError(f'Unclassified assembled molecule {name}; label it explicitly before use')
        matches = list(dict.fromkeys(x for x in [f[6], f[4], f[0]] if x in lengths))
        if len(matches) != 1:
            raise ValueError(f'Ambiguous/missing reference contig for chromosome {name}: {matches}')
        contig = matches[0]
        if int(f[8]) != lengths[contig]:
            raise ValueError(f'Assembly length mismatch: {contig}')
        rows.append(dict(label=name.upper(), contig=contig, length=lengths[contig],
                         autosome=name.isdigit()))
    rows.sort(key=lambda r: (0, int(r['label'])) if r['autosome'] else (1, {'Z':0,'W':1,'X':2,'Y':3}[r['label']]))
    if not rows or not any(r['autosome'] for r in rows):
        raise ValueError('No autosomes resolved from assembly report')
    for key in ['label', 'contig']:
        if len({r[key] for r in rows}) != len(rows):
            raise ValueError(f'Duplicate chromosome {key}')
    return rows


def prepare(a):
    project = Path(a.project).resolve()
    root = project / 'analysis/genomewide_mixing_input'
    ref = Path(a.reference) if a.reference else Path((root / 'metadata/reference_fasta.path').read_text().strip())
    if not ref.is_absolute():
        ref = project / ref
    ref = ref.resolve(strict=True)
    if a.assembly_report:
        report = Path(a.assembly_report).resolve(strict=True)
    else:
        reports = list(ref.parent.glob('*assembly_report.txt'))
        if len(reports) != 1:
            raise ValueError('Provide --assembly-report: expected one report beside reference')
        report = reports[0]
    with pysam.FastaFile(str(ref)) as fa:
        lengths = dict(zip(fa.references, fa.lengths))
    chromosomes = chromosome_map(report, lengths)
    names = [getattr(a, slot + '_sample') for slot in SLOTS]
    if len(set(names)) != 4:
        raise ValueError('Four distinct T2/T3 sample names are required')
    explicit = [getattr(a, slot + '_bam') for slot in SLOTS]
    if any(explicit) and not all(explicit):
        raise ValueError('Supply all four explicit BAM paths, or none')
    if all(explicit):
        selected = [Path(value) for value in explicit]
    else:
        manifest = Path(a.bamlist) if a.bamlist else root / 'metadata/bam_sets_24/all_24_samples.bamlist'
        found = {name: [] for name in names}
        for line in manifest.read_text().splitlines():
            if not line.strip():
                continue
            path = Path(line.strip())
            if not path.is_absolute():
                path = project / path
            with pysam.AlignmentFile(str(path), 'rb') as bam:
                samples = {r.get('SM') for r in bam.header.to_dict().get('RG', [])}
            for name in names:
                if samples == {name}:
                    found[name].append(path)
        if any(len(found[name]) != 1 for name in names):
            raise ValueError(f'Need exactly one BAM per sample, found {found}; supply explicit BAMs if necessary')
        selected = [found[name][0] for name in names]
    selected = [p.resolve(strict=True) for p in selected]
    if len(set(selected)) != 4:
        raise ValueError('Four distinct BAMs are required')
    for name, path in zip(names, selected):
        with pysam.AlignmentFile(str(path), 'rb') as bam:
            bam.check_index()
            h = bam.header.to_dict()
            if h.get('HD', {}).get('SO') != 'coordinate':
                raise ValueError(f'BAM must be coordinate sorted: {path}')
            if {r.get('SM') for r in h.get('RG', [])} != {name}:
                raise ValueError(f'BAM read-group sample mismatch: {name}: {path}')
            if dict(zip(bam.references, bam.lengths)) != lengths:
                raise ValueError(f'BAM/reference dictionary mismatch: {path}')
    samtools = a.samtools or shutil.which('samtools')
    if not samtools:
        candidate = Path.home() / '.local/share/mamba/envs/colias_w_popgen/bin/samtools'
        if candidate.is_file():
            samtools = str(candidate)
    if not samtools:
        raise ValueError('samtools not found; supply --samtools /absolute/path')
    samtools = str(Path(shutil.which(samtools) or samtools).resolve(strict=True))
    version = subprocess.check_output([samtools, '--version'], text=True).splitlines()[0]
    match = re.search(r'(\d+)\.(\d+)', version)
    if not match or tuple(map(int, match.groups())) < (1, 13):
        raise ValueError('samtools >=1.13 required for consistent overlapping-mate handling')
    run = Path(a.run_dir).resolve() if a.run_dir else project / 'analysis/chromosome_coverage_bars' / datetime.now(timezone.utc).strftime('run_%Y%m%dT%H%M%S%fZ')
    if run.exists():
        raise ValueError(f'Run already exists; use scan/plot to resume: {run}')
    run.mkdir(parents=True)
    (run / 'logs').mkdir()
    (run / 'chromosomes').mkdir()
    cfg = dict(version=VERSION, reference=stamp(ref), assembly_report=stamp(report),
               bams=[stamp(p) for p in selected], samples=names, slots=SLOTS, chromosomes=chromosomes,
               samtools=samtools, samtools_version=version, mapq=a.mapq, baseq=a.baseq,
               window=a.window, exclude_flags='0xF04', remove_mate_overlap=True,
               denominator='All A/C/G/T reference bases, including zero-depth positions',
               normalization='Per-sample median autosomal window mean depth; windows >=90% ACGT')
    cfg['signature'] = hashlib.sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest()
    write_json(run / 'config.json', cfg)
    script = run / 'coverage_bars.py'
    shutil.copyfile(Path(__file__).resolve(), script)
    print('PREFLIGHT_PASS\t' + str(run), flush=True)
    for slot, name, bam in zip(SLOTS, names, selected):
        print(f'SAMPLE\t{slot}\t{name}\t{bam}', flush=True)
    print('CHROMOSOMES\t' + ','.join(r['label'] for r in chromosomes), flush=True)
    print('W_PRESENT\t' + str(any(r['label'] == 'W' for r in chromosomes)), flush=True)
    if a.submit:
        def submit(stage, options):
            command = shlex.join([sys.executable, str(script), stage, '--run', str(run)])
            if stage == 'scan':
                command += ' --index "${SLURM_ARRAY_TASK_ID}"'
            job = ['sbatch', '--parsable', '--job-name=cov_' + stage,
                   '--cpus-per-task=1', '--mem=4G', '--time=02:00:00',
                   '--output=' + str(run / 'logs' / (stage + '.%A_%a.out')),
                   '--error=' + str(run / 'logs' / (stage + '.%A_%a.err'))] + options + ['--wrap', 'set -eu\n' + command]
            result = subprocess.run(job, text=True, capture_output=True)
            if result.returncode:
                raise RuntimeError(result.stderr)
            identifier = result.stdout.strip().split(';')[0]
            if not identifier.isdigit():
                raise RuntimeError(f'Unexpected sbatch response: {result.stdout}')
            with (run / 'submissions.jsonl').open('a') as handle:
                handle.write(json.dumps(dict(stage=stage, job=identifier, command=job)) + '\n')
            print(stage.upper() + '_JOB\t' + identifier, flush=True)
            return identifier
        array = submit('scan', ['--partition=' + a.partition, f'--array=0-{len(chromosomes)-1}%{a.max_parallel}'])
        submit('plot', ['--partition=' + a.plot_partition, '--dependency=afterok:' + array])
    print('RUN_DIRECTORY=' + str(run))


def load(run):
    run = Path(run).resolve()
    cfg = json.loads((run / 'config.json').read_text())
    return run, cfg


def scan(a):
    run, cfg = load(a.run)
    for entry in [cfg['reference'], *cfg['bams']]:
        if stamp(entry['path']) != entry:
            raise ValueError('Input changed since preflight: ' + entry['path'])
    row = cfg['chromosomes'][a.index]
    dest = run / 'chromosomes' / f'{a.index:03d}.tsv'
    done = dest.with_suffix('.complete.json')
    if done.exists() and dest.exists():
        if json.loads(done.read_text())['signature'] == cfg['signature']:
            print('REUSE\t' + row['label'])
            return
    with pysam.FastaFile(cfg['reference']['path']) as fa:
        seq = fa.fetch(row['contig']).upper()
    width = cfg['window']
    starts = np.arange(0, len(seq), width)
    ends = np.minimum(starts + width, len(seq))
    valid = np.isin(np.frombuffer(seq.encode('ascii'), dtype='S1'), [b'A', b'C', b'G', b'T'])
    bases = np.add.reduceat(valid.astype(np.int64), starts)
    sums = np.zeros((len(starts), len(cfg["slots"])), dtype=np.int64)
    covered = np.zeros_like(sums)
    command = [cfg['samtools'], 'depth', '-aa', '-s', '-q', str(cfg['baseq']),
               '-Q', str(cfg['mapq']), '-G', cfg['exclude_flags'], '-r', row['contig'],
               *[b['path'] for b in cfg['bams']]]
    previous = 0
    with (run / 'logs' / f'depth.{a.index:03d}.stderr').open('w') as err:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=err, text=True)
        try:
            for line in process.stdout:
                fields = line.split()
                if len(fields) != 2 + len(cfg["slots"]) or fields[0] != row['contig']:
                    raise ValueError('Unexpected samtools depth record')
                pos = int(fields[1])
                depths = list(map(int, fields[2:]))
                if pos != previous + 1 or pos > len(seq):
                    raise ValueError('Missing, duplicate, or out-of-range depth positions')
                previous = pos
                if valid[pos - 1]:
                    index = (pos - 1) // width
                    for sample_index, depth in enumerate(depths):
                        sums[index, sample_index] += depth
                        covered[index, sample_index] += depth > 0
            if process.wait() != 0:
                raise RuntimeError('samtools depth failed; inspect depth stderr')
        finally:
            process.stdout.close()
            if process.poll() is None:
                process.kill()
                process.wait()
    if previous != row['length']:
        raise ValueError(f'Depth ended at {previous}; expected {row["length"]}')
    table = pd.DataFrame(dict(chrom=row['contig'], chromosome=row['label'], start=starts,
                              end=ends, acgt_bases=bases))
    for j, slot in enumerate(cfg['slots']):
        table[slot + '_sum'] = sums[:, j]
        table[slot + '_covered_bases'] = covered[:, j]
    tmp = dest.with_suffix('.partial')
    table.to_csv(tmp, sep='\t', index=False)
    tmp.replace(dest)
    write_json(done, dict(signature=cfg['signature'], chromosome=row, command=command))
    print('COMPLETE\t' + row['label'], flush=True)


def plot(a):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    run, cfg = load(a.run)
    z_rows = [
        row for row in cfg['chromosomes']
        if row['contig'] == 'NC_140409.1'
    ]
    if len(z_rows) != 1:
        raise RuntimeError("Expected project Z contig NC_140409.1")
    z_rows[0]['label'] = 'Z'
    z_rows[0]['autosome'] = False
    write_json(run / 'mirrored_plot_chromosome_map.json',
               cfg['chromosomes'])
    tables = []
    for i, row in enumerate(cfg['chromosomes']):
        path = run / 'chromosomes' / f'{i:03d}.tsv'
        completion = json.loads(path.with_suffix('.complete.json').read_text())
        if completion['signature'] != cfg['signature']:
            raise ValueError(f'Wrong-run depth output: {path}')
        t = pd.read_csv(path, sep='\t', dtype={'chromosome':str})
        t['autosome'] = row['autosome']
        t['chromosome'] = row['label']
        tables.append(t)
    table = pd.concat(tables, ignore_index=True)
    eligible = table.autosome & (table.acgt_bases / (table.end - table.start) >= .9)
    if eligible.sum() < 5:
        raise ValueError('Fewer than five usable autosomal normalization windows')
    baseline = {}
    for side in cfg['slots']:
        table[side + '_mean_depth'] = table[side + '_sum'] / table.acgt_bases.replace(0, np.nan)
        baseline[side] = float(table.loc[eligible, side + '_mean_depth'].median())
        if not np.isfinite(baseline[side]) or baseline[side] <= 0:
            raise ValueError('Invalid autosomal baseline: ' + side)
        table[side + '_relative_depth'] = table[side + '_mean_depth'] / baseline[side]
    rows = []
    for info in cfg['chromosomes']:
        subset = table.loc[table.chrom == info['contig']]
        bases = int(subset.acgt_bases.sum())
        if not bases:
            raise ValueError('No ACGT bases: ' + info['contig'])
        row = dict(info, acgt_bases=bases)
        for side in cfg['slots']:
            depth = float(subset[side + '_sum'].sum()) / bases
            row[side + '_mean_depth'] = depth
            row[side + '_relative_depth'] = depth / baseline[side]
            row[side + '_covered_fraction'] = float(subset[side + '_covered_bases'].sum()) / bases
        rows.append(row)
    summary = pd.DataFrame(rows)
    table.to_csv(run / 'coverage_windows.tsv', sep='\t', index=False)
    summary.to_csv(run / 'chromosome_coverage.tsv', sep='\t', index=False)
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':8, 'pdf.fonttype':42,
                         'axes.linewidth':.6, 'xtick.labelsize':7, 'ytick.labelsize':8})
    fig, axes = plt.subplots(1, 2, sharey=True, figsize=(3.4, max(3.5, .15*len(rows)+1.0)))
    fig.subplots_adjust(left=.15, right=.97, top=.79, bottom=.15, wspace=.25)
    y = np.arange(len(rows))
    maximum = max(1.6, float(summary[[slot + '_relative_depth' for slot in cfg['slots']]].to_numpy().max()) * 1.08)
    for ax, side, colours, title in zip(axes, ['left','right'],
            [('#65429A','#9981C3'), ('#E87519','#F6AD60')], ['Left','Right']):
        for suffix, offset, colour in zip(['t2','t3'], [-.19,.19], colours):
            ax.barh(y + offset, summary[side + '_' + suffix + '_relative_depth'],
                    color=colour, height=.34, linewidth=0, label=suffix.upper())
        ax.legend(loc='lower center', bbox_to_anchor=(.5,1.015), ncol=2,
                  frameon=False, fontsize=6, handlelength=1, columnspacing=.8)
        for level in [.5, 1, 1.5]:
            ax.axvline(level, color='#444444', linewidth=.5, linestyle='--', alpha=.6, zorder=3)
        ax.set(xlim=(0, maximum), ylim=(len(rows)-.4, -.6), title=title)
        ax.set_title(title, pad=26)
        ax.set_xticks([0, .5, 1, 1.5])
        ax.set_yticks(y)
        ax.set_yticklabels(summary.label)
        ax.tick_params(axis='y', length=0)
        for spine in ['top','right','left']:
            ax.spines[spine].set_visible(False)
    # Shared chromosome positions; labels in the central gutter.
    for ax in axes:
        ax.tick_params(axis='y', left=False, right=False,
                       labelleft=False, labelright=False)

    # Left zero is at the central edge, increasing towards the left.
    axes[0].set_xlim(maximum, 0)
    axes[1].set_xlim(0, maximum)

    left_box = axes[0].get_position()
    right_box = axes[1].get_position()
    middle = (left_box.x1 + right_box.x0) / 2

    for position, label in zip(y, summary.label):
        figure_y = fig.transFigure.inverted().transform(
            axes[0].transData.transform((0, position))
        )[1]
        fig.text(middle, figure_y, str(label),
                 ha='center', va='center', fontsize=8)

    # Save numerical Z coverage for interpretation.
    z_summary = summary.loc[summary['label'].eq('Z')].copy()
    z_summary.to_csv(run / 'Z_relative_coverage.tsv',
                     sep='\t', index=False)
    print("\\nZ RELATIVE COVERAGE", flush=True)
    for slot, sample in zip(cfg['slots'], cfg['samples']):
        value = float(z_summary.iloc[0][slot + '_relative_depth'])
        print(f"{slot}\\t{sample}\\t{value:.6f}", flush=True)
    fig.suptitle('Sequencing coverage', y=.97, fontsize=10, fontweight='bold')
    fig.text(.56, .035, 'Relative depth (autosomal baseline = 1)', ha='center', fontsize=8)
    fig.savefig(run / 'chromosome_coverage_mirrored.pdf')
    fig.savefig(run / 'chromosome_coverage_mirrored.png', dpi=600)
    plt.close(fig)
    write_json(run / 'COMPLETE.json', dict(version=VERSION, signature=cfg['signature'],
               baseline=baseline, chromosomes=len(rows), samples=cfg['samples'],
               W_present=any(r['label']=='W' for r in rows),
               interpretation='Relative coverage; not absolute ploidy. No GC/mappability correction or mixing correction.',
               pdf=str(run / 'chromosome_coverage_mirrored.pdf')))
    print('PDF\t' + str(run / 'chromosome_coverage_mirrored.pdf'))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    q = sub.add_parser('prepare')
    q.add_argument('--project', default='.')
    q.add_argument('--reference')
    q.add_argument('--assembly-report')
    q.add_argument('--bamlist')
    for slot, default in zip(SLOTS, ['TN2L_S1','TN3L_S3','TN2R_S2','TN3R_S4']):
        flag = slot.replace('_', '-')
        q.add_argument('--' + flag + '-bam')
        q.add_argument('--' + flag + '-sample', default=default)
    q.add_argument('--samtools')
    q.add_argument('--window', type=int, default=100000)
    q.add_argument('--mapq', type=int, default=20)
    q.add_argument('--baseq', type=int, default=20)
    q.add_argument('--run-dir')
    q.add_argument('--submit', action='store_true')
    q.add_argument('--partition', default='cpu')
    q.add_argument('--plot-partition', default='nano')
    q.add_argument('--max-parallel', type=int, default=8)
    q = sub.add_parser('scan')
    q.add_argument('--run', required=True)
    q.add_argument('--index', type=int, required=True)
    q = sub.add_parser('plot')
    q.add_argument('--run', required=True)
    a = p.parse_args()
    if a.command == 'prepare' and (a.window <= 0 or a.max_parallel <= 0 or min(a.mapq, a.baseq) < 0):
        p.error('Window/concurrency must be positive; quality thresholds nonnegative')
    {'prepare':prepare, 'scan':scan, 'plot':plot}[a.command](a)


if __name__ == '__main__':
    main()
