#!/usr/bin/env python3
"""Plot an indexed-coverage-selected control against itself, or two real siblings."""
import argparse
from datetime import datetime, timezone
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import tempfile
import pysam

PROJECT = Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling')
SIBLINGS = ['SAMEA3594590', 'SAMEA3594614']  # F2 offspring 35 and 100
BCFTOOLS = '/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin/bcftools'


def stamp(path):
    path = Path(path)
    st = path.stat()
    return dict(path=str(path.resolve()), bytes=st.st_size, mtime_ns=st.st_mtime_ns)


def module(project):
    path = project / 'pairwise-af-hex45/plot_pairwise_AF.py'
    spec = importlib.util.spec_from_file_location('af_self', path)
    af = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(af)
    af.colourmap('batlowW_r')
    return af


def run_pipeline(first, second):
    a = subprocess.Popen(first, stdout=subprocess.PIPE)
    try:
        b = subprocess.Popen(second, stdin=a.stdout)
        a.stdout.close()
        second_code = b.wait()
        first_code = a.wait()
        if first_code or second_code:
            raise RuntimeError(f'BCF preparation failed: concat={first_code}, subset={second_code}')
    finally:
        if a.poll() is None:
            a.terminate()
            a.wait()


def sources(project, mode):
    original = project / 'analysis/figures/pairwise_AF_cli/TN_merged_control_sites.AF.audit.json'
    tn = json.loads(original.read_text())
    total = tn['records_retained']
    ndj = project / 'analysis/NDJ_control_comparison_v1/figures/NDJ_control_sites.AF.audit.json'
    if ndj.exists():
        total = max(total, json.loads(ndj.read_text())['records_retained'])
    vmax = 10 ** max(1, math.ceil(math.log10(total)))
    if mode == 'self':
        vcf = Path(tn['input']['path'])
        if stamp(vcf) != tn['input']:
            raise RuntimeError('TN input changed since the plotted audit')
        with pysam.VariantFile(str(vcf)) as handle:
            available = set(handle.header.samples)
        bamlist = project / 'analysis/genomewide_mixing_input/metadata/bam_sets_24/controls.bamlist'
        controls = [Path(line.strip()) for line in bamlist.read_text().splitlines() if line.strip()]
        if len(controls) != 18:
            raise RuntimeError(f'Expected 18 controls; found {len(controls)}')
        ranking = []
        for bam in controls:
            if not bam.is_absolute():
                bam = project / bam
            with pysam.AlignmentFile(str(bam), 'rb') as handle:
                handle.check_index()
                names = {rg.get('SM') for rg in handle.header.to_dict().get('RG', [])}
                if len(names) != 1 or None in names:
                    raise RuntimeError(f'Expected one sample name: {bam}')
                name = next(iter(names))
                if name not in available:
                    raise RuntimeError(f'Control missing from VCF: {name}')
                ranking.append(dict(sample=name, mapped_alignments=handle.mapped, bam=str(bam)))
        ranking.sort(key=lambda row: (-row['mapped_alignments'], row['sample']))
        if ranking[0]['mapped_alignments'] <= 0:
            raise RuntimeError('No mapped control reads')
        name = ranking[0]['sample']
        return dict(samples=[name, name], files=[str(vcf)], vmax=vmax,
                    ranking=ranking, selection='Largest indexed mapped-alignment count among 18 controls; coverage proxy, not purity/quality certification')
    folder = project / 'analysis/genomewide_mixing_corrected_v1/full_sib_validation/variant_calls/by_chromosome'
    files = []
    for chromosome in range(1, 22):
        matches = list(folder.glob(f'chr{chromosome}_*/full_sib_14.joint_at_existing_control_sites.bcf'))
        if len(matches) != 1 or matches[0].stat().st_size == 0:
            raise RuntimeError(f'Expected one nonempty sibling BCF for chromosome {chromosome}: {matches}')
        with pysam.VariantFile(str(matches[0])) as handle:
            if not set(SIBLINGS).issubset(handle.header.samples):
                raise RuntimeError(f'Sibling names absent from {matches[0]}')
        files.append(str(matches[0]))
    return dict(samples=SIBLINGS, files=files, vmax=vmax,
                selection='Known F2 offspring 35 and 100; real read counts without artificial mixing')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project', type=Path, default=PROJECT)
    parser.add_argument('--mode', required=True, choices=['self', 'siblings'])
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    project = args.project.resolve()
    af = module(project)
    selected = sources(project, args.mode)
    print('SELECTED\t' + args.mode + '\t' + ','.join(selected['samples']), flush=True)
    if args.mode == 'self':
        print('SELECTION_PROXY\tmapped_alignments=' + str(selected['ranking'][0]['mapped_alignments']), flush=True)
    root = project / 'analysis/figures/pairwise_AF_cli/control_benchmarks'
    prefix = root / ('best_coverage_control_self_hex45_batlowW' if args.mode == 'self'
                     else 'full_sibs_35_vs_100_hex45_batlowW')
    for suffix in ['AF.pdf', 'AF.png', 'AF.audit.json', 'AF_bins.tsv.gz']:
        if Path(str(prefix) + '.' + suffix).exists():
            raise RuntimeError(f'Output already exists: {prefix}')
    if args.mode == 'siblings':
        subprocess.run([BCFTOOLS, '--version'], stdout=subprocess.DEVNULL, check=True)
    if args.check:
        print('PREFLIGHT_PASS')
        return
    root.mkdir(parents=True, exist_ok=True)
    lock = root / (args.mode + '.lock')
    lock.mkdir()  # Protect preparation/publication against simultaneous reruns.
    try:
        source_stamps = [stamp(p) for p in selected['files']]
        vcf = Path(selected['files'][0])
        if args.mode == 'siblings':
            pairdir = root / 'siblings_35_100_input'
            if pairdir.exists():
                stored = json.loads((pairdir / 'sources.json').read_text())
                if stored['sources'] != source_stamps or stored['samples'] != SIBLINGS:
                    raise RuntimeError('Cached sibling input has different provenance')
                vcf = pairdir / 'siblings.bcf'
                if stamp(vcf) != stored['output']:
                    raise RuntimeError('Cached sibling BCF changed')
            else:
                with tempfile.TemporaryDirectory(prefix='.siblings-', dir=root) as tmp:
                    staged = Path(tmp) / 'input'
                    staged.mkdir()
                    subset = staged / 'siblings.bcf'
                    run_pipeline([BCFTOOLS, 'concat', '-Ou', *selected['files']],
                                 [BCFTOOLS, 'view', '-s', ','.join(SIBLINGS), '-Ob', '-o', str(subset)])
                    subprocess.run([BCFTOOLS, 'index', str(subset)], check=True)
                    if [stamp(p) for p in selected['files']] != source_stamps:
                        raise RuntimeError('Sibling inputs changed during preparation')
                    staged.rename(pairdir)
                vcf = pairdir / 'siblings.bcf'
                (pairdir / 'sources.json').write_text(json.dumps(dict(
                    sources=source_stamps, samples=SIBLINGS, output=stamp(vcf)), indent=2) + '\n')
        selected.update(source_fingerprints=source_stamps,
                        created_utc=datetime.now(timezone.utc).isoformat())
        Path(str(prefix) + '.selection.json').write_text(json.dumps(selected, indent=2) + '\n')
        title = ('Control self-plot (identical reads on both axes)' if args.mode == 'self'
                 else 'Full siblings: offspring 35 versus 100 (no added mixing)')
        argv = ['--vcf', str(vcf), '--left-sample', selected['samples'][0],
                '--right-sample', selected['samples'][1], '--output-prefix', str(prefix),
                '--bins', '45', '--bin-shape', 'hex', '--cmap', 'batlowW_r',
                '--min-depth', '10', '--filter', 'all', '--scale', 'count',
                '--vmin', '1', '--vmax', str(selected['vmax']), '--title', title]
        if args.mode == 'self':
            argv.append('--allow-self')
        result = af.main(argv)
        if result:
            raise RuntimeError('Plotting failed; see preceding error')
    finally:
        lock.rmdir()


if __name__ == '__main__':
    main()
