#!/usr/bin/env python3
"""Add NDJ to the 18-control cohort at existing control-discovered sites."""
import argparse
from contextlib import contextmanager
import gzip
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile

import pysam

NAMES = ['NDJ_left', 'NDJ_right']
BCF_NAME = 'NDJ_controls.joint_at_control_sites.bcf'


def require(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f'Missing or empty file: {path}')
    return path


def stamp(path):
    path = require(path).resolve()
    stat = path.stat()
    return {'path': str(path), 'bytes': stat.st_size, 'mtime_ns': stat.st_mtime_ns}


def unchanged(item):
    if stamp(item['path']) != item:
        raise RuntimeError(f"Input changed: {item['path']}")


def run(cmd, **kwargs):
    print('COMMAND\t' + shlex.join(map(str, cmd)), flush=True)
    return subprocess.run(list(map(str, cmd)), check=True, **kwargs)


def pipeline(first, second):
    print('PIPELINE\t' + shlex.join(first) + ' | ' + shlex.join(second), flush=True)
    producer = subprocess.Popen(first, stdout=subprocess.PIPE)
    try:
        consumer = subprocess.Popen(second, stdin=producer.stdout)
    except BaseException:
        producer.stdout.close()
        producer.terminate()
        producer.wait()
        raise
    producer.stdout.close()
    try:
        consumer_status = consumer.wait()
        producer_status = producer.wait()
    except BaseException:
        for process in (consumer, producer):
            if process.poll() is None:
                process.terminate()
        consumer.wait()
        producer.wait()
        raise
    if producer_status or consumer_status:
        raise RuntimeError(f'Calling pipeline failed: mpileup={producer_status}, call={consumer_status}')


@contextmanager
def lock(path):
    path.mkdir()  # Atomic, fails on a second active submission.
    try:
        yield
    finally:
        path.rmdir()


def dump(path, data):
    path.write_text(json.dumps(data, indent=2) + '\n')


def targets(path, contig):
    """Actual project format: CHROM, POS, REF,ALT (THREE columns)."""
    previous = 0
    with gzip.open(path, 'rt') as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 3:
                raise RuntimeError(f'{path}:{number}: expected three columns CHROM POS REF,ALT')
            chrom, position, allele_text = fields
            position = int(position)
            alleles = allele_text.split(',')
            if (chrom != contig or position <= previous or len(alleles) != 2
                    or any(a not in {'A', 'C', 'G', 'T'} for a in alleles)
                    or alleles[0] == alleles[1]):
                raise RuntimeError(f'{path}:{number}: invalid, duplicate or unsorted target')
            previous = position
            yield position, alleles[0], alleles[1]


def discover(project):
    source = project / 'analysis/genomewide_mixing_input'
    reference = Path(require(source / 'metadata/reference_fasta.path').read_text().strip())
    if not reference.is_absolute():
        reference = project / reference
    require(reference)
    require(str(reference) + '.fai')
    control_file = require(source / 'metadata/bam_sets_24/controls.bamlist')
    bams = []
    for line in control_file.read_text().splitlines():
        if not line.strip():
            continue
        bam = Path(line.strip())
        if not bam.is_absolute():
            bam = project / bam
        bams.append(require(bam).resolve())
    if len(bams) != 18 or len(set(bams)) != 18:
        raise RuntimeError(f'Expected 18 distinct control BAMs; found {len(bams)}')
    chromosomes = []
    with open(str(reference) + '.fai') as handle:
        refdict = {f[0]: int(f[1]) for f in (line.split('\t') for line in handle)}
    for chrom in range(1, 22):
        matches = list((source / 'variant_calls/joint_24/by_chromosome').glob(f'chr{chrom}_*'))
        matches = [p for p in matches if p.is_dir()]
        if len(matches) != 1:
            raise RuntimeError(f'Expected one source directory for chromosome {chrom}; found {matches}')
        directory = matches[0]
        contig = directory.name.split('_', 1)[1]
        if contig not in refdict:
            raise RuntimeError(f'Contig absent from active reference: {contig}')
        allele_file = require(directory / 'controls.discovery.alleles.tsv.gz')
        indexes = [Path(str(allele_file) + suffix) for suffix in ('.tbi', '.csi')]
        if not any(p.is_file() and p.stat().st_size > 0 for p in indexes):
            raise RuntimeError(f'Missing target-list tabix index: {allele_file}; do not start calling until this is resolved')
        generator = targets(allele_file, contig)
        try:
            if next(generator, None) is None:
                raise RuntimeError(f'Empty target list: {allele_file}')
        finally:
            generator.close()
        chromosomes.append({'chromosome': chrom, 'contig': contig, 'alleles': stamp(allele_file)})
    return reference, refdict, control_file, bams, chromosomes


def prepare(args):
    reference, refdict, control_file, controls, chromosomes = discover(args.project)
    print(f'CONTROL_INPUTS_PASS\tbams={len(controls)}\tchromosomes={len(chromosomes)}', flush=True)
    if args.stage == 'preflight':
        # NDJ BAMs are future dependency outputs: do not require them at submission.
        run([args.bcftools, '--version'], stdout=subprocess.DEVNULL)
        run([args.samtools, '--version'], stdout=subprocess.DEVNULL)
        require(args.project / 'pairwise-af/plot_pairwise_AF.py')
        return
    ndj = [args.root / f'alignment/{name}/{name}.dedup.bam' for name in NAMES]
    for bam in ndj:
        complete = json.loads(require(bam.parent / 'alignment.complete.json').read_text())
        if complete['sample'] != bam.parent.name or Path(complete['reference']).resolve() != reference.resolve():
            raise RuntimeError(f'NDJ alignment provenance mismatch: {bam}')
    bamlist = controls + ndj
    samples = []
    bam_inputs = []
    for bam in bamlist:
        require(bam)
        run([args.samtools, 'quickcheck', '-v', bam])
        with pysam.AlignmentFile(str(bam), 'rb') as handle:
            if dict(zip(handle.references, handle.lengths)) != refdict:
                raise RuntimeError(f'Reference dictionary mismatch: {bam}')
            groups = handle.header.to_dict().get('RG', [])
            labels = {g.get('SM') for g in groups}
            if not groups or len(labels) != 1 or None in labels or '' in labels:
                raise RuntimeError(f'Expected one nonempty SM per BAM: {bam}: {labels}')
            samples.append(next(iter(labels)))
            handle.check_index()
            if handle.mapped == 0:
                raise RuntimeError(f'No mapped alignments in {bam}')
        bam_inputs.append(stamp(bam))
    if len(set(samples)) != 20 or samples[-2:] != NAMES or any(s.startswith('TN_') for s in samples[:18]):
        raise RuntimeError(f'Expected 18 unique controls and {NAMES}; found {samples}')
    config = {'reference': stamp(reference), 'fai': stamp(str(reference) + '.fai'),
              'control_bamlist': stamp(control_file), 'samples': samples, 'bams': bam_inputs,
              'chromosomes': chromosomes, 'software': {
                  'bcftools': run([args.bcftools, '--version'], capture_output=True, text=True).stdout,
                  'pysam': pysam.__version__}}
    args.root.mkdir(parents=True, exist_ok=True)
    with lock(args.root / '.prepare.lock'):
        final = args.root / 'metadata'
        if final.exists():
            previous = json.loads(require(final / 'cohort.json').read_text())
            if previous != config:
                raise RuntimeError('Existing cohort metadata differs; refusing overwrite')
            print('REUSE\tcohort metadata')
            return
        with tempfile.TemporaryDirectory(prefix='.prepare.', dir=args.root) as tmp:
            output = Path(tmp) / 'metadata'
            output.mkdir()
            dump(output / 'cohort.json', config)
            (output / 'all_20_samples.bamlist').write_text(''.join(str(p) + '\n' for p in bamlist))
            (output / 'all_20_samples.samples').write_text(''.join(s + '\n' for s in samples))
            output.rename(final)
    print('PREPARE_COMPLETE\tcontrols=18\tNDJ=2\ttotal=20', flush=True)


def get_config(args):
    config = json.loads(require(args.root / 'metadata/cohort.json').read_text())
    for item in [config['reference'], config['fai'], config['control_bamlist'], *config['bams']]:
        unchanged(item)
    expected = ''.join(b['path'] + '\n' for b in config['bams'])
    if require(args.root / 'metadata/all_20_samples.bamlist').read_text() != expected:
        raise RuntimeError('BAM list differs from validated cohort metadata')
    return config


def validate_bcf(path, alleles, contig, samples):
    generator = targets(alleles, contig)
    target = next(generator, None)
    total_targets = int(target is not None)
    rows = 0
    covered = {name: 0 for name in samples[-2:]}
    last_position = 0
    digest = hashlib.sha256()
    try:
        with pysam.VariantFile(str(path)) as vcf:
            if list(vcf.header.samples) != samples:
                raise RuntimeError('BCF sample names/order differ from validated BAM cohort')
            if 'AD' not in vcf.header.formats:
                raise RuntimeError('Missing FORMAT/AD in BCF')
            for record in vcf:
                if record.contig != contig or record.pos <= last_position:
                    raise RuntimeError('Unexpected contig or duplicate/unsorted emitted position')
                last_position = record.pos
                while target is not None and target[0] < record.pos:
                    target = next(generator, None)
                    total_targets += int(target is not None)
                if target is None or (record.pos, record.ref, record.alts) != (target[0], target[1], (target[2],)):
                    raise RuntimeError(f'Emitted allele absent from target list: {record.contig}:{record.pos}:{record.ref}:{record.alts}')
                digest.update(f'{contig}\t{record.pos}\t{record.ref}\t{record.alts[0]}\n'.encode())
                for name in covered:
                    ad = record.samples[name].get('AD')
                    if ad and len(ad) == 2 and all(v is not None and v >= 0 for v in ad) and sum(ad) > 0:
                        covered[name] += 1
                rows += 1
        for _ in generator:
            total_targets += 1
    finally:
        generator.close()
    if not rows or not all(covered.values()):
        raise RuntimeError(f'No records or no usable NDJ AD: records={rows}, NDJ_covered={covered}')
    return {'contig': contig, 'target_sites': total_targets, 'emitted_sites': rows,
            'omitted_target_sites': total_targets - rows, 'NDJ_nonzero_AD_sites': covered,
            'emitted_variant_keys_sha256': digest.hexdigest(),
            'omission_policy': 'Exact emitted allele subset required; missing target records are reported, never converted to 0/0'}


def call_chromosome(args):
    config = get_config(args)
    task = int(os.environ['SLURM_ARRAY_TASK_ID'])
    if not 0 <= task < 21:
        raise RuntimeError(f'Invalid chromosome array index: {task}')
    row = config['chromosomes'][task]
    unchanged(row['alleles'])
    parent = args.root / 'variant_calls/by_chromosome'
    parent.mkdir(parents=True, exist_ok=True)
    name = f"chr{row['chromosome']}_{row['contig']}"
    final = parent / name
    if final.exists():
        raise RuntimeError(f'Chromosome output already exists; preserved: {final}')
    with lock(parent / ('.' + name + '.lock')):
        with tempfile.TemporaryDirectory(prefix='.' + name + '.', dir=parent) as tmp:
            output = Path(tmp) / 'publish'
            output.mkdir()
            bcf = output / BCF_NAME
            first = [str(args.bcftools), 'mpileup', '-Ou', '--threads', '1', '-f', config['reference']['path'],
                     '-r', row['contig'], '-T', row['alleles']['path'], '-q', '20', '-Q', '20', '-d', '1000',
                     '-a', 'FORMAT/AD,FORMAT/DP,FORMAT/SP', '-b', str(args.root / 'metadata/all_20_samples.bamlist')]
            second = [str(args.bcftools), 'call', '-Ob', '--threads', '1', '-m', '-A', '-C', 'alleles',
                      '-T', row['alleles']['path'], '-G', '-', '-o', str(bcf)]
            pipeline(first, second)
            run([args.bcftools, 'index', '-f', bcf])
            audit = validate_bcf(bcf, row['alleles']['path'], row['contig'], config['samples'])
            indexed = int(run([args.bcftools, 'index', '-n', bcf], capture_output=True, text=True).stdout.strip())
            if indexed != audit['emitted_sites']:
                raise RuntimeError('BCF index record count mismatch')
            audit.update({'chromosome': row['chromosome'], 'samples': config['samples'],
                          'mpileup_command': first, 'call_command': second,
                          'bcf_bytes': bcf.stat().st_size})
            unchanged(row['alleles'])
            dump(output / 'call.complete.json', audit)
            output.rename(final)
    print(f"CALL_COMPLETE\tchromosome={row['chromosome']}\temitted={audit['emitted_sites']}\tomitted={audit['omitted_target_sites']}", flush=True)


def finalize(args):
    config = get_config(args)
    bcf_files, audits = [], []
    for row in config['chromosomes']:
        folder = args.root / f"variant_calls/by_chromosome/chr{row['chromosome']}_{row['contig']}"
        bcf = require(folder / BCF_NAME)
        require(str(bcf) + '.csi')
        audit = json.loads(require(folder / 'call.complete.json').read_text())
        if audit['samples'] != config['samples'] or audit['bcf_bytes'] != bcf.stat().st_size:
            raise RuntimeError(f'Changed chromosome output: {bcf}')
        bcf_files.append(bcf)
        audits.append(audit)
    final = args.root / 'variant_calls/genomewide'
    if final.exists():
        raise RuntimeError(f'Final output already exists; preserved: {final}')
    with lock(args.root / '.finalize.lock'):
        with tempfile.TemporaryDirectory(prefix='.finalize.', dir=args.root) as tmp:
            output = Path(tmp) / 'publish'
            output.mkdir()
            filelist = Path(tmp) / 'chromosomes.list'
            filelist.write_text(''.join(str(p) + '\n' for p in bcf_files))
            bcf = output / BCF_NAME
            vcf = output / 'NDJ_controls.joint_at_control_sites.vcf.gz'
            run([args.bcftools, 'concat', '-f', filelist, '-Ob', '--threads', '3', '-o', bcf])
            run([args.bcftools, 'index', '-f', bcf])
            run([args.bcftools, 'view', '-Oz', '--threads', '3', '-o', vcf, bcf])
            run([args.bcftools, 'index', '-f', '-t', vcf])
            expected = sum(a['emitted_sites'] for a in audits)
            for path in (bcf, vcf):
                observed_samples = run([args.bcftools, 'query', '-l', path], capture_output=True, text=True).stdout.splitlines()
                count = int(run([args.bcftools, 'index', '-n', path], capture_output=True, text=True).stdout.strip())
                if observed_samples != config['samples'] or count != expected:
                    raise RuntimeError(f'Final count/sample mismatch: {path}')
            dump(output / 'genomewide.complete.json', {'samples': config['samples'], 'chromosomes': 21,
                 'emitted_sites': expected, 'target_sites': sum(a['target_sites'] for a in audits),
                 'omitted_target_sites': sum(a['omitted_target_sites'] for a in audits), 'per_chromosome': audits,
                 'ploidy_note': 'Default diploid calling, as in earlier workflow. Z GTs are not sex-specific; AF uses AD.'})
            output.rename(final)
    print(f'GENOMEWIDE_COMPLETE\tsamples=20\tsites={expected}\tVCF={final / vcf.name}', flush=True)


def plot(args):
    final = args.root / 'variant_calls/genomewide'
    require(final / 'genomewide.complete.json')
    run([sys.executable, args.project / 'pairwise-af/plot_pairwise_AF.py', '--vcf', final / BCF_NAME,
         '--left-sample', NAMES[0], '--right-sample', NAMES[1], '--min-depth', '10', '--bins', '81',
         '--title', 'NDJ Left versus Right at control-discovered SNPs',
         '--output-prefix', args.root / 'figures/NDJ_control_sites'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['preflight', 'prepare', 'call', 'finalize', 'plot'])
    parser.add_argument('--project', type=Path, default=Path('/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling'))
    parser.add_argument('--bcftools', type=Path, default=Path('/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin/bcftools'))
    parser.add_argument('--samtools', type=Path, default=Path('/gpfs/automountdir/gpfs/homes/CCAS/home/jjhanly/.local/share/mamba/envs/colias_w_popgen/bin/samtools'))
    args = parser.parse_args()
    args.project = args.project.resolve()
    args.root = args.project / 'analysis/NDJ_control_comparison_v1'
    try:
        {'preflight': prepare, 'prepare': prepare, 'call': call_chromosome,
         'finalize': finalize, 'plot': plot}[args.stage](args)
    except Exception as error:
        print(f'ERROR\t{args.stage}\t{error}', file=sys.stderr, flush=True)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
