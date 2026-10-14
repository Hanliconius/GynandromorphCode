import gzip
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pysam
import pytest

BASE = Path(__file__).resolve().parents[1]


def module(filename):
    spec = importlib.util.spec_from_file_location(filename.split('.')[0], BASE / 'scripts' / filename)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


driver = module('218_NDJ_controls.py')
submitter = module('218_submit_NDJ_finish.py')
SAMPLES = ['control'] + driver.NAMES


def alleles(tmp_path, text='chr1\t10\tA,G\nchr1\t20\tC,T\nchr1\t30\tG,A\n'):
    path = tmp_path / 'alleles.tsv.gz'
    with gzip.open(path, 'wt') as f:
        f.write(text)
    return path


def bcf(tmp_path, rows=((10, 'A', 'G'), (30, 'G', 'A')), samples=SAMPLES):
    header = pysam.VariantHeader()
    header.contigs.add('chr1', length=100)
    header.formats.add('GT', 1, 'String', 'Genotype')
    header.formats.add('AD', 'R', 'Integer', 'Allelic depth')
    for sample in samples:
        header.add_sample(sample)
    path = tmp_path / 'calls.bcf'
    with pysam.VariantFile(str(path), 'wb', header=header) as output:
        for position, ref, alt in rows:
            r = output.new_record(contig='chr1', start=position - 1, alleles=(ref, alt))
            for sample in samples:
                r.samples[sample]['GT'] = (0, 1)
                r.samples[sample]['AD'] = (10, 10)
            output.write(r)
    return path


def test_missing_targets_reported_not_failed(tmp_path):
    audit = driver.validate_bcf(bcf(tmp_path), alleles(tmp_path), 'chr1', SAMPLES)
    assert audit['target_sites'] == 3
    assert audit['emitted_sites'] == 2
    assert audit['omitted_target_sites'] == 1


def test_missing_at_both_ends(tmp_path):
    audit = driver.validate_bcf(bcf(tmp_path, [(20, 'C', 'T')]), alleles(tmp_path), 'chr1', SAMPLES)
    assert audit['target_sites'] == 3
    assert audit['omitted_target_sites'] == 2


@pytest.mark.parametrize('rows', [[(10, 'A', 'T')], [(11, 'A', 'G')], [(10, 'A', 'G'), (10, 'A', 'G')]])
def test_wrong_or_duplicate_emitted_alleles_fail(tmp_path, rows):
    with pytest.raises(RuntimeError):
        driver.validate_bcf(bcf(tmp_path, rows), alleles(tmp_path), 'chr1', SAMPLES)


def test_sample_order_checked(tmp_path):
    with pytest.raises(RuntimeError, match='sample names/order'):
        driver.validate_bcf(bcf(tmp_path, samples=list(reversed(SAMPLES))), alleles(tmp_path), 'chr1', SAMPLES)


@pytest.mark.parametrize('text', ['chr1\t10\tA\tG\n', 'chr1\t10\tA,G\nchr1\t10\tA,G\n',
                                  'chr1\t20\tA,G\nchr1\t10\tA,G\n', 'wrong\t10\tA,G\n'])
def test_wrong_target_formats_fail(tmp_path, text):
    with pytest.raises(RuntimeError):
        list(driver.targets(alleles(tmp_path, text), 'chr1'))


@pytest.mark.parametrize('producer_rc,consumer_rc', [(1, 0), (0, 1)])
def test_pipeline_checks_both_exit_codes(producer_rc, consumer_rc):
    first = [sys.executable, '-c', f'import sys; print("data"); sys.exit({producer_rc})']
    second = [sys.executable, '-c', f'import sys; sys.stdin.read(); sys.exit({consumer_rc})']
    with pytest.raises(RuntimeError, match='Calling pipeline failed'):
        driver.pipeline(first, second)


def test_pipeline_passes():
    driver.pipeline([sys.executable, '-c', 'print("data")'],
                    [sys.executable, '-c', 'import sys; assert sys.stdin.read().strip() == "data"'])


def test_lock_not_removed_if_not_owned(tmp_path):
    path = tmp_path / 'busy'
    path.mkdir()
    with pytest.raises(FileExistsError):
        with driver.lock(path):
            pass
    assert path.is_dir()


def test_changed_inputs_fail(tmp_path):
    path = tmp_path / 'input'
    path.write_text('abc')
    recorded = driver.stamp(path)
    path.write_text('abcdef')
    with pytest.raises(RuntimeError, match='Input changed'):
        driver.unchanged(recorded)


def setup_submission(tmp_path):
    python = tmp_path / 'pairwise-af/.venv-af/bin/python'
    python.parent.mkdir(parents=True)
    python.touch()
    script = tmp_path / 'scripts/218_NDJ_controls.py'
    script.parent.mkdir()
    script.touch()


def test_submission_whole_array_dependencies(tmp_path, monkeypatch):
    setup_submission(tmp_path)
    calls = []
    def fake_run(cmd, **kw):
        calls.append(cmd)
        if cmd[0] == 'sbatch':
            return SimpleNamespace(stdout=str(100 + sum(c[0] == 'sbatch' for c in calls)) + ';cluster\n', stderr='', returncode=0)
        return SimpleNamespace(stdout='', stderr='', returncode=0)
    monkeypatch.setattr(submitter.subprocess, 'run', fake_run)
    monkeypatch.setattr(sys, 'argv', ['submit', '--alignment-job', '138290948', '--project', str(tmp_path)])
    submitter.main()
    batch = [c for c in calls if c[0] == 'sbatch']
    assert [c[c.index('--dependency') + 1] for c in batch] == ['afterok:138290948', 'afterok:101', 'afterok:102', 'afterok:103']
    assert batch[1][batch[1].index('--array') + 1] == '0-20%21'
    assert batch[0][batch[0].index('--partition') + 1] == 'nano'
    assert len((tmp_path / 'analysis/NDJ_control_comparison_v1/finish_chain.submissions.jsonl').read_text().splitlines()) == 4
    with pytest.raises(FileExistsError):
        submitter.main()
    assert len([c for c in calls if c[0] == 'sbatch']) == 4


def test_partial_submission_stops_and_records_error(tmp_path, monkeypatch):
    setup_submission(tmp_path)
    count = []
    def fake_run(cmd, **kw):
        if cmd[0] == 'sbatch':
            count.append(cmd)
            return SimpleNamespace(stdout='101\n' if len(count) == 1 else '',
                                   stderr='' if len(count) == 1 else 'test rejection', returncode=0 if len(count) == 1 else 1)
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(submitter.subprocess, 'run', fake_run)
    monkeypatch.setattr(sys, 'argv', ['submit', '--alignment-job', '138290948', '--project', str(tmp_path)])
    with pytest.raises(RuntimeError, match='test rejection'):
        submitter.main()
    assert len(count) == 2
    records = [json.loads(line) for line in (tmp_path / 'analysis/NDJ_control_comparison_v1/finish_chain.submissions.jsonl').read_text().splitlines()]
    assert records[0]['submission_response'] == '101'
    assert records[1]['returncode'] == 1


def test_real_bcftools_call_20_samples(tmp_path):
    """Exercise actual bundled bcftools commands, using files instead of a shell pipe."""
    import pysam.bcftools as bc
    reference = tmp_path / 'reference.fa'
    reference.write_text('>chr1\n' + 'A' * 200 + '\n')
    pysam.faidx(str(reference))
    samples = [f'control_{i:02d}' for i in range(18)] + driver.NAMES
    bams = []
    for index, name in enumerate(samples):
        path = tmp_path / (name + '.bam')
        header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
                  'SQ': [{'SN': 'chr1', 'LN': 200}], 'RG': [{'ID': name, 'SM': name}]}
        with pysam.AlignmentFile(str(path), 'wb', header=header) as handle:
            for n in range(20):
                read = pysam.AlignedSegment()
                read.query_name = f'{name}_{n}'
                sequence = list('A' * 90)
                if n < 10 or index == 19:
                    sequence[49] = 'G'  # Reference POS=60, start0=10.
                read.query_sequence = ''.join(sequence)
                read.flag = 0
                read.reference_id = 0
                read.reference_start = 10
                read.mapping_quality = 60
                read.cigar = [(0, 90)]
                read.query_qualities = pysam.qualitystring_to_array('I' * 90)
                read.set_tag('RG', name)
                handle.write(read)
        pysam.index(str(path))
        bams.append(path)
    bamlist = tmp_path / 'bams.list'
    bamlist.write_text(''.join(str(p) + '\n' for p in bams))
    raw_targets = tmp_path / 'targets.tsv'
    raw_targets.write_text('chr1\t60\tA,G\nchr1\t120\tA,G\n')
    compressed = tmp_path / 'targets.tsv.gz'
    pysam.tabix_compress(str(raw_targets), str(compressed))
    pysam.tabix_index(str(compressed), seq_col=0, start_col=1, end_col=1)
    raw, output = tmp_path / 'pileup.bcf', tmp_path / 'calls.bcf'
    bc.mpileup('-Ou', '--threads', '1', '-f', str(reference), '-r', 'chr1', '-T', str(compressed),
               '-q', '20', '-Q', '20', '-d', '1000', '-a', 'FORMAT/AD,FORMAT/DP,FORMAT/SP',
               '-b', str(bamlist), '-o', str(raw), catch_stdout=False)
    bc.call('-Ob', '--threads', '1', '-m', '-A', '-C', 'alleles', '-T', str(compressed),
            '-G', '-', '-o', str(output), str(raw), catch_stdout=False)
    bc.index('-f', str(output), catch_stdout=False)
    audit = driver.validate_bcf(output, compressed, 'chr1', samples)
    assert audit['emitted_sites'] == 1
    assert audit['omitted_target_sites'] == 1
    assert audit['NDJ_nonzero_AD_sites'] == {'NDJ_left': 1, 'NDJ_right': 1}
    assert int(bc.index('-n', str(output)).strip()) == 1
    compressed_vcf = tmp_path / 'calls.vcf.gz'
    bc.view('-Oz', '-o', str(compressed_vcf), str(output), catch_stdout=False)
    bc.index('-f', '-t', str(compressed_vcf), catch_stdout=False)
    assert int(bc.index('-n', str(compressed_vcf)).strip()) == 1
