import importlib.util
import json
from pathlib import Path
import shutil
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('regen', ROOT / 'scripts/219_regenerate_TN_NDJ_AF.py')
regen = importlib.util.module_from_spec(spec)
spec.loader.exec_module(regen)


def prepare_project(tmp_path):
    project = tmp_path / 'project with spaces'
    package = project / 'pairwise-af'
    package.mkdir(parents=True)
    shutil.copyfile(ROOT / 'pairwise-af/plot_pairwise_AF.py', package / 'plot_pairwise_AF.py')
    vcf = project / 'input.vcf'
    shutil.copyfile(ROOT / 'pairwise-af/examples/tiny.vcf', vcf)
    tool = regen.load_tool(project)
    for _, path in regen.AUDITS:
        prefix = project / path.removesuffix('.AF.audit.json')
        assert tool.main(['--vcf', str(vcf), '--left-sample', 'Left', '--right-sample', 'Right',
                          '--output-prefix', str(prefix), '--bin-shape', 'square', '--bins', '81',
                          '--filter', 'pass', '--contig', 'chr1', '--mixing-left', '.2',
                          '--mixing-right', '.2']) == 0
    return project, tool


def test_full_regeneration_preserves_selection_and_shared_scale(tmp_path):
    project, tool = prepare_project(tmp_path)
    plan = regen.make_plan(project, tool)
    assert plan['vmin'] == 1
    assert plan['vmax'] >= max(e['retained'] for e in plan['entries'])
    path = project / 'plan.json'
    path.write_text(json.dumps(plan))
    for index, entry in enumerate(plan['entries']):
        args = tool.parser().parse_args(entry['argv'])
        assert args.bin_shape == 'hex' and args.bins == 65
        assert args.filter == 'pass' and args.contig == ['chr1']
        assert args.mixing_left == .2 and args.mixing_right == .2
        assert args.vmax == plan['vmax']
        regen.run_task(path, index)
        assert Path(entry['output_prefix'] + '.AF.pdf').is_file()
    with pytest.raises(RuntimeError, match='already exists'):
        regen.make_plan(project, tool)


def test_changed_source_refused_before_jobs(tmp_path):
    project, tool = prepare_project(tmp_path)
    with (project / 'input.vcf').open('a') as handle:
        handle.write('\n')
    with pytest.raises(RuntimeError, match='Input changed'):
        regen.make_plan(project, tool)
