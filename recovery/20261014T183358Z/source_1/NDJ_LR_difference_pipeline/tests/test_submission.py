import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ndj_scan.cli import submit_chain


def test_submission_chain_uses_afterok_and_arrays(tmp_path, monkeypatch):
    run = tmp_path / "run"
    (run / "logs").mkdir(parents=True)
    (run / "code").mkdir()
    (run / "code" / "run_ndj_scan.py").write_text("print('worker')\n")
    calls = []

    class Result:
        returncode = 0
        stderr = ""
        stdout = ""

    def fake_run(command, **kwargs):
        calls.append(command)
        result = Result()
        if command[0] == "sbatch":
            result.stdout = str(9000 + sum(call[0] == "sbatch" for call in calls)) + "\n"
        return result

    monkeypatch.setattr("ndj_scan.cli.subprocess.run", fake_run)
    args = type("Args", (), {"parallel": 7})()
    config = {"chromosomes": [{"contig": "a"}, {"contig": "b"}]}
    submit_chain(run, config, args)
    submitted = [call for call in calls if call[0] == "sbatch"]
    assert len(submitted) == 6
    assert "--array=0-1%7" in submitted[0]
    assert "--dependency=afterok:9001" in submitted[1]
    assert "--dependency=afterok:9005" in submitted[5]
    rows = [json.loads(line) for line in (run / "submissions.jsonl").read_text().splitlines()]
    assert [row["stage"] for row in rows] == ["scan", "background", "score", "shortlist", "qc", "report"]
