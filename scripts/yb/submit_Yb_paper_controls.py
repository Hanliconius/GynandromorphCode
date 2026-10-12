#!/usr/bin/env python3
"""Align ten accession-verified local paper controls to both prepared references."""
import argparse,csv,json,shlex,shutil,subprocess,sys,uuid
from pathlib import Path
from datetime import datetime,timezone
ROOT=Path("/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling")
CONTROLS=[
("ERR2298257","present","melpomene","amaryllis"),
("ERR2298260","present","melpomene","amaryllis"),
("ERR2298208","present","melpomene","rosina"),
("ERR2298221","present","melpomene","rosina"),
("ERR2298263","present","timareta","thelxinoe"),
("ERR2298265","present","timareta","thelxinoe"),
("ERR2298201","absent","melpomene","melpomene"),
("ERR2298216","absent","melpomene","melpomene"),
("ERR2298235","absent","timareta","florencia"),
("ERR2298236","absent","timareta","florencia")]
def require(p):
    p=Path(p)
    if not p.is_absolute() or not p.is_file() or not p.stat().st_size:
        raise RuntimeError(f"Missing nonempty absolute input: {p}")
    return str(p)
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--focal-run",type=Path,required=True)
    a=p.parse_args()
    parent=a.focal_run
    require(parent/"run.json")
    old=json.loads((parent/"run.json").read_text())
    executor=Path(require(parent/"rebuild_Yb_local_both_references.py"))
    samples=[]
    for accession,phenotype,species,subspecies in CONTROLS:
        directory=ROOT/"analysis/genomewide_mixing_input/fastq"/accession
        samples.append(dict(sample=accession,library_id=accession,
            r1=require(directory/(accession+"_1.fastq.gz")),
            r2=require(directory/(accession+"_2.fastq.gz")),
            yellow_bar=phenotype,species=species,subspecies=subspecies,
            source="Livraghi2021_fig8_data1",side="control",leg="control"))
    for name,path in old["tools"].items():
        require(path)
    for ref in old["references"]:
        fasta=parent/ref["name"]/"reference/genome.fa"
        for suffix in ("",".fai",".amb",".ann",".bwt",".pac",".sa"):
            require(str(fasta)+suffix)
    runid=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"_"+uuid.uuid4().hex[:8]
    base=ROOT/"analysis/yb_control_both_references"/runid
    base.mkdir(parents=True)
    (base/"logs").mkdir()
    shutil.copy2(executor,base/executor.name)
    for ref in old["references"]:
        (base/ref["name"]).mkdir()
        (base/ref["name"]/"reference").symlink_to(parent/ref["name"]/"reference",target_is_directory=True)
    config={**old,"run_dir":str(base),"code_dir":str(base),
        "tasks":[dict(reference=ref,**sample) for ref in old["references"] for sample in samples],
        "focal_run":str(parent),"manifest_source":"explicit accession-named local FASTQ pairs",
        "scope":"control alignment and raw coverage only; normalized overlays follow",
        "excluded_alias":"ERR1307111 is not substituted for ERR2298208"}
    configpath=base/"run.json"
    configpath.write_text(json.dumps(config,indent=2)+"\n")
    with (base/"controls.tsv").open("w") as h:
        w=csv.DictWriter(h,fieldnames=list(samples[0]),delimiter="\t")
        w.writeheader();w.writerows(samples)
    worker=base/"worker.sh"
    worker.write_text("#!/usr/bin/env bash\nset -euo pipefail\nexec "+
        shlex.join([sys.executable,str(base/executor.name),"--config",str(configpath),"--stage","align"])+
        ' --task "$SLURM_ARRAY_TASK_ID"\n')
    command=["sbatch","--parsable","--partition=cpu","--time=2-00:00:00",
        "--cpus-per-task=16","--mem=64G","--array=0-19%8","--job-name=YbControls",
        "--export=ALL",f"--chdir={base}",f"--output={base}/logs/align.%A_%a.out",
        f"--error={base}/logs/align.%A_%a.err",str(worker)]
    job=subprocess.check_output(command,text=True).strip().split(";")[0]
    (base/"jobs.json").write_text(json.dumps({"align":job})+"\n")
    (base.parent/"latest_run.path").write_text(str(base)+"\n")
    print(f"ALIGN_ARRAY\t{job}")
    print(f"RUN_DIR\t{base}")
    print("CONTROLS\t6 Yb+; 4 Yb-; 20 alignment tasks, maximum 8 concurrent")
    print("INPUTS\tlocal accession-named pairs; existing indexes; no downloads")
if __name__=="__main__":
    try:main()
    except Exception as e:
        print(f"ERROR: {e}",file=sys.stderr);sys.exit(2)
