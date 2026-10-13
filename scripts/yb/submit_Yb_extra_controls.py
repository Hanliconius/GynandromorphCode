#!/usr/bin/env python3
"""Download nine extra Yb controls, then align each to both prepared references."""
import argparse,csv,gzip,hashlib,io,json,shlex,shutil,subprocess,sys,uuid
from pathlib import Path
from datetime import datetime,timezone
from urllib.request import urlopen
from urllib.parse import urlencode
ROOT=Path("/gpfs/automountdir/gpfs/scratch/martinlab/joe/GynandroCalling")
SAMPLES=[(s,"BWPK","present","user-confirmed; region_panel_metadata dorsal=1 ventral=1") for s in
    ("SRR11243344","SRR11243353","SRR11243364","SRR11243375","SRR11243376")]
SAMPLES += [(s,"rosina","present","Livraghi2021 Figure8 source data1") for s in ("ERR2298205","ERR2298206")]
SAMPLES += [(s,"melpomene","absent","Livraghi2021 Figure8 source data1") for s in ("ERR2298204","ERR2298214")]
def require(p):
    p=Path(p)
    if not p.is_absolute() or not p.is_file() or not p.stat().st_size:raise RuntimeError("Missing "+str(p))
    return p
def digest(path):
    h=hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b""):h.update(chunk)
    return h.hexdigest()
def download(cfg,task):
    e=cfg["samples"][task];accession=e["sample"];directory=Path(cfg["fastq_root"])/accession
    directory.mkdir(parents=True,exist_ok=True)
    query=urlencode(dict(accession=accession,result="read_run",format="tsv",
        fields="run_accession,sample_accession,library_layout,fastq_ftp,fastq_md5,fastq_bytes"))
    with urlopen("https://www.ebi.ac.uk/ena/portal/api/filereport?"+query,timeout=120) as r:
        rows=list(csv.DictReader(io.StringIO(r.read().decode()),delimiter="\t"))
    rows=[r for r in rows if r.get("run_accession")==accession]
    if len(rows)!=1 or rows[0]["library_layout"]!="PAIRED":raise RuntimeError("ENA did not resolve one paired run: "+accession)
    row=rows[0];urls=row["fastq_ftp"].split(";");md5=row["fastq_md5"].split(";");sizes=row["fastq_bytes"].split(";")
    if not len(urls)==len(md5)==len(sizes):raise RuntimeError("ENA checksum/URL mismatch")
    resolved=[]
    for mate in (1,2):
        hits=[i for i,u in enumerate(urls) if u.endswith(f"_{mate}.fastq.gz") or u.endswith(f"_{mate}.fq.gz")]
        if len(hits)!=1:raise RuntimeError("No unique ENA paired FASTQ for "+accession+"; review SRA fallback")
        i=hits[0];url=urls[i];target=directory/f"{accession}_{mate}.fastq.gz"
        valid=target.is_file() and target.stat().st_size==int(sizes[i]) and digest(target)==md5[i]
        if not valid:
            partial=target.with_name(target.name+".part")
            if not url.startswith("https://"):url="https://"+url.removeprefix("ftp://")
            subprocess.run(["curl","--fail","--location","--retry","4","--retry-delay","5",
                "--silent","--show-error","--output",str(partial),url],check=True)
            if partial.stat().st_size!=int(sizes[i]) or digest(partial)!=md5[i]:
                partial.unlink();raise RuntimeError("Checksum failed: "+accession)
            with gzip.open(partial,"rb") as f:
                for chunk in iter(lambda:f.read(8*1024*1024),b""):pass
            partial.replace(target)
        resolved.append(dict(mate=mate,path=str(target),ENA_url=urls[i],md5=md5[i],bytes=int(sizes[i])))
    (directory/"download.complete.json").write_text(json.dumps({"sample":e,"ENA":row,"files":resolved},indent=2)+"\n")
    print("FASTQ_READY\t"+accession)
def submit(focal):
    old=json.loads(require(focal/"run.json").read_text())
    require(focal/"rebuild_Yb_local_both_references.py")
    if not shutil.which("curl"):raise RuntimeError("curl unavailable")
    for ref in old["references"]:
        for ext in ("",".fai",".amb",".ann",".bwt",".pac",".sa"):
            require(str(focal/ref["name"]/"reference/genome.fa")+ext)
    base=ROOT/"analysis/yb_extra_controls"/(datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")+"_"+uuid.uuid4().hex[:8])
    (base/"logs").mkdir(parents=True)
    fastq=ROOT/"analysis/yb_extra_control_fastq"
    samples=[dict(sample=s,library_id=s,subspecies=race,species="unknown" if race=="BWPK" else "melpomene",
        yellow_bar=yb,phenotype_source=source,r1=str(fastq/s/f"{s}_1.fastq.gz"),r2=str(fastq/s/f"{s}_2.fastq.gz"),
        role="control",side="control",leg="control") for s,race,yb,source in SAMPLES]
    code=base/Path(__file__).name;shutil.copy2(__file__,code)
    aligner=base/"rebuild_Yb_local_both_references.py";shutil.copy2(focal/aligner.name,aligner)
    for ref in old["references"]:
        (base/ref["name"]).mkdir()
        (base/ref["name"]/"reference").symlink_to(focal/ref["name"]/"reference",target_is_directory=True)
    cfg={**old,"run_dir":str(base),"code_dir":str(base),"samples":samples,"fastq_root":str(fastq),
         "focal_run":str(focal),"tasks":[dict(reference=ref,**e) for ref in old["references"] for e in samples]}
    config=base/"run.json";config.write_text(json.dumps(cfg,indent=2)+"\n")
    worker=base/"worker.sh"
    worker.write_text("#!/usr/bin/env bash\nset -euo pipefail\n"+
        'if [[ "$1" == download ]]; then\n exec '+shlex.join([sys.executable,str(code),"--config",str(config)])+' --task "$SLURM_ARRAY_TASK_ID"\n'+
        'else\n exec '+shlex.join([sys.executable,str(aligner),"--config",str(config),"--stage","align"])+' --task "$SLURM_ARRAY_TASK_ID"\nfi\n')
    jobs={}
    try:
        for stage,array,cpus,mem,limit in [("download","0-8%9",1,"4G","12:00:00"),("align","0-17%18",16,"64G","2-00:00:00")]:
            cmd=["sbatch","--parsable","--partition=cpu",f"--job-name=YbExtra_{stage}",f"--array={array}",
                f"--cpus-per-task={cpus}",f"--mem={mem}",f"--time={limit}",f"--chdir={base}",
                f"--output={base}/logs/{stage}.%A_%a.out",f"--error={base}/logs/{stage}.%A_%a.err"]
            if stage=="align":cmd.append("--dependency=afterok:"+jobs["download"])
            jobs[stage]=subprocess.check_output(cmd+[str(worker),stage],text=True).strip().split(";")[0]
            (base/"jobs.json").write_text(json.dumps(jobs,indent=2)+"\n");print(stage.upper()+"_ARRAY\t"+jobs[stage],flush=True)
    except Exception:
        for j in jobs.values():subprocess.run(["scancel",j],check=False)
        raise
    print("RUN_DIR\t"+str(base))
    print("CONTROLS\t5 BWPK Yb+; 2 extra rosina Yb+; 2 extra melpomene Yb-; 18 alignments")
def main():
    p=argparse.ArgumentParser();p.add_argument("--focal-run",type=Path);p.add_argument("--config",type=Path);p.add_argument("--task",type=int);a=p.parse_args()
    if a.focal_run:submit(a.focal_run)
    else:download(json.loads(a.config.read_text()),a.task)
if __name__=="__main__":
    try:main()
    except Exception as e:print("ERROR: "+str(e),file=sys.stderr);sys.exit(2)
