from pathlib import Path
import importlib.util
import json

project = Path.cwd()
tool_path = project / "pairwise-af-hex20/plot_pairwise_AF.py"
spec = importlib.util.spec_from_file_location("af_tool", tool_path)
tool = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tool)

root = project / "analysis/figures/pairwise_AF_cli"
old_path = root / "TN_merged_control_sites_hex20_batlowW.AF.audit.json"
old = json.loads(old_path.read_text())
params = old["parameters"]
prefix = root / "TN_merged_control_sites_hex40_batlowW"
new_path = Path(str(prefix) + ".AF.audit.json")

if tool.fingerprint(Path(old["input"]["path"])) != old["input"]:
    raise RuntimeError("Input VCF metadata changed since the original plot.")

argv = [
    "--vcf", old["input"]["path"],
    "--output-prefix", str(prefix),
    "--bins", "40",
    "--bin-shape", "hex",
    "--cmap", "batlowW_r",
    "--scale", "count",
]
for key in (
    "left_sample", "right_sample", "min_depth",
    "max_depth", "filter", "title", "vmin", "vmax",
):
    value = params.get(key)
    if value is not None:
        argv.extend(["--" + key.replace("_", "-"), str(value)])

for contig in params.get("contig") or []:
    argv.extend(["--contig", contig])

if params.get("allow_self"):
    argv.append("--allow-self")

if params.get("sites"):
    whitelist = old.get("whitelist")
    if not whitelist:
        raise RuntimeError("Missing original site-whitelist provenance.")
    sites = Path(whitelist["path"])
    if tool.sha256(sites) != whitelist["sha256"]:
        raise RuntimeError("Site whitelist changed.")
    argv.extend(["--sites", str(sites)])

tool.check_args(tool.parser().parse_args(argv))
if tool.main(argv):
    raise RuntimeError("AF plotting failed.")

new = json.loads(new_path.read_text())
for key in (
    "records_read", "records_retained",
    "excluded_first_reason", "retained_by_contig",
):
    if old[key] != new[key]:
        raise RuntimeError(f"Site selection differs: {key}")

print(f"SELECTION_MATCH\tretained={new['records_retained']:,}")
print(str(prefix) + ".AF.pdf")
