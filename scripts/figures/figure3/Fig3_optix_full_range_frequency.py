#!/usr/bin/env python3
"""Full-range optix frequency genoplot; comparative selection independent of TN."""
import argparse
import gzip
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import Fig3_optix_BA099_both_callable as base

CHROM, START, END = "Hmr18", 1_060_000, 1_240_000
KEY = base.KEY
MODULES = ("Band", "Dennis", "Ray")
plt.rcParams.update({"pdf.fonttype":42, "ps.fonttype":42, "font.size":8})

def cohorts():
    manifest = base.load_manifest()
    band_table = pd.read_csv(base.source_path("Band"), sep="\t",
                             usecols=["sample", "plot_group"])
    band_table["phenotype"] = band_table["plot_group"].map(
        {"Band absent (B=0)":0, "Band present (B=1)":1})
    status = band_table.dropna(subset=["phenotype"])[["sample","phenotype"]].drop_duplicates()
    if status.groupby("sample")["phenotype"].nunique().gt(1).any():
        raise RuntimeError("Conflicting Band phenotype assignments")
    answer = {}
    for module in MODULES:
        table = manifest.loc[manifest.BA_pattern.eq(module)].copy()
        if module == "Band":
            table = table.drop(columns="phenotype").merge(status, on="sample",
                         how="left", validate="many_to_one")
        # The focal specimen never contributes to comparative frequencies or selection.
        table = table.loc[~table["sample"].str.match(r"^TN[_-]?(Left|Right)", case=False)]
        table = table.loc[table.phenotype.isin([0,1]), ["sample","phenotype"]].drop_duplicates()
        if table.groupby("sample").phenotype.nunique().gt(1).any():
            raise RuntimeError(f"{module}: conflicting comparative phenotype")
        table = table.drop_duplicates("sample")
        if set(table.phenotype) != {0,1}:
            raise RuntimeError(f"{module}: both comparative groups required")
        answer[module] = table
    return answer

def read_calls(path, selected, samples):
    wanted = {tuple(row) for row in selected[KEY].itertuples(index=False, name=None)}
    calls, seen, indices = {}, set(), {}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("##"):continue
            if line.startswith("#CHROM"):
                header = line.rstrip().split("\t")
                lookup = {name:i for i,name in enumerate(header)}
                missing = sorted(set(samples)-set(lookup))
                if missing:raise RuntimeError("Comparative samples absent from VCF: "+", ".join(missing))
                indices = {name:lookup[name] for name in samples}
                continue
            if line.startswith("#"):continue
            if not indices:raise RuntimeError("VCF sample header not found")
            fields = line.rstrip().split("\t")
            key = (fields[0], int(fields[1]), fields[3], fields[4])
            if key not in wanted:continue
            if key in seen:raise RuntimeError(f"Duplicate selected VCF variant: {key}")
            seen.add(key)
            fmt = {name:i for i,name in enumerate(fields[8].split(":"))}
            if not {"GT","DP","GQ"} <= set(fmt):
                raise RuntimeError(f"Selected VCF site lacks GT/DP/GQ: {key}; FORMAT={fields[8]}")
            for sample,index in indices.items():
                value = fields[index].split(":")
                dosage = np.nan
                try:
                    if float(value[fmt["DP"]]) >= 2 and float(value[fmt["GQ"]]) >= 20:
                        dosage = base.alt_dosage(value[fmt["GT"]])
                except (ValueError, IndexError):pass
                calls[(key,sample)] = dosage
    return calls, seen

def build_values(selected, recalled, groups, calls):
    records, site_records = [], []
    for module in MODULES:
        cohort = groups[module]
        sites = selected.loc[selected.module.eq(module)].merge(recalled, on=KEY,
                   how="left", validate="one_to_one").copy()
        for _,site in sites.iterrows():
            key = tuple(site[k] for k in KEY)
            group_values = {}
            for phenotype in (0,1):
                members = cohort.loc[cohort.phenotype.eq(phenotype), "sample"].tolist()
                values = np.array([calls.get((key,sample),np.nan) for sample in members])
                if not bool(site.pattern_is_alt):values = 1-values
                finite = np.isfinite(values)
                group_values[phenotype] = (float(np.mean(values[finite])) if finite.any() else np.nan,
                                          int(finite.sum()),len(members))
            # Comparative callability applies before any focal information is displayed.
            eligible = all(n/total >= .8 for _,n,total in group_values.values())
            site_records.append({**{k:site[k] for k in selected.columns},
                                 "call_rate_absent":group_values[0][1]/group_values[0][2],
                                 "call_rate_present":group_values[1][1]/group_values[1][2],
                                 "comparative_callable":eligible})
            if not eligible:continue
            common = {k:site[k] for k in KEY}
            common.update(module=module, balanced_accuracy=site.balanced_accuracy,
                          selected_orientation=site.selected_orientation)
            for phenotype,(freq,n,total) in group_values.items():
                records.append({**common,"row":("present" if phenotype else "absent"),
                                "frequency":freq,"called":n,"total":total,
                                "pair_posterior":np.nan,"source":"comparative DP2 GQ20 diploid allele frequency"})
            confidence = pd.notna(site.get("both_recalled_callable")) and bool(site.get("both_recalled_callable"))
            for side in ("TN_Left","TN_Right"):
                freq = base.pattern_dosage(site.get(side+"_GT"),site.pattern_is_alt) if confidence else np.nan
                records.append({**common,"row":side,"frequency":freq,"called":int(np.isfinite(freq)),
                                "total":1,"pair_posterior":site.get("pair_posterior"),
                                "source":"reciprocal-mixture posterior>=0.95; no original GT fallback"})
    return pd.DataFrame(records),pd.DataFrame(site_records)

def plot(values,out):
    cmap = LinearSegmentedColormap.from_list("pattern_frequency",
           [(0,"#D7EDF5"),(.5,"#ECCB63"),(1,"#662A83")])
    cmap.set_bad("#C8C8C8")
    fig,axes = plt.subplots(3,1,figsize=(12,7.5),sharex=True)
    fig.subplots_adjust(left=.19,right=.88,top=.91,bottom=.13,hspace=.48)
    for axis,module in zip(axes,MODULES):
        rows = (["absent","TN_Left","TN_Right","present"] if module=="Band"
                else ["present","TN_Left","TN_Right","absent"])
        table=values.loc[values.module.eq(module)]
        count=table[KEY].drop_duplicates().shape[0]
        for y,row in enumerate(rows):
            subset=table.loc[table.row.eq(row)]
            axis.axhline(y,color="#eeeeee",lw=.6,zorder=0)
            axis.scatter(subset.pos/1e6,np.full(len(subset),y),c=subset.frequency.to_numpy(),
                         cmap=cmap,vmin=0,vmax=1,marker="|",s=180,linewidths=1.2,plotnonfinite=True)
        axis.set_yticks(range(4),[f"{module} {r}" if r in ("present","absent")
                                  else r+" adjusted" for r in rows])
        axis.set_ylim(3.5,-.5)
        axis.set_title(f"{module}: {count} comparative-selected SNPs",loc="left",fontsize=10)
        axis.set_xlim(START/1e6,END/1e6)
        for spine in axis.spines.values():spine.set_visible(False)
        axis.tick_params(axis="y",length=0)
    axes[-1].set_xlabel("Position on Hmr18 (Mb); SNP ticks at actual genomic positions")
    bar=fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0,1),cmap=cmap),
                    cax=fig.add_axes([.915,.24,.013,.5]))
    bar.set_ticks([0,.5,1]);bar.set_label("Pattern-associated allele frequency / focal dosage")
    fig.suptitle("optix: full-range B / Dennis / Ray frequency genoplot",x=.19,ha="left",fontsize=13)
    fig.text(.19,.045,"Selection: comparative BA > 0.99; both phenotype groups >=80% callable (DP>=2, GQ>=20).\n"
             "Focal display: adjusted reciprocal-mixture posterior >=0.95; grey = unavailable. Focal state does not select sites.",
             fontsize=8,ha="left")
    pdf=out/"optix_full_range_BDR_frequency_adjusted_LR.pdf"
    fig.savefig(pdf);fig.savefig(pdf.with_suffix(".png"),dpi=250);plt.close(fig)
    return pdf

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--vcf",type=Path,default=base.PROJECT /
       "external/krista_genotype_plots_original/Genotype_Plots/Optix_Hmr18-938195-1239421_biallelic.vcf")
    parser.add_argument("--output-dir",type=Path,default=base.PROJECT/"analysis/figure3/optix_full_range_frequency")
    args=parser.parse_args()
    for p in (args.vcf,base.MANIFEST_PATH,base.BAND_BA_PATH,base.DR_BA_PATH,base.LR_PATH,base.source_path("Band")):
        if not p.is_file() or not p.stat().st_size:raise RuntimeError(f"Missing required input: {p}")
    selected=base.load_ba_sites()
    selected=selected.loc[selected.chrom.eq(CHROM)&selected.pos.between(START,END)].copy()
    selected=selected.loc[selected.ref.str.len().eq(1)&selected.alt.str.len().eq(1)]
    if selected.empty:raise RuntimeError("No full-range BA>0.99 SNPs")
    groups=cohorts()
    samples=sorted(set().union(*(set(g["sample"]) for g in groups.values())))
    calls,seen=read_calls(args.vcf,selected,samples)
    recalled=base.load_recalled_pairs()
    values,sites=build_values(selected,recalled,groups,calls)
    if values.empty:raise RuntimeError("No sites pass comparative callability")
    args.output_dir.mkdir(parents=True,exist_ok=True)
    sites.to_csv(args.output_dir/"candidate_sites_and_callability.tsv",sep="\t",index=False)
    values.to_csv(args.output_dir/"frequency_plot_values.tsv.gz",sep="\t",index=False)
    pdf=plot(values,args.output_dir)
    audit={"range":[CHROM,START,END],"vcf":str(args.vcf.resolve()),
           "selection":"comparative BA>0.99; DP2 GQ20; >=80% callable in each phenotype group",
           "focal_selection":False,"focal_missing_sites_retained":True,"phase_inferred":False,
           "BA_sources":[str(base.BAND_BA_PATH),str(base.DR_BA_PATH)],
           "comparative_groups":{m:g.to_dict("records") for m,g in groups.items()},
           "candidate_module_counts":selected.groupby("module").size().to_dict(),
           "display_module_counts":values.groupby("module").apply(lambda t:len(t[KEY].drop_duplicates())).to_dict(),
           "candidate_variant_keys_absent_from_vcf":len(set(map(tuple,selected[KEY].values))-seen),
           "frequency_definition":"mean diploid pattern-associated genotype dosage across QC-passing individuals",
           "focal_definition":"adjusted diploid genotype dosage 0/0.5/1; not raw mixed-read AF",
           "pdf":str(pdf),"palette":{"absent":"#D7EDF5","midpoint":"#ECCB63","present":"#662A83","missing":"#C8C8C8"}}
    (args.output_dir/"audit.json").write_text(json.dumps(audit,indent=2)+"\n")
    for module in MODULES:
        n=values.loc[values.module.eq(module),KEY].drop_duplicates().shape[0]
        print(f"SITES\t{module}\t{n}")
    print(f"PDF\t{pdf}")
    print(f"AUDIT\t{args.output_dir/'audit.json'}")

if __name__=="__main__":main()
