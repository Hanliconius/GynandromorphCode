#!/usr/bin/env python3
"""Exploratory assignment of local fragments using control-selected SNPs."""
import argparse,json,sys
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
ANCHORS={"Hmel2.5":("Hmel215003o",1404112,1405075),"HmRos":("NC_140403.1",1633657,1635225)}
def usable(r):return not (r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_duplicate or r.is_qcfail) and r.mapping_quality>=30
def key(r):return (r.get_tag("RG") if r.has_tag("RG") else "",r.query_name)
def observations(path,contig,start,end,sites=None):
    import pysam
    # One observation per fragment per position; overlapping mates count once.
    calls=defaultdict(dict);disagree=defaultdict(set);overlap=set()
    with pysam.AlignmentFile(str(path),"rb") as bam:
        for r in bam.fetch(contig,start,end):
            if not usable(r):continue
            k=key(r)
            for q,pos in r.get_aligned_pairs(matches_only=True):
                if pos is None or q is None or not start<=pos<end:continue
                if sites is not None and pos not in sites:continue
                if r.query_qualities is None or r.query_qualities[q]<30:continue
                base=r.query_sequence[q].upper()
                if base not in "ACGT" or pos in disagree[k]:continue
                if pos in calls[k] and calls[k][pos]!=base:
                    calls[k].pop(pos);disagree[k].add(pos)
                else:calls[k][pos]=base
    return calls
def assign(base,ref):
    import pysam
    cfg=json.loads((base/"run.json").read_text());contig,lo,hi=ANCHORS[ref]
    manifest=json.loads((base/"allele_evidence"/ref/"evidence_manifest.json").read_text())
    entries={e["sample"]:e for e in cfg["tasks"] if e["reference"]["name"]==ref}
    paths={e["sample"]:Path(e["bam"]) for e in manifest["samples"]}
    controls={race:[s for s,e in entries.items() if e["role"]=="control" and e.get("species")=="melpomene" and e.get("subspecies")==race and s in paths] for race in ("rosina","melpomene")}
    if any(len(ids)<2 for ids in controls.values()):raise RuntimeError("Need at least 2 controls per race")
    start,end=lo-2000,hi+2000
    counts={}
    for ids in controls.values():
        for sample in ids:
            calls=observations(paths[sample],contig,start,end)
            c=defaultdict(Counter)
            for per_fragment in calls.values():
                for pos,b in per_fragment.items():c[pos][b]+=1
            counts[sample]=c
    candidates=sorted(set().union(*(set(c) for c in counts.values())))
    sites={};rows=[]
    for pos in candidates:
        consensus={};valid=True;evidence={}
        for race,ids in controls.items():
            alleles=[]
            for sample in ids:
                c=counts[sample].get(pos,Counter());depth=sum(c.values())
                if depth<8:valid=False;break
                allele,n=c.most_common(1)[0]
                if n/depth<.9:valid=False;break
                alleles.append(allele);evidence[sample]={"counts":dict(c),"depth":depth}
            if not valid or len(set(alleles))!=1:valid=False;break
            consensus[race]=alleles[0]
        if valid and consensus["rosina"]!=consensus["melpomene"]:
            sites[pos]=consensus
            rows.append(dict(reference=ref,contig=contig,pos1=pos+1,rosina=consensus["rosina"],melpomene=consensus["melpomene"],
                inside_primer_span=lo<=pos<hi,evidence=json.dumps(evidence)))
    out=base/"SNP_assignment"/ref;out.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows,columns=["reference","contig","pos1","rosina","melpomene","inside_primer_span","evidence"]).to_csv(out/"control_selected_SNPs.tsv",sep="\t",index=False)
    summaries=[]
    for sample,path in paths.items():
        calls=observations(path,contig,start,end,sites)
        labels={};votes={}
        for k,per_fragment in calls.items():
            v=Counter()
            for pos,b in per_fragment.items():
                if b==sites[pos]["rosina"]:v["rosina"]+=1
                elif b==sites[pos]["melpomene"]:v["melpomene"]+=1
                else:v["other"]+=1
            if v["other"] or (v["rosina"] and v["melpomene"]):label="conflicting"
            elif v["rosina"]:label="rosina"
            elif v["melpomene"]:label="melpomene"
            else:label="uninformative"
            labels[k]=label;votes[k]=dict(v)
        # Propagate fragment label to both mates; no side-of-origin assignment.
        seen=set();touch=set();writers={}
        with pysam.AlignmentFile(str(path),"rb") as src:
            for label in ("rosina","melpomene","conflicting","uninformative"):
                writers[label]=pysam.AlignmentFile(str(out/(sample+"."+label+".bam")),"wb",template=src)
            try:
                for r in src.fetch(until_eof=True):
                    if not usable(r):continue
                    k=key(r);seen.add(k);label=labels.get(k,"uninformative")
                    if r.reference_name==contig and r.reference_start<hi and r.reference_end is not None and r.reference_end>lo:touch.add(k)
                    r.set_tag("YA",label,value_type="Z")
                    writers[label].write(r)
            finally:
                for w in writers.values():w.close()
        for label in writers:pysam.index(str(out/(sample+"."+label+".bam")))
        tally=Counter(labels.get(k,"uninformative") for k in seen)
        near=Counter(labels.get(k,"uninformative") for k in touch)
        n=tally["rosina"]+tally["melpomene"]
        summaries.append(dict(reference=ref,sample=sample,role=entries[sample]["role"],
            side=entries[sample]["side"],leg=entries[sample]["leg"],race=entries[sample].get("subspecies","focal"),
            informative_fragments=n,rosina=tally["rosina"],melpomene=tally["melpomene"],
            conflicting=tally["conflicting"],uninformative=tally["uninformative"],
            rosina_fraction=tally["rosina"]/n if n else np.nan,
            assigned_fragments_touching_primer_span=near["rosina"]+near["melpomene"]))
        if labels:
            pd.DataFrame([dict(read_group=k[0],query_name=k[1],assignment=v,votes=json.dumps(votes[k])) for k,v in labels.items()]).to_csv(out/(sample+".fragment_assignments.tsv.gz"),sep="\t",index=False)
    pd.DataFrame(summaries).to_csv(out/"fragment_summary.tsv",sep="\t",index=False)
    (out/"audit.json").write_text(json.dumps({"reference":ref,"selected_SNPs":len(sites),
        "within_primer_span":sum(lo<=p<hi for p in sites),"thresholds":{"MAPQ":30,"baseQ":30,"control_depth_per_individual":8,"consensus_fraction":.9},
        "focal_used_for_marker_selection":False,"markers":"race-associated local SNPs, not validated causal/phase markers",
        "side_of_origin_assigned":False,"indel_assignment_performed":False,
        "status":"exploratory; control-selection samples are not independent validation; awaiting extra controls"},indent=2)+"\n")
    print("SNP_ASSIGNMENT_COMPLETE",ref,len(sites),out,sep="\t")
def report(base):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    tables=[pd.read_csv(base/"SNP_assignment"/r/"fragment_summary.tsv",sep="\t") for r in ANCHORS]
    summary=pd.concat(tables,ignore_index=True)
    corrected=[]
    for ref in ANCHORS:
        for leg in ("T2","T3"):
            t=summary.loc[(summary.reference==ref)&(summary.role=="focal")&(summary.leg==leg)].set_index("side")
            L=t.loc["Left","rosina_fraction"];R=t.loc["Right","rosina_fraction"]
            for side,value,other in (("Left",L,R),("Right",R,L)):
                corrected.append(dict(reference=ref,leg=leg,side=side,observed=value,
                    assumed20pct_unmixed=(.8*value-.2*other)/.6))
    pd.DataFrame(corrected).to_csv(base/"SNP_assignment"/"focal_marker_fraction_assumed20pct.tsv",sep="\t",index=False)
    out=base/"Yb_d_control_selected_SNP_fragment_assignment_both_references.pdf"
    with PdfPages(out) as pdf:
        fig,axes=plt.subplots(1,2,figsize=(12,5),sharey=True)
        fig.subplots_adjust(bottom=.35,top=.8,wspace=.2)
        for ax,(ref,t) in zip(axes,zip(ANCHORS,tables)):
            for i,row in enumerate(t.itertuples()):
                color="#1764AB" if row.side=="Left" else "#D04B35" if row.side=="Right" else "#662A83" if row.race=="rosina" else "#555555"
                if row.informative_fragments:ax.bar(i,row.rosina_fraction,color=color)
                ax.text(i,1.025,str(row.informative_fragments),ha="center",fontsize=7)
            ax.set_xticks(range(len(t)),t["sample"],rotation=65,ha="right",fontsize=7)
            ax.set_ylim(0,1.12);ax.set_title(ref);ax.spines[["top","right"]].set_visible(False)
        axes[0].set_ylabel("Rosina-supporting fraction of assigned fragments")
        fig.suptitle("Yb local SNP assignment | numbers above bars = informative fragments")
        fig.text(.1,.03,"Markers selected from controls, without focal differences. Each fragment counted once; conflicting evidence remains unassigned.\n"
            "These are local race-associated SNP alleles, not yet phased to the CRE insertion or assigned to an anatomical side.\n"
            "Blank bars = no informative fragments. Additional controls will provide independent marker validation.\n"
            "Exploratory 20% corrected marker fractions are in TSV; correction assumes equal local total dosage and no allele-specific mapping bias.",fontsize=8)
        pdf.savefig(fig);plt.close(fig)
    print("PDF\t"+str(out))
def main():
    p=argparse.ArgumentParser();p.add_argument("--run-dir",required=True,type=Path)
    p.add_argument("--stage",choices=["assign","report"],required=True);p.add_argument("--task",type=int,default=0);a=p.parse_args()
    if a.stage=="assign":assign(a.run_dir,list(ANCHORS)[a.task])
    else:report(a.run_dir)
if __name__=="__main__":
    try:main()
    except Exception as e:print("ERROR: "+str(e),file=sys.stderr);sys.exit(2)
