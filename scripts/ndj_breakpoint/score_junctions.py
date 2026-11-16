#!/usr/bin/env python3
"""Score verified oriented alternate junctions, never a reference breakpoint guess."""
import argparse
import csv


def assess(left, right, inserted, observed, cap=50):
    left,right,inserted,observed=[s.upper() for s in [left,right,inserted,observed]]
    if not left or not right or not observed or any(set(s)-set('ACGT') for s in [left,right,inserted,observed]):
        raise ValueError('Require nonempty unambiguous flanks and assembled junction')
    candidates=[]
    for k in range(min(cap,len(left),len(right))+1):
        if k and (inserted or left[-k:]!=right[:k]):continue
        if observed == left+inserted+right[k:]:candidates.append(k)
    if not candidates:raise ValueError('Junction does not reconcile with supplied oriented flanks/insertion')
    k=max(candidates)
    return dict(microhomology_bp=k,microhomology=left[-k:] if k else '',
        insertion_bp=len(inserted),insertion_sequence=inserted,
        homology_at_search_limit=k==min(cap,len(left),len(right)),
        interpretation='compatible junction signature; not a unique repair-mechanism assignment')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();rows=[]
    with open(a.input) as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if not r['support_evidence'].strip():raise ValueError('Read/assembly support evidence required')
            rows.append(dict(junction_id=r['junction_id'],support_evidence=r['support_evidence'],
                **assess(r['left_flank'],r['right_flank'],r['inserted'],r['junction_sequence'])))
    if not rows:raise ValueError('No junctions')
    with open(a.output,'w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)


if __name__=='__main__':main()
