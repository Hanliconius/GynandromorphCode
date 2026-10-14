#!/usr/bin/env python3
"""Write synthetic layout-demo data to stdout; not a biological calibration."""
import random
import sys

rng = random.Random(421)
print('##fileformat=VCFv4.2')
print('##contig=<ID=synthetic,length=100000>')
print('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Synthetic allelic depths">')
print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tRight\tLeft')
pos = 0
for left in (0, 0.5, 1):
    for right in (0, 0.5, 1):
        ql, qr = .8 * left + .2 * right, .2 * left + .8 * right
        for replicate in range(1200):
            pos += 1
            nl, nr = rng.randint(15, 65), rng.randint(15, 65)
            al = sum(rng.random() < ql for _ in range(nl))
            ar = sum(rng.random() < qr for _ in range(nr))
            print(f'synthetic\t{pos}\t.\tA\tG\t60\tPASS\t.\tAD\t{nr-ar},{ar}\t{nl-al},{al}')
