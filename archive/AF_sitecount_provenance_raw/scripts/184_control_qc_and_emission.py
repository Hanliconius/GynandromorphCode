import argparse
import csv
import gzip
import json
import math
import os
import time
from collections import Counter
from pathlib import Path

CONTROLS = '''ERR3790122 ERR3790123 ERR2298257 ERR2298260 ERR3790118
ERR2298201 ERR2298265 ERR2298235 ERR2298236 ERR2298216 ERR1307113
ERR1307112 ERR1307111 ERR3790124 ERR2298208 ERR2298263 ERR2298221
SRR11243339'''.split()

FOCAL = {
    'TN_Left_T2', 'TN_Left_T3', 'TN_Right_T2', 'TN_Right_T3',
    'TN_Left_merged', 'TN_Right_merged'
}
ROOT = Path('analysis/genomewide_mixing_corrected_v1')


def classify(cell):
    fields = cell.split(':')
    if len(fields) != 5:
        raise ValueError('Expected GT:AD:DP:PL:SP')
    gt, ad, dp, pl, sp = fields
    gt = gt.replace('|', '/')
    state = {'0/0': 0, '0/1': 1, '1/0': 1, '1/1': 2}.get(gt)
    if state is None:
        return 'missing_or_non_diploid_GT', None
    try:
        depths = [int(x) for x in ad.split(',')]
        depth = int(dp)
        likelihoods = [int(x) for x in pl.split(',')]
    except ValueError:
        return 'missing_or_invalid_AD_DP_PL', None
    if (len(depths) != 2 or len(likelihoods) != 3 or
            min(depths + likelihoods + [depth]) < 0):
        return 'invalid_AD_DP_PL_shape_or_range', None
    total = sum(depths)
    if total < 2 or depth < 2:
        return 'depth_below_2', None
    best = min(likelihoods)
    if likelihoods[state] != best:
        return 'GT_disagrees_with_best_PL', None
    gap = sorted(likelihoods)[1] - best
    if gap < 20:
        return 'PL_gap_below_20', None
    return 'PASS', (state, total, depths[1], int(depth != total))


def scan(chrom):
    if not 1 <= chrom <= 20:
        raise ValueError('Control calibration uses autosomes 1-20 only')
    source = ROOT / 'inputs' / f'chr{chrom}_NC_{140388+chrom}.1.qc.tsv.gz'
    out = ROOT / 'control_qc'
    out.mkdir(parents=True, exist_ok=True)
    target = out / f'chr{chrom}.hist.tsv.gz'
    audit = out / f'chr{chrom}.audit.json'
    signature = [source.stat().st_size, source.stat().st_mtime_ns]

    if audit.exists() and target.exists():
        previous = json.loads(audit.read_text())
        if previous.get('source_signature') == signature:
            print(f'REUSE chromosome={chrom}', flush=True)
            return
        raise RuntimeError('Existing QC has a different input; no overwrite')

    hist = Counter()
    qc = {sample: Counter() for sample in CONTROLS}
    sites = Counter()
    count = 0

    with gzip.open(source, 'rt') as handle:
        header = next(handle).rstrip('\n').split('\t')
        if header[:6] != ['#CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER']:
            raise ValueError('Unexpected input header')
        if len(header) != 30 or set(header[6:]) != set(CONTROLS) | FOCAL:
            raise ValueError('Unexpected 24-sample identities')
        indices = [(s, header.index(s)) for s in CONTROLS]

        for count, line in enumerate(handle, 1):
            row = line.rstrip('\n').split('\t')
            if len(row) != 30:
                raise ValueError(f'Invalid row {count}')
            ref, alt = row[2:4]
            if not ref or not alt or any(x not in 'ACGT' for x in ref+alt):
                sites['non_sequence_or_multiallelic_excluded'] += 1
                continue
            kind = 'snp' if len(ref) == len(alt) == 1 else 'other_biallelic'
            sites[kind] += 1

            for sample, index in indices:
                reason, value = classify(row[index])
                qc[sample][reason] += 1
                if value is not None:
                    state, total, alt_depth, dp_diff = value
                    hist[kind, state, total, alt_depth] += 1
                    qc[sample][f'pass_GT_{state}'] += 1
                    qc[sample]['pass_DP_ne_AD_sum'] += dp_diff

            if count % 250000 == 0:
                print(f'SCAN chromosome={chrom} rows={count}', flush=True)

    with source.with_name(source.name+'.complete.tsv').open() as handle:
        expected = int(dict(csv.reader(handle, delimiter='\t'))['rows'])
    if count != expected:
        raise ValueError(f'Row count {count} != {expected}')

    temp = str(target) + f'.{os.getpid()}.new'
    with gzip.open(temp, 'wt') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['variant_type', 'GT_dosage', 'AD_total', 'AD_alt', 'count'])
        for key, multiplicity in sorted(hist.items()):
            writer.writerow([*key, multiplicity])
    os.replace(temp, target)

    result = dict(
        chromosome=chrom, rows=count, sites=dict(sites),
        source_signature=signature, controls=qc,
        rule='DP>=2; sumAD>=2; GT at unique PL minimum; PL gap>=20',
        focal_used=False, allele_balance_filter=False
    )
    temp = str(audit) + f'.{os.getpid()}.new'
    Path(temp).write_text(json.dumps(result, indent=2)+'\n')
    os.replace(temp, audit)
    print(f'COMPLETE chromosome={chrom} rows={count} bins={len(hist)}', flush=True)


def fit():
    import numpy as np
    from scipy.optimize import minimize_scalar
    from scipy.special import betaln, gammaln

    pooled = Counter()
    qc = {sample: Counter() for sample in CONTROLS}
    for chrom in range(1, 21):
        folder = ROOT / 'control_qc'
        audit = json.loads((folder / f'chr{chrom}.audit.json').read_text())
        for sample in CONTROLS:
            qc[sample].update(audit['controls'][sample])
        with gzip.open(folder / f'chr{chrom}.hist.tsv.gz', 'rt') as handle:
            for row in csv.DictReader(handle, delimiter='\t'):
                key = (
                    row['variant_type'], int(row['GT_dosage']),
                    int(row['AD_total']), int(row['AD_alt'])
                )
                pooled[key] += int(row['count'])

    out = ROOT / 'control_qc'
    summaries = []
    for scope in ['all_biallelic', 'snps_only']:
        selected = [
            (g, n, k, c) for (kind, g, n, k), c in pooled.items()
            if scope == 'all_biallelic' or kind == 'snp'
        ]
        hom_n = sum(n*c for g, n, k, c in selected if g != 1)
        hom_error = sum(
            (k if g == 0 else n-k)*c
            for g, n, k, c in selected if g != 1
        )
        het = [(n, k, c) for g, n, k, c in selected if g == 1]
        if hom_n <= 0 or sum(c for n, k, c in het) < 50:
            raise ValueError(f'Insufficient QC controls: {scope}')

        n, k, weight = np.asarray(het, dtype=float).T
        raw_error = hom_error / hom_n

        # Bounds retained from the original emission-calibration procedure.
        error = float(np.clip(raw_error, 0.0005, 0.05))
        mean = float(np.sum(k*weight) / np.sum(n*weight))
        target = float(np.clip((mean-error)/(1-2*error), 0.01, 0.99))
        bias = math.log(target/(1-target))
        expected = error + (1-2*error)*target
        constant = gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)

        def objective(log_phi):
            phi = math.exp(log_phi)
            a, b = expected*phi, (1-expected)*phi
            ll = constant + betaln(k+a, n-k+b)-betaln(a, b)
            return -float(np.sum(weight*ll))

        opt = minimize_scalar(
            objective, bounds=(math.log(2), math.log(5000)),
            method='bounded', options={'xatol': 1e-5}
        )
        if not opt.success or not math.isfinite(opt.fun):
            raise RuntimeError(f'Emission fit failed: {scope}: {opt.message}')

        phi = math.exp(opt.x)
        summaries.append(dict(
            scope=scope, raw_error_rate=raw_error,
            error_rate=error, heterozygote_mean_AF=mean,
            bias_log_odds=bias, concentration=phi,
            heterozygotes=int(sum(weight)),
            error_clipped=raw_error != error,
            concentration_near_bound=phi < 2.02 or phi > 4950,
            selected_for_founder_fit=False
        ))

    with (out / 'emission_calibration.diagnostic.json').open('w') as handle:
        json.dump(summaries, handle, indent=2)

    reasons = sorted({r for counts in qc.values() for r in counts})
    with (out / 'control_QC_summary.tsv').open('w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['sample', *reasons])
        for sample in CONTROLS:
            writer.writerow([sample, *[qc[sample][r] for r in reasons]])

    print(json.dumps(summaries, indent=2), flush=True)
    print('ENDPOINT=CONTROL_QC_AND_EMISSION_DIAGNOSTICS_READY', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=['scan', 'fit'])
    parser.add_argument('--chromosome', type=int)
    args = parser.parse_args()
    started = time.monotonic()
    print(f'START mode={args.mode}', flush=True)
    if args.mode == 'scan':
        scan(args.chromosome)
    else:
        fit()
    print(f'END elapsed_seconds={time.monotonic()-started:.1f}', flush=True)
