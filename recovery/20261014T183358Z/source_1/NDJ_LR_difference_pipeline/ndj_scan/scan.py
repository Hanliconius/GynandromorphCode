"""Parallel chromosome extraction, empirical ranking and regional summaries."""
from array import array
from collections import Counter
import gzip
import json
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
import pysam
from .common import *


def scan(run, index):
    run, cfg = load_config(run)
    info = cfg['chromosomes'][index]
    chrom, length = info['contig'], info['length']
    out = chrom_dir(run, index)
    if (out / 'scan.complete.json').exists():
        print(f'SCAN_REUSE\t{chrom}', flush=True)
        return
    with pysam.FastaFile(cfg['reference']) as fa:
        sequence = fa.fetch(chrom).upper()
    columns = {k: array('q') for k in ['pos', 'left_ref', 'left_alt', 'right_ref', 'right_alt']}
    refs, alts = [], []
    counts = Counter()
    oldpos, seen = 0, set()
    with pysam.VariantFile(cfg['vcf']) as vcf:
        vcf.subset_samples([cfg['left'], cfg['right']])
        for rec in vcf.fetch(chrom):
            counts['read'] += 1
            if rec.pos < oldpos:
                raise ValueError(f'Unsorted input at {chrom}:{rec.pos}')
            if rec.pos != oldpos:
                seen.clear()
            key = (rec.ref, rec.alts)
            if key in seen:
                raise ValueError(f'Duplicate allele record at {chrom}:{rec.pos}')
            seen.add(key)
            oldpos = rec.pos
            if rec.alts and any(a in ('<NON_REF>', '<*>') for a in rec.alts):
                raise ValueError('gVCF detected; use the finalized joint BCF')
            if (len(rec.ref) != 1 or rec.ref not in 'ACGT' or not rec.alts or
                len(rec.alts) != 1 or rec.alts[0] not in 'ACGT' or len(rec.alts[0]) != 1 or rec.alts[0] == rec.ref):
                counts['not_biallelic_snp'] += 1
                continue
            if sequence[rec.pos - 1] != rec.ref:
                raise ValueError(f'REF mismatch at {chrom}:{rec.pos}')
            if cfg['filter'] == 'pass' and list(rec.filter) != ['PASS']:
                counts['non_PASS'] += 1
                continue
            a, b = rec.samples[cfg['left']].get('AD'), rec.samples[cfg['right']].get('AD')
            if not valid_ad(a) or not valid_ad(b):
                counts['missing_or_invalid_AD'] += 1
                continue
            if min(sum(a), sum(b)) < cfg['min_depth']:
                counts['low_AD_depth'] += 1
                continue
            if cfg['max_depth'] and max(sum(a), sum(b)) > cfg['max_depth']:
                counts['high_AD_depth'] += 1
                continue
            columns['pos'].append(rec.pos)
            for k, val in zip(['left_ref', 'left_alt', 'right_ref', 'right_alt'], (*a, *b)):
                columns[k].append(val)
            refs.append(rec.ref)
            alts.append(rec.alts[0])
            counts['retained'] += 1
    d = {k: np.asarray(v) for k, v in columns.items()}
    d.update(ref=np.asarray(refs, dtype='U1'), alt=np.asarray(alts, dtype='U1'))
    d['left_total'] = d['left_ref'] + d['left_alt']
    d['right_total'] = d['right_ref'] + d['right_alt']
    np.savez_compressed(out / 'counts.npz', **d)
    nl, nr, l, r = fractions(d)
    hex_count = np.bincount(hex_ids(l, r), minlength=len(hex_geometry()[0]))
    np.save(out / 'AF_hex_counts.npy', hex_count)
    # Full-site plots use every record. Only empirical background estimation is sampled.
    priority = priority_hash(d['pos'], index, cfg['seed'])
    take = np.argsort(priority)[:cfg['background_sites']]
    np.savez_compressed(out / 'background.npz', **{k: d[k][take] for k in COUNT_KEYS},
                        priority=priority[take], chrom_index=np.full(len(take), index))
    windows = bam_depth(cfg, chrom, length, sequence, out)
    windows.to_csv(out / 'depth.tsv', sep='\t', index=False)
    write_json(out / 'scan.complete.json', dict(contig=chrom, counts=dict(counts),
               eligible_sites=int(len(l)), full_AF_sites=int(hex_count.sum())))
    print(f'SCAN_COMPLETE\t{chrom}\tretained={len(l):,}', flush=True)


def bam_depth(cfg, chrom, length, sequence, out):
    width = cfg['window']
    n = (length + width - 1) // width
    totals = np.zeros((n, 2), dtype=np.int64)
    covered = np.zeros((n, 2), dtype=np.int64)
    command = [cfg['samtools'], 'depth', '-aa', '-s', '-q', str(cfg['baseq']),
               '-Q', str(cfg['mapq']), '-G', '0xF04', '-r', chrom, cfg['bam_left'], cfg['bam_right']]
    (out / 'depth.command.json').write_text(json.dumps(command))
    previous = 0
    with (out / 'depth.stderr.log').open('w') as err:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=err, text=True)
        try:
            for line in process.stdout:
                fields = line.rstrip().split('\t')
                if len(fields) != 4 or fields[0] != chrom:
                    raise ValueError('Unexpected samtools depth output')
                pos, dl, dr = map(int, fields[1:])
                if pos != previous + 1:
                    raise ValueError('samtools depth -aa did not emit consecutive bases')
                previous = pos
                w = (pos - 1) // width
                totals[w] += (dl, dr)
                covered[w] += (dl > 0, dr > 0)
            if process.wait() != 0:
                raise RuntimeError(f'samtools depth failed; inspect {out}/depth.stderr.log')
        except BaseException:
            process.kill()
            process.wait()
            raise
    if previous != length:
        raise RuntimeError(f'Depth ended at {previous}; expected {length} for {chrom}')
    rows = []
    for w in range(n):
        start, end = w * width, min((w + 1)*width, length)
        seq = sequence[start:end]
        acgt = sum(seq.count(b) for b in 'ACGT')
        rows.append(dict(chrom=chrom, start=start, end=end, window=w,
            left_depth=totals[w, 0] / (end-start), right_depth=totals[w, 1] / (end-start),
            left_covered=covered[w, 0] / (end-start), right_covered=covered[w, 1] / (end-start),
            gc=(seq.count('G') + seq.count('C'))/acgt if acgt else np.nan,
            acgt_fraction=acgt/(end-start)))
    return pd.DataFrame(rows)


def calibrate(run):
    run, cfg = load_config(run)
    backgrounds = []
    depths = []
    for i, c in enumerate(cfg['chromosomes']):
        folder = chrom_dir(run, i)
        if not (folder / 'scan.complete.json').is_file():
            raise RuntimeError(f'Missing scan completion: {c["contig"]}')
        if c['background']:
            with np.load(folder / 'background.npz', allow_pickle=False) as f:
                backgrounds.append({k: f[k] for k in f.files})
            depths.append(pd.read_csv(folder / 'depth.tsv', sep='\t'))
    d = {k: np.concatenate([x[k] for x in backgrounds]) for k in backgrounds[0]}
    idx = np.argsort(d['priority'])[:cfg['background_sites']]
    d = {k: v[idx] for k, v in d.items()}
    np.savez_compressed(run / 'background.npz', **d)
    depth = pd.concat(depths, ignore_index=True)
    use = (depth.left_covered >= .5) & (depth.right_covered >= .5) & (depth.acgt_fraction >= .9)
    if use.sum() < 5:
        raise RuntimeError('Too few adequately covered windows for depth normalization')
    norm_l = float(depth.loc[use, 'left_depth'].median())
    norm_r = float(depth.loc[use, 'right_depth'].median())
    if min(norm_l, norm_r) <= 0:
        raise RuntimeError('Zero median autosomal depth')
    fits = {}
    for i, c in enumerate(cfg['chromosomes']):
        take = d['chrom_index'] != i  # independent chromosome background for every score
        fit = fit_background({k: v[take] for k, v in d.items()})
        fit['excluded_contig'] = c['contig']
        fits[str(i)] = fit
    write_json(run / 'background.json', dict(version=VERSION, fits=fits,
        depth_normalizers=dict(left=norm_l, right=norm_r),
        background_contigs=[x['contig'] for x in cfg['chromosomes'] if x['background']],
        sites=int(len(idx)), seed=cfg['seed'],
        warning='Empirical genome contains biology and technical variation. Ranks are not p-values. Global depth normalization cannot identify whole-genome ploidy changes.'))
    print(f'BACKGROUND_COMPLETE\tsites={len(idx):,}', flush=True)


def score(run, index):
    run, cfg = load_config(run)
    c = cfg['chromosomes'][index]
    out = chrom_dir(run, index)
    if (out / 'score.complete.json').exists():
        print(f'SCORE_REUSE\t{c["contig"]}', flush=True)
        return
    bg = json.loads((run / 'background.json').read_text())
    rho = bg['fits'][str(index)]['rho']
    with np.load(out / 'counts.npz', allow_pickle=False) as f:
        d = {k: f[k] for k in f.files}
    with np.load(run / 'background.npz', allow_pickle=False) as f:
        b = {k: f[k][f['chrom_index'] != index] for k in f.files}
    nl, nr, l, r = fractions(d)
    delta, z, _, _, _ = statistics(d, rho)
    _, bz, _, _, _ = statistics(b, rho)
    tail, bg_n = empirical_tail(z, stratum(d), bz, stratum(b))
    table = pd.DataFrame({k: d[k] for k in COUNT_KEYS})
    table.insert(0, 'chrom', c['contig'])
    table['left_total'], table['right_total'] = nl, nr
    table['left_AF'], table['right_AF'] = l, r
    table['delta_AF'] = delta
    table['abs_delta_AF'] = np.abs(delta)
    table['heterozygosity_contrast'] = 2*r*(1-r) - 2*l*(1-l)
    table['background_variance_z'] = z
    table['abs_background_z'] = np.abs(z)
    table['background_tail_rank'] = tail
    table['background_stratum_n'] = bg_n
    table['window'] = (table.pos - 1) // cfg['window']
    table['balanced_either'] = ((l >= .2) & (l <= .8)) | ((r >= .2) & (r <= .8))
    table['balanced_left'] = (l >= .2) & (l <= .8)
    table['balanced_right'] = (r >= .2) & (r <= .8)
    table['near_fixed_left'] = (l <= .05) | (l >= .95)
    table['near_fixed_right'] = (r <= .05) | (r >= .95)
    table['left_balanced_right_fixed'] = table.balanced_left & table.near_fixed_right
    table['right_balanced_left_fixed'] = table.balanced_right & table.near_fixed_left
    # Canonical exact allele table, no GT calls used. Gzip TSV readable outside Python.
    table.to_csv(out / 'paired_sites.tsv.gz', sep='\t', index=False, compression='gzip', float_format='%.8g')
    eligible = table.loc[table.abs_delta_AF >= cfg['candidate_delta']]
    picks = []
    for metric in ['abs_background_z', 'abs_delta_AF']:
        picks.append(eligible.nlargest(cfg['candidates'], metric))
    candidate = pd.concat(picks).drop_duplicates(['chrom', 'pos', 'ref', 'alt'])
    candidate.to_csv(out / 'site_shortlist.tsv', sep='\t', index=False)
    summary = table.groupby('window').agg(sites=('pos', 'size'),
        mean_abs_delta=('abs_delta_AF', 'mean'), median_abs_delta=('abs_delta_AF', 'median'),
        mean_signed_delta=('delta_AF', 'mean'),
        median_heterozygosity_contrast=('heterozygosity_contrast', 'median'),
        left_balanced_fraction=('balanced_left', 'mean'), right_balanced_fraction=('balanced_right', 'mean'),
        left_balanced_right_fixed_fraction=('left_balanced_right_fixed', 'mean'),
        right_balanced_left_fixed_fraction=('right_balanced_left_fixed', 'mean'))
    het = table.loc[table.balanced_either].groupby('window').agg(
        balanced_either_sites=('pos', 'size'), balanced_mean_abs_delta=('abs_delta_AF', 'mean'),
        balanced_mean_signed_delta=('delta_AF', 'mean'))
    dep = pd.read_csv(out / 'depth.tsv', sep='\t').merge(summary, on='window', how='left', validate='one_to_one')
    dep = dep.merge(het, on='window', how='left', validate='one_to_one')
    dep['sites'] = dep['sites'].fillna(0).astype(int)
    dep['balanced_either_sites'] = dep['balanced_either_sites'].fillna(0).astype(int)
    dep['left_depth_relative'] = dep.left_depth / bg['depth_normalizers']['left']
    dep['right_depth_relative'] = dep.right_depth / bg['depth_normalizers']['right']
    good = (dep.left_depth > 0) & (dep.right_depth > 0)
    dep['log2_LR_depth'] = np.nan
    dep.loc[good, 'log2_LR_depth'] = np.log2(dep.loc[good, 'left_depth_relative'] / dep.loc[good, 'right_depth_relative'])
    dep['usable_depth'] = (dep.left_covered >= .5) & (dep.right_covered >= .5) & (dep.acgt_fraction >= .9)
    dep['depth_shift_flag'] = dep.usable_depth & (dep.log2_LR_depth.abs() >= cfg['depth_shift'])
    dep['contig_background'] = c['background']
    dep.to_csv(out / 'windows.tsv', sep='\t', index=False)
    write_json(out / 'score.complete.json', dict(contig=c['contig'], sites=len(table), rho=rho,
        candidates=len(candidate), windows=len(dep)))
    print(f'SCORE_COMPLETE\t{c["contig"]}\tsites={len(table):,}', flush=True)


def shortlist(run):
    run, cfg = load_config(run)
    tables, windows = [], []
    for i, c in enumerate(cfg['chromosomes']):
        out = chrom_dir(run, i)
        if not (out / 'score.complete.json').exists():
            raise RuntimeError(f'Missing score completion {c["contig"]}')
        tables.append(pd.read_csv(out / 'site_shortlist.tsv', sep='\t'))
        windows.append(pd.read_csv(out / 'windows.tsv', sep='\t'))
    table = pd.concat(tables, ignore_index=True)
    chosen = pd.concat([table.nlargest(cfg['candidates']//2, 'abs_background_z'),
                        table.nlargest(cfg['candidates']//2, 'abs_delta_AF')])
    chosen = chosen.drop_duplicates(['chrom', 'pos', 'ref', 'alt']).sort_values('abs_background_z', ascending=False)
    chosen.insert(0, 'candidate_id', [f'S{i+1:03d}' for i in range(len(chosen))])
    chosen.to_csv(run / 'candidates.tsv', sep='\t', index=False)
    win = pd.concat(windows, ignore_index=True)
    baseline = win.loc[win.contig_background & (win.sites >= cfg['min_window_sites']), 'mean_abs_delta'].dropna()
    if len(baseline) < 5:
        raise RuntimeError('Too few SNP-supported windows for regional background')
    threshold = float(baseline.quantile(.99))
    win['AF_shift_flag'] = (win.sites >= cfg['min_window_sites']) & (win.mean_abs_delta > threshold)
    win.to_csv(run / 'windows.tsv', sep='\t', index=False)
    ranked = win.loc[win.sites >= cfg['min_window_sites']].nlargest(20, 'mean_abs_delta')
    ranked.to_csv(run / 'top_AF_windows.tsv', sep='\t', index=False)
    win.loc[win.usable_depth].assign(abs_depth_shift=lambda x: x.log2_LR_depth.abs()).nlargest(20, 'abs_depth_shift').to_csv(run / 'top_depth_windows.tsv', sep='\t', index=False)
    # BED is 0-based, half-open. Merge only adjacent windows with identical signal class/direction.
    regions = []
    for row in win.itertuples():
        if not (row.AF_shift_flag or row.depth_shift_flag):
            continue
        signal = ('AF' if row.AF_shift_flag else '') + ('+depth_L' if row.depth_shift_flag and row.log2_LR_depth > 0 else '+depth_R' if row.depth_shift_flag else '')
        if regions and regions[-1]['chrom'] == row.chrom and regions[-1]['end'] == row.start and regions[-1]['signal'] == signal:
            regions[-1]['end'] = row.end
            regions[-1]['windows'] += 1
        else:
            regions.append(dict(chrom=row.chrom, start=row.start, end=row.end, signal=signal, windows=1))
    pd.DataFrame(regions, columns=['chrom','start','end','signal','windows']).to_csv(run / 'candidate_regions.bed', sep='\t', index=False, header=False)
    write_json(run / 'shortlist.complete.json', dict(sites=len(chosen), regions=len(regions),
        AF_window_threshold=threshold, threshold_description='99th percentile of observed background windows, not a significance cutoff'))
    print(f'SHORTLIST_COMPLETE\tsites={len(chosen)}\tregions={len(regions)}', flush=True)
