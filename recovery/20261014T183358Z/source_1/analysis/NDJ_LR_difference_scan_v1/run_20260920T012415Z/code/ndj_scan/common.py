"""Shared count semantics, statistics and immutable run configuration."""
import hashlib
import json
from pathlib import Path
import numpy as np

VERSION = '1.0.0'
COUNT_KEYS = ['pos', 'ref', 'alt', 'left_ref', 'left_alt', 'right_ref', 'right_alt']


def fingerprint(path):
    p = Path(path).resolve()
    s = p.stat()
    if not s.st_size:
        raise ValueError(f'Empty input: {p}')
    return dict(path=str(p), bytes=s.st_size, mtime_ns=s.st_mtime_ns)


def write_json(path, value):
    path = Path(path)
    temp = path.with_name(path.name + '.partial')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temp.replace(path)


def load_config(run):
    run = Path(run).resolve()
    cfg = json.loads((run / 'config.json').read_text())
    for old in cfg['inputs']:
        if fingerprint(old['path']) != old:
            raise RuntimeError(f"Input changed since preflight: {old['path']}")
    return run, cfg


def chrom_dir(run, index):
    p = Path(run) / 'chromosomes' / f'{index:02d}'
    p.mkdir(parents=True, exist_ok=True)
    return p


def valid_ad(ad):
    return (ad is not None and len(ad) == 2 and
            all(isinstance(x, int) and x >= 0 for x in ad))


def fractions(d):
    nl = d['left_ref'] + d['left_alt']
    nr = d['right_ref'] + d['right_alt']
    l = d['left_alt'] / nl
    r = d['right_alt'] / nr
    return nl, nr, l, r


def statistics(d, rho=0.):
    nl, nr, l, r = fractions(d)
    p = (d['left_alt'] + d['right_alt'] + .5) / (nl + nr + 1.)
    v0 = p * (1 - p) * (1 / nl + 1 / nr)
    vr = p * (1 - p) * (2 - 1 / nl - 1 / nr)
    delta = l - r
    return delta, delta / np.sqrt(v0 + rho * vr), p, v0, vr


def fit_background(d):
    """Moment-based variance inflation, NOT a noise-only or calibrated null fit."""
    delta, _, p, v0, vr = statistics(d)
    use = (p >= .1) & (p <= .9)
    n = int(use.sum())
    if n < 100:
        raise RuntimeError(f'Only {n} intermediate-AF background sites; need >=100')
    raw = float(np.sum(delta[use] ** 2 - v0[use]) / np.sum(vr[use]))
    return dict(rho=float(np.clip(raw, 0, .99)), raw_rho=raw,
                intermediate_sites=n, sampled_sites=len(delta),
                interpretation='Method-of-moments descriptive beta-binomial variance inflation; real L/R differences also inflate it. Not a fitted error-only null or biological p-value.')


def stratum(d):
    nl, nr, _, _ = fractions(d)
    _, _, p, _, _ = statistics(d)
    af = np.digitize(np.minimum(p, 1 - p), [.02, .1, .25, .4])
    depth = np.digitize(np.minimum(nl, nr), [20, 40, 80])
    return af * 4 + depth


def empirical_tail(z, groups, background_z, background_groups):
    """Leave-chromosome-out observed-background tail rank, never a p-value."""
    result = np.full(len(z), np.nan)
    sizes = np.zeros(len(z), dtype=int)
    for g in np.unique(groups):
        b = np.sort(np.abs(background_z[background_groups == g]))
        mask = groups == g
        sizes[mask] = len(b)
        if len(b) >= 100:
            result[mask] = (1 + len(b) - np.searchsorted(b, np.abs(z[mask]), side='left')) / (len(b) + 1)
    return result, sizes


def hex_geometry(nx=45):
    ny = int(nx / np.sqrt(3))
    sx, sy = (1 + 2e-9) / nx, 1 / ny
    centres = np.concatenate([
        np.column_stack([np.repeat(np.arange(nx + 1), ny + 1), np.tile(np.arange(ny + 1), nx + 1)]),
        np.column_stack([np.repeat(np.arange(nx) + .5, ny), np.tile(np.arange(ny) + .5, nx)])])
    centres *= [sx, sy]
    centres[:, 0] -= 1e-9
    return centres, sx, sy


def hex_ids(l, r, nx=45):
    ny = int(nx / np.sqrt(3))
    x, y = (l + 1e-9) / ((1 + 2e-9) / nx), r * ny
    ix, iy = np.rint(x).astype(int), np.rint(y).astype(int)
    jx, jy = np.floor(x).astype(int), np.floor(y).astype(int)
    first = (x - ix)**2 + 3*(y - iy)**2 < (x - jx - .5)**2 + 3*(y - jy - .5)**2
    return np.where(first, ix * (ny + 1) + iy,
                    (nx + 1)*(ny + 1) + jx * ny + jy)


def priority_hash(pos, chrom_index, seed):
    # Stable SplitMix64 priorities: uniform reservoir selection independent of AF.
    with np.errstate(over='ignore'):
        z = np.asarray(pos, dtype=np.uint64) + np.uint64(seed + chrom_index * 1000000007)
        z += np.uint64(0x9e3779b97f4a7c15)
        z = (z ^ (z >> 30)) * np.uint64(0xbf58476d1ce4e5b9)
        z = (z ^ (z >> 27)) * np.uint64(0x94d049bb133111eb)
        return z ^ (z >> 31)
