"""Transparent sequence screens; motifs indicate potential, not observed structures."""
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict, deque
import math
import re


def revcomp(s):
    return s.translate(str.maketrans('ACGTacgt', 'TGCAtgca'))[::-1]


class Track:
    """Union intervals, with exact O(log n) overlap-length queries."""
    def __init__(self, intervals=()):
        merged = []
        for a, b in sorted(intervals):
            if b <= a:
                raise ValueError('Invalid interval')
            if merged and a <= merged[-1][1]:
                merged[-1][1] = max(b, merged[-1][1])
            else:
                merged.append([a, b])
        self.starts = [a for a, b in merged]
        self.ends = [b for a, b in merged]
        self.cumulative = [0]
        for a, b in merged:
            self.cumulative.append(self.cumulative[-1] + b-a)

    def covered_before(self, x):
        i = bisect_right(self.ends, x)
        extra = max(0, x-self.starts[i]) if i < len(self.starts) else 0
        return self.cumulative[i] + extra

    def coverage(self, a, b):
        return self.covered_before(b)-self.covered_before(a)

    def intervals(self, a, b):
        i, j = bisect_right(self.ends, a), bisect_left(self.starts, b)
        return [(max(a, x), min(b, y)) for x, y in zip(self.starts[i:j], self.ends[i:j])]


def regex_intervals(s, pattern, overlap=False):
    expression = '(?=(' + pattern + '))' if overlap else pattern
    for m in re.finditer(expression, s):
        yield (m.start(), m.start()+len(m.group(1))) if overlap else m.span()


def repeat_pairs(s, arm=10, max_gap=100):
    """Exact seed-pair screen: nearest upstream eligible arm for each right arm.

    10-bp arms; 0..100-bp gap; direct, inverted and mirror orientation.
    Bounded rolling dictionary avoids whole-genome k-mer storage. It deliberately
    does not enumerate all nested repeats or attempt imperfect repeat alignment.
    """
    hits = defaultdict(list)
    previous = defaultdict(deque)
    queue = deque()
    for j in range(arm, len(s)-arm+1):
        i = j-arm
        word = s[i:i+arm]
        if set(word) <= set('ACGT'):
            previous[word].append(i)
            queue.append((i, word))
        cutoff = j-2*arm-max_gap
        while queue and queue[0][0] < cutoff:
            old, key = queue.popleft()
            previous[key].popleft()
            if not previous[key]:
                del previous[key]
        right = s[j:j+arm]
        if not set(right) <= set('ACGT'):
            continue
        for name, key in [('direct10', right), ('inverted10', revcomp(right)), ('mirror10', right[::-1])]:
            positions = previous.get(key)
            if positions:
                start = positions[-1]
                # upstream arm ends <= j; gap = j - (start+arm)
                if j-start-arm <= max_gap:
                    hits[name].append((start, j+arm))
                    if name == 'mirror10' and j-start-arm <= 8:
                        pur = sum(x in 'AG' for x in right)/arm
                        if pur >= .9 or pur <= .1:
                            hits['triplex_proxy'].append((start, j+arm))
        # Collapse periodically, bounding storage in very repetitive chromosomes.
        if j % 100000 == 0:
            for name in list(hits):
                t = Track(hits[name])
                hits[name] = list(zip(t.starts, t.ends))
    return hits


def scan(s):
    upper = s.upper()
    patterns = {
        'ambiguous': '[^ACGT]+',
        'softmask': '[acgt]+',
        'homopolymer8': r'A{8,}|C{8,}|G{8,}|T{8,}',
        'AT_run20': '[AT]{20,}',
        'purine_pyrimidine20': '[AG]{20,}|[CT]{20,}',
        'RY_alternating12_proxy': '(?:[AG][CT]){6,}[AG]?|(?:[CT][AG]){6,}[CT]?',
        'telomere_TTAGG_3': '(?:TTAGG){3,}|(?:CCTAA){3,}',
        'G4_canonical': r'G{3,7}[ACGT]{1,7}G{3,7}[ACGT]{1,7}G{3,7}[ACGT]{1,7}G{3,7}|C{3,7}[ACGT]{1,7}C{3,7}[ACGT]{1,7}C{3,7}[ACGT]{1,7}C{3,7}',
    }
    tracks = {k: Track(regex_intervals(s if k == 'softmask' else upper, p,
                                     overlap=k == 'G4_canonical'))
              for k, p in patterns.items()}
    tandem = []
    for period in range(2, 7):
        # >=3 copies and >=12 bp. Primitive period not inferred; union avoids double counting.
        copies = max(3, math.ceil(12/period))
        pattern = r'([ACGT]{%d})\1{%d,}' % (period, copies-1)
        tandem.extend(m.span() for m in re.finditer(pattern, upper))
    tracks['tandem2_6'] = Track(tandem)
    pairs = repeat_pairs(upper)
    for name in ['direct10', 'inverted10', 'mirror10', 'triplex_proxy']:
        tracks[name] = Track(pairs.get(name, []))
    return tracks


def score(s, tracks, start, end):
    seq = s[start:end].upper()
    n = len(seq)
    if n != end-start or n == 0:
        raise ValueError('Truncated/empty interval')
    counts = Counter(seq)
    valid = sum(counts[b] for b in 'ACGT')
    gc = (counts['G']+counts['C'])/valid if valid else None
    entropy = -sum((counts[b]/valid)*math.log2(counts[b]/valid)
                   for b in 'ACGT' if counts[b]) if valid else None
    result = {'GC': gc, 'valid_fraction': valid/n, 'base_entropy': entropy}
    result.update({key+'_fraction': track.coverage(start, end)/n for key, track in tracks.items()})
    return result
