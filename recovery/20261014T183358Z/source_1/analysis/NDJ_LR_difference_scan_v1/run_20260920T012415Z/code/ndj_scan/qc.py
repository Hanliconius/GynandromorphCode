"""Read-level candidate diagnostics. No claim of independent validation."""
from collections import defaultdict
import json
import numpy as np
import pandas as pd
import pysam
from scipy.stats import fisher_exact
from .common import *


def nearby_indel(read, target, flank=10):
    pos = read.reference_start
    for op, size in read.cigartuples or []:
        if op in (1, 2) and pos - flank <= target <= pos + (size if op == 2 else 0) + flank:
            return True
        if op in (0, 2, 3, 7, 8):
            pos += size
    return False


def inspect_reads(bam, chrom, pos, ref, alt, cfg):
    target = pos - 1
    observed = defaultdict(list)
    raw, lowmap = 0, 0
    for read in bam.fetch(chrom, target, target + 1):
        if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate or read.is_qcfail:
            continue
        raw += 1
        if read.mapping_quality < cfg['mapq']:
            lowmap += 1
            continue
        qp = next((q for q, r in read.get_aligned_pairs() if r == target), None)
        if qp is None or read.query_qualities is None:
            continue
        bq = int(read.query_qualities[qp])
        if bq < cfg['baseq']:
            continue
        base = read.query_sequence[qp].upper()
        key = (read.get_tag('RG') if read.has_tag('RG') else '', read.query_name)
        observed[key].append(dict(base=base, bq=bq, mq=read.mapping_quality,
            strand=int(read.is_reverse), start=read.reference_start,
            clipped=any(op == 4 for op, n in read.cigartuples or []),
            near_indel=nearby_indel(read, target),
            end_distance=min(qp, read.query_length - 1 - qp)))
    fragments, discordant = [], 0
    for values in observed.values():
        if len({x['base'] for x in values}) > 1:
            discordant += 1  # discard discordant overlapping mates
            continue
        fragments.append(max(values, key=lambda x: x['bq']))
    result = dict(raw_primary_reads=raw, low_MAPQ_reads=lowmap,
                  discordant_overlap_fragments=discordant,
                  other_base_fragments=sum(x['base'] not in (ref, alt) for x in fragments))
    for allele, label in [(ref,'ref'), (alt,'alt')]:
        values = [x for x in fragments if x['base'] == allele]
        n = len(values)
        result.update({label+'_fragments':n, label+'_forward':sum(x['strand']==0 for x in values),
            label+'_reverse':sum(x['strand']==1 for x in values),
            label+'_unique_starts':len({(x['start'],x['strand']) for x in values})})
        for metric in ['bq','mq','clipped','near_indel','end_distance']:
            result[label+'_'+metric+'_mean'] = float(np.mean([x[metric] for x in values])) if n else np.nan
    n = result['ref_fragments'] + result['alt_fragments']
    result['readcheck_AF'] = result['alt_fragments']/n if n else np.nan
    result['strand_fisher_p'] = float(fisher_exact([[result['ref_forward'], result['ref_reverse']],
                                                   [result['alt_forward'], result['alt_reverse']]])[1])
    return result


def qc(run, index):
    run, cfg = load_config(run)
    c = cfg['chromosomes'][index]
    out = chrom_dir(run,index)
    table = pd.read_csv(run/'candidates.tsv', sep='\t')
    table = table.loc[table.chrom == c['contig']].copy()
    rows = []
    with pysam.AlignmentFile(cfg['bam_left'],'rb') as left, pysam.AlignmentFile(cfg['bam_right'],'rb') as right, pysam.FastaFile(cfg['reference']) as ref:
        for row in table.to_dict('records'):
            record = dict(row)
            flags = []
            for side,bam in [('left',left), ('right',right)]:
                values = inspect_reads(bam,row['chrom'],int(row['pos']),row['ref'],row['alt'],cfg)
                record.update({side+'_'+k:v for k,v in values.items()})
                if values['ref_fragments']+values['alt_fragments'] < cfg['min_depth']:
                    flags.append(side+':low_readcheck_depth')
                if values['raw_primary_reads'] and values['low_MAPQ_reads']/values['raw_primary_reads'] > .2:
                    flags.append(side+':low_MAPQ_fraction')
                if values['strand_fisher_p'] < .01:
                    flags.append(side+':allelic_strand_imbalance')
                for allele in ['ref','alt']:
                    n = values[allele+'_fragments']
                    if n >= 3 and values[allele+'_unique_starts'] < 3:
                        flags.append(side+':'+allele+'_few_starts')
                    if n and values[allele+'_clipped_mean'] > .2:
                        flags.append(side+':'+allele+'_clipping')
                    if n and values[allele+'_near_indel_mean'] > .2:
                        flags.append(side+':'+allele+'_near_indel')
                if np.isfinite(values['readcheck_AF']) and abs(values['readcheck_AF']-row[side+'_AF'])>.15:
                    flags.append(side+':AD_readcheck_disagreement')
            seq = ref.fetch(row['chrom'],max(0,int(row['pos'])-11),int(row['pos'])+10).upper()
            record['reference_context'] = seq
            record['low_complexity_context'] = len(set(seq) & set('ACGT')) <= 2
            if record['low_complexity_context']:
                flags.append('low_complexity_context')
            record['repeat_annotation'] = 'NOT_ASSESSED'
            if cfg.get('repeats'):
                overlap = False
                with pysam.TabixFile(cfg['repeats']) as bed:
                    if row['chrom'] in bed.contigs:
                        overlap = any(True for _ in bed.fetch(row['chrom'],int(row['pos'])-1,int(row['pos'])))
                record['repeat_annotation'] = 'overlap' if overlap else 'no_overlap'
                if overlap:
                    flags.append('repeat_overlap')
            record['fisher_read_sampling_p'] = float(fisher_exact([[row['left_ref'],row['left_alt']],
                                                                  [row['right_ref'],row['right_alt']]])[1])
            record['QC_flags'] = ';'.join(flags) if flags else 'no_threshold_flags'
            record['interpretation'] = 'provisional_single_library_candidate'
            rows.append(record)
    # Empty task still emits a table with a header and completion record.
    pd.DataFrame(rows, columns=None if rows else list(table.columns)+['QC_flags','interpretation']).to_csv(out/'candidate_read_QC.tsv',sep='\t',index=False)
    write_json(out/'qc.complete.json',dict(contig=c['contig'],sites=len(rows)))
    print(f'QC_COMPLETE\t{c["contig"]}\tsites={len(rows)}',flush=True)
