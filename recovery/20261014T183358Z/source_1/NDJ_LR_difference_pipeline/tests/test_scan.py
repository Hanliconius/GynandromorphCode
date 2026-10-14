import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
import pytest
from ndj_scan.common import *


def counts(l,r,n=100):
    l=np.asarray(l);r=np.asarray(r)
    return dict(left_ref=n-l,left_alt=l,right_ref=n-r,right_alt=r)


def test_missing_AD_is_not_zero():
    assert not valid_ad(None)
    assert not valid_ad((0,None))
    assert not valid_ad((2,3,0))
    assert valid_ad((0,0))


def test_side_swap_reverses_signed_score_only():
    a=counts([20,50],[80,50]);b=counts([80,50],[20,50])
    da,za,*_=statistics(a,.1);db,zb,*_=statistics(b,.1)
    np.testing.assert_allclose(da,-db);np.testing.assert_allclose(za,-zb)


def test_REF_ALT_flip_keeps_absolute_difference():
    a=counts([20,50],[80,30]);b=counts([80,50],[20,70])
    np.testing.assert_allclose(np.abs(statistics(a)[0]),np.abs(statistics(b)[0]))


def test_overdispersion_reduces_score():
    assert abs(statistics(counts([20],[80]),.1)[1][0])<abs(statistics(counts([20],[80]),0)[1][0])


def test_background_is_inflated_by_true_difference():
    rng=np.random.default_rng(1)
    a=counts(rng.binomial(100,.5,5000),rng.binomial(100,.5,5000))
    b=counts(rng.binomial(100,.8,5000),rng.binomial(100,.2,5000))
    assert fit_background(a)['rho']<.01
    assert fit_background(b)['rho']>.3


def test_empirical_tail_sparse_is_unresolved():
    rank,n=empirical_tail(np.array([10.,1.]),np.array([0,1]),np.arange(200)/100,np.zeros(200))
    assert rank[0]==1/201 and np.isnan(rank[1]) and n[1]==0


def test_hex_conservation_and_edges():
    rng=np.random.default_rng(3)
    l=np.r_[rng.random(1000),0,1,0,1];r=np.r_[rng.random(1000),0,1,1,0]
    ids=hex_ids(l,r);assert ids.min()>=0 and ids.max()<len(hex_geometry()[0])
    assert np.bincount(ids).sum()==len(l)


def test_hex_matches_matplotlib():
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rng=np.random.default_rng(5);l=rng.random(10000);r=rng.random(10000)
    fig,ax=plt.subplots();h=ax.hexbin(l,r,gridsize=45,extent=(0,1,0,1),mincnt=None)
    np.testing.assert_array_equal(np.bincount(hex_ids(l,r),minlength=len(hex_geometry()[0])),h.get_array())
    plt.close(fig)


def test_changed_input_fails(tmp_path):
    p=tmp_path/'in';p.write_text('hello');s=fingerprint(p)
    write_json(tmp_path/'config.json',dict(inputs=[s]))
    p.write_text('different')
    with pytest.raises(RuntimeError,match='Input changed'):load_config(tmp_path)


def test_read_overlap_collapsed_and_conflicts_dropped(tmp_path):
    import pysam
    from ndj_scan.qc import inspect_reads
    p=tmp_path/'reads.bam'
    hdr={'HD':{'SO':'coordinate'},'SQ':[{'SN':'chr','LN':1000}]}
    with pysam.AlignmentFile(str(p),'wb',header=hdr) as out:
        for name,base in [('same','A'),('same','A'),('conflict','A'),('conflict','G'),('unique','G')]:
            r=pysam.AlignedSegment();r.query_name=name;r.query_sequence=base*100;r.flag=0
            r.reference_id=0;r.reference_start=0;r.cigar=[(0,100)];r.mapping_quality=60
            r.query_qualities=pysam.qualitystring_to_array('I'*100);out.write(r)
    pysam.index(str(p))
    with pysam.AlignmentFile(str(p),'rb') as f:
        v=inspect_reads(f,'chr',50,'A','G',dict(mapq=20,baseq=20))
    assert v['ref_fragments']==1 and v['alt_fragments']==1 and v['discordant_overlap_fragments']==1
