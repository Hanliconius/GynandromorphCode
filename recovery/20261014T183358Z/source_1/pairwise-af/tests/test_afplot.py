import gzip
import json
from pathlib import Path

import numpy as np
import pysam
import pytest

import plot_pairwise_AF as af

FIXTURE = Path(__file__).resolve().parents[1] / "examples" / "tiny.vcf"


def arguments(tmp_path, *extra):
    return af.parser().parse_args([
        "--vcf", str(FIXTURE), "--left-sample", "Left", "--right-sample", "Right",
        "--output-prefix", str(tmp_path / "result"), "--progress-every", "0", "--bin-shape", "square", "--bins", "81", *extra])


def altered(tmp_path, transformation):
    path = tmp_path / "altered.vcf"
    path.write_text(transformation(FIXTURE.read_text()))
    return path


def test_exact_bin_counts_sample_order_and_ad_denominator(tmp_path):
    args = arguments(tmp_path)
    hist, audit = af.accumulate(args)
    # Header order is Right, Decoy, Left. FORMAT/DP deliberately disagrees with sum(AD).
    expected = np.histogram2d([0, .2, .5, 1, .5, .7], [1, .8, 0, 0, .5, .3],
                              bins=(np.linspace(0, 1, 82), np.linspace(0, 1, 82)))[0]
    np.testing.assert_array_equal(hist, expected)
    assert audit["records_read"] == 16
    assert audit["records_retained"] == 6
    assert audit["excluded_first_reason"] == {
        "missing_or_invalid_ad": 5, "below_min_depth": 2, "not_biallelic_snp": 3}


@pytest.mark.parametrize("extra, expected", [
    (["--filter", "pass"], 4), (["--max-depth", "20"], 5),
    (["--min-depth", "1"], 7), (["--min-depth", "40"], 1),
])
def test_filters(tmp_path, extra, expected):
    hist, _ = af.accumulate(arguments(tmp_path, *extra))
    assert hist.sum() == expected


def test_reversing_samples_transposes_histogram(tmp_path):
    args = arguments(tmp_path)
    original, _ = af.accumulate(args)
    args.left_sample, args.right_sample = args.right_sample, args.left_sample
    reverse, _ = af.accumulate(args)
    np.testing.assert_array_equal(reverse, original.T)


@pytest.mark.parametrize("mode", ["vcf.gz", "bcf"])
def test_bgzf_and_bcf_no_index(tmp_path, mode):
    path = tmp_path / ("converted." + mode)
    if mode == "vcf.gz":
        pysam.tabix_compress(str(FIXTURE), str(path))
    else:
        with pysam.VariantFile(str(FIXTURE)) as source:
            with pysam.VariantFile(str(path), "wb", header=source.header) as output:
                for record in source:
                    output.write(record)
    plain, _ = af.accumulate(arguments(tmp_path))
    converted, _ = af.accumulate(arguments(tmp_path, "--vcf", str(path)))
    np.testing.assert_array_equal(plain, converted)


def test_whitelist_exact_alleles_and_missing_sites(tmp_path):
    path = tmp_path / "sites.tsv.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("CHROM\tPOS\tREF\tALT\nchr1\t2\tA\tG\nchr1\t3\tA\tT\nchr1\t900\tA\tG\n")
    targets = af.load_sites(path)
    hist, audit = af.accumulate(arguments(tmp_path), targets)
    assert hist.sum() == 1
    assert audit["sites_not_matched_in_selected_contigs"] == 2


@pytest.mark.parametrize("content", ["chr1\t0\tA\tG\n", "chr1\t1\tA,G\n",
                                        "chr1\t1\tA\tG\nchr1\t1\tA\tG\n", ""])
def test_bad_whitelists(tmp_path, content):
    path = tmp_path / "sites.tsv"
    path.write_text(content)
    with pytest.raises(ValueError):
        af.load_sites(path)


@pytest.mark.parametrize("extra", [["--left-sample", "absent"], ["--contig", "missing"]])
def test_wrong_names(tmp_path, extra):
    with pytest.raises(ValueError):
        af.accumulate(arguments(tmp_path, *extra))


@pytest.mark.parametrize("extra", [["--left-sample", "Right"], ["--min-depth", "0"],
    ["--max-depth", "5"], ["--vmin", "1"], ["--vmin", "5", "--vmax", "1"],
    ["--mixing-left", ".2"], ["--mixing-left", ".2", "--mixing-right", "nan"], ["--bins", "1"]])
def test_invalid_options(tmp_path, extra):
    with pytest.raises(ValueError):
        af.check_args(arguments(tmp_path, *extra))


def test_gvcf_fails_even_after_usable_records(tmp_path):
    path = altered(tmp_path, lambda s: s + "chr1\t17\t.\tA\tG,<NON_REF>\t60\tPASS\t.\tGT:AD\t0/1:5,5,0\t0/1:5,5,0\t0/1:5,5,0\n")
    args = arguments(tmp_path, "--vcf", str(path))
    with pytest.raises(ValueError, match="gVCF"):
        af.accumulate(args)


def test_duplicate_fails(tmp_path):
    path = altered(tmp_path, lambda s: s + s.splitlines()[-1] + "\n")
    with pytest.raises(ValueError, match="Duplicate"):
        af.accumulate(arguments(tmp_path, "--vcf", str(path)))


def test_unsorted_fails(tmp_path):
    path = altered(tmp_path, lambda s: s + s.splitlines()[7] + "\n")
    with pytest.raises(ValueError, match="Unsorted"):
        af.accumulate(arguments(tmp_path, "--vcf", str(path)))


def test_no_ad_header(tmp_path):
    path = altered(tmp_path, lambda s: "\n".join(l for l in s.splitlines() if not l.startswith("##FORMAT=<ID=AD")) + "\n")
    with pytest.raises(ValueError, match="FORMAT/AD"):
        af.accumulate(arguments(tmp_path, "--vcf", str(path)))


def test_bad_ad_declaration(tmp_path):
    path = altered(tmp_path, lambda s: s.replace("AD,Number=R", "AD,Number=1"))
    with pytest.raises(ValueError, match="Number=R"):
        af.accumulate(arguments(tmp_path, "--vcf", str(path)))


def test_empty_after_filtering(tmp_path):
    with pytest.raises(ValueError, match="No usable"):
        af.accumulate(arguments(tmp_path, "--min-depth", "10000"))


def cli_args(tmp_path):
    return ["--vcf", str(FIXTURE), "--left-sample", "Left", "--right-sample", "Right",
            "--output-prefix", str(tmp_path / "result"), "--bin-shape", "square", "--bins", "81"]


def test_end_to_end_and_no_overwrite(tmp_path):
    assert af.main(cli_args(tmp_path)) == 0
    audit_file = tmp_path / "result.AF.audit.json"
    audit = json.loads(audit_file.read_text())
    assert audit["records_retained"] == 6
    with gzip.open(tmp_path / "result.AF_bins.tsv.gz", "rt") as handle:
        next(handle)
        rows = [line.split("\t") for line in handle]
    assert len(rows) == 81 * 81
    assert sum(int(row[6]) for row in rows) == 6
    assert sum(float(row[7]) for row in rows) == pytest.approx(1)
    for filename, meta in audit["outputs"].items():
        assert af.sha256(Path(filename)) == meta["sha256"]
    before = audit_file.read_bytes()
    assert af.main(cli_args(tmp_path)) == 1
    assert before == audit_file.read_bytes()


def test_validation_only_no_outputs(tmp_path):
    assert af.main(cli_args(tmp_path) + ["--validate-only"]) == 0
    assert not list(tmp_path.iterdir())


def test_failed_run_no_outputs(tmp_path):
    assert af.main(cli_args(tmp_path) + ["--min-depth", "9999"]) == 1
    assert not list(tmp_path.iterdir())


def test_fixed_proportion_scale_and_overlay(tmp_path):
    assert af.main(cli_args(tmp_path) + ["--scale", "proportion", "--vmin", "0.01", "--vmax", ".1",
        "--mixing-left", ".197", "--mixing-right", ".2024"]) == 0
    audit = json.loads((tmp_path / "result.AF.audit.json").read_text())
    assert audit["plot"]["colour_limits_clip_positive_bins"]


def test_format_order_not_assumed(tmp_path):
    lines = FIXTURE.read_text().splitlines()
    for index, line in enumerate(lines):
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        fields[8] = ":".join(reversed(fields[8].split(":")))
        for sample_index in range(9, 12):
            fields[sample_index] = ":".join(reversed(fields[sample_index].split(":")))
        lines[index] = "\t".join(fields)
    path = altered(tmp_path, lambda s: "\n".join(lines) + "\n")
    original, _ = af.accumulate(arguments(tmp_path))
    reordered, _ = af.accumulate(arguments(tmp_path, "--vcf", str(path)))
    np.testing.assert_array_equal(original, reordered)


def test_truncated_text_record_fails(tmp_path):
    path = altered(tmp_path, lambda s: s + "chr1\t17\t.\tA\tG\t60\tPASS\t.\tGT:AD\t0/1:5,5\n")
    with pytest.raises((OSError, ValueError)):
        af.accumulate(arguments(tmp_path, "--vcf", str(path)))


def test_display_keeps_singletons_visible(tmp_path):
    hist, _ = af.accumulate(arguments(tmp_path))
    metadata = af.draw(hist, arguments(tmp_path), tmp_path)
    assert "white reserved for zero" in metadata["colormap"]


@pytest.mark.parametrize("nx", [2, 17, 65, 81])
def test_streaming_hex_matches_matplotlib_including_boundaries(nx):
    import matplotlib.pyplot as plt
    rng = np.random.default_rng(192)
    x = np.r_[rng.random(10000), 0, 1, 0, 1, .5, .2, .8]
    y = np.r_[rng.random(10000), 0, 0, 1, 1, .5, .8, .2]
    expected_centres, _, _ = af.hex_geometry(nx)
    counts = np.zeros(len(expected_centres), dtype=int)
    for left, right in zip(x, y):
        counts[af.hex_index(left, right, nx)] += 1
    fig, ax = plt.subplots()
    actual = ax.hexbin(x, y, gridsize=nx, extent=(0, 1, 0, 1), mincnt=0)
    np.testing.assert_array_equal(counts, actual.get_array())
    np.testing.assert_allclose(expected_centres, actual.get_offsets(), atol=1e-15)
    assert counts.sum() == len(x)
    plt.close(fig)


def test_hex_end_to_end_and_filter_invariance(tmp_path):
    args = arguments(tmp_path, "--bin-shape", "hex", "--bins", "65")
    hist, audit = af.accumulate(args)
    square, old = af.accumulate(arguments(tmp_path))
    assert audit == old
    assert hist.sum() == square.sum() == 6
    centres, _, _ = af.hex_geometry(65)
    assert len(hist) == len(centres) == 4913
    assert af.main(cli_args(tmp_path) + ["--bin-shape", "hex", "--bins", "65"]) == 0
    import csv
    with gzip.open(tmp_path / "result.AF_bins.tsv.gz", "rt") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 4913
    assert sum(int(r["count"]) for r in rows) == 6
    saved = json.loads((tmp_path / "result.AF.audit.json").read_text())
    assert saved["plot"]["bin_shape"] == "hex"
    assert saved["plot"]["bin_edges"] is None


def test_new_defaults_and_white_zero():
    from cmcrameri import cm
    args = af.parser().parse_args(["--vcf", "x.vcf", "--left-sample", "L",
                                  "--right-sample", "R", "--output-prefix", "out"])
    assert (args.bins, args.bin_shape, args.cmap) == (65, "hex", "batlowW_r")
    cmap = af.colourmap("batlowW_r")
    np.testing.assert_array_equal(cmap(np.ma.masked), [1, 1, 1, 1])
    assert not np.allclose(cmap(0.), [1, 1, 1, 1])
    np.testing.assert_allclose(cmap(0.), cm.batlowW_r(.04))
    np.testing.assert_allclose(cmap(1.), cm.batlowW_r(1.))
