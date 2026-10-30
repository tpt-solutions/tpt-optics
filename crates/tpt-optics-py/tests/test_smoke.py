"""Smoke tests for the tpt_optics Python extension.

Build the extension first, for example with
    maturin develop --release   (in crates/tpt-optics-py)
or copy the cdylib built by `cargo build -p tpt-optics-py --features extension-module`
to a directory on sys.path as tpt_optics.pyd (Windows) or tpt_optics.so.
"""

import json

import pytest

tpt = pytest.importorskip("tpt_optics")


def test_glass_index_matches_the_sellmeier_model():
    assert abs(tpt.glass_index("N-BK7", 632.8) - 1.5151) < 1e-3


def test_unknown_glass_raises_key_error():
    with pytest.raises(KeyError):
        tpt.glass_index("NOT-A-GLASS", 587.6)


def test_singlet_report_gives_the_focal_length_and_a_focused_spot():
    report = json.loads(tpt.singlet_report(100.0, -100.0, 5.0, 12.5, "N-BK7"))
    assert 0.09 < report["efl"] < 0.1
    assert report["fields"][0]["rms_spot"] < 5e-4


def test_achromat_ranking_is_sorted_by_spot():
    rows = tpt.rank_achromat_pairs(100.0, ["bk7"], ["f2", "n-sf11"])
    assert len(rows) == 2
    assert rows[0][2] <= rows[1][2]
