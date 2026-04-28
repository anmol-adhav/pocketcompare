"""
Integration tests — fetch real structures and run comparisons.
Marked with @pytest.mark.integration so they can be skipped offline:
  pytest -m "not integration"
"""
import pytest
from pocketcompare import fetch, compare


@pytest.mark.integration
def test_homolog_2zcp_5iys(tmp_path):
    p1 = fetch.fetch_pdb("2ZCP", out_dir=tmp_path)
    p2 = fetch.fetch_pdb("5IYS", out_dir=tmp_path)
    pairs = [("A:172","A:173"), ("A:176","A:177"), ("A:225","A:237")]
    report = compare.compare(p1, p2, pairs, name1="2ZCP", name2="5IYS")
    assert report.global_rmsd > 5.0,   "Global RMSD should be large for distant homologs"
    assert report.pocket_rmsd < report.global_rmsd, "Pocket should be more conserved"
    assert 0 < report.conservation_ratio < 1.0
    assert report.recommendation != ""


@pytest.mark.integration
def test_mutant_self_compare(tmp_path):
    """Comparing a structure to itself should give near-zero RMSD and UNCHANGED pocket."""
    p = fetch.fetch_pdb("2ZCP", out_dir=tmp_path)
    pocket = ["A:172", "A:176", "A:225", "A:228"]
    report = compare.compare_mutant(p, p, pocket, wt_label="WT", mut_label="Same")
    assert report.global_rmsd < 0.01
    assert report.pocket_rmsd < 0.01
    assert report.effect == "UNCHANGED"
