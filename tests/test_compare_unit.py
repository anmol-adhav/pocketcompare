"""Unit tests — no network, no PDB files needed."""
import numpy as np
import pytest
from pocketcompare.compare import (
    _radius, _mean_pairwise, _global_rmsd, HomologReport, MutantReport
)


def test_radius_centred():
    # 4 points equidistant from origin
    coords = np.array([[1,0,0],[0,1,0],[-1,0,0],[0,-1,0]], dtype=float)
    assert abs(_radius(coords) - 1.0) < 1e-6


def test_mean_pairwise_triangle():
    # equilateral triangle with side = 1
    coords = np.array([[0,0,0],[1,0,0],[0.5, 0.866,0]], dtype=float)
    d = _mean_pairwise(coords)
    assert abs(d - 1.0) < 0.01


def test_mean_pairwise_single():
    coords = np.array([[1.0, 2.0, 3.0]])
    assert _mean_pairwise(coords) == 0.0


def test_homolog_report_fields():
    r = HomologReport(
        structure1="A", structure2="B",
        global_rmsd=8.0, pocket_rmsd=4.0,
        conservation_ratio=0.5,
        recommendation="ACTIVITY first, then test FLUX",
        reasoning="test"
    )
    assert r.conservation_ratio == 0.5
    assert "FLUX" in r.recommendation


def test_mutant_report_effect():
    r = MutantReport(
        wt_label="WT", mutant_label="MUT",
        global_rmsd=0.3, pocket_rmsd=1.2,
        wt_pocket_radius=8.0, mut_pocket_radius=9.0,
        radius_change=1.0,
        wt_mean_ca_dist=12.0, mut_mean_ca_dist=13.0,
        ca_dist_change=1.0,
        effect="EXPANDED",
        engineering_impact="Pocket expanded."
    )
    assert r.effect == "EXPANDED"
    assert r.radius_change > 0
