"""
Core comparison engine.

Two modes
---------
compare()         Homolog mode: two different proteins
                  -> ACTIVITY vs FLUX recommendation

compare_mutant()  Mutant mode: WT vs engineered variant of same protein
                  -> pocket EXPANDED / CONTRACTED / UNCHANGED
                  -> steric hindrance relief assessment
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser, Superimposer


# ── internal helpers ──────────────────────────────────────────────────────────

def _load(path: Path):
    return PDBParser(QUIET=True).get_structure("s", str(path))


def _get_ca(structure, chain_id: str, resnum: int):
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for res in chain:
                    if res.id[1] == resnum and "CA" in res:
                        return res["CA"]
        break
    return None


def _all_ca(path: Path) -> np.ndarray:
    coords = []
    for model in _load(path):
        for chain in model:
            for res in chain:
                if "CA" in res:
                    coords.append(res["CA"].coord)
        break
    return np.array(coords, dtype=float)


def _global_rmsd(path1: Path, path2: Path) -> float:
    C1, C2 = _all_ca(path1), _all_ca(path2)
    n = min(len(C1), len(C2))
    P = C1[:n] - C1[:n].mean(0)
    Q = C2[:n] - C2[:n].mean(0)
    V, _, Wt = np.linalg.svd(P.T @ Q)
    R = V @ np.diag([1., 1., np.sign(np.linalg.det(V @ Wt))]) @ Wt
    d = P - Q @ R.T
    return float(np.sqrt((d * d).sum() / n))


def _superimpose(s_ref, s_mob, pairs):
    a1, a2 = [], []
    for r1, r2 in pairs:
        c1, n1 = r1.split(":")
        c2, n2 = r2.split(":")
        at1 = _get_ca(s_ref, c1, int(n1))
        at2 = _get_ca(s_mob, c2, int(n2))
        if at1 and at2:
            a1.append(at1)
            a2.append(at2)
    if len(a1) < 3:
        raise ValueError(
            f"Only {len(a1)} matched pairs found (need >=3). "
            "Check residue IDs and chain letters."
        )
    sup = Superimposer()
    sup.set_atoms(a1, a2)
    sup.apply(list(s_mob.get_atoms()))
    ref_coords = np.array([a.coord for a in a1])
    mob_coords = np.array([a.coord for a in a2])
    return float(sup.rms), ref_coords, mob_coords


def _radius(coords: np.ndarray) -> float:
    return float(np.linalg.norm(coords - coords.mean(0), axis=1).mean())


def _mean_pairwise(coords: np.ndarray) -> float:
    n = len(coords)
    if n < 2:
        return 0.0
    dists = [np.linalg.norm(coords[i] - coords[j])
             for i in range(n) for j in range(i+1, n)]
    return float(np.mean(dists))


# ── report dataclasses ────────────────────────────────────────────────────────

@dataclass
class HomologReport:
    structure1: str
    structure2: str
    global_rmsd: float
    pocket_rmsd: float
    conservation_ratio: float
    recommendation: str
    reasoning: str


@dataclass
class MutantReport:
    wt_label: str
    mutant_label: str
    global_rmsd: float
    pocket_rmsd: float
    wt_pocket_radius: float
    mut_pocket_radius: float
    radius_change: float
    wt_mean_ca_dist: float
    mut_mean_ca_dist: float
    ca_dist_change: float
    effect: str                 # EXPANDED | CONTRACTED | UNCHANGED
    engineering_impact: str


# ── public API ────────────────────────────────────────────────────────────────

def compare(
    path1: str | Path,
    path2: str | Path,
    pairs: list[tuple[str, str]],
    name1: str = "structure1",
    name2: str = "structure2",
) -> HomologReport:
    """
    Compare two homologous proteins and recommend ACTIVITY or FLUX engineering.

    Parameters
    ----------
    path1, path2  : PDB file paths
    pairs         : equivalent active-site residue pairs
                    e.g. [("A:172","A:173"), ("A:176","A:177")]
    name1, name2  : display labels
    """
    path1, path2 = Path(path1), Path(path2)
    pocket_rmsd, _, _ = _superimpose(_load(path1), _load(path2), pairs)
    global_r = _global_rmsd(path1, path2)
    ratio = round(pocket_rmsd / global_r, 3) if global_r > 0 else 1.0

    if ratio < 0.5:
        rec = "ACTIVITY"
        reason = (
            f"Pocket RMSD ({pocket_rmsd:.2f} A) < half global RMSD ({global_r:.2f} A). "
            "Active site geometry is well conserved. Single-residue mutations "
            "at the binding pocket are likely to improve kcat or Km."
        )
    elif ratio < 0.8:
        rec = "ACTIVITY first, then test FLUX"
        reason = (
            f"Moderate conservation (pocket {pocket_rmsd:.2f} A / global {global_r:.2f} A, "
            f"ratio {ratio}). Start with activity mutations, then measure whole-cell flux."
        )
    else:
        rec = "FLUX"
        reason = (
            f"Pocket RMSD ({pocket_rmsd:.2f} A) ~ global RMSD ({global_r:.2f} A), "
            f"ratio {ratio}. Active site as diverged as fold. Focus on pathway-level "
            "changes: overexpress upstream enzymes, remove competing branches."
        )

    return HomologReport(
        structure1=name1, structure2=name2,
        global_rmsd=round(global_r, 3),
        pocket_rmsd=round(pocket_rmsd, 3),
        conservation_ratio=ratio,
        recommendation=rec, reasoning=reason,
    )


def compare_mutant(
    wt_path: str | Path,
    mut_path: str | Path,
    pocket_res_ids: list[str],
    wt_label: str = "WT",
    mut_label: str = "Mutant",
) -> MutantReport:
    """
    Compare WT vs engineered mutant of the SAME protein.

    Parameters
    ----------
    wt_path, mut_path   : PDB file paths
    pocket_res_ids      : pocket residues in WT numbering
                          e.g. ["A:172", "A:176", "A:225", "A:228"]
                          (assumed same numbering in mutant for point mutants)
    wt_label, mut_label : display labels
    """
    wt_path, mut_path = Path(wt_path), Path(mut_path)
    pairs = [(r, r) for r in pocket_res_ids]
    pocket_rmsd, wt_coords, mut_coords = _superimpose(
        _load(wt_path), _load(mut_path), pairs
    )
    global_r    = _global_rmsd(wt_path, mut_path)
    wt_rad      = _radius(wt_coords)
    mut_rad     = _radius(mut_coords)
    rad_change  = mut_rad - wt_rad
    wt_dist     = _mean_pairwise(wt_coords)
    mut_dist    = _mean_pairwise(mut_coords)
    dist_change = mut_dist - wt_dist

    if abs(rad_change) < 0.3:
        effect = "UNCHANGED"
        impact = (
            f"Pocket geometry essentially unchanged (radius delta {rad_change:+.2f} A). "
            "Mutation unlikely to act via steric change. "
            "Check electrostatics or run MD simulation for dynamic effects."
        )
    elif rad_change > 0:
        effect = "EXPANDED"
        impact = (
            f"Pocket radius increased by {rad_change:+.2f} A, "
            f"mean inter-residue distance changed by {dist_change:+.2f} A. "
            "Binding cavity is larger — reduced steric hindrance improves "
            "substrate access and product release. Good candidate for whole-cell flux assay."
        )
    else:
        effect = "CONTRACTED"
        impact = (
            f"Pocket radius decreased by {abs(rad_change):.2f} A, "
            f"mean inter-residue distance changed by {dist_change:+.2f} A. "
            "Tighter binding cavity — may improve Km (better substrate positioning) "
            "but could restrict product release (lower kcat). "
            "Run in vitro kinetics before whole-cell assay."
        )

    return MutantReport(
        wt_label=wt_label, mutant_label=mut_label,
        global_rmsd=round(global_r, 3),
        pocket_rmsd=round(pocket_rmsd, 3),
        wt_pocket_radius=round(wt_rad, 3),
        mut_pocket_radius=round(mut_rad, 3),
        radius_change=round(rad_change, 3),
        wt_mean_ca_dist=round(wt_dist, 3),
        mut_mean_ca_dist=round(mut_dist, 3),
        ca_dist_change=round(dist_change, 3),
        effect=effect,
        engineering_impact=impact,
    )
