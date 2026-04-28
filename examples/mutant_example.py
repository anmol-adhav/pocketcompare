"""
Example: compare WT (2ZCP) vs a mutant model.

Replace mutant_model.pdb with your AlphaFold/Rosetta mutant structure.
The pocket residues should be the key active-site residues you mutated around.
"""
from pocketcompare import fetch, compare
from pathlib import Path

wt_path = fetch.fetch_pdb("2ZCP")

# Replace this with your actual mutant model path
mutant_path = Path("mutant_model.pdb")
if not mutant_path.exists():
    print("No mutant_model.pdb found — using WT as stand-in for demonstration.")
    mutant_path = wt_path

pocket_residues = ["A:172", "A:176", "A:225", "A:228"]

report = compare.compare_mutant(
    wt_path, mutant_path,
    pocket_res_ids=pocket_residues,
    wt_label="CrtM_WT",
    mut_label="CrtM_Mutant",
)

print(f"Global RMSD        : {report.global_rmsd} A")
print(f"Pocket RMSD        : {report.pocket_rmsd} A")
print(f"Pocket radius WT   : {report.wt_pocket_radius} A")
print(f"Pocket radius MUT  : {report.mut_pocket_radius} A ({report.radius_change:+.3f})")
print(f"Effect             : {report.effect}")
print(f"Impact             : {report.engineering_impact}")
