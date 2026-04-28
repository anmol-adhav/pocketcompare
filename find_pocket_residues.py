"""
Utility: print all ASP and GLU residues in a structure.
Useful for identifying DXXXD motif residues to use as --pocket or --pairs input.

Usage: python find_pocket_residues.py 2ZCP
       python find_pocket_residues.py my_model.pdb
"""
import sys
from pathlib import Path
from Bio.PDB import PDBParser
from pocketcompare import fetch

def find_charged_residues(pdb_id_or_path: str, resnames=("ASP", "GLU", "HIS", "LYS", "ARG")):
    src = pdb_id_or_path
    if len(src) == 4 and src.isalnum():
        path = fetch.fetch_pdb(src)
    else:
        path = Path(src)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("x", str(path))

    print(f"\nCharged residues in {src}:")
    print(f"{'Chain:ResNum':<14} {'ResName':<8}")
    print("-" * 22)
    for model in structure:
        for chain in model:
            for res in chain:
                if res.resname in resnames and "CA" in res:
                    print(f"{chain.id}:{res.id[1]:<12} {res.resname}")
        break


if __name__ == "__main__":
    src = sys.argv[1] if len(sys.argv) > 1 else "2ZCP"
    find_charged_residues(src)
