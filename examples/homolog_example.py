"""
Example: compare two terpene synthase homologs (2ZCP vs 5IYS)
and get an ACTIVITY vs FLUX recommendation.
"""
from pocketcompare import fetch, compare

p1 = fetch.fetch_pdb("2ZCP")
p2 = fetch.fetch_pdb("5IYS")

# DXXXD motif residues — equivalent active-site aspartates
pairs = [
    ("A:172", "A:173"),
    ("A:176", "A:177"),
    ("A:225", "A:237"),
]

report = compare.compare(p1, p2, pairs, name1="CrtM (2ZCP)", name2="PhySyn (5IYS)")

print(f"Global RMSD         : {report.global_rmsd} A")
print(f"Pocket RMSD         : {report.pocket_rmsd} A")
print(f"Conservation ratio  : {report.conservation_ratio}")
print(f"Recommendation      : {report.recommendation}")
print(f"Reasoning           : {report.reasoning}")
