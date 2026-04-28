from __future__ import annotations
import argparse
from pathlib import Path
from . import fetch
from . import compare as cmp


def _resolve(src: str):
    if len(src) == 4 and src.isalnum():
        return fetch.fetch_pdb(src), src.upper()
    return fetch.load_local(src), Path(src).stem


def _parse_pairs(s: str):
    result = []
    for entry in s.split(";"):
        parts = [x.strip() for x in entry.split(",")]
        if len(parts) == 2:
            result.append((parts[0], parts[1]))
    return result


def _parse_residues(s: str):
    return [r.strip() for r in s.split(";") if r.strip()]


def _fmt_homolog(r) -> str:
    W = 58
    sep = "=" * W
    div = "-" * W
    return "\n".join([
        "",
        sep,
        "  PocketCompare - Homolog Comparison",
        sep,
        "  Structure 1         : " + r.structure1,
        "  Structure 2         : " + r.structure2,
        div,
        "  Global RMSD         : " + str(r.global_rmsd) + " A",
        "  Pocket RMSD         : " + str(r.pocket_rmsd) + " A",
        "  Conservation ratio  : " + str(r.conservation_ratio),
        div,
        "  RECOMMENDATION      : " + r.recommendation,
        div,
        "  " + r.reasoning,
        sep,
        "",
    ])


def _fmt_mutant(r) -> str:
    W = 58
    sep = "=" * W
    div = "-" * W
    arrow = "UP" if r.radius_change > 0 else ("DOWN" if r.radius_change < 0 else "=")
    r_line = "  Pocket radius MUT   : " + str(r.mut_pocket_radius) + " A  [" + arrow + " " + format(r.radius_change, "+.3f") + " A]"
    d_line = "  Mean CA dist MUT    : " + str(r.mut_mean_ca_dist) + " A  [" + arrow + " " + format(r.ca_dist_change, "+.3f") + " A]"
    return "\n".join([
        "",
        sep,
        "  PocketCompare - WT vs Mutant Analysis",
        sep,
        "  WT                  : " + r.wt_label,
        "  Mutant              : " + r.mutant_label,
        div,
        "  Global RMSD         : " + str(r.global_rmsd) + " A  (low = fold preserved)",
        "  Pocket RMSD         : " + str(r.pocket_rmsd) + " A",
        div,
        "  Pocket radius WT    : " + str(r.wt_pocket_radius) + " A",
        r_line,
        "  Mean CA dist WT     : " + str(r.wt_mean_ca_dist) + " A",
        d_line,
        div,
        "  POCKET EFFECT       : " + r.effect,
        div,
        "  " + r.engineering_impact,
        sep,
        "",
    ])


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="pocketcompare",
        description="Compare binding pockets for protein engineering decisions.",
    )
    sub = p.add_subparsers(dest="mode", required=True)

    hom = sub.add_parser("homolog", help="Compare two different proteins")
    hom.add_argument("--pdb1",     help="PDB ID or local file, structure 1")
    hom.add_argument("--pdb2",     help="PDB ID or local file, structure 2")
    hom.add_argument("--uniprot1", help="UniProt ID (AlphaFold) structure 1")
    hom.add_argument("--uniprot2", help="UniProt ID (AlphaFold) structure 2")
    hom.add_argument("--pairs", required=True, help="Equivalent residue pairs: A:172,A:173;A:176,A:177")
    hom.add_argument("--out", help="Save report to text file")

    mut = sub.add_parser("mutant", help="Compare WT vs mutant of same protein")
    mut.add_argument("--wt",  required=True, help="WT: PDB ID or local file")
    mut.add_argument("--mut", required=True, help="Mutant: PDB ID or local file")
    mut.add_argument("--pocket", required=True, help="Pocket residues: A:172;A:176;A:225;A:228")
    mut.add_argument("--wt-label",  default="WT",     dest="wt_label")
    mut.add_argument("--mut-label", default="Mutant", dest="mut_label")
    mut.add_argument("--out", help="Save report to text file")

    return p


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.mode == "homolog":
        src1 = args.pdb1 or args.uniprot1
        src2 = args.pdb2 or args.uniprot2
        if not src1 or not src2:
            parser.error("Provide --pdb1/--uniprot1 and --pdb2/--uniprot2")
        if args.uniprot1:
            path1, lbl1 = fetch.fetch_alphafold(args.uniprot1), "AF-" + args.uniprot1
        else:
            path1, lbl1 = _resolve(src1)
        if args.uniprot2:
            path2, lbl2 = fetch.fetch_alphafold(args.uniprot2), "AF-" + args.uniprot2
        else:
            path2, lbl2 = _resolve(src2)
        pairs = _parse_pairs(args.pairs)
        if not pairs:
            parser.error("No valid pairs. Format: A:172,A:173;A:176,A:177")
        report = cmp.compare(path1, path2, pairs, name1=lbl1, name2=lbl2)
        output = _fmt_homolog(report)

    else:
        path_wt,  lbl_wt  = _resolve(args.wt)
        path_mut, lbl_mut = _resolve(args.mut)
        pocket = _parse_residues(args.pocket)
        if not pocket:
            parser.error("No valid residues. Format: A:172;A:176;A:225")
        report = cmp.compare_mutant(
            path_wt, path_mut, pocket,
            wt_label=args.wt_label,
            mut_label=args.mut_label,
        )
        output = _fmt_mutant(report)

    print(output)
    if args.out:
        Path(args.out).write_text(output)
        print("Saved to " + args.out)


if __name__ == "__main__":
    main()
