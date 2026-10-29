"""Per-residue solvent accessibility over the binding footprint.

Run on the UNBOUND receptor. The point is to establish that the footprint
residues were solvent-exposed before any ligand was introduced, independently of
the docking result — a docked pose at a buried site would be an artifact, and
this is what rules that out.

Optionally takes the buried surface area of the complex and reports what
fraction of the available footprint it covers.

    python sasa_footprint.py --receptor 2N0A.pdb --chain A \
        --residues 68-72,75,91-95 --bsa 138
"""
import argparse

from Bio.PDB import PDBParser, ShrakeRupley


def parse_residues(spec):
    """'68-72,75,91-95' -> [68, 69, 70, 71, 72, 75, 91, 92, 93, 94, 95]"""
    out = []
    for part in spec.split(","):
        if "-" in part:
            lo, hi = part.split("-")
            out.extend(range(int(lo), int(hi) + 1))
        else:
            out.append(int(part))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", required=True, help="unbound receptor PDB")
    ap.add_argument("--chain", default="A")
    ap.add_argument("--residues", required=True,
                    help="footprint, e.g. 68-72,75,91-95")
    ap.add_argument("--bsa", type=float,
                    help="buried surface area of the complex, A^2")
    ap.add_argument("--model", type=int, default=0,
                    help="NMR model index; 2N0A is an ensemble")
    a = ap.parse_args()

    footprint = parse_residues(a.residues)

    model = PDBParser(QUIET=True).get_structure("rec", a.receptor)[a.model]
    ShrakeRupley().compute(model, level="R")

    total = 0.0
    print(f"{'resi':<6}{'resn':<6}{'SASA':>8}")
    for res in model[a.chain]:
        if res.id[1] in footprint:
            print(f"{res.id[1]:<6}{res.get_resname():<6}{res.sasa:>8.1f}")
            total += res.sasa

    print(f"\ntotal exposed footprint area = {total:.1f} A^2")
    if a.bsa:
        print(f"complex buries {a.bsa:.0f} A^2, "
              f"{a.bsa / total * 100:.0f}% of the available footprint")


if __name__ == "__main__":
    main()
