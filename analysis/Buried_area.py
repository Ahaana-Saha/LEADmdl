"""Buried surface area of the docked complex.

BSA = (SASA_receptor + SASA_ligand - SASA_complex) / 2

The halving is because the formula counts the interface once from each side.
Both components are measured in their docked coordinate frames and simply
concatenated — no superposition, because the pose is already in the receptor's
frame.

The PDBQT is converted to PDB by truncating at column 66, which drops the
AutoDock atom type and partial charge. Biopython then infers elements from atom
names; check the iron is read as Fe and not something else, since a wrong
element means a wrong radius.

    python buried_area.py --receptor 2N0A.pdb --ligand webina_out.pdbqt
"""
import argparse
import os
import tempfile

from Bio.PDB import PDBParser, ShrakeRupley


def pdbqt_to_pdb(src, dst):
    """First pose only, as HETATM records truncated to standard PDB width."""
    n = 0
    with open(dst, "w") as out:
        for line in open(src):
            if line.startswith(("ATOM", "HETATM")):
                out.write("HETATM" + line[6:66].rstrip() + "\n")
                n += 1
            elif line.startswith("ENDMDL"):
                break
        out.write("END\n")
    return n


def first_model_to_pdb(src, dst):
    """Strip an NMR ensemble down to its first model."""
    with open(dst, "w") as out:
        for line in open(src):
            if line.startswith(("ATOM", "HETATM")):
                out.write(line)
            elif line.startswith("ENDMDL"):
                break
        out.write("END\n")


def total_sasa(parser, path, name):
    model = parser.get_structure(name, path)[0]
    ShrakeRupley().compute(model, level="A")
    return sum(a.sasa for a in model.get_atoms())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--ligand", required=True, help="docked .pdbqt")
    a = ap.parse_args()

    parser = PDBParser(QUIET=True)
    tmp = tempfile.mkdtemp()
    lig_pdb = os.path.join(tmp, "lig.pdb")
    rec_pdb = os.path.join(tmp, "rec.pdb")
    cpx_pdb = os.path.join(tmp, "complex.pdb")

    n = pdbqt_to_pdb(a.ligand, lig_pdb)
    print(f"ligand atoms: {n}")
    first_model_to_pdb(a.receptor, rec_pdb)

    with open(cpx_pdb, "w") as c:
        for path in (rec_pdb, lig_pdb):
            for line in open(path):
                if line.startswith(("ATOM", "HETATM")):
                    c.write(line)
        c.write("END\n")

    sasa_rec = total_sasa(parser, rec_pdb, "rec")
    sasa_lig = total_sasa(parser, lig_pdb, "lig")
    sasa_cpx = total_sasa(parser, cpx_pdb, "cpx")

    print(f"\nreceptor alone : {sasa_rec:10.1f} A^2")
    print(f"ligand alone   : {sasa_lig:10.1f} A^2")
    print(f"complex        : {sasa_cpx:10.1f} A^2")
    print(f"\nburied surface area = "
          f"{(sasa_rec + sasa_lig - sasa_cpx) / 2:.1f} A^2")


if __name__ == "__main__":
    main()
