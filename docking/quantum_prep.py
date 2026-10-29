"""Metal-aware ligand preparation for AutoDock Vina.
 
Automated PDBQT preparation discards the iron outright, leaving the pyrrole
nitrogens unbonded and the axial hydroxide floating free. Docking under those
conditions does not converge to a consistent site across independent runs, and
the tool reports no error.

Routing through mol2 preserves the Fe-N coordination. This script writes ab
initio Mulliken charges from a converged DFT density over the mol2's Gasteiger
defaults, then converts to PDBQT with --partialcharge none so Open Babel does
not silently recompute them on write.

Verify the output before trusting it: the iron must appear inside the ROOT
block alongside the four pyrrole nitrogens.

    python quantum_prep.py --mol2 ligand.mol2 --chk gas_fe3.chk -o ligand.pdbqt
"""
import argparse
import shutil
import subprocess
import sys

import pyscf
from pyscf import dft

# Must match the basis the checkpoint was computed in. from_chk() rebuilds the
# density from the checkpoint's MO coefficients, which live in that basis — a
# mismatch gives a dimension error at best and wrong charges at worst.
BASIS = {"Fe": "def2-svp", "default": "6-31g*"}


def parse_mol2_atoms(path):
    """Return (element, x, y, z) per atom, taking the element from the SYBYL
    type column rather than the atom name: 'C.ar' -> 'C'. Atom names are free
    text and a name like 'HN2' would parse as element 'HN'."""
    atoms = []
    in_block = False
    for line in open(path):
        if line.startswith("@<TRIPOS>ATOM"):
            in_block = True
            continue
        if line.startswith("@<TRIPOS>"):
            in_block = False
            continue
        if in_block and line.strip():
            p = line.split()
            if len(p) >= 6:
                atoms.append((p[5].split(".")[0], p[2], p[3], p[4]))
    return atoms


def mulliken_charges(atoms, chk, charge, spin):
    spec = "; ".join(f"{el} {x} {y} {z}" for el, x, y, z in atoms)
    mol = pyscf.M(atom=spec, basis=BASIS, charge=charge, spin=spin)
    mf = dft.UKS(mol)
    mf.xc = "b3lyp"
    _, charges = mf.mulliken_pop(dm=mf.from_chk(chk))
    return charges


def write_charges(src, dst, charges):
    """Replace column 9 (partial charge) of each ATOM record."""
    out, in_block, i = [], False, 0
    for line in open(src):
        if line.startswith("@<TRIPOS>ATOM"):
            in_block = True
            out.append(line)
            continue
        if line.startswith("@<TRIPOS>"):
            in_block = False
        if in_block and line.strip():
            p = line.split()
            if len(p) >= 9 and i < len(charges):
                out.append(
                    "%7s %-8s %10.4f %10.4f %10.4f %-6s %4s %-8s %10.6f\n"
                    % (p[0], p[1], float(p[2]), float(p[3]), float(p[4]),
                       p[5], p[6], p[7], charges[i])
                )
                i += 1
                continue
        out.append(line)
    open(dst, "w").writelines(out)
    return i


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mol2", required=True)
    ap.add_argument("--chk", required=True, help="converged DFT checkpoint")
    ap.add_argument("-o", "--out", required=True, help="output .pdbqt")
    ap.add_argument("--charge", type=int, default=0)
    ap.add_argument("--spin", type=int, default=1,
                    help="unpaired electrons (PySCF convention), not 2S+1")
    a = ap.parse_args()

    if shutil.which("obabel") is None:
        sys.exit("obabel not found on PATH")

    atoms = parse_mol2_atoms(a.mol2)
    print(f"{len(atoms)} atoms in {a.mol2}")

    charges = mulliken_charges(atoms, a.chk, a.charge, a.spin)
    if len(charges) != len(atoms):
        sys.exit(f"charge/atom mismatch: {len(charges)} vs {len(atoms)} — "
                 "is the checkpoint from this structure and basis?")

    tmp = a.out.rsplit(".", 1)[0] + "_quantum.mol2"
    n = write_charges(a.mol2, tmp, charges)
    print(f"wrote {n} Mulliken charges to {tmp}")

    # --partialcharge none: without it Open Babel recomputes Gasteiger on write
    # and discards everything above.
    subprocess.run(
        ["obabel", tmp, "-O", a.out, "--partialcharge", "none"], check=True
    )
    print(f"wrote {a.out}")
    print("verify: iron present, inside ROOT with the four pyrrole nitrogens")


if __name__ == "__main__":
    main()
