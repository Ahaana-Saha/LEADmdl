"""Interfacial contact mapping between a docked ligand and a fibril receptor.

Reports every receptor residue with at least one atom within the cutoff of any
ligand atom — contacting residues, not atom-atom contacts. The chain summary is
the point: a pose localized to one chain is terminal, a pose bridging chains is
not.

Both files are used in their own coordinate frames. The docked pose comes out of
Vina in the receptor's frame already, so no superposition is applied and none
should be.

    python contact_map.py --receptor 2N0A.pdb --ligand webina_out.pdbqt
"""
import argparse

import numpy as np
from Bio.PDB import PDBParser


def read_ligand(path):
    """First pose only — stop at the first ENDMDL."""
    coords = []
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")):
            coords.append([float(line[30:38]),
                           float(line[38:46]),
                           float(line[46:54])])
        elif line.startswith("ENDMDL"):
            break
    return np.array(coords)


def contacting_residues(model, lig, cutoff):
    hits = []
    for chain in model:
        for res in chain:
            for atom in res:
                if np.min(np.linalg.norm(lig - atom.coord, axis=1)) < cutoff:
                    hits.append((chain.id, res.id[1], res.get_resname()))
                    break
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--ligand", required=True)
    ap.add_argument("--cutoff", type=float, default=4.5)
    ap.add_argument("--model", type=int, default=0,
                    help="NMR model index; 2N0A is an ensemble, first by default")
    a = ap.parse_args()

    lig = read_ligand(a.ligand)
    if lig.size == 0:
        raise SystemExit(f"no atoms read from {a.ligand}")
    print(f"ligand atoms read: {len(lig)}\n")

    model = PDBParser(QUIET=True).get_structure("rec", a.receptor)[a.model]
    hits = contacting_residues(model, lig, a.cutoff)

    print(f"{'chain':<8}{'resi':<8}{'resn':<8}")
    for chain, resi, resn in sorted(hits):
        print(f"{chain:<8}{resi:<8}{resn:<8}")

    chains = sorted({c for c, _, _ in hits})
    print(f"\ncontacting residues: {len(hits)}  (cutoff {a.cutoff} A)")
    print(f"chains involved: {chains}")
    print(">> bridges multiple chains" if len(chains) > 1 else ">> single chain")


if __name__ == "__main__":
    main()
