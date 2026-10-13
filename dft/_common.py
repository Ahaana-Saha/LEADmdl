"""Shared setup for LEADmdl DFT runs.

All single-point UB3LYP, def2-SVP on Fe and 6-31G* elsewhere, density fitting.
Open-shell iron needs convergence aids to reach an SCF solution at all, but a
level shift biases the result toward the starting occupation. Every run therefore
converges with the aids on, then releases them and re-converges from the
checkpoint. The released energy is the reported one.
"""
import os
import sys

from pyscf import gto
from gpu4pyscf import dft as dft_mod

BASIS = {"Fe": "def2-svp", "default": "6-31g*"}
N_ATOMS = 60
MAX_MEMORY = 12000          # MB; keep near 70% of available RAM


def build_mol(xyz, charge, spin):
    """spin is the number of UNPAIRED electrons (PySCF convention), not 2S+1."""
    mol = gto.M(atom=xyz, charge=charge, spin=spin, basis=BASIS,
                max_memory=MAX_MEMORY, verbose=4)
    assert mol.natm == N_ATOMS, f"expected {N_ATOMS} atoms, got {mol.natm}"
    print("natm", mol.natm, "nelec", mol.nelec)
    return mol


def make_mf(mol, chkfile, eps=None):
    mf = dft_mod.UKS(mol).density_fit().to_gpu()
    if eps is not None:
        mf = mf.PCM()
        mf.with_solvent.method = "IEF-PCM"
        mf.with_solvent.eps = eps
    # verbose/stdout must be set AFTER to_gpu() and PCM() — both return new objects
    mf.verbose = 4
    mf.stdout = sys.stdout
    mf.xc = "b3lyp"
    mf.chkfile = chkfile
    mf.conv_tol = 1e-5
    mf.max_cycle = 300
    return mf


def run(mf, chkfile, seed=None, label="E"):
    """Two-stage convergence: aids on, then released."""
    mf.level_shift, mf.damp = 0.4, 0.5

    dm = None
    src = chkfile if os.path.exists(chkfile) else seed
    if src and os.path.exists(src):
        try:
            dm = mf.from_chk(src)
            print("seeded from", os.path.basename(src))
        except Exception as exc:                     # stale or mismatched chk
            print("fresh start:", type(exc).__name__)

    mf.kernel(dm0=dm)
    print("stage 1 converged:", mf.converged)

    mf.level_shift, mf.damp = 0.0, 0.0
    energy = mf.kernel(dm0=mf.from_chk(chkfile))

    print("converged:", mf.converged)
    print(f"{label} =", energy)
    print("spin_square:", mf.to_cpu().spin_square())
    return energy

"""Fe(II) at one spin state and dielectric. Use after fe2_scan has established
which state is the ground state; seed from the same spin at the other dielectric.
"""
import argparse

from _common import build_mol, make_mf, run

p = argparse.ArgumentParser()
p.add_argument("--eps", type=float, required=True)
p.add_argument("--spin", type=int, default=2)        # triplet ground state
p.add_argument("--xyz", default="../structures/leadmdl_fixed.xyz")
p.add_argument("--seed", default="fe2_eps78.3553_s2.chk")
p.add_argument("--out")
a = p.parse_args()

out = a.out or f"fe2_eps{a.eps:g}_s{a.spin}.chk"
mol = build_mol(a.xyz, charge=-1, spin=a.spin)
mf = make_mf(mol, out, eps=a.eps)
run(mf, out, seed=a.seed, label=f"E_fe2_eps{a.eps:g}_s{a.spin}")

"""Fe(II) spin scan at one dielectric: S = 0, 1, 2 (spin = 0, 2, 4).

The S = 0 solution comes back broken-symmetry rather than closed-shell
(<S^2> ~ 0.9 against an ideal 0.0), so the triplet is the lowest pure-spin
state. Check the printed <S^2> for every state before using any of these
energies — that is the point of the scan.
"""
import argparse

from _common import build_mol, make_mf, run

p = argparse.ArgumentParser()
p.add_argument("--eps", type=float, required=True)
p.add_argument("--xyz", default="../structures/leadmdl_fixed.xyz")
p.add_argument("--seed", default="gas_fe3.chk")
a = p.parse_args()

results = {}
seed = a.seed
for spin in (0, 2, 4):
    out = f"fe2_eps{a.eps:g}_s{spin}.chk"
    mol = build_mol(a.xyz, charge=-1, spin=spin)
    mf = make_mf(mol, out, eps=a.eps)
    results[spin] = run(mf, out, seed=seed, label=f"E_fe2_eps{a.eps:g}_s{spin}")
    seed = out                      # each state seeds the next

HA2KCAL = 627.5095
ref = results[2]                    # triplet
print("\nrelative to the S = 1 triplet, kcal/mol")
for spin, energy in sorted(results.items()):
    print(f"  spin={spin}: {(energy - ref) * HA2KCAL:+.2f}")

"""Fe(III) under IEF-PCM. --spin 1 is the doublet ground state; --spin 5 is the
sextet, run to verify the low-spin assignment rather than assume it.

The sextet does not converge from a doublet density — omit --seed for it.
"""
import argparse

from _common import build_mol, make_mf, run

p = argparse.ArgumentParser()
p.add_argument("--eps", type=float, required=True)   # 78.3553 water, 4.0 buried
p.add_argument("--spin", type=int, default=1)        # unpaired electrons
p.add_argument("--xyz", default="../structures/leadmdl_fixed.xyz")
p.add_argument("--seed", default="gas_fe3.chk")
p.add_argument("--out")
a = p.parse_args()

out = a.out or f"fe3_eps{a.eps:g}_s{a.spin}.chk"
mol = build_mol(a.xyz, charge=0, spin=a.spin)
mf = make_mf(mol, out, eps=a.eps)
run(mf, out, seed=a.seed, label=f"E_fe3_eps{a.eps:g}_s{a.spin}")

"""Gas-phase Fe(III), low-spin doublet. Run this first — it seeds everything else."""
import argparse

from _common import build_mol, make_mf, run

p = argparse.ArgumentParser()
p.add_argument("--xyz", default="../structures/leadmdl_fixed.xyz")
p.add_argument("--out", default="gas_fe3.chk")
a = p.parse_args()

mol = build_mol(a.xyz, charge=0, spin=1)        # doublet, S = 1/2
mf = make_mf(mol, a.out)
run(mf, a.out, label="E_gas_fe3")

"""Reduction potential from the two single-point energies.

V = [E(Fe3+) - E(Fe2+)] * 27.2114 - 4.44

The absolute values are not quantitative: 6-31G* carries no diffuse functions
and so underbinds the added electron on the anion, the reduced state is not
geometry-relaxed, and the computed product retains the axial hydroxide where
experiment is expected to release it. These errors act in the same direction at
both dielectrics and largely cancel in the difference, so the inter-environment
shift is the more reliable quantity.
"""
import argparse

HARTREE_TO_EV = 27.2114
SHE_ABSOLUTE = 4.44

p = argparse.ArgumentParser()
p.add_argument("--e-ox", type=float, required=True)     # Fe(III), Hartree
p.add_argument("--e-red", type=float, required=True)    # Fe(II), Hartree
a = p.parse_args()

delta = a.e_ox - a.e_red
print(f"dE   = {delta:.7f} Ha = {delta * HARTREE_TO_EV:.3f} eV")
print(f"V    = {delta * HARTREE_TO_EV - SHE_ABSOLUTE:.2f} V vs SHE")

# DFT

Single-point UB3LYP, def2-SVP on Fe and 6-31G* elsewhere, density fitting,
IEF-PCM. Geometry is the MMFF94 pre-optimized structure — **not** relaxed at the
DFT level.

Note that PySCF's `spin` is the number of unpaired electrons, not the
multiplicity: `spin=1` is a doublet, `spin=4` a quintet.

## Order

```bash
python gas_fe3.py                                  # seeds everything below

python fe3_pcm.py --eps 78.3553                    # Fe(III) doublet, water
python fe3_pcm.py --eps 4.0
python fe3_pcm.py --eps 78.3553 --spin 5 --seed "" # sextet check, no seed

python fe2_scan.py --eps 78.3553                   # S = 0, 1, 2
python fe2_pcm.py --eps 4.0 --spin 2               # triplet, buried interface

python potential.py --e-ox <E> --e-red <E>
```

The sextet will not converge from a doublet density — pass `--seed ""` and let it
start from the atomic guess. It takes far longer than the seeded runs.

## Results

| state | ε = 78.4 | ε = 4 | ⟨S²⟩ (ideal) |
|---|---|---|---|
| Fe(III) S = ½ | −3179.2554565 | −3179.2417665 | — |
| Fe(III) S = 5/2 | −3179.1400709 | — | 9.04 (8.75) |
| Fe(II) S = 0 (BS) | −3179.3833322 | — | 0.91 (0.00) |
| Fe(II) S = 1 | −3179.3781417 | −3179.3518833 | 2.19 (2.00) |
| Fe(II) S = 2 | −3179.3174093 | — | 6.23 (6.00) |

Ground states: Fe(III) doublet, with the sextet 0.1154 Ha higher. Fe(II) triplet,
with the quintet 38.1 kcal/mol higher. The S = 0 solution is broken-symmetry, not
a closed-shell singlet — see `yamaguchi.py`.

Potentials: **−1.10 V** at ε = 78.4, **−1.44 V** at ε = 4. Yamaguchi-projected
singlet gives −0.86 V, bracketing the aqueous value to 0.24 V.

Gas-phase Fe(III): −3179.2193614 Ha, dipole 5.84 D, α gap 2.529 eV, β gap 2.393 eV.


"""Yamaguchi spin projection for the broken-symmetry Fe(II) singlet.

A UKS calculation asked for S = 0 returns an antiferromagnetically coupled
solution with <S^2> near 1, not a closed-shell singlet. Its energy is
variationally below the triplet but is not a spin eigenstate. This projects a
pure singlet energy from the BS/triplet pair.

Yamaguchi, Takahara & Fueno (1986), in Applied Quantum Chemistry, 155-184.
"""
import argparse

p = argparse.ArgumentParser()
p.add_argument("--e-bs", type=float, required=True)
p.add_argument("--s2-bs", type=float, required=True)
p.add_argument("--e-triplet", type=float, required=True)
p.add_argument("--s2-triplet", type=float, required=True)
a = p.parse_args()

ratio = a.s2_bs / (a.s2_triplet - a.s2_bs)
e_singlet = a.e_bs - ratio * (a.e_triplet - a.e_bs)
print(f"projected singlet = {e_singlet:.7f} Ha")
