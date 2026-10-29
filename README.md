# LEADmdl

A de novo neutral iron(III) porphyrin designed to bind the ⁶⁸GAVVT⁷² steric zipper at the
α-synuclein fibril growth terminus.

Every characterized metalloporphyrin that binds amyloidogenic protein carries a net charge,
which either bars CNS entry or directs binding away from the aggregation interface. LEADmdl
is net-neutral by design: a desymmetrized two-winged macrocycle with a perpendicular axial
hydroxide on Fe(III). The target is the zipper interface monomers pair through — buried
throughout a mature fibril, exposed only at the growth terminus, which means occupancy
scales with fibril number rather than total protein.

**Nothing here has been synthesized.** All results are computational predictions.

---

## Key results

| Quantity | Value |
|---|---|
| Formula / MW | C₃₃H₁₇F₃FeN₅O, 612.37 Da |
| Net charge | 0 |
| Fe(III) ground state | low-spin doublet (sextet 0.1154 Ha higher) |
| Fe(II) ground state | S = 1 triplet (quintet 38.1 kcal/mol higher) |
| α / β HOMO–LUMO gap | 2.529 eV / 2.393 eV |
| Net dipole | 5.84 D |
| Reduction potential | −1.10 V (ε = 78.4), −1.44 V (ε = 4) |
| Docking ΔG | −7.95 kcal/mol (positional control: −4.7) |
| Contacting residues at terminus | 12, single chain |
| Val70 bidentate clip | donates 2.77 Å, accepts 3.15 Å |
| Predicted logBB | +0.17 |
| Fibril mass at 10 μM | −88.4% by 37.5 h |

---

## Repository structure

```
structures/    input geometry
dft/           PySCF / gpu4pyscf run scripts, one per spin and dielectric
docking/       metal-aware ligand preparation
analysis/      contact mapping, solvent accessibility, buried surface area
kinetics/      Finke–Watzky capping model, Cohen–Knowles–Buell half-times
```

PLIP interaction typing and the PyMOL coordinate measurements were run with standard tools and are not scripted here.

Solvent accessibility, PLIP interaction typing, and the PyMOL coordinate measurements were run with standard tools and are not scripted here.

---

## Metal-aware ligand preparation

Automated PDBQT preparation discards the iron outright, leaving the pyrrole nitrogens
unbonded and the axial hydroxide floating free. Docking under those conditions does not
converge to a consistent site across independent runs — and the tool reports no error.

`docking/quantum_prep.py` routes preparation through mol2, which preserves the Fe–N
coordination, and writes ab initio Mulliken charges from the converged DFT density in place
of Gasteiger defaults. Any docking study of a metalloporphyrin hits this failure; the fix
transfers.

---

## Descriptors

Computed with standard tools — Ertl TPSA, Wildman–Crippen logP, Clark logBB regression — so
no wrapper is included here. The one non-standard step is a chelation correction to polar
surface area: descriptor models are parameterized on neutral organics and treat the four
macrocycle nitrogens as hydrogen-bond capable, which they are not once coordinated to iron.
Subtracting their contribution — 2 × 15.79 Å² (pyrrole NH) + 2 × 12.89 Å² (aromatic N) =
57.36 Å² — gives a corrected TPSA of 44.02 Å² from a base value of 101.38 Å². Experimental
logPow was substituted for the permanently ionized benchmarks, where computed values are
not meaningful.

---

## Environment

```
python            3.11
pyscf             2.7
gpu4pyscf-cuda12x 
biopython         1.88
numpy
```

Open Babel 3.x, AutoDock Vina 1.2.x, PLIP 2.x, PyMOL 3.0. DFT ran on a single NVIDIA T4 under CUDA 13.4.

---

## Reproducing the DFT

All runs use UB3LYP with def2-SVP on Fe and 6-31G* elsewhere, density fitting, IEF-PCM,
level shift 0.4 and damping 0.5, conv_tol 1e-5. Each script converges with those aids on,
then releases them and re-converges from the checkpoint — the released energy is the
reported one.

```
python dft/gas_fe3.py                              # neutral, spin=1
python dft/fe3_pcm.py --eps 78.3553
python dft/fe3_pcm.py --eps 4.0
python dft/fe3_pcm.py --eps 78.3553 --spin 5 --seed ""   # sextet check
python dft/fe2_scan.py --eps 78.3553               # spin = 0, 2, 4
python dft/fe2_pcm.py --eps 4.0 --spin 2
```

Spin-state assignments were verified, not assumed. ⟨S²⟩ is reported for every run, and
`dft/README.md` carries the full energy table.

---

## Known limitations

- All DFT results are **single-point** evaluations at the MMFF94 pre-optimized geometry.
  The structure was not relaxed at the DFT level. Geometry optimization at GFN2-xTB is the
  next step.
- AutoDock Vina's scoring function is not parameterized for transition metals; the pose was
  validated by contact mapping, solvent accessibility and interaction typing rather than by
  score alone.
- Descriptor models are parameterized on neutral organics. The chelation correction above is
  chemically motivated and uniformly applied, but not independently validated.
- Absolute reduction potentials are not quantitative: no diffuse functions, no relaxation of
  the reduced state, and the computed product state retains the axial hydroxide. The
  inter-environment difference (0.34 V) is the more reliable quantity.
- The fibril mass simulation uses apparent first-order rates rather than the literal
  bimolecular elongation constant, which saturates the system before 37.5 h. The capping
  ratio carries the inhibitor effect; the absolute rates set the timescale. See
  `kinetics/fw_capping.py`.
- Molecular dynamics has not yet been run; the docking result assumes a rigid receptor.
- Ligand self-association is unassessed.

---

## Contact

<name> — <email>
