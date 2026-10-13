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


