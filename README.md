# LEADmdl
A de novo neutral iron(III) porphyrin designed to bind the ⁶⁸GAVVT⁷² steric zipper at the α-synuclein fibril growth terminus.
Every characterized metalloporphyrin that binds amyloidogenic protein carries a net charge, which either bars CNS entry or directs binding away from the aggregation interface. LEADmdl is net-neutral by design: a desymmetrized two-winged macrocycle with a perpendicular axial hydroxide on Fe(III). The target is the zipper interface monomers pair through — buried throughout a mature fibril, exposed only at the growth terminus, which means occupancy scales with fibril number rather than total protein.
This repository contains the full computational pipeline: DFT electronic structure, interfacial docking, interaction analysis, descriptor profiling, and aggregation kinetics.
Nothing here has been synthesized. All results are computational predictions.
**Key results**
Quantity	Value

Formula / MW	C₃₃H₁₇F₃FeN₅O, 612.37 Da
Net charge	0
Fe(III) ground state	low-spin doublet (sextet 0.1154 Ha higher)
Fe(II) ground state	S = 1 triplet (quintet 38.1 kcal/mol higher)
α / β HOMO–LUMO gap	2.529 eV / 2.393 eV
Net dipole	5.84 D
Reduction potential	−1.10 V (ε = 78.4), −1.44 V (ε = 4)
Docking ΔG	−7.95 kcal/mol (control: −4.7)
Contacts at terminus	12, terminal chain
Val70 bidentate clip	donates 2.77 Å, accepts 3.15 Å
Predicted logBB	+0.17
Fibril mass at 10 μM	−88.4% by 37.5 h
