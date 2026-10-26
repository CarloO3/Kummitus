# Kummitus — DOF from Poles

![Kummitus](Kummitus.png)

**Kummitus** is a Mathematica package for computing the physical degrees of freedom of higher-spin field theories directly from the poles of the saturated propagator. It is built on top of the [xAct](http://www.xact.es/) tensor algebra system (xTensor, xPert, xCoba, xTras).

Copyright © 2026, Carlo Marzo. Released under the GNU General Public License v3.0 (see [LICENSE](LICENSE)). If you use Kummitus in published work, please [cite it](#how-to-cite).

---

## Installation

1. Install [xAct](http://www.xact.es/) (xTensor, xPert, xCoba, xTras) so that ``Needs["xAct`xTensor`"]`` works in your Mathematica.
2. Clone this repository and keep its files together: `Kummitus.wl` loads `ProjectorsLibrary.m` and `SystemSolver.wl` from its own folder.
3. Load the package by its path:

```mathematica
Get["/path/to/Kummitus/Kummitus.wl"]
```

## Quick start

A massive Fierz–Pauli graviton: one pole at `mg^2` with five positive residues.

```mathematica
Get["/path/to/Kummitus/Kummitus.wl"];
ConstantSymbols[mg];

model = 1/2 dd[-c][H[-a, -b]] dd[c][H[a, b]] - dd[-a][H[a, b]] dd[c][H[-c, -b]] +
   dd[-a][H[a, b]] dd[-b][H[c, -c]] - 1/2 dd[-a][H[b, -b]] dd[a][H[c, -c]] -
   mg^2/2 (H[-a, -b] H[a, b] - H[a, -a] H[b, -b]);

SourceDeclare[model];
CompSatProp[model];
masses                    (* {mg^2} *)
ComputeMassiveStates[]    (* {{mg^2, {2, 2, 2, 1, 1}}} *)
```

Single letters `a`…`t` (except `g`, `q`) are the tensor indices of the package: use multi-letter names for your own parameters and variables.

---

## How to cite

If Kummitus contributes to published work, please cite the Kummitus paper and the paper introducing the method:

- C. Marzo, *Kummitus: a light-weight toolbox for counting DOF in perturbative QFT*, [arXiv:2603.22451](https://arxiv.org/abs/2603.22451) (2026).
- C. Marzo, *Ghost and tachyon free propagation up to spin 3 in Lorentz invariant field theories*, [Phys. Rev. D 105, 065017 (2022)](https://doi.org/10.1103/PhysRevD.105.065017), [arXiv:2108.11982](https://arxiv.org/abs/2108.11982).

```bibtex
@article{Marzo:2026yjg,
    author = "Marzo, Carlo",
    title = "{Kummitus: a light-weight toolbox for counting DOF in perturbative QFT}",
    eprint = "2603.22451",
    archivePrefix = "arXiv",
    primaryClass = "hep-th",
    month = "3",
    year = "2026"
}

@article{Marzo:2021esg,
    author = "Marzo, Carlo",
    title = "{Ghost and tachyon free propagation up to spin 3 in Lorentz invariant field theories}",
    eprint = "2108.11982",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "MI-TH-1941",
    doi = "10.1103/PhysRevD.105.065017",
    journal = "Phys. Rev. D",
    volume = "105",
    number = "6",
    pages = "065017",
    year = "2022"
}
```

The same references are in [CITATION.cff](CITATION.cff), which GitHub shows as "Cite this repository".

---

## Overview

Given a free-field action written in terms of tensors of rank 0–3, Kummitus automates the full chain:

1. Vary the action twice to extract the kinetic operator matrices per spin-parity sector.
2. Determine gauge invariances and fix them to obtain irreducible blocks.
3. Find the propagator poles (candidate masses).
4. Construct the saturated propagator `J* Δ J`.
5. Compute residues at each pole and read off the physical spectrum.

---

## Files

| File | Role |
|---|---|
| `Kummitus.wl` | Main package — loads xAct, defines manifold/fields/sources, all core routines |
| `SystemSolver.wl` | Algebraic solver utilities (called internally by `Kummitus.wl`) |
| `ProjectorsLibrary.m` | Pre-computed spin-parity projectors `P[J,P,{i,j}]` for all field combinations |
| `FieldConfiguration.m` | (Optional) legacy field configuration helpers, not loaded by default |
| `Models.nb` | Example notebook with all the models of the Kummitus paper (arXiv:2603.22451), from Maxwell to massive spin-3; open it from this folder, it loads `Kummitus.wl` from its own directory |
| `MultiGravitonModels.wl` | Example models with several copies of a field: bigravity, graviton chains and rings, vector multiplets, mixed scalars; each with its expected and computed spectrum |
| `KalbRamondModels.wl` | Example models for the antisymmetric rank-2 (Kalb–Ramond) field |
| `LICENSE` | GNU General Public License v3.0 |
| `CITATION.cff` | Citation metadata (papers to cite) |

---

## Built-in Geometry

The package defines a fixed 4D Minkowski background on load:

| Symbol | Meaning |
|---|---|
| `M` | 4D manifold with Lorentz-signature indices `{a,b,c,...,t}` (excluding `g`,`q`) |
| `g[-a,-b]` | Flat metric with PD as fiducial derivative; `MetricRule1/2` give Lorentz components `diag(1,-1,-1,-1)` |
| `q[-a]` | Momentum 4-vector |
| `qq` | Momentum squared: `q[a] q[-a] = qq` |
| `qm` | Momentum magnitude: `qm^2 = qq` |
| `ω`, `κ`, `ε` | Constant symbols for energy, spatial momentum, perturbation parameter |
| `dimM` | Spacetime dimension (set to 4) |

---

## Fields and Sources

The package pre-defines the following tensor fields and their first-order perturbations:

| Field | Tensor | Perturbation | Source | Conj. source |
|---|---|---|---|---|
| Rank-3 `A` | `A[-a,-b,-c]` | `εAp[LI[1],-a,-b,-c]` | `AJ[a,b,c]` | `AJc[a,b,c]` |
| Rank-2 symmetric `H` | `H[-a,-b]` | `εh[LI[1],-a,-b]` | `HJ[a,b]` | `HJc[a,b]` |
| Vector `V` | `V[-a]` | `εVp[LI[1],-a]` | `VJ[a]` | `VJc[a]` |
| Scalar `S` | `S[]` | `εSp[LI[1]]` | `SJ` | `SJc` |

The symmetry type of the rank-3 field `A` is controlled by the string variable `rank3SymmetryType`, which must be set **before** loading the package:

```mathematica
rank3SymmetryType = "TotallySymmetric";  (* default *)
(* Options: "Generic", "Sym23", "Antisym23", "TotallySymmetric", "TotallyAntisymmetric" *)
```

### Declaring more fields (copies)

Further fields, including several copies of the same type, are declared after loading with

```mathematica
DeclareField[V2, "vector"]                         (* defines V2, εV2p, V2J, V2Jc *)
DeclareField[B, "rank3", "Symmetry" -> "Sym23"]    (* each rank-3 field may have its own symmetry *)
```

Supported types: `"rank3"`, `"rank2sym"`, `"rank2anti"` (antisymmetric, Kalb–Ramond), `"vector"`, `"scalar"`. Default names are `εXp` for the perturbation and `XJ`, `XJc` for the sources; the options `"Perturbation"`, `"Source"`, `"ConjugateSource"` override them. The built-in `A`, `H`, `V`, `S` are declared the same way at load time. `FieldInstances[]` lists the registry.

The antisymmetric rank-2 field has no built-in instance; declare one with `DeclareField[B, "rank2anti"]`. It lives in the 1⁺ (projector label 4) and 1⁻ (label 9) sectors and mixes with all other fields through the projectors in `ProjectorsLibrary.m`. The library sets the phase `XX3 = π`, which the antisymmetric mixing projectors require: with `XX3 = 0` the vector–antisymmetric products in 1⁻ do not close. `Kummitus.wl` therefore leaves `XX3` untouched (`XX1 = XX2 = XX4 = 0`).

Copies share the projectors: in each J^P sector the basis is the list of pairs `{field, label}`, ordered by rank (highest first) and then by declaration order. Gauge fixing removes, in order of preference, higher rank, then later-declared copies, then higher labels.

---

## Public Routines

### Setup

**`SourceDeclare[action]`**
> Detects which fields appear in `action` and initialises the source component rules. Call this once before `CompSatProp`. Expensive step, computed only once.

```mathematica
SourceDeclare[myAction]
```

### Core computation

**`CompSatProp[action]`**
> Main workhorse. Given a quadratic action, it:
> - Detects active fields and builds representation tables per spin-parity sector.
> - Varies the action twice to obtain the kinetic matrices `a(J^P)`.
> - Computes source gauge constraints (null spaces of each sector matrix).
> - Fixes gauges (removes redundant components, preferring scalar > vector > rank-2 > rank-3).
> - Extracts candidate mass poles from `det[a(J^P)] = 0`.
> - Constructs and stores the full saturated propagator.
>
> After this call, `masses` holds the list of poles and `prop` holds the full propagator.

```mathematica
CompSatProp[myAction]
```

### Spectrum extraction

**`ComputeMasslessStates[]`**
> Analyses the propagator near `qq → 0` (lightlike frame). Applies source constraints, solves for dependent source components, expands the saturated propagator as a Laurent series in `qq`, and diagonalises the residue matrix.
> Returns `masslessData = {0, {residues...}}`.
> Stores the propagator (as a Laurent series) in `LatestPropMassless`.

**`ComputeMassiveStates[]`**
> Same analysis at each non-zero mass pole `qq = m²` (rest frame).
> Returns `massiveData = {{m²₁, {residues...}}, {m²₂, ...}, ...}`.
> Stores propagators in `LatestPropMassive`.

### Helicity decomposition

**`TensorToHelicity[]`**
> Expresses the tensor source components in terms of helicity eigenstates `sv[h,i]`, `sh[h,i]`, `sa[h,i]` (and conjugates `svc`, `shc`, `sac`) for the vector, rank-2, and rank-3 fields respectively.
> The result, combined with `masslessSolution`, gives the propagator directly in helicity language.
> Only the fields present in the last `CompSatProp` are included. The built-in `V`, `H`, `A`, `S` keep the untagged names `sv[h,i]`, `sh[h,i]`, `sa[h,i]`, `ss[h,i]`; any other field carries its name as a third argument, e.g. `sv[h,i,V2]`, `sh[h,i,H2]`. Antisymmetric rank-2 fields use `sb[h,i,B]` (helicity ±1 from B_{i3}, B_{0i}; helicity 0 from B_12, B_03).

### Accessors

| Symbol | Content |
|---|---|
| `masses` | List of all poles `qq = m²` (including `0` if massless) |
| `masslessData` | `{0, {non-zero residues}}` after `ComputeMasslessStates[]` |
| `massiveData` | `{{m², {residues}}, ...}` after `ComputeMassiveStates[]` |
| `LatestPropMassless` | Saturated propagator Laurent-expanded near `qq = 0` |
| `LatestPropMassive` | List of propagators expanded near each massive pole |
| `matricesOutput` | List of raw `a(J^P)` matrices (for external use, e.g. Survey) |

### Utilities

| Function | Purpose |
|---|---|
| `dd[∂][expr]` | Alias for xAct's `PD[∂][expr]` (flat partial derivative) |
| `ConstantSymbols[sym]` | Wrapper for xAct's `DefConstantSymbol` |
| `PrintK[...]` | Styled print with grayed-out static text and highlighted expressions |

---

## Global Control Variables

Set these **before** loading the package (or before calling `CompSatProp`):

| Variable | Default | Effect |
|---|---|---|
| `rank3SymmetryType` | `"TotallySymmetric"` | Index symmetry of the rank-3 field `A` |
| `DivideEtImpera` | `True` | If `True`, sums propagator term-by-term before Laurent expansion (faster for large expressions) |

---

## How-to: Typical Workflow

The model files `MultiGravitonModels.wl` and `KalbRamondModels.wl` follow this pattern; load one after the package and pass any of its Lagrangians to the steps below.

### Step 0 — Configure and load

```mathematica
(* Optional: set before Get[] *)
rank3SymmetryType = "Antisym23";   (* or "TotallySymmetric", etc. *)
DivideEtImpera = True;

Get["/path/to/Kummitus/Kummitus.wl"]
```

### Step 1 — Write the action

Build a quadratic action using the pre-defined tensors `A`, `H`, `V`, `S` (or fields added with `DeclareField`) and the flat derivative `dd`. Example: a Proca field and a massive scalar.

```mathematica
ConstantSymbols /@ {mV, mS};
model = -1/2 (dd[-a][V[-b]] dd[a][V[b]] - dd[-a][V[-b]] dd[b][V[a]]) + mV^2/2 V[-a] V[a] +
   1/2 dd[-a][S[]] dd[a][S[]] - mS^2/2 S[]^2;
```

Two copies of the vector with kinetic mixing (a "dark photon"):

```mathematica
ConstantSymbols /@ {mV, ee};
DeclareField[V2, "vector"];
maxwellMix[X_, Y_] := -1/2 (dd[-a][X[-b]] dd[a][Y[b]] - dd[-a][X[-b]] dd[b][Y[a]]);
model = maxwellMix[V, V] + maxwellMix[V2, V2] + 2 ee maxwellMix[V, V2] + mV^2/2 V2[-a] V2[a];
```

### Step 2 — Declare sources

```mathematica
SourceDeclare[model]
```

This detects which fields appear and prepares component rules. It is slow but only needs to run once per field content.

### Step 3 — Compute the saturated propagator

```mathematica
CompSatProp[model]
```

The package prints the spin-parity sector matrices, source constraints, gauge-fixed matrices, and the candidate mass spectrum. Inspect `masses` to see the poles:

```mathematica
masses   (* e.g. {mV^2, mS^2} *)
```

### Step 4 — Extract the massless spectrum

```mathematica
ComputeMasslessStates[]
masslessData   (* {0, {residue₁, residue₂, ...}} *)
```

The non-zero residues correspond to physical massless degrees of freedom. Access the full propagator:

```mathematica
LatestPropMassless
```

### Step 5 — Extract the massive spectrum

```mathematica
ComputeMassiveStates[]
massiveData   (* {{m², {residues}}, ...} *)
```

### Step 6 — Helicity basis (massless case)

```mathematica
helRules = TensorToHelicity[]
LatestPropMassless /. helRules
```

This rewrites the propagator in terms of helicity eigenstates `sv[h,i]`, `sh[h,i]`, `sa[h,i]`.

### Compact pattern (notebook style)

For a scan over multiple actions or parameters, the notebook uses a pattern like:

```mathematica
SourceDeclare[model];
CompSatProp[model];

masslessDataX = If[MemberQ[masses, 0], ComputeMasslessStates[], {}];
massiveDataX  = If[DeleteCases[masses, 0] =!= {}, ComputeMassiveStates[], {}];

Print["Massless poles and residues: ", masslessDataX]
Print["Massive poles and residues: ", massiveDataX]
```

---

## Known limitations

- Spacetime dimension is fixed to D = 4, with signature (+,−,−,−).
- The projector operators are parity preserving: Lagrangians with parity-violating terms (ε tensor at quadratic order) are not supported.
- `ComputeMasslessStates[]` returns nothing when the massless pole is of higher order (e.g. four-derivative theories); inspect `LatestPropMassless` in that case.

---

## Dependencies

- Mathematica 12+
- xAct suite: `xTensor`, `xPert`, `xCoba`, `xTras`
  Available at [http://www.xact.es/](http://www.xact.es/)

---

## License

Kummitus is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License, version 3, as published by the Free Software Foundation. It is distributed without any warranty; see [LICENSE](LICENSE) for the full text.
