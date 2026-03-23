# Kummitus — DOF from Poles

![Kummitus](Kummitus.png)

**Kummitus** is a Mathematica package for computing the physical degrees of freedom of higher-spin field theories directly from the poles of the saturated propagator. It is built on top of the [xAct](http://www.xact.es/) tensor algebra system (xTensor, xPert, xCoba, xTras).

Copyright © 2026, Carlo Marzo — General Public License.

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
| `FieldConfiguration.m` | (Optional) legacy field configuration helpers |
| `TheNotebook.nb` | Example notebook — illustrative computations and workflow |

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
| `DivideEtImpera` | `"True"` | If `True`, sums propagator term-by-term before Laurent expansion (faster for large expressions) |

---

## How-to: Typical Workflow

The examples in `TheNotebook.nb` follow this pattern.

### Step 0 — Configure and load

```mathematica
(* Optional: set before Get[] *)
rank3SymmetryType = "Antisym23";   (* or "TotallySymmetric", etc. *)
DivideEtImpera = True;

Get["/path/to/Kummitus.wl"]
```

### Step 1 — Write the action

Build a quadratic action using the pre-defined tensors `A`, `H`, `V`, `S` and the flat derivative `dd`. Example (Fierz–Pauli-like for a rank-2 field `H`):

```mathematica
model = dd[-a][H[-b,-c]] dd[a][H[b,c]] - (1/2) dd[-a][H[b,b]] dd[a][H[c,c]]
```

For a more complex mixed action involving a rank-3 field `A` and a scalar `S`:

```mathematica
max = dd[-a][A[-b,-c,-d]] dd[a][A[b,c,d]] + m2 A[-a,-b,-c] A[a,b,c] + ...
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
masses   (* e.g. {0, m2} *)
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

## Dependencies

- Mathematica 12+
- xAct suite: `xTensor`, `xPert`, `xCoba`, `xTras`
  Available at [http://www.xact.es/](http://www.xact.es/)
