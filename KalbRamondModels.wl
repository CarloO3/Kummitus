(* ::Package:: *)

(* ============================================================
   KalbRamondModels.wl -- antisymmetric rank-2 (Kalb-Ramond)
   test models for Kummitus.

   Usage (after loading the package):
     Get["/path/to/Kummitus/Kummitus.wl"];
     Get["/path/to/Kummitus/KalbRamondModels.wl"];
     SourceDeclare[KRmassive]; CompSatProp[KRmassive];
     masses
     ComputeMasslessStates[]
     ComputeMassiveStates[]

   Conventions: mostly-minus metric diag(1,-1,-1,-1).
   B lives in the 1+ (label 4) and 1- (label 9) sectors; the mixing
   projectors require XX3 = Pi (set in ProjectorsLibrary.m).
   Residues quoted below are the values Kummitus returned on 2026-10-14.
   ============================================================ *)

If[!xAct`xTensor`ConstantSymbolQ[#], ConstantSymbols[#]] & /@ {mB, cc};

(* Antisymmetric fields. Skipped if already declared. *)
If[!KeyExistsQ[FieldInstances[], B], DeclareField[B, "rank2anti"]];
If[!KeyExistsQ[FieldInstances[], B2], DeclareField[B2, "rank2anti"]];

(* (1/12) H_abc H^abc with H = dB *)
KR[X_] := 1/4 dd[-a][X[-b, -c]] dd[a][X[b, c]] + 1/2 dd[-a][X[-b, -c]] dd[b][X[c, a]];

(* Field strength of a vector, indices down and up *)
FSdn[X_] := dd[-a][X[-b]] - dd[-b][X[-a]];
FSup[X_] := dd[a][X[b]] - dd[b][X[a]];

(* Linearized Einstein-Hilbert Lagrangian *)
EH[X_] := 1/2 dd[-c][X[-a, -b]] dd[c][X[a, b]] - dd[-a][X[a, b]] dd[c][X[-c, -b]] +
    dd[-a][X[a, b]] dd[-b][X[c, -c]] - 1/2 dd[-a][X[b, -b]] dd[a][X[c, -c]];


(* ------------------------------------------------------------
   1. Massless Kalb-Ramond. Gauge invariance dB = d(xi) removes
   the 1- part; 1 dof in 4D (dual to a massless scalar).
   Expected: masses = {0}; 1 residue. Got: {0, {4}}.
   ------------------------------------------------------------ *)

KRmassless = KR[B];


(* ------------------------------------------------------------
   2. Massive Kalb-Ramond. Pole in 1+; 1- is non-dynamical.
   3 dof (massive 2-form = axial vector, dual to Proca).
   Expected: masses = {mB^2}; 3 residues. Got: {{mB^2, {4, 4, 4}}}.
   ------------------------------------------------------------ *)

KRmassive = KR[B] - mB^2/4 B[-a, -b] B[a, b];


(* ------------------------------------------------------------
   3. Stueckelberg form: mass term on B - dV/mB. Gauge symmetry
   dB = dxi, dV = mB xi restores the vector gauge invariance.
   Tests the B-V mixing projectors (1- labels 8, 9).
   Expected: same as model 2; the 0 pole is pure gauge.
   Got: masses {mB^2, 0}; massless {0, {}}; massive {{mB^2, {4, 4, 4}}}.
   ------------------------------------------------------------ *)

KRstueckelberg = KR[B] -
    mB^2/4 (B[-a, -b] - FSdn[V]/mB) (B[a, b] - FSup[V]/mB);


(* ------------------------------------------------------------
   4. Two copies with a mass term on the difference.
   B+ = (B + B2)/Sqrt[2]: massless KR, 1 dof.
   B- = (B - B2)/Sqrt[2]: massive KR with M^2 = 2 mB^2, 3 dof.
   Expected: masses = {0, 2 mB^2}. Got: {0, {4}}, {{2 mB^2, {4, 4, 4}}}.
   ------------------------------------------------------------ *)

KRtwoCopies = KR[B] + KR[B2] -
    mB^2/4 (B[-a, -b] - B2[-a, -b]) (B[a, b] - B2[a, b]);


(* ------------------------------------------------------------
   5. Graviton + massive KR with derivative mixing
   cc (d^a H_ab)(d_c B^cb), which lives in the 1- sector
   (H label 7, B label 9).
   The coupling is NOT invariant under linearized diffeomorphisms
   dH_ab = d_a xi_b + d_b xi_a: the 1- part of H is pure gauge in EH,
   so any H-B mixing there breaks the gauge symmetry. Spectrum not
   known a priori; use this model to look at what the breaking does.
   Got (cc generic): masses {0, mB^2};
     massive {{mB^2, {4, 4, 4}}};
     massless: Null (ComputeMasslessStates found higher-order poles or
     did not simplify; inspect LatestPropMassless).
     1- matrix: {{0, -cc qq/2}, {-cc qq/2, -mB^2/2}} in basis {H7, B9}.
   ------------------------------------------------------------ *)

KRplusH = EH[H] + KR[B] - mB^2/4 B[-a, -b] B[a, b] +
    cc dd[-a][H[a, b]] dd[c][B[-c, -b]];
