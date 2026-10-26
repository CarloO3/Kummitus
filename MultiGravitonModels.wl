(* ::Package:: *)

(* ============================================================
   MultiGravitonModels.wl -- linear multi-graviton test models
   for Kummitus.

   Usage (after loading the package):
     Get["/path/to/Kummitus/Kummitus.wl"];
     Get["/path/to/Kummitus/MultiGravitonModels.wl"];
     SourceDeclare[bigravity]; CompSatProp[bigravity];
     masses
     ComputeMasslessStates[]
     ComputeMassiveStates[]

   Conventions: mostly-minus metric diag(1,-1,-1,-1).
   "Got" lines are the results Kummitus returned on 2026-10-14.
   ============================================================ *)

If[!xAct`xTensor`ConstantSymbolQ[mg], ConstantSymbols[mg]];

(* Extra copies of the rank-2 symmetric field (H is built in). Skipped if already declared. *)
If[!KeyExistsQ[FieldInstances[], H2], DeclareField[H2, "rank2sym"]];
If[!KeyExistsQ[FieldInstances[], H3], DeclareField[H3, "rank2sym"]];

(* Linearized Einstein-Hilbert (massless Fierz-Pauli) Lagrangian *)
EH[X_] := 1/2 dd[-c][X[-a, -b]] dd[c][X[a, b]] - dd[-a][X[a, b]] dd[c][X[-c, -b]] +
    dd[-a][X[a, b]] dd[-b][X[c, -c]] - 1/2 dd[-a][X[b, -b]] dd[a][X[c, -c]];

(* Fierz-Pauli mass structure for the difference X - Y *)
FPdiff[X_, Y_] := (X[-a, -b] - Y[-a, -b]) (X[a, b] - Y[a, b]) -
    (X[a, -a] - Y[a, -a]) (X[b, -b] - Y[b, -b]);


(* ------------------------------------------------------------
   1. Bigravity: two gravitons, mass term on the difference.
   h+ = (H + H2)/Sqrt[2]: massless graviton, 2 dof.
   h- = (H - H2)/Sqrt[2]: massive Fierz-Pauli graviton, M^2 = 2 mg^2, 5 dof.
   Expected: masses = {0, 2 mg^2}; 2 residues at qq = 0, 5 at qq = 2 mg^2.
   Got (38 s): {0, {2, 1}}; {{2 mg^2, {2, 2, 2, 1, 1}}}.
   ------------------------------------------------------------ *)

bigravity = EH[H] + EH[H2] - mg^2/2 FPdiff[H, H2];


(* ------------------------------------------------------------
   2. Three-graviton open chain (deconstructed dimension):
   H - H2 - H3 with equal nearest-neighbour Fierz-Pauli masses.
   Mass matrix = mg^2 x Laplacian of the 3-site path, eigenvalues
   M^2 = 4 mg^2 Sin[k Pi/6]^2, k = 0, 1, 2  ->  {0, mg^2, 3 mg^2}.
   Eigenvectors: (1,1,1)/Sqrt[3], (1,0,-1)/Sqrt[2], (1,-2,1)/Sqrt[6].
   Expected: masses = {0, mg^2, 3 mg^2}; 2 residues at qq = 0,
   5 at each massive pole. Only the diagonal diffeomorphism survives.
   Got (129 s): {0, {2, 1}}; {{mg^2, {2, 2, 2, 1, 1}}, {3 mg^2, {2, 2, 2, 1, 1}}}.
   ------------------------------------------------------------ *)

gravitonChain3 = EH[H] + EH[H2] + EH[H3] -
    mg^2/2 (FPdiff[H, H2] + FPdiff[H2, H3]);


(* ============================================================
   OTHER MULTI-FIELD MODELS
   Copies of vectors, scalars and gravitons, with mass and kinetic
   mixing. Each has a spectrum known in closed form.
   ============================================================ *)

If[!xAct`xTensor`ConstantSymbolQ[#], ConstantSymbols[#]] & /@ {mv, ee, al};

If[!KeyExistsQ[FieldInstances[], V2], DeclareField[V2, "vector"]];
If[!KeyExistsQ[FieldInstances[], V3], DeclareField[V3, "vector"]];
If[!KeyExistsQ[FieldInstances[], S2], DeclareField[S2, "scalar"]];

(* -1/4 F_ab F^ab, and the mixed term -1/4 F(X)_ab F(Y)^ab (so 2 ee maxwellMix = -ee/2 F.F2) *)
maxwell[X_] := -1/2 (dd[-a][X[-b]] dd[a][X[b]] - dd[-a][X[-b]] dd[b][X[a]]);
maxwellMix[X_, Y_] := -1/2 (dd[-a][X[-b]] dd[a][Y[b]] - dd[-a][X[-b]] dd[b][Y[a]]);


(* ------------------------------------------------------------
   3. Proca triplet with a non-diagonal mass matrix
   M^2 = mv^2 {{2,-1,0},{-1,2,-1},{0,-1,2}} on (V, V2, V3).
   Eigenvalues: mv^2 (2 - Sqrt[2]), 2 mv^2, mv^2 (2 + Sqrt[2]).
   Expected: three massive poles, 3 residues each; no massless pole.
   Got (3 s): {1, 1, 1} at each of the three poles.
   ------------------------------------------------------------ *)

procaTriplet = maxwell[V] + maxwell[V2] + maxwell[V3] +
    mv^2/2 (2 V[-a] V[a] + 2 V2[-a] V2[a] + 2 V3[-a] V3[a] - 2 V[-a] V2[a] - 2 V2[-a] V3[a]);


(* ------------------------------------------------------------
   4. Kinetically mixed photons ("dark photon"): -ee/2 F.F2, and a
   mass for V2 only. det(q^2 K - M) = q^2 ((1 - ee^2) q^2 - mv^2), so
   Expected: masses {0, mv^2/(1 - ee^2)} (needs ee^2 < 1);
   2 residues at 0, 3 at the massive pole.
   Got (1 s): {0, {1, 1}}; 3 residues (1 + ee^2)/(1 - ee^2) at mv^2/(1 - ee^2).
   ------------------------------------------------------------ *)

darkPhoton = maxwell[V] + maxwell[V2] + 2 ee maxwellMix[V, V2] + mv^2/2 V2[-a] V2[a];


(* ------------------------------------------------------------
   5. Stueckelberg Proca with copies: V2 and S2 with the gauge symmetry
   dV2 = d xi, dS2 = mv xi.  mv^2/2 (V2 - dS2/mv)^2.
   Expected: masses {mv^2} (plus a pure-gauge 0 candidate, if any),
   3 residues at mv^2.
   Got: masses {mv^2, 0}; {0, {}} (pure gauge); {{mv^2, {1, 1, 1}}}.
   ------------------------------------------------------------ *)

stueckelbergCopies = maxwell[V2] +
    mv^2/2 (V2[-a] - dd[-a][S2[]]/mv) (V2[a] - dd[a][S2[]]/mv);


(* ------------------------------------------------------------
   6. Bigravity with different Planck masses: kinetic weights 1 and al^2.
   det = q^2 (al^2 q^2 - mg^2 (1 + al^2)), so
   Expected: masses {0, mg^2 (1 + al^2)/al^2};
   2 residues at 0, 5 at the massive pole.
   Got (50 s): {0, {4, 2}/(1 + al^2)};
     {(1 + 1/al^2) mg^2, (1 + al^4)/(al^2 (1 + al^2)) {1, 1, 2, 2, 2}}.
     At al = 1 this reproduces bigravity.
   ------------------------------------------------------------ *)

bigravityAsym = EH[H] + al^2 EH[H2] - mg^2/2 FPdiff[H, H2];


(* ------------------------------------------------------------
   7. Three-graviton ring (closed chain): adds the H3 - H link.
   Mass matrix = mg^2 x Laplacian of the 3-cycle: eigenvalues 0, 3, 3.
   Degenerate massive pole: two massive gravitons at 3 mg^2.
   Expected: masses {0, 3 mg^2}; 2 residues at 0, 10 at 3 mg^2.
   Got (113 s): {0, {2, 1}}; {{3 mg^2, {2, 2, 2, 2, 2, 2, 1, 1, 1, 1}}}.
   ------------------------------------------------------------ *)

gravitonRing3 = EH[H] + EH[H2] + EH[H3] -
    mg^2/2 (FPdiff[H, H2] + FPdiff[H2, H3] + FPdiff[H3, H]);


(* ------------------------------------------------------------
   8. Two mass-mixed scalars.
   Mass matrix mv^2 {{1,-1},{-1,1}}: eigenvalues 0, 2 mv^2.
   Expected: masses {0, 2 mv^2}; 1 residue at each pole.
   Got (2026-10-15, after the single-term fix): {0, {1}}; {{2 mv^2, {1}}}.
   ------------------------------------------------------------ *)

scalarPair = 1/2 dd[-a][S[]] dd[a][S[]] + 1/2 dd[-a][S2[]] dd[a][S2[]] -
    mv^2/2 (S[] - S2[])^2;
