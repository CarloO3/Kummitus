(* ::Package:: *)

(* Projector Library for Spectral Analysis *)
(* Symmetry: Symmetric in last 2 indices (Sym23) *)
(* Extracted from ProjectorsLibrary.nb *)

(* ============================================== *)
(* USAGE:                                         *)
(* Load AFTER setting up xAct and defining M.     *)
(*                                                *)
(* << "/path/to/ProjectorsLibrary.m"              *)
(*                                                *)
(* Required pre-defined symbols:                  *)
(*   - q[-a]: momentum vector                     *)
(*   - qm, qq: momentum magnitude symbols         *)
(*   - delta[a,b]: Kronecker delta                *)
(*   - dimM: manifold dimension                   *)
(* ============================================== *)

(* === Momentum magnitude relations === *)
qm /: qm^2 := qq;
qm /: qm^4 := qq^2;
qm /: qm^6 := qq^3;
qm /: qm^8 := qq^4;

(* === Momentum contraction rules === *)
q /: q[a_] q[-a_] := qq;
q /: q[-a_] q[a_] := qq;
(*
(* === Control variables for mixed operators === *)
DefConstantSymbol[#]&/@{XX1,XX2,XX3,XX4};
*)
(* === Transverse and Longitudinal Projectors === *)
PT[x_, y_] := delta[x, y] - q[x] q[y]/qq;
PL[x_, y_] := q[x] q[y]/qq;

(* === Default: undefined projectors return 0 === *)
P[x_, p_, {y_, z_}][{m___}] := P[x, p, {y, z}][m];
P[x_, p_, {y_, z_}][m___] := 0;


(* ============================================== *)
(* PROJECTOR DEFINITIONS BY SECTOR                *)
(* ============================================== *)

(*Projectors definitions, in InputForm (1912.01023 paper`s typos taken care of) (sym->2,3)*)

(*This should allow to use list oriented tools on the arguments*)
P[x_,p_,{y_,z_}][{m___}]:=P[x,p,{y,z}][m];
P[x_,p_,{y_,z_}][m___]:=0;
(*Sector0-*)
P[0,m,{1,1}][c_,a_,b_,d_,e_,f_]:=-(1/6) PT[a,f] PT[b,e] PT[c,d]+1/6 PT[a,e] PT[b,f] PT[c,d]+1/6 PT[a,f] PT[b,d] PT[c,e]-1/6 PT[a,d] PT[b,f] PT[c,e]-1/6 PT[a,e] PT[b,d] PT[c,f]+1/6 PT[a,d] PT[b,e] PT[c,f]
(*Sector0+*)
P[0,p,{1,1}][c_,a_,b_,d_,e_,f_]:=(PL[c,f] PT[a,b] PT[d,e])/(3 (-1+dimM))+(PL[b,f] PT[a,c] PT[d,e])/(3 (-1+dimM))+(PL[a,f] PT[b,c] PT[d,e])/(3 (-1+dimM))+(PL[c,e] PT[a,b] PT[d,f])/(3 (-1+dimM))+(PL[b,e] PT[a,c] PT[d,f])/(3 (-1+dimM))+(PL[a,e] PT[b,c] PT[d,f])/(3 (-1+dimM))+(PL[c,d] PT[a,b] PT[e,f])/(3 (-1+dimM))+(PL[b,d] PT[a,c] PT[e,f])/(3 (-1+dimM))+(PL[a,d] PT[b,c] PT[e,f])/(3 (-1+dimM))

P[0,p,{1,2}][c_,a_,b_,d_,e_,f_]:=-((PL[c,f] PT[a,b] PT[d,e])/(3 Sqrt[2] (-1+dimM)))-(PL[b,f] PT[a,c] PT[d,e])/(3 Sqrt[2] (-1+dimM))-(PL[a,f] PT[b,c] PT[d,e])/(3 Sqrt[2] (-1+dimM))-(PL[c,e] PT[a,b] PT[d,f])/(3 Sqrt[2] (-1+dimM))-(PL[b,e] PT[a,c] PT[d,f])/(3 Sqrt[2] (-1+dimM))-(PL[a,e] PT[b,c] PT[d,f])/(3 Sqrt[2] (-1+dimM))+(Sqrt[2] PL[c,d] PT[a,b] PT[e,f])/(3 (-1+dimM))+(Sqrt[2] PL[b,d] PT[a,c] PT[e,f])/(3 (-1+dimM))+(Sqrt[2] PL[a,d] PT[b,c] PT[e,f])/(3 (-1+dimM))

P[0,p,{1,3}][c_,a_,b_,d_,e_,f_]:=-((PL[c,f] PT[a,b] PT[d,e])/(Sqrt[6] (-1+dimM)))-(PL[b,f] PT[a,c] PT[d,e])/(Sqrt[6] (-1+dimM))-(PL[a,f] PT[b,c] PT[d,e])/(Sqrt[6] (-1+dimM))+(PL[c,e] PT[a,b] PT[d,f])/(Sqrt[6] (-1+dimM))+(PL[b,e] PT[a,c] PT[d,f])/(Sqrt[6] (-1+dimM))+(PL[a,e] PT[b,c] PT[d,f])/(Sqrt[6] (-1+dimM))

P[0,p,{1,4}][c_,a_,b_,d_,e_,f_]:=1/(3 Sqrt[3]) Sqrt[1/(-1+dimM)] (PL[c,f] PL[d,e] PT[a,b]+PL[c,e] PL[d,f] PT[a,b]+PL[c,d] PL[e,f] PT[a,b]+PL[b,f] PL[d,e] PT[a,c]+PL[b,e] PL[d,f] PT[a,c]+PL[b,d] PL[e,f] PT[a,c]+PL[a,f] PL[d,e] PT[b,c]+PL[a,e] PL[d,f] PT[b,c]+PL[a,d] PL[e,f] PT[b,c])


P[0,p,{1,5}][c_,a_,b_,e_,f_]:=(PT[b,c] PT[f,e] q[a])/(Sqrt[3] (-1+dimM) qm)+(PT[a,c] PT[f,e] q[b])/(Sqrt[3] (-1+dimM) qm)+(PT[a,b] PT[f,e] q[c])/(Sqrt[3] (-1+dimM) qm)

P[0,p,{1,6}][c_,a_,b_,e_,f_]:=(PL[f,e] PT[b,c] q[a])/(Sqrt[3] Sqrt[-1+dimM] qm)+(PL[f,e] PT[a,c] q[b])/(Sqrt[3] Sqrt[-1+dimM] qm)+(PL[f,e] PT[a,b] q[c])/(Sqrt[3] Sqrt[-1+dimM] qm)

P[0,p,{2,1}][c_,a_,b_,d_,e_,f_]:=(Sqrt[2] PL[c,f] PT[a,b] PT[d,e])/(3 (-1+dimM))-(PL[b,f] PT[a,c] PT[d,e])/(3 Sqrt[2] (-1+dimM))-(PL[a,f] PT[b,c] PT[d,e])/(3 Sqrt[2] (-1+dimM))+(Sqrt[2] PL[c,e] PT[a,b] PT[d,f])/(3 (-1+dimM))-(PL[b,e] PT[a,c] PT[d,f])/(3 Sqrt[2] (-1+dimM))-(PL[a,e] PT[b,c] PT[d,f])/(3 Sqrt[2] (-1+dimM))+(Sqrt[2] PL[c,d] PT[a,b] PT[e,f])/(3 (-1+dimM))-(PL[b,d] PT[a,c] PT[e,f])/(3 Sqrt[2] (-1+dimM))-(PL[a,d] PT[b,c] PT[e,f])/(3 Sqrt[2] (-1+dimM))

P[0,p,{2,2}][c_,a_,b_,d_,e_,f_]:=-((PL[c,f] PT[a,b] PT[d,e])/(3 (-1+dimM)))+(PL[b,f] PT[a,c] PT[d,e])/(6 (-1+dimM))+(PL[a,f] PT[b,c] PT[d,e])/(6 (-1+dimM))-(PL[c,e] PT[a,b] PT[d,f])/(3 (-1+dimM))+(PL[b,e] PT[a,c] PT[d,f])/(6 (-1+dimM))+(PL[a,e] PT[b,c] PT[d,f])/(6 (-1+dimM))+(2 PL[c,d] PT[a,b] PT[e,f])/(3 (-1+dimM))-(PL[b,d] PT[a,c] PT[e,f])/(3 (-1+dimM))-(PL[a,d] PT[b,c] PT[e,f])/(3 (-1+dimM))

P[0,p,{2,3}][c_,a_,b_,d_,e_,f_]:=-((PL[c,f] PT[a,b] PT[d,e])/(Sqrt[3] (-1+dimM)))+(PL[b,f] PT[a,c] PT[d,e])/(2 Sqrt[3] (-1+dimM))+(PL[a,f] PT[b,c] PT[d,e])/(2 Sqrt[3] (-1+dimM))+(PL[c,e] PT[a,b] PT[d,f])/(Sqrt[3] (-1+dimM))-(PL[b,e] PT[a,c] PT[d,f])/(2 Sqrt[3] (-1+dimM))-(PL[a,e] PT[b,c] PT[d,f])/(2 Sqrt[3] (-1+dimM))

P[0,p,{2,4}][c_,a_,b_,d_,e_,f_]:=1/3 Sqrt[2/3] Sqrt[1/(-1+dimM)] PL[c,f] PL[d,e] PT[a,b]+1/3 Sqrt[2/3] Sqrt[1/(-1+dimM)] PL[c,e] PL[d,f] PT[a,b]+1/3 Sqrt[2/3] Sqrt[1/(-1+dimM)] PL[c,d] PL[e,f] PT[a,b]-(Sqrt[1/(-1+dimM)] PL[b,f] PL[d,e] PT[a,c])/(3 Sqrt[6])-(Sqrt[1/(-1+dimM)] PL[b,e] PL[d,f] PT[a,c])/(3 Sqrt[6])-(Sqrt[1/(-1+dimM)] PL[b,d] PL[e,f] PT[a,c])/(3 Sqrt[6])-(Sqrt[1/(-1+dimM)] PL[a,f] PL[d,e] PT[b,c])/(3 Sqrt[6])-(Sqrt[1/(-1+dimM)] PL[a,e] PL[d,f] PT[b,c])/(3 Sqrt[6])-(Sqrt[1/(-1+dimM)] PL[a,d] PL[e,f] PT[b,c])/(3 Sqrt[6])

P[0,p,{2,5}][c_,a_,b_,e_,f_]:=-((PT[b,c] PT[e,f] q[a])/(Sqrt[6] (-1+dimM) qm))-(PT[a,c] PT[e,f] q[b])/(Sqrt[6] (-1+dimM) qm)+(Sqrt[2/3] PT[a,b] PT[e,f] q[c])/((-1+dimM) qm)

P[0,p,{2,6}][c_,a_,b_,e_,f_]:=-((PL[e,f] PT[b,c] q[a])/(Sqrt[6] Sqrt[-1+dimM] qm))-(PL[e,f] PT[a,c] q[b])/(Sqrt[6] Sqrt[-1+dimM] qm)+(Sqrt[2/3] PL[e,f] PT[a,b] q[c])/(Sqrt[-1+dimM] qm)

P[0,p,{3,1}][c_,a_,b_,d_,e_,f_]:=-((PL[b,f] PT[a,c] PT[d,e])/(Sqrt[6] (-1+dimM)))+(PL[a,f] PT[b,c] PT[d,e])/(Sqrt[6] (-1+dimM))-(PL[b,e] PT[a,c] PT[d,f])/(Sqrt[6] (-1+dimM))+(PL[a,e] PT[b,c] PT[d,f])/(Sqrt[6] (-1+dimM))-(PL[b,d] PT[a,c] PT[e,f])/(Sqrt[6] (-1+dimM))+(PL[a,d] PT[b,c] PT[e,f])/(Sqrt[6] (-1+dimM))

P[0,p,{3,2}][c_,a_,b_,d_,e_,f_]:=(PL[b,f] PT[a,c] PT[d,e])/(2 Sqrt[3] (-1+dimM))-(PL[a,f] PT[b,c] PT[d,e])/(2 Sqrt[3] (-1+dimM))+(PL[b,e] PT[a,c] PT[d,f])/(2 Sqrt[3] (-1+dimM))-(PL[a,e] PT[b,c] PT[d,f])/(2 Sqrt[3] (-1+dimM))-(PL[b,d] PT[a,c] PT[e,f])/(Sqrt[3] (-1+dimM))+(PL[a,d] PT[b,c] PT[e,f])/(Sqrt[3] (-1+dimM))

P[0,p,{3,3}][c_,a_,b_,d_,e_,f_]:=(PL[b,f] PT[a,c] PT[d,e])/(2 (-1+dimM))-(PL[a,f] PT[b,c] PT[d,e])/(2 (-1+dimM))-(PL[b,e] PT[a,c] PT[d,f])/(2 (-1+dimM))+(PL[a,e] PT[b,c] PT[d,f])/(2 (-1+dimM))

P[0,p,{3,4}][c_,a_,b_,d_,e_,f_]:=-((PL[b,f] PL[d,e] PT[a,c])/(3 Sqrt[2] Sqrt[-1+dimM]))-(PL[b,e] PL[d,f] PT[a,c])/(3 Sqrt[2] Sqrt[-1+dimM])-(PL[b,d] PL[e,f] PT[a,c])/(3 Sqrt[2] Sqrt[-1+dimM])+(PL[a,f] PL[d,e] PT[b,c])/(3 Sqrt[2] Sqrt[-1+dimM])+(PL[a,e] PL[d,f] PT[b,c])/(3 Sqrt[2] Sqrt[-1+dimM])+(PL[a,d] PL[e,f] PT[b,c])/(3 Sqrt[2] Sqrt[-1+dimM])

P[0,p,{3,5}][c_,a_,b_,e_,f_]:=(PT[b,c] PT[e,f] q[a])/(3 Sqrt[2] (-1+dimM) qm)+(Sqrt[2] PT[b,c] PT[e,f] q[a])/(3 (-1+dimM) qm)-(PT[a,c] PT[e,f] q[b])/(3 Sqrt[2] (-1+dimM) qm)-(Sqrt[2] PT[a,c] PT[e,f] q[b])/(3 (-1+dimM) qm)

P[0,p,{3,6}][c_,a_,b_,e_,f_]:=(PL[e,f] PT[b,c] q[a])/(3 Sqrt[2] Sqrt[-1+dimM] qm)+(Sqrt[2] PL[e,f] PT[b,c] q[a])/(3 Sqrt[-1+dimM] qm)-(PL[e,f] PT[a,c] q[b])/(3 Sqrt[2] Sqrt[-1+dimM] qm)-(Sqrt[2] PL[e,f] PT[a,c] q[b])/(3 Sqrt[-1+dimM] qm)

P[0,p,{4,1}][c_,a_,b_,d_,e_,f_]:=1/(3 Sqrt[3]) Sqrt[1/(-1+dimM)] (PL[a,f] PL[b,c] PT[d,e]+PL[a,b] PL[c,f] PT[d,e]+PL[a,e] PL[b,c] PT[d,f]+PL[a,b] PL[c,e] PT[d,f]+PL[a,d] PL[b,c] PT[e,f]+PL[a,b] PL[c,d] PT[e,f]+PL[a,c] (PL[b,f] PT[d,e]+PL[b,e] PT[d,f]+PL[b,d] PT[e,f]))

P[0,p,{4,2}][c_,a_,b_,d_,e_,f_]:=-((Sqrt[1/(-1+dimM)] PL[a,f] PL[b,c] PT[d,e])/(3 Sqrt[6]))-(Sqrt[1/(-1+dimM)] PL[a,c] PL[b,f] PT[d,e])/(3 Sqrt[6])-(Sqrt[1/(-1+dimM)] PL[a,b] PL[c,f] PT[d,e])/(3 Sqrt[6])-(Sqrt[1/(-1+dimM)] PL[a,e] PL[b,c] PT[d,f])/(3 Sqrt[6])-(Sqrt[1/(-1+dimM)] PL[a,c] PL[b,e] PT[d,f])/(3 Sqrt[6])-(Sqrt[1/(-1+dimM)] PL[a,b] PL[c,e] PT[d,f])/(3 Sqrt[6])+1/3 Sqrt[2/3] Sqrt[1/(-1+dimM)] PL[a,d] PL[b,c] PT[e,f]+1/3 Sqrt[2/3] Sqrt[1/(-1+dimM)] PL[a,c] PL[b,d] PT[e,f]+1/3 Sqrt[2/3] Sqrt[1/(-1+dimM)] PL[a,b] PL[c,d] PT[e,f]

P[0,p,{4,3}][c_,a_,b_,d_,e_,f_]:=-((PL[a,f] PL[b,c] PT[d,e])/(3 Sqrt[2] Sqrt[-1+dimM]))-(PL[a,c] PL[b,f] PT[d,e])/(3 Sqrt[2] Sqrt[-1+dimM])-(PL[a,b] PL[c,f] PT[d,e])/(3 Sqrt[2] Sqrt[-1+dimM])+(PL[a,e] PL[b,c] PT[d,f])/(3 Sqrt[2] Sqrt[-1+dimM])+(PL[a,c] PL[b,e] PT[d,f])/(3 Sqrt[2] Sqrt[-1+dimM])+(PL[a,b] PL[c,e] PT[d,f])/(3 Sqrt[2] Sqrt[-1+dimM])

P[0,p,{4,4}][c_,a_,b_,d_,e_,f_]:=1/6 PL[a,f] PL[b,e] PL[c,d]+1/6 PL[a,e] PL[b,f] PL[c,d]+1/6 PL[a,f] PL[b,d] PL[c,e]+1/6 PL[a,d] PL[b,f] PL[c,e]+1/6 PL[a,e] PL[b,d] PL[c,f]+1/6 PL[a,d] PL[b,e] PL[c,f]

P[0,p,{4,5}][c_,a_,b_,e_,f_]:=(PL[b,c] PT[e,f] q[a])/(3 Sqrt[-1+dimM] qm)+(PL[a,c] PT[e,f] q[b])/(3 Sqrt[-1+dimM] qm)+(PL[a,b] PT[e,f] q[c])/(3 Sqrt[-1+dimM] qm)

P[0,p,{4,6}][c_,a_,b_,e_,f_]:=(PL[b,c] PL[e,f] q[a])/(3 qm)+(PL[a,c] PL[e,f] q[b])/(3 qm)+(PL[a,b] PL[e,f] q[c])/(3 qm)

P[0,p,{5,1}][a_,b_,d_,e_,f_]:=(PT[a,b] PT[e,f] q[d])/(Sqrt[3] (-1+dimM) qm)+(PT[a,b] PT[d,f] q[e])/(Sqrt[3] (-1+dimM) qm)+(PT[a,b] PT[d,e] q[f])/(Sqrt[3] (-1+dimM) qm)

P[0,p,{5,2}][a_,b_,d_,e_,f_]:=(Sqrt[2/3] PT[a,b] PT[e,f] q[d])/((-1+dimM) qm)-(PT[a,b] PT[d,f] q[e])/(Sqrt[6] (-1+dimM) qm)-(PT[a,b] PT[d,e] q[f])/(Sqrt[6] (-1+dimM) qm)

P[0,p,{5,3}][a_,b_,d_,e_,f_]:=(PT[a,b] PT[d,f] q[e])/(3 Sqrt[2] (-1+dimM) qm)+(Sqrt[2] PT[a,b] PT[d,f] q[e])/(3 (-1+dimM) qm)-(PT[a,b] PT[d,e] q[f])/(3 Sqrt[2] (-1+dimM) qm)-(Sqrt[2] PT[a,b] PT[d,e] q[f])/(3 (-1+dimM) qm)

P[0,p,{5,4}][a_,b_,d_,e_,f_]:=(PL[e,f] PT[a,b] q[d])/(3 Sqrt[-1+dimM] qm)+(PL[d,f] PT[a,b] q[e])/(3 Sqrt[-1+dimM] qm)+(PL[d,e] PT[a,b] q[f])/(3 Sqrt[-1+dimM] qm)

P[0,p,{5,5}][a_,b_,e_,f_]:=1/(-1+dimM) PT[a,b] PT[e,f]

P[0,p,{5,6}][a_,b_,e_,f_]:=(PL[e,f] PT[a,b])/Sqrt[-1+dimM]

P[0,p,{6,1}][a_,b_,d_,e_,f_]:=(PL[a,b] PT[e,f] q[d])/(Sqrt[3] Sqrt[-1+dimM] qm)+(PL[a,b] PT[d,f] q[e])/(Sqrt[3] Sqrt[-1+dimM] qm)+(PL[a,b] PT[d,e] q[f])/(Sqrt[3] Sqrt[-1+dimM] qm)

P[0,p,{6,2}][a_,b_,d_,e_,f_]:=(Sqrt[2/3] PL[a,b] PT[e,f] q[d])/(Sqrt[-1+dimM] qm)-(PL[a,b] PT[d,f] q[e])/(Sqrt[6] Sqrt[-1+dimM] qm)-(PL[a,b] PT[d,e] q[f])/(Sqrt[6] Sqrt[-1+dimM] qm)

P[0,p,{6,3}][a_,b_,d_,e_,f_]:=(PL[a,b] PT[d,f] q[e])/(3 Sqrt[2] Sqrt[-1+dimM] qm)+(Sqrt[2] PL[a,b] PT[d,f] q[e])/(3 Sqrt[-1+dimM] qm)-(PL[a,b] PT[d,e] q[f])/(3 Sqrt[2] Sqrt[-1+dimM] qm)-(Sqrt[2] PL[a,b] PT[d,e] q[f])/(3 Sqrt[-1+dimM] qm)

P[0,p,{6,4}][a_,b_,d_,e_,f_]:=(PL[a,b] PL[e,f] q[d])/(3 qm)+(PL[a,b] PL[d,f] q[e])/(3 qm)+(PL[a,b] PL[d,e] q[f])/(3 qm)

P[0,p,{6,5}][a_,b_,e_,f_]:=(PL[a,b] PT[e,f])/Sqrt[-1+dimM]

P[0,p,{6,6}][a_,b_,e_,f_]:=PL[a,b] PL[e,f]
(*Sector1-*)
P[1,m,{1,1}][c_,a_,b_,d_,e_,f_]:=(PT[a,f] PT[b,c] PT[d,e])/(3 (1+dimM))+(PT[a,c] PT[b,f] PT[d,e])/(3 (1+dimM))+(PT[a,b] PT[c,f] PT[d,e])/(3 (1+dimM))+(PT[a,e] PT[b,c] PT[d,f])/(3 (1+dimM))+(PT[a,c] PT[b,e] PT[d,f])/(3 (1+dimM))+(PT[a,b] PT[c,e] PT[d,f])/(3 (1+dimM))+(PT[a,d] PT[b,c] PT[e,f])/(3 (1+dimM))+(PT[a,c] PT[b,d] PT[e,f])/(3 (1+dimM))+(PT[a,b] PT[c,d] PT[e,f])/(3 (1+dimM))
P[1,m,{1,2}][c_,a_,b_,d_,e_,f_]:=1/(3 Sqrt[2] (-2+dimM) (1+dimM))Sqrt[-2-dimM+dimM^2](PT[a,f] PT[b,c] PT[d,e]+PT[a,c] PT[b,f] PT[d,e]+PT[a,b] PT[c,f] PT[d,e]+PT[a,e] PT[b,c] PT[d,f]+PT[a,c] PT[b,e] PT[d,f]+PT[a,b] PT[c,e] PT[d,f]-2 PT[a,d] PT[b,c] PT[e,f]-2 PT[a,c] PT[b,d] PT[e,f]-2 PT[a,b] PT[c,d] PT[e,f])
P[1,m,{1,3}][c_,a_,b_,d_,e_,f_]:=(Sqrt[1/(-2-dimM+dimM^2)] PT[a,f] PT[b,c] PT[d,e])/Sqrt[6]+(Sqrt[1/(-2-dimM+dimM^2)] PT[a,c] PT[b,f] PT[d,e])/Sqrt[6]+(Sqrt[1/(-2-dimM+dimM^2)] PT[a,b] PT[c,f] PT[d,e])/Sqrt[6]-(Sqrt[1/(-2-dimM+dimM^2)] PT[a,e] PT[b,c] PT[d,f])/Sqrt[6]-(Sqrt[1/(-2-dimM+dimM^2)] PT[a,c] PT[b,e] PT[d,f])/Sqrt[6]-(Sqrt[1/(-2-dimM+dimM^2)] PT[a,b] PT[c,e] PT[d,f])/Sqrt[6]

P[1,m,{1,4}][c_,a_,b_,d_,e_,f_]:=(PL[e,f] PT[a,d] PT[b,c])/(3 Sqrt[1+dimM])+(PL[d,f] PT[a,e] PT[b,c])/(3 Sqrt[1+dimM])+(PL[d,e] PT[a,f] PT[b,c])/(3 Sqrt[1+dimM])+(PL[e,f] PT[a,c] PT[b,d])/(3 Sqrt[1+dimM])+(PL[d,f] PT[a,c] PT[b,e])/(3 Sqrt[1+dimM])+(PL[d,e] PT[a,c] PT[b,f])/(3 Sqrt[1+dimM])+(PL[e,f] PT[a,b] PT[c,d])/(3 Sqrt[1+dimM])+(PL[d,f] PT[a,b] PT[c,e])/(3 Sqrt[1+dimM])+(PL[d,e] PT[a,b] PT[c,f])/(3 Sqrt[1+dimM])

P[1,m,{1,5}][c_,a_,b_,d_,e_,f_]:=-((Sqrt[2] PL[e,f] PT[a,d] PT[b,c])/(3 Sqrt[1+dimM]))+(PL[d,f] PT[a,e] PT[b,c])/(3 Sqrt[2] Sqrt[1+dimM])+(PL[d,e] PT[a,f] PT[b,c])/(3 Sqrt[2] Sqrt[1+dimM])-(Sqrt[2] PL[e,f] PT[a,c] PT[b,d])/(3 Sqrt[1+dimM])+(PL[d,f] PT[a,c] PT[b,e])/(3 Sqrt[2] Sqrt[1+dimM])+(PL[d,e] PT[a,c] PT[b,f])/(3 Sqrt[2] Sqrt[1+dimM])-(Sqrt[2] PL[e,f] PT[a,b] PT[c,d])/(3 Sqrt[1+dimM])+(PL[d,f] PT[a,b] PT[c,e])/(3 Sqrt[2] Sqrt[1+dimM])+(PL[d,e] PT[a,b] PT[c,f])/(3 Sqrt[2] Sqrt[1+dimM])

P[1,m,{1,6}][c_,a_,b_,d_,e_,f_]:=-((Sqrt[1/(1+dimM)] PL[d,f] PT[a,e] PT[b,c])/Sqrt[6])+(Sqrt[1/(1+dimM)] PL[d,e] PT[a,f] PT[b,c])/Sqrt[6]-(Sqrt[1/(1+dimM)] PL[d,f] PT[a,c] PT[b,e])/Sqrt[6]+(Sqrt[1/(1+dimM)] PL[d,e] PT[a,c] PT[b,f])/Sqrt[6]-(Sqrt[1/(1+dimM)] PL[d,f] PT[a,b] PT[c,e])/Sqrt[6]+(Sqrt[1/(1+dimM)] PL[d,e] PT[a,b] PT[c,f])/Sqrt[6]

P[1,m,{1,7}][c_,a_,b_,e_,f_]:=(Sqrt[1/(1+dimM)] PT[a,f] PT[b,c] q[e])/(Sqrt[6] qm)+(Sqrt[1/(1+dimM)] PT[a,c] PT[b,f] q[e])/(Sqrt[6] qm)+(Sqrt[1/(1+dimM)] PT[a,b] PT[c,f] q[e])/(Sqrt[6] qm)+(Sqrt[1/(1+dimM)] PT[a,e] PT[b,c] q[f])/(Sqrt[6] qm)+(Sqrt[1/(1+dimM)] PT[a,c] PT[b,e] q[f])/(Sqrt[6] qm)+(Sqrt[1/(1+dimM)] PT[a,b] PT[c,e] q[f])/(Sqrt[6] qm)

P[1,m,{2,1}][c_,a_,b_,d_,e_,f_]:=1/(3 Sqrt[2] (-2+dimM) (1+dimM))Sqrt[-2-dimM+dimM^2] (PT[a,f] PT[b,c] PT[d,e]+PT[a,c] PT[b,f] PT[d,e]-2 PT[a,b] PT[c,f] PT[d,e]+PT[a,e] PT[b,c] PT[d,f]+PT[a,c] PT[b,e] PT[d,f]-2 PT[a,b] PT[c,e] PT[d,f]+PT[a,d] PT[b,c] PT[e,f]+PT[a,c] PT[b,d] PT[e,f]-2 PT[a,b] PT[c,d] PT[e,f])

P[1,m,{2,2}][c_,a_,b_,d_,e_,f_]:=(PT[a,f] PT[b,c] PT[d,e])/(6 (-2+dimM))+(PT[a,c] PT[b,f] PT[d,e])/(6 (-2+dimM))-(PT[a,b] PT[c,f] PT[d,e])/(3 (-2+dimM))+(PT[a,e] PT[b,c] PT[d,f])/(6 (-2+dimM))+(PT[a,c] PT[b,e] PT[d,f])/(6 (-2+dimM))-(PT[a,b] PT[c,e] PT[d,f])/(3 (-2+dimM))-(PT[a,d] PT[b,c] PT[e,f])/(3 (-2+dimM))-(PT[a,c] PT[b,d] PT[e,f])/(3 (-2+dimM))+(2 PT[a,b] PT[c,d] PT[e,f])/(3 (-2+dimM))

P[1,m,{2,3}][c_,a_,b_,d_,e_,f_]:=(PT[a,f] PT[b,c] PT[d,e])/(2 Sqrt[3] (-2+dimM))+(PT[a,c] PT[b,f] PT[d,e])/(2 Sqrt[3] (-2+dimM))-(PT[a,b] PT[c,f] PT[d,e])/(Sqrt[3] (-2+dimM))-(PT[a,e] PT[b,c] PT[d,f])/(2 Sqrt[3] (-2+dimM))-(PT[a,c] PT[b,e] PT[d,f])/(2 Sqrt[3] (-2+dimM))+(PT[a,b] PT[c,e] PT[d,f])/(Sqrt[3] (-2+dimM))

P[1,m,{2,4}][c_,a_,b_,d_,e_,f_]:=(PL[e,f] PT[a,d] PT[b,c])/(3 Sqrt[2] Sqrt[-2+dimM])+(PL[d,f] PT[a,e] PT[b,c])/(3 Sqrt[2] Sqrt[-2+dimM])+(PL[d,e] PT[a,f] PT[b,c])/(3 Sqrt[2] Sqrt[-2+dimM])+(PL[e,f] PT[a,c] PT[b,d])/(3 Sqrt[2] Sqrt[-2+dimM])+(PL[d,f] PT[a,c] PT[b,e])/(3 Sqrt[2] Sqrt[-2+dimM])+(PL[d,e] PT[a,c] PT[b,f])/(3 Sqrt[2] Sqrt[-2+dimM])-(Sqrt[2] PL[e,f] PT[a,b] PT[c,d])/(3 Sqrt[-2+dimM])-(Sqrt[2] PL[d,f] PT[a,b] PT[c,e])/(3 Sqrt[-2+dimM])-(Sqrt[2] PL[d,e] PT[a,b] PT[c,f])/(3 Sqrt[-2+dimM])

P[1,m,{2,5}][c_,a_,b_,d_,e_,f_]:=-((PL[e,f] PT[a,d] PT[b,c])/(3 Sqrt[-2+dimM]))+(PL[d,f] PT[a,e] PT[b,c])/(6 Sqrt[-2+dimM])+(PL[d,e] PT[a,f] PT[b,c])/(6 Sqrt[-2+dimM])-(PL[e,f] PT[a,c] PT[b,d])/(3 Sqrt[-2+dimM])+(PL[d,f] PT[a,c] PT[b,e])/(6 Sqrt[-2+dimM])+(PL[d,e] PT[a,c] PT[b,f])/(6 Sqrt[-2+dimM])+(2 PL[e,f] PT[a,b] PT[c,d])/(3 Sqrt[-2+dimM])-(PL[d,f] PT[a,b] PT[c,e])/(3 Sqrt[-2+dimM])-(PL[d,e] PT[a,b] PT[c,f])/(3 Sqrt[-2+dimM])

P[1,m,{2,6}][c_,a_,b_,d_,e_,f_]:=-((PL[d,f] PT[a,e] PT[b,c])/(2 Sqrt[3] Sqrt[-2+dimM]))+(PL[d,e] PT[a,f] PT[b,c])/(2 Sqrt[3] Sqrt[-2+dimM])-(PL[d,f] PT[a,c] PT[b,e])/(2 Sqrt[3] Sqrt[-2+dimM])+(PL[d,e] PT[a,c] PT[b,f])/(2 Sqrt[3] Sqrt[-2+dimM])+(PL[d,f] PT[a,b] PT[c,e])/(Sqrt[3] Sqrt[-2+dimM])-(PL[d,e] PT[a,b] PT[c,f])/(Sqrt[3] Sqrt[-2+dimM])

P[1,m,{2,7}][c_,a_,b_,e_,f_]:=(PT[a,f] PT[b,c] q[e])/(2 Sqrt[3] Sqrt[-2+dimM] qm)+(PT[a,c] PT[b,f] q[e])/(2 Sqrt[3] Sqrt[-2+dimM] qm)-(PT[a,b] PT[c,f] q[e])/(Sqrt[3] Sqrt[-2+dimM] qm)+(PT[a,e] PT[b,c] q[f])/(2 Sqrt[3] Sqrt[-2+dimM] qm)+(PT[a,c] PT[b,e] q[f])/(2 Sqrt[3] Sqrt[-2+dimM] qm)-(PT[a,b] PT[c,e] q[f])/(Sqrt[3] Sqrt[-2+dimM] qm)

P[1,m,{3,1}][c_,a_,b_,d_,e_,f_]:=-((Sqrt[1/(-2-dimM+dimM^2)] PT[a,f] PT[b,c] PT[d,e])/Sqrt[6])+(Sqrt[1/(-2-dimM+dimM^2)] PT[a,c] PT[b,f] PT[d,e])/Sqrt[6]-(Sqrt[1/(-2-dimM+dimM^2)] PT[a,e] PT[b,c] PT[d,f])/Sqrt[6]+(Sqrt[1/(-2-dimM+dimM^2)] PT[a,c] PT[b,e] PT[d,f])/Sqrt[6]-(Sqrt[1/(-2-dimM+dimM^2)] PT[a,d] PT[b,c] PT[e,f])/Sqrt[6]+(Sqrt[1/(-2-dimM+dimM^2)] PT[a,c] PT[b,d] PT[e,f])/Sqrt[6]

P[1,m,{3,2}][c_,a_,b_,d_,e_,f_]:=-((PT[a,f] PT[b,c] PT[d,e])/(2 Sqrt[3] (-2+dimM)))+(PT[a,c] PT[b,f] PT[d,e])/(2 Sqrt[3] (-2+dimM))-(PT[a,e] PT[b,c] PT[d,f])/(2 Sqrt[3] (-2+dimM))+(PT[a,c] PT[b,e] PT[d,f])/(2 Sqrt[3] (-2+dimM))+(PT[a,d] PT[b,c] PT[e,f])/(Sqrt[3] (-2+dimM))-(PT[a,c] PT[b,d] PT[e,f])/(Sqrt[3] (-2+dimM))

P[1,m,{3,3}][c_,a_,b_,d_,e_,f_]:=-((PT[a,f] PT[b,c] PT[d,e])/(2 (-2+dimM)))+(PT[a,c] PT[b,f] PT[d,e])/(2 (-2+dimM))+(PT[a,e] PT[b,c] PT[d,f])/(2 (-2+dimM))-(PT[a,c] PT[b,e] PT[d,f])/(2 (-2+dimM))

P[1,m,{3,4}][c_,a_,b_,d_,e_,f_]:=-((PL[e,f] PT[a,d] PT[b,c])/(Sqrt[6] Sqrt[-2+dimM]))-(PL[d,f] PT[a,e] PT[b,c])/(Sqrt[6] Sqrt[-2+dimM])-(PL[d,e] PT[a,f] PT[b,c])/(Sqrt[6] Sqrt[-2+dimM])+(PL[e,f] PT[a,c] PT[b,d])/(Sqrt[6] Sqrt[-2+dimM])+(PL[d,f] PT[a,c] PT[b,e])/(Sqrt[6] Sqrt[-2+dimM])+(PL[d,e] PT[a,c] PT[b,f])/(Sqrt[6] Sqrt[-2+dimM])

P[1,m,{3,5}][c_,a_,b_,d_,e_,f_]:=(PL[e,f] PT[a,d] PT[b,c])/(Sqrt[3] Sqrt[-2+dimM])-(PL[d,f] PT[a,e] PT[b,c])/(2 Sqrt[3] Sqrt[-2+dimM])-(PL[d,e] PT[a,f] PT[b,c])/(2 Sqrt[3] Sqrt[-2+dimM])-(PL[e,f] PT[a,c] PT[b,d])/(Sqrt[3] Sqrt[-2+dimM])+(PL[d,f] PT[a,c] PT[b,e])/(2 Sqrt[3] Sqrt[-2+dimM])+(PL[d,e] PT[a,c] PT[b,f])/(2 Sqrt[3] Sqrt[-2+dimM])

P[1,m,{3,6}][c_,a_,b_,d_,e_,f_]:=(PL[d,f] PT[a,e] PT[b,c])/(2 Sqrt[-2+dimM])-(PL[d,e] PT[a,f] PT[b,c])/(2 Sqrt[-2+dimM])-(PL[d,f] PT[a,c] PT[b,e])/(2 Sqrt[-2+dimM])+(PL[d,e] PT[a,c] PT[b,f])/(2 Sqrt[-2+dimM])

P[1,m,{3,7}][c_,a_,b_,e_,f_]:=-((PT[a,f] PT[b,c] q[e])/(2 Sqrt[-2+dimM] qm))+(PT[a,c] PT[b,f] q[e])/(2 Sqrt[-2+dimM] qm)-(PT[a,e] PT[b,c] q[f])/(2 Sqrt[-2+dimM] qm)+(PT[a,c] PT[b,e] q[f])/(2 Sqrt[-2+dimM] qm)

P[1,m,{4,1}][c_,a_,b_,d_,e_,f_]:=(PL[b,c] PT[a,f] PT[d,e])/(3 Sqrt[1+dimM])+(PL[a,c] PT[b,f] PT[d,e])/(3 Sqrt[1+dimM])+(PL[a,b] PT[c,f] PT[d,e])/(3 Sqrt[1+dimM])+(PL[b,c] PT[a,e] PT[d,f])/(3 Sqrt[1+dimM])+(PL[a,c] PT[b,e] PT[d,f])/(3 Sqrt[1+dimM])+(PL[a,b] PT[c,e] PT[d,f])/(3 Sqrt[1+dimM])+(PL[b,c] PT[a,d] PT[e,f])/(3 Sqrt[1+dimM])+(PL[a,c] PT[b,d] PT[e,f])/(3 Sqrt[1+dimM])+(PL[a,b] PT[c,d] PT[e,f])/(3 Sqrt[1+dimM])

P[1,m,{4,2}][c_,a_,b_,d_,e_,f_]:=(PL[b,c] PT[a,f] PT[d,e])/(3 Sqrt[2] Sqrt[-2+dimM])+(PL[a,c] PT[b,f] PT[d,e])/(3 Sqrt[2] Sqrt[-2+dimM])+(PL[a,b] PT[c,f] PT[d,e])/(3 Sqrt[2] Sqrt[-2+dimM])+(PL[b,c] PT[a,e] PT[d,f])/(3 Sqrt[2] Sqrt[-2+dimM])+(PL[a,c] PT[b,e] PT[d,f])/(3 Sqrt[2] Sqrt[-2+dimM])+(PL[a,b] PT[c,e] PT[d,f])/(3 Sqrt[2] Sqrt[-2+dimM])-(Sqrt[2] PL[b,c] PT[a,d] PT[e,f])/(3 Sqrt[-2+dimM])-(Sqrt[2] PL[a,c] PT[b,d] PT[e,f])/(3 Sqrt[-2+dimM])-(Sqrt[2] PL[a,b] PT[c,d] PT[e,f])/(3 Sqrt[-2+dimM])

P[1,m,{4,3}][c_,a_,b_,d_,e_,f_]:=(PL[b,c] PT[a,f] PT[d,e])/(Sqrt[6] Sqrt[-2+dimM])+(PL[a,c] PT[b,f] PT[d,e])/(Sqrt[6] Sqrt[-2+dimM])+(PL[a,b] PT[c,f] PT[d,e])/(Sqrt[6] Sqrt[-2+dimM])-(PL[b,c] PT[a,e] PT[d,f])/(Sqrt[6] Sqrt[-2+dimM])-(PL[a,c] PT[b,e] PT[d,f])/(Sqrt[6] Sqrt[-2+dimM])-(PL[a,b] PT[c,e] PT[d,f])/(Sqrt[6] Sqrt[-2+dimM])

P[1,m,{4,4}][c_,a_,b_,d_,e_,f_]:=1/6 PL[b,f] PL[c,e] PT[a,d]+1/6 PL[b,e] PL[c,f] PT[a,d]+1/6 PL[b,f] PL[c,d] PT[a,e]+1/6 PL[b,d] PL[c,f] PT[a,e]+1/6 PL[b,e] PL[c,d] PT[a,f]+1/6 PL[b,d] PL[c,e] PT[a,f]+1/6 PL[a,f] PL[c,e] PT[b,d]+1/6 PL[a,e] PL[c,f] PT[b,d]+1/6 PL[a,f] PL[c,d] PT[b,e]+1/6 PL[a,d] PL[c,f] PT[b,e]+1/6 PL[a,e] PL[c,d] PT[b,f]+1/6 PL[a,d] PL[c,e] PT[b,f]+1/6 PL[a,f] PL[b,e] PT[c,d]+1/6 PL[a,e] PL[b,f] PT[c,d]+1/6 PL[a,f] PL[b,d] PT[c,e]+1/6 PL[a,d] PL[b,f] PT[c,e]+1/6 PL[a,e] PL[b,d] PT[c,f]+1/6 PL[a,d] PL[b,e] PT[c,f]

P[1,m,{4,5}][c_,a_,b_,d_,e_,f_]:=-((PL[b,f] PL[c,e] PT[a,d])/(3 Sqrt[2]))-(PL[b,e] PL[c,f] PT[a,d])/(3 Sqrt[2])+(PL[b,f] PL[c,d] PT[a,e])/(6 Sqrt[2])+(PL[b,d] PL[c,f] PT[a,e])/(6 Sqrt[2])+(PL[b,e] PL[c,d] PT[a,f])/(6 Sqrt[2])+(PL[b,d] PL[c,e] PT[a,f])/(6 Sqrt[2])-(PL[a,f] PL[c,e] PT[b,d])/(3 Sqrt[2])-(PL[a,e] PL[c,f] PT[b,d])/(3 Sqrt[2])+(PL[a,f] PL[c,d] PT[b,e])/(6 Sqrt[2])+(PL[a,d] PL[c,f] PT[b,e])/(6 Sqrt[2])+(PL[a,e] PL[c,d] PT[b,f])/(6 Sqrt[2])+(PL[a,d] PL[c,e] PT[b,f])/(6 Sqrt[2])-(PL[a,f] PL[b,e] PT[c,d])/(3 Sqrt[2])-(PL[a,e] PL[b,f] PT[c,d])/(3 Sqrt[2])+(PL[a,f] PL[b,d] PT[c,e])/(6 Sqrt[2])+(PL[a,d] PL[b,f] PT[c,e])/(6 Sqrt[2])+(PL[a,e] PL[b,d] PT[c,f])/(6 Sqrt[2])+(PL[a,d] PL[b,e] PT[c,f])/(6 Sqrt[2])

P[1,m,{4,6}][c_,a_,b_,d_,e_,f_]:=-((PL[b,f] PL[c,d] PT[a,e])/(2 Sqrt[6]))-(PL[b,d] PL[c,f] PT[a,e])/(2 Sqrt[6])+(PL[b,e] PL[c,d] PT[a,f])/(2 Sqrt[6])+(PL[b,d] PL[c,e] PT[a,f])/(2 Sqrt[6])-(PL[a,f] PL[c,d] PT[b,e])/(2 Sqrt[6])-(PL[a,d] PL[c,f] PT[b,e])/(2 Sqrt[6])+(PL[a,e] PL[c,d] PT[b,f])/(2 Sqrt[6])+(PL[a,d] PL[c,e] PT[b,f])/(2 Sqrt[6])-(PL[a,f] PL[b,d] PT[c,e])/(2 Sqrt[6])-(PL[a,d] PL[b,f] PT[c,e])/(2 Sqrt[6])+(PL[a,e] PL[b,d] PT[c,f])/(2 Sqrt[6])+(PL[a,d] PL[b,e] PT[c,f])/(2 Sqrt[6])

P[1,m,{4,7}][c_,a_,b_,e_,f_]:=(PL[b,c] PT[a,f] q[e])/(Sqrt[6] qm)+(PL[a,c] PT[b,f] q[e])/(Sqrt[6] qm)+(PL[a,b] PT[c,f] q[e])/(Sqrt[6] qm)+(PL[b,c] PT[a,e] q[f])/(Sqrt[6] qm)+(PL[a,c] PT[b,e] q[f])/(Sqrt[6] qm)+(PL[a,b] PT[c,e] q[f])/(Sqrt[6] qm)

P[1,m,{5,1}][c_,a_,b_,d_,e_,f_]:=(PL[b,c] PT[a,f] PT[d,e])/(3 Sqrt[2] Sqrt[1+dimM])+(PL[a,c] PT[b,f] PT[d,e])/(3 Sqrt[2] Sqrt[1+dimM])-(Sqrt[2] PL[a,b] PT[c,f] PT[d,e])/(3 Sqrt[1+dimM])+(PL[b,c] PT[a,e] PT[d,f])/(3 Sqrt[2] Sqrt[1+dimM])+(PL[a,c] PT[b,e] PT[d,f])/(3 Sqrt[2] Sqrt[1+dimM])-(Sqrt[2] PL[a,b] PT[c,e] PT[d,f])/(3 Sqrt[1+dimM])+(PL[b,c] PT[a,d] PT[e,f])/(3 Sqrt[2] Sqrt[1+dimM])+(PL[a,c] PT[b,d] PT[e,f])/(3 Sqrt[2] Sqrt[1+dimM])-(Sqrt[2] PL[a,b] PT[c,d] PT[e,f])/(3 Sqrt[1+dimM])

P[1,m,{5,2}][c_,a_,b_,d_,e_,f_]:=(PL[b,c] PT[a,f] PT[d,e])/(6 Sqrt[-2+dimM])+(PL[a,c] PT[b,f] PT[d,e])/(6 Sqrt[-2+dimM])-(PL[a,b] PT[c,f] PT[d,e])/(3 Sqrt[-2+dimM])+(PL[b,c] PT[a,e] PT[d,f])/(6 Sqrt[-2+dimM])+(PL[a,c] PT[b,e] PT[d,f])/(6 Sqrt[-2+dimM])-(PL[a,b] PT[c,e] PT[d,f])/(3 Sqrt[-2+dimM])-(PL[b,c] PT[a,d] PT[e,f])/(3 Sqrt[-2+dimM])-(PL[a,c] PT[b,d] PT[e,f])/(3 Sqrt[-2+dimM])+(2 PL[a,b] PT[c,d] PT[e,f])/(3 Sqrt[-2+dimM])

P[1,m,{5,3}][c_,a_,b_,d_,e_,f_]:=(PL[b,c] PT[a,f] PT[d,e])/(2 Sqrt[3] Sqrt[-2+dimM])+(PL[a,c] PT[b,f] PT[d,e])/(2 Sqrt[3] Sqrt[-2+dimM])-(PL[a,b] PT[c,f] PT[d,e])/(Sqrt[3] Sqrt[-2+dimM])-(PL[b,c] PT[a,e] PT[d,f])/(2 Sqrt[3] Sqrt[-2+dimM])-(PL[a,c] PT[b,e] PT[d,f])/(2 Sqrt[3] Sqrt[-2+dimM])+(PL[a,b] PT[c,e] PT[d,f])/(Sqrt[3] Sqrt[-2+dimM])

P[1,m,{5,4}][c_,a_,b_,d_,e_,f_]:=(PL[b,f] PL[c,e] PT[a,d])/(6 Sqrt[2])+(PL[b,e] PL[c,f] PT[a,d])/(6 Sqrt[2])+(PL[b,f] PL[c,d] PT[a,e])/(6 Sqrt[2])+(PL[b,d] PL[c,f] PT[a,e])/(6 Sqrt[2])+(PL[b,e] PL[c,d] PT[a,f])/(6 Sqrt[2])+(PL[b,d] PL[c,e] PT[a,f])/(6 Sqrt[2])+(PL[a,f] PL[c,e] PT[b,d])/(6 Sqrt[2])+(PL[a,e] PL[c,f] PT[b,d])/(6 Sqrt[2])+(PL[a,f] PL[c,d] PT[b,e])/(6 Sqrt[2])+(PL[a,d] PL[c,f] PT[b,e])/(6 Sqrt[2])+(PL[a,e] PL[c,d] PT[b,f])/(6 Sqrt[2])+(PL[a,d] PL[c,e] PT[b,f])/(6 Sqrt[2])-(PL[a,f] PL[b,e] PT[c,d])/(3 Sqrt[2])-(PL[a,e] PL[b,f] PT[c,d])/(3 Sqrt[2])-(PL[a,f] PL[b,d] PT[c,e])/(3 Sqrt[2])-(PL[a,d] PL[b,f] PT[c,e])/(3 Sqrt[2])-(PL[a,e] PL[b,d] PT[c,f])/(3 Sqrt[2])-(PL[a,d] PL[b,e] PT[c,f])/(3 Sqrt[2])

P[1,m,{5,5}][c_,a_,b_,d_,e_,f_]:=-(1/6) PL[b,f] PL[c,e] PT[a,d]-1/6 PL[b,e] PL[c,f] PT[a,d]+1/3 PL[b,f] PL[c,d] PT[a,e]-1/6 PL[b,d] PL[c,f] PT[a,e]+1/3 PL[b,e] PL[c,d] PT[a,f]-1/6 PL[b,d] PL[c,e] PT[a,f]-1/6 PL[a,f] PL[c,e] PT[b,d]-1/6 PL[a,e] PL[c,f] PT[b,d]+1/3 PL[a,f] PL[c,d] PT[b,e]-1/6 PL[a,d] PL[c,f] PT[b,e]+1/3 PL[a,e] PL[c,d] PT[b,f]-1/6 PL[a,d] PL[c,e] PT[b,f]+1/3 PL[a,f] PL[b,e] PT[c,d]+1/3 PL[a,e] PL[b,f] PT[c,d]-1/6 PL[a,f] PL[b,d] PT[c,e]-1/6 PL[a,d] PL[b,f] PT[c,e]-1/6 PL[a,e] PL[b,d] PT[c,f]-1/6 PL[a,d] PL[b,e] PT[c,f]

P[1,m,{5,6}][c_,a_,b_,d_,e_,f_]:=-((PL[b,f] PL[c,e] PT[a,d])/(2 Sqrt[3]))+(PL[b,e] PL[c,f] PT[a,d])/(2 Sqrt[3])-(PL[b,f] PL[c,d] PT[a,e])/(2 Sqrt[3])+(PL[b,e] PL[c,d] PT[a,f])/(2 Sqrt[3])-(PL[a,f] PL[c,e] PT[b,d])/(2 Sqrt[3])+(PL[a,e] PL[c,f] PT[b,d])/(2 Sqrt[3])-(PL[a,f] PL[c,d] PT[b,e])/(2 Sqrt[3])+(PL[a,e] PL[c,d] PT[b,f])/(2 Sqrt[3])+(PL[a,f] PL[b,d] PT[c,e])/(2 Sqrt[3])+(PL[a,d] PL[b,f] PT[c,e])/(2 Sqrt[3])-(PL[a,e] PL[b,d] PT[c,f])/(2 Sqrt[3])-(PL[a,d] PL[b,e] PT[c,f])/(2 Sqrt[3])

P[1,m,{5,7}][c_,a_,b_,e_,f_]:=(PL[b,c] PT[a,f] q[e])/(2 Sqrt[3] qm)+(PL[a,c] PT[b,f] q[e])/(2 Sqrt[3] qm)-(PL[a,b] PT[c,f] q[e])/(Sqrt[3] qm)+(PL[b,c] PT[a,e] q[f])/(2 Sqrt[3] qm)+(PL[a,c] PT[b,e] q[f])/(2 Sqrt[3] qm)-(PL[a,b] PT[c,e] q[f])/(Sqrt[3] qm)

P[1,m,{6,1}][c_,a_,b_,d_,e_,f_]:=-((Sqrt[1/(1+dimM)] PL[b,c] PT[a,f] PT[d,e])/Sqrt[6])+(Sqrt[1/(1+dimM)] PL[a,c] PT[b,f] PT[d,e])/Sqrt[6]-(Sqrt[1/(1+dimM)] PL[b,c] PT[a,e] PT[d,f])/Sqrt[6]+(Sqrt[1/(1+dimM)] PL[a,c] PT[b,e] PT[d,f])/Sqrt[6]-(Sqrt[1/(1+dimM)] PL[b,c] PT[a,d] PT[e,f])/Sqrt[6]+(Sqrt[1/(1+dimM)] PL[a,c] PT[b,d] PT[e,f])/Sqrt[6]

P[1,m,{6,2}][c_,a_,b_,d_,e_,f_]:=-((PL[b,c] PT[a,f] PT[d,e])/(2 Sqrt[3] Sqrt[-2+dimM]))+(PL[a,c] PT[b,f] PT[d,e])/(2 Sqrt[3] Sqrt[-2+dimM])-(PL[b,c] PT[a,e] PT[d,f])/(2 Sqrt[3] Sqrt[-2+dimM])+(PL[a,c] PT[b,e] PT[d,f])/(2 Sqrt[3] Sqrt[-2+dimM])+(PL[b,c] PT[a,d] PT[e,f])/(Sqrt[3] Sqrt[-2+dimM])-(PL[a,c] PT[b,d] PT[e,f])/(Sqrt[3] Sqrt[-2+dimM])

P[1,m,{6,3}][c_,a_,b_,d_,e_,f_]:=-((PL[b,c] PT[a,f] PT[d,e])/(2 Sqrt[-2+dimM]))+(PL[a,c] PT[b,f] PT[d,e])/(2 Sqrt[-2+dimM])+(PL[b,c] PT[a,e] PT[d,f])/(2 Sqrt[-2+dimM])-(PL[a,c] PT[b,e] PT[d,f])/(2 Sqrt[-2+dimM])

P[1,m,{6,4}][c_,a_,b_,d_,e_,f_]:=-((PL[b,f] PL[c,e] PT[a,d])/(2 Sqrt[6]))-(PL[b,e] PL[c,f] PT[a,d])/(2 Sqrt[6])-(PL[b,f] PL[c,d] PT[a,e])/(2 Sqrt[6])-(PL[b,d] PL[c,f] PT[a,e])/(2 Sqrt[6])-(PL[b,e] PL[c,d] PT[a,f])/(2 Sqrt[6])-(PL[b,d] PL[c,e] PT[a,f])/(2 Sqrt[6])+(PL[a,f] PL[c,e] PT[b,d])/(2 Sqrt[6])+(PL[a,e] PL[c,f] PT[b,d])/(2 Sqrt[6])+(PL[a,f] PL[c,d] PT[b,e])/(2 Sqrt[6])+(PL[a,d] PL[c,f] PT[b,e])/(2 Sqrt[6])+(PL[a,e] PL[c,d] PT[b,f])/(2 Sqrt[6])+(PL[a,d] PL[c,e] PT[b,f])/(2 Sqrt[6])

P[1,m,{6,5}][c_,a_,b_,d_,e_,f_]:=(PL[b,f] PL[c,e] PT[a,d])/(2 Sqrt[3])+(PL[b,e] PL[c,f] PT[a,d])/(2 Sqrt[3])-(PL[b,f] PL[c,d] PT[a,e])/(2 Sqrt[3])-(PL[b,e] PL[c,d] PT[a,f])/(2 Sqrt[3])-(PL[a,f] PL[c,e] PT[b,d])/(2 Sqrt[3])-(PL[a,e] PL[c,f] PT[b,d])/(2 Sqrt[3])+(PL[a,f] PL[c,d] PT[b,e])/(2 Sqrt[3])+(PL[a,e] PL[c,d] PT[b,f])/(2 Sqrt[3])+(PL[a,f] PL[b,d] PT[c,e])/(2 Sqrt[3])-(PL[a,d] PL[b,f] PT[c,e])/(2 Sqrt[3])+(PL[a,e] PL[b,d] PT[c,f])/(2 Sqrt[3])-(PL[a,d] PL[b,e] PT[c,f])/(2 Sqrt[3])

P[1,m,{6,6}][c_,a_,b_,d_,e_,f_]:=1/6 PL[b,f] PL[c,e] PT[a,d]-1/6 PL[b,e] PL[c,f] PT[a,d]+1/3 PL[b,f] PL[c,d] PT[a,e]+1/6 PL[b,d] PL[c,f] PT[a,e]-1/3 PL[b,e] PL[c,d] PT[a,f]-1/6 PL[b,d] PL[c,e] PT[a,f]-1/6 PL[a,f] PL[c,e] PT[b,d]+1/6 PL[a,e] PL[c,f] PT[b,d]-1/3 PL[a,f] PL[c,d] PT[b,e]-1/6 PL[a,d] PL[c,f] PT[b,e]+1/3 PL[a,e] PL[c,d] PT[b,f]+1/6 PL[a,d] PL[c,e] PT[b,f]-1/3 PL[a,f] PL[b,e] PT[c,d]+1/3 PL[a,e] PL[b,f] PT[c,d]-1/6 PL[a,f] PL[b,d] PT[c,e]+1/6 PL[a,d] PL[b,f] PT[c,e]+1/6 PL[a,e] PL[b,d] PT[c,f]-1/6 PL[a,d] PL[b,e] PT[c,f]

P[1,m,{6,7}][c_,a_,b_,e_,f_]:=-((PL[b,c] PT[a,f] q[e])/(2 qm))+(PL[a,c] PT[b,f] q[e])/(2 qm)-(PL[b,c] PT[a,e] q[f])/(2 qm)+(PL[a,c] PT[b,e] q[f])/(2 qm)

P[1,m,{7,1}][a_,b_,d_,e_,f_]:=(Sqrt[1/(1+dimM)] PT[b,f] PT[d,e] q[a])/(Sqrt[6] qm)+(Sqrt[1/(1+dimM)] PT[b,e] PT[d,f] q[a])/(Sqrt[6] qm)+(Sqrt[1/(1+dimM)] PT[b,d] PT[e,f] q[a])/(Sqrt[6] qm)+(Sqrt[1/(1+dimM)] PT[a,f] PT[d,e] q[b])/(Sqrt[6] qm)+(Sqrt[1/(1+dimM)] PT[a,e] PT[d,f] q[b])/(Sqrt[6] qm)+(Sqrt[1/(1+dimM)] PT[a,d] PT[e,f] q[b])/(Sqrt[6] qm)

P[1,m,{7,2}][a_,b_,d_,e_,f_]:=(PT[b,f] PT[d,e] q[a])/(2 Sqrt[3] Sqrt[-2+dimM] qm)+(PT[b,e] PT[d,f] q[a])/(2 Sqrt[3] Sqrt[-2+dimM] qm)-(PT[b,d] PT[e,f] q[a])/(Sqrt[3] Sqrt[-2+dimM] qm)+(PT[a,f] PT[d,e] q[b])/(2 Sqrt[3] Sqrt[-2+dimM] qm)+(PT[a,e] PT[d,f] q[b])/(2 Sqrt[3] Sqrt[-2+dimM] qm)-(PT[a,d] PT[e,f] q[b])/(Sqrt[3] Sqrt[-2+dimM] qm)

P[1,m,{7,3}][a_,b_,d_,e_,f_]:=(PT[b,f] PT[d,e] q[a])/(2 Sqrt[-2+dimM] qm)-(PT[b,e] PT[d,f] q[a])/(2 Sqrt[-2+dimM] qm)+(PT[a,f] PT[d,e] q[b])/(2 Sqrt[-2+dimM] qm)-(PT[a,e] PT[d,f] q[b])/(2 Sqrt[-2+dimM] qm)

P[1,m,{7,4}][a_,b_,d_,e_,f_]:=(PL[e,f] PT[b,d] q[a])/(Sqrt[6] qm)+(PL[d,f] PT[b,e] q[a])/(Sqrt[6] qm)+(PL[d,e] PT[b,f] q[a])/(Sqrt[6] qm)+(PL[e,f] PT[a,d] q[b])/(Sqrt[6] qm)+(PL[d,f] PT[a,e] q[b])/(Sqrt[6] qm)+(PL[d,e] PT[a,f] q[b])/(Sqrt[6] qm)

P[1,m,{7,5}][a_,b_,d_,e_,f_]:=-((PL[e,f] PT[b,d] q[a])/(Sqrt[3] qm))+(PL[d,f] PT[b,e] q[a])/(2 Sqrt[3] qm)+(PL[d,e] PT[b,f] q[a])/(2 Sqrt[3] qm)-(PL[e,f] PT[a,d] q[b])/(Sqrt[3] qm)+(PL[d,f] PT[a,e] q[b])/(2 Sqrt[3] qm)+(PL[d,e] PT[a,f] q[b])/(2 Sqrt[3] qm)

P[1,m,{7,6}][a_,b_,d_,e_,f_]:=-((PL[d,f] PT[b,e] q[a])/(2 qm))+(PL[d,e] PT[b,f] q[a])/(2 qm)-(PL[d,f] PT[a,e] q[b])/(2 qm)+(PL[d,e] PT[a,f] q[b])/(2 qm)

P[1,m,{7,7}][a_,b_,e_,f_]:=1/2 (PL[b,f] PT[a,e]+PL[b,e] PT[a,f]+PL[a,f] PT[b,e]+PL[a,e] PT[b,f])

(*Sector1+*)
P[1,p,{1,1}][c_,a_,b_,d_,e_,f_]:=1/4 PL[b,f] PT[a,e] PT[c,d]+1/4 PL[b,e] PT[a,f] PT[c,d]+1/4 PL[a,f] PT[b,e] PT[c,d]+1/4 PL[a,e] PT[b,f] PT[c,d]-1/4 PL[b,f] PT[a,d] PT[c,e]-1/4 PL[a,f] PT[b,d] PT[c,e]-1/4 PL[b,e] PT[a,d] PT[c,f]-1/4 PL[a,e] PT[b,d] PT[c,f]
P[1,p,{1,2}][c_,a_,b_,d_,e_,f_]:=(PL[b,f] PT[a,e] PT[c,d])/(4 Sqrt[3])-(PL[b,e] PT[a,f] PT[c,d])/(4 Sqrt[3])+(PL[a,f] PT[b,e] PT[c,d])/(4 Sqrt[3])-(PL[a,e] PT[b,f] PT[c,d])/(4 Sqrt[3])-(PL[b,f] PT[a,d] PT[c,e])/(4 Sqrt[3])-(PL[b,d] PT[a,f] PT[c,e])/(2 Sqrt[3])-(PL[a,f] PT[b,d] PT[c,e])/(4 Sqrt[3])-(PL[a,d] PT[b,f] PT[c,e])/(2 Sqrt[3])+(PL[b,e] PT[a,d] PT[c,f])/(4 Sqrt[3])+(PL[b,d] PT[a,e] PT[c,f])/(2 Sqrt[3])+(PL[a,e] PT[b,d] PT[c,f])/(4 Sqrt[3])+(PL[a,d] PT[b,e] PT[c,f])/(2 Sqrt[3])
P[1,p,{1,3}][c_,a_,b_,d_,e_,f_]:=(PL[b,f] PT[a,e] PT[c,d])/(2 Sqrt[6])-(PL[b,e] PT[a,f] PT[c,d])/(2 Sqrt[6])+(PL[a,f] PT[b,e] PT[c,d])/(2 Sqrt[6])-(PL[a,e] PT[b,f] PT[c,d])/(2 Sqrt[6])-(PL[b,f] PT[a,d] PT[c,e])/(2 Sqrt[6])+(PL[b,d] PT[a,f] PT[c,e])/(2 Sqrt[6])-(PL[a,f] PT[b,d] PT[c,e])/(2 Sqrt[6])+(PL[a,d] PT[b,f] PT[c,e])/(2 Sqrt[6])+(PL[b,e] PT[a,d] PT[c,f])/(2 Sqrt[6])-(PL[b,d] PT[a,e] PT[c,f])/(2 Sqrt[6])+(PL[a,e] PT[b,d] PT[c,f])/(2 Sqrt[6])-(PL[a,d] PT[b,e] PT[c,f])/(2 Sqrt[6])
P[1,p,{2,1}][c_,a_,b_,d_,e_,f_]:=(PL[c,f] PT[a,e] PT[b,d])/(2 Sqrt[3])+(PL[c,e] PT[a,f] PT[b,d])/(2 Sqrt[3])-(PL[c,f] PT[a,d] PT[b,e])/(2 Sqrt[3])-(PL[c,e] PT[a,d] PT[b,f])/(2 Sqrt[3])+(PL[b,f] PT[a,e] PT[c,d])/(4 Sqrt[3])+(PL[b,e] PT[a,f] PT[c,d])/(4 Sqrt[3])-(PL[a,f] PT[b,e] PT[c,d])/(4 Sqrt[3])-(PL[a,e] PT[b,f] PT[c,d])/(4 Sqrt[3])-(PL[b,f] PT[a,d] PT[c,e])/(4 Sqrt[3])+(PL[a,f] PT[b,d] PT[c,e])/(4 Sqrt[3])-(PL[b,e] PT[a,d] PT[c,f])/(4 Sqrt[3])+(PL[a,e] PT[b,d] PT[c,f])/(4 Sqrt[3])
P[1,p,{2,2}][c_,a_,b_,d_,e_,f_]:=1/6 PL[c,f] PT[a,e] PT[b,d]-1/6 PL[c,e] PT[a,f] PT[b,d]-1/6 PL[c,f] PT[a,d] PT[b,e]-1/3 PL[c,d] PT[a,f] PT[b,e]+1/6 PL[c,e] PT[a,d] PT[b,f]+1/3 PL[c,d] PT[a,e] PT[b,f]+1/12 PL[b,f] PT[a,e] PT[c,d]-1/12 PL[b,e] PT[a,f] PT[c,d]-1/12 PL[a,f] PT[b,e] PT[c,d]+1/12 PL[a,e] PT[b,f] PT[c,d]-1/12 PL[b,f] PT[a,d] PT[c,e]-1/6 PL[b,d] PT[a,f] PT[c,e]+1/12 PL[a,f] PT[b,d] PT[c,e]+1/6 PL[a,d] PT[b,f] PT[c,e]+1/12 PL[b,e] PT[a,d] PT[c,f]+1/6 PL[b,d] PT[a,e] PT[c,f]-1/12 PL[a,e] PT[b,d] PT[c,f]-1/6 PL[a,d] PT[b,e] PT[c,f]
P[1,p,{2,3}][c_,a_,b_,d_,e_,f_]:=(PL[c,f] PT[a,e] PT[b,d])/(3 Sqrt[2])-(PL[c,e] PT[a,f] PT[b,d])/(3 Sqrt[2])-(PL[c,f] PT[a,d] PT[b,e])/(3 Sqrt[2])+(PL[c,d] PT[a,f] PT[b,e])/(3 Sqrt[2])+(PL[c,e] PT[a,d] PT[b,f])/(3 Sqrt[2])-(PL[c,d] PT[a,e] PT[b,f])/(3 Sqrt[2])+(PL[b,f] PT[a,e] PT[c,d])/(6 Sqrt[2])-(PL[b,e] PT[a,f] PT[c,d])/(6 Sqrt[2])-(PL[a,f] PT[b,e] PT[c,d])/(6 Sqrt[2])+(PL[a,e] PT[b,f] PT[c,d])/(6 Sqrt[2])-(PL[b,f] PT[a,d] PT[c,e])/(6 Sqrt[2])+(PL[b,d] PT[a,f] PT[c,e])/(6 Sqrt[2])+(PL[a,f] PT[b,d] PT[c,e])/(6 Sqrt[2])-(PL[a,d] PT[b,f] PT[c,e])/(6 Sqrt[2])+(PL[b,e] PT[a,d] PT[c,f])/(6 Sqrt[2])-(PL[b,d] PT[a,e] PT[c,f])/(6 Sqrt[2])-(PL[a,e] PT[b,d] PT[c,f])/(6 Sqrt[2])+(PL[a,d] PT[b,e] PT[c,f])/(6 Sqrt[2])
P[1,p,{3,1}][c_,a_,b_,d_,e_,f_]:=-((PL[c,f] PT[a,e] PT[b,d])/(2 Sqrt[6]))-(PL[c,e] PT[a,f] PT[b,d])/(2 Sqrt[6])+(PL[c,f] PT[a,d] PT[b,e])/(2 Sqrt[6])+(PL[c,e] PT[a,d] PT[b,f])/(2 Sqrt[6])+(PL[b,f] PT[a,e] PT[c,d])/(2 Sqrt[6])+(PL[b,e] PT[a,f] PT[c,d])/(2 Sqrt[6])-(PL[a,f] PT[b,e] PT[c,d])/(2 Sqrt[6])-(PL[a,e] PT[b,f] PT[c,d])/(2 Sqrt[6])-(PL[b,f] PT[a,d] PT[c,e])/(2 Sqrt[6])+(PL[a,f] PT[b,d] PT[c,e])/(2 Sqrt[6])-(PL[b,e] PT[a,d] PT[c,f])/(2 Sqrt[6])+(PL[a,e] PT[b,d] PT[c,f])/(2 Sqrt[6])
P[1,p,{3,2}][c_,a_,b_,d_,e_,f_]:=-((PL[c,f] PT[a,e] PT[b,d])/(6 Sqrt[2]))+(PL[c,e] PT[a,f] PT[b,d])/(6 Sqrt[2])+(PL[c,f] PT[a,d] PT[b,e])/(6 Sqrt[2])+(PL[c,d] PT[a,f] PT[b,e])/(3 Sqrt[2])-(PL[c,e] PT[a,d] PT[b,f])/(6 Sqrt[2])-(PL[c,d] PT[a,e] PT[b,f])/(3 Sqrt[2])+(PL[b,f] PT[a,e] PT[c,d])/(6 Sqrt[2])-(PL[b,e] PT[a,f] PT[c,d])/(6 Sqrt[2])-(PL[a,f] PT[b,e] PT[c,d])/(6 Sqrt[2])+(PL[a,e] PT[b,f] PT[c,d])/(6 Sqrt[2])-(PL[b,f] PT[a,d] PT[c,e])/(6 Sqrt[2])-(PL[b,d] PT[a,f] PT[c,e])/(3 Sqrt[2])+(PL[a,f] PT[b,d] PT[c,e])/(6 Sqrt[2])+(PL[a,d] PT[b,f] PT[c,e])/(3 Sqrt[2])+(PL[b,e] PT[a,d] PT[c,f])/(6 Sqrt[2])+(PL[b,d] PT[a,e] PT[c,f])/(3 Sqrt[2])-(PL[a,e] PT[b,d] PT[c,f])/(6 Sqrt[2])-(PL[a,d] PT[b,e] PT[c,f])/(3 Sqrt[2])
P[1,p,{3,3}][c_,a_,b_,d_,e_,f_]:=-(1/6) PL[c,f] PT[a,e] PT[b,d]+1/6 PL[c,e] PT[a,f] PT[b,d]+1/6 PL[c,f] PT[a,d] PT[b,e]-1/6 PL[c,d] PT[a,f] PT[b,e]-1/6 PL[c,e] PT[a,d] PT[b,f]+1/6 PL[c,d] PT[a,e] PT[b,f]+1/6 PL[b,f] PT[a,e] PT[c,d]-1/6 PL[b,e] PT[a,f] PT[c,d]-1/6 PL[a,f] PT[b,e] PT[c,d]+1/6 PL[a,e] PT[b,f] PT[c,d]-1/6 PL[b,f] PT[a,d] PT[c,e]+1/6 PL[b,d] PT[a,f] PT[c,e]+1/6 PL[a,f] PT[b,d] PT[c,e]-1/6 PL[a,d] PT[b,f] PT[c,e]+1/6 PL[b,e] PT[a,d] PT[c,f]-1/6 PL[b,d] PT[a,e] PT[c,f]-1/6 PL[a,e] PT[b,d] PT[c,f]+1/6 PL[a,d] PT[b,e] PT[c,f]
(*Sector2-*)
P[2,m,{1,1}][c_,a_,b_,d_,e_,f_]:=1/3 PT[a,f] PT[b,e] PT[c,d]+1/3 PT[a,e] PT[b,f] PT[c,d]-1/6 PT[a,f] PT[b,d] PT[c,e]-1/6 PT[a,d] PT[b,f] PT[c,e]-1/6 PT[a,e] PT[b,d] PT[c,f]-1/6 PT[a,d] PT[b,e] PT[c,f]-(PT[a,f] PT[b,c] PT[d,e])/(6 (-2+dimM))-(PT[a,c] PT[b,f] PT[d,e])/(6 (-2+dimM))+(PT[a,b] PT[c,f] PT[d,e])/(3 (-2+dimM))-(PT[a,e] PT[b,c] PT[d,f])/(6 (-2+dimM))-(PT[a,c] PT[b,e] PT[d,f])/(6 (-2+dimM))+(PT[a,b] PT[c,e] PT[d,f])/(3 (-2+dimM))+(PT[a,d] PT[b,c] PT[e,f])/(3 (-2+dimM))+(PT[a,c] PT[b,d] PT[e,f])/(3 (-2+dimM))-(2 PT[a,b] PT[c,d] PT[e,f])/(3 (-2+dimM))
P[2,m,{1,2}][c_,a_,b_,d_,e_,f_]:=(PT[a,f] PT[b,d] PT[c,e])/(2 Sqrt[3])+(PT[a,d] PT[b,f] PT[c,e])/(2 Sqrt[3])-(PT[a,e] PT[b,d] PT[c,f])/(2 Sqrt[3])-(PT[a,d] PT[b,e] PT[c,f])/(2 Sqrt[3])-(PT[a,f] PT[b,c] PT[d,e])/(2 Sqrt[3] (-2+dimM))-(PT[a,c] PT[b,f] PT[d,e])/(2 Sqrt[3] (-2+dimM))+(PT[a,b] PT[c,f] PT[d,e])/(Sqrt[3] (-2+dimM))+(PT[a,e] PT[b,c] PT[d,f])/(2 Sqrt[3] (-2+dimM))+(PT[a,c] PT[b,e] PT[d,f])/(2 Sqrt[3] (-2+dimM))-(PT[a,b] PT[c,e] PT[d,f])/(Sqrt[3] (-2+dimM))
P[2,m,{2,1}][c_,a_,b_,d_,e_,f_]:=-((PT[a,f] PT[b,d] PT[c,e])/(2 Sqrt[3]))+(PT[a,d] PT[b,f] PT[c,e])/(2 Sqrt[3])-(PT[a,e] PT[b,d] PT[c,f])/(2 Sqrt[3])+(PT[a,d] PT[b,e] PT[c,f])/(2 Sqrt[3])+(PT[a,f] PT[b,c] PT[d,e])/(2 Sqrt[3] (-2+dimM))-(PT[a,c] PT[b,f] PT[d,e])/(2 Sqrt[3] (-2+dimM))+(PT[a,e] PT[b,c] PT[d,f])/(2 Sqrt[3] (-2+dimM))-(PT[a,c] PT[b,e] PT[d,f])/(2 Sqrt[3] (-2+dimM))-(PT[a,d] PT[b,c] PT[e,f])/(Sqrt[3] (-2+dimM))+(PT[a,c] PT[b,d] PT[e,f])/(Sqrt[3] (-2+dimM))
P[2,m,{2,2}][c_,a_,b_,d_,e_,f_]:=-(1/3) PT[a,f] PT[b,e] PT[c,d]+1/3 PT[a,e] PT[b,f] PT[c,d]-1/6 PT[a,f] PT[b,d] PT[c,e]+1/6 PT[a,d] PT[b,f] PT[c,e]+1/6 PT[a,e] PT[b,d] PT[c,f]-1/6 PT[a,d] PT[b,e] PT[c,f]+(PT[a,f] PT[b,c] PT[d,e])/(2 (-2+dimM))-(PT[a,c] PT[b,f] PT[d,e])/(2 (-2+dimM))-(PT[a,e] PT[b,c] PT[d,f])/(2 (-2+dimM))+(PT[a,c] PT[b,e] PT[d,f])/(2 (-2+dimM))

P[2,p,{1,1}][c_,a_,b_,d_,e_,f_]:=1/6 PL[c,f] PT[a,e] PT[b,d]+1/6 PL[c,e] PT[a,f] PT[b,d]+1/6 PL[c,f] PT[a,d] PT[b,e]+1/6 PL[c,d] PT[a,f] PT[b,e]+1/6 PL[c,e] PT[a,d] PT[b,f]+1/6 PL[c,d] PT[a,e] PT[b,f]+1/6 PL[b,f] PT[a,e] PT[c,d]+1/6 PL[b,e] PT[a,f] PT[c,d]+1/6 PL[a,f] PT[b,e] PT[c,d]+1/6 PL[a,e] PT[b,f] PT[c,d]+1/6 PL[b,f] PT[a,d] PT[c,e]+1/6 PL[b,d] PT[a,f] PT[c,e]+1/6 PL[a,f] PT[b,d] PT[c,e]+1/6 PL[a,d] PT[b,f] PT[c,e]+1/6 PL[b,e] PT[a,d] PT[c,f]+1/6 PL[b,d] PT[a,e] PT[c,f]+1/6 PL[a,e] PT[b,d] PT[c,f]+1/6 PL[a,d] PT[b,e] PT[c,f]-(PL[c,f] PT[a,b] PT[d,e])/(3 (-1+dimM))-(PL[b,f] PT[a,c] PT[d,e])/(3 (-1+dimM))-(PL[a,f] PT[b,c] PT[d,e])/(3 (-1+dimM))-(PL[c,e] PT[a,b] PT[d,f])/(3 (-1+dimM))-(PL[b,e] PT[a,c] PT[d,f])/(3 (-1+dimM))-(PL[a,e] PT[b,c] PT[d,f])/(3 (-1+dimM))-(PL[c,d] PT[a,b] PT[e,f])/(3 (-1+dimM))-(PL[b,d] PT[a,c] PT[e,f])/(3 (-1+dimM))-(PL[a,d] PT[b,c] PT[e,f])/(3 (-1+dimM))
P[2,p,{1,2}][c_,a_,b_,d_,e_,f_]:=-((PL[c,f] PT[a,e] PT[b,d])/(6 Sqrt[2]))-(PL[c,e] PT[a,f] PT[b,d])/(6 Sqrt[2])-(PL[c,f] PT[a,d] PT[b,e])/(6 Sqrt[2])+(PL[c,d] PT[a,f] PT[b,e])/(3 Sqrt[2])-(PL[c,e] PT[a,d] PT[b,f])/(6 Sqrt[2])+(PL[c,d] PT[a,e] PT[b,f])/(3 Sqrt[2])-(PL[b,f] PT[a,e] PT[c,d])/(6 Sqrt[2])-(PL[b,e] PT[a,f] PT[c,d])/(6 Sqrt[2])-(PL[a,f] PT[b,e] PT[c,d])/(6 Sqrt[2])-(PL[a,e] PT[b,f] PT[c,d])/(6 Sqrt[2])-(PL[b,f] PT[a,d] PT[c,e])/(6 Sqrt[2])+(PL[b,d] PT[a,f] PT[c,e])/(3 Sqrt[2])-(PL[a,f] PT[b,d] PT[c,e])/(6 Sqrt[2])+(PL[a,d] PT[b,f] PT[c,e])/(3 Sqrt[2])-(PL[b,e] PT[a,d] PT[c,f])/(6 Sqrt[2])+(PL[b,d] PT[a,e] PT[c,f])/(3 Sqrt[2])-(PL[a,e] PT[b,d] PT[c,f])/(6 Sqrt[2])+(PL[a,d] PT[b,e] PT[c,f])/(3 Sqrt[2])+(PL[c,f] PT[a,b] PT[d,e])/(3 Sqrt[2] (-1+dimM))+(PL[b,f] PT[a,c] PT[d,e])/(3 Sqrt[2] (-1+dimM))+(PL[a,f] PT[b,c] PT[d,e])/(3 Sqrt[2] (-1+dimM))+(PL[c,e] PT[a,b] PT[d,f])/(3 Sqrt[2] (-1+dimM))+(PL[b,e] PT[a,c] PT[d,f])/(3 Sqrt[2] (-1+dimM))+(PL[a,e] PT[b,c] PT[d,f])/(3 Sqrt[2] (-1+dimM))-(Sqrt[2] PL[c,d] PT[a,b] PT[e,f])/(3 (-1+dimM))-(Sqrt[2] PL[b,d] PT[a,c] PT[e,f])/(3 (-1+dimM))-(Sqrt[2] PL[a,d] PT[b,c] PT[e,f])/(3 (-1+dimM))
P[2,p,{1,3}][c_,a_,b_,d_,e_,f_]:=-((PL[c,f] PT[a,e] PT[b,d])/(2 Sqrt[6]))+(PL[c,e] PT[a,f] PT[b,d])/(2 Sqrt[6])-(PL[c,f] PT[a,d] PT[b,e])/(2 Sqrt[6])+(PL[c,e] PT[a,d] PT[b,f])/(2 Sqrt[6])-(PL[b,f] PT[a,e] PT[c,d])/(2 Sqrt[6])+(PL[b,e] PT[a,f] PT[c,d])/(2 Sqrt[6])-(PL[a,f] PT[b,e] PT[c,d])/(2 Sqrt[6])+(PL[a,e] PT[b,f] PT[c,d])/(2 Sqrt[6])-(PL[b,f] PT[a,d] PT[c,e])/(2 Sqrt[6])-(PL[a,f] PT[b,d] PT[c,e])/(2 Sqrt[6])+(PL[b,e] PT[a,d] PT[c,f])/(2 Sqrt[6])+(PL[a,e] PT[b,d] PT[c,f])/(2 Sqrt[6])+(PL[c,f] PT[a,b] PT[d,e])/(Sqrt[6] (-1+dimM))+(PL[b,f] PT[a,c] PT[d,e])/(Sqrt[6] (-1+dimM))+(PL[a,f] PT[b,c] PT[d,e])/(Sqrt[6] (-1+dimM))-(PL[c,e] PT[a,b] PT[d,f])/(Sqrt[6] (-1+dimM))-(PL[b,e] PT[a,c] PT[d,f])/(Sqrt[6] (-1+dimM))-(PL[a,e] PT[b,c] PT[d,f])/(Sqrt[6] (-1+dimM))
P[2,p,{1,4}][c_,a_,b_,e_,f_]:=(PT[b,f] PT[c,e] q[a])/(2 Sqrt[3] qm)+(PT[b,e] PT[c,f] q[a])/(2 Sqrt[3] qm)-(PT[b,c] PT[e,f] q[a])/(Sqrt[3] (-1+dimM) qm)+(PT[a,f] PT[c,e] q[b])/(2 Sqrt[3] qm)+(PT[a,e] PT[c,f] q[b])/(2 Sqrt[3] qm)-(PT[a,c] PT[e,f] q[b])/(Sqrt[3] (-1+dimM) qm)+(PT[a,f] PT[b,e] q[c])/(2 Sqrt[3] qm)+(PT[a,e] PT[b,f] q[c])/(2 Sqrt[3] qm)-(PT[a,b] PT[e,f] q[c])/(Sqrt[3] (-1+dimM) qm)
P[2,p,{2,1}][c_,a_,b_,d_,e_,f_]:=(PL[c,f] PT[a,e] PT[b,d])/(3 Sqrt[2])+(PL[c,e] PT[a,f] PT[b,d])/(3 Sqrt[2])+(PL[c,f] PT[a,d] PT[b,e])/(3 Sqrt[2])+(PL[c,d] PT[a,f] PT[b,e])/(3 Sqrt[2])+(PL[c,e] PT[a,d] PT[b,f])/(3 Sqrt[2])+(PL[c,d] PT[a,e] PT[b,f])/(3 Sqrt[2])-(PL[b,f] PT[a,e] PT[c,d])/(6 Sqrt[2])-(PL[b,e] PT[a,f] PT[c,d])/(6 Sqrt[2])-(PL[a,f] PT[b,e] PT[c,d])/(6 Sqrt[2])-(PL[a,e] PT[b,f] PT[c,d])/(6 Sqrt[2])-(PL[b,f] PT[a,d] PT[c,e])/(6 Sqrt[2])-(PL[b,d] PT[a,f] PT[c,e])/(6 Sqrt[2])-(PL[a,f] PT[b,d] PT[c,e])/(6 Sqrt[2])-(PL[a,d] PT[b,f] PT[c,e])/(6 Sqrt[2])-(PL[b,e] PT[a,d] PT[c,f])/(6 Sqrt[2])-(PL[b,d] PT[a,e] PT[c,f])/(6 Sqrt[2])-(PL[a,e] PT[b,d] PT[c,f])/(6 Sqrt[2])-(PL[a,d] PT[b,e] PT[c,f])/(6 Sqrt[2])-(Sqrt[2] PL[c,f] PT[a,b] PT[d,e])/(3 (-1+dimM))+(PL[b,f] PT[a,c] PT[d,e])/(3 Sqrt[2] (-1+dimM))+(PL[a,f] PT[b,c] PT[d,e])/(3 Sqrt[2] (-1+dimM))-(Sqrt[2] PL[c,e] PT[a,b] PT[d,f])/(3 (-1+dimM))+(PL[b,e] PT[a,c] PT[d,f])/(3 Sqrt[2] (-1+dimM))+(PL[a,e] PT[b,c] PT[d,f])/(3 Sqrt[2] (-1+dimM))-(Sqrt[2] PL[c,d] PT[a,b] PT[e,f])/(3 (-1+dimM))+(PL[b,d] PT[a,c] PT[e,f])/(3 Sqrt[2] (-1+dimM))+(PL[a,d] PT[b,c] PT[e,f])/(3 Sqrt[2] (-1+dimM))
P[2,p,{2,2}][c_,a_,b_,d_,e_,f_]:=-(1/6) PL[c,f] PT[a,e] PT[b,d]-1/6 PL[c,e] PT[a,f] PT[b,d]-1/6 PL[c,f] PT[a,d] PT[b,e]+1/3 PL[c,d] PT[a,f] PT[b,e]-1/6 PL[c,e] PT[a,d] PT[b,f]+1/3 PL[c,d] PT[a,e] PT[b,f]+1/12 PL[b,f] PT[a,e] PT[c,d]+1/12 PL[b,e] PT[a,f] PT[c,d]+1/12 PL[a,f] PT[b,e] PT[c,d]+1/12 PL[a,e] PT[b,f] PT[c,d]+1/12 PL[b,f] PT[a,d] PT[c,e]-1/6 PL[b,d] PT[a,f] PT[c,e]+1/12 PL[a,f] PT[b,d] PT[c,e]-1/6 PL[a,d] PT[b,f] PT[c,e]+1/12 PL[b,e] PT[a,d] PT[c,f]-1/6 PL[b,d] PT[a,e] PT[c,f]+1/12 PL[a,e] PT[b,d] PT[c,f]-1/6 PL[a,d] PT[b,e] PT[c,f]+(PL[c,f] PT[a,b] PT[d,e])/(3 (-1+dimM))-(PL[b,f] PT[a,c] PT[d,e])/(6 (-1+dimM))-(PL[a,f] PT[b,c] PT[d,e])/(6 (-1+dimM))+(PL[c,e] PT[a,b] PT[d,f])/(3 (-1+dimM))-(PL[b,e] PT[a,c] PT[d,f])/(6 (-1+dimM))-(PL[a,e] PT[b,c] PT[d,f])/(6 (-1+dimM))-(2 PL[c,d] PT[a,b] PT[e,f])/(3 (-1+dimM))+(PL[b,d] PT[a,c] PT[e,f])/(3 (-1+dimM))+(PL[a,d] PT[b,c] PT[e,f])/(3 (-1+dimM))
P[2,p,{2,3}][c_,a_,b_,d_,e_,f_]:=-((PL[c,f] PT[a,e] PT[b,d])/(2 Sqrt[3]))+(PL[c,e] PT[a,f] PT[b,d])/(2 Sqrt[3])-(PL[c,f] PT[a,d] PT[b,e])/(2 Sqrt[3])+(PL[c,e] PT[a,d] PT[b,f])/(2 Sqrt[3])+(PL[b,f] PT[a,e] PT[c,d])/(4 Sqrt[3])-(PL[b,e] PT[a,f] PT[c,d])/(4 Sqrt[3])+(PL[a,f] PT[b,e] PT[c,d])/(4 Sqrt[3])-(PL[a,e] PT[b,f] PT[c,d])/(4 Sqrt[3])+(PL[b,f] PT[a,d] PT[c,e])/(4 Sqrt[3])+(PL[a,f] PT[b,d] PT[c,e])/(4 Sqrt[3])-(PL[b,e] PT[a,d] PT[c,f])/(4 Sqrt[3])-(PL[a,e] PT[b,d] PT[c,f])/(4 Sqrt[3])+(PL[c,f] PT[a,b] PT[d,e])/(Sqrt[3] (-1+dimM))-(PL[b,f] PT[a,c] PT[d,e])/(2 Sqrt[3] (-1+dimM))-(PL[a,f] PT[b,c] PT[d,e])/(2 Sqrt[3] (-1+dimM))-(PL[c,e] PT[a,b] PT[d,f])/(Sqrt[3] (-1+dimM))+(PL[b,e] PT[a,c] PT[d,f])/(2 Sqrt[3] (-1+dimM))+(PL[a,e] PT[b,c] PT[d,f])/(2 Sqrt[3] (-1+dimM))
P[2,p,{2,4}][c_,a_,b_,e_,f_]:=-((PT[b,f] PT[c,e] q[a])/(2 Sqrt[6] qm))-(PT[b,e] PT[c,f] q[a])/(2 Sqrt[6] qm)+(PT[b,c] PT[e,f] q[a])/(Sqrt[6] (-1+dimM) qm)-(PT[a,f] PT[c,e] q[b])/(2 Sqrt[6] qm)-(PT[a,e] PT[c,f] q[b])/(2 Sqrt[6] qm)+(PT[a,c] PT[e,f] q[b])/(Sqrt[6] (-1+dimM) qm)+(PT[a,f] PT[b,e] q[c])/(Sqrt[6] qm)+(PT[a,e] PT[b,f] q[c])/(Sqrt[6] qm)-(Sqrt[2/3] PT[a,b] PT[e,f] q[c])/((-1+dimM) qm)
P[2,p,{3,1}][c_,a_,b_,d_,e_,f_]:=-((PL[b,f] PT[a,e] PT[c,d])/(2 Sqrt[6]))-(PL[b,e] PT[a,f] PT[c,d])/(2 Sqrt[6])+(PL[a,f] PT[b,e] PT[c,d])/(2 Sqrt[6])+(PL[a,e] PT[b,f] PT[c,d])/(2 Sqrt[6])-(PL[b,f] PT[a,d] PT[c,e])/(2 Sqrt[6])-(PL[b,d] PT[a,f] PT[c,e])/(2 Sqrt[6])+(PL[a,f] PT[b,d] PT[c,e])/(2 Sqrt[6])+(PL[a,d] PT[b,f] PT[c,e])/(2 Sqrt[6])-(PL[b,e] PT[a,d] PT[c,f])/(2 Sqrt[6])-(PL[b,d] PT[a,e] PT[c,f])/(2 Sqrt[6])+(PL[a,e] PT[b,d] PT[c,f])/(2 Sqrt[6])+(PL[a,d] PT[b,e] PT[c,f])/(2 Sqrt[6])+(PL[b,f] PT[a,c] PT[d,e])/(Sqrt[6] (-1+dimM))-(PL[a,f] PT[b,c] PT[d,e])/(Sqrt[6] (-1+dimM))+(PL[b,e] PT[a,c] PT[d,f])/(Sqrt[6] (-1+dimM))-(PL[a,e] PT[b,c] PT[d,f])/(Sqrt[6] (-1+dimM))+(PL[b,d] PT[a,c] PT[e,f])/(Sqrt[6] (-1+dimM))-(PL[a,d] PT[b,c] PT[e,f])/(Sqrt[6] (-1+dimM))
P[2,p,{3,2}][c_,a_,b_,d_,e_,f_]:=(PL[b,f] PT[a,e] PT[c,d])/(4 Sqrt[3])+(PL[b,e] PT[a,f] PT[c,d])/(4 Sqrt[3])-(PL[a,f] PT[b,e] PT[c,d])/(4 Sqrt[3])-(PL[a,e] PT[b,f] PT[c,d])/(4 Sqrt[3])+(PL[b,f] PT[a,d] PT[c,e])/(4 Sqrt[3])-(PL[b,d] PT[a,f] PT[c,e])/(2 Sqrt[3])-(PL[a,f] PT[b,d] PT[c,e])/(4 Sqrt[3])+(PL[a,d] PT[b,f] PT[c,e])/(2 Sqrt[3])+(PL[b,e] PT[a,d] PT[c,f])/(4 Sqrt[3])-(PL[b,d] PT[a,e] PT[c,f])/(2 Sqrt[3])-(PL[a,e] PT[b,d] PT[c,f])/(4 Sqrt[3])+(PL[a,d] PT[b,e] PT[c,f])/(2 Sqrt[3])-(PL[b,f] PT[a,c] PT[d,e])/(2 Sqrt[3] (-1+dimM))+(PL[a,f] PT[b,c] PT[d,e])/(2 Sqrt[3] (-1+dimM))-(PL[b,e] PT[a,c] PT[d,f])/(2 Sqrt[3] (-1+dimM))+(PL[a,e] PT[b,c] PT[d,f])/(2 Sqrt[3] (-1+dimM))+(PL[b,d] PT[a,c] PT[e,f])/(Sqrt[3] (-1+dimM))-(PL[a,d] PT[b,c] PT[e,f])/(Sqrt[3] (-1+dimM))
P[2,p,{3,3}][c_,a_,b_,d_,e_,f_]:=1/4 PL[b,f] PT[a,e] PT[c,d]-1/4 PL[b,e] PT[a,f] PT[c,d]-1/4 PL[a,f] PT[b,e] PT[c,d]+1/4 PL[a,e] PT[b,f] PT[c,d]+1/4 PL[b,f] PT[a,d] PT[c,e]-1/4 PL[a,f] PT[b,d] PT[c,e]-1/4 PL[b,e] PT[a,d] PT[c,f]+1/4 PL[a,e] PT[b,d] PT[c,f]-(PL[b,f] PT[a,c] PT[d,e])/(2 (-1+dimM))+(PL[a,f] PT[b,c] PT[d,e])/(2 (-1+dimM))+(PL[b,e] PT[a,c] PT[d,f])/(2 (-1+dimM))-(PL[a,e] PT[b,c] PT[d,f])/(2 (-1+dimM))
P[2,p,{3,4}][c_,a_,b_,e_,f_]:=(PT[b,f] PT[c,e] q[a])/(2 Sqrt[2] qm)+(PT[b,e] PT[c,f] q[a])/(2 Sqrt[2] qm)-(PT[b,c] PT[e,f] q[a])/(3 Sqrt[2] (-1+dimM) qm)-(Sqrt[2] PT[b,c] PT[e,f] q[a])/(3 (-1+dimM) qm)-(PT[a,f] PT[c,e] q[b])/(2 Sqrt[2] qm)-(PT[a,e] PT[c,f] q[b])/(2 Sqrt[2] qm)+(PT[a,c] PT[e,f] q[b])/(3 Sqrt[2] (-1+dimM) qm)+(Sqrt[2] PT[a,c] PT[e,f] q[b])/(3 (-1+dimM) qm)
P[2,p,{4,1}][a_,b_,d_,e_,f_]:=(PT[a,f] PT[b,e] q[d])/(2 Sqrt[3] qm)+(PT[a,e] PT[b,f] q[d])/(2 Sqrt[3] qm)-(PT[a,b] PT[e,f] q[d])/(Sqrt[3] (-1+dimM) qm)+(PT[a,f] PT[b,d] q[e])/(2 Sqrt[3] qm)+(PT[a,d] PT[b,f] q[e])/(2 Sqrt[3] qm)-(PT[a,b] PT[d,f] q[e])/(Sqrt[3] (-1+dimM) qm)+(PT[a,e] PT[b,d] q[f])/(2 Sqrt[3] qm)+(PT[a,d] PT[b,e] q[f])/(2 Sqrt[3] qm)-(PT[a,b] PT[d,e] q[f])/(Sqrt[3] (-1+dimM) qm)
P[2,p,{4,2}][a_,b_,d_,e_,f_]:=(PT[a,f] PT[b,e] q[d])/(Sqrt[6] qm)+(PT[a,e] PT[b,f] q[d])/(Sqrt[6] qm)-(Sqrt[2/3] PT[a,b] PT[e,f] q[d])/((-1+dimM) qm)-(PT[a,f] PT[b,d] q[e])/(2 Sqrt[6] qm)-(PT[a,d] PT[b,f] q[e])/(2 Sqrt[6] qm)+(PT[a,b] PT[d,f] q[e])/(Sqrt[6] (-1+dimM) qm)-(PT[a,e] PT[b,d] q[f])/(2 Sqrt[6] qm)-(PT[a,d] PT[b,e] q[f])/(2 Sqrt[6] qm)+(PT[a,b] PT[d,e] q[f])/(Sqrt[6] (-1+dimM) qm)
P[2,p,{4,3}][a_,b_,d_,e_,f_]:=(PT[a,f] PT[b,d] q[e])/(2 Sqrt[2] qm)+(PT[a,d] PT[b,f] q[e])/(2 Sqrt[2] qm)-(PT[a,b] PT[d,f] q[e])/(3 Sqrt[2] (-1+dimM) qm)-(Sqrt[2] PT[a,b] PT[d,f] q[e])/(3 (-1+dimM) qm)-(PT[a,e] PT[b,d] q[f])/(2 Sqrt[2] qm)-(PT[a,d] PT[b,e] q[f])/(2 Sqrt[2] qm)+(PT[a,b] PT[d,e] q[f])/(3 Sqrt[2] (-1+dimM) qm)+(Sqrt[2] PT[a,b] PT[d,e] q[f])/(3 (-1+dimM) qm)
P[2,p,{4,4}][a_,b_,e_,f_]:=1/2 (PT[a,f] PT[b,e]+PT[a,e] PT[b,f])-1/(-1+dimM) PT[a,b] PT[e,f]

P[3,m,{1,1}][c_,a_,b_,d_,e_,f_]:=1/6 PT[a,f] PT[b,e] PT[c,d]+1/6 PT[a,e] PT[b,f] PT[c,d]+1/6 PT[a,f] PT[b,d] PT[c,e]+1/6 PT[a,d] PT[b,f] PT[c,e]+1/6 PT[a,e] PT[b,d] PT[c,f]+1/6 PT[a,d] PT[b,e] PT[c,f]-(PT[a,f] PT[b,c] PT[d,e])/(3 (1+dimM))-(PT[a,c] PT[b,f] PT[d,e])/(3 (1+dimM))-(PT[a,b] PT[c,f] PT[d,e])/(3 (1+dimM))-(PT[a,e] PT[b,c] PT[d,f])/(3 (1+dimM))-(PT[a,c] PT[b,e] PT[d,f])/(3 (1+dimM))-(PT[a,b] PT[c,e] PT[d,f])/(3 (1+dimM))-(PT[a,d] PT[b,c] PT[e,f])/(3 (1+dimM))-(PT[a,c] PT[b,d] PT[e,f])/(3 (1+dimM))-(PT[a,b] PT[c,d] PT[e,f])/(3 (1+dimM))

(*r3r1*)
P[0,p,{1,7}][a_,b_,c_,m_]:=(E^((-I)*XX1+I*XX2)*(qq*g[b,c]*q[a]+qq*g[a,c]*q[b]+(qq*g[a,b]-3*q[a]*q[b])*q[c])*q[m])/(3*qq^2);
P[0,p,{2,7}][a_,b_,c_,m_]:=(E^((-I)*XX1+I*XX2)*(2*g[b,c]*q[a]-g[a,c]*q[b]-g[a,b]*q[c])*q[m])/(3*Sqrt[2]*qq);
P[0,p,{3,7}][a_,b_,c_,m_]:=(E^((-I)*XX1+I*XX2)*(g[a,c]*q[b]-g[a,b]*q[c])*q[m])/(Sqrt[6]*qq);
P[0,p,{4,7}][a_,b_,c_,m_]:=(E^((-I)*XX1+I*XX2)*q[a]*q[b]*q[c]*q[m])/qq^2;
(*r2symr1*)
P[0,p,{5,7}][a_,b_,m_]:=(E^((-I)*XX1+I*XX2)*(qq*g[a,b]-q[a]*q[b])*q[m])/(Sqrt[3]*qq^(3/2));
P[0,p,{6,7}][a_,b_,m_]:=(E^((-I)*XX1+I*XX2)*q[a]*q[b]*q[m])/qq^(3/2);
(*r1r1*)
P[0,p,{7,7}][a_,b_]:=PL[a,b];
(*r1r3*)
P[0,p,{7,1}][a_,b_,c_,m_]:=(E^(I*XX1-I*XX2)*q[a]*(qq*g[c,m]*q[b]+qq*g[b,m]*q[c]+(qq*g[b,c]-3*q[b]*q[c])*q[m]))/(3*qq^2);
P[0,p,{7,2}][a_,b_,c_,m_]:=(E^(I*XX1-I*XX2)*q[a]*(2*g[c,m]*q[b]-g[b,m]*q[c]-g[b,c]*q[m]))/(3*Sqrt[2]*qq);
P[0,p,{7,3}][a_,b_,c_,m_]:=(E^(I*XX1-I*XX2)*q[a]*(g[b,m]*q[c]-g[b,c]*q[m]))/(Sqrt[6]*qq);
P[0,p,{7,4}][a_,b_,c_,m_]:=(E^(I*XX1-I*XX2)*q[a]*q[b]*q[c]*q[m])/qq^2;
(*r1r2sym*)
P[0,p,{7,5}][a_,b_,m_]:=(E^(I*XX1-I*XX2)*q[a]*(qq*g[b,m]-q[b]*q[m]))/(Sqrt[3]*qq^(3/2));

P[0,p,{7,6}][a_,b_,m_]:=(E^(I*XX1-I*XX2)*q[a]*q[b]*q[m])/qq^(3/2);
(*r3r1*)
(*{x,y,z,t,r,p,q,w,a,b,n,m,l,k,j,f,o,z,c}*)
P[1,m,{1,8}][a_,b_,c_,m_]:=(qq^2*g[a,b]*g[c,m]-qq*g[c,m]*q[a]*q[b]-qq*g[b,m]*q[a]*q[c]+qq*g[a,m]*(qq*g[b,c]-q[b]*q[c])-qq*g[b,c]*q[a]*q[m]-qq*g[a,b]*q[c]*q[m]+3*q[a]*q[b]*q[c]*q[m]+qq*g[a,c]*(qq*g[b,m]-q[b]*q[m]))/(Sqrt[15]*E^(I*XX3)*qq^2);
P[1,m,{2,8}][a_,b_,c_,m_]:=-1/2*(-(qq*g[a,b]*g[c,m])+g[c,m]*q[a]*q[b]+g[b,m]*q[a]*q[c]+2*g[a,m]*(qq*g[b,c]-q[b]*q[c])-2*g[b,c]*q[a]*q[m]+g[a,b]*q[c]*q[m]+g[a,c]*(-(qq*g[b,m])+q[b]*q[m]))/(Sqrt[3]*E^(I*XX3)*qq);
P[1,m,{3,8}][a_,b_,c_,m_]:=(q[a]*(-(g[c,m]*q[b])+g[b,m]*q[c])+g[a,c]*(-(qq*g[b,m])+q[b]*q[m])+g[a,b]*(qq*g[c,m]-q[c]*q[m]))/(2*E^(I*XX3)*qq);
P[1,m,{4,8}][a_,b_,c_,m_]:=(qq*g[c,m]*q[a]*q[b]+q[c]*(qq*g[b,m]*q[a]+q[b]*(qq*g[a,m]-3*q[a]*q[m])))/(Sqrt[3]*E^(I*XX3)*qq^2);
P[1,m,{5,8}][a_,b_,c_,m_]:=(g[c,m]*q[a]*q[b]+g[b,m]*q[a]*q[c]-2*g[a,m]*q[b]*q[c])/(Sqrt[6]*E^(I*XX3)*qq);
P[1,m,{6,8}][a_,b_,c_,m_]:=(q[a]*(g[c,m]*q[b]-g[b,m]*q[c]))/(Sqrt[2]*E^(I*XX3)*qq);
P[1,m,{7,8}][a_,b_,m_]:=(qq*g[b,m]*q[a]+q[b]*(qq*g[a,m]-2*q[a]*q[m]))/(Sqrt[2]*E^(I*XX3)*qq^(3/2));
(*r1r1*)
P[1,m,{8,8}][a_,b_]:=PT[a,b];
(*{x,y,z,t,r,p,q,w,a,b,n,m,l,k,j,f,o,z,c}*)
(*r1r3*)
P[1,m,{8,1}][a_,b_,c_,m_]:=(E^(I*XX3)*(qq^2*g[a,b]*g[c,m]-qq*g[c,m]*q[a]*q[b]-qq*g[b,m]*q[a]*q[c]+qq*g[a,m]*(qq*g[b,c]-q[b]*q[c])-qq*g[b,c]*q[a]*q[m]-qq*g[a,b]*q[c]*q[m]+3*q[a]*q[b]*q[c]*q[m]+qq*g[a,c]*(qq*g[b,m]-q[b]*q[m])))/(Sqrt[15]*qq^2);
P[1,m,{8,2}][a_,b_,c_,m_]:=-1/2*(E^(I*XX3)*(2*qq*g[a,b]*g[c,m]-2*g[c,m]*q[a]*q[b]+g[b,m]*q[a]*q[c]+g[a,m]*(-(qq*g[b,c])+q[b]*q[c])+g[b,c]*q[a]*q[m]-2*g[a,b]*q[c]*q[m]+g[a,c]*(-(qq*g[b,m])+q[b]*q[m])))/(Sqrt[3]*qq);
P[1,m,{8,3}][a_,b_,c_,m_]:=(E^(I*XX3)*(g[a,m]*(qq*g[b,c]-q[b]*q[c])+q[a]*(g[b,m]*q[c]-g[b,c]*q[m])+g[a,c]*(-(qq*g[b,m])+q[b]*q[m])))/(2*qq);
P[1,m,{8,4}][a_,b_,c_,m_]:=(E^(I*XX3)*(qq*g[a,m]*q[b]*q[c]+(qq*g[a,c]*q[b]+(qq*g[a,b]-3*q[a]*q[b])*q[c])*q[m]))/(Sqrt[3]*qq^2);
P[1,m,{8,5}][a_,b_,c_,m_]:=(E^(I*XX3)*(g[a,m]*q[b]*q[c]+(g[a,c]*q[b]-2*g[a,b]*q[c])*q[m]))/(Sqrt[6]*qq);
P[1,m,{8,6}][a_,b_,c_,m_]:=(E^(I*XX3)*q[b]*(g[a,m]*q[c]-g[a,c]*q[m]))/(Sqrt[2]*qq);
(*r1r2*)
P[1,m,{8,7}][a_,b_,m_]:=(E^(I*XX3)*(qq*g[a,m]*q[b]+(qq*g[a,b]-2*q[a]*q[b])*q[m]))/(Sqrt[2]*qq^(3/2))
(*r3r0*)
P[0,p,{1,8}][a_,b_,c_]:=(E^(I*XX2)*(qq*g[b,c]*q[a]+qq*g[a,c]*q[b]+(qq*g[a,b]-3*q[a]*q[b])*q[c]))/(3*qq^(3/2));
P[0,p,{2,8}][a_,b_,c_]:=(E^(I*XX2)*(2*g[b,c]*q[a]-g[a,c]*q[b]-g[a,b]*q[c]))/(3*Sqrt[2]*Sqrt[qq]);
P[0,p,{3,8}][a_,b_,c_]:=(E^(I*XX2)*(g[a,c]*q[b]-g[a,b]*q[c]))/(Sqrt[6]*Sqrt[qq]);
P[0,p,{4,8}][a_,b_,c_]:=(E^(I*XX2)*q[a]*q[b]*q[c])/qq^(3/2);
(*r2symr0*)
P[0,p,{5,8}][a_,b_]:=(E^(I*XX2)*(qq*g[a,b]-q[a]*q[b]))/(Sqrt[3]*qq);
P[0,p,{6,8}][a_,b_]:=(E^(I*XX2)*q[a]*q[b])/qq
(*r1r0*)
P[0,p,{7,8}][a_]:=(E^(I*XX1)*q[a])/Sqrt[qq]
(*r0r0*)
P[0,p,{8,8}][]:=1
(*r3r0*)
P[0,p,{8,1}][a_,b_,c_]:=(qq*g[b,c]*q[a]+qq*g[a,c]*q[b]+(qq*g[a,b]-3*q[a]*q[b])*q[c])/(3*E^(I*XX2)*qq^(3/2));
P[0,p,{8,2}][a_,b_,c_]:=(2*g[b,c]*q[a]-g[a,c]*q[b]-g[a,b]*q[c])/(3*Sqrt[2]*E^(I*XX2)*Sqrt[qq]);
P[0,p,{8,3}][a_,b_,c_]:=(g[a,c]*q[b]-g[a,b]*q[c])/(Sqrt[6]*E^(I*XX2)*Sqrt[qq]);
P[0,p,{8,4}][a_,b_,c_]:=(q[a]*q[b]*q[c])/(E^(I*XX2)*qq^(3/2));
(*r2symr0*)
P[0,p,{8,5}][a_,b_]:=(qq*g[a,b]-q[a]*q[b])/(Sqrt[3]*E^(I*XX2)*qq);
P[0,p,{8,6}][a_,b_]:=(q[a]*q[b])/(E^(I*XX2)*qq);
(*r1r0*)
P[0,p,{8,7}][a_]:=q[a]/(E^(I*XX1)*Sqrt[qq]);

(* === End of Projector Library === *)
