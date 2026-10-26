(* ::Package:: *)

BeginPackage["Kummitus`"];

Print @ Row[{
  Style["Package ", GrayLevel[0.3]],
  Style["Kummitus", Darker[Blue], Bold],
  Style[" -- dof from poles!
CopyRight \[Copyright] 2026, Carlo Marzo, under the General Public License.", GrayLevel[0.3]]
}];

Print[" "];
Print @ Row[{
  Style["Powered by ", GrayLevel[0.3]],
  Style["xAct", Black, Bold],
  Style[" (xTensor, xPert, xCoba, xTras)", GrayLevel[0.3]]
}];
Print[" "];
 
 Block[{Print = (Null &)},

Needs["xAct`xTensor`"];
Print[" "];
Print[" "];
Needs["xAct`xPert`"];
Print[" "];
Print[" "];
Needs["xAct`xCoba`"];
Print[" "];
Print[" "];
Needs["xAct`xTras`"];
Print[" "];
Print[" "];

];

(*xAct`$CVVerbose = True;*)

Off[xAct`xTensor`VarD::nouse];
Off[xAct`xTensor`DefMetric::error];
Off[Solve::svars];


pkgDir = DirectoryName[$InputFileName];

(*Get[FileNameJoin[{pkgDir, "FieldConfiguration.m"}]];
*)
Get[FileNameJoin[{pkgDir, "SystemSolver.wl"}]];


(* Public symbols + usage *)
dd::usage = "Flat Derivative - Access to XAct's PD fiducial derivative.";
ConstantSymbols::usage = "Access XAct' DefConstantSymbol.";
Canonicalize::usage = "Access XAct' CollectTensors.";
DeclareField::usage = "DeclareField[X, type, opts] declares a new field X of the given type (\"rank3\", \"rank2sym\", \"rank2anti\", \"vector\", \"scalar\"), together with its perturbation \[Epsilon]Xp and sources XJ, XJc. Several fields of the same type may be declared. Options: \"Symmetry\" (rank-3 only, default rank3SymmetryType), \"Perturbation\", \"Source\", \"ConjugateSource\".";
FieldInstances::usage = "FieldInstances[] returns the registry of declared fields.";


SourceDeclare::usage = "It defines all the necessary rules to go from tensor to the minimal set of independent components.";
CompSatProp::usage = "Compute the saturated propagator near the appropriate poles (Massive and not).";
ComputeMassiveStates::usage = "Compute the residues around the non-zero poles.";
ComputeMasslessStates::usage = "Compute the residues around the zero poles.";
matricesOutput::usage = "Access the SPO matrices for external use. Like Survey.";
PrintK::usage = "Print with Style";


masses::usage = "Access the propagator's poles of the last successfull run.";
masslessData::usage = "Access the residues of the saturated propagator in the limit qq -> 0.";
massiveData::usage = "Access the masses and residues of the saturated propagator in the limit qq -> 0.";
LatestPropMassive::usage = "Component form of the Saturated Propagator in Light-Like frame.";
LatestPropMassless::usage = "Component form of the Saturated Propagator in Time-Like frame.";

masslessSolution::usage = "...";
masslessSolutionConj::usage = "...";
TensorToHelicity::usage = "TensorToHelicity[] returns rules expressing the source components of the active fields in terms of helicity amplitudes, reduced by the massless source constraints. Apply them to LatestPropMassless. Built-in fields use sv[h,i], sh[h,i], sa[h,i], ss[h,i]; other fields carry their name as a third argument, e.g. sv[h,i,V2]; Kalb-Ramond fields use sb[h,i,B].";


sv::usage = "Helicity amplitude.";
svc::usage = "Conjugate helicity amplitude.";
sh::usage = "Helicity amplitude.";
shc::usage = "Conjugate helicity amplitude.";
sa::usage = "Helicity amplitude.";
sac::usage = "Conjugate helicity amplitude.";
ss::usage = "Helicity amplitude (scalar).";
ssc::usage = "Conjugate helicity amplitude (scalar).";
sb::usage = "Helicity amplitude (antisymmetric rank-2).";
sbc::usage = "Conjugate helicity amplitude (antisymmetric rank-2).";



(*Manifold definitions:*)


(*SetAttributes[PrintK, HoldAll];
*)
(* there must be a smarter way of doing this... *)

PrintK[expr_] := Print[Style[expr, GrayLevel[0.3]]];

PrintK[expr1_,term1_] := Print[Style[expr1, GrayLevel[0.3]],term1];
PrintK[expr1_,term1_,expr2_] := Print[Style[expr1, GrayLevel[0.3]], term1, Style[expr2, GrayLevel[0.3]]];
PrintK[expr1_,term1_,expr2_,term2_] := Print[Style[expr1, GrayLevel[0.3]], term1, Style[expr2, GrayLevel[0.3]], term2];
PrintK[expr1_,term1_,expr2_,term2_, expr3_] := Print[Style[expr1, GrayLevel[0.3]],term1, Style[expr2, GrayLevel[0.3]], term2, Style[expr3, GrayLevel[0.3]]];
PrintK[expr1_,term1_,expr2_,term2_, expr3_,term3_] := Print[Style[expr1, GrayLevel[0.3]],term1, Style[expr2, GrayLevel[0.3]], term2, Style[expr3, GrayLevel[0.3]], term3];
PrintK[expr1_,term1_,expr2_,term2_, expr3_,term3_, expr4_] := Print[Style[expr1, GrayLevel[0.3]],term1, Style[expr2, GrayLevel[0.3]], term2, Style[expr3, GrayLevel[0.3]], term3, Style[expr4, GrayLevel[0.3]]];
PrintK[expr1_,term1_,expr2_,term2_, expr3_,term3_, expr4_,term4_] := Print[Style[expr1, GrayLevel[0.3]],term1, Style[expr2, GrayLevel[0.3]], term2, Style[expr3, GrayLevel[0.3]], term3, Style[expr4, GrayLevel[0.3]], term4];

(* ... or more arguments: behave exactly like Print *)
(*PrintK[x_, y_, rest___] := Print[x, y, rest];*)



Block[{Print = (Null &)},

(* 4D manifold - exclude g and q from indices since we use them for metric and momentum *)
DefManifold[M, 4, DeleteCases[IndexRange[a, t], g | q]];

(* Flat metric with partial derivative PD *)
DefMetric[-1, g[-a, -b], PD, FlatMetric -> True];

(* Momentum vector *)
DefTensor[q[-a], M];

(* Momentum squared: q^2 = q_a q^a *)
DefConstantSymbol[qq];
q /: q[a_] q[-a_] := qq;
q /: q[-a_] q[a_] := qq;

(* Momentum magnitude *)
DefConstantSymbol[qm];
qm /: qm^2 := qq;
qm /: qm^4 := qq^2;
qm /: qm^6 := qq^3;
qm /: qm^8 := qq^4;

(* Dimension symbol *)
DefConstantSymbol[dimM];
dimM = 4;

DefBasis[Lorentz, TangentM, {0, 1, 2, 3}];
MetricRule1 = MetricInBasis[g, -Lorentz, {1, -1, -1, -1}] // Flatten;
MetricRule2 = MetricInBasis[g, Lorentz, {1, -1, -1, -1}] // Flatten;
MomentumRule1 = ChangeComponents[q[-{a, Lorentz}], q[{a, Lorentz}]];

(*Energy and Momenta*)
DefConstantSymbol[\[Epsilon]];
DefConstantSymbol[\[Omega]];
DefConstantSymbol[\[Kappa]];

$Assumptions={\[Omega]>0,Element[\[Kappa],Reals]};


(* === Control variables for mixed operators === *)
DefConstantSymbol[#]&/@{XX1,XX2,XX3,XX4};


];



Block[{Print = (Null &)},

projectorsLibraryPath = 
  FileNameJoin[{DirectoryName[$InputFileName], "ProjectorsLibrary.m"}];
];

PrintK["============================================================================"]

(* Alternative: use absolute path *)
(* projectorsLibraryPath = "/path/to/ProjectorsLibrary.m"; *)

(* Load the library *)
If[FileExistsQ[projectorsLibraryPath],
    Get[projectorsLibraryPath];
    PrintK["Loaded projector library from: "];
    PrintK[projectorsLibraryPath];
,
    PrintK["WARNING: Projector library not found"];
    PrintK["Please ensure ProjectorsLibrary.m is in the same directory as this notebook."];
    PrintK["Or set projectorsLibraryPath to the correct path."];
];

XX1=XX2=(*XX3=*)XX4=0; (*These are control-variables used to check the mixed projectors operators. I used them in my original computation. Never appear in physical objects. Choosing 0.*)
(*XX3 = \[Pi] required by brute force calculation of mixed operators for antisymmetric field *)




If[!ValueQ[Global`DivideEtImpera],
  Global`DivideEtImpera = True;
];

If[!ValueQ[Global`rank3SymmetryType],
  Global`rank3SymmetryType = "TotallySymmetric";
];


PrintK["Defining fields. Rank-3 field A is declared: ", Global`rank3SymmetryType];
PrintK["You can change the value of the string variable 'rank3SymmetryType' before starting the package. You can choose among 'Generic', 'Sym23' , 'Antisym23', 'TotallySymmetric', 'TotallyAntisymmetric'  ."];
PrintK["Use H for rank-2 totally symmetric tensor, V for vector and S for scalar. More fields, including copies of the same type and the antisymmetric \"rank2anti\" (Kalb-Ramond), are added with DeclareField[name, type]."];
PrintK["The source of 'field' is 'field name'J. The complex conjugate source is 'field name'Jc."];
PrintK["============================================================================"]

(*BEGIN FIELD DEFINITION*)

Begin["`Private`"]

(* === Field-instance registry ===
   Every field in the theory is an instance of one of the supported types. Several instances
   of the same type (copies) may coexist: they share the projectors and differ only by name.
   name -> <|"type", "rank", "symmetry", "pert", "J", "Jc", "order"|> *)
$FieldInstances = <||>;

$FieldTypeRank = <|"rank3" -> 3, "rank2sym" -> 2, "rank2anti" -> 2, "vector" -> 1, "scalar" -> 0|>;

$Rank3Symmetries = {"Generic", "Sym23", "Antisym23", "TotallySymmetric", "TotallyAntisymmetric"};

(* xAct symmetry specification on the given index slots *)
symmetrySpec["rank3", "Sym23", {_, i2_, i3_}] := Symmetric[{i2, i3}];
symmetrySpec["rank3", "Antisym23", {_, i2_, i3_}] := Antisymmetric[{i2, i3}];
symmetrySpec["rank3", "TotallySymmetric", inds_] := Symmetric[inds];
symmetrySpec["rank3", "TotallyAntisymmetric", inds_] := Antisymmetric[inds];
symmetrySpec["rank2sym", _, inds_] := Symmetric[inds];
symmetrySpec["rank2anti", _, inds_] := Antisymmetric[inds];
symmetrySpec[_, _, _] := Sequence[];

Options[DeclareField] = {
    "Symmetry" :> Global`rank3SymmetryType,
    "Perturbation" -> Automatic,
    "Source" -> Automatic,
    "ConjugateSource" -> Automatic
};

DeclareField::type = "Unknown field type `1`. Supported types: `2`.";
DeclareField::exists = "Field `1` is already declared.";
DeclareField::sym = "Unknown rank-3 symmetry `1`; using \"Generic\". Supported: `2`.";

DeclareField[name_Symbol, type_String, OptionsPattern[]] := Module[
    {rank, sym, pert, src, srcc, inds},

    If[!KeyExistsQ[$FieldTypeRank, type],
        Message[DeclareField::type, type, Keys[$FieldTypeRank]];
        Return[$Failed]
    ];
    If[KeyExistsQ[$FieldInstances, name],
        Message[DeclareField::exists, name];
        Return[$Failed]
    ];

    rank = $FieldTypeRank[type];
    sym = If[type === "rank3", OptionValue["Symmetry"], None];
    If[type === "rank3" && !MemberQ[$Rank3Symmetries, sym],
        Message[DeclareField::sym, sym, $Rank3Symmetries];
        sym = "Generic";
    ];

    (* Default names: field X -> perturbation \[Epsilon]Xp, sources XJ and XJc *)
    pert = Replace[OptionValue["Perturbation"], Automatic :> Symbol["\[Epsilon]" <> SymbolName[name] <> "p"]];
    src  = Replace[OptionValue["Source"], Automatic :> Symbol[SymbolName[name] <> "J"]];
    srcc = Replace[OptionValue["ConjugateSource"], Automatic :> Symbol[SymbolName[name] <> "Jc"]];

    inds = Take[{-a, -b, -c}, rank];

    Block[{Print = (Null &)},
        DefTensor[name @@ inds, M, symmetrySpec[type, sym, inds]];
        DefTensorPerturbation[pert @@ Prepend[inds, LI[1]], name @@ inds, M, symmetrySpec[type, sym, inds]];
        (* Scalar sources are constants; the others carry upper indices *)
        If[rank === 0,
            DefConstantSymbol[src];
            DefConstantSymbol[srcc];,
            DefTensor[src @@ (-inds), M, symmetrySpec[type, sym, -inds]];
            DefTensor[srcc @@ (-inds), M, symmetrySpec[type, sym, -inds]];
        ];
    ];

    $FieldInstances[name] = <|
        "type" -> type, "rank" -> rank, "symmetry" -> sym,
        "pert" -> pert, "J" -> src, "Jc" -> srcc,
        "order" -> Length[$FieldInstances] + 1
    |>;
    name
];

FieldInstances[] := $FieldInstances;

End[] (* `Private` *)

(* Built-in fields. The rank-3 symmetry follows rank3SymmetryType. *)
DeclareField[A, "rank3"];
DeclareField[H, "rank2sym", "Perturbation" -> \[Epsilon]h];
DeclareField[V, "vector"];
DeclareField[S, "scalar"];

(*END FIELD DEFINITION*)

Begin["`Private`"]

ConstantSymbols[input_]:= Block[{Print = (Null &)},DefConstantSymbol[input];];

Canonicalize[input_]:= CollectTensors[input]


dd[input_][input2_]:=PD[input][input2];

(* === Field-instance accessors === *)

fieldRank[f_] := $FieldInstances[f]["rank"];
fieldPert[f_] := $FieldInstances[f]["pert"];
fieldJ[f_] := $FieldInstances[f]["J"];
fieldJc[f_] := $FieldInstances[f]["Jc"];

(* Canonical field order: higher rank first, then declaration order *)
sortFields[fields_List] := SortBy[fields, {-fieldRank[#] &, $FieldInstances[#]["order"] &}];

(* Fields appearing in an action, in canonical order *)
detectFields[act_] := sortFields @ Select[Keys[$FieldInstances], !FreeQ[act, #] &];

(* Fields in the theory currently analysed. Set by SourceDeclare / CompSatProp. *)
activeFields = {};

toPerturbationRules[] := Table[
    With[{F = fld, ep = fieldPert[fld]}, F[inds___] :> ep[LI[1], inds]],
    {fld, Keys[$FieldInstances]}
];

(* J -> Jc for every declared field *)
srcToConjRules[] := Table[
    With[{j = fieldJ[fld], jc = fieldJc[fld]},
        If[fieldRank[fld] === 0, j -> jc, j[ind__] :> jc[ind]]
    ],
    {fld, Keys[$FieldInstances]}
];


(* === Spin-parity content ===
   Labels of the projector library carried by each field type, per J^P sector. Labels are only
   unique within a type; a basis element of a sector is therefore the pair {field, label}. *)

$RepLabels = <|
    {"rank3", "Generic"} -> <|
        {3, m} -> {1}, {2, p} -> {1, 2, 3}, {2, m} -> {1, 2}, {1, p} -> {1, 2, 3},
        {1, m} -> {1, 2, 3, 4, 5, 6}, {0, p} -> {1, 2, 3, 4}, {0, m} -> {1}|>,
    {"rank3", "Sym23"} -> <|
        {3, m} -> {1}, {2, p} -> {1, 2}, {2, m} -> {1}, {1, p} -> {1},
        {1, m} -> {1, 2, 4, 5}, {0, p} -> {1, 2, 4}|>,
    {"rank3", "Antisym23"} -> <|
        {2, p} -> {3}, {2, m} -> {2}, {1, p} -> {2, 3},
        {1, m} -> {3, 6}, {0, p} -> {3}, {0, m} -> {1}|>,
    {"rank3", "TotallySymmetric"} -> <|
        {3, m} -> {1}, {2, p} -> {1}, {2, m} -> {1}, {1, p} -> {1},
        {1, m} -> {1, 4}, {0, p} -> {1, 4}|>,
    {"rank3", "TotallyAntisymmetric"} -> <|
        {1, p} -> {3}, {0, m} -> {1}|>,
    {"rank2sym", None} -> <|{2, p} -> {4}, {1, m} -> {7}, {0, p} -> {5, 6}|>,
    {"rank2anti", None} -> <|{1, p} -> {4}, {1, m} -> {9}|>,
    {"vector", None} -> <|{1, m} -> {8}, {0, p} -> {7}|>,
    {"scalar", None} -> <|{0, p} -> {8}|>
|>;

fieldLabels[f_, jp_] := Lookup[
    Lookup[$RepLabels, Key[{$FieldInstances[f]["type"], $FieldInstances[f]["symmetry"]}], <||>],
    Key[jp], {}];

(* List of all sectors *)
listRep = {{3, m}, {2, p}, {2, m}, {1, p}, {1, m}, {0, p}, {0, m}};

(* Basis of a sector: {field, label} pairs, fields in canonical order *)
listC[spin_, parity_] := Flatten[
    Table[{fld, lab}, {fld, activeFields}, {lab, fieldLabels[fld, {spin, parity}]}], 1];

(* Position of a basis element in listC[spin, parity] *)
car[spin_, parity_, lab_] := FirstPosition[listC[spin, parity], lab, {0}, {1}][[1]];

(* Get active sectors (non-empty listC) *)
activeSectors[] := Select[listRep, listC @@ # =!= {} &];

(* Dimension of spin-j representation *)
dimRep[j_] := 2 j + 1;

printConfiguration[] := Module[{},
    PrintK["=== Field Content Configuration ==="];
    PrintK["Active fields: ", activeFields];
    Do[
        PrintK["  ", jp[[1]], If[jp[[2]] === p, "+", "-"], ": basis ", listC @@ jp],
        {jp, activeSectors[]}
    ];
];


(* === Index bookkeeping ===
   Free indices shared by projectors, kinetic blocks and sources. For the pair (row field fi,
   column field fj) the row takes the first rank(fi) indices, the column the next rank(fj). *)

$idx = {m, n, f, r, s, t};

rowInds[fi_] := $idx[[;; fieldRank[fi]]];
colInds[fi_, fj_] := $idx[[fieldRank[fi] + 1 ;; fieldRank[fi] + fieldRank[fj]]];

(* Source / conjugate source with the given upper indices *)
sourceOf[f_, inds_List] := If[fieldRank[f] === 0, fieldJ[f], fieldJ[f] @@ inds];
conjSourceOf[f_, inds_List] := If[fieldRank[f] === 0, fieldJc[f], fieldJc[f] @@ inds];

(* Projector between basis elements {fi, li} and {fj, lj}, all indices down *)
projector[spin_, parity_, {fi_, li_}, {fj_, lj_}] :=
    P[spin, parity, {li, lj}] @@ (-Join[rowInds[fi], colInds[fi, fj]]);


(* === Source component rules ===
   Expensive (ChangeComponents), so computed once per field and cached. Each entry is a list
   of rule sets, applied one after the other by ToComponents. *)

$SourceRules = <||>;

componentRules[src_, 0] := {};

componentRules[src_, 1] := {ChangeComponents[src[-{a, Lorentz}], src[{a, Lorentz}]]};

componentRules[src_, 2] := {
    ChangeComponents[src[{-a, -Lorentz}, {-b, -Lorentz}], src[{a, Lorentz}, {b, Lorentz}]] //
    ReplaceRepeated[#, MomentumRule1] & //
    ReplaceRepeated[#, MetricRule1] &//
    ReplaceRepeated[#, MetricRule2] & // FullSimplify};

componentRules[src_, 3] := Flatten[Table[
    ChangeComponents[
        src[x1*{a, Lorentz}, x2*{b, Lorentz}, x3*{c, Lorentz}],
        src[{a, Lorentz}, {b, Lorentz}, {c, Lorentz}]] //
    ReplaceRepeated[#, MomentumRule1] & //
    ReplaceRepeated[#, MetricRule1] & //
    ReplaceRepeated[#, MetricRule2] & // FullSimplify,
    {x1, {-1, +1}}, {x2, {-1, +1}}, {x3, {-1, +1}}
]];

initializeSourceComponentRules[] := initializeSourceComponentRules[activeFields];

initializeSourceComponentRules[fields_List] := Module[{initialized = {}},
    Do[
        If[!KeyExistsQ[$SourceRules, fld],
            $SourceRules[fld] = Join[
                componentRules[fieldJ[fld], fieldRank[fld]],
                componentRules[fieldJc[fld], fieldRank[fld]]];
            AppendTo[initialized, fld];
        ],
        {fld, fields}
    ];
    If[Length[initialized] > 0,
        PrintK["Source component rules initialized for: ", initialized];,
        PrintK["Source component rules already available for: ", fields];
    ];
];


(* ToComponents: applies index conversion rules (lower -> upper) *)
ToComponents[expr_] := Module[{result},
    result = expr;

    (* Apply source rules, lowest rank first *)
    Do[
        result = result // ReplaceRepeated[#, ruleSet] &,
        {fld, Reverse[activeFields]}, {ruleSet, Lookup[$SourceRules, fld, {}]}
    ];

    (* Apply momentum and metric rules *)
    result = result // ReplaceRepeated[#, MomentumRule1] &;
    result = result // ReplaceRepeated[#, MetricRule1] &;
    result = result // ReplaceRepeated[#, MetricRule2] &;

    result
];

TensorToComponents[expr_] := expr // DummyToBasis[Lorentz] // TraceBasisDummy // ToComponents;

VarAction[ac_, x_] := VarD[x, PD][ac]

ToMomentumSpaceRules = {
    PD[s_][delta[x_, y_]] :> (-I * q[s]) * delta[x, y],
    PD[s_][q[x_]] :> 0,
    delta[x__, z__] :> 1
};


SourceC[x_, y_, {}] := {};

SourceC[spin_, parity_, nullInp_] := Module[{prova, null, labels},
    (* Handle complex conjugate *)
    null = nullInp // Conjugate // ReplaceAll[Conjugate[x_] :> x];
    
    labels = listC[spin, parity];
    If[labels === {}, Return[{}]];
    
    prova = Table[
        Sum[
            null[[car[spin, parity, j]]] * projector[spin, parity, i, j] *
                sourceOf[j[[1]], colInds[i[[1]], j[[1]]]],
            {j, labels}
        ],
        {i, labels}
    ];
    
    prova = prova // ContractMetric // ReplaceAll[dimM -> 4] // FullSimplify;
    prova
];

GaugeC[spin_, parity_, null_] := Module[{prova, labels},
    
    labels = listC[spin, parity];
    If[labels === {}, Return[{}]];
    
    prova = Table[
        Sum[
            null[[car[spin, parity, i]]] * projector[spin, parity, i, j] *
                sourceOf[j[[1]], colInds[i[[1]], j[[1]]]],
            {i, labels}
        ],
        {j, labels}
    ];
    
    prova = prova // ContractMetric // ReplaceAll[dimM -> 4] // FullSimplify;
    prova
];

(* Compute null space for a given sector *)
computeNullSpace[spin_, parity_] := Module[{mat, ns},
    mat = aMatrix[spin, parity];
    If[mat === {} || mat === Null, Return[{}]];
    ns = NullSpace[FullSimplify[mat /. qm -> Sqrt[qq]]] // FullSimplify;
    ns
];

(* Gather all null spaces *)
gatherNullSpaces[] := Module[{},
    PrintK["=== Computing Null Spaces ==="];
    Do[
        {spin, parity} = sector;
        parityStr = If[parity === p, "+", "-"];
        NS[spin, parity] = computeNullSpace[spin, parity];
        If[NS[spin, parity] =!= {},
            PrintK["NS[", spin, "^", parityStr, "]: ", Length[NS[spin, parity]], " null vector(s)"];
        ],
        {sector, activeSectors[]}
    ];
];

(* Gather source constraints (vincoli) for all sectors *)
gatherSourceConstraints[] := Module[{sectorList, vincoliList},
    PrintK["\n=== Computing Source Constraints (Vincoli) ==="];
    
    Do[
        {spin, parity} = sector;
        parityStr = If[parity === p, "+", "-"];
        
        If[NS[spin, parity] =!= {} && NS[spin, parity] =!= Null,
            sourceConstraints[spin, parity] = Table[
                SourceC[spin, parity, nullvec] // ReplaceAll[qm -> Sqrt[qq]] // FullSimplify,
                {nullvec, NS[spin, parity]}
            ];
            PrintK["Sector ", spin, "^", parityStr, ": ", Length[sourceConstraints[spin, parity]], " constraint(s)"];
        ,
            sourceConstraints[spin, parity] = {};
        ],
        {sector, activeSectors[]}
    ];
    
    (* Flatten all constraints into one list *)
    vincoliList = Flatten[Table[
        sourceConstraints[sector[[1]], sector[[2]]],
        {sector, activeSectors[]}
    ]];
    vincoliList
];

(* One-step: gather everything *)
gatherAllConstraints[] := Module[{},
    gatherNullSpaces[];
    vincoli = gatherSourceConstraints[];
    
    PrintK["\n=== Summary ==="];
    PrintK["Total source constraints (vincoli): ", Length[vincoli]];
    
    vincoli
];

(* Least preferred first: higher rank, then later-declared copy, then higher label *)
buildRemovalOrder[spin_, parity_] := SortBy[
    listC[spin, parity],
    {-fieldRank[#[[1]]] &, -$FieldInstances[#[[1]]]["order"] &, -#[[2]] &}
];

fixGauge[spin_, parity_] := Module[
    {mat, fullLabels, removalOrder, targetRank, currentLabels, currentIndices, 
     testLabels, testIndices, testMat, testRank, removed, parityStr},
    
    parityStr = If[parity === p, "+", "-"];
    mat = aMatrix[spin, parity];
    
    (* Handle empty or null matrices *)
    If[mat === {} || mat === Null || mat === 0,
        listCpost[spin, parity] = {};
        sectorPost[{spin, parity}] = {};
        PrintK["Sector ", spin, "^", parityStr, ": empty or null matrix"];
        Return[{}];
    ];
    
    fullLabels = listC[spin, parity];
    targetRank = MatrixRank[mat /. qm -> Sqrt[qq] // FullSimplify];
    
    (* If rank is 0, sector is completely gauged away *)
    If[targetRank === 0,
        listCpost[spin, parity] = {};
        sectorPost[{spin, parity}] = {};
        PrintK["Sector ", spin, "^", parityStr, ": completely gauged away (rank 0)"];
        Return[{}];
    ];
    
    (* If already full rank, no gauge fixing needed *)
    If[targetRank === Length[fullLabels],
        listCpost[spin, parity] = fullLabels;
        sectorPost[{spin, parity}] = mat;
        PrintK["Sector ", spin, "^", parityStr, ": already full rank (", targetRank, ")"];
        Return[mat];
    ];
    
    (* Build removal order: least preferred first *)
    removalOrder = buildRemovalOrder[spin, parity];
    
    (* Start with all labels *)
    currentLabels = fullLabels;
    removed = {};
    
    (* Try to remove each label in order *)
    Do[
        (* Skip if label not in current set *)
        If[Not[MemberQ[currentLabels, labelToRemove]], Continue[]];
        
        (* Test removing this label *)
        testLabels = DeleteCases[currentLabels, labelToRemove];
        testIndices = Map[Position[fullLabels, #][[1, 1]] &, testLabels];
        
        (* Extract submatrix *)
        testMat = mat[[testIndices, testIndices]];
        testRank = MatrixRank[testMat /. qm -> Sqrt[qq] // FullSimplify];
        
        (* If rank is preserved, accept the removal *)
        If[testRank === targetRank,
            currentLabels = testLabels;
            AppendTo[removed, labelToRemove];
        ];
        
        (* Stop if we've reached the target size *)
        If[Length[currentLabels] === targetRank, Break[]];
        
        , {labelToRemove, removalOrder}
    ];
    
    (* Build final submatrix *)
    currentIndices = Map[Position[fullLabels, #][[1, 1]] &, currentLabels];
    sectorPost[{spin, parity}] = mat[[currentIndices, currentIndices]];
    listCpost[spin, parity] = currentLabels;
    
    PrintK["Sector ", spin, "^", parityStr, ": ", Length[fullLabels], " -> ", Length[currentLabels], 
          " (removed: ", removed, ")"];
    
    sectorPost[{spin, parity}]
];

fixAllGauges[] := Module[{activeSecs, nonEmptySectors},
    PrintK["=== Fixing Gauge (Removing Redundancies) ==="];
    PrintK["Preference: keep Scalar > Vector > Rank-2 > Rank-3\n"];
    
    activeSecs = activeSectors[];
    
    Do[
        {spin, parity} = sector;
        fixGauge[spin, parity];
        , {sector, activeSecs}
    ];
    
    (* Update listRep to only include non-empty sectors *)
    nonEmptySectors = Select[activeSecs, listCpost[#[[1]], #[[2]]] =!= {} &];
    listRepPost = nonEmptySectors;
    
    PrintK["\n=== Summary ==="];
    PrintK["Active sectors after gauge fixing: ", listRepPost];
    
    listRepPost
];

extractPoles[] := Module[{dets, massRulesList, allMassRules, masses},
    PrintK["=== Extracting (Candidate) Poles from Gauge-Fixed Matrices ===\n"];
    
    (* Compute determinant for each sector *)
    Do[
        {spin, parity} = sector;
        parityStr = If[parity === p, "+", "-"];
        
        mat = sectorPost[{spin, parity}];
        
        If[mat === {} || mat === Null || Length[mat] === 0,
            detSector[spin, parity] = 1;  (* No contribution *)
            PrintK["Sector ", spin, "^", parityStr, ": empty (det = 1)"];
            Continue[];
        ];
        
        (* For 1x1 matrix, det is just the element *)
        If[Dimensions[mat] === {1, 1},
            detSector[spin, parity] = mat[[1, 1]] // FullSimplify;,
            detSector[spin, parity] = Det[mat] // FullSimplify;
        ];
        
        PrintK["Sector ", spin, "^", parityStr, ": det = ", detSector[spin, parity]];
        
        , {sector, listRepPost}
    ];
    
    (* Solve for poles in each sector *)
    PrintK["\n=== Solving Det = 0 for each sector ===\n"];
    
    massRulesList = {};
    Do[
        {spin, parity} = sector;
        parityStr = If[parity === p, "+", "-"];
        
        det = detSector[spin, parity];
        
        (* Skip trivial determinants *)
        If[det === 1 || det === 0 || FreeQ[det, qq],
            PrintK["Sector ", spin, "^", parityStr, ": no qq dependence, skipping"];
            Continue[];
        ];
        
        (* Solve Det = 0 for qq *)
        sols = Quiet[Solve[det == 0, qq]];
        
        If[sols === {} || sols === {{}},
            PrintK["Sector ", spin, "^", parityStr, ": no solutions found"];,
            PrintK["Sector ", spin, "^", parityStr, ": ", Length[sols], " solution(s)"];
            Do[PrintK["  qq -> ", qq /. sol // FullSimplify], {sol, sols}];
            AppendTo[massRulesList, sols];
        ];
        
        , {sector, listRepPost}
    ];
    
    (* Gather all mass rules *)
    allMassRules = Flatten[massRulesList];
    
    (* Extract unique masses *)
    masses = DeleteDuplicates[(qq /. #) & /@ allMassRules] // FullSimplify;
    
    PrintK["\n=== Mass Spectrum ==="];
    PrintK["Unique poles (qq = m^2):"];
    Do[PrintK["  m^2 = ", mass], {mass, masses}];
    
    masses
];

carp[SPIN_, PARITY_, num_] := Module[{pos},
    pos = Position[listCpost[SPIN, PARITY], num];
    If[pos === {}, 
        PrintK["Warning: label ", num, " not in listCpost[", SPIN, ",", PARITY, "]"];
        Return[0]
    ];
    pos[[1, 1]]
];

saturatedPropagator[SPIN_, PARITY_, subs_:{}] := Module[
    {spin = SPIN, parity = PARITY, mat, bm1, labels, prova, term, parityStr},
    
    parityStr = If[parity === p, "+", "-"];
    
    (* Get the gauge-fixed matrix *)
    mat = sectorPost[{spin, parity}];
    labels = listCpost[spin, parity];
    
    (* Handle empty sector *)
    If[mat === {} || mat === Null || labels === {},
        PrintK["Sector ", spin, "^", parityStr, ": empty, no propagator"];
        Return[0]
    ];
    
    (* Compute inverse *)
    bm1 = Inverse[
        Expand[mat] // ReplaceRepeated[#, subs] & // ReplaceAll[qm -> Sqrt[qq]]
    ];
    
    (* Build saturated propagator: sum over all basis pairs, Jc_i P_ij J_j *)
    prova = Sum[
        bm1[[carp[spin, parity, i], carp[spin, parity, j]]] *
            conjSourceOf[i[[1]], rowInds[i[[1]]]] * projector[spin, parity, i, j] *
            sourceOf[j[[1]], colInds[i[[1]], j[[1]]]],
        {i, labels}, {j, labels}
    ];
    
    (* Simplify *)
    prova = prova // ReplaceAll[dimM -> 4] // CollectTensors // ContractMetric // ReplaceAll[dimM -> 4];
    
    prova
];

computeAllPropagators[subs_:{}] := Module[{},
    PrintK["=== Computing Saturated Propagators ===\n"];
    
    Do[
        {spin, parity} = sector;
        parityStr = If[parity === p, "+", "-"];
        
        PrintK["Computing propagator for sector ", spin, "^", parityStr, "..."];
        propagator[spin, parity] = saturatedPropagator[spin, parity, subs];
        PrintK["  Done.\n"];
        
        , {sector, listRepPost}
    ];
    
    PrintK["=== All propagators computed ==="];
];

fullPropagator[subs_:{}] := Module[{total},
    computeAllPropagators[subs];
    
    total = Sum[
        propagator[sector[[1]], sector[[2]]],
        {sector, listRepPost}
    ];
    
    PrintK["\n=== Full Propagator ==="];
    total
];

lightlikeFrameRules = {
    q[{1, Lorentz}] -> 0, 
    q[{2, Lorentz}] -> 0,
    qq -> q[{0, Lorentz}]^2 - q[{3, Lorentz}]^2
};

momentumRenameRules = {
    q[{0, Lorentz}] :> \[Omega], 
    q[{3, Lorentz}] :> \[Kappa]
};

(* Independent components of a source, with upper Lorentz indices *)
sourceComponents[src_, 0] := {src};

sourceComponents[src_, 1] := Table[src[{k1, Lorentz}], {k1, 0, 3}];

(* For an antisymmetric source, Simplification gives zeros and signs: drop them *)
sourceComponents[src_, 2] := Table[src[{k1, Lorentz}, {k2, Lorentz}], {k1, 0, 3}, {k2, 0, 3}] //
    Simplification // Flatten // ReplaceAll[-x__ :> x] //
    DeleteDuplicates // DeleteCases[#, 0] &;

sourceComponents[src_, 3] := Table[src[{k1, Lorentz}, {k2, Lorentz}, {k3, Lorentz}],
    {k1, 0, 3}, {k2, 0, 3}, {k3, 0, 3}] //
    Simplification // Flatten // ReplaceAll[-x__ :> x] //
    DeleteDuplicates // DeleteCases[#, 0] &;

buildSourceComponentLists[] := Module[{},
    allSourceVars = Flatten[sourceComponents[fieldJ[#], fieldRank[#]] & /@ activeFields];
    allSourceVarsConj = Flatten[sourceComponents[fieldJc[#], fieldRank[#]] & /@ activeFields];
];

convertConstraintsToLightlike[vincoliList_] := Module[{eq1, eq2},
    eq1 = Table[
        ComponentArray[vincoliList[[ii]]] // Expand, 
        {ii, 1, Length[vincoliList]}
    ] // Flatten // DeleteDuplicates//ToComponents // 
    TraceBasisDummy // ToCanonical // ToComponents //
    ReplaceAll[Scalar[x_] :> x] //
    ReplaceAll[lightlikeFrameRules] //
    ReplaceAll[momentumRenameRules] // 
    DeleteDuplicates //
    DeleteCases[#, 0] &;
    
    (* Convert to equation form *)
    eq2 = (# == 0) & /@ eq1;
    eq2 = DeleteCases[eq2, True];
    
    PrintK["Constraints in lightlike frame: ", Length[eq2], " equations"];
    eq2
];

isSingularAtMasslessLimit[sol_] := Module[{expr, testLimit, singular},
    If[sol === {} || sol === {{}}, Return[True]];  (* Empty is "bad" *)
    
    singular = False;
   
     Do[
        expr = FullSimplify[rule[[2]]];  (* RHS of the rule *)
        
        (* Check if denominator contains (\[Omega] - \[Kappa]) *)
(*        If[Not[FreeQ[Denominator[Together[expr]], \[Omega] - \[Kappa]]],
            singular = True; Break[];
        ];
*)      
  
        (* Also try limit test *)
       (* testLimit = Quiet[Limit[expr, \[Omega] -> \[Kappa]]];
        If[testLimit === Indeterminate || testLimit === ComplexInfinity || 
           testLimit === -Infinity || testLimit === Infinity ||
           MatchQ[testLimit, _DirectedInfinity],
            singular = True; Break[];
        ];
       *) 
        
        testLimit = Quiet[Limit[expr, \[Omega] -> \[Kappa]]];

		If[!FreeQ[testLimit, Indeterminate | Infinity | -Infinity | ComplexInfinity | DirectedInfinity], singular = True;
           Break[];
          ];

        
        , {rule, sol}
    ];
    
    singular
];

(* our first attempt

robustSolveMassless[equations_, variables_, maxAttempts_:100] := Module[
    {eqs, vars, nVars, sol, nDependent, subsets, subset, attempt, found, k},
    
    (* Equations are already in form expr == 0, just remove trivial ones *)
    eqs = DeleteCases[equations, True];
    vars = variables;
    nVars = Length[vars];
    
    If[eqs === {},
        PrintK["No constraints - all variables are free."];
        Return[{{}, vars, 0}]  (* {solution, freeVars, nDependent} *)
    ];
    
    (* First attempt: solve for all variables *)
    PrintK["Attempting initial Solve..."];
    sol = Quiet[Solve[eqs, vars]];
    
    If[sol =!= {} && sol =!= {{}} && !isSingularAtMasslessLimit[sol[[1]]],
        (* First attempt is OK *)
        nDependent = Length[sol[[1]]];
        PrintK["Initial solution is OK. Dependent variables: ", nDependent];
        PrintK["Physical DOF: ", nVars - nDependent];
        Return[{sol[[1]], Complement[vars, sol[[1]][[All, 1]]], nDependent}]
    ];
    
    PrintK["Initial solution is BAD (empty or singular). Searching for OK solution..."];
    
    (* Determine starting size: try to get from initial solve *)
    If[sol =!= {} && sol =!= {{}},
        nDependent = Length[sol[[1]]];,
        nDependent = Min[Length[eqs], nVars - 1];  (* Estimate *)
    ];
    
    (* Search from size nDependent down to 1 *)
    found = False;
    attempt = 0;
    
    Do[
        PrintK["Trying subsets of size ", k, "..."];
        subsets = Subsets[vars, {k}];
        
        Do[
            attempt++;
            If[attempt > maxAttempts,
                PrintK["Max attempts reached. Returning best effort."];
                Break[];
            ];
            
            sol = Quiet[Solve[eqs, subset]];
            
            (* Check: non-empty? *)
            If[sol === {} || sol === {{}}, 
                Continue[];  (* This subset doesn't work, try next *)
            ];
            
            (* Check: non-singular? *)
            If[!isSingularAtMasslessLimit[sol[[1]]],
                (* Found OK solution! *)
                found = True;
                PrintK["Found OK solution at attempt ", attempt];
                PrintK["Dependent variables (", k, "): ", subset];
                PrintK["Physical DOF: ", nVars - k];
                Break[];
            ];
            
            , {subset, subsets}
        ];
        
        If[found, Break[]];
        
        , {k, nDependent, 1, -1}
    ];
    
    If[!found,
        PrintK["WARNING: Could not find non-singular solution!"];
        Return[{{}, vars, 0}]
    ];
    
    {sol[[1]], Complement[vars, subset], k}
];

*)
robustSolveMassless[equations_, variables_] := Module[
    {eqs, vars, nVars, sol, nDependent, attempt = 0,
     candidates, newCandidates, bestSol, bestSubset={}, bestK = 0,
     remaining, tried, subset},
    
    eqs = DeleteCases[equations, True|0];
    vars = variables;
    
    
     sol=SolveSystemFull[eqs,vars];
    
    bestSol = sol["Solution"][[1]];
    
    bestSubset = bestSol[[All,1]];
    
   
    {bestSol, Complement[vars, bestSubset], bestK}
]

analyzeMasslessSpectrum[] := Module[{eqsLightlike, result},
    
    PrintK["=== Massless Spectrum Analysis ===\n"];
    
    (* Build source component lists *)
    buildSourceComponentLists[];
    
    (* Convert constraints to lightlike frame *)
    PrintK["\nConverting constraints to lightlike frame..."];
    eqsLightlike = convertConstraintsToLightlike[vincoli];
    
    (* Find non-singular solution *)
    PrintK["\nSolving constraints..."];
    result = robustSolveMassless[eqsLightlike, allSourceVars];
    
    {masslessSolution, masslessFreeVars, masslessNDependent} = result;
    
    (* Compute conjugate solution *)
    PrintK["\nComputing conjugate solution..."];
    masslessSolutionConj = conjugateRules[masslessSolution];
    
   (* PrintK["\n=== Massless Analysis Summary ==="];
    PrintK["Total source components: ", Length[allSourceVars]];
    PrintK["Constrained (dependent): ", masslessNDependent];
    PrintK["Free (physical DOF): ", Length[masslessFreeVars]];
    PrintK["Free variables: ", masslessFreeVars];
  *)  
    (* Return both solutions *)
{    masslessSolution, masslessSolutionConj, masslessFreeVars, masslessNDependent}
];

conjugateSolution[expr_] := 
    Conjugate[expr // Expand] // 
    FullSimplify[#, Assumptions -> {\[Omega] > 0, \[Kappa] > 0}] & //
    ComplexExpand // 
    ReplaceAll[srcToConjRules[]];

(* Convert a list of rules to conjugate rules *)
conjugateRules[rules_List] := Module[{newRules},
    newRules = Table[
        Rule[
            rule[[1]] /. srcToConjRules[],
            conjugateSolution[rule[[2]]]
        ],
        {rule, rules}
    ];
    newRules
];

applyMasslessSolution[prop_, sol_, solConj_] := Module[{propLightlike, propSolved},
    
    (* Convert propagator to lightlike frame *)
    propLightlike = prop // DummyToBasis[Lorentz] // TraceBasisDummy // ToComponents //
        ReplaceAll[Scalar[x_] :> x] //CollectTensors//
        ReplaceAll[lightlikeFrameRules] //
        ReplaceAll[momentumRenameRules];
    
    (* Apply both solutions: sol for J sources, solConj for Jc sources *)
    propSolved = propLightlike // ReplaceAll[sol] // ReplaceAll[solConj];
    
    propSolved
];

restFrameRules = {
    q[{1, Lorentz}] -> 0, 
    q[{2, Lorentz}] -> 0,
    q[{3, Lorentz}] -> 0,
    qq -> \[Omega]^2
};

restFrameMomentumRenameRules = {
    q[{0, Lorentz}] :> \[Omega]
};

convertConstraintsToRestFrame[vincoliList_] := Module[{eq1, eq2},
    eq1 = Table[
        ComponentArray[vincoliList[[ii]]] // Expand, 
        {ii, 1, Length[vincoliList]}
    ] // Flatten // DeleteDuplicates // 
    TraceBasisDummy // ToCanonical // ToComponents //
    ReplaceAll[Scalar[x_] :> x] //
    ReplaceAll[restFrameRules] //
    ReplaceAll[restFrameMomentumRenameRules] // 
    DeleteDuplicates //
    DeleteCases[#, 0] &;
    
    (* Convert to equation form *)
    eq2 = (# == 0) & /@ eq1;
    eq2 = DeleteCases[eq2, True];
    
    PrintK["Constraints in rest frame: ", Length[eq2], " equations"];
    eq2
];

solveMassiveConstraints[equations_, variables_] := Module[
    {eqs, vars, nVars, sol, nDependent, freeVars},
    
    (* Equations are already in form expr == 0, just remove trivial ones *)
    eqs = DeleteCases[equations, True];
    vars = variables;
    nVars = Length[vars];
    
    If[eqs === {},
        PrintK["No constraints - all variables are free."];
        Return[{{}, vars, 0}]
    ];
    
    PrintK["Solving constraints..."];
    sol = Quiet[Solve[eqs, vars]];
    
    If[sol === {} || sol === {{}},
        PrintK["No solution found or empty solution."];
        Return[{{}, vars, 0}]
    ];
    
    nDependent = Length[sol[[1]]];
    freeVars = Complement[vars, sol[[1]][[All, 1]]];
    
    PrintK["Solution found."];
    PrintK["Dependent variables: ", nDependent];
    PrintK["Physical DOF: ", nVars - nDependent];
    
    {sol[[1]], freeVars, nDependent}
];

analyzeMassiveSpectrum[] := Module[{eqsRestFrame, result},
    
    PrintK["=== Massive Spectrum Analysis ===\n"];
    
    (* Rebuild the source component lists for the current fields (a previous model may have left others) *)
    buildSourceComponentLists[];
    
    (* Convert constraints to rest frame *)
    PrintK["\nConverting constraints to rest frame..."];
    eqsRestFrame = convertConstraintsToRestFrame[vincoli];
    
    (* Solve constraints *)
    PrintK["\nSolving constraints..."];
    result = solveMassiveConstraints[eqsRestFrame, allSourceVars];
    
    {massiveSolution, massiveFreeVars, massiveNDependent} = result;
    
    (* Compute conjugate solution *)
    PrintK["\nComputing conjugate solution..."];
    massiveSolutionConj = conjugateRules[massiveSolution];
    
    (*PrintK["\n=== Massive Analysis Summary ==="];
    PrintK["Total source components: ", Length[allSourceVars]];
    PrintK["Constrained (dependent): ", massiveNDependent];
    PrintK["Free (physical DOF): ", Length[massiveFreeVars]];
    PrintK["Free variables: ", massiveFreeVars];
    *)
    (* Return both solutions *)
    {massiveSolution, massiveSolutionConj, massiveFreeVars, massiveNDependent}
];

applyMassiveSolution[prop_, sol_, solConj_] := Module[{propRestFrame, propSolved},
    
    (* Convert propagator to rest frame *)
    propRestFrame = prop // DummyToBasis[Lorentz] // TraceBasisDummy // ToComponents //
        ReplaceAll[Scalar[x_] :> x] //CollectTensors//
        ReplaceAll[restFrameRules] //
        ReplaceAll[restFrameMomentumRenameRules];
    
    (* Apply both solutions: sol for J sources, solConj for Jc sources *)
    propSolved = propRestFrame // ReplaceAll[sol] // ReplaceAll[solConj];
    
    propSolved
];

analyzeAtMassPole[massSq_] := Module[{eqsAtPole, result},
    
    PrintK["=== Analysis at mass pole qq = ", massSq, " ===\n"];
    
    (* Convert constraints to rest frame and substitute specific mass *)
    eqsAtPole = convertConstraintsToRestFrame[vincoli] // 
        ReplaceAll[\[Omega]^2 -> massSq] //
        ReplaceAll[\[Omega] -> Sqrt[massSq]];
    
    (* Solve *)
    result = solveMassiveConstraints[eqsAtPole, allSourceVars];
    
    result
];



(* T[fi, fj]: second variation of the action w.r.t. the perturbations of fi (row, LI[2]) and
   fj (column, LI[3]), in momentum space. Free indices follow rowInds / colInds. *)
computeKineticBlocks[actionPert_] := Module[{},
    Clear[TBlock];
    Do[
        TBlock[fi, fj] = VarAction[
            VarAction[actionPert // ReplaceDummies, fieldPert[fi] @@ Prepend[-rowInds[fi], LI[2]]],
            fieldPert[fj] @@ Prepend[-colInds[fi, fj], LI[3]]] // CollectTensors;
        TBlock[fi, fj] = TBlock[fi, fj] //. ToMomentumSpaceRules,
        {fi, activeFields}, {fj, activeFields}
    ];
];

amatrix2[spin_, parity_] := Module[{prova, labels},
    labels = listC[spin, parity];
    
    If[labels === {},
        PrintK["No labels for sector ", spin, "^", If[parity === p, "+", "-"]];
        Return[{}]
    ];
    
    prova = Table[
        projector[spin, parity, i, j] * TBlock[i[[1]], j[[1]]],
        {i, labels},
        {j, labels}
    ];
    
    (* Normalize by dimension and simplify *)
    prova = prova / dimRep[spin] // ContractMetric // ReplaceAll[dimM -> 4] // FullSimplify;
    prova
];


(*Source declaration in component terms take too long. Do it once for all.*)
SourceDeclare[action_]:=
Block[{Print = (Null &)},
activeFields = detectFields[action];
initializeSourceComponentRules[];
];


CompSatProp[action_]:=Module[{actionPert},

(* Fields not yet seen by SourceDeclare get their component rules here *)
activeFields = detectFields[action];
initializeSourceComponentRules[];

(* View configuration *)
printConfiguration[];


(* Apply to action *)

actionPert = action /. toPerturbationRules[];

computeKineticBlocks[actionPert];

PrintK["\n=== Computing matrices a(\!\(\*SuperscriptBox[\(J\), \(P\)]\)) ==="];

matricesOutput={};

Do[
    {spin, parity} = sector;
    parityStr = If[parity === p, "+", "-"];
    PrintK["\nSector ", spin, "^", parityStr, ":"];
    aMatrix[spin, parity] = amatrix2[spin, parity];
    
    AppendTo[matricesOutput,aMatrix[spin, parity]];
    
    PrintK[MatrixForm[aMatrix[spin, parity]]],
    {sector, activeSectors[]}
];

vincoli = gatherAllConstraints[];
PrintK["\n=== Source Constraints ==="];
Do[PrintK["vincoli[[", k, "]] = ", vincoli[[k]]], {k, Length[vincoli]}];

listRepPost = fixAllGauges[];
PrintK["\n=== Gauge-Fixed Matrices ==="];
Do[
    {spin, parity} = sector;
    parityStr = If[parity === p, "+", "-"];
    PrintK["\nSector ", spin, "^", parityStr, " (labels: ", listCpost[spin, parity], "):"];
    PrintK[MatrixForm[sectorPost[{spin, parity}]]];
    , {sector, listRepPost}
];

masses = extractPoles[];
massesNotZero=DeleteCases[masses,0];

vincoli = Table[
    ToBasis[Lorentz][vincoli[[ii]] // ContractMetric // Simplification // Expand], 
    {ii, 1, Length[vincoli]}
];

PrintK["Source Gauge Constraints converted to Lorentz basis: ", Length[vincoli], " constraints"];

(*listRepPost = fixAllGauges[];
masses = extractPoles[];
*)
computeAllPropagators[];
(* --- For massless poles (qq = 0): --- *)
(* {masslessSolution, masslessSolutionConj, masslessFreeVars, masslessNDependent} = analyzeMasslessSpectrum[]; *)
(*{masslessSolution, masslessSolutionConj, masslessFreeVars, masslessNDependent} = analyzeMasslessSpectrum[]; 
*)(* --- For massive poles (qq = m^2 > 0): --- *)
(*{massiveSolution, massiveSolutionConj, massiveFreeVars, massiveNDependent} = analyzeMassiveSpectrum[]; 
*)
(* --- Or analyze at a specific mass pole: --- *)
(* analyzeAtMassPole[MM^2]; *)

prop=fullPropagator[];

(* Step 1: Constraints *)
(*vincoli = gatherAllConstraints[];
*)
(* Step 2: Convert to Lorentz basis *)
(*vincoli = Table[
    ToBasis[Lorentz][vincoli[[ii]] // ContractMetric // Simplification // Expand], 
    {ii, 1, Length[vincoli]}
];
*)
(* Step 3: Gauge fixing *)
(*listRepPost = fixAllGauges[];
*)
(* Step 4: Mass poles *)
(*masses = extractPoles[];
*)
If[masses=={}, PrintK["No poles detected. Empty (tree-level) spectrum. "];,

(* Step 5: Propagators *)
(*computeAllPropagators[];
*)
(* Step 6a: Massless analysis *)
If[MemberQ[masses, 0],
    PrintK["=== Analyzing massless sector ==="];
    {masslessSolution, masslessSolutionConj, masslessFreeVars, masslessNDependent} = 
        Simplify[analyzeMasslessSpectrum[]];
PrintK[masslessSolution];
];


(* Step 6b: Massive analysis *)
massivePoles = DeleteCases[masses, 0];
If[massivePoles =!= {},
    PrintK["=== Analyzing massive sectors ==="];
    {massiveSolution, massiveSolutionConj, massiveFreeVars, massiveNDependent} = 
        Simplify[analyzeMassiveSpectrum[]];
];


(* Step 7: Apply to propagators *)
    
    (* Massless *)
    If[MemberQ[masses, 0],
        propMasslessLimit = 
            applyMasslessSolution[prop, masslessSolution, masslessSolutionConj];
    ];
    
    (* Massive *)
    If[massivePoles =!= {},
        propMassiveLimit = 
            applyMassiveSolution[prop, massiveSolution, massiveSolutionConj];
    ];
    
    (*this only consider one mass...?*)

];
];


(*Now we compute the spectrum
The massless case is a one-time thing, so the propagator can be stored in LatestPropMassless and printed out upon request. 
The massive case might be a sequence of computations. So we need to store it in a list.

*)



(******************************************************************************************************************)
(******************************************************************************************************************)
(******************************************************************************************************************)
(*********************************************MASSLESS CASE*********************************************************)
(******************************************************************************************************************)
(******************************************************************************************************************)
(******************************************************************************************************************)



ComputeMasslessStates[]:=Module[{step0, step1,Mat},

basisH=DeleteCases[masslessFreeVars/.{},0];
basisHConj=basisH/.srcToConjRules[];

(*Print[{masslessSolution, masslessSolutionConj}];
*)
LatestPropMassless={};

If[Global`DivideEtImpera,

(*DivideEtImpera*)

propMasslessLimitList=If[Head[propMasslessLimit]===Plus, List@@propMasslessLimit, {propMasslessLimit}]; (*a single term is not split into its factors*)

step0=Plus@@Table[Series[propMasslessLimitList[[elemm]]//ReplaceAll[\[Kappa]:>Sqrt[\[Omega]^2-qq]]//Simplification,{qq,0,-1}]//FullSimplify[#]&,{elemm,Length[propMasslessLimitList]}];
LatestPropMassless={step0,qq->0};

If[SeriesCoefficient[step0,-2]=!=0||SeriesCoefficient[step0,-3]=!=0, PrintK["Higher poles present or failure to simplify to zero. Call function 'LatestPropMassless' to access the propagator."];,

step1 = Table[Coefficient[FullSimplify[SeriesCoefficient[step0,-1], Assumptions -> \[Omega]>0]//Expand,source1*source1C],{source1,basisH},{source1C,basisHConj}];

Mat=step1;

{vals, vecs} = Eigensystem[Mat];      (* vecs = eigenvectors *)

masslessData ={0,DeleteCases[vals,0]}]

,
step0=Series[propMasslessLimit//ReplaceAll[\[Kappa]:>Sqrt[\[Omega]^2-qq]]//Simplification,{qq,0,-1}]//FullSimplify[#]&;
LatestPropMassless={step0,qq->0};



If[SeriesCoefficient[step0,-2]=!=0||SeriesCoefficient[step0,-3]=!=0, PrintK["Higher poles present or failure to simplify to zero. Call function 'LatestPropMassless' to access the propagator."];,

step1 = Table[Coefficient[FullSimplify[SeriesCoefficient[step0,-1], Assumptions -> \[Omega]>0]//Expand,source1*source1C],{source1,basisH},{source1C,basisHConj}];

Mat=step1;

{vals, vecs} = Eigensystem[Mat];      (* vecs = eigenvectors *)

masslessData ={0,DeleteCases[vals,0]}]

]]



(******************************************************************************************************************)
(******************************************************************************************************************)
(******************************************************************************************************************)
(*********************************************MASSIVE CASE*********************************************************)
(******************************************************************************************************************)
(******************************************************************************************************************)
(******************************************************************************************************************)


ComputeMassiveStates[]:=Module[{step0, step1,Mat},

basisH=DeleteCases[massiveFreeVars(*/.AJ[x__]:>0*),0];
basisHConj=basisH/.srcToConjRules[];

LatestPropMassive={};
massiveData={};

If[Global`DivideEtImpera,

(*DivideEtImpera*)
propMassiveLimitList=If[Head[propMassiveLimit]===Plus, List@@propMassiveLimit, {propMassiveLimit}]; (*a single term is not split into its factors*)

Table[

PrintK["Computing propagator in the limit q^2 -> ", maxes];

step0=Plus@@Table[Series[propMassiveLimitList[[elemm]]//ReplaceAll[\[Omega]:>Sqrt[qq]]//Simplification,{qq,maxes,-1}],{elemm,Length[propMassiveLimitList]}];
AppendTo[LatestPropMassive,{step0,qq->maxes}];

If[SeriesCoefficient[step0,-2]=!=0||SeriesCoefficient[step0,-3]=!=0, PrintK["Higher poles present or failure to simplify to zero. Call function 'LatestPropMassive' to access the propagator."];,

step1 = Table[Coefficient[FullSimplify[SeriesCoefficient[step0,-1], Assumptions -> \[Omega]>0]//Expand,source1*source1C],{source1,basisH},{source1C,basisHConj}];

Mat=step1;

{vals, vecs} = Eigensystem[Mat];      (* vecs = eigenvectors *)

AppendTo[massiveData,{maxes,DeleteCases[vals,0]}];

];(*endif*)

,{maxes,massesNotZero}];

massiveData
,

Table[
PrintK["Computing propagator in the limit q^2 -> ", maxes];

step0=Series[propMassiveLimit//ReplaceAll[\[Omega]:>Sqrt[qq]]//Simplification,{qq,maxes,-1}];
AppendTo[LatestPropMassive,{step0,qq->maxes}];

If[SeriesCoefficient[step0,-2]=!=0||SeriesCoefficient[step0,-3]=!=0, PrintK["Higher poles present or failure to simplify to zero. Call function 'LatestPropMassive' to access the propagator."];,

step1 = Table[Coefficient[FullSimplify[SeriesCoefficient[step0,-1], Assumptions -> \[Omega]>0]//Expand,source1*source1C],{source1,basisH},{source1C,basisHConj}];

Mat=step1;

{vals, vecs} = Eigensystem[Mat];      (* vecs = eigenvectors *)

AppendTo[massiveData,{maxes,DeleteCases[vals,0]}];
];(*endif*)

,{maxes,massesNotZero}]; 

massiveData]
]



(*Now the addendum to get the propagator in terms of helicity states when massless propagation is present. This bit is not relevant for computing
the unitarity conditions. Yet, it has an educative scope which might turn useful here and there.*)

(*For all type of fields we use, we consider the connection between tensor components and helicity eigenstates. *)

(* Helicity templates, per field type. They are written for the prototype sources
   (VJ, HJ, AJ, SJ, protoBJ) and renamed for each declared field by helicityRules. *)

helicityTemplate["vector", _] := {sv[-1, 1] == I*VJ[{1, Lorentz}] + VJ[{2, Lorentz}], sv[1, 2] == (-I)*VJ[{1, Lorentz}] + VJ[{2, Lorentz}], sv[0, 3] == VJ[{3, Lorentz}], sv[0, 4] == VJ[{0, Lorentz}], svc[-1, 1] == (-I)*VJc[{1, Lorentz}] + VJc[{2, Lorentz}], 
 svc[1, 2] == I*VJc[{1, Lorentz}] + VJc[{2, Lorentz}], svc[0, 3] == VJc[{3, Lorentz}], svc[0, 4] == VJc[{0, Lorentz}]};

helicityTemplate["rank2sym", _] := {sh[-2, 1] == -HJ[{1, Lorentz}, {1, Lorentz}] + (2*I)*HJ[{1, Lorentz}, {2, Lorentz}] + HJ[{2, Lorentz}, {2, Lorentz}], 
 sh[2, 2] == -HJ[{1, Lorentz}, {1, Lorentz}] - (2*I)*HJ[{1, Lorentz}, {2, Lorentz}] + HJ[{2, Lorentz}, {2, Lorentz}], sh[-1, 3] == I*HJ[{1, Lorentz}, {3, Lorentz}] + HJ[{2, Lorentz}, {3, Lorentz}], 
 sh[-1, 4] == I*HJ[{0, Lorentz}, {1, Lorentz}] + HJ[{0, Lorentz}, {2, Lorentz}], sh[1, 5] == (-I)*HJ[{1, Lorentz}, {3, Lorentz}] + HJ[{2, Lorentz}, {3, Lorentz}], 
 sh[1, 6] == (-I)*HJ[{0, Lorentz}, {1, Lorentz}] + HJ[{0, Lorentz}, {2, Lorentz}], sh[0, 7] == HJ[{3, Lorentz}, {3, Lorentz}], sh[0, 8] == HJ[{1, Lorentz}, {1, Lorentz}] + HJ[{2, Lorentz}, {2, Lorentz}], 
 sh[0, 9] == HJ[{0, Lorentz}, {3, Lorentz}], sh[0, 10] == HJ[{0, Lorentz}, {0, Lorentz}], shc[-2, 1] == -HJc[{1, Lorentz}, {1, Lorentz}] - (2*I)*HJc[{1, Lorentz}, {2, Lorentz}] + HJc[{2, Lorentz}, {2, Lorentz}], 
 shc[2, 2] == -HJc[{1, Lorentz}, {1, Lorentz}] + (2*I)*HJc[{1, Lorentz}, {2, Lorentz}] + HJc[{2, Lorentz}, {2, Lorentz}], shc[-1, 3] == (-I)*HJc[{1, Lorentz}, {3, Lorentz}] + HJc[{2, Lorentz}, {3, Lorentz}], 
 shc[-1, 4] == (-I)*HJc[{0, Lorentz}, {1, Lorentz}] + HJc[{0, Lorentz}, {2, Lorentz}], shc[1, 5] == I*HJc[{1, Lorentz}, {3, Lorentz}] + HJc[{2, Lorentz}, {3, Lorentz}], 
 shc[1, 6] == I*HJc[{0, Lorentz}, {1, Lorentz}] + HJc[{0, Lorentz}, {2, Lorentz}], shc[0, 7] == HJc[{3, Lorentz}, {3, Lorentz}], shc[0, 8] == HJc[{1, Lorentz}, {1, Lorentz}] + HJc[{2, Lorentz}, {2, Lorentz}], 
 shc[0, 9] == HJc[{0, Lorentz}, {3, Lorentz}], shc[0, 10] == HJc[{0, Lorentz}, {0, Lorentz}]};



helicityTemplate["rank3", sym_] := 
 Switch[sym,
        
        "Generic",
        {sa[-3, 1] == (-I)*AJ[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] - AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] - AJ[{1, Lorentz}, {2, Lorentz}, {1, Lorentz}] + I*AJ[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] - 
   AJ[{2, Lorentz}, {1, Lorentz}, {1, Lorentz}] + I*AJ[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}] + I*AJ[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[3, 2] == I*AJ[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] - AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] - AJ[{1, Lorentz}, {2, Lorentz}, {1, Lorentz}] - I*AJ[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] - 
   AJ[{2, Lorentz}, {1, Lorentz}, {1, Lorentz}] - I*AJ[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}] - I*AJ[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[-2, 3] == -AJ[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] + I*AJ[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] + I*AJ[{2, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sa[-2, 4] == -AJ[{3, Lorentz}, {1, Lorentz}, {1, Lorentz}] + I*AJ[{3, Lorentz}, {1, Lorentz}, {2, Lorentz}] + I*AJ[{3, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJ[{3, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[-2, 5] == -AJ[{1, Lorentz}, {3, Lorentz}, {1, Lorentz}] + I*AJ[{1, Lorentz}, {3, Lorentz}, {2, Lorentz}] + I*AJ[{2, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {3, Lorentz}, {2, Lorentz}], 
 sa[-2, 6] == -AJ[{1, Lorentz}, {0, Lorentz}, {1, Lorentz}] + I*AJ[{1, Lorentz}, {0, Lorentz}, {2, Lorentz}] + I*AJ[{2, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sa[-2, 7] == -AJ[{0, Lorentz}, {1, Lorentz}, {1, Lorentz}] + I*AJ[{0, Lorentz}, {1, Lorentz}, {2, Lorentz}] + I*AJ[{0, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[-2, 8] == -AJ[{1, Lorentz}, {1, Lorentz}, {0, Lorentz}] + I*AJ[{1, Lorentz}, {2, Lorentz}, {0, Lorentz}] + I*AJ[{2, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {0, Lorentz}], 
 sa[2, 9] == -AJ[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] - I*AJ[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] - I*AJ[{2, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sa[2, 10] == -AJ[{3, Lorentz}, {1, Lorentz}, {1, Lorentz}] - I*AJ[{3, Lorentz}, {1, Lorentz}, {2, Lorentz}] - I*AJ[{3, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJ[{3, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[2, 11] == -AJ[{1, Lorentz}, {3, Lorentz}, {1, Lorentz}] - I*AJ[{1, Lorentz}, {3, Lorentz}, {2, Lorentz}] - I*AJ[{2, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {3, Lorentz}, {2, Lorentz}], 
 sa[2, 12] == -AJ[{1, Lorentz}, {0, Lorentz}, {1, Lorentz}] - I*AJ[{1, Lorentz}, {0, Lorentz}, {2, Lorentz}] - I*AJ[{2, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sa[2, 13] == -AJ[{0, Lorentz}, {1, Lorentz}, {1, Lorentz}] - I*AJ[{0, Lorentz}, {1, Lorentz}, {2, Lorentz}] - I*AJ[{0, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[2, 14] == -AJ[{1, Lorentz}, {1, Lorentz}, {0, Lorentz}] - I*AJ[{1, Lorentz}, {2, Lorentz}, {0, Lorentz}] - I*AJ[{2, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {0, Lorentz}], 
 sa[-1, 15] == I*AJ[{3, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{3, Lorentz}, {2, Lorentz}, {3, Lorentz}], sa[-1, 16] == I*AJ[{1, Lorentz}, {3, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {3, Lorentz}, {3, Lorentz}], 
 sa[-1, 17] == I*AJ[{1, Lorentz}, {0, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {0, Lorentz}, {3, Lorentz}], sa[-1, 18] == I*AJ[{0, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sa[-1, 19] == I*AJ[{3, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJ[{3, Lorentz}, {3, Lorentz}, {2, Lorentz}], sa[-1, 20] == I*AJ[{3, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{3, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sa[-1, 21] == I*AJ[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + I*AJ[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[-1, 22] == I*AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] - I*AJ[{1, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}] - AJ[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}], 
 sa[-1, 23] == I*AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJ[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] - I*AJ[{2, Lorentz}, {1, Lorentz}, {1, Lorentz}] - AJ[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}], 
 sa[-1, 24] == I*AJ[{0, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJ[{0, Lorentz}, {3, Lorentz}, {2, Lorentz}], sa[-1, 25] == I*AJ[{0, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{0, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sa[-1, 26] == I*AJ[{3, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJ[{3, Lorentz}, {2, Lorentz}, {0, Lorentz}], sa[-1, 27] == I*AJ[{1, Lorentz}, {3, Lorentz}, {0, Lorentz}] + AJ[{2, Lorentz}, {3, Lorentz}, {0, Lorentz}], 
 sa[-1, 28] == I*AJ[{1, Lorentz}, {0, Lorentz}, {0, Lorentz}] + AJ[{2, Lorentz}, {0, Lorentz}, {0, Lorentz}], sa[-1, 29] == I*AJ[{0, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {0, Lorentz}], 
 sa[1, 30] == (-I)*AJ[{3, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{3, Lorentz}, {2, Lorentz}, {3, Lorentz}], sa[1, 31] == (-I)*AJ[{1, Lorentz}, {3, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {3, Lorentz}, {3, Lorentz}], 
 sa[1, 32] == (-I)*AJ[{1, Lorentz}, {0, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {0, Lorentz}, {3, Lorentz}], sa[1, 33] == (-I)*AJ[{0, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sa[1, 34] == (-I)*AJ[{3, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJ[{3, Lorentz}, {3, Lorentz}, {2, Lorentz}], sa[1, 35] == (-I)*AJ[{3, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{3, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sa[1, 36] == (-I)*AJ[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] - I*AJ[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[1, 37] == (-I)*AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + I*AJ[{1, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}] - AJ[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}], 
 sa[1, 38] == (-I)*AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJ[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] + I*AJ[{2, Lorentz}, {1, Lorentz}, {1, Lorentz}] - AJ[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}], 
 sa[1, 39] == (-I)*AJ[{0, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJ[{0, Lorentz}, {3, Lorentz}, {2, Lorentz}], sa[1, 40] == (-I)*AJ[{0, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{0, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sa[1, 41] == (-I)*AJ[{3, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJ[{3, Lorentz}, {2, Lorentz}, {0, Lorentz}], sa[1, 42] == (-I)*AJ[{1, Lorentz}, {3, Lorentz}, {0, Lorentz}] + AJ[{2, Lorentz}, {3, Lorentz}, {0, Lorentz}], 
 sa[1, 43] == (-I)*AJ[{1, Lorentz}, {0, Lorentz}, {0, Lorentz}] + AJ[{2, Lorentz}, {0, Lorentz}, {0, Lorentz}], sa[1, 44] == (-I)*AJ[{0, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {0, Lorentz}], 
 sa[0, 45] == AJ[{3, Lorentz}, {3, Lorentz}, {3, Lorentz}], sa[0, 46] == AJ[{3, Lorentz}, {0, Lorentz}, {3, Lorentz}], sa[0, 47] == AJ[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sa[0, 48] == -AJ[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {1, Lorentz}, {3, Lorentz}], sa[0, 49] == AJ[{0, Lorentz}, {3, Lorentz}, {3, Lorentz}], sa[0, 50] == AJ[{0, Lorentz}, {0, Lorentz}, {3, Lorentz}], 
 sa[0, 51] == AJ[{3, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJ[{3, Lorentz}, {2, Lorentz}, {2, Lorentz}], sa[0, 52] == AJ[{3, Lorentz}, {1, Lorentz}, {2, Lorentz}] - AJ[{3, Lorentz}, {2, Lorentz}, {1, Lorentz}], 
 sa[0, 53] == AJ[{1, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {3, Lorentz}, {2, Lorentz}], sa[0, 54] == AJ[{1, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sa[0, 55] == AJ[{1, Lorentz}, {3, Lorentz}, {2, Lorentz}] - AJ[{2, Lorentz}, {3, Lorentz}, {1, Lorentz}], sa[0, 56] == AJ[{1, Lorentz}, {0, Lorentz}, {2, Lorentz}] - AJ[{2, Lorentz}, {0, Lorentz}, {1, Lorentz}], 
 sa[0, 57] == AJ[{0, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {2, Lorentz}], sa[0, 58] == AJ[{0, Lorentz}, {1, Lorentz}, {2, Lorentz}] - AJ[{0, Lorentz}, {2, Lorentz}, {1, Lorentz}], 
 sa[0, 59] == AJ[{3, Lorentz}, {3, Lorentz}, {0, Lorentz}], sa[0, 60] == AJ[{3, Lorentz}, {0, Lorentz}, {0, Lorentz}], sa[0, 61] == AJ[{1, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {0, Lorentz}], 
 sa[0, 62] == -AJ[{1, Lorentz}, {2, Lorentz}, {0, Lorentz}] + AJ[{2, Lorentz}, {1, Lorentz}, {0, Lorentz}], sa[0, 63] == AJ[{0, Lorentz}, {3, Lorentz}, {0, Lorentz}], sa[0, 64] == AJ[{0, Lorentz}, {0, Lorentz}, {0, Lorentz}], 
 sac[-3, 1] == I*AJc[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] - AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] - AJc[{1, Lorentz}, {2, Lorentz}, {1, Lorentz}] - I*AJc[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] - 
   AJc[{2, Lorentz}, {1, Lorentz}, {1, Lorentz}] - I*AJc[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}] - I*AJc[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[3, 2] == (-I)*AJc[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] - AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] - AJc[{1, Lorentz}, {2, Lorentz}, {1, Lorentz}] + I*AJc[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] - 
   AJc[{2, Lorentz}, {1, Lorentz}, {1, Lorentz}] + I*AJc[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}] + I*AJc[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[-2, 3] == -AJc[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] - I*AJc[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] - I*AJc[{2, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sac[-2, 4] == -AJc[{3, Lorentz}, {1, Lorentz}, {1, Lorentz}] - I*AJc[{3, Lorentz}, {1, Lorentz}, {2, Lorentz}] - I*AJc[{3, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJc[{3, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[-2, 5] == -AJc[{1, Lorentz}, {3, Lorentz}, {1, Lorentz}] - I*AJc[{1, Lorentz}, {3, Lorentz}, {2, Lorentz}] - I*AJc[{2, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {3, Lorentz}, {2, Lorentz}], 
 sac[-2, 6] == -AJc[{1, Lorentz}, {0, Lorentz}, {1, Lorentz}] - I*AJc[{1, Lorentz}, {0, Lorentz}, {2, Lorentz}] - I*AJc[{2, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sac[-2, 7] == -AJc[{0, Lorentz}, {1, Lorentz}, {1, Lorentz}] - I*AJc[{0, Lorentz}, {1, Lorentz}, {2, Lorentz}] - I*AJc[{0, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[-2, 8] == -AJc[{1, Lorentz}, {1, Lorentz}, {0, Lorentz}] - I*AJc[{1, Lorentz}, {2, Lorentz}, {0, Lorentz}] - I*AJc[{2, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {0, Lorentz}], 
 sac[2, 9] == -AJc[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] + I*AJc[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] + I*AJc[{2, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sac[2, 10] == -AJc[{3, Lorentz}, {1, Lorentz}, {1, Lorentz}] + I*AJc[{3, Lorentz}, {1, Lorentz}, {2, Lorentz}] + I*AJc[{3, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJc[{3, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[2, 11] == -AJc[{1, Lorentz}, {3, Lorentz}, {1, Lorentz}] + I*AJc[{1, Lorentz}, {3, Lorentz}, {2, Lorentz}] + I*AJc[{2, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {3, Lorentz}, {2, Lorentz}], 
 sac[2, 12] == -AJc[{1, Lorentz}, {0, Lorentz}, {1, Lorentz}] + I*AJc[{1, Lorentz}, {0, Lorentz}, {2, Lorentz}] + I*AJc[{2, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sac[2, 13] == -AJc[{0, Lorentz}, {1, Lorentz}, {1, Lorentz}] + I*AJc[{0, Lorentz}, {1, Lorentz}, {2, Lorentz}] + I*AJc[{0, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[2, 14] == -AJc[{1, Lorentz}, {1, Lorentz}, {0, Lorentz}] + I*AJc[{1, Lorentz}, {2, Lorentz}, {0, Lorentz}] + I*AJc[{2, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {0, Lorentz}], 
 sac[-1, 15] == (-I)*AJc[{3, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{3, Lorentz}, {2, Lorentz}, {3, Lorentz}], sac[-1, 16] == (-I)*AJc[{1, Lorentz}, {3, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {3, Lorentz}, {3, Lorentz}], 
 sac[-1, 17] == (-I)*AJc[{1, Lorentz}, {0, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {0, Lorentz}, {3, Lorentz}], sac[-1, 18] == (-I)*AJc[{0, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sac[-1, 19] == (-I)*AJc[{3, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJc[{3, Lorentz}, {3, Lorentz}, {2, Lorentz}], sac[-1, 20] == (-I)*AJc[{3, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{3, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sac[-1, 21] == (-I)*AJc[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] - I*AJc[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[-1, 22] == (-I)*AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + I*AJc[{1, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}] - AJc[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}], 
 sac[-1, 23] == (-I)*AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJc[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] + I*AJc[{2, Lorentz}, {1, Lorentz}, {1, Lorentz}] - AJc[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}], 
 sac[-1, 24] == (-I)*AJc[{0, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJc[{0, Lorentz}, {3, Lorentz}, {2, Lorentz}], sac[-1, 25] == (-I)*AJc[{0, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{0, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sac[-1, 26] == (-I)*AJc[{3, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJc[{3, Lorentz}, {2, Lorentz}, {0, Lorentz}], sac[-1, 27] == (-I)*AJc[{1, Lorentz}, {3, Lorentz}, {0, Lorentz}] + AJc[{2, Lorentz}, {3, Lorentz}, {0, Lorentz}], 
 sac[-1, 28] == (-I)*AJc[{1, Lorentz}, {0, Lorentz}, {0, Lorentz}] + AJc[{2, Lorentz}, {0, Lorentz}, {0, Lorentz}], sac[-1, 29] == (-I)*AJc[{0, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {0, Lorentz}], 
 sac[1, 30] == I*AJc[{3, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{3, Lorentz}, {2, Lorentz}, {3, Lorentz}], sac[1, 31] == I*AJc[{1, Lorentz}, {3, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {3, Lorentz}, {3, Lorentz}], 
 sac[1, 32] == I*AJc[{1, Lorentz}, {0, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {0, Lorentz}, {3, Lorentz}], sac[1, 33] == I*AJc[{0, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sac[1, 34] == I*AJc[{3, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJc[{3, Lorentz}, {3, Lorentz}, {2, Lorentz}], sac[1, 35] == I*AJc[{3, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{3, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sac[1, 36] == I*AJc[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + I*AJc[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[1, 37] == I*AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] - I*AJc[{1, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}] - AJc[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}], 
 sac[1, 38] == I*AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJc[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] - I*AJc[{2, Lorentz}, {1, Lorentz}, {1, Lorentz}] - AJc[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}], 
 sac[1, 39] == I*AJc[{0, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJc[{0, Lorentz}, {3, Lorentz}, {2, Lorentz}], sac[1, 40] == I*AJc[{0, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{0, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sac[1, 41] == I*AJc[{3, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJc[{3, Lorentz}, {2, Lorentz}, {0, Lorentz}], sac[1, 42] == I*AJc[{1, Lorentz}, {3, Lorentz}, {0, Lorentz}] + AJc[{2, Lorentz}, {3, Lorentz}, {0, Lorentz}], 
 sac[1, 43] == I*AJc[{1, Lorentz}, {0, Lorentz}, {0, Lorentz}] + AJc[{2, Lorentz}, {0, Lorentz}, {0, Lorentz}], sac[1, 44] == I*AJc[{0, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {0, Lorentz}], 
 sac[0, 45] == AJc[{3, Lorentz}, {3, Lorentz}, {3, Lorentz}], sac[0, 46] == AJc[{3, Lorentz}, {0, Lorentz}, {3, Lorentz}], sac[0, 47] == AJc[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sac[0, 48] == -AJc[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {1, Lorentz}, {3, Lorentz}], sac[0, 49] == AJc[{0, Lorentz}, {3, Lorentz}, {3, Lorentz}], sac[0, 50] == AJc[{0, Lorentz}, {0, Lorentz}, {3, Lorentz}], 
 sac[0, 51] == AJc[{3, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJc[{3, Lorentz}, {2, Lorentz}, {2, Lorentz}], sac[0, 52] == AJc[{3, Lorentz}, {1, Lorentz}, {2, Lorentz}] - AJc[{3, Lorentz}, {2, Lorentz}, {1, Lorentz}], 
 sac[0, 53] == AJc[{1, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {3, Lorentz}, {2, Lorentz}], sac[0, 54] == AJc[{1, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sac[0, 55] == AJc[{1, Lorentz}, {3, Lorentz}, {2, Lorentz}] - AJc[{2, Lorentz}, {3, Lorentz}, {1, Lorentz}], sac[0, 56] == AJc[{1, Lorentz}, {0, Lorentz}, {2, Lorentz}] - AJc[{2, Lorentz}, {0, Lorentz}, {1, Lorentz}], 
 sac[0, 57] == AJc[{0, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {2, Lorentz}], sac[0, 58] == AJc[{0, Lorentz}, {1, Lorentz}, {2, Lorentz}] - AJc[{0, Lorentz}, {2, Lorentz}, {1, Lorentz}], 
 sac[0, 59] == AJc[{3, Lorentz}, {3, Lorentz}, {0, Lorentz}], sac[0, 60] == AJc[{3, Lorentz}, {0, Lorentz}, {0, Lorentz}], sac[0, 61] == AJc[{1, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {0, Lorentz}], 
 sac[0, 62] == -AJc[{1, Lorentz}, {2, Lorentz}, {0, Lorentz}] + AJc[{2, Lorentz}, {1, Lorentz}, {0, Lorentz}], sac[0, 63] == AJc[{0, Lorentz}, {3, Lorentz}, {0, Lorentz}], sac[0, 64] == AJc[{0, Lorentz}, {0, Lorentz}, {0, Lorentz}]},
        
        "Sym23",
        {sa[-3, 1] == (-I)*AJ[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] - 2*AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + I*AJ[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] - AJ[{2, Lorentz}, {1, Lorentz}, {1, Lorentz}] + 
   (2*I)*AJ[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[3, 2] == I*AJ[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] - 2*AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] - I*AJ[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] - AJ[{2, Lorentz}, {1, Lorentz}, {1, Lorentz}] - 
   (2*I)*AJ[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[-2, 3] == -AJ[{3, Lorentz}, {1, Lorentz}, {1, Lorentz}] + (2*I)*AJ[{3, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJ[{3, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[-2, 4] == -AJ[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] + I*AJ[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] + I*AJ[{2, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sa[-2, 5] == -AJ[{0, Lorentz}, {1, Lorentz}, {1, Lorentz}] + (2*I)*AJ[{0, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[-2, 6] == -AJ[{1, Lorentz}, {0, Lorentz}, {1, Lorentz}] + I*AJ[{1, Lorentz}, {0, Lorentz}, {2, Lorentz}] + I*AJ[{2, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sa[2, 7] == -AJ[{3, Lorentz}, {1, Lorentz}, {1, Lorentz}] - (2*I)*AJ[{3, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJ[{3, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[2, 8] == -AJ[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] - I*AJ[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] - I*AJ[{2, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sa[2, 9] == -AJ[{0, Lorentz}, {1, Lorentz}, {1, Lorentz}] - (2*I)*AJ[{0, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[2, 10] == -AJ[{1, Lorentz}, {0, Lorentz}, {1, Lorentz}] - I*AJ[{1, Lorentz}, {0, Lorentz}, {2, Lorentz}] - I*AJ[{2, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sa[-1, 11] == I*AJ[{1, Lorentz}, {3, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {3, Lorentz}, {3, Lorentz}], sa[-1, 12] == I*AJ[{3, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{3, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sa[-1, 13] == I*AJ[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + I*AJ[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[-1, 14] == I*AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJ[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] - I*AJ[{2, Lorentz}, {1, Lorentz}, {1, Lorentz}] - AJ[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}], 
 sa[-1, 15] == I*AJ[{0, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {3, Lorentz}], sa[-1, 16] == I*AJ[{3, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{3, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sa[-1, 17] == I*AJ[{1, Lorentz}, {0, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {0, Lorentz}, {3, Lorentz}], sa[-1, 18] == I*AJ[{1, Lorentz}, {0, Lorentz}, {0, Lorentz}] + AJ[{2, Lorentz}, {0, Lorentz}, {0, Lorentz}], 
 sa[-1, 19] == I*AJ[{0, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{0, Lorentz}, {0, Lorentz}, {2, Lorentz}], sa[1, 20] == (-I)*AJ[{1, Lorentz}, {3, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {3, Lorentz}, {3, Lorentz}], 
 sa[1, 21] == (-I)*AJ[{3, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{3, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sa[1, 22] == (-I)*AJ[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] - I*AJ[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[1, 23] == (-I)*AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJ[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] + I*AJ[{2, Lorentz}, {1, Lorentz}, {1, Lorentz}] - AJ[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}], 
 sa[1, 24] == (-I)*AJ[{0, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {3, Lorentz}], sa[1, 25] == (-I)*AJ[{3, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{3, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sa[1, 26] == (-I)*AJ[{1, Lorentz}, {0, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {0, Lorentz}, {3, Lorentz}], sa[1, 27] == (-I)*AJ[{1, Lorentz}, {0, Lorentz}, {0, Lorentz}] + AJ[{2, Lorentz}, {0, Lorentz}, {0, Lorentz}], 
 sa[1, 28] == (-I)*AJ[{0, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{0, Lorentz}, {0, Lorentz}, {2, Lorentz}], sa[0, 29] == AJ[{3, Lorentz}, {3, Lorentz}, {3, Lorentz}], sa[0, 30] == AJ[{0, Lorentz}, {3, Lorentz}, {3, Lorentz}], 
 sa[0, 31] == AJ[{3, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJ[{3, Lorentz}, {2, Lorentz}, {2, Lorentz}], sa[0, 32] == AJ[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sa[0, 33] == AJ[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] - AJ[{2, Lorentz}, {1, Lorentz}, {3, Lorentz}], sa[0, 34] == AJ[{0, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[0, 35] == AJ[{3, Lorentz}, {0, Lorentz}, {3, Lorentz}], sa[0, 36] == AJ[{3, Lorentz}, {0, Lorentz}, {0, Lorentz}], sa[0, 37] == AJ[{1, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sa[0, 38] == -AJ[{1, Lorentz}, {0, Lorentz}, {2, Lorentz}] + AJ[{2, Lorentz}, {0, Lorentz}, {1, Lorentz}], sa[0, 39] == AJ[{0, Lorentz}, {0, Lorentz}, {3, Lorentz}], sa[0, 40] == AJ[{0, Lorentz}, {0, Lorentz}, {0, Lorentz}], 
 sac[-3, 1] == I*AJc[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] - 2*AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] - I*AJc[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] - AJc[{2, Lorentz}, {1, Lorentz}, {1, Lorentz}] - 
   (2*I)*AJc[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[3, 2] == (-I)*AJc[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] - 2*AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + I*AJc[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] - AJc[{2, Lorentz}, {1, Lorentz}, {1, Lorentz}] + 
   (2*I)*AJc[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[-2, 3] == -AJc[{3, Lorentz}, {1, Lorentz}, {1, Lorentz}] - (2*I)*AJc[{3, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJc[{3, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[-2, 4] == -AJc[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] - I*AJc[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] - I*AJc[{2, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sac[-2, 5] == -AJc[{0, Lorentz}, {1, Lorentz}, {1, Lorentz}] - (2*I)*AJc[{0, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[-2, 6] == -AJc[{1, Lorentz}, {0, Lorentz}, {1, Lorentz}] - I*AJc[{1, Lorentz}, {0, Lorentz}, {2, Lorentz}] - I*AJc[{2, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sac[2, 7] == -AJc[{3, Lorentz}, {1, Lorentz}, {1, Lorentz}] + (2*I)*AJc[{3, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJc[{3, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[2, 8] == -AJc[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] + I*AJc[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] + I*AJc[{2, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sac[2, 9] == -AJc[{0, Lorentz}, {1, Lorentz}, {1, Lorentz}] + (2*I)*AJc[{0, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[2, 10] == -AJc[{1, Lorentz}, {0, Lorentz}, {1, Lorentz}] + I*AJc[{1, Lorentz}, {0, Lorentz}, {2, Lorentz}] + I*AJc[{2, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sac[-1, 11] == (-I)*AJc[{1, Lorentz}, {3, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {3, Lorentz}, {3, Lorentz}], sac[-1, 12] == (-I)*AJc[{3, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{3, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sac[-1, 13] == (-I)*AJc[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] - I*AJc[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[-1, 14] == (-I)*AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJc[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] + I*AJc[{2, Lorentz}, {1, Lorentz}, {1, Lorentz}] - AJc[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}], 
 sac[-1, 15] == (-I)*AJc[{0, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {3, Lorentz}], sac[-1, 16] == (-I)*AJc[{3, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{3, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sac[-1, 17] == (-I)*AJc[{1, Lorentz}, {0, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {0, Lorentz}, {3, Lorentz}], sac[-1, 18] == (-I)*AJc[{1, Lorentz}, {0, Lorentz}, {0, Lorentz}] + AJc[{2, Lorentz}, {0, Lorentz}, {0, Lorentz}], 
 sac[-1, 19] == (-I)*AJc[{0, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{0, Lorentz}, {0, Lorentz}, {2, Lorentz}], sac[1, 20] == I*AJc[{1, Lorentz}, {3, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {3, Lorentz}, {3, Lorentz}], 
 sac[1, 21] == I*AJc[{3, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{3, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sac[1, 22] == I*AJc[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + I*AJc[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[1, 23] == I*AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJc[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] - I*AJc[{2, Lorentz}, {1, Lorentz}, {1, Lorentz}] - AJc[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}], 
 sac[1, 24] == I*AJc[{0, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {3, Lorentz}], sac[1, 25] == I*AJc[{3, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{3, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sac[1, 26] == I*AJc[{1, Lorentz}, {0, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {0, Lorentz}, {3, Lorentz}], sac[1, 27] == I*AJc[{1, Lorentz}, {0, Lorentz}, {0, Lorentz}] + AJc[{2, Lorentz}, {0, Lorentz}, {0, Lorentz}], 
 sac[1, 28] == I*AJc[{0, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{0, Lorentz}, {0, Lorentz}, {2, Lorentz}], sac[0, 29] == AJc[{3, Lorentz}, {3, Lorentz}, {3, Lorentz}], 
 sac[0, 30] == AJc[{0, Lorentz}, {3, Lorentz}, {3, Lorentz}], sac[0, 31] == AJc[{3, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJc[{3, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[0, 32] == AJc[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], sac[0, 33] == AJc[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] - AJc[{2, Lorentz}, {1, Lorentz}, {3, Lorentz}], 
 sac[0, 34] == AJc[{0, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {2, Lorentz}], sac[0, 35] == AJc[{3, Lorentz}, {0, Lorentz}, {3, Lorentz}], sac[0, 36] == AJc[{3, Lorentz}, {0, Lorentz}, {0, Lorentz}], 
 sac[0, 37] == AJc[{1, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {0, Lorentz}, {2, Lorentz}], sac[0, 38] == -AJc[{1, Lorentz}, {0, Lorentz}, {2, Lorentz}] + AJc[{2, Lorentz}, {0, Lorentz}, {1, Lorentz}], 
 sac[0, 39] == AJc[{0, Lorentz}, {0, Lorentz}, {3, Lorentz}], sac[0, 40] == AJc[{0, Lorentz}, {0, Lorentz}, {0, Lorentz}]},
        
        "Antisym23",
        {sa[-2, 1] == -AJ[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] + I*AJ[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] + I*AJ[{2, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sa[-2, 2] == -AJ[{1, Lorentz}, {0, Lorentz}, {1, Lorentz}] + I*AJ[{1, Lorentz}, {0, Lorentz}, {2, Lorentz}] + I*AJ[{2, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sa[2, 3] == -AJ[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] - I*AJ[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] - I*AJ[{2, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sa[2, 4] == -AJ[{1, Lorentz}, {0, Lorentz}, {1, Lorentz}] - I*AJ[{1, Lorentz}, {0, Lorentz}, {2, Lorentz}] - I*AJ[{2, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sa[-1, 5] == I*AJ[{3, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{3, Lorentz}, {2, Lorentz}, {3, Lorentz}], sa[-1, 6] == I*AJ[{1, Lorentz}, {0, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {0, Lorentz}, {3, Lorentz}], 
 sa[-1, 7] == I*AJ[{0, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {3, Lorentz}], sa[-1, 8] == I*AJ[{3, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{3, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sa[-1, 9] == I*AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJ[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}], sa[-1, 10] == I*AJ[{0, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{0, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sa[1, 11] == (-I)*AJ[{3, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{3, Lorentz}, {2, Lorentz}, {3, Lorentz}], sa[1, 12] == (-I)*AJ[{1, Lorentz}, {0, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {0, Lorentz}, {3, Lorentz}], 
 sa[1, 13] == (-I)*AJ[{0, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {3, Lorentz}], sa[1, 14] == (-I)*AJ[{3, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{3, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sa[1, 15] == (-I)*AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJ[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}], sa[1, 16] == (-I)*AJ[{0, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{0, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sa[0, 17] == AJ[{3, Lorentz}, {0, Lorentz}, {3, Lorentz}], sa[0, 18] == AJ[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sa[0, 19] == -AJ[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {1, Lorentz}, {3, Lorentz}], sa[0, 20] == AJ[{0, Lorentz}, {0, Lorentz}, {3, Lorentz}], sa[0, 21] == AJ[{3, Lorentz}, {1, Lorentz}, {2, Lorentz}], 
 sa[0, 22] == AJ[{1, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {0, Lorentz}, {2, Lorentz}], sa[0, 23] == AJ[{1, Lorentz}, {0, Lorentz}, {2, Lorentz}] - AJ[{2, Lorentz}, {0, Lorentz}, {1, Lorentz}], 
 sa[0, 24] == AJ[{0, Lorentz}, {1, Lorentz}, {2, Lorentz}], sac[-2, 1] == -AJc[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] - I*AJc[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] - I*AJc[{2, Lorentz}, {1, Lorentz}, {3, Lorentz}] + 
   AJc[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], sac[-2, 2] == -AJc[{1, Lorentz}, {0, Lorentz}, {1, Lorentz}] - I*AJc[{1, Lorentz}, {0, Lorentz}, {2, Lorentz}] - I*AJc[{2, Lorentz}, {0, Lorentz}, {1, Lorentz}] + 
   AJc[{2, Lorentz}, {0, Lorentz}, {2, Lorentz}], sac[2, 3] == -AJc[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] + I*AJc[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] + I*AJc[{2, Lorentz}, {1, Lorentz}, {3, Lorentz}] + 
   AJc[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], sac[2, 4] == -AJc[{1, Lorentz}, {0, Lorentz}, {1, Lorentz}] + I*AJc[{1, Lorentz}, {0, Lorentz}, {2, Lorentz}] + I*AJc[{2, Lorentz}, {0, Lorentz}, {1, Lorentz}] + 
   AJc[{2, Lorentz}, {0, Lorentz}, {2, Lorentz}], sac[-1, 5] == (-I)*AJc[{3, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{3, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sac[-1, 6] == (-I)*AJc[{1, Lorentz}, {0, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {0, Lorentz}, {3, Lorentz}], sac[-1, 7] == (-I)*AJc[{0, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sac[-1, 8] == (-I)*AJc[{3, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{3, Lorentz}, {0, Lorentz}, {2, Lorentz}], sac[-1, 9] == (-I)*AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJc[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}], 
 sac[-1, 10] == (-I)*AJc[{0, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{0, Lorentz}, {0, Lorentz}, {2, Lorentz}], sac[1, 11] == I*AJc[{3, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{3, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sac[1, 12] == I*AJc[{1, Lorentz}, {0, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {0, Lorentz}, {3, Lorentz}], sac[1, 13] == I*AJc[{0, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sac[1, 14] == I*AJc[{3, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{3, Lorentz}, {0, Lorentz}, {2, Lorentz}], sac[1, 15] == I*AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJc[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}], 
 sac[1, 16] == I*AJc[{0, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{0, Lorentz}, {0, Lorentz}, {2, Lorentz}], sac[0, 17] == AJc[{3, Lorentz}, {0, Lorentz}, {3, Lorentz}], 
 sac[0, 18] == AJc[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], sac[0, 19] == -AJc[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {1, Lorentz}, {3, Lorentz}], 
 sac[0, 20] == AJc[{0, Lorentz}, {0, Lorentz}, {3, Lorentz}], sac[0, 21] == AJc[{3, Lorentz}, {1, Lorentz}, {2, Lorentz}], sac[0, 22] == AJc[{1, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sac[0, 23] == AJc[{1, Lorentz}, {0, Lorentz}, {2, Lorentz}] - AJc[{2, Lorentz}, {0, Lorentz}, {1, Lorentz}], sac[0, 24] == AJc[{0, Lorentz}, {1, Lorentz}, {2, Lorentz}]},
        
        "TotallySymmetric",
        {sa[-3, 1] == (-I)*AJ[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] - 3*AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + (3*I)*AJ[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[3, 2] == I*AJ[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] - 3*AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] - (3*I)*AJ[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[-2, 3] == -AJ[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] + (2*I)*AJ[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sa[-2, 4] == -AJ[{0, Lorentz}, {1, Lorentz}, {1, Lorentz}] + (2*I)*AJ[{0, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[2, 5] == -AJ[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] - (2*I)*AJ[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sa[2, 6] == -AJ[{0, Lorentz}, {1, Lorentz}, {1, Lorentz}] - (2*I)*AJ[{0, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[-1, 7] == I*AJ[{1, Lorentz}, {3, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {3, Lorentz}, {3, Lorentz}], 
 sa[-1, 8] == I*AJ[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + I*AJ[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[-1, 9] == I*AJ[{0, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {3, Lorentz}], sa[-1, 10] == I*AJ[{0, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{0, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sa[1, 11] == (-I)*AJ[{1, Lorentz}, {3, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {3, Lorentz}, {3, Lorentz}], 
 sa[1, 12] == (-I)*AJ[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] - I*AJ[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[1, 13] == (-I)*AJ[{0, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {3, Lorentz}], sa[1, 14] == (-I)*AJ[{0, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{0, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sa[0, 15] == AJ[{3, Lorentz}, {3, Lorentz}, {3, Lorentz}], sa[0, 16] == AJ[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], sa[0, 17] == AJ[{0, Lorentz}, {3, Lorentz}, {3, Lorentz}], 
 sa[0, 18] == AJ[{0, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {2, Lorentz}], sa[0, 19] == AJ[{0, Lorentz}, {0, Lorentz}, {3, Lorentz}], sa[0, 20] == AJ[{0, Lorentz}, {0, Lorentz}, {0, Lorentz}], 
 sac[-3, 1] == I*AJc[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] - 3*AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] - (3*I)*AJc[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[3, 2] == (-I)*AJc[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] - 3*AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + (3*I)*AJc[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[-2, 3] == -AJc[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] - (2*I)*AJc[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sac[-2, 4] == -AJc[{0, Lorentz}, {1, Lorentz}, {1, Lorentz}] - (2*I)*AJc[{0, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[2, 5] == -AJc[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] + (2*I)*AJc[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sac[2, 6] == -AJc[{0, Lorentz}, {1, Lorentz}, {1, Lorentz}] + (2*I)*AJc[{0, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[-1, 7] == (-I)*AJc[{1, Lorentz}, {3, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {3, Lorentz}, {3, Lorentz}], 
 sac[-1, 8] == (-I)*AJc[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] - I*AJc[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[-1, 9] == (-I)*AJc[{0, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {3, Lorentz}], sac[-1, 10] == (-I)*AJc[{0, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{0, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sac[1, 11] == I*AJc[{1, Lorentz}, {3, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {3, Lorentz}, {3, Lorentz}], 
 sac[1, 12] == I*AJc[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + I*AJc[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[1, 13] == I*AJc[{0, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {3, Lorentz}], sac[1, 14] == I*AJc[{0, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{0, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sac[0, 15] == AJc[{3, Lorentz}, {3, Lorentz}, {3, Lorentz}], sac[0, 16] == AJc[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], sac[0, 17] == AJc[{0, Lorentz}, {3, Lorentz}, {3, Lorentz}], 
 sac[0, 18] == AJc[{0, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {2, Lorentz}], sac[0, 19] == AJc[{0, Lorentz}, {0, Lorentz}, {3, Lorentz}], sac[0, 20] == AJc[{0, Lorentz}, {0, Lorentz}, {0, Lorentz}]},
        
        "TotallyAntisymmetric",
       {sa[-1, 1] == I*AJ[{0, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {3, Lorentz}], sa[1, 2] == (-I)*AJ[{0, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sa[0, 3] == AJ[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}], sa[0, 4] == AJ[{0, Lorentz}, {1, Lorentz}, {2, Lorentz}], sac[-1, 1] == (-I)*AJc[{0, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sac[1, 2] == I*AJc[{0, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {3, Lorentz}], sac[0, 3] == AJc[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}], sac[0, 4] == AJc[{0, Lorentz}, {1, Lorentz}, {2, Lorentz}]},
        
        _,
        Print["Unknown symmetry type: ", sym, ". Using Generic."];
        {sa[-3, 1] == (-I)*AJ[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] - AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] - AJ[{1, Lorentz}, {2, Lorentz}, {1, Lorentz}] + I*AJ[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] - 
   AJ[{2, Lorentz}, {1, Lorentz}, {1, Lorentz}] + I*AJ[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}] + I*AJ[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[3, 2] == I*AJ[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] - AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] - AJ[{1, Lorentz}, {2, Lorentz}, {1, Lorentz}] - I*AJ[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] - 
   AJ[{2, Lorentz}, {1, Lorentz}, {1, Lorentz}] - I*AJ[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}] - I*AJ[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[-2, 3] == -AJ[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] + I*AJ[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] + I*AJ[{2, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sa[-2, 4] == -AJ[{3, Lorentz}, {1, Lorentz}, {1, Lorentz}] + I*AJ[{3, Lorentz}, {1, Lorentz}, {2, Lorentz}] + I*AJ[{3, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJ[{3, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[-2, 5] == -AJ[{1, Lorentz}, {3, Lorentz}, {1, Lorentz}] + I*AJ[{1, Lorentz}, {3, Lorentz}, {2, Lorentz}] + I*AJ[{2, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {3, Lorentz}, {2, Lorentz}], 
 sa[-2, 6] == -AJ[{1, Lorentz}, {0, Lorentz}, {1, Lorentz}] + I*AJ[{1, Lorentz}, {0, Lorentz}, {2, Lorentz}] + I*AJ[{2, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sa[-2, 7] == -AJ[{0, Lorentz}, {1, Lorentz}, {1, Lorentz}] + I*AJ[{0, Lorentz}, {1, Lorentz}, {2, Lorentz}] + I*AJ[{0, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[-2, 8] == -AJ[{1, Lorentz}, {1, Lorentz}, {0, Lorentz}] + I*AJ[{1, Lorentz}, {2, Lorentz}, {0, Lorentz}] + I*AJ[{2, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {0, Lorentz}], 
 sa[2, 9] == -AJ[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] - I*AJ[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] - I*AJ[{2, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sa[2, 10] == -AJ[{3, Lorentz}, {1, Lorentz}, {1, Lorentz}] - I*AJ[{3, Lorentz}, {1, Lorentz}, {2, Lorentz}] - I*AJ[{3, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJ[{3, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[2, 11] == -AJ[{1, Lorentz}, {3, Lorentz}, {1, Lorentz}] - I*AJ[{1, Lorentz}, {3, Lorentz}, {2, Lorentz}] - I*AJ[{2, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {3, Lorentz}, {2, Lorentz}], 
 sa[2, 12] == -AJ[{1, Lorentz}, {0, Lorentz}, {1, Lorentz}] - I*AJ[{1, Lorentz}, {0, Lorentz}, {2, Lorentz}] - I*AJ[{2, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sa[2, 13] == -AJ[{0, Lorentz}, {1, Lorentz}, {1, Lorentz}] - I*AJ[{0, Lorentz}, {1, Lorentz}, {2, Lorentz}] - I*AJ[{0, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[2, 14] == -AJ[{1, Lorentz}, {1, Lorentz}, {0, Lorentz}] - I*AJ[{1, Lorentz}, {2, Lorentz}, {0, Lorentz}] - I*AJ[{2, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {0, Lorentz}], 
 sa[-1, 15] == I*AJ[{3, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{3, Lorentz}, {2, Lorentz}, {3, Lorentz}], sa[-1, 16] == I*AJ[{1, Lorentz}, {3, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {3, Lorentz}, {3, Lorentz}], 
 sa[-1, 17] == I*AJ[{1, Lorentz}, {0, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {0, Lorentz}, {3, Lorentz}], sa[-1, 18] == I*AJ[{0, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sa[-1, 19] == I*AJ[{3, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJ[{3, Lorentz}, {3, Lorentz}, {2, Lorentz}], sa[-1, 20] == I*AJ[{3, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{3, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sa[-1, 21] == I*AJ[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + I*AJ[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[-1, 22] == I*AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] - I*AJ[{1, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}] - AJ[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}], 
 sa[-1, 23] == I*AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJ[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] - I*AJ[{2, Lorentz}, {1, Lorentz}, {1, Lorentz}] - AJ[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}], 
 sa[-1, 24] == I*AJ[{0, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJ[{0, Lorentz}, {3, Lorentz}, {2, Lorentz}], sa[-1, 25] == I*AJ[{0, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{0, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sa[-1, 26] == I*AJ[{3, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJ[{3, Lorentz}, {2, Lorentz}, {0, Lorentz}], sa[-1, 27] == I*AJ[{1, Lorentz}, {3, Lorentz}, {0, Lorentz}] + AJ[{2, Lorentz}, {3, Lorentz}, {0, Lorentz}], 
 sa[-1, 28] == I*AJ[{1, Lorentz}, {0, Lorentz}, {0, Lorentz}] + AJ[{2, Lorentz}, {0, Lorentz}, {0, Lorentz}], sa[-1, 29] == I*AJ[{0, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {0, Lorentz}], 
 sa[1, 30] == (-I)*AJ[{3, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{3, Lorentz}, {2, Lorentz}, {3, Lorentz}], sa[1, 31] == (-I)*AJ[{1, Lorentz}, {3, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {3, Lorentz}, {3, Lorentz}], 
 sa[1, 32] == (-I)*AJ[{1, Lorentz}, {0, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {0, Lorentz}, {3, Lorentz}], sa[1, 33] == (-I)*AJ[{0, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sa[1, 34] == (-I)*AJ[{3, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJ[{3, Lorentz}, {3, Lorentz}, {2, Lorentz}], sa[1, 35] == (-I)*AJ[{3, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{3, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sa[1, 36] == (-I)*AJ[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] - I*AJ[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sa[1, 37] == (-I)*AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + I*AJ[{1, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}] - AJ[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}], 
 sa[1, 38] == (-I)*AJ[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJ[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] + I*AJ[{2, Lorentz}, {1, Lorentz}, {1, Lorentz}] - AJ[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}], 
 sa[1, 39] == (-I)*AJ[{0, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJ[{0, Lorentz}, {3, Lorentz}, {2, Lorentz}], sa[1, 40] == (-I)*AJ[{0, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{0, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sa[1, 41] == (-I)*AJ[{3, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJ[{3, Lorentz}, {2, Lorentz}, {0, Lorentz}], sa[1, 42] == (-I)*AJ[{1, Lorentz}, {3, Lorentz}, {0, Lorentz}] + AJ[{2, Lorentz}, {3, Lorentz}, {0, Lorentz}], 
 sa[1, 43] == (-I)*AJ[{1, Lorentz}, {0, Lorentz}, {0, Lorentz}] + AJ[{2, Lorentz}, {0, Lorentz}, {0, Lorentz}], sa[1, 44] == (-I)*AJ[{0, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {0, Lorentz}], 
 sa[0, 45] == AJ[{3, Lorentz}, {3, Lorentz}, {3, Lorentz}], sa[0, 46] == AJ[{3, Lorentz}, {0, Lorentz}, {3, Lorentz}], sa[0, 47] == AJ[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sa[0, 48] == -AJ[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] + AJ[{2, Lorentz}, {1, Lorentz}, {3, Lorentz}], sa[0, 49] == AJ[{0, Lorentz}, {3, Lorentz}, {3, Lorentz}], sa[0, 50] == AJ[{0, Lorentz}, {0, Lorentz}, {3, Lorentz}], 
 sa[0, 51] == AJ[{3, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJ[{3, Lorentz}, {2, Lorentz}, {2, Lorentz}], sa[0, 52] == AJ[{3, Lorentz}, {1, Lorentz}, {2, Lorentz}] - AJ[{3, Lorentz}, {2, Lorentz}, {1, Lorentz}], 
 sa[0, 53] == AJ[{1, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {3, Lorentz}, {2, Lorentz}], sa[0, 54] == AJ[{1, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJ[{2, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sa[0, 55] == AJ[{1, Lorentz}, {3, Lorentz}, {2, Lorentz}] - AJ[{2, Lorentz}, {3, Lorentz}, {1, Lorentz}], sa[0, 56] == AJ[{1, Lorentz}, {0, Lorentz}, {2, Lorentz}] - AJ[{2, Lorentz}, {0, Lorentz}, {1, Lorentz}], 
 sa[0, 57] == AJ[{0, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJ[{0, Lorentz}, {2, Lorentz}, {2, Lorentz}], sa[0, 58] == AJ[{0, Lorentz}, {1, Lorentz}, {2, Lorentz}] - AJ[{0, Lorentz}, {2, Lorentz}, {1, Lorentz}], 
 sa[0, 59] == AJ[{3, Lorentz}, {3, Lorentz}, {0, Lorentz}], sa[0, 60] == AJ[{3, Lorentz}, {0, Lorentz}, {0, Lorentz}], sa[0, 61] == AJ[{1, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJ[{2, Lorentz}, {2, Lorentz}, {0, Lorentz}], 
 sa[0, 62] == -AJ[{1, Lorentz}, {2, Lorentz}, {0, Lorentz}] + AJ[{2, Lorentz}, {1, Lorentz}, {0, Lorentz}], sa[0, 63] == AJ[{0, Lorentz}, {3, Lorentz}, {0, Lorentz}], sa[0, 64] == AJ[{0, Lorentz}, {0, Lorentz}, {0, Lorentz}], 
 sac[-3, 1] == I*AJc[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] - AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] - AJc[{1, Lorentz}, {2, Lorentz}, {1, Lorentz}] - I*AJc[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] - 
   AJc[{2, Lorentz}, {1, Lorentz}, {1, Lorentz}] - I*AJc[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}] - I*AJc[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[3, 2] == (-I)*AJc[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] - AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] - AJc[{1, Lorentz}, {2, Lorentz}, {1, Lorentz}] + I*AJc[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] - 
   AJc[{2, Lorentz}, {1, Lorentz}, {1, Lorentz}] + I*AJc[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}] + I*AJc[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[-2, 3] == -AJc[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] - I*AJc[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] - I*AJc[{2, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sac[-2, 4] == -AJc[{3, Lorentz}, {1, Lorentz}, {1, Lorentz}] - I*AJc[{3, Lorentz}, {1, Lorentz}, {2, Lorentz}] - I*AJc[{3, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJc[{3, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[-2, 5] == -AJc[{1, Lorentz}, {3, Lorentz}, {1, Lorentz}] - I*AJc[{1, Lorentz}, {3, Lorentz}, {2, Lorentz}] - I*AJc[{2, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {3, Lorentz}, {2, Lorentz}], 
 sac[-2, 6] == -AJc[{1, Lorentz}, {0, Lorentz}, {1, Lorentz}] - I*AJc[{1, Lorentz}, {0, Lorentz}, {2, Lorentz}] - I*AJc[{2, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sac[-2, 7] == -AJc[{0, Lorentz}, {1, Lorentz}, {1, Lorentz}] - I*AJc[{0, Lorentz}, {1, Lorentz}, {2, Lorentz}] - I*AJc[{0, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[-2, 8] == -AJc[{1, Lorentz}, {1, Lorentz}, {0, Lorentz}] - I*AJc[{1, Lorentz}, {2, Lorentz}, {0, Lorentz}] - I*AJc[{2, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {0, Lorentz}], 
 sac[2, 9] == -AJc[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] + I*AJc[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] + I*AJc[{2, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sac[2, 10] == -AJc[{3, Lorentz}, {1, Lorentz}, {1, Lorentz}] + I*AJc[{3, Lorentz}, {1, Lorentz}, {2, Lorentz}] + I*AJc[{3, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJc[{3, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[2, 11] == -AJc[{1, Lorentz}, {3, Lorentz}, {1, Lorentz}] + I*AJc[{1, Lorentz}, {3, Lorentz}, {2, Lorentz}] + I*AJc[{2, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {3, Lorentz}, {2, Lorentz}], 
 sac[2, 12] == -AJc[{1, Lorentz}, {0, Lorentz}, {1, Lorentz}] + I*AJc[{1, Lorentz}, {0, Lorentz}, {2, Lorentz}] + I*AJc[{2, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sac[2, 13] == -AJc[{0, Lorentz}, {1, Lorentz}, {1, Lorentz}] + I*AJc[{0, Lorentz}, {1, Lorentz}, {2, Lorentz}] + I*AJc[{0, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[2, 14] == -AJc[{1, Lorentz}, {1, Lorentz}, {0, Lorentz}] + I*AJc[{1, Lorentz}, {2, Lorentz}, {0, Lorentz}] + I*AJc[{2, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {0, Lorentz}], 
 sac[-1, 15] == (-I)*AJc[{3, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{3, Lorentz}, {2, Lorentz}, {3, Lorentz}], sac[-1, 16] == (-I)*AJc[{1, Lorentz}, {3, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {3, Lorentz}, {3, Lorentz}], 
 sac[-1, 17] == (-I)*AJc[{1, Lorentz}, {0, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {0, Lorentz}, {3, Lorentz}], sac[-1, 18] == (-I)*AJc[{0, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sac[-1, 19] == (-I)*AJc[{3, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJc[{3, Lorentz}, {3, Lorentz}, {2, Lorentz}], sac[-1, 20] == (-I)*AJc[{3, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{3, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sac[-1, 21] == (-I)*AJc[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] - I*AJc[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[-1, 22] == (-I)*AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + I*AJc[{1, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}] - AJc[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}], 
 sac[-1, 23] == (-I)*AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJc[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] + I*AJc[{2, Lorentz}, {1, Lorentz}, {1, Lorentz}] - AJc[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}], 
 sac[-1, 24] == (-I)*AJc[{0, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJc[{0, Lorentz}, {3, Lorentz}, {2, Lorentz}], sac[-1, 25] == (-I)*AJc[{0, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{0, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sac[-1, 26] == (-I)*AJc[{3, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJc[{3, Lorentz}, {2, Lorentz}, {0, Lorentz}], sac[-1, 27] == (-I)*AJc[{1, Lorentz}, {3, Lorentz}, {0, Lorentz}] + AJc[{2, Lorentz}, {3, Lorentz}, {0, Lorentz}], 
 sac[-1, 28] == (-I)*AJc[{1, Lorentz}, {0, Lorentz}, {0, Lorentz}] + AJc[{2, Lorentz}, {0, Lorentz}, {0, Lorentz}], sac[-1, 29] == (-I)*AJc[{0, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {0, Lorentz}], 
 sac[1, 30] == I*AJc[{3, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{3, Lorentz}, {2, Lorentz}, {3, Lorentz}], sac[1, 31] == I*AJc[{1, Lorentz}, {3, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {3, Lorentz}, {3, Lorentz}], 
 sac[1, 32] == I*AJc[{1, Lorentz}, {0, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {0, Lorentz}, {3, Lorentz}], sac[1, 33] == I*AJc[{0, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sac[1, 34] == I*AJc[{3, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJc[{3, Lorentz}, {3, Lorentz}, {2, Lorentz}], sac[1, 35] == I*AJc[{3, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{3, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sac[1, 36] == I*AJc[{1, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + I*AJc[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {2, Lorentz}], 
 sac[1, 37] == I*AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] - I*AJc[{1, Lorentz}, {2, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {1, Lorentz}, {2, Lorentz}] - AJc[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}], 
 sac[1, 38] == I*AJc[{1, Lorentz}, {1, Lorentz}, {2, Lorentz}] + AJc[{1, Lorentz}, {2, Lorentz}, {2, Lorentz}] - I*AJc[{2, Lorentz}, {1, Lorentz}, {1, Lorentz}] - AJc[{2, Lorentz}, {2, Lorentz}, {1, Lorentz}], 
 sac[1, 39] == I*AJc[{0, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJc[{0, Lorentz}, {3, Lorentz}, {2, Lorentz}], sac[1, 40] == I*AJc[{0, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{0, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sac[1, 41] == I*AJc[{3, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJc[{3, Lorentz}, {2, Lorentz}, {0, Lorentz}], sac[1, 42] == I*AJc[{1, Lorentz}, {3, Lorentz}, {0, Lorentz}] + AJc[{2, Lorentz}, {3, Lorentz}, {0, Lorentz}], 
 sac[1, 43] == I*AJc[{1, Lorentz}, {0, Lorentz}, {0, Lorentz}] + AJc[{2, Lorentz}, {0, Lorentz}, {0, Lorentz}], sac[1, 44] == I*AJc[{0, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {0, Lorentz}], 
 sac[0, 45] == AJc[{3, Lorentz}, {3, Lorentz}, {3, Lorentz}], sac[0, 46] == AJc[{3, Lorentz}, {0, Lorentz}, {3, Lorentz}], sac[0, 47] == AJc[{1, Lorentz}, {1, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {3, Lorentz}], 
 sac[0, 48] == -AJc[{1, Lorentz}, {2, Lorentz}, {3, Lorentz}] + AJc[{2, Lorentz}, {1, Lorentz}, {3, Lorentz}], sac[0, 49] == AJc[{0, Lorentz}, {3, Lorentz}, {3, Lorentz}], sac[0, 50] == AJc[{0, Lorentz}, {0, Lorentz}, {3, Lorentz}], 
 sac[0, 51] == AJc[{3, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJc[{3, Lorentz}, {2, Lorentz}, {2, Lorentz}], sac[0, 52] == AJc[{3, Lorentz}, {1, Lorentz}, {2, Lorentz}] - AJc[{3, Lorentz}, {2, Lorentz}, {1, Lorentz}], 
 sac[0, 53] == AJc[{1, Lorentz}, {3, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {3, Lorentz}, {2, Lorentz}], sac[0, 54] == AJc[{1, Lorentz}, {0, Lorentz}, {1, Lorentz}] + AJc[{2, Lorentz}, {0, Lorentz}, {2, Lorentz}], 
 sac[0, 55] == AJc[{1, Lorentz}, {3, Lorentz}, {2, Lorentz}] - AJc[{2, Lorentz}, {3, Lorentz}, {1, Lorentz}], sac[0, 56] == AJc[{1, Lorentz}, {0, Lorentz}, {2, Lorentz}] - AJc[{2, Lorentz}, {0, Lorentz}, {1, Lorentz}], 
 sac[0, 57] == AJc[{0, Lorentz}, {1, Lorentz}, {1, Lorentz}] + AJc[{0, Lorentz}, {2, Lorentz}, {2, Lorentz}], sac[0, 58] == AJc[{0, Lorentz}, {1, Lorentz}, {2, Lorentz}] - AJc[{0, Lorentz}, {2, Lorentz}, {1, Lorentz}], 
 sac[0, 59] == AJc[{3, Lorentz}, {3, Lorentz}, {0, Lorentz}], sac[0, 60] == AJc[{3, Lorentz}, {0, Lorentz}, {0, Lorentz}], sac[0, 61] == AJc[{1, Lorentz}, {1, Lorentz}, {0, Lorentz}] + AJc[{2, Lorentz}, {2, Lorentz}, {0, Lorentz}], 
 sac[0, 62] == -AJc[{1, Lorentz}, {2, Lorentz}, {0, Lorentz}] + AJc[{2, Lorentz}, {1, Lorentz}, {0, Lorentz}], sac[0, 63] == AJc[{0, Lorentz}, {3, Lorentz}, {0, Lorentz}], sac[0, 64] == AJc[{0, Lorentz}, {0, Lorentz}, {0, Lorentz}]}
    ];
    
helicityTemplate["scalar", _] := {ss[0, 1] == SJ, ssc[0, 1] == SJc};

(* B_{i3} and B_{0i} (i = 1, 2) rotate like a vector index about the momentum axis; B_12 and B_03 are invariant *)
helicityTemplate["rank2anti", _] := {
    sb[-1, 1] == I*protoBJ[{1, Lorentz}, {3, Lorentz}] + protoBJ[{2, Lorentz}, {3, Lorentz}],
    sb[-1, 2] == I*protoBJ[{0, Lorentz}, {1, Lorentz}] + protoBJ[{0, Lorentz}, {2, Lorentz}],
    sb[1, 3] == (-I)*protoBJ[{1, Lorentz}, {3, Lorentz}] + protoBJ[{2, Lorentz}, {3, Lorentz}],
    sb[1, 4] == (-I)*protoBJ[{0, Lorentz}, {1, Lorentz}] + protoBJ[{0, Lorentz}, {2, Lorentz}],
    sb[0, 5] == protoBJ[{1, Lorentz}, {2, Lorentz}],
    sb[0, 6] == protoBJ[{0, Lorentz}, {3, Lorentz}],
    sbc[-1, 1] == (-I)*protoBJc[{1, Lorentz}, {3, Lorentz}] + protoBJc[{2, Lorentz}, {3, Lorentz}],
    sbc[-1, 2] == (-I)*protoBJc[{0, Lorentz}, {1, Lorentz}] + protoBJc[{0, Lorentz}, {2, Lorentz}],
    sbc[1, 3] == I*protoBJc[{1, Lorentz}, {3, Lorentz}] + protoBJc[{2, Lorentz}, {3, Lorentz}],
    sbc[1, 4] == I*protoBJc[{0, Lorentz}, {1, Lorentz}] + protoBJc[{0, Lorentz}, {2, Lorentz}],
    sbc[0, 5] == protoBJc[{1, Lorentz}, {2, Lorentz}],
    sbc[0, 6] == protoBJc[{0, Lorentz}, {3, Lorentz}]
};

(* type -> {prototype source, prototype conjugate source, built-in field keeping the untagged names} *)
$helicityPrototype = <|
    "vector" -> {VJ, VJc, V}, "rank2sym" -> {HJ, HJc, H}, "rank3" -> {AJ, AJc, A},
    "scalar" -> {SJ, SJc, S}, "rank2anti" -> {protoBJ, protoBJc, None}
|>;

(* Helicity definitions for one field: its own sources; copies get the field name as a third argument *)
helicityRules[fld_] := Module[{type, proto, rules},
    type = $FieldInstances[fld]["type"];
    proto = $helicityPrototype[type];
    rules = helicityTemplate[type, $FieldInstances[fld]["symmetry"]] /.
        {proto[[1]] -> fieldJ[fld], proto[[2]] -> fieldJc[fld]};
    If[fld === proto[[3]],
        rules,
        rules /. (hd : (sv | svc | sh | shc | sa | sac | ss | ssc | sb | sbc))[h_, i_] :> hd[h, i, fld]
    ]
];

TensorToHelicity[] := Module[{},

helicityGenericRules = Flatten[helicityRules /@ activeFields];
(*These undergo a further reduction due to source constraints*)

Reduce[FullSimplify[helicityGenericRules/.masslessSolution/.masslessSolutionConj/.\[Kappa]->\[Omega]],Join[allSourceVars,allSourceVarsConj]]//ToRules

]




End[] (* `Private` *)
EndPackage[]
