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
TensorToHelicity::usage = "...";


sv::usage = "Helicity amplitude.";
svc::usage = "Conjugate helicity amplitude.";
sh::usage = "Helicity amplitude.";
shc::usage = "Conjugate helicity amplitude.";
sa::usage = "Helicity amplitude.";
sac::usage = "Conjugate helicity amplitude.";



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

XX1=XX2=XX3=XX4=0; (*These are control-variables used to check the mixed projectors operators. I used them in my original computation. Never appear in physical objects. Choosing 0.*)


If[!ValueQ[Global`DivideEtImpera],
  Global`DivideEtImpera = "True";
];

If[!ValueQ[Global`rank3SymmetryType],
  Global`rank3SymmetryType = "TotallySymmetric";
];


PrintK["Defining fields. Rank-3 field A is declared: ", Global`rank3SymmetryType];
PrintK["You can change the value of the string variable 'rank3SymmetryType' before starting the package. You can choose among 'Generic', 'Sym23' , 'Antisym23', 'TotallySymmetric', 'TotallyAntisymmetric'  ."];
PrintK["Use H for rank-2 totally symmetric tensor, V for vector and S for scalar. "];
PrintK["The source of 'field' is 'field name'J. The complex conjugate source is 'field name'Jc."];
PrintK["============================================================================"]

Block[{Print = (Null &)},

(* Setup function for rank-3 tensor based on symmetry type *)
    Switch[Global`rank3SymmetryType,
        
        "Generic",
        DefTensor[A[-a, -b, -c], M];
        DefTensorPerturbation[\[Epsilon]Ap[LI[1], -a, -b, -c], A[-a, -b, -c], M];
        Print["Rank-3 tensor A: Generic (no symmetry)"];,
        
        "Sym23",
        DefTensor[A[-a, -b, -c], M, Symmetric[{-b, -c}]];
        DefTensorPerturbation[\[Epsilon]Ap[LI[1], -a, -b, -c], A[-a, -b, -c], M, Symmetric[{-b, -c}]];
        Print["Rank-3 tensor A: Symmetric in last 2 indices"];,
        
        "Antisym23",
        DefTensor[A[-a, -b, -c], M, Antisymmetric[{-b, -c}]];
        DefTensorPerturbation[\[Epsilon]Ap[LI[1], -a, -b, -c], A[-a, -b, -c], M, Antisymmetric[{-b, -c}]];
        Print["Rank-3 tensor A: Antisymmetric in last 2 indices"];,
        
        "TotallySymmetric",
        DefTensor[A[-a, -b, -c], M, Symmetric[{-a, -b, -c}]];
        DefTensorPerturbation[\[Epsilon]Ap[LI[1], -a, -b, -c], A[-a, -b, -c], M, Symmetric[{-a, -b, -c}]];
        Print["Rank-3 tensor A: Totally symmetric"];,
        
        "TotallyAntisymmetric",
        DefTensor[A[-a, -b, -c], M, Antisymmetric[{-a, -b, -c}]];
        DefTensorPerturbation[\[Epsilon]Ap[LI[1], -a, -b, -c], A[-a, -b, -c], M, Antisymmetric[{-a, -b, -c}]];
        Print["Rank-3 tensor A: Totally antisymmetric"];,
        
        _,
        Print["Unknown symmetry type: ", symType, ". Using Generic."];
        DefTensor[A[-a, -b, -c], M];
        DefTensorPerturbation[\[Epsilon]Ap[LI[1], -a, -b, -c], A[-a, -b, -c], M];
    ];
    
    
DefTensor[H[-a, -b], M, Symmetric[{-a, -b}]];
DefTensorPerturbation[\[Epsilon]h[LI[1], -a, -b], H[-a, -b], M, Symmetric[{-a, -b}]];
DefTensor[V[-a], M];
DefTensorPerturbation[\[Epsilon]Vp[LI[1], -a], V[-a], M];

DefTensor[S[], M];
DefTensorPerturbation[\[Epsilon]Sp[LI[1]], S[], M];

];


(*BEGIN SOURCE DEFINITION*)


(* Setup function for rank-3 source tensors based on symmetry type *)

Block[{Print = (Null &)},


    Switch[Global`rank3SymmetryType,
        
        "Generic",
        DefTensor[AJ[a, b, c], M];
        DefTensor[AJc[a, b, c], M];,
        
        "Sym23",
        DefTensor[AJ[a, b, c], M, Symmetric[{b, c}]];
        DefTensor[AJc[a, b, c], M, Symmetric[{b, c}]];,
        
        "Antisym23",
        DefTensor[AJ[a, b, c], M, Antisymmetric[{b, c}]];
        DefTensor[AJc[a, b, c], M, Antisymmetric[{b, c}]];,
        
        "TotallySymmetric",
        DefTensor[AJ[a, b, c], M, Symmetric[{a, b, c}]];
        DefTensor[AJc[a, b, c], M, Symmetric[{a, b, c}]];,
        
        "TotallyAntisymmetric",
        DefTensor[AJ[a, b, c], M, Antisymmetric[{a, b, c}]];
        DefTensor[AJc[a, b, c], M, Antisymmetric[{a, b, c}]];,
        
        _,
        DefTensor[AJ[a, b, c], M];
        DefTensor[AJc[a, b, c], M];
   ];
   
(* Other sources - fixed symmetries *)
DefTensor[HJ[a, b], M, Symmetric[{a, b}]];
DefTensor[HJc[a, b], M, Symmetric[{a, b}]];
DefTensor[VJ[a], M];
DefTensor[VJc[a], M];
DefConstantSymbol[SJ];
DefConstantSymbol[SJc];

 ];



(*END SOURCE DEFINITION*)

Begin["`Private`"]

ConstantSymbols[input_]:= Block[{Print = (Null &)},DefConstantSymbol[input];];

Canonicalize[input_]:= CollectTensors[input]


dd[input_][input2_]:=PD[input][input2];

(*******************Questo era un pacchetto separato, vediamo se funziona*********************************************)
(*********************************************************************************************************************)
(*********************************************************************************************************************)
(*********************************************************************************************************************)
(*********************************************************************************************************************)
(*********************************************************************************************************************)
(*********************************************************************************************************************)
(*********************************************************************************************************************)


(* === Internal data structures === *)

(* Field type registry: name -> <|"rank" -> n, "symmetry" -> sym, "slot" -> slotIndex|> *)
$FieldTypeRegistry = <||>;

(* Representation content: {fieldType, {J, parity}} -> {labels} *)
$RepContentTable = <||>;

(* Current configuration (built by BuildRepresentationTables) *)
$CurrentListRep = {};
$CurrentListC = <||>;
$CurrentListF = <||>;
$CurrentFields = {};

(* Slot mapping: rank -> slot index *)
$SlotFromRank = <|3 -> 1, 2 -> 2, 1 -> 3, 0 -> 4|>;
$SlotSymbol = <|1 -> i3, 2 -> i2, 3 -> i1, 4 -> i0|>;
$NumSlots = 4;  (* Current maximum: rank 0,1,2,3 *)

(* === Field Type Registry === *)

DefineFieldType[name_String, rank_Integer, symmetry_] := Module[{},
    If[rank < 0 || rank > 3,
        Message[DefineFieldType::badrank, rank];
        Return[$Failed]
    ];
    $FieldTypeRegistry[name] = <|
        "rank" -> rank,
        "symmetry" -> symmetry,
        "slot" -> $SlotFromRank[rank]
    |>;
    name
];

DefineFieldType::badrank = "Rank `1` is not supported. Currently supported: 0, 1, 2, 3.";

GetFieldType[name_String] := Lookup[$FieldTypeRegistry, name, Missing["NotFound", name]];

ListFieldTypes[] := Keys[$FieldTypeRegistry];

(* === Representation Content === *)

SetRepContent[fieldType_String, {j_Integer, parity_}, labels_List] := Module[{},
    If[!KeyExistsQ[$FieldTypeRegistry, fieldType],
        Message[SetRepContent::unknown, fieldType];
        Return[$Failed]
    ];
    $RepContentTable[{fieldType, {j, parity}}] = labels;
];

SetRepContent::unknown = "Field type `1` is not registered. Use DefineFieldType first.";

GetRepContent[fieldType_String] := Module[{keys, result},
    keys = Select[Keys[$RepContentTable], #[[1]] === fieldType &];
    result = <||>;
    Do[
        result[key[[2]]] = $RepContentTable[key],
        {key, keys}
    ];
    result
];

(* === Build Representation Tables === *)

BuildRepresentationTables[fieldList_List] := Module[
    {allJP, fieldSlots, jp, slot, labels, combined},
    
    (* Validate all fields exist *)
    Do[
        If[!KeyExistsQ[$FieldTypeRegistry, f],
            Message[BuildRepresentationTables::unknown, f];
            Return[$Failed]
        ],
        {f, fieldList}
    ];
    
    (* Store current field configuration *)
    $CurrentFields = fieldList;
    
    (* Get all J^P sectors that appear in any of the specified fields *)
    allJP = DeleteDuplicates @ Flatten[
        Table[
            Keys[GetRepContent[f]],
            {f, fieldList}
        ],
        1
    ];
    
    (* Sort: by spin descending, then by parity (m before p for same spin) *)
    allJP = SortBy[allJP, {-#[[1]] &, #[[2]] === p &}];
    
    $CurrentListRep = allJP;
    
    (* Build listC and listF for each J^P sector *)
    $CurrentListC = <||>;
    $CurrentListF = <||>;
    
    Do[
        (* Initialize empty slots *)
        combined = Table[{}, {$NumSlots}];
        
        (* Fill in labels from each field that contributes to this J^P *)
        Do[
            slot = $FieldTypeRegistry[f]["slot"];
            labels = Lookup[$RepContentTable, Key[{f, jp}], {}];
            combined[[slot]] = Join[combined[[slot]], labels],
            {f, fieldList}
        ];
        
        (* Store listC *)
        $CurrentListC[jp] = combined;
        
        (* Store listF for each slot *)
        Do[
            $CurrentListF[{jp, $SlotSymbol[s]}] = combined[[s]],
            {s, 1, $NumSlots}
        ],
        
        {jp, allJP}
    ];
    
    (* Return summary *)
    <|
        "fields" -> fieldList,
        "sectors" -> allJP,
        "listC" -> $CurrentListC
    |>
];

BuildRepresentationTables::unknown = "Field type `1` is not registered.";

(* === Accessors === *)

GetListRep[] := $CurrentListRep;

GetListC[j_Integer, parity_] := Lookup[$CurrentListC, Key[{j, parity}], {}];
GetListC[{j_Integer, parity_}] := GetListC[j, parity];

GetListF[j_Integer, parity_, slot_] := Lookup[$CurrentListF, Key[{{j, parity}, slot}], {}];
GetListF[{j_Integer, parity_}, slot_] := GetListF[j, parity, slot];

(* === Utilities === *)

ClearConfiguration[] := Module[{},
    $CurrentListRep = {};
    $CurrentListC = <||>;
    $CurrentListF = <||>;
    $CurrentFields = {};
];

PrintConfiguration[] := Module[{jp},
    PrintK["=== Current Field Configuration ==="];
    PrintK["Fields: ", $CurrentFields];
    PrintK[""];
    PrintK["J^P Sectors: ", $CurrentListRep];
    PrintK[""];
    PrintK["Representation Tables:"];
    Do[
        PrintK["  sector[", jp[[1]], ",", jp[[2]], "] = ", $CurrentListC[jp]],
        {jp, $CurrentListRep}
    ];
];

(* === Initialize built-in field types and representation content === *)

(* Register the standard field types *)
DefineFieldType["rank3", 3, Generic];
DefineFieldType["rank2sym", 2, Symmetric];
DefineFieldType["rank2anti", 2, Antisymmetric];
DefineFieldType["vector", 1, None];
DefineFieldType["scalar", 0, None];

(* Set representation content for rank-3 generic tensor *)
SetRepContent["rank3", {3, m}, {1}];
SetRepContent["rank3", {2, p}, {1, 2, 3}];
SetRepContent["rank3", {2, m}, {1, 2}];
SetRepContent["rank3", {1, p}, {1, 2, 3}];
SetRepContent["rank3", {1, m}, {1, 2, 3, 4, 5, 6}];
SetRepContent["rank3", {0, p}, {1, 2, 3, 4}];
SetRepContent["rank3", {0, m}, {1}];

(* Set representation content for rank-2 symmetric tensor *)
SetRepContent["rank2sym", {2, p}, {4}];
SetRepContent["rank2sym", {1, m}, {7}];
SetRepContent["rank2sym", {0, p}, {5, 6}];

(* Set representation content for vector *)
SetRepContent["vector", {1, m}, {8}];
SetRepContent["vector", {0, p}, {7}];

(* Set representation content for scalar *)
SetRepContent["scalar", {0, p}, {8}];

(* Note: rank2anti content to be added when projectors are available *)


(*********************************************************************************************************************)
(*********************************************************************************************************************)
(*********************************************************************************************************************)
(*********************************************************************************************************************)
(*********************************************************************************************************************)
(*********************************************************************************************************************)
(*********************************************************************************************************************)
(*********************************************************************************************************************)
(*********************************************************************************************************************)
(*********************************************************************************************************************)




(* Placeholder - will be populated later *)
(* Global variable to track active fields in the theory *)
(* Set by BuildRepresentationTables or manually *)
activeFieldTypes = {};

initializeSourceComponentRules[] := initializeSourceComponentRules[activeFieldTypes];

initializeSourceComponentRules[fieldTypes_List] := Module[{initialized = {}},
    
    (* Initialize empty rules *)
    VJrules = {}; VJcrules = {};
    HJrules = {}; HJcrules = {};
    tabAJrules = {}; tabAJcrules = {};
    
    (* Vector source rules - only if vector field is active *)
    If[MemberQ[fieldTypes, "vector"],
        VJrules = ChangeComponents[VJ[-{a, Lorentz}], VJ[{a, Lorentz}]];
        VJcrules = ChangeComponents[VJc[-{a, Lorentz}], VJc[{a, Lorentz}]];
        AppendTo[initialized, "vector (VJ, VJc)"];
    ];
    
    (* Rank-2 source rules - only if rank2sym field is active *)
    If[MemberQ[fieldTypes, "rank2sym"],
        HJrules = ChangeComponents[HJ[{-a, -Lorentz}, {-b, -Lorentz}], 
            HJ[{a, Lorentz}, {b, Lorentz}]] // 
            ReplaceRepeated[#, MomentumRule1] & // 
            ReplaceRepeated[#, MetricRule1] & // 
            ReplaceRepeated[#, MetricRule2] & // FullSimplify;
        
        HJcrules = ChangeComponents[HJc[{-a, -Lorentz}, {-b, -Lorentz}], 
            HJc[{a, Lorentz}, {b, Lorentz}]] // 
            ReplaceRepeated[#, MomentumRule1] & // 
            ReplaceRepeated[#, MetricRule1] & // 
            ReplaceRepeated[#, MetricRule2] & // FullSimplify;
        AppendTo[initialized, "rank2sym (HJ, HJc)"];
    ];
    
    (* Rank-3 source rules - only if rank3 field is active *)
    If[MemberQ[fieldTypes, "rank3"],
        tabAJrules = Flatten[Table[
            ChangeComponents[
                AJ[x1*{a, Lorentz}, x2*{b, Lorentz}, x3*{c, Lorentz}], 
                AJ[{a, Lorentz}, {b, Lorentz}, {c, Lorentz}]] // 
            ReplaceRepeated[#, MomentumRule1] & // 
            ReplaceRepeated[#, MetricRule1] & // 
            ReplaceRepeated[#, MetricRule2] & // FullSimplify,
            {x1, {-1, +1}}, {x2, {-1, +1}}, {x3, {-1, +1}}
        ]];
        
        tabAJcrules = Flatten[Table[
            ChangeComponents[
                AJc[x1*{a, Lorentz}, x2*{b, Lorentz}, x3*{c, Lorentz}], 
                AJc[{a, Lorentz}, {b, Lorentz}, {c, Lorentz}]] // 
            ReplaceRepeated[#, MomentumRule1] & // 
            ReplaceRepeated[#, MetricRule1] & // 
            ReplaceRepeated[#, MetricRule2] & // FullSimplify,
            {x1, {-1, +1}}, {x2, {-1, +1}}, {x3, {-1, +1}}
        ]];
        AppendTo[initialized, "rank3 (AJ, AJc)"];
    ];
    
    (* Scalar sources are just constants - no component rules needed *)
    If[MemberQ[fieldTypes, "scalar"],
        AppendTo[initialized, "scalar (SJ, SJc - constants)"];
    ];
    
    If[Length[initialized] > 0,
        PrintK["Source component rules initialized for: ", StringRiffle[initialized, ", "]];,
        PrintK["No source component rules initialized (no active fields)."];
    ];
];


(* ToComponents: applies index conversion rules (lower -> upper) *)
ToComponents[expr_] := Module[{result},
    result = expr;
    
    (* Apply vector source rules *)
    result = result // ReplaceRepeated[#, VJrules] &;
    result = result // ReplaceRepeated[#, VJcrules] &;
    
    (* Apply rank-2 source rules *)
    result = result // ReplaceRepeated[#, HJrules] &;
    result = result // ReplaceRepeated[#, HJcrules] &;
    
    (* Apply rank-3 source rules *)
    Do[result = result // ReplaceRepeated[#, tabAJrules[[k]]] &, {k, Length[tabAJrules]}];
    Do[result = result // ReplaceRepeated[#, tabAJcrules[[k]]] &, {k, Length[tabAJcrules]}];
    
    (* Apply momentum and metric rules *)
    result = result // ReplaceRepeated[#, MomentumRule1] &;
    result = result // ReplaceRepeated[#, MetricRule1] &;
    result = result // ReplaceRepeated[#, MetricRule2] &;
    
    result
];

TensorToComponents[expr_] := expr // DummyToBasis[Lorentz] // TraceBasisDummy // ToComponents;


(* Check if a field appears in an expression *)
containsField[expr_, field_] := Not[FreeQ[expr, field]];

(* Detect fields in action *)
detectFields[act_] := Module[{fields = {}},
    If[containsField[act, A], AppendTo[fields, "rank3"]];
    If[containsField[act, H], AppendTo[fields, "rank2sym"]];
    If[containsField[act, V], AppendTo[fields, "vector"]];
    If[containsField[act, S], AppendTo[fields, "scalar"]];
    fields
];

(* Define the substitution rules *)
ToPerturbations = {
    A[inds__] :> \[Epsilon]Ap[LI[1], inds],
    H[inds__] :> \[Epsilon]h[LI[1], inds],
    V[inds__] :> \[Epsilon]Vp[LI[1], inds],
    S[] :> \[Epsilon]Sp[LI[1]]
};

VarAction[ac_, x_] := VarD[x, PD][ac]

ToMomentumSpaceRules = {
    PD[s_][delta[x_, y_]] :> (-I * q[s]) * delta[x, y],
    PD[s_][q[x_]] :> 0,
    delta[x__, z__] :> 1
};

(* car[spin, parity, label] returns position of label in listC[spin, parity] *)
car[SPIN_, PARITY_, num_] := Position[listC[SPIN, PARITY], num][[1, 1]];


SourceC[x_, y_, {}] := {};

SourceC[SPIN_, PARITY_, nullInp_] := Module[{spin = SPIN, parity = PARITY, prova, null, labels, elem},
    (* Handle complex conjugate *)
    null = nullInp // Conjugate // ReplaceAll[Conjugate[x_] :> x];
    
    labels = listC[spin, parity];
    If[labels === {}, Return[{}]];
    
    prova = Table[
        Sum[
            elem = Which[
                (* i in rank3, j in rank3 *)
                MemberQ[listF[spin, parity, i3], i] && MemberQ[listF[spin, parity, i3], j],
                null[[car[spin, parity, j]]] * P[spin, parity, {i, j}][-m, -n, -f, -r, -s, -t] * AJ[r, s, t],
                
                (* i in rank2, j in rank3 *)
                MemberQ[listF[spin, parity, i2], i] && MemberQ[listF[spin, parity, i3], j],
                null[[car[spin, parity, j]]] * P[spin, parity, {i, j}][-m, -n, -f, -r, -s] * AJ[f, r, s],
                
                (* i in vector, j in rank3 *)
                MemberQ[listF[spin, parity, i1], i] && MemberQ[listF[spin, parity, i3], j],
                null[[car[spin, parity, j]]] * P[spin, parity, {i, j}][-m, -n, -f, -r] * AJ[n, f, r],
                
                (* i in scalar, j in rank3 *)
                MemberQ[listF[spin, parity, i0], i] && MemberQ[listF[spin, parity, i3], j],
                null[[car[spin, parity, j]]] * P[spin, parity, {i, j}][-m, -n, -f] * AJ[m, n, f],
                
                (* i in rank3, j in rank2 *)
                MemberQ[listF[spin, parity, i3], i] && MemberQ[listF[spin, parity, i2], j],
                null[[car[spin, parity, j]]] * P[spin, parity, {i, j}][-m, -n, -f, -r, -s] * HJ[r, s],
                
                (* i in rank2, j in rank2 *)
                MemberQ[listF[spin, parity, i2], i] && MemberQ[listF[spin, parity, i2], j],
                null[[car[spin, parity, j]]] * P[spin, parity, {i, j}][-m, -n, -f, -r] * HJ[f, r],
                
                (* i in vector, j in rank2 *)
                MemberQ[listF[spin, parity, i1], i] && MemberQ[listF[spin, parity, i2], j],
                null[[car[spin, parity, j]]] * P[spin, parity, {i, j}][-m, -n, -f] * HJ[n, f],
                
                (* i in scalar, j in rank2 *)
                MemberQ[listF[spin, parity, i0], i] && MemberQ[listF[spin, parity, i2], j],
                null[[car[spin, parity, j]]] * P[spin, parity, {i, j}][-m, -n] * HJ[m, n],
                
                (* i in rank3, j in vector *)
                MemberQ[listF[spin, parity, i3], i] && MemberQ[listF[spin, parity, i1], j],
                null[[car[spin, parity, j]]] * P[spin, parity, {i, j}][-m, -n, -f, -r] * VJ[r],
                
                (* i in rank2, j in vector *)
                MemberQ[listF[spin, parity, i2], i] && MemberQ[listF[spin, parity, i1], j],
                null[[car[spin, parity, j]]] * P[spin, parity, {i, j}][-m, -n, -f] * VJ[f],
                
                (* i in vector, j in vector *)
                MemberQ[listF[spin, parity, i1], i] && MemberQ[listF[spin, parity, i1], j],
                null[[car[spin, parity, j]]] * P[spin, parity, {i, j}][-m, -n] * VJ[n],
                
                (* i in scalar, j in vector *)
                MemberQ[listF[spin, parity, i0], i] && MemberQ[listF[spin, parity, i1], j],
                null[[car[spin, parity, j]]] * P[spin, parity, {i, j}][-m] * VJ[m],
                
                (* i in rank3, j in scalar *)
                MemberQ[listF[spin, parity, i3], i] && MemberQ[listF[spin, parity, i0], j],
                null[[car[spin, parity, j]]] * P[spin, parity, {i, j}][-m, -n, -f] * SJ,
                
                (* i in rank2, j in scalar *)
                MemberQ[listF[spin, parity, i2], i] && MemberQ[listF[spin, parity, i0], j],
                null[[car[spin, parity, j]]] * P[spin, parity, {i, j}][-m, -n] * SJ,
                
                (* i in vector, j in scalar *)
                MemberQ[listF[spin, parity, i1], i] && MemberQ[listF[spin, parity, i0], j],
                null[[car[spin, parity, j]]] * P[spin, parity, {i, j}][-m] * SJ,
                
                (* i in scalar, j in scalar *)
                MemberQ[listF[spin, parity, i0], i] && MemberQ[listF[spin, parity, i0], j],
                null[[car[spin, parity, j]]] * P[spin, parity, {i, j}][] * SJ,
                
                (* Default *)
                True, 0
            ];
            elem,
            {j, labels}
        ],
        {i, labels}
    ];
    
    prova = prova // ContractMetric // ReplaceAll[dimM -> 4] // FullSimplify;
    prova
];

GaugeC[SPIN_, PARITY_, null_] := Module[{spin = SPIN, parity = PARITY, prova, labels, elem},
    
    labels = listC[spin, parity];
    If[labels === {}, Return[{}]];
    
    prova = Table[
        Sum[
            elem = Which[
                (* i in rank3, j in rank3 *)
                MemberQ[listF[spin, parity, i3], i] && MemberQ[listF[spin, parity, i3], j],
                null[[car[spin, parity, i]]] * P[spin, parity, {i, j}][-m, -n, -f, -r, -s, -t] * AJ[r, s, t],
                
                (* i in rank2, j in rank3 *)
                MemberQ[listF[spin, parity, i2], i] && MemberQ[listF[spin, parity, i3], j],
                null[[car[spin, parity, i]]] * P[spin, parity, {i, j}][-m, -n, -f, -r, -s] * AJ[f, r, s],
                
                (* i in vector, j in rank3 *)
                MemberQ[listF[spin, parity, i1], i] && MemberQ[listF[spin, parity, i3], j],
                null[[car[spin, parity, i]]] * P[spin, parity, {i, j}][-m, -n, -f, -r] * AJ[n, f, r],
                
                (* i in scalar, j in rank3 *)
                MemberQ[listF[spin, parity, i0], i] && MemberQ[listF[spin, parity, i3], j],
                null[[car[spin, parity, i]]] * P[spin, parity, {i, j}][-m, -n, -f] * AJ[m, n, f],
                
                (* i in rank3, j in rank2 *)
                MemberQ[listF[spin, parity, i3], i] && MemberQ[listF[spin, parity, i2], j],
                null[[car[spin, parity, i]]] * P[spin, parity, {i, j}][-m, -n, -f, -r, -s] * HJ[r, s],
                
                (* i in rank2, j in rank2 *)
                MemberQ[listF[spin, parity, i2], i] && MemberQ[listF[spin, parity, i2], j],
                null[[car[spin, parity, i]]] * P[spin, parity, {i, j}][-m, -n, -f, -r] * HJ[f, r],
                
                (* i in vector, j in rank2 *)
                MemberQ[listF[spin, parity, i1], i] && MemberQ[listF[spin, parity, i2], j],
                null[[car[spin, parity, i]]] * P[spin, parity, {i, j}][-m, -n, -f] * HJ[n, f],
                
                (* i in scalar, j in rank2 *)
                MemberQ[listF[spin, parity, i0], i] && MemberQ[listF[spin, parity, i2], j],
                null[[car[spin, parity, i]]] * P[spin, parity, {i, j}][-m, -n] * HJ[m, n],
                
                (* i in rank3, j in vector *)
                MemberQ[listF[spin, parity, i3], i] && MemberQ[listF[spin, parity, i1], j],
                null[[car[spin, parity, i]]] * P[spin, parity, {i, j}][-m, -n, -f, -r] * VJ[r],
                
                (* i in rank2, j in vector *)
                MemberQ[listF[spin, parity, i2], i] && MemberQ[listF[spin, parity, i1], j],
                null[[car[spin, parity, i]]] * P[spin, parity, {i, j}][-m, -n, -f] * VJ[f],
                
                (* i in vector, j in vector *)
                MemberQ[listF[spin, parity, i1], i] && MemberQ[listF[spin, parity, i1], j],
                null[[car[spin, parity, i]]] * P[spin, parity, {i, j}][-m, -n] * VJ[n],
                
                (* i in scalar, j in vector *)
                MemberQ[listF[spin, parity, i0], i] && MemberQ[listF[spin, parity, i1], j],
                null[[car[spin, parity, i]]] * P[spin, parity, {i, j}][-m] * VJ[m],
                
                (* i in rank3, j in scalar *)
                MemberQ[listF[spin, parity, i3], i] && MemberQ[listF[spin, parity, i0], j],
                null[[car[spin, parity, i]]] * P[spin, parity, {i, j}][-m, -n, -f] * SJ,
                
                (* i in rank2, j in scalar *)
                MemberQ[listF[spin, parity, i2], i] && MemberQ[listF[spin, parity, i0], j],
                null[[car[spin, parity, i]]] * P[spin, parity, {i, j}][-m, -n] * SJ,
                
                (* i in vector, j in scalar *)
                MemberQ[listF[spin, parity, i1], i] && MemberQ[listF[spin, parity, i0], j],
                null[[car[spin, parity, i]]] * P[spin, parity, {i, j}][-m] * SJ,
                
                (* i in scalar, j in scalar *)
                MemberQ[listF[spin, parity, i0], i] && MemberQ[listF[spin, parity, i0], j],
                null[[car[spin, parity, i]]] * P[spin, parity, {i, j}][] * SJ,
                
                (* Default *)
                True, 0
            ];
            elem,
            {i, labels}
        ],
        {j, labels}
    ];
    
    prova = prova // ContractMetric // ReplaceAll[dimM -> 4] // FullSimplify;
    prova
];

defineSourcesForFields[fields_List] := Module[{},
    PrintK["=== Defining Source Tensors ==="];
    
    If[MemberQ[fields, "rank3"],
        PrintK["  AJ[a,b,c] - source for rank-3 tensor A"];
    ];
    If[MemberQ[fields, "rank2sym"],
        PrintK["  HJ[a,b] - source for rank-2 tensor H"];
    ];
    If[MemberQ[fields, "vector"],
        PrintK["  VJ[a] - source for vector V"];
    ];
    If[MemberQ[fields, "scalar"],
        PrintK["  SJ - source for scalar S"];
    ];
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

buildRemovalOrder[spin_, parity_] := Module[{order},
    order = Join[
        Reverse[Sort[listF[spin, parity, i3]]],  (* rank-3: remove first *)
        Reverse[Sort[listF[spin, parity, i2]]],  (* rank-2: remove second *)
        Reverse[Sort[listF[spin, parity, i1]]],  (* vector: remove third *)
        Reverse[Sort[listF[spin, parity, i0]]]   (* scalar: remove last *)
    ];
    order
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
    
    (* Build saturated propagator: sum over all label pairs *)
    prova = Sum[
        term = Which[
            (* i in rank3, j in rank3 *)
            MemberQ[listF[spin, parity, i3], i] && MemberQ[listF[spin, parity, i3], j],
            bm1[[carp[spin, parity, i], carp[spin, parity, j]]] *
                AJc[m, n, f] * P[spin, parity, {i, j}][-m, -n, -f, -r, -s, -t] * AJ[r, s, t],
            
            (* i in rank2, j in rank3 *)
            MemberQ[listF[spin, parity, i2], i] && MemberQ[listF[spin, parity, i3], j],
            bm1[[carp[spin, parity, i], carp[spin, parity, j]]] *
                HJc[m, n] * P[spin, parity, {i, j}][-m, -n, -f, -r, -s] * AJ[f, r, s],
            
            (* i in vector, j in rank3 *)
            MemberQ[listF[spin, parity, i1], i] && MemberQ[listF[spin, parity, i3], j],
            bm1[[carp[spin, parity, i], carp[spin, parity, j]]] *
                VJc[m] * P[spin, parity, {i, j}][-m, -n, -f, -r] * AJ[n, f, r],
            
            (* i in scalar, j in rank3 *)
            MemberQ[listF[spin, parity, i0], i] && MemberQ[listF[spin, parity, i3], j],
            bm1[[carp[spin, parity, i], carp[spin, parity, j]]] *
                SJc * P[spin, parity, {i, j}][-m, -n, -f] * AJ[m, n, f],
            
            (* i in rank3, j in rank2 *)
            MemberQ[listF[spin, parity, i3], i] && MemberQ[listF[spin, parity, i2], j],
            bm1[[carp[spin, parity, i], carp[spin, parity, j]]] *
                AJc[m, n, f] * P[spin, parity, {i, j}][-m, -n, -f, -r, -s] * HJ[r, s],
            
            (* i in rank2, j in rank2 *)
            MemberQ[listF[spin, parity, i2], i] && MemberQ[listF[spin, parity, i2], j],
            bm1[[carp[spin, parity, i], carp[spin, parity, j]]] *
                HJc[m, n] * P[spin, parity, {i, j}][-m, -n, -f, -r] * HJ[f, r],
            
            (* i in vector, j in rank2 *)
            MemberQ[listF[spin, parity, i1], i] && MemberQ[listF[spin, parity, i2], j],
            bm1[[carp[spin, parity, i], carp[spin, parity, j]]] *
                VJc[m] * P[spin, parity, {i, j}][-m, -n, -f] * HJ[n, f],
            
            (* i in scalar, j in rank2 *)
            MemberQ[listF[spin, parity, i0], i] && MemberQ[listF[spin, parity, i2], j],
            bm1[[carp[spin, parity, i], carp[spin, parity, j]]] *
                SJc * P[spin, parity, {i, j}][-m, -n] * HJ[m, n],
            
            (* i in rank3, j in vector *)
            MemberQ[listF[spin, parity, i3], i] && MemberQ[listF[spin, parity, i1], j],
            bm1[[carp[spin, parity, i], carp[spin, parity, j]]] *
                AJc[m, n, f] * P[spin, parity, {i, j}][-m, -n, -f, -r] * VJ[r],
            
            (* i in rank2, j in vector *)
            MemberQ[listF[spin, parity, i2], i] && MemberQ[listF[spin, parity, i1], j],
            bm1[[carp[spin, parity, i], carp[spin, parity, j]]] *
                HJc[m, n] * P[spin, parity, {i, j}][-m, -n, -f] * VJ[f],
            
            (* i in vector, j in vector *)
            MemberQ[listF[spin, parity, i1], i] && MemberQ[listF[spin, parity, i1], j],
            bm1[[carp[spin, parity, i], carp[spin, parity, j]]] *
                VJc[m] * P[spin, parity, {i, j}][-m, -n] * VJ[n],
            
            (* i in scalar, j in vector *)
            MemberQ[listF[spin, parity, i0], i] && MemberQ[listF[spin, parity, i1], j],
            bm1[[carp[spin, parity, i], carp[spin, parity, j]]] *
                SJc * P[spin, parity, {i, j}][-m] * VJ[m],
            
            (* i in rank3, j in scalar *)
            MemberQ[listF[spin, parity, i3], i] && MemberQ[listF[spin, parity, i0], j],
            bm1[[carp[spin, parity, i], carp[spin, parity, j]]] *
                AJc[m, n, f] * P[spin, parity, {i, j}][-m, -n, -f] * SJ,
            
            (* i in rank2, j in scalar *)
            MemberQ[listF[spin, parity, i2], i] && MemberQ[listF[spin, parity, i0], j],
            bm1[[carp[spin, parity, i], carp[spin, parity, j]]] *
                HJc[m, n] * P[spin, parity, {i, j}][-m, -n] * SJ,
            
            (* i in vector, j in scalar *)
            MemberQ[listF[spin, parity, i1], i] && MemberQ[listF[spin, parity, i0], j],
            bm1[[carp[spin, parity, i], carp[spin, parity, j]]] *
                VJc[m] * P[spin, parity, {i, j}][-m] * SJ,
            
            (* i in scalar, j in scalar *)
            MemberQ[listF[spin, parity, i0], i] && MemberQ[listF[spin, parity, i0], j],
            bm1[[carp[spin, parity, i], carp[spin, parity, j]]] *
                SJc * P[spin, parity, {i, j}][] * SJ,
            
            (* Default *)
            True, 0
        ];
        term,
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

buildSourceComponentLists[] := Module[{},
    
    (* Scalar: 1 component *)
    SJInd = {SJ};
    SJcInd = {SJc};
    
    (* Vector: 4 components *)
    VJInd = Table[VJ[{m, Lorentz}], {m, 0, 3}] // Flatten // DeleteDuplicates;
    VJcInd = Table[VJc[{m, Lorentz}], {m, 0, 3}] // Flatten // DeleteDuplicates;
    
    (* Rank-2 symmetric: 10 independent components *)
    HJInd = Table[HJ[{m, Lorentz}, {n, Lorentz}], {m, 0, 3}, {n, 0, 3}] // 
        Simplification // Flatten // DeleteDuplicates;
    HJcInd = Table[HJc[{m, Lorentz}, {n, Lorentz}], {m, 0, 3}, {n, 0, 3}] // 
        Simplification // Flatten // DeleteDuplicates;
    
    (* Rank-3 symmetric: 20 independent components *)
    AJInd = Table[AJ[{m, Lorentz}, {n, Lorentz}, {r, Lorentz}], 
        {m, 0, 3}, {n, 0, 3}, {r, 0, 3}] // 
        Simplification // Flatten // ReplaceAll[-x__ :> x] // 
        DeleteDuplicates // DeleteCases[#, 0] &;
    AJcInd = Table[AJc[{m, Lorentz}, {n, Lorentz}, {r, Lorentz}], 
        {m, 0, 3}, {n, 0, 3}, {r, 0, 3}] // 
        Simplification // Flatten // ReplaceAll[-x__ :> x] // 
        DeleteDuplicates // DeleteCases[#, 0] &;
    
    (* All source variables *)
    (* allSourceVars = Join[SJInd, VJInd, HJInd, AJInd];
     allSourceVarsConj = Join[SJcInd, VJcInd, HJcInd, AJcInd];*)
    
    (*
    allSourceVars = AJInd;
    allSourceVarsConj = AJcInd;
    *)
    
    allSourceVars={};
    allSourceVarsConj={};
 
 If[MemberQ[activeFieldTypes,"rank3"], AppendTo[allSourceVars,AJInd]; AppendTo[allSourceVarsConj,AJcInd];];
 If[MemberQ[activeFieldTypes,"rank2sym"], AppendTo[allSourceVars,HJInd]; AppendTo[allSourceVarsConj,HJcInd];];
 If[MemberQ[activeFieldTypes,"vector"], AppendTo[allSourceVars,VJInd]; AppendTo[allSourceVarsConj,VJcInd];];
 If[MemberQ[activeFieldTypes,"scalar"], AppendTo[allSourceVars,SJInd]; AppendTo[allSourceVarsConj,SJcInd];];
 
 allSourceVars=Flatten[allSourceVars];
 allSourceVarsConj=Flatten[allSourceVarsConj];
 
 
(*    PrintK[allSourceVars];
    PrintK[allSourceVarsConj];
*)
    
    
 (*   allSourceVars = Join[SJInd, VJInd, HJInd,AJInd];
    allSourceVarsConj = Join[SJcInd, VJcInd, HJcInd,AJcInd];
 *)
(*    
    PrintK["Source components built:"];
    PrintK["  Scalar: ", Length[SJInd]];
    PrintK["  Vector: ", Length[VJInd]];
    PrintK["  Rank-2: ", Length[HJInd]];
    PrintK["  Rank-3: ", Length[AJInd]];
    PrintK["  Total: ", Length[allSourceVars]];
*)

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
    ReplaceAll[{
        AJ[ind__] :> AJc[ind], 
        HJ[ind__] :> HJc[ind], 
        VJ[ind__] :> VJc[ind], 
        SJ -> SJc
    }];

(* Convert a list of rules to conjugate rules *)
conjugateRules[rules_List] := Module[{newRules},
    newRules = Table[
        Rule[
            rule[[1]] /. {AJ[ind__] :> AJc[ind], HJ[ind__] :> HJc[ind], 
                          VJ[ind__] :> VJc[ind], SJ -> SJc},
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
    
    (* Build source component lists if not already done *)
    If[!ValueQ[allSourceVars] || allSourceVars === {},
        buildSourceComponentLists[];
    ];
    
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



(*Source declaration in component terms take too long. Do it once for all.*)
SourceDeclare[action_]:=
Block[{Print = (Null &)},
activeFieldTypes = detectFields[action](*{"rank2sym"(*, "vector", "scalar"*)}*);
initializeSourceComponentRules[];
];


CompSatProp[action_]:=Module[{actionPert},

activeFieldTypesSecondTime = detectFields[action];  

If[activeFieldTypesSecondTime =!= activeFieldTypes,

activeFieldTypes = detectFields[action];
initializeSourceComponentRules[];
];


BuildRepresentationTables[activeFieldTypes];

(* View configuration *)
PrintConfiguration[];


(* Apply to action *)

actionPert = action /. ToPerturbations;

TAA = VarAction[
    VarAction[actionPert // ReplaceDummies, \[Epsilon]Ap[LI[2], -m, -n, -f]], 
    \[Epsilon]Ap[LI[3], -r, -s, -t]] // CollectTensors;

THA = VarAction[
    VarAction[actionPert // ReplaceDummies, \[Epsilon]h[LI[2], -m, -n]], 
    \[Epsilon]Ap[LI[3], -f, -r, -s]] // CollectTensors;

TVA = VarAction[
    VarAction[actionPert // ReplaceDummies, \[Epsilon]Vp[LI[2], -m]], 
    \[Epsilon]Ap[LI[3], -n, -f, -r]] // CollectTensors;

TSA = VarAction[
    VarAction[actionPert // ReplaceDummies, \[Epsilon]Sp[LI[2]]], 
    \[Epsilon]Ap[LI[3], -m, -n, -f]] // CollectTensors;

TAH = VarAction[
    VarAction[actionPert // ReplaceDummies, \[Epsilon]Ap[LI[2], -m, -n, -f]], 
    \[Epsilon]h[LI[3], -r, -s]] // CollectTensors;

THH = VarAction[
    VarAction[actionPert // ReplaceDummies, \[Epsilon]h[LI[2], -m, -n]], 
    \[Epsilon]h[LI[3], -f, -r]] // CollectTensors;

TVH = VarAction[
    VarAction[actionPert // ReplaceDummies, \[Epsilon]Vp[LI[2], -m]], 
    \[Epsilon]h[LI[3], -n, -f]] // CollectTensors;

TSH = VarAction[
    VarAction[actionPert // ReplaceDummies, \[Epsilon]Sp[LI[2]]], 
    \[Epsilon]h[LI[3], -m, -n]] // CollectTensors;

TAV = VarAction[
    VarAction[actionPert // ReplaceDummies, \[Epsilon]Ap[LI[2], -m, -n, -f]], 
    \[Epsilon]Vp[LI[3], -r]] // CollectTensors;

THV = VarAction[
    VarAction[actionPert // ReplaceDummies, \[Epsilon]h[LI[2], -m, -n]], 
    \[Epsilon]Vp[LI[3], -f]] // CollectTensors;

TVV = VarAction[
    VarAction[actionPert // ReplaceDummies, \[Epsilon]Vp[LI[2], -m]], 
    \[Epsilon]Vp[LI[3], -n]] // CollectTensors;

TSV = VarAction[
    VarAction[actionPert // ReplaceDummies, \[Epsilon]Sp[LI[2]]], 
    \[Epsilon]Vp[LI[3], -m]] // CollectTensors;

TAS = VarAction[
    VarAction[actionPert // ReplaceDummies, \[Epsilon]Ap[LI[2], -m, -n, -f]], 
    \[Epsilon]Sp[LI[3]]] // CollectTensors;

THS = VarAction[
    VarAction[actionPert // ReplaceDummies, \[Epsilon]h[LI[2], -m, -n]], 
    \[Epsilon]Sp[LI[3]]] // CollectTensors;

TVS = VarAction[
    VarAction[actionPert // ReplaceDummies, \[Epsilon]Vp[LI[2], -m]], 
    \[Epsilon]Sp[LI[3]]] // CollectTensors;

TSS = VarAction[
    VarAction[actionPert // ReplaceDummies, \[Epsilon]Sp[LI[2]]], 
    \[Epsilon]Sp[LI[3]]] // CollectTensors;


(* Apply to all T-matrices *)
TAA = TAA //. ToMomentumSpaceRules;
THA = THA //. ToMomentumSpaceRules;
TVA = TVA //. ToMomentumSpaceRules;
TSA = TSA //. ToMomentumSpaceRules;

TAH = TAH //. ToMomentumSpaceRules;
THH = THH //. ToMomentumSpaceRules;
TVH = TVH //. ToMomentumSpaceRules;
TSH = TSH //. ToMomentumSpaceRules;

TAV = TAV //. ToMomentumSpaceRules;
THV = THV //. ToMomentumSpaceRules;
TVV = TVV //. ToMomentumSpaceRules;
TSV = TSV //. ToMomentumSpaceRules;

TAS = TAS //. ToMomentumSpaceRules;
THS = THS //. ToMomentumSpaceRules;
TVS = TVS //. ToMomentumSpaceRules;
TSS = TSS //. ToMomentumSpaceRules;

(* Initialize all listF to empty *)
initializeListF[] := Module[{},
    (* 0- sector *)
    listF[0, m, i3] = {}; listF[0, m, i2] = {}; listF[0, m, i1] = {}; listF[0, m, i0] = {};
    (* 0+ sector *)
    listF[0, p, i3] = {}; listF[0, p, i2] = {}; listF[0, p, i1] = {}; listF[0, p, i0] = {};
    (* 1- sector *)
    listF[1, m, i3] = {}; listF[1, m, i2] = {}; listF[1, m, i1] = {}; listF[1, m, i0] = {};
    (* 1+ sector *)
    listF[1, p, i3] = {}; listF[1, p, i2] = {}; listF[1, p, i1] = {}; listF[1, p, i0] = {};
    (* 2- sector *)
    listF[2, m, i3] = {}; listF[2, m, i2] = {}; listF[2, m, i1] = {}; listF[2, m, i0] = {};
    (* 2+ sector *)
    listF[2, p, i3] = {}; listF[2, p, i2] = {}; listF[2, p, i1] = {}; listF[2, p, i0] = {};
    (* 3- sector *)
    listF[3, m, i3] = {}; listF[3, m, i2] = {}; listF[3, m, i1] = {}; listF[3, m, i0] = {};
];

(* Populate listF based on detected fields AND rank-3 symmetry type *)
populateListF[fields_List] := Module[{},
    
    (* Initialize all to empty *)
    initializeListF[];
    
    (* Rank-3 tensor labels - depend on symmetry type *)
    If[MemberQ[fields, "rank3"],
        Switch[Global`rank3SymmetryType,
            
            "Generic",
            (* All sectors, full labels *)
            listF[3, m, i3] = {1};
            listF[2, p, i3] = {1, 2, 3};
            listF[2, m, i3] = {1, 2};
            listF[1, p, i3] = {1, 2, 3};
            listF[1, m, i3] = {1, 2, 3, 4, 5, 6};
            listF[0, p, i3] = {1, 2, 3, 4};
            listF[0, m, i3] = {1};,
            
            "Sym23",
            (* Symmetric in last 2 indices *)
            listF[3, m, i3] = {1};
            listF[2, p, i3] = {1, 2};
            listF[2, m, i3] = {1};
            listF[1, p, i3] = {1};
            listF[1, m, i3] = {1, 2, 4, 5};
            listF[0, p, i3] = {1, 2, 4};
            listF[0, m, i3] = {};,
            
            "Antisym23",
            (* Antisymmetric in last 2 indices *)
            listF[3, m, i3] = {};
            listF[2, p, i3] = {3};
            listF[2, m, i3] = {2};
            listF[1, p, i3] = {2, 3};
            listF[1, m, i3] = {3, 6};
            listF[0, p, i3] = {3};
            listF[0, m, i3] = {1};,
            
            "TotallySymmetric",
            (* Fully symmetric *)
            listF[3, m, i3] = {1};
            listF[2, p, i3] = {1};
            listF[2, m, i3] = {1};
            listF[1, p, i3] = {1};
            listF[1, m, i3] = {1, 4};
            listF[0, p, i3] = {1, 4};
            listF[0, m, i3] = {};,
            
            "TotallyAntisymmetric",
            (* Fully antisymmetric *)

            listF[3, m, i3] = {};
            listF[2, p, i3] = {};
            listF[2, m, i3] = {};
            listF[1, p, i3] = {3};  
            listF[1, m, i3] = {};  
            listF[0, p, i3] = {};  
            listF[0, m, i3] = {1};  
            ,
            _,
            (* Default to Generic *)
            listF[3, m, i3] = {1};
            listF[2, p, i3] = {1, 2, 3};
            listF[2, m, i3] = {1, 2};
            listF[1, p, i3] = {1, 2, 3};
            listF[1, m, i3] = {1, 2, 3, 4, 5, 6};
            listF[0, p, i3] = {1, 2, 3, 4};
            listF[0, m, i3] = {1};
        ];
    ];
    
    (* Rank-2 symmetric tensor labels - same for all cases *)
    If[MemberQ[fields, "rank2sym"],
        listF[0, p, i2] = {5, 6};
        listF[1, m, i2] = {7};
        listF[2, p, i2] = {4};
    ];
    
    (* Vector labels - same for all cases *)
    If[MemberQ[fields, "vector"],
        listF[0, p, i1] = {7};
        listF[1, m, i1] = {8};
    ];
    
    (* Scalar labels - same for all cases *)
    If[MemberQ[fields, "scalar"],
        listF[0, p, i0] = {8};
    ];
    
    (* Print summary *)
    PrintK["=== Field Content Configuration ==="];
    PrintK["Detected fields: ", fields];
    PrintK["Rank-3 symmetry type: ", Global`rank3SymmetryType];
    PrintK[""];
    PrintK["Active sectors:"];
    If[listC[3, m] =!= {}, PrintK["  3-: labels ", listC[3, m]]];
    If[listC[2, p] =!= {}, PrintK["  2+: labels ", listC[2, p]]];
    If[listC[2, m] =!= {}, PrintK["  2-: labels ", listC[2, m]]];
    If[listC[1, p] =!= {}, PrintK["  1+: labels ", listC[1, p]]];
    If[listC[1, m] =!= {}, PrintK["  1-: labels ", listC[1, m]]];
    If[listC[0, p] =!= {}, PrintK["  0+: labels ", listC[0, p]]];
    If[listC[0, m] =!= {}, PrintK["  0-: labels ", listC[0, m]]];
];

(* One-step setup: detect fields and populate *)
setupFieldContent[act_] := Module[{fields},
    fields = detectFields[act];
    populateListF[fields];
    fields
];

(* Dimension of spin-j representation *)
dimRep[j_] := 2 j + 1;

(* Flatten to get full list of labels for a sector *)
listC[spin_, parity_] := Flatten[{listF[spin, parity, i3], listF[spin, parity, i2], listF[spin, parity, i1], listF[spin, parity, i0]}];

(* List of all sectors *)
listRep = {{3, m}, {2, p}, {2, m}, {1, p}, {1, m}, {0, p}, {0, m}};

(* Get active sectors (non-empty listC) *)
activeSectors[] := Select[listRep, listC[#[[1]], #[[2]]] =!= {} &];

amatrix2[SPIN_, PARITY_] := Module[{spin = SPIN, parity = PARITY, prova, labels, elem},
    labels = listC[spin, parity];
    
    If[labels === {},
        PrintK["No labels for sector ", spin, "^", If[parity === p, "+", "-"]];
        Return[{}]
    ];
    
    prova = Table[
        elem = Which[
            (* rank3-rank3: 6 indices *)
            MemberQ[listF[spin, parity, i3], i] && MemberQ[listF[spin, parity, i3], j],
            P[spin, parity, {i, j}][-m, -n, -f, -r, -s, -t] * TAA,
            
            (* rank2-rank3: 5 indices *)
            MemberQ[listF[spin, parity, i2], i] && MemberQ[listF[spin, parity, i3], j],
            P[spin, parity, {i, j}][-m, -n, -f, -r, -s] * THA,
            
            (* vector-rank3: 4 indices *)
            MemberQ[listF[spin, parity, i1], i] && MemberQ[listF[spin, parity, i3], j],
            P[spin, parity, {i, j}][-m, -n, -f, -r] * TVA,
            
            (* scalar-rank3: 3 indices *)
            MemberQ[listF[spin, parity, i0], i] && MemberQ[listF[spin, parity, i3], j],
            P[spin, parity, {i, j}][-m, -n, -f] * TSA,
            
            (* rank3-rank2: 5 indices *)
            MemberQ[listF[spin, parity, i3], i] && MemberQ[listF[spin, parity, i2], j],
            P[spin, parity, {i, j}][-m, -n, -f, -r, -s] * TAH,
            
            (* rank2-rank2: 4 indices *)
            MemberQ[listF[spin, parity, i2], i] && MemberQ[listF[spin, parity, i2], j],
            P[spin, parity, {i, j}][-m, -n, -f, -r] * THH,
            
            (* vector-rank2: 3 indices *)
            MemberQ[listF[spin, parity, i1], i] && MemberQ[listF[spin, parity, i2], j],
            P[spin, parity, {i, j}][-m, -n, -f] * TVH,
            
            (* scalar-rank2: 2 indices *)
            MemberQ[listF[spin, parity, i0], i] && MemberQ[listF[spin, parity, i2], j],
            P[spin, parity, {i, j}][-m, -n] * TSH,
            
            (* rank3-vector: 4 indices *)
            MemberQ[listF[spin, parity, i3], i] && MemberQ[listF[spin, parity, i1], j],
            P[spin, parity, {i, j}][-m, -n, -f, -r] * TAV,
            
            (* rank2-vector: 3 indices *)
            MemberQ[listF[spin, parity, i2], i] && MemberQ[listF[spin, parity, i1], j],
            P[spin, parity, {i, j}][-m, -n, -f] * THV,
            
            (* vector-vector: 2 indices *)
            MemberQ[listF[spin, parity, i1], i] && MemberQ[listF[spin, parity, i1], j],
            P[spin, parity, {i, j}][-m, -n] * TVV,
            
            (* scalar-vector: 1 index *)
            MemberQ[listF[spin, parity, i0], i] && MemberQ[listF[spin, parity, i1], j],
            P[spin, parity, {i, j}][-m] * TSV,
            
            (* rank3-scalar: 3 indices *)
            MemberQ[listF[spin, parity, i3], i] && MemberQ[listF[spin, parity, i0], j],
            P[spin, parity, {i, j}][-m, -n, -f] * TAS,
            
            (* rank2-scalar: 2 indices *)
            MemberQ[listF[spin, parity, i2], i] && MemberQ[listF[spin, parity, i0], j],
            P[spin, parity, {i, j}][-m, -n] * THS,
            
            (* vector-scalar: 1 index *)
            MemberQ[listF[spin, parity, i1], i] && MemberQ[listF[spin, parity, i0], j],
            P[spin, parity, {i, j}][-m] * TVS,
            
            (* scalar-scalar: 0 indices *)
            MemberQ[listF[spin, parity, i0], i] && MemberQ[listF[spin, parity, i0], j],
            P[spin, parity, {i, j}][] * TSS,
            
            (* Default: no match *)
            True,
            PrintK["No match for i=", i, ", j=", j]; 0
        ];
        elem,
        {i, labels},
        {j, labels}
    ];
    
    (* Normalize by dimension and simplify *)
    prova = prova / dimRep[spin] // ContractMetric // ReplaceAll[dimM -> 4] // FullSimplify;
    prova
];

detectedFields = setupFieldContent[action];

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
basisHConj=basisH/.{AJ[x__]:>AJc[x],HJ[x__]:>HJc[x],VJ[x__]:>VJc[x],SJ->SJc};

(*Print[{masslessSolution, masslessSolutionConj}];
*)
LatestPropMassless={};

If[Global`DivideEtImpera,

(*DivideEtImpera*)

propMasslessLimitList=List@@propMasslessLimit;

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
basisHConj=basisH/.{AJ[x__]:>AJc[x],HJ[x__]:>HJc[x],VJ[x__]:>VJc[x],SJ->SJc};

LatestPropMassive={};
massiveData={};

If[Global`DivideEtImpera,

(*DivideEtImpera*)
propMassiveLimitList=List@@propMassiveLimit;

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

TensorToHelicity[]:=Module[{},


helicityGenericRulesPart1={sv[-1, 1] == I*VJ[{1, Lorentz}] + VJ[{2, Lorentz}], sv[1, 2] == (-I)*VJ[{1, Lorentz}] + VJ[{2, Lorentz}], sv[0, 3] == VJ[{3, Lorentz}], sv[0, 4] == VJ[{0, Lorentz}], svc[-1, 1] == (-I)*VJc[{1, Lorentz}] + VJc[{2, Lorentz}], 
 svc[1, 2] == I*VJc[{1, Lorentz}] + VJc[{2, Lorentz}], svc[0, 3] == VJc[{3, Lorentz}], svc[0, 4] == VJc[{0, Lorentz}], 
 sh[-2, 1] == -HJ[{1, Lorentz}, {1, Lorentz}] + (2*I)*HJ[{1, Lorentz}, {2, Lorentz}] + HJ[{2, Lorentz}, {2, Lorentz}], 
 sh[2, 2] == -HJ[{1, Lorentz}, {1, Lorentz}] - (2*I)*HJ[{1, Lorentz}, {2, Lorentz}] + HJ[{2, Lorentz}, {2, Lorentz}], sh[-1, 3] == I*HJ[{1, Lorentz}, {3, Lorentz}] + HJ[{2, Lorentz}, {3, Lorentz}], 
 sh[-1, 4] == I*HJ[{0, Lorentz}, {1, Lorentz}] + HJ[{0, Lorentz}, {2, Lorentz}], sh[1, 5] == (-I)*HJ[{1, Lorentz}, {3, Lorentz}] + HJ[{2, Lorentz}, {3, Lorentz}], 
 sh[1, 6] == (-I)*HJ[{0, Lorentz}, {1, Lorentz}] + HJ[{0, Lorentz}, {2, Lorentz}], sh[0, 7] == HJ[{3, Lorentz}, {3, Lorentz}], sh[0, 8] == HJ[{1, Lorentz}, {1, Lorentz}] + HJ[{2, Lorentz}, {2, Lorentz}], 
 sh[0, 9] == HJ[{0, Lorentz}, {3, Lorentz}], sh[0, 10] == HJ[{0, Lorentz}, {0, Lorentz}], shc[-2, 1] == -HJc[{1, Lorentz}, {1, Lorentz}] - (2*I)*HJc[{1, Lorentz}, {2, Lorentz}] + HJc[{2, Lorentz}, {2, Lorentz}], 
 shc[2, 2] == -HJc[{1, Lorentz}, {1, Lorentz}] + (2*I)*HJc[{1, Lorentz}, {2, Lorentz}] + HJc[{2, Lorentz}, {2, Lorentz}], shc[-1, 3] == (-I)*HJc[{1, Lorentz}, {3, Lorentz}] + HJc[{2, Lorentz}, {3, Lorentz}], 
 shc[-1, 4] == (-I)*HJc[{0, Lorentz}, {1, Lorentz}] + HJc[{0, Lorentz}, {2, Lorentz}], shc[1, 5] == I*HJc[{1, Lorentz}, {3, Lorentz}] + HJc[{2, Lorentz}, {3, Lorentz}], 
 shc[1, 6] == I*HJc[{0, Lorentz}, {1, Lorentz}] + HJc[{0, Lorentz}, {2, Lorentz}], shc[0, 7] == HJc[{3, Lorentz}, {3, Lorentz}], shc[0, 8] == HJc[{1, Lorentz}, {1, Lorentz}] + HJc[{2, Lorentz}, {2, Lorentz}], 
 shc[0, 9] == HJc[{0, Lorentz}, {3, Lorentz}], shc[0, 10] == HJc[{0, Lorentz}, {0, Lorentz}]};


helicityGenericRulesPart2 = 
 Switch[Global`rank3SymmetryType,
        
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
        Print["Unknown symmetry type: ", symType, ". Using Generic."];
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
    
helicityGenericRules = Join[helicityGenericRulesPart1,helicityGenericRulesPart2];
(*These undergo a further reduction due to source constraints*)

Reduce[FullSimplify[helicityGenericRules/.masslessSolution/.masslessSolutionConj/.\[Kappa]->\[Omega]],Join[allSourceVars,allSourceVarsConj]]//ToRules

]




End[] (* `Private` *)
EndPackage[]
