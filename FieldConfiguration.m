(* ::Package:: *)

(* FieldConfiguration.m *)
(* Module for defining field content and building representation tables *)
(* Part of the Spectral Analysis package for quadratic Lagrangians *)

BeginPackage["SpectralAnalysis`FieldConfiguration`"];

(* === Public symbols === *)

DefineFieldType::usage = "DefineFieldType[name, rank, symmetry] registers a new field type.";
GetFieldType::usage = "GetFieldType[name] returns the field type specification.";
ListFieldTypes::usage = "ListFieldTypes[] lists all registered field types.";

SetRepContent::usage = "SetRepContent[fieldType, jpSector, labels] sets the representation content.";
GetRepContent::usage = "GetRepContent[fieldType] returns all J^P content for a field type.";

BuildRepresentationTables::usage = "BuildRepresentationTables[fieldList] builds listRep, listC, listF from the specified fields.";
GetListRep::usage = "GetListRep[] returns the current listRep.";
GetListC::usage = "GetListC[J, parity] returns listC for the given J^P sector.";
GetListF::usage = "GetListF[J, parity, fieldSlot] returns listF for the given J^P sector and field slot.";

ClearConfiguration::usage = "ClearConfiguration[] resets all configuration to empty state.";
PrintConfiguration::usage = "PrintConfiguration[] displays the current configuration.";

(* Slot identifiers *)
i3::usage = "i3 represents the rank-3 field slot (slot 1).";
i2::usage = "i2 represents the rank-2 field slot (slot 2).";
i1::usage = "i1 represents the rank-1/vector field slot (slot 3).";
i0::usage = "i0 represents the rank-0/scalar field slot (slot 4).";

(* Parity symbols *)
(*
p::usage = "p represents positive parity.";
m::usage = "m represents negative parity.";
*)
Begin["`Private`"];

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
    Print["=== Current Field Configuration ==="];
    Print["Fields: ", $CurrentFields];
    Print[""];
    Print["J^P Sectors (listRep): ", $CurrentListRep];
    Print[""];
    Print["Representation Tables:"];
    Do[
        Print["  listC[", jp[[1]], ",", jp[[2]], "] = ", $CurrentListC[jp]],
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

End[];
EndPackage[];
