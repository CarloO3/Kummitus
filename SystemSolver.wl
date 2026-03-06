(* ::Package:: *)

(* ============================================================
   SystemSolver.wl \[LongDash] Modular routines for solving symbolic
   linear systems with automatic simplification and
   singularity analysis.

   Load with:
     Get["SystemSolver.wl"]
   or
     << SystemSolver`

   Main entry point:
     result = SolveSystemFull[rawEqs, vars, opts]

   Individual routines are also exported for pipeline use.
   ============================================================ *)

BeginPackage["SystemSolver`"];

(* ---- Public symbols ---- *)

SimplifyEquations::usage =
  "SimplifyEquations[rawEqs, vars] takes a list of equations (expr == 0) \
and returns a reduced list: denominators cleared, parameter-only factors \
stripped, duplicates and proportional equations removed.";

AnalyzeLinearity::usage =
  "AnalyzeLinearity[eqs, vars] returns an Association with keys \
\"Linear\" (True/False), \"CoefficientMatrix\", \"Rank\", \"FreeCount\".";

CheckSingularRules::usage =
  "CheckSingularRules[sol, replaceRule, assumps] checks each rule in a \
solution branch for singularities under the substitution replaceRule \
(e.g. \[Kappa] -> \[Omega]).  Returns a list of singular rules.";

ChooseNonSingularPivots::usage =
  "ChooseNonSingularPivots[coeffMatrix, vars, rk, replaceRule] \
identifies dependent variables whose solution is regular under the \
given substitution.  Returns {depVars, freeVars}.";

SolveReduced::usage =
  "SolveReduced[eqs, vars, assumps] solves the (already simplified) \
equations for vars.  Returns the solution list.";

SolveSystemFull::usage =
  "SolveSystemFull[rawEqs, vars, opts] is the main driver.  Options: \
\"Assumptions\" -> {}, \"SingularityCheck\" -> (\[Kappa] -> \[Omega]), \
\"Verbose\" -> True.  Returns an Association with all results.";

Begin["`Private`"];

(* ============================================================
   1. SimplifyEquations
   ============================================================ *)

stripParamFactors[expr_, vars_List] := Module[
  {fac, keep},
  fac = FactorList[expr];
  keep = Select[fac, !FreeQ[#[[1]], Alternatives @@ vars] &];
  If[keep === {},
    expr,
    Times @@ ((#[[1]]^#[[2]]) & /@ keep) // Expand
  ]
];

removeDuplicates[exprList_List, vars_List] := Module[
  {unique = {}, dominated, ratio},
  Do[
    If[e === 0, Continue[]];
    dominated = False;
    Do[
      ratio = Cancel[Together[e / u]];
      If[FreeQ[ratio, Alternatives @@ vars],
        dominated = True; Break[]
      ],
      {u, unique}
    ];
    If[!dominated, AppendTo[unique, e]],
    {e, exprList}
  ];
  unique
];

SimplifyEquations[rawEqs_List, vars_List] := Module[
  {numerators, cleaned, distinct,
   nTrivial, nDuplicates},

  (* Extract LHS, clear denominators *)
  numerators = Together /@ (rawEqs /. (lhs_ == 0) :> lhs);
  numerators = Numerator /@ numerators;
  numerators = Expand /@ numerators;

  (* Strip multiplicative factors that are free of the unknowns *)
  cleaned = stripParamFactors[#, vars] & /@ numerators;

  (* Remove trivial (0==0) and proportional duplicates *)
  nTrivial = Count[cleaned, 0];
  distinct = removeDuplicates[cleaned, vars];
  nDuplicates = Length[cleaned] - nTrivial - Length[distinct];

  <|
    "Equations"  -> distinct,
    "InputCount" -> Length[rawEqs],
    "TrivialCount" -> nTrivial,
    "DuplicateCount" -> nDuplicates,
    "DistinctCount" -> Length[distinct]
  |>
];

(* ============================================================
   2. AnalyzeLinearity
   ============================================================ *)

AnalyzeLinearity[eqs_List, vars_List] := Module[
  {isLinear, coeffMat, rk},

  isLinear = And @@ (Function[eq,
    And @@ (Exponent[eq, #] <= 1 & /@ vars) &&
    FreeQ[eq /. Thread[vars -> 0], Alternatives @@ vars]
  ] /@ eqs);

  If[isLinear,
    coeffMat = Outer[Function[{eq, v}, Coefficient[eq, v]], eqs, vars];
    coeffMat = Together /@ coeffMat // Simplify;
    rk = MatrixRank[coeffMat];
    <|
      "Linear" -> True,
      "CoefficientMatrix" -> coeffMat,
      "Rank" -> rk,
      "FreeCount" -> Length[vars] - rk
    |>,
    (* Nonlinear *)
    <|
      "Linear" -> False,
      "CoefficientMatrix" -> None,
      "Rank" -> None,
      "FreeCount" -> None
    |>
  ]
];

(* ============================================================
   3. CheckSingularRules
   ============================================================ *)

CheckSingularRules[sol_List, replaceRule_Rule, assumps_List: {}] := Module[
  {singulars = {}, den, denAt, lim, param, val},

  param = replaceRule[[1]];
  val   = replaceRule[[2]];

  Do[
    den = Denominator[Together[rule[[2]]]];
    denAt = den /. replaceRule // Simplify;
    If[TrueQ[denAt == 0],
      lim = Quiet[Limit[rule[[2]], param -> val, Assumptions -> assumps]];
      If[!FreeQ[lim, DirectedInfinity | Indeterminate | ComplexInfinity],
        AppendTo[singulars, rule]
      ]
    ],
    {rule, sol}
  ];

  singulars
];

(* ============================================================
   4. ChooseNonSingularPivots
   ============================================================ *)

ChooseNonSingularPivots[coeffMatrix_, vars_List, rk_Integer,
                         replaceRule_Rule] := Module[
  {coeffAt, rkAt, rref, pivotCols = {}, rrefGen, depVars, freeVars},

  coeffAt = coeffMatrix /. replaceRule // Simplify;
  rkAt = MatrixRank[coeffAt];

  (* Row-reduce at the special point to find stable pivots *)
  rref = RowReduce[coeffAt];
  Do[
    Do[
      If[rref[[i, j]] =!= 0 && TrueQ[Simplify[rref[[i, j]] != 0]],
        AppendTo[pivotCols, j]; Break[]
      ],
      {j, Length[vars]}
    ],
    {i, Length[rref]}
  ];
  pivotCols = DeleteDuplicates[pivotCols];

  (* If the rank drops at the special point, augment from generic pivots *)
  If[Length[pivotCols] < rk,
    rrefGen = RowReduce[coeffMatrix // Simplify];
    Do[
      Do[
        If[rrefGen[[i, j]] =!= 0 && TrueQ[Simplify[rrefGen[[i, j]] != 0]],
          AppendTo[pivotCols, j]; Break[]
        ],
        {j, Length[vars]}
      ],
      {i, Length[rrefGen]}
    ];
    pivotCols = DeleteDuplicates[pivotCols];
    pivotCols = Take[pivotCols, Min[rk, Length[pivotCols]]];
  ];

  depVars  = vars[[pivotCols]];
  freeVars = Complement[vars, depVars];

  <|
    "DependentVars" -> depVars,
    "FreeVars" -> freeVars,
    "PivotColumns" -> pivotCols,
    "RankAtPoint" -> rkAt
  |>
];

(* ============================================================
   5. SolveReduced
   ============================================================ *)

SolveReduced[eqs_List, vars_List, assumps_List: {}] := Module[
  {eqsToSolve, sol},

  eqsToSolve = (# == 0) & /@ eqs;
  sol = Solve[eqsToSolve, vars];
  sol = Simplify[#, Assumptions -> assumps] & /@ sol;
  sol
];

(* ============================================================
   6. SolveSystemFull \[LongDash] main driver
   ============================================================ *)

Options[SolveSystemFull] = {
  "Assumptions" -> {},
  "SingularityCheck" -> None,
  "Verbose" -> True
};

SolveSystemFull[rawEqs_List, vars_List, OptionsPattern[]] := Module[
  {assumps, singRule, verbose,
   simpResult, eqs,
   linResult, coeffMat, rk,
   sol, singulars, pivotResult, depVars,
   verified, residuals, allTrue,
   result},

  assumps  = OptionValue["Assumptions"];
  singRule = OptionValue["SingularityCheck"];
  verbose  = OptionValue["Verbose"];

  (* --- Step 1: Simplify --- *)
  simpResult = SimplifyEquations[rawEqs, vars];
  eqs = simpResult["Equations"];

  If[verbose,
    Print["=== Simplification ==="];
    Print["  Input equations    : ", simpResult["InputCount"]];
    Print["  Trivially zero     : ", simpResult["TrivialCount"]];
    Print["  Duplicates removed : ", simpResult["DuplicateCount"]];
    Print["  Distinct equations : ", simpResult["DistinctCount"]];
    Print[""];
  ];

  If[Length[eqs] == 0,
    If[verbose, Print["All equations trivially satisfied."]];
    Return[<|"Solution" -> {{}}, "Singular" -> False,
              "Simplification" -> simpResult, "Linearity" -> None|>]
  ];

  (* --- Step 2: Linearity analysis --- *)
  linResult = AnalyzeLinearity[eqs, vars];
  coeffMat  = linResult["CoefficientMatrix"];
  rk        = linResult["Rank"];

  If[verbose,
    Print["=== Linearity ==="];
    Print["  Linear  : ", linResult["Linear"]];
    If[linResult["Linear"],
      Print["  Rank    : ", rk, "  (out of ", Length[vars], ")"];
      Print["  Independent    : ", linResult["FreeCount"]];
    ];
    Print[""];
  ];

  (* --- Step 3: Initial solve --- *)
  If[verbose, Print["=== Solving ==="]];
  sol = SolveReduced[eqs, vars, assumps];
  If[verbose,
    Print["  Branches: ", Length[sol]];
    Print[""];
  ];

  (* --- Step 4: Singularity check --- *)
  result = <|
    "Simplification" -> simpResult,
    "Linearity" -> linResult,
    "Solution" -> sol,
    "Singular" -> False,
    "SingularRules" -> {},
    "PivotInfo" -> None
  |>;

  If[singRule =!= None && Length[sol] > 0,
    If[verbose, Print["=== Singularity check: ", singRule, " ==="]];

    singulars = CheckSingularRules[sol[[1]], singRule, assumps];

    If[Length[singulars] > 0,
      If[verbose,
        Print["  SINGULAR rules:"];
        Do[Print["    ", r[[1]]], {r, singulars}];
        Print[""];
      ];
      result["Singular"] = True;
      result["SingularRules"] = singulars;

      (* --- Step 4b: Re-solve --- *)
      If[linResult["Linear"],
        If[verbose, Print["  Choosing non-singular pivots ..."]];
        pivotResult = ChooseNonSingularPivots[coeffMat, vars, rk, singRule];
        depVars = pivotResult["DependentVars"];

        If[verbose,
          Print["  Rank at point : ", pivotResult["RankAtPoint"]];
          Print["  Dependent vars: ", depVars];
          Print["  Free vars     : ", pivotResult["FreeVars"]];
          Print[""];
          Print["  Re-solving ..."];
        ];

        sol = SolveReduced[eqs, depVars, assumps];

        (* Verify non-singularity *)
        singulars = If[Length[sol] > 0,
          CheckSingularRules[sol[[1]], singRule, assumps], {}];

        If[verbose,
          If[Length[singulars] > 0,
            Print["  WARNING: still singular in: ", First /@ singulars],
            Print["  OK \[LongDash] solution now regular."]
          ];
          Print[""];
        ];

        result["Solution"] = sol;
        result["PivotInfo"] = pivotResult;
        result["Singular"] = Length[singulars] > 0;
        result["SingularRules"] = singulars;
      ,
        If[verbose,
          Print["  System is nonlinear; automatic pivot selection not available."];
          Print[""];
        ];
      ];
    ,
      If[verbose,
        Print["  All rules regular."];
        Print[""];
      ];
    ];
  ];

  (* --- Step 5: Verification --- *)
  If[verbose, Print["=== Verification ==="]];
  Do[
    residuals = ((# == 0) & /@ eqs) /. sol[[i]] //
                Simplify[#, Assumptions -> assumps] &;
    allTrue = And @@ residuals;
    If[verbose,
      Print["  Branch ", i, ": ",
        If[TrueQ[allTrue], "VERIFIED", "CHECK NEEDED"]]
    ],
    {i, Length[sol]}
  ];
  If[verbose, Print[""]];

  (* --- Step 6: Display final solution --- *)
  If[verbose,
    Print["=== Solution ==="];
    Do[
      Print["--- Branch ", i, " ---"];
      Module[{branch = sol[[i]], determined, free},
        determined = First /@ branch;
        free = Complement[vars, determined];
        If[free =!= {}, Print["  Independent: ", free]];
        Do[Print["  ", rule], {rule, branch}];
      ];
      Print[""],
      {i, Length[sol]}
    ];
  ];

  result
];

End[];       (* `Private` *)
EndPackage[];
