BeginPackage["LINVARIANT`ISODISTORT`",{"LINVARIANT`Structure`","LINVARIANT`GroupTheory`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
ShowIsoModes                 ::usage = "ShowIsoModes[PosVec]"
GetIsoBasis                  ::usage = "GetIsoBasis[grp0, Wyckoff0]"
ISODISTORT                   ::usage = "ISODISTORT[R0, pos0, pos, IsoMatrix, label]"
ImposeMode                   ::usage = "ImposeMode[Wyckoff0, IsoMatrix, modeset, s]"
GetIRvector                  ::usage = "GetIRvector[id, pos]"
GetBasisField                ::usage = "GetBasisField[id, BasisMatrix, pos]"
GetOpMatrix                  ::usage = "GetOpMatrix[SymOpFile, pos, IsoMatrix, Modes]"
GetIsoVars                   ::usage = "GetIsoVars[IsoDispModes]"
GetIsoTransformRules         ::usage = "GetIsoDispTransformRules[OpMatrix, IsoDispModes, TranType]"
GetIsoStrainTransformRules   ::usage = "GetIsoStrainTransformRules[GridSymFile]"
Epsilon2Field                ::usage = "Epsilon2Field[strain]"
Field2Epsilon                ::usage = "Field2Epsilon[field]"
GetEpsilonijRule             ::usage = "GetEpsilonijRule[symfile]"
GetTransformationRules       ::usage = "GetTransformationRules[spg0, OpDispModes, OpDispMatrix, SpinModes, OpSpinMatrix]"
GetInvariants                ::usage = "GetInvariants[seeds, order, AllModes, OpMatrix, GridSymFile]"
ImposeDW                     ::usage = "ImposeDW[Wyckoff0, IsoMatrix, modeset, {Nx, Ny, Nz}]"
ImposeIsoStrainVariedDspMode ::usage = "ImposeIsoStrainVariedDspMode[Wyckoff0, IsoMatrix, modeset, LV]"
ShowInvariantTable           ::usage = "ShowInvariantTable[TableData]"
Jij                          ::usage = "Jij[r0, MeshDim]"
FieldCode2var                ::usage = "FieldCode2var[code, varstr]"
GetSiteInt                   ::usage = "GetSiteInt[spg0, field]"
GetSiteCluster               ::usage = "GetSiteCluster[spg0, PosVec]"
InvariantEOnSite             ::usage = "InvariantEOnSite[xyz, expr]"

(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
mIso
Iso
Epsilon
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)
(*Options[ImportIsodistortCIF]    = {Fractional->False, CorrectLabels->True, Tolerance->10^-6}*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
ShowIsoModes[PosVec_] := Module[{StructData},
  StructData = Table[{ElementData[i[[3]], "IconColor"], 
                      Sphere[i[[1]], QuantityMagnitude[ElementData[i[[3]], "AtomicRadius"], "Angstroms"]/2], 
                      Black, 
                      Text[i[[3]]<>ToString@Position[PosVec, i][[1, 1]], 
                      i[[1]]]}, {i, PosVec}];

  ArrowData = Table[{Green, 
                     Arrowheads[0.03], 
                     Arrow@Tube[{i[[1]], 
                     i[[1]] + i[[2]]}]}, {i, PosVec}];

  Print[Graphics3D[{StructData, ArrowData}, 
                   ImageSize -> 500, 
                   Axes -> True, 
                   AxesLabel -> (Style[#, Bold, 64]&/@{"a", "b", "c"}), 
                   ViewPoint -> {0, 0, \[Infinity]}]];
]

ISODISTORT[R0_, pos0_, pos_, IsoMatrix_, label_] := Module[{imode, Amp, NN},
  Amp = Table[Flatten[R0.PbcDiff[#] & /@ (Transpose[pos][[1]] - Transpose[pos0][[1]])].Normalize@Normal[IsoMatrix[[;; , imode]]], {imode, Length@label}];
  NN = Table[1/Norm@Flatten[R0.# & /@ Partition[Normal[IsoMatrix][[;; , imode]], 3]], {imode, Length@label}];
  Return[{Range[Length@label], label, NN, Chop[Amp, 2 10^-4]}\[Transpose]]
]

ImposeMode[Wyckoff0_, IsoMatrix_, modeset_, s_] := Module[{mode, id, Amp, pos},
  pos = Wyckoff0;
  Do[id = mode[[1]]; Amp = s mode[[3]] mode[[4]];
     pos = Transpose[{#[[1]] & /@ pos + Amp If[IntegerQ[id], # & /@ Partition[IsoMatrix[[;; , id]] // Normal, 3], Print["mode not exist!"]], First@StringCases[#,RegularExpression["[[:upper:]][[:lower:]]*"]] & /@ Transpose[pos][[2]]}], {mode, modeset}];
  Return[pos]
]

(*Jij[spg0_, OptionsPattern[{"Cluster" -> False}]] := Module[{tij, i, j}, 
  tij = DeleteDuplicates[DeleteCases[Flatten[Table[Total[Rationalize[#1[[2]]].{Subscript[Iso, 1, #1[[1]][[1]], #1[[1]][[2]], #1[[1]][[3]]], Subscript[Iso, 2, #1[[1]][[1]], #1[[1]][[2]], #1[[1]][[3]]], Subscript[Iso, 3, #1[[1]][[1]], #1[[1]][[2]], #1[[1]][[3]]]} Rationalize[#2[[2]]].{Subscript[Iso, 1, #2[[1]][[1]], #2[[1]][[2]], #2[[1]][[3]]], Subscript[Iso, 2, #2[[1]][[1]], #2[[1]][[2]], #2[[1]][[3]]], Subscript[Iso, 3, #2[[1]][[1]], #2[[1]][[2]], #2[[1]][[3]]]} & @@@
  ({Flatten[GetSiteCluster[spg0, {{{0, 0, 0}, v1}}, "spin" -> False], 1], Flatten[GetSiteCluster[spg0, {{neighbor, v2}}, "spin" -> False], 1]}\[Transpose])], {v1, IdentityMatrix[3]}, {v2, IdentityMatrix[3]}, {neighbor, Flatten[GetUpTo3rdNeighbors[], 1]}]], 0], #1 === -#2 || #1 === #2 &];
  Return[Table[Expand[tij[[j]]/(GCD @@ Table[tij[[j]][[i]] /. {Subscript[Iso, i_, x_, y_, z_] -> 1}, {i, Length@tij[[j]]}])], {j, Length@tij}]]
]*)

Jij[spg0_, fielddef_?ListQ, OptionsPattern[{"OnSite" -> False, "IsoTranRules"->{}}]] := Module[{tij, i, j, field, v0, v1, neighbor},
  tij = Expand[DeleteDuplicates[DeleteCases[Flatten[
           Table[
                 field = {{fielddef[[1,1]], {{0, 0, 0}, v0}, fielddef[[1,3]]}, {fielddef[[2,1]], {neighbor, v1}, fielddef[[2,3]]}};
                 GetSiteInt[spg0, field, "IsoTranRules"->OptionValue["IsoTranRules"]], {v0, fielddef[[1,2]]}, {v1, fielddef[[2,2]]}, {neighbor, Flatten[GetUpTo3rdNeighbors["OnSite" -> OptionValue["OnSite"]], 1]}
                 ]
        ], 0], #1 === -#2 || #1 === #2 &]];
  Return[Expand[#/GCD @@ (If[MatchQ[#, Plus[_, __]], Level[#, {1}], Level[#, {0}]] /. Thread[Variables[#] -> ConstantArray[1, Length[Variables[#]]]])] & /@ tij]
]

Epsilon2Field[strain_] := Module[{}, 
  {{{0, 0, 0}, IdentityMatrix[3][[strain[[2]]]]}, {IdentityMatrix[3][[strain[[3]]]], IdentityMatrix[3][[strain[[2]]]]}, {{0, 0, 0}, IdentityMatrix[3][[strain[[3]]]]}, {IdentityMatrix[3][[strain[[2]]]], IdentityMatrix[3][[strain[[3]]]]}}]

Field2Epsilon[field_] := Module[{}, 
  1/2 (Sum[Normalize[field[[2]][[1]]][[j]] Normalize[field[[1]][[2]]][[i]] Subscript[Epsilon, i, j], {j, 3}, {i, 3}] + Sum[Normalize[field[[4]][[1]]][[j]] Normalize[field[[3]][[2]]][[i]] Subscript[Epsilon, i, j], {j, 3}, {i, 3}])/.{Subscript[Epsilon, 2, 1] -> Subscript[Epsilon, 1, 2], Subscript[Epsilon, 3, 1] -> Subscript[Epsilon, 1, 3], Subscript[Epsilon, 3, 2] -> Subscript[Epsilon, 2, 3]}]

(*Epsilon2Field[strain_] := Module[{},
  {{{0, 0, 0}, IdentityMatrix[3][[strain[[2]]]]}, {IdentityMatrix[3][[strain[[3]]]], IdentityMatrix[3][[strain[[2]]]]}}
]

Field2Epsilon[field_] := Module[{},
  Sum[Normalize[field[[2]][[1]]][[j]] Normalize[field[[1]][[2]]][[i]] Subscript[Epsilon, i, j], {j, 3}, {i, 3}] /. {Subscript[Epsilon, 2, 1] -> Subscript[Epsilon, 1, 2], Subscript[Epsilon, 3, 1] -> Subscript[Epsilon, 1, 3], Subscript[Epsilon, 3, 2] -> Subscript[Epsilon, 2, 3]}
]*)

GetEpsilonijRule[spg0_] := Module[{tij, i},
  strains = DeleteDuplicates@Flatten[SparseArray[{{i_, j_} /; i == j -> Subscript[Epsilon, i, j], {i_, j_} /; i < j -> Subscript[Epsilon, i, j], {i_, j_} /; i > j -> Subscript[Epsilon, j, i]}, {3, 3}] // Normal];
  Thread[strains -> #] & /@ (Table[Field2Epsilon[#] & /@ GetSiteCluster[spg0, Epsilon2Field[ep]], {ep, strains}]\[Transpose])
]

GetIsoVars[IsoDispModes_] := Module[{VarString, Var},
  VarString = {#1 & @@@ IsoDispModes, StringReplace[First@StringCases[#2, RegularExpression["[[:upper:]]*\\d[+-]"]], {"+" -> "Plus", "-" -> "Minus"}] & @@@ IsoDispModes, StringPart[#2, -2] & @@@ IsoDispModes}\[Transpose];
  Var = {#1, Subscript[ToExpression[#2], ToExpression[#3]]} & @@@ VarString;
  Return[Var]
]

GetIRvector[id_, IsoMatrix_, pos_] := Module[{IRvector},
  IRvector = {If[IntegerQ[id], # & /@ Partition[IsoMatrix[[;; , id]] // Normal, 3], Print["mode not exist!"]], # & /@ (pos\[Transpose][[2]])}\[Transpose];
  Return[IRvector]
]

GetBasisField[id_, BasisMatrix_, pos_, OptionsPattern["spin" -> False]] := Module[{BasisField},
  BasisField = Which[ListQ[id], GetBasisField[#, BasisMatrix, pos] & /@ id, IntegerQ[id], {pos\[Transpose][[2]], {pos\[Transpose][[1]], # & /@ Partition[BasisMatrix[[;; , id]], 3]}\[Transpose], ConstantArray[If[OptionValue["spin"],1,2], Length[pos]]}\[Transpose]];
  Return[BasisField]
]

GetOpMatrix[grp_, pos_, IsoMatrix_, Modes_, OptionsPattern["spin" -> False]] := Module[{op2, ir1, ir2, mat, OpMatrix, AllIRvectors, AllTransformedIRvectors, IsoDim},

  IsoDim = Length@Modes;
  OpMatrix = If[IsoDim != 0,
  AllIRvectors = Flatten[GetIRvector[#1, IsoMatrix, pos]\[Transpose][[1]]] & @@@ Modes;
  AllTransformedIRvectors = Transpose[ParallelTable[SymmetryOpVectorField[grp, pos, GetIRvector[id, IsoMatrix, pos], "spin" -> OptionValue["spin"]], {id, IsoDim}, DistributedContexts -> {"LINVARIANT`ISODISTORT`Private`"}]];
  ParallelTable[mat=Table[Flatten[op2[[ir1]]\[Transpose][[1]]].AllIRvectors[[ir2]], {ir1, Range@IsoDim}, {ir2, Range@IsoDim}]; SparseArray[Normalize[#] & /@ mat], {op2, AllTransformedIRvectors}, DistributedContexts -> {"LINVARIANT`ISODISTORT`Private`"}],
  Table[{}, {Length@grp}]
  ];
  Return[OpMatrix]
]

GetOpMatrixHeavy[grp_, pos_, BasisMatrix_, BasisLabels_, OptionsPattern["spin" -> False]] := Module[{op, op2, ir1, ir2, mat, OpMatrix, AllIRvectors, BasisField, TransformedBasisField, BasisDim},
  Which[AssociationQ[grp], GetOpMatrix[#, pos, BasisMatrix, BasisLabels] & /@ Keys[grp],
   ListQ[grp] && ! MatrixQ[grp], GetOpMatrix[#, pos, BasisMatrix, BasisLabels] & /@ grp,
   MatrixQ[grp], GetOpMatrix[M42xyzStr[grp], pos, BasisMatrix, BasisLabels],
   StringQ[grp],
   BasisDim = Length@BasisLabels;
   OpMatrix = If[BasisDim != 0,
     BasisField = GetBasisField[Range[BasisDim], BasisMatrix, pos, "spin" -> OptionValue["spin"]];
     TransformedBasisField = ParallelTable[SymmetryOpBasisField[grp, field], {field, BasisField}];
     mat = ParallelTable[Flatten[TransformedBasisField[[i]]\[Transpose][[2]]\[Transpose][[2]]].Flatten[BasisField[[j]]\[Transpose][[2]]\[Transpose][[2]]], {i, Range@BasisDim}, {j, Range@BasisDim}]; 
     SparseArray[Normalize[#] & /@ mat], Table[{}, {Length@grp}]]; 
   Return[OpMatrix]]
]

GetIsoTransformRules[OpMatrix_, IsoModes_, TranType_] := Module[{IsoDim, IsoVars, SpinDispRules, rules, i, var},
  IsoDim = Length@IsoModes;
  rules = If[IsoDim != 0,
  IsoVars = Which[TranType == "disp", Subscript[Iso, ToExpression[#1]] &@@@ IsoModes, TranType == "spin", Subscript[mIso, ToExpression[#1]] &@@@ IsoModes];
  SpinDispRules = IsoVars[[#1[[1]]]] -> #2 IsoVars[[#1[[2]]]] & @@@ Drop[ArrayRules[OpMatrix], -1];
  Table[First@DeleteDuplicates[Keys[SpinDispRules[[#]]] & /@ i] -> Total[Values[SpinDispRules[[#]]] & /@ i], {i, Table[Flatten[Position[Keys@SpinDispRules, var]], {var, Which[TranType == "disp", Subscript[Iso, #] & /@ Range[IsoDim], TranType == "spin", Subscript[mIso, #] & /@ Range[IsoDim]]}]}], {}];
  Return[rules]
]

GetIsoStrainTransformRules[spg0_] := Module[{StrainRules},
  StrainRules = GetEpsilonijRule[spg0];
  Return[StrainRules]
]

GetTransformationRules[spg0_, OpDispModes_, OpDispMatrix_, SpinModes_, OpSpinMatrix_] := Module[{},
  Join[#1, #2, #3] & @@@ ({GetIsoTransformRules[#, OpDispModes, "disp"] & /@ OpDispMatrix, 
    GetIsoTransformRules[#, OpSpinModes, "spin"] & /@ OpSpinMatrix, 
    GetIsoStrainTransformRules[spg0]}\[Transpose])
]

GetInvariants[seeds_, order_, DispModes_, OpDispMatrix_, SpinModes_, OpSpinMatrix_, spg0_, OptionsPattern[{"strain"->False}]] := Module[{monomials, invariant, TransformRules, i, j, ss, strains},
  strains = DeleteDuplicates@Flatten[SparseArray[{{i_, j_} /; i == j -> Subscript[Epsilon, i, j], {i_, j_} /; i < j -> Subscript[Epsilon, i, j], {i_, j_} /; i > j -> Subscript[Epsilon, j, i]}, {3, 3}] // Normal];
  CoeffNorm = Thread[seeds -> ConstantArray[1, Length@seeds]];
  monomials = MonomialList[Total[seeds]^order];
  monomials = If[OptionValue["strain"], Flatten[Table[# ss & /@ monomials, {ss, strains}]], monomials];
  TransformRules = Dispatch[Join[#1, #2, #3] &@@@ ({GetIsoTransformRules[#, DispModes, "disp"] & /@ OpDispMatrix, GetIsoTransformRules[#, SpinModes, "spin"] & /@ OpSpinMatrix, GetIsoStrainTransformRules[spg0]}\[Transpose])];
  
  invariant = Rationalize[DeleteDuplicates[DeleteCases[Union[Expand[Total[(monomials /. # & /@ TransformRules)]]], i_/;i==0], (#1 -#2 == 0 || #1 + #2 == 0) &]];

  Return[Expand[#/GCD @@ (If[MatchQ[#, Plus[_, __]], Level[#, {1}], Level[#, {0}]] /. Thread[Variables[#] -> ConstantArray[1, Length[Variables[#]]]])] &/@ invariant]

  (*Return[Table[invariant[[i]]/(GCD @@ (Table[invariant[[i]][[j]], {j, Length[invariant[[i]]]}] /. Thread[Variables[invariant[[i]]] -> ConstantArray[1, Length[Variables[invariant[[i]]]]]])), {i, Length@invariant}]]*)
]

ImposeDW[Wyckoff0_, IsoMatrix_, modeset_, {Nx_, Ny_, Nz_}] := Module[{mode, id, Amp, pos, s, ix, iy, iz, Superpos},
  Superpos = Table[{#1 + {ix, iy, iz}, #2} & @@@ Wyckoff0, {ix, 0, Nx - 1}, {iy, 0, Ny - 1}, {iz, 0, Nz - 1}];
  Do[pos = Superpos[[ix]][[iy]][[iz]];
     s = Cos[2 Pi {1/Nx, 1/Ny, 1/Nz}.{ix, iy, iz}];
     Do[id = mode[[1]]; Amp = s mode[[3]] mode[[4]];
        pos = Transpose[{#[[1]] & /@ pos + Amp If[IntegerQ[id], # & /@ Partition[IsoMatrix[[;; , id]] // Normal, 3], Print["mode not exist!"]], First@StringCases[#, RegularExpression["[[:upper:]][[:lower:]]*"]] & /@ Transpose[pos][[2]]}], {mode, modeset}
        ];
     Superpos[[ix]][[iy]][[iz]] = pos, {ix, Nx}, {iy, Ny}, {iz, Nz}
     ];
  Return[Superpos]
]

ImposeIsoStrainVariedDspMode[Wyckoff0_, IsoMatrix_, modeset_, LV_] := Module[{mode, modesetnew, id, Amp, pos, NN},
  pos = Wyckoff0;
  Do[id = mode[[1]];
     NN = 1/Norm@Flatten[LV.# & /@ Partition[Normal[IsoMatrix][[;; , id]], 3]];
     Amp = NN mode[[4]];
     pos = Transpose[{#[[1]] & /@ pos + Amp If[IntegerQ[id] , # & /@ Partition[IsoMatrix[[;; , id]] // Normal, 3], Print["mode not exist!"]],
                      First@StringCases[#, RegularExpression["[[:upper:]][[:lower:]]*"]] & /@ Transpose[pos][[2]]}],
   {mode, modeset}];
  Return[pos]
]

ShowInvariantTable[TableData_, param_, OptionsPattern["FontSize" -> 12]] := Module[{m, n},
  Print[Rotate[Grid[Table[Style[Rotate[# // Expand, -270 Degree], Black, Bold, OptionValue["FontSize"]] & /@ (Flatten[Table[{Flatten[{"param", param[[n]]}], Prepend[TableData[[n]], n]}, {n, Length@TableData}], 1][[m]]), {m, 2 Length@TableData}], Alignment -> Left, Frame -> All], 270 Degree]]]

GetIsoBasis[grp0_, Wyckoff0_, kpoint_: {{0, 0, 0}, "\[CapitalGamma]"}, ct0_: {}, OptionsPattern[{"spin" -> False}]] := Module[{w, p, imode, WyckoffSites, ZeroModeMatrix, SymmetryAdaptedBasis, Sites, IsoDispModes, IsoDispModeMatrix, OpDispMatrix, ct, ProjMat, g, basis, grpk, m, n, lp},
  grpk = GetGroupK[grp0, kpoint[[1]]];
  ct = If[Length[ct0] == 0, GetCharacters[grpk, "kpoint" -> kpoint[[2]]], ct0];
  g = Length[grpk];
  WyckoffSites = GetWyckoffImages[grp0, {#}] & /@ Wyckoff0;
  SymmetryAdaptedBasis = Table[
    Sites = Map[{#[[1]], #[[2]]} &, WyckoffSites[[w]]];
    IsoDispModes = MapIndexed[{First[#2], #1, 1, 0} &, Flatten[Table[Subscript[#[[2]], xyz], {xyz, 3}] & /@ Sites]];
    IsoDispModeMatrix = IdentityMatrix[3 Length[Sites]];
    OpDispMatrix = GetOpMatrix[grpk, Sites, IsoDispModeMatrix, IsoDispModes, "spin" -> OptionValue["spin"]];
    ProjMat = Table[lp = ct[[2]][[p]][[1]]; Sum[{m, n} = First@Position[ct[[1]], Keys[grpk][[i]]]; lp/g Conjugate[ct[[2]][[p, m]]] Rationalize@Normal[OpDispMatrix[[i]]], {i, g}], {p, Length[ct[[2]]]}];
    basis = Table[Complement[Orthogonalize[ProjMat[[p]].# & /@ (IsoDispModeMatrix\[Transpose])], {Table[0, {i, 3 Length[Sites]}]}], {p, Length[ct[[2]]]}];
    If[#2 == {}, Unevaluated[Sequence[]], {Table[Sites[[1]][[2]] <> " " <> #1 <> "-" <> ToString[Length@#2] <> "(" <> ToString[imode] <> ")", {imode, Length@#2}], #2}] & @@@ Thread[ct[[3]] -> basis], {w, Length@WyckoffSites}];
  IsoDispModes = MapIndexed[{First[#2], #1, 1, 0.} &, Flatten[Flatten[SymmetryAdaptedBasis, 1]\[Transpose][[1]]]];
  IsoDispModeMatrix = # & /@ Transpose[Fold[ArrayFlatten[{{#1, 0}, {0, #2}}] &, Flatten[#, 1] & /@ (#\[Transpose][[2]] & /@ SymmetryAdaptedBasis)]];
  Return[{IsoDispModes, IsoDispModeMatrix}]
]

FieldCode2var[code_, varstr_, OptionsPattern[{"rec"->False}]] := Module[{x1, x2, x3},
  {x1, x2, x3} = If[OptionValue["rec"], Mod[code[[1]], 1], code[[1]]];
  If[ListQ[code[[2]]], 
     Sign[code[[2]].{1, 2, 3}] Subscript[ToExpression[varstr], Abs[code[[2]].{1, 2, 3}], x1, x2, x3],
     Which[MatchQ[code[[2]], _Subscript], Construct@@Join[{Head[code[[2]]]}, Level[code[[2]], {1}], code[[1]]], 
           MatchQ[code[[2]], _Times], Level[code[[2]], {1}][[1]] Construct@@Join[{Subscript}, Level[code[[2]], {2}], code[[1]]],
           True, Print["expression not in the right form!"]
       ]
     ]
]

GetSiteCluster[spg0_, PosVec_, OptionsPattern[{"spin"->False, "IsoTranRules"->{}}]] := Block[{xyzStrData, xyzRotTranData, xyzTranslation, xyzRotData, field, newpos, newvec, difftable, diff, posmap, i, j, det},
  xyzStrData = xyzStrData = Keys[spg0];
  xyzRotTranData = Table[ToExpression["{" <> xyzStrData[[i]] <> "}"], {i, Length[xyzStrData]}];
  xyzTranslation = xyzRotTranData /. {ToExpression["x"] -> 0, ToExpression["y"] -> 0, ToExpression["z"] -> 0};
  xyzRotData = xyzRotTranData - xyzTranslation;

  If[OptionValue["IsoTranRules"]=={},
    det = If[OptionValue["spin"], 1, 2];
    newvec = Table[Det[xyz2Rot[op]]^det op /. Thread[ToExpression[{"x", "y", "z"}] -> #] & /@ (PosVec\[Transpose][[2]]), {op, xyzRotData}];
    newpos = Table[op /. Thread[ToExpression[{"x", "y", "z"}] -> #] & /@ (PosVec\[Transpose][[1]]), {op, xyzRotData}];,
    newvec = Table[PosVec\[Transpose][[2]] /. op, {op, OptionValue["IsoTranRules"]}];
    newpos = Table[op /. Thread[ToExpression[{"x", "y", "z"}] -> #] & /@ (PosVec\[Transpose][[1]]), {op, xyzRotData}];
  ];
 
  Return[Rationalize[{#[[1]], #[[2]]}]\[Transpose] &/@ ({newpos, newvec}\[Transpose])]
]

GetSiteInt[spg0_, field_, OptionsPattern[{"IsoTranRules"->{}}]] := Module[{tij, factor, OpMatrix},
  tij = Total@Fold[Times, Table[FieldCode2var[#, posvec[[1]]] & /@ Flatten[GetSiteCluster[spg0, {posvec[[2]]}, "spin" -> posvec[[3]], "IsoTranRules"->OptionValue["IsoTranRules"]], 1], {posvec, field}]];
  factor = GCD @@ (If[MatchQ[tij, Plus[_, __]], Level[tij, {1}], Level[tij, {0}]] /. Thread[Variables[tij] -> ConstantArray[1, Length[Variables[tij]]]]);
  Return[If[factor == 0, 0, Expand[tij/factor]]]
] 

InvariantEOnSite[expr_] := Module[{X, i, dx, dy, dz},
  expr /. {Subscript[X_, i_, dx_, dy_, dz_] :> Subscript[X, i, ToExpression["ix"] + dx, ToExpression["iy"] + dy, ToExpression["iz"] + dz]}
]
(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
