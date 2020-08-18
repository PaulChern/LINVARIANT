BeginPackage["LINVARIANT`INVARIANT`",{"LINVARIANT`Structure`","LINVARIANT`GroupTheory`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
ShowIsoModes                 ::usage = "ShowIsoModes[PosVec]"
GetIsoBasis                  ::usage = "GetIsoBasis[grp0, Wyckoff0]"
GetSymAdaptedBasis           ::usage = "GetSymAdaptedBasis[grp0, pos, kpoint, ftype, ct0]"
ISODISTORT                   ::usage = "ISODISTORT[R0, pos0, pos, IsoMatrix, label]"
ImposeMode                   ::usage = "ImposeMode[Wyckoff0, IsoMatrix, modeset, s]"
GetIRvector                  ::usage = "GetIRvector[id, pos]"
GetBasisField                ::usage = "GetBasisField[id, BasisMatrix, pos]"
GetOpMatrix                  ::usage = "GetOpMatrix[SymOpFile, pos, IsoMatrix, Modes]"
GetMatrixRep                 ::usage = "GetMatrixRep[SymOpFile, pos, IsoMatrix, Modes]"
SymmetryOpBasisField         ::usage = "SymmetryOpBasisField[grp, BasisField]"
MonomialTransform            ::usage = "MonomialTransform[monomials, rules]"
GetIsoVars                   ::usage = "GetIsoVars[IsoDispModes]"
GetIsoTransformRules         ::usage = "GetIsoDispTransformRules[OpMatrix, IsoDispModes, TranType]"
GetIsoStrainTransformRules   ::usage = "GetIsoStrainTransformRules[GridSymFile]"
Epsilon2Field                ::usage = "Epsilon2Field[strain]"
Field2Epsilon                ::usage = "Field2Epsilon[field]"
GetEpsilonijRule             ::usage = "GetEpsilonijRule[symfile]"
GetTransformationRules       ::usage = "GetTransformationRules[spg0, OpMatrix]"
GetInvariants                ::usage = "GetInvariants[seeds, order, AllModes, OpMatrix, GridSymFile]"
NumberCommonDivisor          ::usage = "NumberCommonDivisor[NumList]"
GetConstantFactor            ::usage = "GetConstantFactor[expr]"
SimplifyCommonFactor         ::usage = "SimplifyCommonFactor[expr]"
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

ISODISTORT[R0_, pos0_, pos_, IsoMatrix_, label_] := Module[{imode, Amp, NN, posmatched},
  posmatched = Transpose[{PosMatchTo[pos0\[Transpose][[1]], pos\[Transpose][[1]], 0.01][[2]], Transpose[pos0][[2]]}];
  Amp = Table[Flatten[R0.PbcDiff[#] & /@ (Transpose[posmatched][[1]] - Transpose[pos0][[1]])].Normalize@Normal[IsoMatrix[[;; , imode]]], {imode, Length@label}];
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
  Which[Length@Dimensions[fielddef] == 2,
        tij = Expand[DeleteDuplicates[DeleteCases[Flatten[
           Table[
                 field = {{fielddef[[1,1]], {{0, 0, 0}, v0}, fielddef[[1,3]]}, {fielddef[[2,1]], {neighbor, v1}, fielddef[[2,3]]}};
                 GetSiteInt[spg0, field, "IsoTranRules"->OptionValue["IsoTranRules"]], {v0, fielddef[[1,2]]}, {v1, fielddef[[2,2]]}, {neighbor, Flatten[GetUpTo3rdNeighbors["OnSite" -> OptionValue["OnSite"]], 1]}
                 ]
        ], 0], #1 === -#2 || #1 === #2 &]];
        Return[Expand[#/GCD @@ (If[MatchQ[#, Plus[_, __]], Level[#, {1}], Level[#, {0}]] /. Thread[Variables[#] -> ConstantArray[1, Length[Variables[#]]]])] & /@ tij],
        Length@Dimensions[fielddef] > 2,
        Flatten[Jij[spg0, #, "OnSite" -> OptionValue["OnSite"], "IsoTranRules"->OptionValue["IsoTranRules"]] &/@ fielddef]
      ]
]

Epsilon2Field[strain_] := Module[{}, 
  {{{0, 0, 0}, IdentityMatrix[3][[strain[[2]]]]}, {IdentityMatrix[3][[strain[[3]]]], IdentityMatrix[3][[strain[[2]]]]}, {{0, 0, 0}, IdentityMatrix[3][[strain[[3]]]]}, {IdentityMatrix[3][[strain[[2]]]], IdentityMatrix[3][[strain[[3]]]]}}
]

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
  Thread[strains -> #] & /@ (Table[Field2Epsilon[#] & /@ GetSiteCluster[spg0, Epsilon2Field[ep], "disp"], {ep, strains}]\[Transpose])
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

GetBasisField[id_, BasisMatrix_, BasisLabels_, pos_, ftype_] := Module[{BasisLabelsupdn, BasisField, latt, sites, i, j, lm, p, p0, dim},
  {latt, sites} = pos;
  BasisField = Which[ListQ[id], 
                     GetBasisField[#, BasisMatrix, BasisLabels, pos, ftype] & /@ id, 
                     IntegerQ[id], 
                     If[ftype=="disp"||ftype=="spin",
                        {sites\[Transpose][[2]], 
                        {sites\[Transpose][[1]], Partition[Normal[BasisMatrix][[;; , id]], 3]}\[Transpose], 
                        ConstantArray[ftype, Length[sites]]}\[Transpose],
                        (*ftype=="spinor"||"orbital"*)
                        BasisLabelsupdn = Join[BasisLabels, BasisLabels];
                        sites = Join[sites, sites];
                        dim = Total@Flatten[Table[2 # + 1 &/@ lm, {lm, BasisLabelsupdn}]];
                        If[dim != Length[BasisMatrix], 
                           Print["Error: Basis Matrix dimension mismatch!!!"];Abort[]];
                        p0 = 1;
                        Table[{sites[[i,2]],
                               {sites[[i,1]],
                               Table[p = p0; p0 = p0 + 2*lm +1;
                                     BasisMatrix[[p;;p+2*lm,id]], {lm, BasisLabelsupdn[[i]]}]}, 
                               If[i<=Length[sites]/2, "up", "dn"]}, {i, Length@BasisLabelsupdn}]
                
                        ]
                    ];
  Return[BasisField]
]

GetOpMatrix[grp_, pos_, IsoMatrix_, Modes_, ftype_] := Module[{latt, sites, op2, ir1, ir2, mat, OpMatrix, AllIRvectors, AllTransformedIRvectors, NormFactor, IsoDim, i, vec},
  {latt, sites} = pos;
  IsoDim = Length@Modes;

  Which[AssociationQ[grp],
        If[IsoDim != 0,
        {NormFactor, AllIRvectors} = Table[
                    vec=GetIRvector[i, IsoMatrix, sites]\[Transpose][[1]];
                    {1/Norm[Flatten[vec]], vec}, {i, IsoDim}]\[Transpose];
        AllTransformedIRvectors = Transpose[ParallelTable[SymmetryOpVectorField[grp, sites, GetIRvector[id, IsoMatrix, sites], ftype], {id, IsoDim}, DistributedContexts -> {"LINVARIANT`INVARIANT`Private`"}]];
        ParallelTable[
          mat=Table[Rationalize[NormFactor[[ir1]]*NormFactor[[ir2]]*Flatten[op2[[ir1]]\[Transpose][[1]]].Flatten[AllIRvectors[[ir2]]]], {ir1, Range@IsoDim}, {ir2, Range@IsoDim}]; 
          SparseArray[mat], {op2, AllTransformedIRvectors}, DistributedContexts -> {"LINVARIANT`INVARIANT`Private`"}], Table[{{}}, {Length@grp}]],
        ListQ[grp],
        OpMatrix = GetOpMatrix[#, pos, IsoMatrix, Modes, ftype] &/@ grp;
        If[IsoDim != 0,
           Fold[Dot, #] &/@ Tuples[OpMatrix],
           Table[{{}}, Fold[Times, Length[#] &/@ OpMatrix]]]
  ]
]

GetMatrixRep[grp0_, grpt_, pos_, BasisMatrix_, BasisLabels_, ftype_] := Module[{i, j, grp, op, op2, ir1, ir2, mat, OpMatrix, AllIRvectors, BasisField, TransformedBasisField, BasisDim, g, ig, field, NormFactor, latt, sites, cell, basis},
  grp = {grpt, grp0};
  {latt, sites} = pos;
  BasisDim = Length@BasisMatrix;
  cell = GetCellFromGrp[grpt];
  OpMatrix =
    Table[
      If[BasisDim != 0,
         {NormFactor, BasisField} = Table[
         basis=GetBasisField[i, BasisMatrix, BasisLabels, pos, ftype];
         {1/Norm[Flatten[basis\[Transpose][[2]]\[Transpose][[2]]]], basis}, {i,BasisDim}]\[Transpose];
         TransformedBasisField = ParallelTable[SymmetryOpBasisField[g, pos, cell, #, ftype] &/@ BasisField, {g, Keys[grp[[ig]]]}, DistributedContexts -> {"LINVARIANT`INVARIANT`Private`"}];
         ParallelTable[mat = Table[Rationalize[NormFactor[[i]]*NormFactor[[j]]*Flatten[TransformedBasisField[[ig]][[i]]\[Transpose][[2]]\[Transpose][[2]]].Flatten[BasisField[[j]]\[Transpose][[2]]\[Transpose][[2]]]], {i, Range@BasisDim}, {j, Range@BasisDim}];
         SparseArray[mat], {ig, Length[grp[[ig]]]}, DistributedContexts -> {"LINVARIANT`INVARIANT`Private`"}],
         Table[{}, {Length[grp[[ig]]]}]], {ig, 2}];
  If[BasisDim != 0,
     Flatten[Table[OpMatrix[[1,i]].OpMatrix[[2,j]], {i, Length[grpt]}, {j, Length[grp0]}], 1],
     (*Fold[Dot, #] &/@ Tuples[OpMatrix],*)
     Table[{{}}, Fold[Times, Length[#] &/@ OpMatrix]]
    ]
]

MonomialTransform[monomials_, rules_] := Module[{P},
  P = (monomials /. #) & /@ rules;
  Return[P]
]

SymmetryOpBasisField[grp_, pos_, cell_, BasisField_, ftype_] := Module[{latt, sites, site, i, trans, rot, rotL, su2, NewField, NewFieldup, NewFielddn, posvec, difftable, posmap, updn, basis, upbasis, dnbasis, tesseral, dim, dgQ},
  Which[
   AssociationQ[grp],
   SymmetryOpBasisField[#, pos, cell, BasisField, ftype] & /@ Keys[grp],
   ListQ[grp],
   SymmetryOpBasisField[#, pos, cell, BasisField, ftype] & /@ grp,
   StringQ[grp] && Length[Dimensions[BasisField]] == 3,
   SymmetryOpBasisField[grp, pos, cell, #, ftype] & /@ BasisField,
   StringQ[grp] && Length[Dimensions[BasisField]] == 2,
   {latt, sites} = pos;
   {rot, trans, su2} = xyz2RotTSU2[grp];
   If[ftype=="disp"||ftype=="spin",
     Which[
       ftype=="disp",
       posvec={ModCell[N[rot.(#2[[1]])+trans], cell], N[rot.(#2[[2]])]} & @@@ BasisField,
       ftype=="spin",
       posvec={ModCell[N[rot.(#2[[1]])+trans], cell], Det[rot] N[rot.(#2[[2]])]} & @@@ BasisField
       ];
     difftable = DistMatrix[BasisField\[Transpose][[2]]\[Transpose][[1]], posvec\[Transpose][[1]], cell];
     posmap = Position[difftable, x_ /; TrueQ[Chop[x] == 0]];
     NewField = {BasisField[[#1]][[1]], posvec[[#2]], BasisField[[#1]][[3]]} & @@@ posmap,
     (*ftype=="spinor"||"orbital"*)
     dim = Length[BasisField];
     basis = Transpose[Transpose[#][[2]]] &/@ Transpose[Partition[BasisField, dim/2]];
     {upbasis, dnbasis} = Which[
       ftype=="orbital",
       Table[
         site = ModCell[N[rot.#+trans], cell] &/@ basis[[i,1]];
         tesseral = Partition[GetAngularMomentumRep[latt, grp, (Length[#]-1)/2, "Tesseral"].# &/@Flatten[basis[[i,2]], 1], Length@Flatten[basis[[i,2]], 1]/2];
         Transpose[{site, tesseral, {"up","dn"}}], {i, dim/2}]\[Transpose],
       ftype=="spinor",
       If[su2 === Null, Print["Error: For spinor, you need a double group!"];Abort[]];
       Table[
         site = ModCell[N[rot.#+trans], cell] &/@ basis[[i,1]];
         tesseral = Partition[GetAngularMomentumRep[latt, grp, (Length[#]-1)/2, "Tesseral"].# &/@Flatten[basis[[i,2]], 1], Length@Flatten[basis[[i,2]], 1]/2];
         Transpose[{site, Transpose[Transpose[Map[su2.# &, Transpose[#]]] &/@ Transpose[tesseral]], {"up", "dn"}}], {i, dim/2}]\[Transpose]        
     ];

     difftable = DistMatrix[sites\[Transpose][[1]], upbasis\[Transpose][[1]], cell];
     posmap = Position[difftable, x_ /; TrueQ[Chop[x] == 0]];
     NewFieldup = {sites[[#1]][[2]], upbasis[[#2]][[1;;2]], "up"} & @@@ posmap;
     NewFielddn = {sites[[#1]][[2]], dnbasis[[#2]][[1;;2]], "dn"} & @@@ posmap;
     NewField = Join[NewFieldup, NewFielddn];
   ];
  Return[NewField]]
]


GetIsoTransformRules[OpMatrix_, TranType_] := Module[{IsoDim, IsoVars, VarRules, rules, i, var},

  Which[Length[Dimensions@OpMatrix]==2, 
   IsoDim = Length[OpMatrix];
   rules = If[OpMatrix != {{}},
   IsoVars = Which[TranType == "disp", 
                   Subscript[ToExpression["Iso"], #] &/@ Range[IsoDim], 
                   TranType == "spin", 
                   Subscript[ToExpression["mIso"], #] &/@ Range[IsoDim],
                   TranType == "spinor" || TranType == "orbital",
                   Subscript[ToExpression["eIso"], #] &/@ Range[IsoDim]];

   VarRules = IsoVars[[#1[[1]]]] -> #2 IsoVars[[#1[[2]]]] & @@@ Drop[ArrayRules[OpMatrix], -1];

   Table[First@DeleteDuplicates[Keys[VarRules[[#]]] & /@ i] -> Total[Values[VarRules[[#]]] & /@ i], {i, Table[Flatten[Position[Keys@VarRules, var]], {var, IsoVars}]}], {}];
   Return[rules],
   Length[Dimensions@OpMatrix]==3,
   GetIsoTransformRules[#, TranType] &/@ OpMatrix]
]

GetIsoStrainTransformRules[spg0_] := Module[{StrainRules},
  Which[AssociationQ[spg0],
        GetEpsilonijRule[spg0],
        ListQ[spg0],
        GetEpsilonijRule[xyz2Grp[Flatten[GTimes[spg0]], "fast"->True]]
  ]
]

GetTransformationRules[spg0_, OpMatrix_] := Module[{OpDispMatrix, OpSpinMatrix, OpOrbitalMatrix},
  {OpDispMatrix, OpSpinMatrix, OpOrbitalMatrix} = OpMatrix;
  Join[#1, #2, #3, #4] & @@@ ({GetIsoTransformRules[OpDispMatrix, "disp"], 
    GetIsoTransformRules[OpSpinMatrix, "spin"], GetIsoTransformRules[OpOrbitalMatrix, "spinor"],
    GetIsoStrainTransformRules[spg0]}\[Transpose])
]

GetInvariants[seeds_, order_, OpMatrix_, spg0_, OptionsPattern[{"MustInclude"->{}}]] := Module[{fixlist, monomials, invariant, TransformRules, n, i, j, ss, factor, out, factorLCM, factorGCD, OpDispMatrix, OpSpinMatrix, OpOrbitalMatrix},
  {OpDispMatrix, OpSpinMatrix, OpOrbitalMatrix} = OpMatrix;
  out = Table[
  monomials = If[
    OptionValue["MustInclude"]=={}, 
    MonomialList[Total[seeds]^n],
    fixlist = Flatten[Table[MonomialList[Total[#1]^i], {i, #2}] &@@@ OptionValue["MustInclude"]];
    Flatten[Table[# ss & /@ MonomialList[Total[seeds]^n], {ss, fixlist}]]];
  TransformRules = GetTransformationRules[spg0, OpMatrix];
  invariant = Rationalize[DeleteDuplicates[DeleteCases[Union[Expand[Total[MonomialTransform[monomials, TransformRules]]]], i_/;i==0], (#1 -#2 == 0 || #1 + #2 == 0) &]];
  SimplifyCommonFactor[invariant], {n, order}];
  Return[DeleteDuplicates[#, (#1 - #2 == 0 || #1 + #2 == 0 || Expand[#1 + I #2] == 0 || Expand[#1 - I #2] == 0) &] &/@ out]
]

NumberCommonDivisor[NumList_] := Module[{TempList, DenominatorLCM},
 TempList = Which[Head[#] === Integer, #, Head[#] === Times, First@Level[#, {1}], Head[#] === Power, 1, Head[#] === Rational, #, (Head[#] === Complex)&&(Re[#]!=0), Abs[#], (Head[#] === Complex)&&(Re[#]==0), #] &/@ NumList;
 DenominatorLCM = If[MemberQ[TempList, _Rational], LCM @@ (Denominator[#] & /@ Extract[TempList, Position[TempList, _Rational]]), 1];
 Return[{DenominatorLCM, GCD @@ (TempList DenominatorLCM)}]
]

GetConstantFactor[expr_] := Module[{},
  If[ListQ[expr],  GetConstantFactor[#] & /@ expr, Return[(If[MatchQ[Expand@expr, Plus[_, __]], Level[Expand@expr, {1}], Level[Expand@expr, {0}]] /. Thread[Variables[Expand@expr] -> ConstantArray[1, Length[Variables[Expand@expr]]]])]]
]

SimplifyCommonFactor[expr_] := Module[{factorLCM, factorGCD},
  If[ListQ[expr], SimplifyCommonFactor[#] & /@ expr,
     {factorLCM, factorGCD} = NumberCommonDivisor[GetConstantFactor[Expand[expr]]];
     Return[Expand[factorLCM expr/factorGCD]]
   ]
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

GetIsoBasis[grp0_, pos_, ftype_, kpoint_: {{0, 0, 0}, "\[CapitalGamma]"}, ct0_: {}] := Module[{w, p, imode, latt, Wyckoff0, WyckoffSites, ZeroModeMatrix, SymmetryAdaptedBasis, Sites, IsoDispModes, IsoDispModeMatrix, OpDispMatrix, ct, ProjMat, g, basis, grpk, m, n, lp},
  grpk = GetGrpK[grp0, kpoint[[1]]];
  ct = If[Length[ct0] == 0, GetCharacters[grpk, "kpoint" -> kpoint[[2]]], ct0];
  g = Length[grpk];
  {latt, Wyckoff0} = pos;
  WyckoffSites = GetWyckoffImages[grp0, Wyckoff0];
  SymmetryAdaptedBasis = Table[
    Sites = Map[{#[[1]], #[[2]]} &, WyckoffSites[[w]]];
    IsoDispModes = MapIndexed[{First[#2], #1, 1, 0} &, Flatten[Table[Subscript[#[[2]], xyz], {xyz, 3}] & /@ Sites]];
    IsoDispModeMatrix = IdentityMatrix[3 Length[Sites]];
    OpDispMatrix = GetOpMatrix[grpk, {latt, Sites}, IsoDispModeMatrix, IsoDispModes, ftype];
    ProjMat = Table[lp = ct[[2]][[p]][[1]]; Sum[{m, n} = First@Position[ct[[1]], Keys[grpk][[i]]]; lp/g Conjugate[ct[[2]][[p, m]]] Rationalize@Normal[OpDispMatrix[[i]]], {i, g}], {p, Length[ct[[2]]]}];
    basis = Table[Complement[Orthogonalize[ProjMat[[p]].# & /@ (IsoDispModeMatrix\[Transpose])], {Table[0, {i, 3 Length[Sites]}]}], {p, Length[ct[[2]]]}];
    If[#2 == {}, Unevaluated[Sequence[]], {Table[Sites[[1]][[2]] <> " " <> #1 <> "-" <> ToString[Length@#2] <> "(" <> ToString[imode] <> ")", {imode, Length@#2}], #2}] & @@@ Thread[ct[[3]] -> basis], {w, Length@WyckoffSites}];
  IsoDispModes = MapIndexed[{First[#2], #1, 1, 0.} &, Flatten[Flatten[SymmetryAdaptedBasis, 1]\[Transpose][[1]]]];
  IsoDispModeMatrix = # & /@ Transpose[Fold[ArrayFlatten[{{#1, 0}, {0, #2}}] &, Flatten[#, 1] & /@ (#\[Transpose][[2]] & /@ SymmetryAdaptedBasis)]];
  Return[{IsoDispModes, IsoDispModeMatrix}]
]

GetSymAdaptedBasis[grp0_, grpt_, pos_, kpoint_: {{0, 0, 0}, "\[CapitalGamma]"}, ftype_] := Module[{i, j, k, t, l, c1, c2, w, p, imode, WyckoffSites, ZeroModeMatrix, SymmetryAdaptedBasis, Latt, tran, Wyckoff0, Sites, IsoDispModes, IsoDispModeMatrix, OpDispMatrix, character, ireps0, ireps, classes, ProjMat, g0, basis, irepbasis, ireplabel, grpk, m, n, lp, BasisMatrix, BasisModes, grp, ET, cell},
  tran = {ToExpression["\!\(\*SubscriptBox[\(t\), \(1\)]\)"], ToExpression["\!\(\*SubscriptBox[\(t\),  \(2\)]\)"], ToExpression["\!\(\*SubscriptBox[\(t\), \(3\)]\)"]};
  {Latt, Wyckoff0} = pos;
  cell = GetCellFromGrp[grpt];
  grpk = GetGrpK[grp0, kpoint[[1]]];
  classes = GetClasses[grpk];
  {ireps0, character} = GetSpgIreps[grp0, kpoint, Latt, "print"->False];
  g0 = Length[grp0];
  grp = {grpt, grp0};
  ireps = Flatten[Table[ET = Simplify[#[[1]] /. Thread[tran->xyz2RotTSU2[Keys[grpt][[t]]][[2]]]];
                        ET.(#[[i]] /. Thread[tran->{0, 0, 0}]), {t, Length[grpt]}, {i, g0}], 1] &/@ ireps0;
  WyckoffSites = GetWyckoffImages[grp, Wyckoff0];
  basis = Table[Sites = WyckoffSites[[w]];
                IsoDispModes = MapIndexed[{First[#2], #1, 1, 0} &, Flatten[Table[Subscript[#[[2]], xyz], {xyz, 3}] & /@ Sites]];
                IsoDispModeMatrix = IdentityMatrix[3 Length[Sites]];
                OpDispMatrix = GetMatrixRep[grp0, grpt, {Latt, Sites}, IsoDispModeMatrix, IsoDispModes, ftype];
                irepbasis = ProjectOnBasis[ireps, OpDispMatrix, #] & /@ IsoDispModeMatrix;
                irepbasis = Flatten[({#1, Sites[[1, 2]] <> " " <> kpoint[[2]] <> "-" <> #2} & @@@ #) & /@ irepbasis, 1];
                DeleteDuplicates[DeleteCases[irepbasis, {Table[0, {i, 3 Length[Sites]}], __}], #1[[1]] == #2[[1]] || #1[[1]] == -#2[[1]]&], 
                {w, Length@WyckoffSites}];
  ireplabel = Flatten[#\[Transpose][[2]] & /@ basis];
  SymmetryAdaptedBasis = Fold[ArrayFlatten[{{#1, 0}, {0, #2}}] &, #\[Transpose][[1]] & /@ basis];
  BasisMatrix = Table[{c1, c2} = NumberCommonDivisor[SymmetryAdaptedBasis[[i]]]; 
                      c1/c2 SymmetryAdaptedBasis[[i]], {i, Length@SymmetryAdaptedBasis}]\[Transpose];
  BasisModes = Table[{i, ireplabel[[i]], 1/Norm[Flatten[Latt\[Transpose].# & /@ Partition[BasisMatrix[[;; , i]], 3]]], 0}, {i, Length[BasisMatrix\[Transpose]]}];
  Return[{BasisModes, BasisMatrix}]
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

GetSiteCluster[spg0_, PosVec_, ftype_?StringQ, OptionsPattern[{"IsoTranRules"->{}}]] := Block[{xyzStrData, xyzRotTranData, xyzTranslation, xyzRotData, field, newpos, newvec, difftable, diff, posmap, i, j, det},
  det = If[ftype=="spin",1,2];
  xyzStrData = Which[AssociationQ[spg0], Keys[spg0], ListQ[spg0], spg0];
  xyzRotTranData = Table[ToExpression["{" <> xyzStrData[[i]] <> "}"], {i, Length[xyzStrData]}];
  xyzTranslation = xyzRotTranData /. {ToExpression["x"] -> 0, ToExpression["y"] -> 0, ToExpression["z"] -> 0};
  xyzRotData = xyzRotTranData - xyzTranslation;

  If[OptionValue["IsoTranRules"]=={},
    newvec = If[ftype=="spin"||ftype=="disp",
                Table[Det[Expr2Rot[op]]^det op /. Thread[ToExpression[{"x", "y", "z"}] -> #] & /@ (PosVec\[Transpose][[2]]), {op, xyzRotData}],
                Print["Only l=1 type vector field can be transformed without given IsoTranRules"];
                Abort[];
      ];
    newpos = Table[op /. Thread[ToExpression[{"x", "y", "z"}] -> #] & /@ (PosVec\[Transpose][[1]]), {op, xyzRotData}];,
    newvec = Table[PosVec\[Transpose][[2]] /. op, {op, OptionValue["IsoTranRules"]}];
    newpos = Table[op /. Thread[ToExpression[{"x", "y", "z"}] -> #] & /@ (PosVec\[Transpose][[1]]), {op, xyzRotData}];
  ];
 
  Return[Rationalize[{#[[1]], #[[2]]}]\[Transpose] &/@ ({newpos, newvec}\[Transpose])]
]

GetSiteInt[spg0_, field_, OptionsPattern[{"IsoTranRules"->{}}]] := Module[{tij, factor, OpMatrix},
  tij = Total@Fold[Times, Table[FieldCode2var[#, posvec[[1]]] & /@ Flatten[GetSiteCluster[spg0, {posvec[[2]]}, posvec[[3]], "IsoTranRules"->OptionValue["IsoTranRules"]], 1], {posvec, field}]];
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
