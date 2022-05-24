BeginPackage["LINVARIANT`INVARIANT`",{"LINVARIANT`Structure`","LINVARIANT`GroupTheory`", "LINVARIANT`MathematicaPlus`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
ShowIsoModes                 ::usage = "ShowIsoModes[PosVec]"
GetIsoBasis                  ::usage = "GetIsoBasis[grp0, Wyckoff0]"
GetSymAdaptedBasis           ::usage = "GetSymAdaptedBasis[grp0, pos, kpoint, ftype, ct0]"
ISODISTORT                   ::usage = "ISODISTORT[R0, pos0, pos, IsoMatrix, label]"
ImposeMode                   ::usage = "ImposeMode[Wyckoff0, IsoMatrix, modeset, s]"
GetIRvector                  ::usage = "GetIRvector[id, pos]"
GetBasisField                ::usage = "GetBasisField[id, BasisMatrix, pos]"
ExportBasis                  ::usage = "ExportBasis[fname, BasisMatrix, BasisLabel, pos, ftype]"
GetOpMatrix                  ::usage = "GetOpMatrix[SymOpFile, pos, IsoMatrix, Modes]"
GetMatrixRep                 ::usage = "GetMatrixRep[SymOpFile, pos, IsoMatrix, Modes]"
ExportMatrixRep              ::usage = "ExportMatrixRep[SymOpFile, pos, IsoMatrix, Modes]"
SymmetryOpBasisField         ::usage = "SymmetryOpBasisField[grp, BasisField]"
MonomialTransform            ::usage = "MonomialTransform[monomials, rules]"
GetIsoVars                   ::usage = "GetIsoVars[IsoDispModes]"
GetIsoTransformRules         ::usage = "GetIsoDispTransformRules[OpMatrix, IsoDispModes, TranType]"
GetIsoStrainTransformRules   ::usage = "GetIsoStrainTransformRules[GridSymFile]"
Epsilon2Field                ::usage = "Epsilon2Field[strain]"
Field2Epsilon                ::usage = "Field2Epsilon[field]"
GetEpsilonijRule             ::usage = "GetEpsilonijRule[symfile]"
GetTransformationRules       ::usage = "GetTransformationRules[spg0, OpMatrix]"
StrainInvQ                   ::usage = "StrainInvQ[x, e]"
NumPolynomialVar             ::usage = "NumPolynomialVar[x]"
PolynomialOrder              ::usage = " PolynomialOrder[x, vars]"
SortInvByNumVar              ::usage = "SortInvByNumVar[list]"
SortInvariants               ::usage = "SortInvariants[list, vars]"
InvariantCharacter           ::usage = "InvariantCharacter[inv, vars]"
GetInvariants                ::usage = "GetInvariants[seeds, order, AllModes, OpMatrix, GridSymFile]"
GenMonomialList              ::usage = "GenMonomialList[seeds, n]"
NumberCommonDivisor          ::usage = "NumberCommonDivisor[NumList]"
GetConstantFactor            ::usage = "GetConstantFactor[expr]"
SimplifyCommonFactor         ::usage = "SimplifyCommonFactor[expr]"
ImposeDW                     ::usage = "ImposeDW[Wyckoff0, IsoMatrix, modeset, {Nx, Ny, Nz}]"
ImposeIsoStrainVariedDspMode ::usage = "ImposeIsoStrainVariedDspMode[Wyckoff0, IsoMatrix, modeset, LV]"
ShowInvariantTable           ::usage = "ShowInvariantTable[TableData]"
Jij                          ::usage = "Jij[r0, MeshDim]"
FieldCode2var                ::usage = "FieldCode2var[code, varstr]"
Var2field                    ::usage = "Var2field[expr, site]"
GetSiteInt                   ::usage = "GetSiteInt[spg0, field]"
GetSiteCluster               ::usage = "GetSiteCluster[spg0, PosVec]"
GetSiteEpsilon               ::usage = "GetSiteEpsilon[spg0, PosVec]"
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
ShowIsoModes[PosVec_, OptionsPattern[{"imagesize"->300}]] := Module[{StructData, ArrowData},
  StructData = Table[{ElementData[i[[3]], "IconColor"], 
                      Sphere[i[[1]], QuantityMagnitude[ElementData[i[[3]], "AtomicRadius"], "Angstroms"]/2], 
                      Black, 
                      Text[i[[3]]<>ToString@Position[PosVec, i][[1, 1]], 
                      i[[1]]]}, {i, PosVec}];

  ArrowData = Table[{Green, 
                     Arrowheads[0.03], 
                     Arrow@Tube[{i[[1]], 
                     i[[1]] + i[[2]]}]}, {i, PosVec}];

  Graphics3D[{StructData, ArrowData},
              ImagePadding -> {{75, 75}, {75, 75}},
              ImageSize -> OptionValue["imagesize"], 
              Axes -> True, 
              AxesLabel -> (Style[#, Bold, 32]&/@{"a", "b", "c"}), 
              ViewPoint -> {0, 0, \[Infinity]}]
]

ISODISTORT[Lattice_, pos0_, pos_, IsoMatrix_, label_, OptionsPattern[{"round"->10.0^-6, "match"->True}]] := Module[{imode, Amp, NN, posmatched, phonon, basis, origin, posshifted},
  posmatched = If[OptionValue["match"], 
                  Transpose[{PosMatchTo[Lattice, pos0\[Transpose][[1]], pos\[Transpose][[1]]][[2]], Transpose[pos0][[2]]}],
                  pos];
  origin = Total[PbcDiff[#] & /@ (posmatched\[Transpose][[1]] - pos0\[Transpose][[1]])]/Length[pos0];
  posshifted = {#1-origin, #2} &@@@ posmatched;

  Amp = Table[phonon=Flatten[Lattice\[Transpose].PbcDiff[#] & /@ (Transpose[posshifted][[1]] - Transpose[pos0][[1]])];
              basis=Normalize[Flatten[Lattice\[Transpose].# &/@ Partition[Normal[IsoMatrix[[;;, imode]]],3]]];
              phonon.basis, {imode, Length@label}];
  NN = Table[1/Norm@Flatten[Lattice\[Transpose].# & /@ Partition[Normal[IsoMatrix][[;; , imode]], 3]], {imode, Length@label}];
  Return[{Range[Length@label], label, NN, N@Round[Amp, OptionValue["round"]]}\[Transpose]]
]

ImposeMode[latt_, Wyckoff0_, IsoMatrix_, modeset_, s_] := Module[{mode, id, NN, Amp, pos},
  pos = Wyckoff0;
  Do[id = mode[[1]];
     NN = 1/Norm@Flatten[latt\[Transpose].# & /@ Partition[Normal[IsoMatrix][[;; , id]], 3]];
     Amp = s NN mode[[4]];
     pos = Transpose[{#[[1]] & /@ pos + Amp Partition[IsoMatrix[[;; , id]] // Normal, 3], First@StringCases[#,RegularExpression["[[:upper:]][[:lower:]]*"]] & /@ Transpose[pos][[2]]}], {mode, modeset}];
  Return[pos]
]

(*Jij[spg0_, OptionsPattern[{"Cluster" -> False}]] := Module[{tij, i, j}, 
  tij = DeleteDuplicates[DeleteCases[Flatten[Table[Total[Rationalize[#1[[2]]].{Subscript[Iso, 1, #1[[1]][[1]], #1[[1]][[2]], #1[[1]][[3]]], Subscript[Iso, 2, #1[[1]][[1]], #1[[1]][[2]], #1[[1]][[3]]], Subscript[Iso, 3, #1[[1]][[1]], #1[[1]][[2]], #1[[1]][[3]]]} Rationalize[#2[[2]]].{Subscript[Iso, 1, #2[[1]][[1]], #2[[1]][[2]], #2[[1]][[3]]], Subscript[Iso, 2, #2[[1]][[1]], #2[[1]][[2]], #2[[1]][[3]]], Subscript[Iso, 3, #2[[1]][[1]], #2[[1]][[2]], #2[[1]][[3]]]} & @@@
  ({Flatten[GetSiteCluster[spg0, {{{0, 0, 0}, v1}}, "spin" -> False], 1], Flatten[GetSiteCluster[spg0, {{neighbor, v2}}, "spin" -> False], 1]}\[Transpose])], {v1, IdentityMatrix[3]}, {v2, IdentityMatrix[3]}, {neighbor, Flatten[GetUpTo3rdNeighbors[], 1]}]], 0], #1 === -#2 || #1 === #2 &];
  Return[Table[Expand[tij[[j]]/(GCD @@ Table[tij[[j]][[i]] /. {Subscript[Iso, i_, x_, y_, z_] -> 1}, {i, Length@tij[[j]]}])], {j, Length@tij}]]
]*)

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
  Thread[strains -> #] & /@ (Table[Field2Epsilon[#] & /@ GetSiteEpsilon[spg0, Epsilon2Field[ep]], {ep, strains}]\[Transpose])
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

ExportBasis[fname_, BasisMatrix_, BasisLabel_, pos_, ftype_] := Module[{data, basis, dim, i},
  dim = Length[BasisLabel];
  data = Flatten[Table[basis = GetBasisField[i, BasisMatrix, BasisLabel, pos, ftype];
                       Join[{BasisLabel[[i]][[1 ;; 2]]}, {StringRiffle[ToString[NumberForm[N[#], {10, 8}]] & /@ #1], 
                            StringRiffle[ToString[NumberForm[N[#], {6, 8}]] & /@ #2]} & @@@ (basis\[Transpose][[2]])], 
                       {i, dim}], 1];
  Export[fname, data]
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

GetMatrixRep[grp0_, grpt_, pos_, BasisMatrix_, BasisLabels_, ftype_] := Module[{i, j, grp, op, op2, ir1, ir2, mat, OpMatrix, AllIRvectors, BasisField, TransformedBasisField, BasisDim, NumBasis, g, ig, field, NormFactor, latt, slatt, sites, cell, basis},
  grp = {grpt, grp0};
  {latt, sites} = pos;
  (*slatt = Lattice2Symbolic[latt] /. {ToExpression["a"]->1, ToExpression["b"]->2, ToExpression["c"]->4};*)
  slatt = Lattice2Symbolic[latt][[1]];
  {BasisDim, NumBasis} = Dimensions[BasisMatrix];
  cell = GetCellFromGrp[grpt];
  OpMatrix =
    Table[
      If[BasisDim != 0,
         {NormFactor, BasisField} = Table[
         basis=GetBasisField[i, BasisMatrix, BasisLabels, pos, ftype];
         {If[ftype!="orbital",1/ComplexExpand@Norm[Flatten[slatt\[Transpose].# &/@ (basis\[Transpose][[2]]\[Transpose][[2]])]],
                              1/Norm[Flatten[basis\[Transpose][[2]]\[Transpose][[2]]]]], 
          basis}, {i,NumBasis}]\[Transpose];
         TransformedBasisField = ParallelTable[SymmetryOpBasisField[g, pos, cell, #, ftype] &/@ BasisField, {g, Keys[grp[[ig]]]}, DistributedContexts -> {"LINVARIANT`INVARIANT`Private`"}];
         ParallelTable[mat = 
                       If[ftype!="orbital",
                          Table[Simplify@Expand[NormFactor[[i]] NormFactor[[j]] Flatten[slatt\[Transpose].# &/@ (TransformedBasisField[[g]][[i]]\[Transpose][[2]]\[Transpose][[2]])].Flatten[slatt\[Transpose].# &/@ (BasisField[[j]]\[Transpose][[2]]\[Transpose][[2]])]], {i, Range@NumBasis}, {j, Range@NumBasis}],
                          Table[Simplify@Expand[NormFactor[[i]] NormFactor[[j]] Flatten[TransformedBasisField[[g]][[i]]\[Transpose][[2]]\[Transpose][[2]]].Flatten[BasisField[[j]]\[Transpose][[2]]\[Transpose][[2]]]], {i, Range@NumBasis}, {j, Range@NumBasis}]];
         SparseArray[mat], {g, Length[grp[[ig]]]}, DistributedContexts -> {"LINVARIANT`INVARIANT`Private`"}],
         Table[{}, {Length[grp[[ig]]]}]], {ig, 2}];
  If[BasisDim != 0,
     Flatten[Table[Chop[OpMatrix[[1,i]].OpMatrix[[2,j]]], {i, Length[grpt]}, {j, Length[grp0]}], 1],
     (*Fold[Dot, #] &/@ Tuples[OpMatrix],*)
     Table[{{}}, Fold[Times, Length[#] &/@ OpMatrix]]
    ]
]

ExportMatrixRep[grp0_, grpt_, ig_, g_, pos_, BasisMatrix_, BasisLabels_, ftype_, dir_:"~/"] := Module[{i, j, grp, xyz, op, op2, ir1, ir2, mat, OpMatrix, AllIRvectors, BasisField, TransformedBasisField, BasisDim, field, NormFactor, latt, slatt, sites, cell, basis, vecleft, vecright, val},
  grp = {grpt, grp0};
  {latt, sites} = pos;
  (*slatt = Lattice2Symbolic[latt] /. {ToExpression["a"]->1, ToExpression["b"]->2, ToExpression["c"]->4};*)
  slatt = Lattice2Symbolic[latt][[1]]\[Transpose];
  BasisDim = Length@BasisMatrix;
  cell = GetCellFromGrp[grpt];
  xyz = Keys[grp[[ig]]];
  {NormFactor, BasisField} = 
     Table[basis=GetBasisField[i, BasisMatrix, BasisLabels, pos, ftype];
           vecleft=basis\[Transpose][[2]]\[Transpose][[2]];
           {If[ftype!="orbital",1/ComplexExpand@Norm[Flatten[slatt.# &/@ vecleft]],1/Norm[Flatten[vecleft]]], basis}, 
           {i,BasisDim}]\[Transpose];
  TransformedBasisField = Table[SymmetryOpBasisField[xyz[[g]], pos, cell, BasisField[[i]], ftype], {i, BasisDim}];
  mat = Table[If[ftype!="orbital", 
                 vecleft=Flatten[slatt.# &/@ (TransformedBasisField[[i]]\[Transpose][[2]]\[Transpose][[2]])];
                 vecright=Flatten[slatt.# &/@ (BasisField[[j]]\[Transpose][[2]]\[Transpose][[2]])];,
                 vecleft=Flatten[TransformedBasisField[[i]]\[Transpose][[2]]\[Transpose][[2]]];
                 vecright=Flatten[BasisField[[j]]\[Transpose][[2]]\[Transpose][[2]]];];
              Chop@Simplify@Expand[NormFactor[[i]] NormFactor[[j]] vecleft.vecright], {i, Range@BasisDim}, {j, Range@BasisDim}];
  Export[dir<>"/OpDispMat-"<>ToString[ig]<>"."<>ToString[g]<>".dat", SparseArray[mat,{BasisDim,BasisDim}], "ExpressionML"];
  Print["Saving: "<>"OpDispMat-"<>ToString[ig]<>"."<>ToString[g]<>".dat"];
]

GenMonomialList[seeds_, n_] := Module[{seed1, seed2, isolated, inter, intra, nirrep, out, i, j}, 
  If[Flatten[seeds] === seeds, 
     isolated = #^n & /@ Flatten[seeds];
     seed1 = #/First[GetConstantFactor[#]] & /@ Flatten[MonomialList[Total[Flatten[seeds]]^n]];
     intra = Complement[seed1, isolated];
     out = {isolated, intra}, 
     nirrep = Length[seeds];
     isolated = #^n & /@ Flatten[seeds];
     seed1 = #/First[GetConstantFactor[#]] & /@ Flatten[MonomialList[Total[#]^n] & /@ seeds];
     intra = Complement[seed1, isolated];
     out = {isolated, intra};
     Do[seed2 = MonomialList[Total[Flatten[#]]^n] & /@ Subsets[seeds, {i}];
        Do[AppendTo[out, Complement[#/First[GetConstantFactor[#]] & /@ (seed2[[j]]), Flatten[out]]], {j, Length[seed2]}], {i, 2, nirrep}]];
  Return[out]
]

MonomialTransform[monomials_, rules_] := Module[{P},
  P = (monomials /. #) & /@ rules;
  Return[P]
]

SymmetryOpBasisField[grp_, pos_, cell_, BasisField_, ftype_] := Module[{latt, sites, site, i, trans, rot, rotL, su2, NewField, NewFieldup, NewFielddn, posvec, difftable, posmap, updn, basis, upbasis, dnbasis, tesseral, dim, dgQ, originshift},
  originshift = 0.0*{0.25, 0.25, 0.25};
  Which[
   AssociationQ[grp],
   SymmetryOpBasisField[#, pos, cell, BasisField, ftype] & /@ Keys[grp],
   ListQ[grp],
   SymmetryOpBasisField[#, pos, cell, BasisField, ftype] & /@ grp,
   StringQ[grp] && Length[Dimensions[BasisField]] == 3,
   SymmetryOpBasisField[grp, pos, cell, #, ftype] & /@ BasisField,
   StringQ[grp] && Length[Dimensions[BasisField]] == 2,
   {latt, sites} = pos;
   {rot, trans, su2} = xyz2RotTsu2[latt, grp];
   If[ftype=="disp"||ftype=="spin",
     Which[
       ftype=="disp",
       posvec={ModCell[N[rot.(#2[[1]]+originshift)+trans], cell], rot.(#2[[2]])} & @@@ BasisField,
       ftype=="spin",
       posvec={ModCell[N[rot.(#2[[1]]+originshift)+trans], cell], Det[rot] rot.(#2[[2]])} & @@@ BasisField
       ];
     difftable = DistMatrix[pos[[1]], #+originshift&/@(BasisField\[Transpose][[2]]\[Transpose][[1]]), posvec\[Transpose][[1]], cell];
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

     difftable = DistMatrix[pos[[1]], sites\[Transpose][[1]], upbasis\[Transpose][[1]], cell];
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

GetIsoStrainTransformRulesOld[latt_, spg0_] := Module[{StrainRules},
  Which[AssociationQ[spg0],
        GetEpsilonijRule[spg0],
        ListQ[spg0],
        GetEpsilonijRule[xyz2Grp[latt, Flatten[GTimes[latt, spg0]], "fast"->True]]
  ]
]

GetIsoStrainTransformRules[latt_, spg0_] := Module[{slatt, i, j, ig, g, strain0, strain, q, rot, tran, v1, v2},
  g= Length[spg0];
  slatt = Lattice2Symbolic[latt][[1]];
  strain0 = Normal@SparseArray[{{i_, j_} /; i == j -> Subscript[Epsilon, i, j], {i_, j_} /; i < j -> Subscript[Epsilon, i, j], {i_, j_} /; i > j -> Subscript[Epsilon, j, i]}, {3, 3}];
  Table[{rot, tran}=xyz2RotT[Keys[spg0][[ig]]];
        q = Table[v1 = (rot.slatt)[[i]]; 
                  v2 = slatt[[j]];
                  ComplexExpand[v1.v2/(Norm[v1] Norm[v2])], {i, 3}, {j, 3}];
        strain = q\[Transpose].strain0.q;
        Thread[{strain0[[1,1]], strain0[[2,2]], strain0[[3,3]], strain0[[2,3]], strain0[[1,3]], strain0[[1,2]]}
             ->{strain[[1,1]], strain[[2,2]], strain[[3,3]], strain[[2,3]], strain[[1,3]], strain[[1,2]]}], {ig, g}]
]

GetTransformationRules[latt_, spg0_, OpMatrix_] := Module[{OpDispMatrix, OpSpinMatrix, OpOrbitalMatrix},
  {OpDispMatrix, OpSpinMatrix, OpOrbitalMatrix} = OpMatrix;
  Join[#1, #2, #3, #4] & @@@ ({GetIsoTransformRules[OpDispMatrix, "disp"], 
    GetIsoTransformRules[OpSpinMatrix, "spin"], GetIsoTransformRules[OpOrbitalMatrix, "spinor"],
    GetIsoStrainTransformRules[latt, spg0]}\[Transpose])
]

NumPolynomialVar[x_] := Module[{},
  Which[ListQ[x], NumPolynomialVar[#] & /@ x,
        True, If[Head[x] === Plus, If[Head[PowerExpand[Log[x[[1]]]]] === Plus, Length[PowerExpand[Log[x[[1]]]]], 1], 
                                   If[Head[PowerExpand[Log[x]]] === Plus, Length[PowerExpand[Log[x]]], 1]]]
]

PolynomialOrder[x_, vars_, OptionsPattern[{"tot" -> True}]] := Module[{},
  Which[ListQ[x], PolynomialOrder[#, vars] & /@ x, 
        Head[x] === Plus, If[OptionValue["tot"], Total@Exponent[First[Level[x, 1]], vars], Exponent[First[Level[x, 1]], vars]], 
        True, If[OptionValue["tot"], Total@Exponent[x, vars], Exponent[x, vars]]]
]

StrainInvQ[x_, e_] := Module[{}, If[ListQ[x], StrainInvQ[#, e] & /@ x, MemberQ[Level[x, All], e]]]

SortInvByNumVar[list_] := Module[{},
  SortBy[list, NumPolynomialVar[#] &]
]

SortInvariants[list_, vars_] := Module[{},
  SortBy[#, PolynomialOrder[#, vars] &] & /@ Flatten[Values[KeySortBy[GroupBy[#, Sort@InvariantCharacter[#, vars] &], Minus] & /@ Values[KeySortBy[GroupBy[list, NumPolynomialVar[#] &], Plus]]], 1]
]

InvariantCharacter[inv_, vars_] := Module[{v0, v1, sub, s, seeds},
  Which[ListQ[inv], InvariantCharacter[#, vars] & /@ inv,
        True, seeds = If[Head[inv] === Plus, Level[inv, 1], {inv}];
              DeleteDuplicates@Table[v1 = Variables[s];
                                     v0 = Complement[vars, v1];
                                     sub = Join[Thread[v0 -> ConstantArray[0, Length@v0]], Thread[v1 -> ConstantArray[1, Length@v1]]];
                                     vars /. sub, {s, seeds}]]
]

GetInvariants[TRules_, seeds_, order_, OptionsPattern[{"MustInclude"->{{{},{}}}, "Simplify"->True, "round"->10^-6}]] := Module[{fixlist, monomials, invariant, m, n, i, j, ss, factor, out, factorLCM, factorGCD, OpDispMatrix, OpSpinMatrix, OpOrbitalMatrix, ReducedOut, SortedOut, mm, fixseeds},
  out = Table[
  monomials = If[
    OptionValue["MustInclude"]=={{{},{}}}, 
    GenMonomialList[seeds,n],
    fixlist = Flatten[Table[fixseeds=MonomialList[Total[#1]^i];
                            fixseeds=Flatten[GenMonomialList[#1,i]];
                            If[i<n,{ConstantArray[i,Length@fixseeds],fixseeds}\[Transpose],{}], {i, #2}] &@@@ OptionValue["MustInclude"],2];
    Flatten[#]&/@(Table[Expand[# ss[[2]]] & /@ GenMonomialList[seeds,n-ss[[1]]], {ss, fixlist}]\[Transpose])];
  invariant = Table[DeleteDuplicates[DeleteCases[Chop[Expand[Total[MonomialTransform[m, TRules]]], OptionValue["round"]], i_/;i==0], (Chop[#1 -#2]*Chop[#1 + #2] == 0) &], {m, monomials}];
  If[OptionValue["Simplify"], SimplifyCommonFactor[#,OptionValue["round"]] &/@ invariant, invariant], {n, order}];
  ReducedOut = Table[DeleteDuplicates[#, (#1 - #2 == 0 || #1 + #2 == 0 || Expand[#1 + I #2] == 0 || Expand[#1 - I #2] == 0) &] &/@ invariant, {invariant,out}];
  SortedOut = Table[SortInvByNumVar[#] &/@ invariant, {invariant, ReducedOut}];
  Do[mm = PolynomialReduce[#, Flatten@m, Join[Flatten@seeds,Flatten[OptionValue["MustInclude"]\[Transpose][[1]]]]] & /@ Flatten[m];
        If[! (DiagonalMatrixQ[Quiet[mm\[Transpose][[1]]]] || (Flatten@m === {})), 
           Print["Warnning! ", "polynomial may not be independent: ", MatrixForm[mm]]], {m, SortedOut}];
  Return[SortedOut]
]



NumberCommonDivisor[NumList_,prec_:10^-12] := Module[{TempList, DenominatorLCM},
 TempList = Which[Head[#] === Real, Round[#,prec], Head[#] === Integer, #, Head[#] === Times, First@Level[#, {1}], Head[#] === Power, 1, Head[#] === Rational, #, (Head[#] === Complex)&&(Re[#]!=0), Abs[#], (Head[#] === Complex)&&(Re[#]==0), #] &/@ NumList;
 DenominatorLCM = If[MemberQ[TempList, _Rational], LCM @@ (Denominator[#] & /@ Extract[TempList, Position[TempList, _Rational]]), 1];
 DenominatorLCM = If[AllTrue[TempList, Negative], -DenominatorLCM, DenominatorLCM];
 Return[{DenominatorLCM, GCD @@ (TempList DenominatorLCM)}]
]

GetConstantFactor[expr_] := Module[{},
  If[ListQ[expr],  GetConstantFactor[#] & /@ expr, Return[(If[MatchQ[Expand@expr, Plus[_, __]], Level[Expand@expr, {1}], Level[Expand@expr, {0}]] /. Thread[Variables[Expand@expr] -> ConstantArray[1, Length[Variables[Expand@expr]]]])]]
]

SimplifyCommonFactor[expr_,prec_:10^-6] := Module[{factorLCM, factorGCD},
  Which[ListQ[expr], 
        SimplifyCommonFactor[#, prec] & /@ expr,
        expr === 0,
        Return[0],
        True,
        {factorLCM, factorGCD} = NumberCommonDivisor[GetConstantFactor[Expand[expr]],prec prec];
        Return[Expand[expr factorLCM/factorGCD]/. x_Real :> Round[x,prec]]
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

ShowInvariantTable[InvData_, param_, OptionsPattern[{"FontSize" -> 12, "color"->Pink}]] := Module[{colors, bg, i, j, len, newparam, NumOrder, TableData, vars, order},
  NumOrder= Length[InvData];
  bg = Flatten[If[Flatten[#] === {}, ## &[], {Gray,Table[ConstantArray[Hue[OptionValue["color"], 0.5 + i 0.5/Length[#]],Length[#[[i]]]], {i, Length[#]}]}] & /@ InvData];
  newparam = If[param===Null, Table[ConstantArray[0, Length[Flatten[InvData[[i]]]]], {i, NumOrder}], DeleteCases[param,{}]];
  Grid[Flatten[Table[If[Flatten[InvData[[i]]] != {}, vars = Variables[Flatten[InvData[[i]]][[1, 1]]];
                     order=Total@Exponent[Flatten[InvData[[i]]][[1, 1]], vars];
                     Prepend[{newparam, Flatten[#]&/@InvData}\[Transpose][[i]]\[Transpose], 
                             {"#", ToString[ToString["U"]^"("<>ToString[order]<>")",StandardForm]<>" contains                          "<>ToString[Length[Flatten[InvData[[i]]]]]<>" invariants"}], ##&[]], {i, NumOrder}], 1], 
       Background -> {{Yellow, White}, bg}, 
       Alignment -> Left, 
       Spacings -> {1, 0.5}, 
       ItemSize -> Full, 
       ItemStyle -> Directive[FontSize -> OptionValue["FontSize"], Black], 
       Frame -> All]
]

GetIsoBasis[latt_, grp0_, pos_, ftype_, kpoint_: {{0, 0, 0}, "\[CapitalGamma]"}, ct0_: {}] := Module[{w, p, imode, latt, Wyckoff0, WyckoffSites, ZeroModeMatrix, SymmetryAdaptedBasis, Sites, IsoDispModes, IsoDispModeMatrix, OpDispMatrix, ct, ProjMat, g, basis, grpk, m, n, lp},
  grpk = GetGrpK[latt, grp0, kpoint[[1]]];
  ct = If[Length[ct0] == 0, GetCharacters[latt, grpk, "kpoint" -> kpoint[[2]]], ct0];
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

GetSymAdaptedBasis[grp0_, grpt_, pos_, kpoint_: {{0, 0, 0}, "\[CapitalGamma]"}, ftype_] := Module[{i, j, k, t, l, c1, c2, w, p, imode, WyckoffSites, ZeroModeMatrix, SymmetryAdaptedBasis, latt, tran, Wyckoff0, Sites, IsoDispModes, IsoDispModeMatrix, OpDispMatrix, character, ireps0, ireps, classes, ProjMat, g0, basis, irepbasis, ireplabel, grpk, m, n, lp, BasisMatrix, BasisModes, grp, ET, cell},
  tran = {ToExpression["\!\(\*SubscriptBox[\(t\), \(1\)]\)"], ToExpression["\!\(\*SubscriptBox[\(t\),  \(2\)]\)"], ToExpression["\!\(\*SubscriptBox[\(t\), \(3\)]\)"]};
  {latt, Wyckoff0} = pos;
  cell = GetCellFromGrp[grpt];
  grpk = GetGrpK[latt, grp0, kpoint[[1]]];
  classes = GetClasses[latt, grpk];
  {ireps0, character} = GetSpgIreps[latt, grp0, kpoint, "print"->False];
  g0 = Length[grp0];
  grp = {grpt, grp0};
  ireps = Flatten[Table[ET = Simplify[#[[1]] /. Thread[tran->xyz2RotT[Keys[grpt][[t]]][[2]]]];
                        ET.(#[[i]] /. Thread[tran->{0, 0, 0}]), {t, Length[grpt]}, {i, g0}], 1] &/@ ireps0;
  WyckoffSites = GetWyckoffImages[latt, grp, Wyckoff0];
  basis = Table[Sites = WyckoffSites[[w]];
                IsoDispModes = MapIndexed[{First[#2], #1, 1, 0} &, Flatten[Table[Subscript[#[[2]], xyz], {xyz, 3}] & /@ Sites]];
                IsoDispModeMatrix = IdentityMatrix[3 Length[Sites]];
                OpDispMatrix = GetMatrixRep[grp0, grpt, {latt, Sites}, IsoDispModeMatrix, IsoDispModes, ftype];
                irepbasis = ProjectOnBasis[ireps, OpDispMatrix, #] & /@ IsoDispModeMatrix;
                irepbasis = Flatten[({#1, Sites[[1, 2]] <> " " <> kpoint[[2]] <> "-" <> #2} & @@@ #) & /@ irepbasis, 1];
                DeleteDuplicates[DeleteCases[irepbasis, {Table[0, {i, 3 Length[Sites]}], __}], #1[[1]] == #2[[1]] || #1[[1]] == -#2[[1]]&], 
                {w, Length@WyckoffSites}];
  ireplabel = Flatten[#\[Transpose][[2]] & /@ basis];
  SymmetryAdaptedBasis = Fold[ArrayFlatten[{{#1, 0}, {0, #2}}] &, #\[Transpose][[1]] & /@ basis];
  BasisMatrix = Table[{c1, c2} = NumberCommonDivisor[SymmetryAdaptedBasis[[i]]]; 
                      c1/c2 SymmetryAdaptedBasis[[i]], {i, Length@SymmetryAdaptedBasis}]\[Transpose];
  BasisModes = Table[{i, ireplabel[[i]], 1/Norm[Flatten[latt\[Transpose].# & /@ Partition[BasisMatrix[[;; , i]], 3]]], 0}, {i, Length[BasisMatrix\[Transpose]]}];
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

Var2field[expr_, site_] := Module[{},
  Which[MatchQ[Head[expr], Plus], 
  Plus @@ (Var2field[#, site] & /@ Level[expr, 1]),
  MatchQ[Head[expr], Times], 
  Times @@ (Var2field[#, site] & /@ Level[expr, 1]),
  MatchQ[Head[expr], Power],
  Var2field[Level[expr, 1][[1]], site]^(Level[expr,1][[2]]),
  MatchQ[Head[expr], Subscript], 
  Subscript @@ Join[GetSubscriptInfo[expr], site],
  NumberQ[expr], expr]
]

GetSiteEpsilon[spg0_, PosVec_] := Block[{xyzStrData, rot, tran, op, newpos, newvec, cell},
  cell = GetCellFromGrp[spg0];
  xyzStrData = Which[AssociationQ[spg0], Keys[spg0], ListQ[spg0], spg0];
  newvec = Table[rot=xyz2RotT[op][[1]];
                 rot.# & /@ (PosVec\[Transpose][[2]]), {op, xyzStrData}];
  newpos = Table[{rot,tran}=xyz2RotT[op];
                 rot.# & /@ (PosVec\[Transpose][[1]]), {op, xyzStrData}];

  Return[Rationalize[{#[[1]], #[[2]]}]\[Transpose] &/@ ({newpos, newvec}\[Transpose])]
]

GetSiteCluster[spg0_, TRules_, field_] := Block[{xyzStrData, rot, tran, sites, clusters, cell},
  cell = GetCellFromGrp[spg0];
  xyzStrData = Which[AssociationQ[spg0], Keys[spg0], ListQ[spg0], spg0];

  sites = Table[{rot,tran}=xyz2RotT[op]; rot.# & /@ (field\[Transpose][[1]]), {op, xyzStrData}];
  clusters = field\[Transpose][[2]] /.# &/@ TRules;

  Return[Rationalize[{#[[1]], #[[2]]}]\[Transpose] &/@ ({sites, clusters}\[Transpose])]
]

GetCellCellInt[spg0_, TRules_, field_] := Module[{g, FieldTransformed, tij, factor, OpMatrix, f},
  g = Length[spg0];
  FieldTransformed = GetSiteCluster[spg0, TRules, field];
  tij = Sum[Times@@(Var2field[#2, #1] &@@@ FieldTransformed[[i]]), {i,g}];
  Return[SimplifyCommonFactor[tij]]
] 

Jij[spg0_, TRules_, fielddef_?ListQ, OptionsPattern[{"OnSite" -> False, "IsoTranRules"->{}}]] := Module[{tij, i, j, n, field, v0, v1, neighbor, nlist},
  nlist = Flatten[GetUpTo3rdNeighbors["OnSite" -> OptionValue["OnSite"]], 1];
  tij = Expand[DeleteDuplicates[DeleteCases[Flatten[
          Table[field = {{{0, 0, 0}, v0}, {n, v1}};
                GetCellCellInt[spg0, TRules, field], {v0, fielddef[[1]]}, {v1, fielddef[[2]]}, {n, nlist}]], 0], #1 === -#2 || #1 === #2 &]];
  Return[tij]
]


InvariantEOnSite[expr_] := Module[{X, i, dx, dy, dz},
  expr /. {Subscript[X_, i_, dx_, dy_, dz_] :> Subscript[X, i, ToExpression["ix"] + dx, ToExpression["iy"] + dy, ToExpression["iz"] + dz]}
]
(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
