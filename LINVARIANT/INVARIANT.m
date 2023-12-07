BeginPackage["LINVARIANT`INVARIANT`",{"LINVARIANT`Structure`","LINVARIANT`GroupTheory`", "LINVARIANT`MathematicaPlus`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
ShowIsoModes                 ::usage = "ShowIsoModes[PosVec]"
GetIsoBasis                  ::usage = "GetIsoBasis[grp0, Wyckoff0]"
GetSymAdaptedBasis           ::usage = "GetSymAdaptedBasis[grp0, pos, kpoint, ftype, ct0]"
ISODISTORT                   ::usage = "ISODISTORT[R0, pos0, pos, IsoMatrix, label]"
ModeDecomposition            ::usage = "ModeDecomposition[pos0, pos, Basis, label]"
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
GetReciprocalTRules          ::usage = "GetReciprocalTRules[latt, spg0, k]"
JijInvQ                      ::usage = "JijInvQ[x]"
NumPolynomialVar             ::usage = "NumPolynomialVar[x]"
PolynomialOrder              ::usage = " PolynomialOrder[x, vars]"
SortInvByNumVar              ::usage = "SortInvByNumVar[list]"
SortInvariants               ::usage = "SortInvariants[list, vars]"
InvariantCharacter           ::usage = "InvariantCharacter[inv, vars]"
GetInvariants                ::usage = "GetInvariants[seeds, order, AllModes, OpMatrix, GridSymFile]"
GenMonomialList              ::usage = "GenMonomialList[seeds, n]"
ImposeDW                     ::usage = "ImposeDW[Wyckoff0, IsoMatrix, modeset, {Nx, Ny, Nz}]"
ImposeDomains                ::usage = "ImposeDomains[pos0, Basis, mode, grid]"
ImposeIsoStrainVariedDspMode ::usage = "ImposeIsoStrainVariedDspMode[Wyckoff0, IsoMatrix, modeset, LV]"
ShowInvariantTable           ::usage = "ShowInvariantTable[TableData]"
Jij                          ::usage = "Jij[r0, MeshDim]"
GetJijStar                   ::usage = "GetJijStar[expr, nn, vars]"
FieldCode2var                ::usage = "FieldCode2var[code, varstr]"
Var2field                    ::usage = "Var2field[expr, site]"
GetSiteInt                   ::usage = "GetSiteInt[spg0, field]"
GetSiteEpsilon               ::usage = "GetSiteEpsilon[spg0, PosVec]"
InvariantEOnSite             ::usage = "InvariantEOnSite[xyz, expr]"
BasisShift                   ::usage = "BasisShift[pos0, spg0, Basis]"
InvariantFolding             ::usage = "InvariantFolding[model]"
GridSymTransformation        ::usage = "GridSymTransformation[site0, xyz]"
CollectPolynomial            ::usage = "CollectPolynomial[expr]"
HModelFolding                ::usage = "HModelFolding[hmodel, grid]"
FoldHopping                  ::usage = "FoldHopping[inv, grid]"
JijTRules                    ::usage = "JijTRules[spg0, TRules, dist]"
PlotPolynomial               ::usage = "PlotPolynomial[model]"
IncludeAcousticMode          ::usage = "IncludeAcousticMode[pos0, Basis, spg0, tgrp, vars]"
GetInvariantTensor           ::usage = "GetInvariantTensor[T0, spg0, symmetric]"
BasisInSupercell             ::usage "BasisInSupercell[pos0, Basis, vars, NSC, SymmetricQ]"

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

ISODISTORT[Lattice_, pos0_, pos_, IsoMatrix_, label_, OptionsPattern[{"round"->10.0^-6, "match"->True}]] := Module[{imode, Amp, NN, posmatched, phonon, basis},
  posmatched = If[OptionValue["match"], 
                  Transpose[{PosMatchTo[Lattice, pos0\[Transpose][[1]], pos\[Transpose][[1]], "OriginShift" -> True][[2]], Transpose[pos0][[2]]}],
                  pos];

  Amp = Table[phonon=Flatten[Lattice\[Transpose].PbcDiff[#] & /@ (Transpose[posmatched][[1]] - Transpose[pos0][[1]])];
              basis=Normalize[Flatten[Lattice\[Transpose].# &/@ Partition[Normal[IsoMatrix[[;;, imode]]],3]]];
              phonon.basis, {imode, Length@label}];
  NN = Table[1/Norm@Flatten[Lattice\[Transpose].# & /@ Partition[Normal[IsoMatrix][[;; , imode]], 3]], {imode, Length@label}];
  Return[{Range[Length@label], label, NN, Chop[Amp, OptionValue["round"]]}\[Transpose]]
]

ModeDecompositionOld[pos0_, pos_, Basis_, label_, OptionsPattern[{"lattnorm" -> False, "table" -> True, "round" -> 10.0^-6, "match" -> True}]] := Module[{latt0, sites0, latt, sites, eij, ix, iy, iz, ngx, ngy, ngz, block, tab, grid, s, lattnorm},
  {latt0, sites0} = pos0;
  {latt, sites} = pos;
  lattnorm = If[OptionValue["lattnorm"], Norm[latt0], 1];
  grid = Round[(Norm[#] & /@ latt)/(Norm[#] & /@ latt0)];
  {ngx, ngy, ngz} = grid;
  
  eij = Chop@GetStrainTensor[latt0, latt/grid, "iso" -> False];
  
  data = Association@Table[block = Table[If[And@@Thread[{ix, iy, iz} - 1 <= s[[1]] grid < {ix, iy, iz}], 
                                            {Mod[s[[1]] grid, 1], s[[2]]}, 
                                            ## &[]], {s, sites}]; 
                                         {ix, iy, iz} -> ({#1, #2, #3, #4/lattnorm} & @@@ ISODISTORT[pos0[[1]], pos0[[2]], block, Basis, label, 
                                                                                                     "round" -> OptionValue["round"], "match" -> OptionValue["match"]]), {iz, ngz}, {iy, ngy}, {ix, ngx}];
  
  tab = Grid[Join[{{MatrixForm@eij}}, 
                  {Flatten[Table[Grid[{{{ix, iy, iz}}, {MatrixForm[data[{ix, iy, iz}]]}}], {iz, ngz}, {iy, ngy}, {ix, ngx}], 2]}], Dividers -> All];

  If[OptionValue["table"], Print[tab]];

  Return[{data, eij2eta[eij]}]
]

ModeDecomposition[pos0_, pos_, Basis_, vars_, svars_, uvars_, OptionsPattern[{"zyx" -> True, "lattnorm" -> False, "table" -> True, "round" -> 10.0^-6, "upto" -> 4, "ndigits"->2, "match" -> True}]] := Module[{latt0, sites0, supersites0, latt, sites, eij, ix, iy, iz, ngx, ngy, ngz, block, tab, grid, i, s, lattnorm, SuperVars, NumAtom, LattCar2Dir, AcousticBasis, FullBasis, SuperBasis, tabdata, data},
  
  {latt0, sites0} = pos0;
  NumAtom = Length[sites0];
  LattCar2Dir = Inverse[Transpose[latt0]];
  AcousticBasis = Transpose[Flatten[Chop@Table[LattCar2Dir . (Partition[Normalize[#], 3][[i]]), {i, NumAtom}]] & /@ (Flatten[Table[1.0 IdentityMatrix[3], {NumAtom}], 1]\[Transpose])];
  FullBasis = If[uvars === {}, Basis, Transpose@Join[Transpose@Basis, Transpose@AcousticBasis]];
  
  {latt, sites} = pos;
  lattnorm = If[OptionValue["lattnorm"], Norm[latt0], 1];
  grid = Round[(Norm[#] & /@ latt)/(Norm[#] & /@ latt0)];
  {ngx, ngy, ngz} = grid;
  
  supersites0 = Flatten[Flatten[Table[{(#[[1]] + {ix - 1, iy - 1, iz - 1})/grid, #[[2]]}, {iz, 1, ngz}, {iy, 1, ngy}, {ix, 1, ngx}], 2] & /@ (pos0[[2]]), 1];
  
  {SuperBasis, SuperVars} = BasisInSupercell[pos0, FullBasis, Append[vars, uvars], grid, False, "rotate" -> OptionValue["zyx"]];
  
  eij = Chop@GetStrainTensor[latt0, latt/grid, "iso" -> False];
  
  data = {#1, #2, #3, #4/lattnorm} & @@@ ISODISTORT[pos0[[1]], supersites0, sites, SuperBasis, SuperVars, "round" -> OptionValue["round"], "match" -> OptionValue["match"]];
  tabdata= {#1, #2, NumberForm[#3, {10,OptionValue["ndigits"]}], NumberForm[#4, {10,OptionValue["ndigits"]}]} &@@@ data; 
  tab = Grid[{{MatrixForm[eij]}, {Grid[Partition[Grid[#] & /@ GatherBy[tabdata, VarGrid[#[[2]]] &], UpTo[OptionValue["upto"]]], Dividers -> All]}}, 
             Frame -> True, Dividers -> All];
  
  If[OptionValue["table"], Print[tab]];
  Return[{GatherBy[#2 -> #4 & @@@ data, VarGrid[First@#]&], Thread[svars -> eij2eta[eij]]}]
  
]

ImposeMode[latt_, Wyckoff0_, IsoMatrix_, modeset_, s_] := Module[{mode, id, NN, Amp, pos, tmp},
  pos = Wyckoff0;
  Do[id = mode[[1]];
     NN = 1/Norm@Flatten[latt\[Transpose].# & /@ Partition[Normal[IsoMatrix][[;; , id]], 3]];
     Amp = s NN mode[[4]];
     tmp = PosMatchTo[latt, #1 &@@@ Wyckoff0, (#[[1]] & /@ pos) + Amp Partition[IsoMatrix[[;; , id]] // Normal, 3], 
                      "OriginShift" -> True, "SymmetricQ" -> False][[2]];
     pos = Transpose[{tmp, First@StringCases[#,RegularExpression["[[:upper:]][[:lower:]]*"]] & /@ Transpose[pos][[2]]}], {mode, modeset}];
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
  (*slatt = Lattice2Symbolic[latt][[1]];*)
  {BasisDim, NumBasis} = Dimensions[BasisMatrix];
  cell = GetCellFromGrp[grpt];
  OpMatrix =
    Table[
      If[BasisDim != 0,
         {NormFactor, BasisField} = Table[
         basis=GetBasisField[i, BasisMatrix, BasisLabels, pos, ftype];
         {If[ftype!="orbital",1/ComplexExpand@Norm[Flatten[latt\[Transpose].# &/@ (basis\[Transpose][[2]]\[Transpose][[2]])]],
                              1/Norm[Flatten[basis\[Transpose][[2]]\[Transpose][[2]]]]], 
          basis}, {i,NumBasis}]\[Transpose];
         TransformedBasisField = ParallelTable[SymmetryOpBasisField[g, pos, cell, #, ftype] &/@ BasisField, {g, Keys[grp[[ig]]]}, DistributedContexts -> {"LINVARIANT`INVARIANT`Private`"}];
         ParallelTable[mat = 
                       If[ftype!="orbital",
                          Table[Simplify@Expand[NormFactor[[i]] NormFactor[[j]] Flatten[latt\[Transpose].# &/@ (TransformedBasisField[[g]][[i]]\[Transpose][[2]]\[Transpose][[2]])].Flatten[latt\[Transpose].# &/@ (BasisField[[j]]\[Transpose][[2]]\[Transpose][[2]])]], {i, Range@NumBasis}, {j, Range@NumBasis}],
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
     isolated = DeleteDuplicates[#^n & /@ Flatten[seeds]];
     seed1 = DeleteDuplicates[#/First[GetConstantFactor[#]] & /@ Flatten[MonomialList[Total[Flatten[seeds]]^n]]];
     intra = Complement[seed1, isolated];
     out = {isolated, intra}, 
     nirrep = Length[seeds];
     isolated = DeleteDuplicates[#^n & /@ Flatten[seeds]];
     seed1 = DeleteDuplicates[#/First[GetConstantFactor[#]] & /@ Flatten[MonomialList[Total[#]^n] & /@ seeds]];
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

SymmetryOpBasisField[grp_, pos_, cell_, BasisField_, ftype_] := Module[{latt, sites, site, i, trans, rot, rotL, su2, NewField, NewFieldup, NewFielddn, posvec, difftable, posmap, updn, basis, upbasis, dnbasis, tesseral, dim, dgQ, originshift, TRQ},
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
   {rot, trans, su2, TRQ} = xyz2RotTsu2[grp];
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
       If[su2 == Null, Print["Error: For spinor, you need a double group!"];Abort[]];
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

GetReciprocalTRules[latt_, op_, k_] := Module[{out},
  Which[AssociationQ[op],
        GetReciprocalTRules[latt, #, k] & /@ Keys[op],
        ListQ[op],
        GetReciprocalTRules[latt, #, k] & /@ op,
        StringQ[op],
        out = Thread[k -> First[xyz2RotT[GInverse[latt, op]]] . k];
        Return[out]]
]

NumPolynomialVar[x_] := Module[{},
  Which[ListQ[x], NumPolynomialVar[#] & /@ x,
        True, If[Head[x] === Plus, If[Head[PowerExpand[Log[x[[1]]]]] === Plus, Length[PowerExpand[Log[x[[1]]]]], 1], 
                                   If[Head[PowerExpand[Log[x]]] === Plus, Length[PowerExpand[Log[x]]], 1]]]
]

PolynomialOrder[x_, vars_, OptionsPattern[{"tot" -> True}]] := Module[{fullvars},
  fullvars = If[vars === {}, Variables[x], vars];
  Which[ListQ[x], PolynomialOrder[#, fullvars] & /@ x, 
        Head[x] === Plus, If[OptionValue["tot"], Total@Exponent[First[Level[x, 1]], fullvars], Exponent[First[Level[x, 1]], fullvars]], 
        True, If[OptionValue["tot"], Total@Exponent[x, fullvars], Exponent[x, fullvars]]]
]

JijInvQ[x_] := Module[{}, If[ListQ[x], JijInvQ[#] & /@ x, Or@@(Length[Level[#, 1]] > 3 & /@ Variables[x])]]

SortInvByNumVar[list_] := Module[{},
  SortBy[list, NumPolynomialVar[#] &]
]

SortInvariants[list_, vars_] := Module[{},
  SortBy[#, PolynomialOrder[#, vars] &] & /@ Flatten[Values[KeySortBy[GroupBy[#, Sort@InvariantCharacter[#, vars] &], Minus] & /@ Values[KeySortBy[GroupBy[list, NumPolynomialVar[#] &], Plus]]], 1]
]

InvariantCharacter[inv_, vars_] := Module[{v0, v1, sub, s, seeds},
  Which[ListQ[inv], InvariantCharacter[#, Flatten@vars] & /@ inv,
        True, seeds = If[Head[inv] === Plus, Level[inv, 1], {inv}];
              DeleteDuplicates@Table[v1 = Variables[s];
                                     If[MemberQ[VarBare[v1], #], 1, 0] &/@ VarBare[Flatten@vars], {s, seeds}]]
]

GetInvariants[TRules_, seeds_, order_, OptionsPattern[{"check" -> True, "MustInclude"->{{{},{}}}, "Simplify"->True, "eventerms"->False, "round"->10^-6}]] := Module[{fixlist, monomials, invariant, m, n, i, j, ss, factor, out, factorLCM, factorGCD, OpDispMatrix, OpSpinMatrix, OpOrbitalMatrix, ReducedOut, SortedOut, mm, fixseeds, eventerms},
  out = Table[
  monomials = If[
    OptionValue["MustInclude"]=={{{},{}}}, 
    GenMonomialList[seeds,n],
    fixlist = Flatten[Table[fixseeds=MonomialList[Total[#1]^i];
                            fixseeds=Flatten[GenMonomialList[#1,i]];
                            If[i<n,{ConstantArray[i,Length@fixseeds],fixseeds}\[Transpose],{}], {i, #2}] &@@@ OptionValue["MustInclude"],2];
    Flatten[#]&/@(Table[Expand[# ss[[2]]] & /@ GenMonomialList[seeds,n-ss[[1]]], {ss, fixlist}]\[Transpose])];
    invariant = DeleteDuplicates[DeleteCases[Table[Chop[Expand[Total[MonomialTransform[m, TRules]]], OptionValue["round"]], {m, Flatten@monomials}], i_/;i==0], (Chop[#1 -#2]*Chop[#1 + #2] == 0) &];
  invariant = If[OptionValue["Simplify"], SimplifyCommonFactor[#,OptionValue["round"]] &/@ invariant, invariant];
    If[OptionValue["eventerms"], If[AllTrue[PolynomialOrder[#, {}, "tot"->False], EvenQ], #, ##&[]] &/@ invariant, invariant], {n, order}];

     (*ReducedOut = Table[DeleteDuplicates[#, (#1 - #2 == 0 || #1 + #2 == 0 || Expand[#1 + I #2] == 0 || Expand[#1 - I #2] == 0) &] &/@ invariant, {invariant,out}];
  SortedOut = Table[SortInvByNumVar[#] &/@ invariant, {invariant, ReducedOut}];*)

  SortedOut = SortInvByNumVar[#] & /@ out;

  If[OptionValue["check"],
     Do[mm = PolynomialReduce[#, Flatten@m, Join[Flatten@seeds,Flatten[OptionValue["MustInclude"]\[Transpose][[1]]]]] & /@ Flatten[m];
        If[! (DiagonalMatrixQ[Quiet[mm\[Transpose][[1]]]] || (Flatten@m === {})), 
        Print["Warnning! ", "polynomial may not be independent: ", MatrixForm[mm]]], {m, SortedOut}]];
  Return[SortedOut]
]

GetInvariantTensor[T0_, latt_, spg0_, symmetric_ : False] := Module[{TT, U, rank, i, s, contraction, out},
  rank = TensorRank[T0];
  contraction = Table[{2 i, 2 rank + i}, {i, rank}];
  out = Sum[U = Inverse[latt].xyz2RotT[s][[1]].latt;
            TT = T0;
            Do[TT = U\[TensorProduct]TT, {i, rank}];
            TensorContract[TT, contraction], {s, Keys@spg0}];
  If[symmetric, out = Expand[SimplifyTensorCommonFactor[out\[Transpose] + out]]];
  Return[SimplifyTensor[T0, out]]
]

ImposeDW[pos0_, IsoMatrix_, modeset_, {Nx_, Ny_, Nz_}] := Module[{mode, id, Amp, latt0, sites0, latt, sites, s, ix, iy, iz, Superpos},
  {latt0, sites0} = pos0;
  Superpos = Table[{Mod[(#1 + {ix, iy, iz})/{Nx, Ny, Nz},1], #2} & @@@ sites0, {ix, 0, Nx - 1}, {iy, 0, Ny - 1}, {iz, 0, Nz - 1}];
  Do[sites = Superpos[[ix]][[iy]][[iz]];
     s = Cos[2 Pi {1/Nx, 1/Ny, 1/Nz}.{ix, iy, iz}];
     Do[id = mode[[1]]; Amp = s mode[[3]] mode[[4]];
        sites = Transpose[{#[[1]] & /@ sites + Amp If[IntegerQ[id], # & /@ Partition[IsoMatrix[[;; , id]] // Normal, 3], Print["mode not exist!"]], First@StringCases[#, RegularExpression["[[:upper:]][[:lower:]]*"]] & /@ Transpose[sites0][[2]]}], {mode, modeset}
        ];
     Superpos[[ix]][[iy]][[iz]] = sites, {ix, Nx}, {iy, Ny}, {iz, Nz}
     ];
  latt = {Nx latt0[[1]], Ny latt0[[2]], Nz latt0[[3]]};
  Return[{latt, SortBy[Flatten[Superpos, 3], #[[2]]&]}]
]

ImposeDomains[pos0_, Basis_, eij_, mode_, grid_] := Module[{block, k, ngx, ngy, ngz, ix, iy, iz, pos, site, sites, latt},
  {ngx, ngy, ngz} = grid;
  block = Association@Flatten[Table[site = {ix, iy, iz};
                                    pos = ImposeDW[pos0, Basis, mode[site], {1, 1, 1}];
                                    site -> pos[[2]], {iz, ngz}, {iy, ngy}, {ix, ngx}]];
  
  sites = SortBy[Flatten[Table[{(#1 + k - 1)/grid, #2} & @@@ (block[k]), {k, Keys@block}], 1], #[[2]] &];
  latt = (IdentityMatrix[3] + eij).(pos0[[1]]) grid;
  
  Return[{latt, sites}]
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

ShowInvariantTable[InvData_, param_, OptionsPattern[{"FontSize" -> 12, "color"->Pink}]] := Module[{colors, bg, i, j, len, newparam, NumOrder, TableData, vars, order, ind, tmp, count},
  NumOrder= Length[InvData];
  bg = Flatten[If[Flatten[#] === {}, ## &[], {Gray,Table[ConstantArray[Hue[OptionValue["color"], 0.5 + i 0.5/Length[#]],Length[#[[i]]]], {i, Length[#]}]}] & /@ InvData];

  bg = Flatten[If[Flatten[#] === {}, ## &[], {Gray, ConstantArray[Hue[OptionValue["color"], 0.5], Length[#]]}] & /@ InvData];
  newparam = If[param===Null, Table[ConstantArray[0, Length[InvData[[i]]]], {i, NumOrder}], DeleteCases[param,{}]];

  count = 0;
  ind = Table[tmp = count + Range[Length[InvData[[i]]]];
              count = count + Length[InvData[[i]]];
              tmp, {i, NumOrder}];

  Grid[Flatten[Table[If[InvData[[i]] != {}, 
                        vars = Variables[InvData[[i]]];
                        order = i;
                        Prepend[{ind, newparam, Flatten[#]&/@InvData}\[Transpose][[i]]\[Transpose], 
                                {"#", "Parameters", ToString[ToString["U"]^"("<>ToString[order]<>")",StandardForm]<>" contains                          "<>ToString[Length[InvData[[i]]]]<>" invariants"}], ##&[]], {i, NumOrder}], 1], 
       Background -> {{Gray, Yellow, White}, bg}, 
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

GetSiteClusterOld[spg0_, TRules_, field_] := Block[{xyzStrData, rot, tran, sites, clusters, cell},
  cell = GetCellFromGrp[spg0];
  xyzStrData = Which[AssociationQ[spg0], Keys[spg0], ListQ[spg0], spg0];

  sites = Table[{rot,tran}=xyz2RotT[op]; rot.# & /@ (field\[Transpose][[1]]), {op, xyzStrData}];
  clusters = field\[Transpose][[2]] /.# &/@ TRules;

  Return[Rationalize[{#[[1]], #[[2]]}]\[Transpose] &/@ ({sites, clusters}\[Transpose])]
]

JijTRules[pos0_, spg0_, TRules_, Basis_, dist_] := Module[{xyzStrData, i, j, s, v, a, o1, o2, o3, s1, s2, s3, ng, field, out, KeyRule, ValueRule, v1, v2, tgrp, numvar, dim, rules0, rules, uvar, AcousticBasis, FullBasis, FullOpMat, spg}, 
   {dim, numvar} = Dimensions[Basis];
   tgrp = <|"x,y,z" -> {{{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}}, False}|>;

   rules0 = If[TRules ==={},  
               uvar = Subscript[ToExpression["u"], #] & /@ Range[3];
               AcousticBasis = (Normalize[#] & /@ (Flatten[Table[1.0 IdentityMatrix[3], {dim/3}], 1]\[Transpose]))\[Transpose];
               FullBasis = Transpose@Join[Basis\[Transpose], AcousticBasis\[Transpose]];
               FullOpMat = Round[GetMatrixRep[spg0, tgrp, pos0, FullBasis, Range[Length[numvar] + 3], "disp"], 0.01];
               Rationalize@Chop[Join[#1, #2] & @@@ ({GetIsoTransformRules[FullOpMat, "disp"],
                                                     GetIsoStrainTransformRules[pos0[[1]], spg0]}\[Transpose])],
               TRules];

  {spg, rules} = Transpose[If[Norm[xyz2RotT[#1][[2]]] == 0, {#1, #2}, ## &[]] & @@@ ({If[AssociationQ[spg0], Keys[spg0], spg0], rules0}\[Transpose])];

  {ng, dim} = Dimensions[rules];
  field = GetNeighborList[dist, "OnSite" -> True];
  xyzStrData = Which[AssociationQ[spg], Keys[spg], ListQ[spg], spg];

  KeyRule = {Subscript[v_, a_] :> Subscript[v, a, o1, o2, o3]};
  ValueRule = {Subscript[v_, a_] :> Subscript[v, a, s1, s2, s3]};
  
  out = Table[{o1, o2, o3} = s;
              {s1, s2, s3} = xyz2RotT[xyzStrData[[i]]][[1]].s;
              {s1, s2, s3} = GridSymTransformation[s, xyzStrData[[i]]];
              {v1, v2} = Level[rules[[i, j]], 1];
              (v1 /. KeyRule) -> (v2 /. ValueRule), {i, ng}, {j, dim}, {s, field}];
  Return[DeleteDuplicates[Flatten[#]] & /@ out]
]

Jij[pos0_, spg0_, dist_, TRules_, Basis_, fieldlist_?ListQ, OptionsPattern[{"OnSite" -> False}]] := Module[{tij, i, j, n, field, v0, v1, neighbor, nlist, rules, out, sub0, subn, sub, v, a, s1, s2, s3, spg, rules4jij},
  nlist = GetNeighborList[dist, "OnSite" -> OptionValue["OnSite"]];

  {spg, rules} = Transpose[If[Norm[xyz2RotT[#1][[2]]] == 0, {#1, #2}, ## &[]] & @@@ ({Keys[spg0], TRules}\[Transpose])];
  rules4jij = JijTRules[pos0, spg, rules, Basis, dist];
  
  SetSharedVariable[rules4jij];

  tij = ParallelTable[Table[{s1, s2, s3} = n;
                            sub0 = {Subscript[v_, a_] :> Subscript[v, a, 0, 0, 0]};
                            subn = {Subscript[v_, a_] :> Subscript[v, a, s1, s2, s3]};
                            field = {v0 /. sub0, v1 /. subn};
                            SimplifyCommonFactor[Chop[Expand[Sum[Times@@field /. sub, {sub, rules4jij}]]], 10^-6], {v0, fieldlist[[1]]}, {v1, fieldlist[[2]]}], {n, nlist}, 
                      DistributedContexts -> {"LINVARIANT`INVARIANT`Private`"}];
  tij = Expand[DeleteDuplicates[DeleteCases[Flatten[tij], 0], #1 === -#2 || #1 === #2 &]];

  out = DeleteDuplicates[DeleteCases[GetJijStar[#, -1] &/@ tij, 0], 
                                        (#1 - #2 == 0 || #1 + #2 == 0 || Expand[#1 + I #2] == 0 || Expand[#1 - I #2] == 0) &];
  Return[out]
]

GridSymTransformation[site0_, xyz_] := Module[{siteset0, siteset, rot},
  rot = xyz2RotT[xyz][[1]];
  siteset0 = Join[{{0, 0, 0}}, {0, 0, 0} + IdentityMatrix[3], {1, 1, 1} - IdentityMatrix[3], {{1, 1, 1}}];
  siteset = rot . (# + site0) & /@ siteset0;
  First[Sort[{Total[#], #} & /@ siteset]][[2]]
]

GetJijStar[expr_, disp_] := Module[{vars, nn, v, \[Alpha], \[Beta], a, b, i, j, k, dd, t, terms, tmp, out, ix, iy, iz},
  If[expr == 0. || expr === {}, Return[0]];
  vars = DeleteDuplicates[Cases[Variables[expr],  Subscript[v_, a_, ix_, iy_, iz_] -> Subscript[v, a]]];
  nn = If[disp<0, Max[Max[#] & /@ Abs@Cases[Variables[expr], Subscript[v_, a_, ix_, iy_, iz_] -> {ix, iy, iz}]], disp];

  tmp = Chop@Expand@Sum[expr /. {Subscript[v_, \[Alpha]_, i_, j_, k_] -> Subscript[v, \[Alpha], i + dd[[1]], j + dd[[2]], k + dd[[3]]]}, 
                        {dd, Tuples[Range[-nn, nn], {3}]}];
  terms = If[Head[Expand@tmp] === Plus, Level[Expand@tmp, {1}], {tmp}];

  out = Total@Table[If[!And @@ Map[D[t, {#, 1}] === 0 &, Flatten[vars /. {Subscript[v_, a_] -> Subscript[v, a, 0, 0, 0]}]], t, ## &[]], {t, terms}];

  Return[SimplifyCommonFactor[out, 10^-6]]
]

InvariantEOnSite[expr_] := Module[{X, i, dx, dy, dz},
  expr /. {Subscript[X_, i_, dx_, dy_, dz_] :> Subscript[X, i, ToExpression["ix"] + dx, ToExpression["iy"] + dy, ToExpression["iz"] + dz]}
]

BasisShift[pos0_, spg0_, Basis_] := Module[{cell, xyzStrData, subsites, latt, sites, rot, tran, active, bvector, shift},
  {latt, sites} = pos0;
  cell = GetCellFromGrp[spg0];
  xyzStrData = Which[AssociationQ[spg0], Keys[spg0], ListQ[spg0], spg0];
  shift = Table[{rot, tran} = xyz2RotT[op];
                subsites = {#1, ! Norm[#2] == 0} & @@@ Transpose[{sites\[Transpose][[1]], Partition[bvector, 3]}];
                active = Select[{Chop[rot . #1 - Mod[rot . #1, {1, 1, 1}]], #2} & @@@ subsites, #[[2]] &];
                If[MemberQ[active, {{0, 0, 0}, True}], {0, 0, 0}, Round@Mean[active\[Transpose][[1]]]], {op, xyzStrData}, {bvector, Basis\[Transpose]}];
  Return[shift]
]

InvariantFolding[model_] := Module[{inv, tab, coeff, vars, vsub, data, i, basis, p, tot, head, var2one, svars, invlabel},
  {coeff, inv} = CollectPolynomial[model[[1]].FixDoubleCounting[model[[2]]], {}, "round" -> 1.0 10^-8]\[Transpose];
  vars = Select[Variables[inv], Length@# == 5 &];
  svars = Select[Variables[inv], Length@# == 3 &];
  vsub = {DeleteDuplicates[Subscript[#1, #2, #3, #4, #5] -> Subscript[#1, #2] & @@@ vars],
          DeleteDuplicates[Subscript[#1, #2, #3] -> Subscript[#1, #2, #3] & @@@ svars]};
  var2one = # -> 1 & /@ Variables[inv /. Flatten@vsub];
  data = {If[Head[#2]===Plus, First[#2], #2], #1, Tally[MonomialList[#2] /. Flatten@vsub]} & @@@ Transpose[{coeff, inv}];
  basis = DeleteDuplicates[Flatten[data\[Transpose][[3]], 1]\[Transpose][[1]]];
  tab = Table[Flatten[{data[[i, 1]] /. vsub[[2]], data[[i, 3]][[1, 2]], Table[Total[If[First[#] === p, (data[[i, 2]] Last[#]) /. var2one, 0] & /@ (data[[i, 3]])], {p, basis}]}], {i, Length@data}];
  head = Flatten[{"", "coordinates", basis}];
  tot = Total[tab]; tot[[1]] = "";
  Print[Grid[Join[{Flatten@{"", "", Range[Length[head]-2]}}, {head}, {tot}, tab], Background -> {{1 -> Green, 2 -> Green}, {1 -> Pink, 2 -> Pink, 3 -> Pink}}, ItemSize -> Full]];
]

CollectPolynomial[expr_, TRules_:{}, OptionsPattern[{"SymmetricQ"->False, "chop"->1.0*10^-6, "round"->1.0*10^-16, "sort" -> True}]] := Module[{terms, t, p, data, model, out, v, invariant, s, f, tmp, rules},
  terms = If[Head[Chop[Expand@expr]] === Plus, Level[Chop[Expand@expr], {1}], {expr}];
  rules = Dispatch[#] &/@ TRules;

  data = Table[p = t /. {Subscript[__] -> 1};
               f = If[p==0, p=1; t, Rationalize[t/p]];
               tmp = SimplifyCommonFactor[Chop[Expand[Sum[f /. s, {s, rules}]]], 10^-6];
               invariant = If[TRules === {}, 0, GetJijStar[tmp, -1]];
               {Round[p, OptionValue["round"]], f, If[invariant === 0, {}, Sort@invariant]}, {t, terms}];

  model = Transpose[#] & /@ If[TRules ==={}, 
                               GatherBy[data, Abs@First[#] &],
                               Gather[If[OptionValue["SymmetricQ"], Select[data, !(#[[3]] === {}) &], data], 
                                      (#1[[3]] - #2[[3]] === 0 || #1[[3]] + #2[[3]] === 0) &]];
  
  out = If[Chop@Mean[Abs@#1] < OptionValue["chop"], ##&[], {Sign[First@#1] Chop@Mean[Abs@#1], (Sign[First@#1] Sign[#1]).#2}] & @@@ model;
  Return[If[OptionValue["sort"], ReverseSortBy[out, Abs[#[[1]]]&], SortBy[out, #[[2]]&]]]
]

HModelFolding[hmodel_, grid_, OptionsPattern[{"round"->10^-9}]] := Module[{Hexpr, sub, v, a, dx, dy, dz, ix, iy, iz, h, ngx, ngy, ngz, out},
  {ngx, ngy, ngz} = grid;
  out = Table[Hexpr = Sum[sub = Subscript[v_, a_, dx_, dy_, dz_] :> VarOnGrid[Subscript[v, a], PbcDiff[{ix + dx - 1, iy + dy - 1, iz + dz - 1}, grid]];
                          Expand[h\[Transpose][[1]].FixDoubleCounting[h\[Transpose][[2]]]] /. sub, {iz, ngz}, {iy, ngy}, {ix, ngx}];  
              CollectPolynomial[Chop[Hexpr, OptionValue["round"]]], {h, hmodel}];
  Return[out]
]

FoldHopping[inv_, grid_] := Module[{sub, ngx, ngy, ngz},
  If[ListQ[inv],
     FoldHopping[#, grid] & /@ inv,
     {ngx, ngy, ngz} = grid;
     Sum[sub = Subscript[v_, a_, dx_, dy_, dz_] :> VarOnGrid[Subscript[v, a], PbcDiff[{ix + dx - 1, iy + dy - 1, iz + dz - 1}, grid]];
         First@FixDoubleCounting[inv] /. sub, {iz, ngz}, {iy, ngy}, {ix, ngx}]]
]

PlotPolynomial[model_] := Module[{var2pltdata, var, grid, s, v0, v1, dv1, t, ii},
  var2pltdata[var_, c_] := Module[{v, a, i, j, k, out, data, ss},
    If[ListQ[var], var2pltdata[#, c] & /@ var,
       ss = c Sign[var/First@Variables[var]];
       data = First@Variables[var] /. {Subscript[v_, a_, i_, j_, k_] -> {v, a, {i, j, k}}};
       out = {Text[Style[ToString[data[[1]]], White, Bold, 12], data[[3]] + {1, 1, 1}/5], Red, Arrowheads[0.1], Arrow@Tube[{data[[3]] - 0.25 ss IdentityMatrix[3][[data[[2]]]], data[[3]] + 0.25 ss IdentityMatrix[3][[data[[2]]]]}, 0.05]};
     Return[out]]
  ];
  
  grid = {FaceForm[], EdgeForm[Directive[Black, Dashed]], Cuboid[# - {1, 1, 1}, #]} & /@ Join[{{0, 0, 0}}, {0, 0, 0} + IdentityMatrix[3], {1, 1, 1} - IdentityMatrix[3], {{1, 1, 1}}];
  
  Table[s = Sign@First[t]; 
        v0 = Cases[Variables[t[[2]]], Subscript[v_, a_, 0, 0, 0]];
        dv1 = D[t[[2]], #] /. {Subscript[v_, a_, b_] -> 1} &/@ v0;
        v1 = If[Head@# === Plus, Level[#, {1}], {#}] &/@ dv1;
        Table[Graphics3D[{grid, var2pltdata[v0[[ii]], 1], var2pltdata[v1[[ii]], -s]}, 
                          ImageSize -> {200, 200}, Boxed -> False, Background -> Gray], {ii, Length@v0}], {t, model}]
]

IncludeAcousticMode[pos0_, Basis_, spg0_, tgrp_, vars_, OptionsPattern[{"chop"->0.01}]] := Module[{AcousticBasis, FullBasis, FullOpMat, FullTRules, FullTRules4Jij, LattCar2Dir, NumAtom, latt, sites, i},

  {latt, sites} = pos0;
  NumAtom = Length[sites];

  LattCar2Dir = Inverse[Transpose[latt]];

  AcousticBasis = Transpose[Flatten[Chop@Table[LattCar2Dir.(Partition[Normalize[#], 3][[i]]), {i, NumAtom}]] & /@ (Flatten[Table[1.0 IdentityMatrix[3], {NumAtom}], 1]\[Transpose])]; 

  FullBasis = Transpose@Join[Basis\[Transpose], AcousticBasis\[Transpose]]; 
  
  FullOpMat = Chop[GetMatrixRep[spg0, tgrp, pos0, FullBasis, Range[Length[Flatten@vars] + 3], "disp"], OptionValue["chop"]]; 
  
  FullTRules = Rationalize@Chop[Join[#1, #2] & @@@ ({GetIsoTransformRules[FullOpMat, "disp"], GetIsoStrainTransformRules[latt, spg0]}\[Transpose])];
  
  FullTRules4Jij = JijTRules[pos0, spg0, FullTRules, Basis, 1];

  Return[{FullBasis, FullOpMat, FullTRules, FullTRules4Jij}]
]

BasisInSupercell[pos0_, Basis_, vars_, NSC_, SymmetricQ_ : True, OptionsPattern[{"round" -> 1.0*10^-6, "rotate"->False}]] := Module[{Nij, Nx, Ny, Nz, latt, sites, ia, i, j, k, s, g, ppos, spos, p2s, s2p, shift, multiplicity, BasisTemplates, SuperBasis, bm, cell0sites, BoolQ, NumBasis, LengthBasis, SuperVars, NumPpos, NumSpos, NeighborList, LattCar2Dir},
  {LengthBasis, NumBasis} = Dimensions@Basis;
  Nij = DiagonalMatrix[NSC];
  {latt, sites} = pos0;
  LattCar2Dir = Inverse@Transpose[latt];
  {Nx, Ny, Nz} = NSC;
  NeighborList = If[OptionValue["rotate"],
                    Flatten[Table[{k, j, i} - {1, 1, 1}, {i, Nz}, {j, Ny}, {k, Nx}], 2],
                    Flatten[Table[{i, j, k} - {1, 1, 1}, {i, Nx}, {j, Ny}, {k, Nz}], 2]];

   ppos = pos0;
   spos = If[OptionValue["rotate"],
             {Nij.latt, Flatten[Table[{#[[1]] + {k - 1, j - 1, i - 1}, #[[2]]}, {i, Nz}, {j, Ny}, {k, Nx}] & /@ sites, 3]},
             {Nij.latt, Flatten[Table[{#[[1]] + {i - 1, j - 1, k - 1}, #[[2]]}, {i, Nx}, {j, Ny}, {k, Nz}] & /@ sites, 3]}];

  NumPpos = Length[ppos[[2]]];
  NumSpos = Length[spos[[2]]];
  p2s = Association[Thread[NeighborList -> Table[Association[Table[i -> First@First@Position[Rationalize[Norm[ppos[[2]][[i]][[1]] + shift - #1]] & @@@ (spos[[2]]), 0], {i, NumPpos}]], {shift, NeighborList}]]];
  s2p = Association@Flatten[Table[p2s[{ix - 1, iy - 1, iz - 1}][i] -> {{ix - 1, iy - 1, iz - 1}, i}, {iz, Nz}, {iy, Ny}, {ix, Nx}, {i, NumPpos}], 3];

  cell0sites = If[SymmetricQ, Table[DeleteDuplicates[Chop[# ({1, 1, 1} - Sign[Abs[Chop@s]]) + Chop[s] Sign[Abs[Chop@s]]] & /@ Tuples[{0, 1}, {3}]], {s, sites\[Transpose][[1]]}], {#} & /@ (sites\[Transpose][[1]])];
  multiplicity = If[SymmetricQ, 1/Power[2, Count[Chop[Mod[#1, 1]], 0]] & @@@ (spos[[2]]), 1];


  BasisTemplates = Table[multiplicity Partition[Flatten[ConstantArray[latt\[Transpose].#, Times @@ NSC] & /@ Partition[Basis[[;; , i]], 3]], 3], {i, NumBasis}];

  SuperBasis = Table[bm = Table[{g, ia} = s2p[s];
                                BoolQ = Or @@ (Chop@Norm[spos[[2]][[s]][[1]] - Mod[# + shift, NSC]] == 0 & /@ (cell0sites[[ia]]));
                                If[BoolQ, {1, 1, 1}, {0, 0, 0}], {s, NumSpos}];
                     Table[Flatten[bm  BasisTemplates[[i]]], {i, NumBasis}], {shift, Keys@p2s}];
  SuperBasis = Flatten[Table[Chop[LattCar2Dir.(Partition[#, 3][[i]])], {i, NumSpos}]] &/@ Flatten[SuperBasis, 1];

  SuperVars = Flatten[Table[s = Mod[{1, 1, 1} + shift, NSC] - {1, 1, 1};
                            s = PbcDiff[shift, NSC];
                            Subscript@@Join[If[Level[#, Infinity]==={},{#},Level[#, Infinity]], s] & /@ Flatten[vars], {shift, Keys[p2s]}], 1];

  Return[{SuperBasis\[Transpose], SuperVars}]
]


(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
