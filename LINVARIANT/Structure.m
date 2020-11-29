BeginPackage["LINVARIANT`Structure`"]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
ModCell                      ::usage = "ModCell[xyz, cell]"
exportXYZ                    ::usage = "exportXYZ[filepath_, {comment_, vertices_, atomcoordinates_}]"
cif2mcif                     ::usage = "cif2mcif[id, IsoMatrix, pos0]"
SymmetryOpOrigin             ::usage = "SymmetryOpOrigin[file, set]"
SymmetryOpVectorField        ::usage = "SymmetryOpVectorField[file, set]"
GetLatticeVectors            ::usage = "GetLatticeVectors[dat]"
Lattice2Symbolic             ::usage = "Lattice2Symbolic[latt]"
ImportIsodistortCIF          ::usage = "ImportIsodistortCIF[\*StyleBox[\"filename\",\"TI\"]] opens .cif file and imports relevant data."
PosMatchTo                   ::usage = "PosMatchTo[spos, dpos, tol]"
GetStrainTensor              ::usage = "GetStrainTensor[parent, sub]"
GridPbc                      ::usage = "GridPbc[ixyz, Lxyz]"
GridNeighbors                ::usage = "GridNeighbors[r0, MeshDim]"
PbcDiff                      ::usage = "PbcDiff[diff]"
GetUpTo3rdNeighbors          ::usage = "GetUpTo3rdNeighbors[]"
MakeSuperCell                ::usage = "MakeSuperCell[Crys, {Nx, Ny, Nz}]"
GetBZPath                    ::usage = "GetBZPath[struc]"
GetBZHighSymK                ::usage = "GetBZHighSymK[struc]"
AxialRotation                ::usage = "AxialRotation[dir, pos]"
PlotVectorField              ::usage = "PlotVectorField[latt, field]"
xright                       ::usage = "xright[x0, nn, NN]"
xleft                        ::usage = "xleft[x0, nn, NN]"
GetKpath                     ::usage = "GetKpath[klist, div]"
LatticeFromLenAng            ::usage = "LatticeFromLenAng[VecAng]"
SimplifyElementSymbol        ::usage = "SimplifyElementSymbol[ele]"
PlotCrystal                  ::usage = "PlotGrid[latt, pos, dim]"
DistMatrix                   ::usage = "DistMatrix[pos1, pos2]"
pos2index                    ::usage = "pos2index[atoms, pos]"
AtomicFormFactor             ::usage = "AtomicFormFactor[atom, q]"
XRayDiffraction              ::usage = "XRayDiffraction[q, latt, Atoms]"
GetBZKList                   ::usage = "GetBZKList[spg0, kinv]"
PlotStruct                   ::usage = "PlotStruct[pos]"

(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)
Options[ImportIsodistortCIF]    = {Fractional->False, CorrectLabels->True, Tolerance->10^-6}

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
pos2index[latt_, atoms_, pos_] := If[Length[Dimensions[pos]] > 1, 
                              pos2index[atoms, #] &/@ pos,
                              First@First@Position[Rationalize[Chop[Flatten[DistMatrix[latt, atoms\[Transpose][[1]], {pos}, {1,1,1}]]]], 0]]
ModCell[xyz_, cell_] := Module[{},
  Mod[#1, #2] &@@@ ({xyz, cell}\[Transpose])
]

exportXYZ[filepath_, {comment_, vertices_, atomcoordinates_}] := Module[{outputfilestream, outputstring},
  If[Or[Not@StringQ[filepath], Not@StringQ[comment]], Message[exportXYZ::usage]; Abort[]];
  If[Length@vertices != First@Dimensions@atomcoordinates, Message[exportXYZ::mismatch]; Abort[]];
  outputstring = ToString@Length[vertices] <> "\n" <> comment <> "\n" <> MapThread[(#1 <> "\t" <> ExportString[{#2}, "TSV"] <> "" &), {vertices, atomcoordinates}];
  outputfilestream = OpenWrite[filepath];
  WriteString[outputfilestream, outputstring];
  Close[outputfilestream];
Return[outputstring]]

GetIRvector[id_, IsoMatrix_, pos_] := Module[{IRvector},
  IRvector = {If[IntegerQ[id], # & /@ Partition[IsoMatrix[[;; , id]] // Normal, 3], Print["mode not exist!"]], # & /@ (pos\[Transpose][[2]])}\[Transpose];
  Return[IRvector]
]

cif2mcif[file_, id_, IsoMatrix_, pos0_] := Module[{CifTagVec, mcif, cif, vectors},
  CifTagVec = {{}, {"loop_"}, {"_atom_site_moment.label"}, {"_atom_site_moment.crystalaxis_x"}, {"_atom_site_moment.crystalaxis_y"}, {"_atom_site_moment.crystalaxis_z"}, {"_atom_site_moment.symmform"}};
  cif = Import[file, "Table"];
  vectors = {First@StringCases[#2, RegularExpression["[[:upper:]][[:lower:]]*"]] <> If[dd = StringCases[#2, RegularExpression["\\d+"]]; dd === {}, Unevaluated[Sequence[]], dd], #1[[1]], #1[[2]], #1[[3]], "mx,my,mz"} & @@@ GetIRvector[id, IsoMatrix, pos0];
  mcif = Join[cif, CifTagVec, vectors];
  Return[mcif]
  ]

SymmetryOpOrigin[file_, set_] := Module[{xyzStrData, xyzExpData, newsets, op, posmap, diff, i, j},
  xyzStrData = CifImportOperators[file];
  xyzExpData = Table[ToExpression["{" <> xyzStrData[[i]] <> "}"], {i, Length[xyzStrData]}];
  newsets = Table[op /. Thread[ToExpression[{"x", "y", "z"}] -> #1] & @@@set, {op, xyzExpData}];
  posmap = Table[Flatten[Table[diff = Mod[newsets[[op]][[i]] - set[[j]][[1]], 1]; If[diff == {0., 0., 0.}, {j, i, set\[Transpose][[2]][[j]]},Unevaluated[Sequence[]]], {i, Length@set}, {j, Length@set}], 1], {op, Length@xyzExpData}];
  (*posmap = Table[SortBy[op, First], {op, posmap}];*)
  newsets = Table[{newsets[[op]], Transpose[posmap[[op]]][[3]]}\[Transpose], {op, Length@xyzExpData}];
  Return[newsets]
]

(*DistMatrixOld = Compile[{{pos1, _Real, 2}, {pos2, _Real, 2}}, Module[{i, j},
   Return[Table[Norm[Mod[Round[pos1[[i]] - pos2[[j]], 10^-8], 1]], {i, Length@pos1}, {j, Length@pos2}]]
]]

DistMatrix = Compile[{{pos1, _Real, 2}, {pos2, _Real, 2}}, Module[{i, j},
   Return[Table[Norm[Round[Which[# > 0.5, # - 1, # <= -0.5, # + 1, True, #] & /@ (pos1[[i]] - pos2[[j]]), 10^-8]^2], {i, Length@pos1}, {j, Length@pos2}]]
]]*)

DistMatrix = Compile[{{latt, _Real, 2}, {pos1, _Real, 2}, {pos2, _Real, 2}, {cell, _Integer, 1}}, Module[{i, j, n, diff},
  Return[Table[diff = pos1[[i]] - pos2[[j]];
               Norm[Round[latt.Table[Which[diff[[n]] > cell[[n]]/2, diff[[n]] - cell[[n]], diff[[n]] < -cell[[n]]/2, diff[[n]] + cell[[n]], True, diff[[n]]], {n, 3}], 10^-6]], {i, Length@pos1}, {j, Length@pos2}]]
]]

(*DistMatrix[pos1_, pos2_, cell_:{1, 1, 1}] := Module[{i, j},
  Return[Table[Norm[Round[PbcDiff[pos1[[i]] - pos2[[j]], cell], 10^-8]], {i, Length@pos1}, {j, Length@pos2}]]
]
*)

PosMatchTo[latt_, spos_, dpos_, OptionsPattern["shift"->False]] := Module[{difftable, posmap, i, newpos},
  difftable = DistMatrix[latt, spos, dpos, {1,1,1}];
  posmap = First@First@Position[#, x_ /; TrueQ[x == Min[#]]] & /@ difftable;
  newpos = If[OptionValue["shift"], PbcDiff[#]&/@(Table[dpos[[posmap[[i]]]], {i, Length@spos}] - spos) + spos,Table[dpos[[posmap[[i]]]], {i, Length@spos}]];
  Return[{posmap, newpos}]
]

SymmetryOpVectorField[grp_, pos_, vec_, ftype_] := Block[{originshift, xyzStrData, xyzRotTranData, xyzTranslation, xyzRotData, field, newpos, newvec, difftable, diff, posmap, i, j},
  originshift = {0., 0., 0.};
  xyzStrData = Keys[grp];
  xyzRotTranData = Table[ToExpression["{" <> xyzStrData[[i]] <> "}"], {i, Length[xyzStrData]}];
  xyzTranslation = xyzRotTranData  /. {ToExpression["x"] -> 0, ToExpression["y"] -> 0, ToExpression["z"] -> 0} ;
  xyzRotData = xyzRotTranData - xyzTranslation;

  newvec = Which[ftype=="disp",
                 Table[{Det[Expr2Rot[op]]^2*(op /. Thread[ToExpression[{"x", "y", "z"}] -> #1]), #2} & @@@ vec, {op, xyzRotData}],
                 ftype=="spin",
                 Table[{Det[Expr2Rot[op]]*(op /. Thread[ToExpression[{"x", "y", "z"}] -> #1]), #2} & @@@ vec, {op, xyzRotData}]];
  newpos = Table[{(op /. Thread[ToExpression[{"x", "y", "z"}] -> (pos[[i]][[1]]+originshift)]), i}, {op, xyzRotTranData}, {i, Length@pos}];
  
  (*difftable = ParallelTable[DistMatrix[#+originshift&/@(pos\[Transpose][[1]]), newpos[[op]]\[Transpose][[1]]], {op, Length@xyzRotTranData}, DistributedContexts -> {"LINVARIANT`Structure`Private`"}];*)
  difftable = Table[DistMatrix[#+originshift&/@(pos\[Transpose][[1]]), newpos[[op]]\[Transpose][[1]], {1,1,1}], {op, Length@xyzRotTranData}];
  posmap = Table[Position[difftable[[op]], x_ /; TrueQ[x == 0]], {op, Length@xyzRotTranData}];
  field = Table[Table[{newvec[[i]][[posmap[[i]][[j]][[2]]]][[1]], pos\[Transpose][[2]][[j]]}, {j, Length@pos}], {i, Length@posmap}];
  Return[field]
]

SymmetryOpVectorFieldCif[file_, pos_, vec_, OptionsPattern["spin" -> False]] := Block[{originshift, xyzStrData, xyzRotTranData, xyzTranslation, xyzRotData, field, newpos, newvec, difftable, diff, posmap, i, j, axial},
  If[! FileExistsQ[file], Print["Error:" <> file <> " was not found in the working directory!"]; Abort[]];
  CifData = Import[file, "string"] <> "\n";
  originshift = ToExpression[StringReplace[First[StringCases[CifData,RegularExpression["origin=\\W\\W*\\d\\.*\\d*,\\W*\\d\\.*\\d*,\\W*\\d\\.*\\d*\\W"]]], {"origin=" -> "", "(" -> "{", ")" -> "}"}]];
  (*---fix for some strange behaviour---*)
  xyzStrData = CifImportOperators[file];
  xyzRotTranData = Table[ToExpression["{" <> xyzStrData[[i]] <> "}"], {i, Length[xyzStrData]}];
  xyzTranslation = xyzRotTranData  /. {ToExpression["x"] -> 0, ToExpression["y"] -> 0, ToExpression["z"] -> 0} ;
  xyzRotData = xyzRotTranData - xyzTranslation;

  axial = If[OptionValue["spin"], 1, 2];

  newvec = Table[{Det[Expr2Rot[op]]^axial*N[op /. Thread[ToExpression[{"x", "y", "z"}] -> #1]], #2} & @@@ vec, {op, xyzRotData}];
  newpos = Table[{N[op /. Thread[ToExpression[{"x", "y", "z"}] -> (pos[[i]][[1]]+originshift)]], i}, {op, xyzRotTranData}, {i, Length@pos}];

  (*difftable = ParallelTable[DistMatrix[#+originshift&/@(pos\[Transpose][[1]]), newpos[[op]]\[Transpose][[1]]], {op, Length@xyzRotTranData}, DistributedContexts ->   {"LINVARIANT`Structure`Private`"}];*)
  difftable = Table[DistMatrix[#+originshift&/@(pos\[Transpose][[1]]), newpos[[op]]\[Transpose][[1]], {1,1,1}], {op, Length@xyzRotTranData}];
  posmap = Table[Position[difftable[[op]], x_ /; TrueQ[x == 0]], {op, Length@xyzRotTranData}];
  field = Table[Table[{newvec[[i]][[posmap[[i]][[j]][[2]]]][[1]], pos\[Transpose][[2]][[j]]}, {j, Length@pos}], {i, Length@posmap}];
  Return[field]
]

GetUpTo3rdNeighbors[OptionsPattern[{"OnSite" -> False}]] := Module[{FirstNeighborList, SecondNeighborList, ThirdNeighborList},
  FirstNeighborList = Flatten[{-1, 1}\[TensorProduct]IdentityMatrix[3], 1];
  SecondNeighborList = DeleteCases[Total[#] & /@ Subsets[FirstNeighborList, {2}], {0, 0, 0}];
  ThirdNeighborList = DeleteCases[Total[#] & /@ Subsets[FirstNeighborList, {3}], a_ /; MemberQ[FirstNeighborList, a]];
  Return[If[OptionValue["OnSite"], 
            {{{0,0,0}}, FirstNeighborList, SecondNeighborList, ThirdNeighborList}, 
            {FirstNeighborList, SecondNeighborList, ThirdNeighborList}]
  ]
]

GetLatticeVectors[dat_] := Module[{a1, b1, c1, \[Alpha]1, \[Beta]1, \[Gamma]1, bt},
     {a1, b1, c1, \[Alpha]1, \[Beta]1, \[Gamma]1} = dat;
     bt = {{a1, 0, 0},
           {b1 Cos[\[Gamma]1], b1 Sin[\[Gamma]1], 0},
           {c1 Cos[\[Beta]1], c1 (Cos[\[Alpha]1] - Cos[\[Beta]1] Cos[\[Gamma]1])/Sin[\[Gamma]1],
            c1 Sqrt[1 - (Cos[\[Alpha]1]^2 + Cos[\[Beta]1]^2 + Cos[\[Gamma]1]^2) + 2 Cos[\[Alpha]1] Cos[\[Beta]1] Cos[\[Gamma]1]]/Sin[\[Gamma]1]}};
     Return[bt]
]

Lattice2Symbolic[latt_] := Module[{G, aa, bb, cc, \[Alpha], \[Beta], \[Gamma], LatticeVectors, LatticeScales},
  {aa, bb, cc} = Rationalize[Norm[#], 10^-9] & /@ latt;
  G = (latt.(latt\[Transpose]));
  {\[Alpha], \[Beta], \[Gamma]} = Rationalize[ArcCos[G[[#1, #2]]/(Norm[latt[[#1]]] Norm[latt[[#2]]])]/Pi] Pi & @@@ {{2, 3}, {1, 3}, {1, 2}};
  LatticeVectors = GetLatticeVectors[{aa, bb, cc, \[Alpha], \[Beta], \[Gamma]}];
  LatticeScales = Table[ToExpression[FromCharacterCode[96 + i]] -> Simplify@Norm[LatticeVectors[[i]]], {i, Length[LatticeVectors]}];
  LatticeVectors = Table[ToExpression[FromCharacterCode[96 + i]]*Simplify@Normalize[LatticeVectors[[i]]], {i, Length[LatticeVectors]}];
  {x, y, z} = ToExpression["{a, b, c}"] /. LatticeScales;
  If[PossibleZeroQ[y - z], LatticeScales = Delete[LatticeScales, 3];LatticeVectors = LatticeVectors /. {ToExpression["c"] -> ToExpression["b"]}, 
     If[PossibleZeroQ[x - z], LatticeScales = Delete[LatticeScales, 3];
     LatticeVectors = LatticeVectors /. {ToExpression["c"] -> ToExpression["a"]}]];
  If[PossibleZeroQ[x - y], LatticeScales = Delete[LatticeScales, 2];LatticeVectors = LatticeVectors /. {ToExpression["b"] -> ToExpression["a"]}];
  Clear[x, y, z];
  (*LatticeVectors= LatticeVectors /. {ToExpression["a"] -> ToExpression["a1"], ToExpression["b"] -> ToExpression["a2"], ToExpression["c"] -> ToExpression["a3"]};*)
  Return[LatticeVectors]
]

ImportIsodistortCIF[file_,OptionsPattern[]] := Module[{CifData, CifFlags, xyzName, ElemSym, ElemName, SpgName, xyzStrData, xyzExpData, xyzTranslation, Atoms, EveryAtomPos, LatticeVectors, LatticeScales, sol, i, j, k, failed, x,y,z, SpgNumber, aa,bb,cc,alpha,beta,gamma,FracOutText,Wyckoff,IsoDim,IsoDispModeRow,IsoDispModeCol,IsoDispModeVal,IsoDispModeMatrix, IsoDispModesID,IsoDispModesNam,IsoDispModesAmp,IsoDispModes,IsoDispModeNormVal,IsoDeltaCoordID,IsoDeltaCoordLabel,IsoDeltaCoordVal,IsoDeltaCoord, IsoCoordLabel,IsoCoordFormula,IsoCoord, IsoStrainModesID,IsoStrainModesNam,IsoStrainModesAmp,IsoStrainModes, IsoStrainModeRow, IsoStrainModeCol, IsoStrainModeVal},

    If[!FileExistsQ[file], Print["Error:" <> file <> " was not found in the working directory!"]; Abort[] ];
	CifData = Import[file, "string"] <> "\n"; (*--- fix for some strange behaviour ---*)
    CifData = StringReplace[CifData,Thread[DeleteDuplicates[StringCases[CifData, RegularExpression[";"]]] -> ","]];
    CifData = StringReplace[CifData,Thread[DeleteDuplicates[StringCases[CifData,RegularExpression["[[:upper:]].{3,6}(\[)\\d,\\d,\\d(\])"]]] -> ""]];
    IsoDispModes = StringCases[CifData, RegularExpression["([[:upper:]][[:lower:]]*\\W*|[[:upper:]]*\\d*[[:lower:]]*)_\\d+"]]; 
    CifData = StringReplace[CifData, Thread[IsoDispModes -> StringReplace[IsoDispModes, {"*" -> "conj", "'" -> "Prime", "_" -> ""}]]];

	CifData = ImportString[CifData, "CIF"];
	CifFlags = Table[Level[CifData[[i]], 1][[1]], {i, Length[CifData]}];
	ElemSym = If[Position[CifFlags, "_chemical_formula_sum"] != {},"_chemical_formula_sum" /. CifData,""];
	ElemName = If[Position[CifFlags, "_chemical_name_mineral"] != {},"_chemical_name_mineral" /. CifData,""];

    IsoDim = If[Position[CifFlags, "_iso_displacivemode_number"] != {}, "_iso_displacivemode_number" /. CifData, 0];
    IsoDispModeRow = If[Position[CifFlags, "_iso_displacivemodematrix_row"] != {}, "_iso_displacivemodematrix_row" /. CifData, {"None"}];
    IsoDispModeCol = If[Position[CifFlags, "_iso_displacivemodematrix_col"] != {}, "_iso_displacivemodematrix_col" /. CifData, {"None"}];
    IsoDispModeVal = If[Position[CifFlags, "_iso_displacivemodematrix_value"] != {}, Rationalize["_iso_displacivemodematrix_value" /. CifData], {"None"}];
    (*IsoDispModeMatrix = Normal[SparseArray[Thread[Transpose[{IsoDispModeRow, IsoDispModeCol}] -> IsoDispModeVal], {IsoDim, IsoDim}]];*)

    IsoDispModesID = If[Position[CifFlags, "_iso_displacivemode_ID"] != {}, "_iso_displacivemode_ID" /. CifData, {"None"}];
    IsoDispModesNam = If[Position[CifFlags, "_iso_displacivemode_label"] != {}, "_iso_displacivemode_label" /. CifData, {"None"}];
    IsoDispModesAmp = If[Position[CifFlags, "_iso_displacivemode_value"] != {}, "_iso_displacivemode_value" /. CifData, {"None"}];
    IsoDispModeNormVal = If[Position[CifFlags, "_iso_displacivemodenorm_value"] != {},"_iso_displacivemodenorm_value" /. CifData, {"None"}];
    IsoDispModes = Transpose[{IsoDispModesID, IsoDispModesNam, IsoDispModeNormVal, IsoDispModesAmp}];

    IsoStrainModeRow = If[Position[CifFlags, "_iso_strainmodematrix_row"] != {}, "_iso_strainmodematrix_row" /. CifData, {"None"}];
    IsoStrainModeCol = If[Position[CifFlags, "_iso_strainmodematrix_col"] != {}, "_iso_strainmodematrix_col" /. CifData, {"None"}];
    IsoStrainModeVal = If[Position[CifFlags, "_iso_strainmodematrix_value"] != {}, Rationalize["_iso_strainmodematrix_value" /. CifData], {"None"}];

    IsoStrainModesID = If[Position[CifFlags, "_iso_strainmode_ID"] != {}, "_iso_strainmode_ID" /. CifData, {"None"}];
    IsoStrainModesNam = If[Position[CifFlags, "_iso_strainmode_label"] != {}, "_iso_strainmode_label" /. CifData, {"None"}];
    IsoStrainModesAmp = If[Position[CifFlags, "_iso_strainmode_value"] != {}, "_iso_strainmode_value" /. CifData, {"None"}];
    IsoStrainModes = Transpose[{IsoStrainModesID, IsoStrainModesNam, IsoStrainModesAmp}];

    IsoDeltaCoordID = If[Position[CifFlags, "_iso_deltacoordinate_ID"] != {}, "_iso_deltacoordinate_ID" /. CifData, {"None"}];
    IsoDeltaCoordLabel = If[Position[CifFlags, "_iso_deltacoordinate_label"] != {}, "_iso_deltacoordinate_label" /. CifData, {"None"}];
    IsoDeltaCoordVal = If[Position[CifFlags, "_iso_deltacoordinate_value"] != {}, "_iso_deltacoordinate_value" /. CifData, {"None"}];
    IsoDeltaCoord = Transpose[{IsoDeltaCoordID, IsoDeltaCoordLabel, IsoDeltaCoordVal}];

    IsoCoordLabel = If[Position[CifFlags, "_iso_coordinate_label"] != {}, "_iso_coordinate_label" /. CifData, {"None"}];
    IsoCoordFormula = If[Position[CifFlags, "_iso_coordinate_formula"] != {}, "_iso_coordinate_formula" /. CifData, {"None"}];
    IsoCoord = Transpose[{IsoCoordLabel, IsoCoordFormula}];

	(*--- get lattice ---*)
	aa = If[Position[CifFlags, "_cell_length_a"] != {},Rationalize["_cell_length_a" /. CifData], Print["lattice error"]; Abort[]];
	bb = If[Position[CifFlags, "_cell_length_b"] != {},Rationalize["_cell_length_b" /. CifData], Print["lattice error"]; Abort[]];
	cc = If[Position[CifFlags, "_cell_length_c"] != {},Rationalize["_cell_length_c" /. CifData], Print["lattice error"]; Abort[]];
	alpha = If[Position[CifFlags, "_cell_angle_alpha"] != {},Rationalize["_cell_angle_alpha" /. CifData] \[Degree], Print["lattice error"]; Abort[]];
	beta = If[Position[CifFlags, "_cell_angle_beta"] != {},Rationalize["_cell_angle_beta" /. CifData] \[Degree], Print["lattice error"]; Abort[]];
	gamma = If[Position[CifFlags, "_cell_angle_gamma"] != {},Rationalize["_cell_angle_gamma" /. CifData] \[Degree], Print["lattice error"]; Abort[]];
	LatticeVectors = GetLatticeVectors[{aa,bb,cc,alpha,beta,gamma}];
	LatticeScales = Table[ToExpression[FromCharacterCode[96+i]] -> Simplify@Norm[LatticeVectors[[i]]] ,{i,Length[LatticeVectors]}];
	LatticeVectors = Table[ToExpression[FromCharacterCode[96+i]]* Simplify@Normalize[LatticeVectors[[i]]] ,{i,Length[LatticeVectors]}];
    If[Position[CifFlags,"_atom_site_fract_x" | "_atom_site_fract_y" |"_atom_site_fract_z"] != {},
       Wyckoff = {ToExpression[#1], #2} &@@@ Transpose[{Transpose[({"_atom_site_fract_x","_atom_site_fract_y", "_atom_site_fract_z"} /. CifData)],"_atom_site_label" /. CifData}];,
       If[Position[CifFlags,"_atom_site_Cartn_x" | "_atom_site_Cartn_y" |"_atom_site_Cartn_z"] != {},
          Print["cartesian coords found, not implemented"];Abort[], 
          Print["no atom_side_fract"];Abort[]
       ];
    ];
	(*--- get multiplicity and transform to cartesian ---*)
    xyzName = Part[CifFlags, Flatten[Position[StringCount[CifFlags, "xyz"], 1]]];
    xyzStrData = Flatten[xyzName /. CifData];
    xyzExpData = Table[ToExpression["{" <> xyzStrData[[i]] <> "}"], {i, Length[xyzStrData]}];

	(*--- apply multiplicity ---*)
    EveryAtomPos = Flatten[Table[Sort[Union[{Mod[#, 1], Wyckoff[[i, 2]]} & /@ (xyzExpData /. Thread[ToExpression[{"x", "y", "z"}] -> Wyckoff[[i, 1]]]), SameTest->(Norm[Chop[#1[[1]]-#2[[1]]]]==0.&)], Total[#1[[1]]] > Total[#2[[1]]] &], {i, Length[("_atom_site_label" /. CifData)]}], 1];
    Atoms = {#[[1]], ElementLabel = Characters[#[[2]]]; ToUpperCase[ElementLabel[[1]]] <> If[Length[ElementLabel] > 1 && LetterQ[ElementLabel[[2]]], ElementLabel[[2]], ""]} & /@ EveryAtomPos;

	(*--- filter equal length ---*)
	{x,y,z} = {a,b,c} /. LatticeScales;
	If[ PossibleZeroQ[y-z],
		LatticeScales = Delete[LatticeScales,3];
		Atoms = Atoms /.{ c -> b};
		LatticeVectors = LatticeVectors /.{c->b},
	  If[ PossibleZeroQ[x-z],
		LatticeScales = Delete[LatticeScales,3];
		Atoms = Atoms /.{c->a};
		LatticeVectors = LatticeVectors /.{c->a}]
	];
	If[ PossibleZeroQ[x-y],
		LatticeScales = Delete[LatticeScales,2];
		Atoms = Atoms /.{b->a};
		LatticeVectors = LatticeVectors /.{b->a}];
	Clear[x,y,z];

	(*--- space group names ---*)
	SpgNumber = If[Position[CifFlags, "_space_group_IT_number"] != {},"_space_group_IT_number" /. CifData,0];
	If[SpgNumber!=0 && IntegerQ[SpgNumber],
		SpgName = "None"
	,
		If[Position[CifFlags, "_symmetry_space_group_name_H-M"] != {},
			SpgName = "_symmetry_space_group_name_H-M" /. CifData;
			SpgName = StringJoin[DeleteCases[Characters[SpgName], " "]];
			SpgName = StringReplace[SpgName, {"m3m" -> "m-3m"}];
			SpgName = BracketingBar[ToExpression[SpgName]],
			Print["space group name error"];
			SpgName = ""
		]
	];

    FracOutText="fractional";

	If[ValueQ[failed], Print[failed]];
	{{ElemSym, ElemName}, "","", SpgName, SpgNumber, LatticeVectors, Atoms, LatticeScales, FracOutText, Wyckoff, IsoDim, IsoDispModes, {IsoDispModeRow, IsoDispModeCol, IsoDispModeVal}, IsoStrainModes, {IsoStrainModeRow, IsoStrainModeCol, IsoStrainModeVal}, IsoDeltaCoord, IsoCoord}
]

GetStrainTensor[parent_, sub_, OptionsPattern["iso"->True]] := Module[{e, strain},
  strain= If[OptionValue["iso"],
             e = N[sub[[6]] /. sub[[8]]].Inverse[N[parent[[6]] /. parent[[8]]]] - IdentityMatrix[3];
             1/2 (e + Transpose[e] + Transpose[e].e),
             e = sub.Inverse[parent] - IdentityMatrix[3];
             1/2 (e + Transpose[e] + Transpose[e].e)]
]

GridPbc[ixyz_, Lxyz_] := Module[{}, Mod[ixyz - 1, Lxyz] + 1]

GridNeighbors[r0_, MeshDim_] := Module[{r00, list, FirstNeighborList, SecondNeighborList, ThirdNeighborList, Neighbours, a},
  FirstNeighborList = Flatten[{-1, 1}\[TensorProduct]IdentityMatrix[3], 1];
  SecondNeighborList = DeleteCases[Total[#] & /@ Subsets[FirstNeighborList, {2}], {0, 0, 0}];
  ThirdNeighborList = DeleteCases[Total[#] & /@ Subsets[FirstNeighborList, {3}], a_ /; MemberQ[FirstNeighborList, a]];
  r00 = GridPbc[r0, MeshDim];
  Neighbours = Table[GridPbc[r00 + #, MeshDim] & /@ list, {list, {FirstNeighborList, SecondNeighborList, ThirdNeighborList}}];
  Return[Neighbours]
]

PbcDiff[diff_, cell_:{1, 1, 1}] := Module[{halflattice}, 
  halflattice = cell/2;
  Which[-#2 <= #1 < #2, #1, #1 >= #2, #1 - 2 #2, #1 < -#2, #1 + 2 #2] &@@@ ({diff, halflattice}\[Transpose])
]

MakeSuperCell[Crys_, mode_, k_, {Nx_, Ny_, Nz_}] := Module[{R, pos, ix, iy, iz, spos, sR},
  {R, pos} = Crys;
  spos = If[mode != {}, 
            Flatten[Table[{Chop[(#1 + {ix, iy, iz})/{Nx, Ny, Nz} + Re[#3 Exp[I 2 Pi {ix, iy, iz}.k]]], #2}, {ix, 0, Nx - 1}, {iy, 0,   Ny - 1}, {iz, 0, Nz - 1}] & @@@ (Append[pos\[Transpose], mode]\[Transpose]), 3], 
            Flatten[Table[{Chop[(#1 + {ix, iy, iz})/{Nx, Ny, Nz}], #2}, {ix, 0, Nx - 1}, {iy, 0, Ny - 1}, {iz, 0, Nz - 1}] & @@@ pos, 3]];
  sR = DiagonalMatrix[{Nx, Ny, Nz}].R;
  Return[{sR, spos}]
]

GetBZPath[struc_] := Module[{path,plist},
  plist={"fcc","bcc","sc","Honeycomb","SquareLattice"};
  If[struc=="Help",Print["The following paths are implemented: ",plist];Return[],None];
  If[struc == "fcc", path = {{{0, 0, 0}, {0, 0, 1}, {0, 1/2, 1}, {1/2, 1/2, 1/2}, {0, 0, 0},{3/4,3/4,0}, {1, 1, 0}}, {"\[CapitalGamma]", "X", "W", "L", "\[CapitalGamma]","K","X"}}, None];
  If[struc == "bcc", path = {{{0, 0, 0}, {0, 0, 1}, {0, 1/2, 1/2}, {0, 0, 0}, {1/2, 1/2, 1/2}, {0, 1/2, 1/2}}, {"\[CapitalGamma]", "H", "N", "\[CapitalGamma]", "P", "N"}}, None];
  If[struc == "sc", path = {{{0, 0, 0}, {0, 0, 1/2}, {0, 1/2, 1/2}, {1/2, 1/2, 1/2}, {0, 0, 0}, {0, 1/2, 1/2}, {1/2, 1/2, 1/2}}, {"\[CapitalGamma]", "X", "M", "R", "\[CapitalGamma]", "M", "R"}}, None];
  If[struc == "Honeycomb", path = {{{-2/3/Sqrt[3], 0, 0}, {0, 0, 0}, {1/2/Sqrt[3], -1/6, 0}, {2/3/Sqrt[3], 0, 0}}, {"K'", "\[CapitalGamma]", "M", "K"}}, None];
  If[struc == "SquareLattice", path = {{{0, 0, 0}, {1/2, 0, 0}, {1/2, 1/2, 0}, {0, 0, 0}}, {"\[CapitalGamma]", "X", "M", "\[CapitalGamma]"}}, None];
  If[struc == "bcc" || struc == "fcc" || struc == "sc" || struc == "Honeycomb"||struc == "SquareLattice", Return[path], Print["Error : structure ", struc, " not implemented!"]]
]

GetBZHighSymK[struc_] := Module[{},
  Return[DeleteDuplicates[Transpose[GetBZPath[struc]]]]
]

AxialRotation[dir_, pos_] := Module[{R, axis, vec, \[Theta]},
  R = Fold[Dot, Table[RotationMatrix[120 Degree, {1, 1, 1}], dir]];
  axis = 0.1 {0, 0, 1};
  vec = 0.1 {1, 0, 0};
  Return[{Black, Arrowheads[0.03], Arrow[Tube[{pos + # - 0.3 (R.axis)\[Cross]#, pos + # + 0.3 (R.axis)\[Cross]#}, 0.01]]} & /@ Table[RotationMatrix[\[Theta], R.axis].R.vec, {\[Theta], 0, 2 Pi, (2 Pi)/10}]]
]  

PlotVectorField[latt_, field_, OptionsPattern[{"PointSize"->0.01, "PointColor" -> Black, "CellColor" -> Blue, "VectorColor" -> Green, "CoorColor" -> Gray, "TextColor" -> Black, "Text" -> True}]] := Module[{TextTable, VectorTable, PointTable, coor, cell},
  coor = {OptionValue["CoorColor"], Arrowheads[0.03], Arrow[Tube[{-#, #}, 0.005]]} & /@ latt;
  TextTable = Table[Text[Style[StringJoin[Riffle[ToString[#] & /@ Expand[latt.i], ","]], OptionValue["TextColor"]], latt.i + latt.{0.1, 0.1, 0.1}], {i, Flatten[field, 1]}];
  VectorTable = Table[{OptionValue["VectorColor"], Arrowheads[2*OptionValue["PointSize"]], Arrow@Tube[{latt.i[[1]], latt.(i[[1]] + i[[2]])}, 2.5*OptionValue["PointSize"]]}, {i, field}];
  PointTable = {OptionValue["PointColor"], PointSize[OptionValue["PointSize"]], Point[Table[latt.i, {i, Flatten[field, 1]}]]};
  cell = {OptionValue["CellColor"], Thickness[0.008], Line[If[Count[#2 - #1, 0] > 1, {latt.#1, latt.#2}, ## &[]] & @@@ Subsets[Permutations[{0, 0, 0, 1, 1, 1}, {3}], {2}]]};
  Return[If[OptionValue["Text"],
  Graphics3D[Join[VectorTable, PointTable, TextTable, coor, cell]],
  Graphics3D[Join[VectorTable, PointTable, coor, cell]]]]
]  

xright[x0_, nn_, NN_] := (x0 + nn) - Floor[N[x0 + nn - 1]/N[NN]] NN
xleft[x0_, nn_, NN_] := (x0 - nn) - Floor[N[x0 - nn - 1]/N[NN]] NN  

GetKpath[klist_, div_] := Module[{kpath, qpoints, label, distance, xticks, dk, i},
  {qpoints, label} = klist\[Transpose];
  {dk, kpath} = Join[{{0., First@qpoints}}, Flatten[Table[{Norm[N[1/div #2]], N[#1 + i/div #2]}, {i, div}] & @@@ ({qpoints[[1 ;; -   2]], Differences[qpoints]}\[Transpose]), 1]]\[Transpose];
  distance = Accumulate[dk];
  xticks = {Table[distance[[div (i - 1) + 1]], {i, Length[qpoints]}], label}\[Transpose];
  Return[{{distance, kpath}\[Transpose], xticks}]
]

LatticeFromLenAng[VecAng_] := Module[{aa, bb, cc, \[Alpha], \[Beta], \[Gamma], Lattice},
  {aa, bb, cc, \[Alpha], \[Beta], \[Gamma]} = VecAng;
  Lattice = {{aa, 0, 0},
             {bb Cos[\[Gamma]], bb Sin[\[Gamma]], 0},
             {cc Cos[\[Beta]], cc (Cos[\[Alpha]] - Cos[\[Beta]] Cos[\[Gamma]])/Sin[\[Gamma]], cc Sqrt[1 - (Cos[\[Alpha]]^2 + Cos[\[Beta]]^2 + Cos[\[Gamma]]^2) + 2 Cos[\[Alpha]] Cos[\[Beta]] Cos[\[Gamma]]]/Sin[\[Gamma]]}};
  Return[Chop[Lattice]]
]

SimplifyElementSymbol[ele_] := Module[{SimplifyList},
  SimplifyList = Which[StringQ[ele], First@StringCases[ele, RegularExpression["[[:upper:]][[:lower:]]*"]],
                       ListQ[ele], SimplifyElementSymbol[#] &/@ ele,
                       True, Print["Input type wrong!"];Abort[]];
  Return[SimplifyList]
]

PlotCrystal[latt_, pos_, dim_, OptionsPattern[{"vec" -> {}, "shift"->{0,0,0}, "AtomSize" -> 0.02, "ImageSize" -> 500}]] := Module[{xlimit, ylimit, zlimit, ExtVec, PointList, PosSymbols, PosTypes, PosColorCode, PosColorRange, PosColorList, PosList, vec},
  {xlimit, ylimit, zlimit} = If[Divisible[#, 2], {-#/2 + 1, #/2}, {-Floor[#/2], Floor[#/2]}] &/@ dim;
  ExtVec = If[OptionValue["vec"] != {}, Table[{vec[[1]], Arrowheads[1.0*OptionValue["AtomSize"]], Arrow[Tube[{latt.(#1-OptionValue["shift"]), latt.(#1+#2-OptionValue["shift"])}, 2.5*OptionValue["AtomSize"]]]} & @@@ vec[[2]], {vec, OptionValue["vec"]}], {{}}];
  PosSymbols = pos\[Transpose][[2]];
  PosTypes = GroupBy[PosSymbols, StringTake[#, First@First@StringPosition[#, "."] - 1] &];
  PosColorCode = Association[MapIndexed[#1 -> First[#2] &, Keys[PosTypes]]];
  PosColorRange = MinMax[Values[PosColorCode]];
  PosColorList = ColorData["Rainbow"][Rescale[PosColorCode[StringTake[#, First@First@StringPosition[#, "."] - 1]], {1, 2}]] & /@ PosSymbols;
  PointList = {PointSize[OptionValue["AtomSize"]], Point[Flatten[Table[latt.({ix, iy, iz} + #1 - OptionValue["shift"]), {ix, xlimit[[1]], xlimit[[2]]}, {iy, ylimit[[1]], ylimit[[2]]}, {iz, zlimit[[1]], zlimit[[2]]}], 2]]} & @@@ pos;
  PosList = {PosColorList, PointList}\[Transpose];
  Print[Graphics3D[Join[ExtVec,
                        {{Gray, PointSize[Medium],Point[Flatten[Table[latt.{ix, iy, iz}, {ix, xlimit[[1]], xlimit[[2]]}, {iy, ylimit[[1]], ylimit[[2]]}, {iz, zlimit[[1]], zlimit[[2]]}], 2]]}},
                        PosList,
                        {{Blue, Thickness[0.005], Line[If[Count[#2 - #1, 0] >= 1, {latt.Join[#1, {0}], latt.Join[#2, {0}]}, ## &[]] & @@@ Subsets[Permutations[{0, 0, 0, 1, 1, 1}, {2}], {2}]]}}, 
                        {{Blue, Arrowheads[0.03], Arrow[Tube[{latt.{0, 0, 0}, 1.15 latt.# + latt.{0, 0, 0}}, 0.005]]} & /@ (IdentityMatrix[3][[1 ;; 2]])},
                        {{Black, FontSize -> 10, Table[Text[StringJoin[Riffle[ToString[#] & /@ {ix, iy}, ","]], latt.{ix, iy, iz} + latt.{0.2, 0.2, 0.0}], {ix, xlimit[[1]], xlimit[[2]]}, {iy, ylimit[[1]], ylimit[[2]]}, {iz, zlimit[[1]], zlimit[[2]]}]}}],
                   ImageSize -> OptionValue["ImageSize"],
                   ViewPoint -> {0, 0, 100}]];
]

AtomicFormFactor[atom_, q_] := Module[{data, f, i, a, b, c},
  If[ListQ[atom], AtomicFormFactor[#] & /@ atom,
    data = Import[$UserBaseDirectory<>"/Applications/LINVARIANT/dataset/AtomicFormFactorsData.m"];
    c = Values[data[atom]][[-1]];
    {a, b} = Partition[Values[data[atom]][[2 ;; -2]], 2]\[Transpose];
    f = Sum[a[[i]] Exp[-b[[i]] (q/(4 Pi))^2], {i, 1, 4}] + c;
    Return[f];
    ];
]

XRayDiffraction[q_, latt_, Atoms_] := Module[{},
  Chop[1/Length[Atoms] Norm[Total[AtomicFormFactor[SimplifyElementSymbol[#2], Norm[2 Pi Inverse[latt].q]] Exp[-I 2 Pi q.#1] & @@@ Atoms]]^2]
]

GetBZKList[spg0_, kinv_] := Module[{kstar, klist, kx, ky, kz, i},
  kstar = DeleteDuplicates[Sort[DeleteDuplicates[Table[Mod[Values[spg0][[i]].#, 1][[1 ;; 3]], {i, Length[spg0]}]]] & /@ Flatten[Table[{kx, ky, kz, 0}, {kx, -0.5, 0.5, kinv}, {ky, -0.5, 0.5, kinv}, {kz, -0.5, 0.5, kinv}], 2], #1 == #2 &];
  klist = {PbcDiff[First[#]], Length[#]/(N@Length@Flatten[kstar, 1])} & /@ kstar;
  Return[klist]
]

PlotStruct[pos_] := Module[{latt, sites, SiteData, BondData, LabelData, tt, HoppingPairs},
   {latt, sites} = pos;
   SiteData = {ElementData[SimplifyElementSymbol[#2], "IconColor"], Specularity[White, 20], Sphere[latt.#1,QuantityMagnitude[ElementData[SimplifyElementSymbol[#2], "AtomicRadius"], "Angstroms"]/2]} & @@@ sites;
   LabelData = Text[Style[#2, Medium, Bold, Black], latt.#1] &@@@ sites;
   Graphics3D[Join[SiteData, LabelData],
                     ImageSize -> 500, Axes -> True,
                     AxesLabel -> (Style[#, Bold, 64] & /@ {"a", "b", "c"}),
                     ViewPoint -> {0, 0, \[Infinity]}]
]

(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
