BeginPackage["LINVARIANT`Vasp`", {"LINVARIANT`INVARIANT`", "LINVARIANT`Structure`", "LINVARIANT`MathematicaPlus`", "LINVARIANT`Parser`", "LINVARIANT`DFT`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
ImportPOSCAR                 ::usage "ImportPOSCAR[f]"
ExportPOSCAR                 ::usage "ExportPOSCAR[dir, fname, f]"
ExportOPTCELL                ::usage "ExportOPTCELL[fname, character]"
ExportCHG                    ::usage "ExportCHG[dir, poscar, chg]"
ParseVasprunBands            ::usage "ParseVasprunBands[xml]"
ParseVasprunFatBands         ::usage "ParseVasprunFatBands[xml]"
VaspKLabel                   ::usage "VaspKLabel[kpt]"
ImportWannierCHK             ::usage "ImportWannierCHK[chk]"
ParsePROCAR                  ::usage "ParsePROCAR[file]"
ImportKPATH                  ::usage "ImportKPATH[kpath]"
VaspBSPlot                   ::usage "VaspBSPlot[bsxml, klabels]"
BerryPol                     ::usage "BerryPol[vasprun]"
Rot2FFTGrid                  ::usage "Rot2FFTGrid[NGrid, k]"
GetUnk                       ::usage "GetUnk[file, is, ik, ib]"
Unk2Bloch                    ::usage "Unk2Bloch[unkr]"
FunctionInherited            ::usage "FunctionInherited[latt, BaseSpace, fiber, func]"
GetTDM                       ::usage "GetTDM[fwavecar, poscar, skb1, skb2]"
WAVECAR2WANNIER              ::usage "WAVECAR2WANNIER[wavecar, vasprun]"
GetBlochLink                 ::usage "GetBlochLink[file]"
UoptPsik                     ::usage "UoptPsik[ik, iw, UmatOpt, ndimwin, lwindow, wavecar]"
UPsik                        ::usage "UPsik[ik, iw, UmatOpt, Umat, lwindow, ndimwin, wavecar]"
ReadVpot                     ::usage "ReadVpot[file]"
GetWannierFunc               ::usage "GetWannierFunc[iw, R, chk, wavecar_, vasprun]"
BundleUp                     ::usage "BundleUp[x, f, l]"
BundleDn                     ::usage "BundleUp[fx, l]"
GetElasticModuli             ::usage "GetElasticModuli[file, vol]"
GetDielectricTensor          ::usage "GetDielectricTensor[file]"
ReadForceConstants           ::usage "ReadForceConstants[file]"
ReadWannierCentres           ::usage "ReadWannierCentres[file]"
GetVaspTrajectory            ::usage "GetVaspTrajectory[dir]"

(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
ParseVasprunBands[xml_, OptionsPattern[{"fermi"->None}]] := Module[{latt, is, ik, klist, NumKpoint, spinupdn=<||>, fermi, ISPIN},
  latt = Last@ParseXMLData[ParseXML[xml, "varray", {"name", "basis"}], "v"];
  klist = Flatten[ParseXMLData[ParseXML[xml, "varray", {"name", "kpointlist"}], "v"], 1];
  fermi = If[OptionValue["fermi"]===None, ParseFortranNumber@First@First@ParseXML[xml, "i", {"name", "efermi"}], OptionValue["fermi"]];
  NumKpoint = Length[klist];
  ISPIN = ParseFortranNumber@First@First@ParseXML[xml, "i", {{"type", "int"}, {"name", "ISPIN"}}];

  spinupdn["latt"] = latt;
  spinupdn["k"] = klist;
  spinupdn["fermi"] = fermi;
  Do[spinupdn[If[is==1, "up", is==2, "dn"]] = Table[{#1-fermi,#2}&@@@Flatten[ParseXMLData[ParseXML[ParseXML[ParseXML[ParseXML[xml, "projected",  {}], "eigenvalues", {}], "set", {"comment", "spin "<>ToString[is]}], "set", {"comment", "kpoint " <> ToString[ik]}], "r"], 1], {ik, NumKpoint}], {is, ISPIN}];

  Return[spinupdn]
]

ParseVasprunFatBands[xml_, OptionsPattern[{"fermi"->None, "FatBandRange"->None}]] := Module[{latt, is, ik, ib, klist, NumKpoint, spinupdn=<||>, character=<||>, fermi, NBANDS, ISPIN, fat, orbits, ifat0, ifat1},
  latt = Last@ParseXMLData[ParseXML[xml, "varray", {"name", "basis"}], "v"];
  klist = Flatten[ParseXMLData[ParseXML[xml, "varray", {"name", "kpointlist"}], "v"], 1];
  fermi = If[OptionValue["fermi"]===None, ParseFortranNumber@First@First@ParseXML[xml, "i", {"name", "efermi"}], OptionValue["fermi"]];
  NumKpoint = Length[klist];
  NBANDS = ParseFortranNumber@First@First@ParseXML[xml, "i", {{"type", "int"}, {"name", "NBANDS"}}];
  ISPIN = ParseFortranNumber@First@First@ParseXML[xml, "i", {{"type", "int"}, {"name", "ISPIN"}}]; 
  {ifat0, ifat1} = If[OptionValue["FatBandRange"]===None, {Floor[0.75*NBANDS], NBANDS}, OptionValue["FatBandRange"]];

  spinupdn["k"] = klist;
  spinupdn["fermi"] = fermi;
  Do[spinupdn[Which[is==1, "up", is==2, "dn"]] = Table[{#1-fermi,#2}&@@@(Flatten[ParseXMLData[ParseXML[ParseXML[ParseXML[ParseXML[xml, "projected", {}], "eigenvalues", {}], "set", {"comment", "spin "<>ToString[is]}], "set", {"comment", "kpoint " <> ToString[ik]}], "r"], 1][[ifat0;;ifat1]]), {ik, NumKpoint}], {is, ISPIN}];

  fat = Complement[ParseXML[ParseXML[xml, "projected", {}], "array", {}], ParseXML[ParseXML[ParseXML[xml, "projected", {}], "eigenvalues", {}], "array", {}]];
  orbits = Flatten[ParseXML[fat, "field", {}]];
  Print[orbits];

  character["orb"] = orbits;
  Do[character[Which[is==1, "up", is==2, "dn"]] = Table[SparseArray[Flatten[ParseXMLData[ParseXML[ParseXML[ParseXML[fat, "set", {"comment", "spin"<>ToString[is]}], "set", {"comment", "kpoint "<>ToString[ik]}], "set", {"comment", "band "<>ToString[ib]}], "r"], 1]], {ik, NumKpoint}, {ib, ifat0, ifat1}], {is, ISPIN}];

  Return[{spinupdn, character}]
]

ImportPOSCAR[f_,OptionsPattern[{"originshift" -> {0,0,0}}]] := Module[{inp, Latt, xyz, xyzType, EleType, EleNum, TotalNum},
  inp = OpenRead[f];
  Read[inp, String];
  Latt = Quiet@ReadList[inp, Number, RecordLists -> True];
  Latt = Latt[[1, 1]]*Latt[[2 ;; 4]];
  EleType = ReadList[StringToStream@Read[inp, String], Word];
  EleNum = ReadList[StringToStream@Read[inp, String], Number];
  TotalNum = Total[EleNum];
  xyzType = ToLowerCase[Read[inp, String]];
  xyz = ReadList[inp, String, RecordLists -> True][[1 ;; TotalNum]];
  xyz = {Chop@Mod[ParseFortranNumber[First[StringSplit[#1]][[1 ;; 3]]] - OptionValue["originshift"], 1], #2} & @@@ ({xyz, Flatten[Table[#1 <> ToString[i], {i, #2}] & @@@ ({EleType, EleNum}\[Transpose])]}\[Transpose]);
  Close[inp];
  Return[{Latt, xyz}]
]

ExportPOSCAR[dir_, fname_, f_, OptionsPattern[{"head"->"From Mathematica"}]] := Module[{i, pos, EleList, POSCAR, head},
  head = OptionValue["head"];
  EleList = SimplifyElementSymbol[#] & /@ (f[[2]]\[Transpose][[2]]);
  pos = {f[[1]], {f[[2]]\[Transpose][[1]], EleList}\[Transpose]};
  POSCAR = Flatten /@ Join[{{head}, {1}}, 
                           Table[ToString[NumberForm[DecimalForm[N@i], {16, 15}]], {i, #}] & /@ pos[[1]], 
                           Through[{Keys, Values}[Counts[pos[[2]]\[Transpose][[-1]]]]], 
                           {{"Direct"}}, 
                           {Table[ToString[NumberForm[DecimalForm[N@i], {17, 16}]], {i, #[[1]]}], #[[2]]} & /@ (pos[[2]])];
  Export[dir <> "/" <> fname, POSCAR, "Table", "FieldSeparators" -> "    "]
]

ExportOPTCELL[fname_, character_, OptionsPattern[{"constrain"->True}]] := Module[{eta, n, i, j, e, sub, optcell},
  optcell = If[OptionValue["constrain"],
               eta = eta2eij[character];
               Flatten[Riffle[Table[sub = Thread[Table[Subscript[e, i, j], {i, 3}, {j, 3}] . (eta[[n, ;;]]) -> (eta[[n, ;;]])];
                                 Table[Subscript[e, i, j], {i, 3}, {j, 3}] /. sub /. Subscript[e, __] -> 0, {n, 3}], {{}}], 1],
               Flatten[If[Norm[character]==0, Table[0, {3}, {3}, {3}], Table[IdentityMatrix[3], {3}]], 1]];
  Export[fname, optcell, "Table", "FieldSeparators" -> " "];
]

ExportCHG[dir_, poscar_, chg_] := Module[{CGrid, rowdata, data, latt, pos, head},
  {latt, pos} = ImportPOSCAR[poscar];
  head = Join[{{"From Mathematica"}, {"1.0"}}, latt, {DeleteDuplicates[pos\[Transpose][[2]]]}, {Count[pos\[Transpose][[2]], #] & /@ DeleteDuplicates[pos\[Transpose][[2]]]}, {{"Direct"}}, pos\[Transpose][[1]], {{"  "}}];
  CGrid = Dimensions[chg];
  rowdata = Flatten[TensorTranspose[chg, {3, 2, 1}]];
  data = Join[head, {CGrid}, Partition[rowdata, 5], {rowdata[[-Mod[Length[rowdata], 5] ;; -1]]}];
  Export[dir, data, "Table", "FieldSeparators" -> "    "]
]

PlotGrid[latt_, pos_, dim_, OptionsPattern[{"vec" -> {}, "AtomSize" -> 0.02, "ImageSize" -> 500}]] := Module[{Nx, Ny, Nz, ExtVec, PointList, PosSymbols, PosTypes, PosColorCode, PosColorRange, PosColorList, PosList}, 
  {Nx, Ny, Nz} = dim;
  ExtVec = If[OptionValue["vec"] != {}, {OptionValue["vec"][[1]], Arrowheads[0.015], Arrow[Tube[{latt.#1, latt.#2}, 0.02]]} & @@@ OptionValue["vec"][[2]], {{}}];
  PosSymbols = pos\[Transpose][[2]];
  PosTypes = GroupBy[PosSymbols, StringTake[#, First@First@StringPosition[#, "."] - 1] &];
  PosColorCode = Association[MapIndexed[#1 -> First[#2] &, Keys[PosTypes]]];
  PosColorRange = MinMax[Values[PosColorCode]];
  PosColorList = ColorData["Rainbow"][Rescale[PosColorCode[StringTake[#, First@First@StringPosition[#, "."] - 1]], {1, 2}]] & /@ PosSymbols;
  PointList = {PointSize[OptionValue["AtomSize"]], Point[Flatten[Table[latt.({ix, iy, iz} + #1), {ix, -Nx, Nx}, {iy, -Ny, Ny}, {iz, -Nz, Nz}], 2]]} & @@@ pos;
  PosList = {PosColorList, PointList}\[Transpose];
  Print[Graphics3D[Join[ExtVec, 
                        {{Gray, PointSize[Medium],Point[Flatten[Table[latt.{ix, iy, iz}, {ix, -Nx, Nx}, {iy, -Ny, Ny}, {iz, -Nz, Nz}], 2]]}},
                        PosList,
                        {{Blue, Thickness[0.005], Line[If[Count[#2 - #1, 0] >= 1, {latt.Join[#1, {0}], latt.Join[#2, {0}]}, ## &[]] & @@@ Subsets[Permutations[{0, 0, 0, 1, 1, 1}, {2}], {2}]]}, {Blue, Arrowheads[0.03], Arrow[Tube[{latt.{0, 0, 0}, 1.15 latt.# + latt.{0, 0, 0}}, 0.005]]} & /@ (IdentityMatrix[3][[1 ;; 2]]), 
                         {Black, FontSize -> 10, Table[Text[StringJoin[Riffle[ToString[#] & /@ {ix, iy}, ","]], latt.{ix, iy, iz} + latt.{0.2, 0.2, 0.0}], {ix, -Nx, Nx}, {iy, -Ny, Ny}, {iz, -Nz, Nz}]}}],
                   ImageSize -> OptionValue["ImageSize"],
                   ViewPoint -> {0, 0, 100}]];
]

VaspKLabel[kpt_] := Module[{},
  Which[StringContainsQ[kpt, "_"], 
        ToString[Subscript[VaspKLabel[#1], #2], StandardForm] & @@ StringSplit[kpt, "_"], 
        StringContainsQ[kpt, "\\"], ToString@ToExpression["\\[" <> "Capital" <> StringReplace[kpt, "\\" -> ""] <> "]"], 
        True, kpt]
]

ImportKPATH[kpath_] := Module[{KPATH, klist, kpt, i},
  KPATH = ReadList[kpath, String][[5 ;;]];
  klist = Table[kpt = StringSplit[KPATH[[i]]]; {ParseFortranNumber[kpt[[1 ;; 3]]], VaspKLabel[kpt[[5]]]}, {i, Length@KPATH}];
  Prepend[Table[If[klist[[i + 1]] === klist[[i]], ## &[], klist[[i + 1]]], {i, Length[klist] - 1}], klist[[1]]]
]

ImportWannierCHK[chk_] := Module[{wchkfile, NumBands, NumExcludeBands, Lattice, NKPT, KGrid, KPOINTS, NumNNkpt, NumWann, chkpos, DisentanglementQ, omega, lwindow, ndimwin, UMatOpt, UMat, Mmat, WannCenters, WannSpreads, i, ik, iw, ib, jw, in},
  wchkfile = OpenRead[chk];
  SetStreamPosition[wchkfile, 0];
  ReadLine[wchkfile];
  NumBands = ParseFortranNumber[wchkfile];
  NumExcludeBands = ParseFortranNumber[wchkfile];
  Do[ReadLine[wchkfile], {i, NumExcludeBands}];
  Lattice = Partition[ParseFortranNumber[wchkfile], 3];
  ReadLine[wchkfile];
  NKPT = ParseFortranNumber[wchkfile];
  KGrid = ParseFortranNumber[wchkfile];
  KPOINTS = Table[ParseFortranNumber[wchkfile], {i, NKPT}];
  NumNNkpt = ParseFortranNumber[wchkfile];
  NumWann = ParseFortranNumber[wchkfile];
  chkpos = ReadLine[wchkfile];
  DisentanglementQ = If[ParseFortranNumber[wchkfile] == 0, False, True];
  If[DisentanglementQ,
     omega = ParseFortranNumber[wchkfile];
     lwindow = Table[ParseFortranNumber[wchkfile], {ik, NKPT}, {ib, NumBands}];
     ndimwin = Table[ParseFortranNumber[wchkfile], {ik, NKPT}];
     UMatOpt = SparseArray[Table[ParseFortranNumber[wchkfile].{1, I}, {ik, NKPT}, {iw, NumWann}, {ib, NumBands}]]];
  UMat = Table[ParseFortranNumber[wchkfile].{1, I}, {ik, NKPT}, {iw, NumWann}, {jw, NumWann}];
  Mmat = Table[ParseFortranNumber[wchkfile].{1, I}, {ik, NKPT}, {in, NumNNkpt}, {iw, NumWann}, {jw, NumWann}];
  WannCenters = Mod[Inverse[Lattice\[Transpose]].#, 1] & /@ Table[ParseFortranNumber[wchkfile], {iw, NumWann}];
  WannSpreads = Table[ParseFortranNumber[wchkfile], {iw, NumWann}];
  Close[wchkfile];
  If[DisentanglementQ, Return[{UMatOpt, UMat, lwindow, ndimwin, Mmat, WannCenters}], Return[{UMat, lwindow, ndimwin, Mmat, WannCenters}]]
]

VaspBSPlot[bsxml_, hkpts_, OptionsPattern[{"PlotRange" -> {All, All}, "AspectRatio"->1/GoldenRatio}]] := Module[{latt, blatt, kpath, updata, dndata, upplot, dnplot, klabels},
  latt = bsxml["latt"];
  blatt = 1/(2*Pi)*Inverse[latt]\[Transpose];
  klabels = Kpoints2Kpath[hkpts,blatt];
  kpath = {Accumulate[Join[{0}, Norm[blatt.#] & /@ Differences[bsxml["k"]]]], bsxml["k"]}\[Transpose];
  updata = {kpath\[Transpose][[1]], #\[Transpose][[1]]}\[Transpose] & /@ (bsxml["up"]\[Transpose]);
  upplot = ListLinePlot[updata, 
                        PlotStyle -> {{Black, Thick}}, 
                        Joined -> True, 
                        PlotRange -> OptionValue["PlotRange"], 
                        AspectRatio -> OptionValue["AspectRatio"], 
                        Frame -> True, 
                        GridLines -> {{klabels\[Transpose][[2]], ConstantArray[Thick, Length[klabels]]}\[Transpose], Automatic},
                        FrameTicks -> {{Automatic, None}, {klabels[[;;,2;;3]], None}},
                        ImageSize -> Medium];
  If[MemberQ[Keys[bsxml], "dn"], 
     dndata = {kpath\[Transpose][[1]], #\[Transpose][[1]]}\[Transpose] & /@ (bsxml["dn"]\[Transpose]);
     dnplot = ListLinePlot[dndata,  
                           PlotStyle -> {{Black, Dashed, Thin}},  
                           Joined -> True,  
                           PlotRange -> OptionValue["PlotRange"],  
                           AspectRatio -> OptionValue["FigRatio"],  
                           Frame -> True, 
                           GridLines -> {{klabels\[Transpose][[1]], ConstantArray[Thick, Length[klabels]]}\[Transpose], Automatic},
                           FrameTicks -> {{Automatic, None}, {klabels, None}},
                           ImageSize -> Medium];
     Return[Show[{upplot,dnplot}]], 
     Return[upplot]]
]

ParsePROCAR[file_] := Module[{procar, templine, NKPT, NBANDS, NumAtoms, bands, weight, kpt, ene, occ, orbits},
  procar = OpenRead[file];
  SetStreamPosition[procar, 0];
  ReadNonEmptyLine[procar];
  templine = StringSplit@ReadNonEmptyLine[procar];
  NKPT = ParseFortranNumber[templine[[4]]];
  NBANDS = ParseFortranNumber[templine[[8]]];
  NumAtoms = ParseFortranNumber[templine[[12]]];
  bands = Table[templine = ReadNonEmptyLine[procar];
                weight = ParseFortranNumber[StringSplit[StringTrim[StringSplit[StringSplit[templine, ":"][[2]], "weight"][[2]]]][[2]]];
                kpt = ParseFortranNumber[StringPartition[" " <> StringTrim[StringSplit[StringSplit[templine, ":"][[2]], "weight"][[1]]], 11]];
                Table[templine = StringSplit@ReadNonEmptyLine[procar];
                      ene = ParseFortranNumber[templine[[5]]];
                      occ = ParseFortranNumber[templine[[8]]];
                      templine = StringSplit@ReadNonEmptyLine[procar];
                      orbits = templine[[2 ;; -2]];
                      {kpt, ene, occ, weight, 
                       SparseArray[Table[ParseFortranNumber[procar][[2 ;; -2]], {i,NumAtoms + 1}][[1 ;; NumAtoms]]]}, {ib, NBANDS}], {ik, NKPT}];
  Close[procar];
  Return[{orbits, bands}]
]

BerryPol[vasprun_, coord_] := Module[{ZVAL, lattice, sites, pos, units, PolIon, PolEle, Ptot},
  ZVAL = Flatten[ConstantArray[{#2, ParseFortranNumber@#4}, ParseFortranNumber@#1] & @@@ Partition[Flatten@ParseXML[ParseXML[vasprun, "array", {"name", "atomtypes"}], "c", {}], 5], 1]\[Transpose];
  lattice = First@ParseXMLData[ParseXML[ParseXML[vasprun, "structure", {"name", "finalpos"}], "varray", {"name", "basis"}], "v"];
  sites = First@ParseXMLData[ParseXML[ParseXML[vasprun, "structure", {"name", "finalpos"}], "varray", {"name", "positions"}], "v"];
  pos = {lattice, Join[{sites}, ZVAL]\[Transpose]};

  units = 1.6021766 10^3/Det[lattice];
  PolIon =Chop@ParseFortranNumber@StringSplit[First@First@ParseXML[vasprun, "v", {"name", "PION"}]];
  PolEle =Chop@ParseFortranNumber@StringSplit[First@First@ParseXML[vasprun, "v", {"name", "PELC"}]];
  PolIon = Total[#1 #3 & @@@ (Join[{sites}, ZVAL]\[Transpose])];
  (*Print[PolIon, Inverse[(lattice\[Transpose])].PolEle];*)
  Ptot = units lattice\[Transpose].PbcDiff[Mod[PolIon - Inverse[(lattice\[Transpose])].PolEle, 1]];
  Ptot = If[coord == "Cartesian", 
            units lattice\[Transpose].PbcDiff[Mod[PolIon - Inverse[(lattice\[Transpose])].PolEle, 1]], 
            units (Norm[#] &@lattice)*PbcDiff[Mod[PolIon - Inverse[(lattice\[Transpose])].PolEle, 1]]];
  PolEle = units lattice\[Transpose].PbcDiff[Mod[-Inverse[(lattice\[Transpose])].PolEle, 1]];
  PolIon = units lattice\[Transpose].PbcDiff[Mod[PolIon, 1]];
  Return[Chop[{Ptot, PolIon, PolEle}, 10^-4]]
]

ReadVpot[file_] := Module[{data, pot, Nx, Ny, Nz, i, j, k},
  data = ReadList[file, Real];
  {Nx, Ny, Nz} = Rationalize[data[[1 ;; 3]]];
  pot = TensorTranspose[ArrayReshape[data[[4 ;;]], {Nz, Ny, Nx}], {3, 2, 1}];
  Return[pot]
]

Rot2FFTGrid[NGrid_, k_] := Module[{},
  If[#2 >= #1/2, #2 - #1, #2] & @@@ ({NGrid, k}\[Transpose])
]

GetUnk[file_, is_, ik_, ib_, OptionsPattern[{"refine" -> 1, "sc"->{1,1,1}, "shift"->{0,0,0}}]] := Module[{Bohr2Ang, Ry2eV, EleKin, wavecar, LineRecord, NumSpin, ComplexTag, WFPrec, info, NumKPTs, NumBands, Ecut, lattice, NumGrid, Bands, PWCoeff, BinaryShift, temp, NumPW, icoeff, unkG, unkr0, unkr, G, BaseSpace, eikr, Ekin, kvec, En, fn, Phi, rR, Ri, Rj, Rk, ii, jj, kk, shift},
  {Ri, Rj, Rk} = OptionValue["sc"];
  shift = OptionValue["shift"];
  Bohr2Ang = 0.529177249;
  Ry2eV = 13.605826;
  EleKin = 3.8100198740807945;
  wavecar = OpenRead[file, BinaryFormat -> True];
  
  (* Read info *)
  SetStreamPosition[wavecar, 0];
  {LineRecord, NumSpin, ComplexTag} = Rationalize[BinaryReadList[wavecar, "Real64", 3]];
  WFPrec = Which[ComplexTag == 45200, "Complex64", ComplexTag == 45210, "Complex128"];
  SetStreamPosition[wavecar, LineRecord];
  info = BinaryReadList[wavecar, "Real64", 12];
  {NumKPTs, NumBands} = Rationalize[info[[1 ;; 2]]];
  Ecut = info[[3]];
  lattice = Partition[info[[4 ;; 12]], 3];
  NumGrid = OptionValue["refine"] 2 Ceiling[Norm[#]/(2 Pi Bohr2Ang) Sqrt[Ecut/Ry2eV]] + 1 & /@ lattice;

  (* Read Bands *)
  BinaryShift = 2 + (is - 1) NumKPTs (NumBands + 1) + (ik - 1) (NumBands + 1);
  SetStreamPosition[wavecar, BinaryShift LineRecord];
  temp = Chop@BinaryReadList[wavecar, "Real64", 4 + 3 NumBands];
  Bands = {Rationalize[temp[[1]]], temp[[2 ;; 4]], {#1, #3} & @@@ Partition[temp[[5 ;;]], 3]};

  (* Read PW Coefficients *)
  NumPW = Bands[[1]];
  BinaryShift = 2 + (is - 1) NumKPTs (NumBands + 1) + (ik - 1) (NumBands + 1) + ib;
  SetStreamPosition[wavecar, BinaryShift LineRecord];
  PWCoeff = Chop@Normalize@BinaryReadList[wavecar, WFPrec, NumPW];

  Close[wavecar];
  
  kvec = Bands[[2]];
  {En, fn} = Bands[[3, ib]];
  (*Rmesh = Table[Rot2FFTGrid[Norm[#] & /@ lattice, lattice\[Transpose].({i, j, k}/NumGrid)], {i, 0, NumGrid[[1]] - 1}, {j, 0, NumGrid[[2]] - 1}, {k, 0, NumGrid[[3]] - 1}];*)
  (*Rmesh = Table[{i, j, k}/NumGrid, {i, 0, NumGrid[[1]] - 1}, {j, 0, NumGrid[[2]] - 1}, {k, 0, NumGrid[[3]] - 1}];*)
  eikr = Table[Exp[I 2 Pi kvec.({i, j, k}/NumGrid)], {i, 0, NumGrid[[1]] - 1}, {j, 0, NumGrid[[2]] - 1}, {k, 0, NumGrid[[3]] - 1}];
  icoeff = 0;
  unkG = TensorTranspose[Table[G = Rot2FFTGrid[NumGrid, {i, j, k}];
               Ekin = EleKin Norm[2 Pi Inverse[lattice].(kvec + G)]^2;
               If[Ekin < Ecut, icoeff = icoeff + 1; PWCoeff[[icoeff]], 0.0], 
               {k, 0, NumGrid[[3]] - 1}, {j, 0, NumGrid[[2]] - 1}, {i, 0, NumGrid[[1]] - 1}], {3,2,1}];
  unkr0 = InverseFourier[unkG, FourierParameters -> {0, -1}];
  (*Phi = eikr unkr;*)
  BaseSpace = Table[N[Mod[{i, j, k}/NumGrid+shift,1]], {i, 0, NumGrid[[1]] Ri - 1}, {j, 0, NumGrid[[2]] Rj - 1}, {k, 0, NumGrid[[3]] Rk - 1}];
  unkr = Table[{ii, jj, kk} = MapThread[Mod[#1, #2, 1] &, {{i+1, j+1, k+1}, NumGrid}];
               unkr0[[ii,jj,kk]], {i, 0, NumGrid[[1]] Ri - 1}, {j, 0, NumGrid[[2]] Rj - 1}, {k, 0, NumGrid[[3]] Rk - 1}];
  Return[{En, fn, lattice, kvec, BaseSpace, unkr}]
]

BundleUp[x_, f_, l_] := Module[{},
  MapThread[Join[#1, {#2}] &, {x, f}, l]
]

BundleDn[fx_, l_] := Module[{},
  {Map[#[[1;;3]]&, fx, {l}], Map[#[[4]]&, fx, {l}]}
]

Unk2Bloch[unkr_, OptionsPattern[{"coord" -> "Cartesian", "plot" -> False}]] := Module[{coord, plot, bloch, unk, latt, kvec, basespace, En, fn},
  coord = OptionValue["coord"];
  plot = OptionValue["plot"];
  {En, fn, latt, kvec, basespace, unk} = unkr;
  bloch = Which[
    plot && coord === "Direct",
    BundleUp[basespace, MapThread[Sign[Arg[Exp[I 2 Pi kvec.#1] #2]] Abs[Exp[I 2 Pi kvec.#1] #2] &, {basespace, unk}, 3], 3],
    plot && coord === "Cartesian",
    BundleUp[Map[latt\[Transpose].# &,basespace,{3}], MapThread[Sign[Arg[Exp[I 2 Pi kvec.#1] #2]] Abs[Exp[I 2 Pi kvec.#1] #2] &, {basespace, unk}, 3], 3],
    Not@plot && coord === "Direct",
    {basespace, MapThread[Exp[I 2 Pi kvec.#1] #2 &, {basespace, unk}, 3]},
    Not@plot && coord === "Cartesian",
    {Map[latt\[Transpose].# &,basespace,{3}], MapThread[Exp[I 2 Pi kvec.#1] #2 &, {basespace, unk}, 3]}
  ];
  Return[bloch]
]

FunctionInherited[latt_, BaseSpace_, fiber_, func_] := Module[{f, x, ndim},
  ndim = Length[Dimensions[BaseSpace]];
  x = Map[latt\[Transpose] . # &, BaseSpace, {ndim}];
  f = MapThread[func, {x, fiber}, ndim];
  Return[f]
]

GetTDM[fwavecar_, poscar_, skb_, LatticePadding_:{0,0,0}, OptionsPattern[{"refine" -> 1}]] := Module[{is, ik, ib, iskb, site, Bohr2Ang, Ry2eV, EleKin, LineRecord, wavecar, pos, NumSpin, ComplexTag, WFPrec, info, NumKPTs, NumBands, Ecut, lattice, NumGrid, Bands, PWCoeff, BinaryShift, temp, NumPW, icoeff, unkG, unkr, G, eikr, Ekin, DeltaK={0,0,0}, kvec, En, fn, Phi, rR, tdm, i, j, k, Ri, Rj, Rk, p1, p2, p3, singularR},
  Bohr2Ang = 0.529177249;
  Ry2eV = 13.605826;
  EleKin = 3.8100198740807945;
  wavecar = OpenRead[fwavecar, BinaryFormat -> True];
  pos = ImportPOSCAR[poscar];
  {p1, p2, p3} = LatticePadding;

  (*Read info*)
  SetStreamPosition[wavecar, 0];
  {LineRecord, NumSpin, ComplexTag} = Rationalize[BinaryReadList[wavecar, "Real64", 3]];
  WFPrec = Which[ComplexTag == 45200, "Complex64", ComplexTag == 45210, "Complex128"];
  SetStreamPosition[wavecar, LineRecord];
  info = BinaryReadList[wavecar, "Real64", 12];
  {NumKPTs, NumBands} = Rationalize[info[[1 ;; 2]]];
  Ecut = info[[3]];
  lattice = Partition[info[[4 ;; 12]], 3];
  NumGrid = OptionValue["refine"] 2 Ceiling[Norm[#]/(2 Pi Bohr2Ang) Sqrt[Ecut/Ry2eV]] + 1 & /@ lattice;
  singularR=Min[lattice\[Transpose].({1, 1, 1}/NumGrid)]/1000;

  Phi = Table[
    {is, ik, ib} = skb[[iskb]];
    (*Read Bands*)
    BinaryShift = 2 + (is - 1) NumKPTs (NumBands + 1) + (ik - 1) (NumBands + 1);
    SetStreamPosition[wavecar, BinaryShift LineRecord];
    temp = Chop@BinaryReadList[wavecar, "Real64", 4 + 3 NumBands];
    Bands = {Rationalize[temp[[1]]], temp[[2 ;; 4]], {#1, #3} & @@@ Partition[temp[[5 ;;]], 3]};
    (*Read PW Coefficients*)
    NumPW = Bands[[1]];
    BinaryShift = 2 + (is - 1) NumKPTs (NumBands + 1) + (ik - 1) (NumBands + 1) + ib;
    SetStreamPosition[wavecar, BinaryShift LineRecord];
    PWCoeff = Chop@Normalize@BinaryReadList[wavecar, WFPrec, NumPW];
    kvec = Bands[[2]];
    DeltaK = DeltaK + If[iskb==1, -kvec, kvec];
    {En, fn} = Bands[[3, ib]];
    eikr = Table[Exp[I 2 Pi kvec.({i, j, k}/NumGrid)], {i, 0, NumGrid[[1]] - 1}, {j, 0, NumGrid[[2]] - 1}, {k, 0, NumGrid[[3]] - 1}];
    icoeff = 0;
    unkG = TensorTranspose[Table[G = Rot2FFTGrid[NumGrid, {i, j, k}];
                 Ekin = EleKin Norm[2 Pi Inverse[lattice].(kvec + G)]^2;
                 If[Ekin < Ecut, icoeff = icoeff + 1; PWCoeff[[icoeff]], 0.0], {k, 0, NumGrid[[3]] - 1}, {j, 0, NumGrid[[2]] - 1}, {i, 0, NumGrid[[1]] - 1}], {3, 2, 1}];
    unkr = InverseFourier[unkG, FourierParameters -> {0, -1}];
    eikr unkr, {iskb, 2}];

    tdm = Sum[Sum[rR = Rot2FFTGrid[Norm[#]&/@lattice, lattice\[Transpose].({i, j, k}/NumGrid)] 
                     + lattice\[Transpose].{Ri, Rj, Rk} 
                     - Rot2FFTGrid[Norm[#]&/@lattice, #[[1]]];
                  Exp[I 2 Pi DeltaK.{Ri,Rj,Rk}] Phi[[1,i+1,j+1,k+1]]\[Conjugate] Phi[[2,i+1,j+1,k+1]] rR If[Norm[rR]<singularR,0,1/Norm[rR]^3], 
                  {i, 0, NumGrid[[1]]-1}, {j, 0, NumGrid[[2]]-1}, {k, 0, NumGrid[[3]]-1}], 
              {Ri, -p1, p1}, {Rj, -p2, p2}, {Rk, -p3, p3}] &/@ (pos[[2]]);
  Return[tdm]
]

WAVECAR2WANNIER[wavecar_, vasprun_] := Module[{klist, lattice, sites, pos, bin, xml, weight, LineRecord, WFPrec, NumSpin, ComplexTag, BlochInfo, NumKPTs, NumBands},
  bin = OpenRead[wavecar, BinaryFormat -> True];
  xml = Import[dir0 <> "vasprun.xml"];

  (* Read vasprun info *)
  lattice = First@ParseXMLData[ParseXML[ParseXML[vasprun, "structure", {"name", "finalpos"}], "varray", {"name", "basis"}], "v"];
  sites = First@ParseXMLData[ParseXML[ParseXML[vasprun, "structure", {"name", "finalpos"}], "varray", {"name", "positions"}], "v"];
  pos = {lattice, sites};
  klist = First@ParseXMLData[ParseXML[xml, "varray", {"name", "kpointlist"}], "v"];
  weight = Flatten@ParseXMLData[ParseXML[xml, "varray", {"name", "weights"}], "v"];

  (* Read wavecar info *)
  SetStreamPosition[bin, 0];
  {LineRecord, NumSpin, ComplexTag} = Rationalize[BinaryReadList[bin, "Real64", 3]];
  WFPrec = Which[ComplexTag == 45200, "Complex64", ComplexTag == 45210, "Complex128"];
  SetStreamPosition[bin, LineRecord];
  BlochInfo = BinaryReadList[bin, "Real64", 12];
  {NumKPTs, NumBands} = Rationalize[BlochInfo[[1 ;; 2]]];
  Close[bin];

  Return[klist]
]

GetBlochLink[file_] := Module[{BlochLink, wavecar, LineRecord, NumSpin, ComplexTag, WFPrec, info, NumKPTs, NumBands, ib, jb, is, ik, jk},
  wavecar = OpenRead[file, BinaryFormat -> True];
  SetStreamPosition[wavecar, 0];
  {LineRecord, NumSpin, ComplexTag} = Rationalize[BinaryReadList[wavecar, "Real64", 3]];
  WFPrec = Which[ComplexTag == 45200, "Complex64", ComplexTag == 45210, "Complex128"];
  SetStreamPosition[wavecar, LineRecord];
  info = BinaryReadList[wavecar, "Real64", 12];
  {NumKPTs, NumBands} = Rationalize[info[[1 ;; 2]]];
  Close[wavecar];

  OverlapMat = Table[Table[Table[ph1 = GetUnk[dir0 <> "WAVECAR", is, ik-1, ib, "refine" -> 1][[4]];
                                 ph2 = GetUnk[dir0 <> "WAVECAR", is, ik, jb, "refine" -> 1][[4]];
                                 Norm[Total[ph1\[Conjugate] ph2, 3]], {jb, ib, NumBands}], {ib, NumBands}], {is, NumSpin}, {ik, 2, NumKPTS}];
  PermChain = Table[Table[First@First@Position[OverlapMat[[is, ik, ib]], Max[OverlapMat[[is, ik, ib]]]]+ib-1, {ib, NumBands}], {is, NumSpin}, {ik, NumKPTS-1}];
  BlochLink = Table[Prepend[PermuteThrough[PermChain], Range[NumBands]], {is, NumSpin}];
  Return[BlochLink]
]

UoptPsik[ik_, iw_, UmatOpt_, ndimwin_, lwindow_, wavecar_, resolution_] := Module[{i, shift, Psik},
  shift = First@First@Position[lwindow[[ik]], 1] - 1;
  Psik = Total@Table[If[Norm[UmatOpt[[ik, iw, i]]] == 0, ## &[], UmatOpt[[ik, iw, i]] GetUnk[wavecar, 1, ik, i + shift, "refine" -> resolution][[3]]], {i, ndimwin[[ik]]}];
  Return[Psik]
]

UPsik[ik_, iw_, UmatOpt_, Umat_, lwindow_, ndimwin_, wavecar_, resolution_] := Module[{jw, PsikOpt, Psik},
  Psik = Sum[PsikOpt = UoptPsik[ik, jw, UmatOpt, ndimwin, lwindow, wavecar, resolution];
             Umat[[ik]][[iw, jw]] PsikOpt, {jw, 72}];
  Return[Psik]
];

GetWannierFunc[iw_, R_, chk_, wavecar_, vasprun_, OptionsPattern[{"refine" -> 1}]] := Module[{xml, UmatOpt, Umat, lwindow, ndimwin, Mmat, WannCenters, Wann, NumKpoints, klist},
  {UmatOpt, Umat, lwindow, ndimwin, Mmat, WannCenters} = ImportWannierCHK[chk];
  xml = Import[vasprun];
  klist = Flatten[ParseXMLData[ParseXML[xml, "varray", {"name", "kpointlist"}], "v"], 1];
  NumKpoints = Length[ndimwin];
  Wann = Sum[Exp[-I 2 Pi klist[[ik]].R] UPsik[ik,iw,UmatOpt,Umat,lwindow,ndimwin,wavecar,OptionValue["refine"]], {ik, NumKpoints}]/NumKpoints;
  Return[Wann]
]

GetElasticModuli[file_, vol_] := Module[{KBar2eV, mat, out, outcar, tmp, Cijkl},
  KBar2eV = 0.1 vol/160.21766208;
  mat = {};
  outcar = OpenRead[file];
  Find[outcar, "TOTAL ELASTIC MODULI (kBar)"];
  tmp = ReadLine[outcar];
  tmp = ReadLine[outcar];
  Do[AppendTo[mat, ParseFortranNumber[outcar][[2 ;; 7]]], {6}];
  Cijkl = SparseArray[{{i_, i_, i_, i_} -> ToExpression["C1111"], 
                       {i_, i_, j_, j_} /; i != j -> ToExpression["C1122"], 
                       {i_, j_, i_, j_} /; i != j -> ToExpression["C1212"]}, {3, 3, 3, 3}];
  out = DeleteDuplicates@Flatten[Table[Which[i == j && k == l,
                                             Cijkl[[i, j, k, l]] -> mat[[i, k]] KBar2eV,
                                             i == k && j == l && i != j,
                                             Cijkl[[i, j, k, l]] -> mat[[i + 3, i + 3]] KBar2eV,
                                             True, ## &[]], {i, 3}, {j, 3}, {k, 3}, {l, 3}]];
  Close[outcar];
  Return[out]
]

GetDielectricTensor[file_] := Module[{KBar2eV, mat, out, outcar, tmp},
  mat = {};
  outcar = OpenRead[file];
  Find[outcar, "MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in RPA (Hartree))"];
  tmp = ReadLine[outcar];
  Do[AppendTo[mat, ParseFortranNumber[outcar]], {3}];
  Close[outcar]
  Return[mat]
]

ReadForceConstants[file_] := Module[{data, NumAtom0, NumAtom, FC},
  data = ReadList[file, Number, RecordLists -> True];
  {NumAtom0, NumAtom} = data[[1]];
  FC = Association[#1 -> {#2, #3, #4} & @@@ Partition[data[[2 ;;]], 4]];
  FC = Table[FC[{i, j}], {i, NumAtom0}, {j, NumAtom}];
  Return[FC]
]

ReadWannierCentres[file_] := Module[{inp},
  inp = OpenRead[file];
  Read[inp, Number];
  Read[inp, String];
  If[#1 === "X", {ToExpression[#2], ToExpression[#3], ToExpression[#4]}, ## &[]] & @@@ StringSplit[ReadList[inp, String, RecordLists -> False]]
]

GetVaspTrajectory[dir_, SYSTEM_] := Module[{XDATCAR, NSteps, latt0, sites0, trj, sites, latt},
  XDATCAR = OpenRead[dir <> "/XDATCAR"];
  NSteps = Length[FindList[XDATCAR, "Direct configuration="]];
  {latt0, sites0} = ImportPOSCAR[dir <> "/POSCAR"];
  XDATCAR = OpenRead[dir0 <> "/XDATCAR"];
  trj = Table[If[! (Find[XDATCAR, SYSTEM] === EndOfFile),
                 ReadLine[XDATCAR]; 
                 latt = Table[ParseFortranNumber[ReadLine[XDATCAR]], {3}];
                 ReadLine[XDATCAR];
                 ReadLine[XDATCAR];
                 ReadLine[XDATCAR];
                 sites = Table[ParseFortranNumber[ReadLine[XDATCAR]], {Length@sites0}];
                 {latt, {sites, sites0\[Transpose][[2]]}\[Transpose]}], {NSteps}];
  Close[XDATCAR];
  Return[trj]
]

(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
