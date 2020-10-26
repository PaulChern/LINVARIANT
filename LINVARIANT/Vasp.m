BeginPackage["LINVARIANT`Vasp`", {"LINVARIANT`INVARIANT`", "LINVARIANT`Structure`", "LINVARIANT`MathematicaPlus`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
ImportPOSCAR                 ::usage "ImportPOSCAR[f]"
ExportPOSCAR                 ::usage "ExportPOSCAR[dir, fname, f]"
ParseXML                     ::usage "ParseXML[xml, tag, label, level]"
ParseXMLData                 ::usage "ParseXMLData[xml, DataType]"
ParseVasprunBands            ::usage "ParseVasprunBands[xml]"
ParseVasprunFatBands         ::usage "ParseVasprunFatBands[xml]"
Kpoints2Kpath                ::usage "Kpoints2Kpath[kp]"
VaspKLabel                   ::usage "VaspKLabel[kpt]"
ImportWannierCHK             ::usage "ImportWannierCHK[chk]"
ParsePROCAR                  ::usage "ParsePROCAR[file]"
ImportKPATH                  ::usage "ImportKPATH[kpath]"
VaspBSPlot                   ::usage "VaspBSPlot[bsxml, klabels]"
BerryPol                     ::usage "BerryPol[vasprun]"
Rot2FFTGrid                  ::usage "Rot2FFTGrid[NGrid, k]"
GetBlochFunc                 ::usage "GetBlochFunc[file, is, ik, ib]"
GetTDM                       ::usage "GetTDM[fwavecar, poscar, skb1, skb2]"
WAVECAR2WANNIER              ::usage "WAVECAR2WANNIER[wavecar, vasprun]"
GetBlochLink                 ::usage "GetBlochLink[file]"
(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
ParseXML[xml_, tag_, label_, level_: Infinity] := Module[{xmldata, DataType, x},
  xmldata = Cases[xml, XMLElement[tag, Flatten[#1 -> #2 & @@@ If[label==={},label,If[Length@Dimensions[label]==1, {label}, label]]], x_] :> x, level];
  Return[xmldata]
]

ParseXMLData[xml_, DataType_] := Module[{xmldata},
  Flatten[Cases[#, XMLElement[DataType, {}, x_] :> ToExpression[StringSplit[x]]], 1] & /@ xml
]

ParseVasprunBands[xml_, OptionsPattern[{"fermi"->None}]] := Module[{is, ik, klist, NumKpoint, spinupdn=<||>, fermi, ISPIN},
  klist = Flatten[ParseXMLData[ParseXML[xml, "varray", {"name", "kpointlist"}], "v"], 1];
  fermi = If[OptionValue["fermi"]===None, ToExpression@First@First@ParseXML[xml, "i", {"name", "efermi"}], OptionValue["fermi"]];
  NumKpoint = Length[klist];
  ISPIN = ToExpression@First@First@ParseXML[xml, "i", {{"type", "int"}, {"name", "ISPIN"}}];

  spinupdn["k"] = klist;
  spinupdn["fermi"] = fermi;
  Do[spinupdn[If[is==1, "up", is==2, "dn"]] = Table[{#1-fermi,#2}&@@@Flatten[ParseXMLData[ParseXML[ParseXML[ParseXML[ParseXML[xml, "projected",  {}], "eigenvalues", {}], "set", {"comment", "spin "<>ToString[is]}], "set", {"comment", "kpoint " <> ToString[ik]}], "r"], 1], {ik, NumKpoint}], {is, ISPIN}];

  Return[spinupdn]
]

ParseVasprunFatBands[xml_, OptionsPattern[{"fermi"->None, "FatBandRange"->None}]] := Module[{is, ik, ib, klist, NumKpoint, spinupdn=<||>, character=<||>, fermi, NBANDS, ISPIN, fat, orbits, ifat0, ifat1},
  klist = Flatten[ParseXMLData[ParseXML[xml, "varray", {"name", "kpointlist"}], "v"], 1];
  fermi = If[OptionValue["fermi"]===None, ToExpression@First@First@ParseXML[xml, "i", {"name", "efermi"}], OptionValue["fermi"]];
  NumKpoint = Length[klist];
  NBANDS = ToExpression@First@First@ParseXML[xml, "i", {{"type", "int"}, {"name", "NBANDS"}}];
  ISPIN = ToExpression@First@First@ParseXML[xml, "i", {{"type", "int"}, {"name", "ISPIN"}}]; 
  {ifat0, ifat1} = If[OptionValue["FatBandRange"]===None, {Floor[0.75*NBANDS], NBANDS}, OptionValue["FatBandRange"]];

  spinupdn["k"] = klist;
  spinupdn["fermi"] = fermi;
  Do[spinupdn[If[is==1, "up", is==2, "dn"]] = Table[{#1-fermi,#2}&@@@(Flatten[ParseXMLData[ParseXML[ParseXML[ParseXML[ParseXML[xml, "projected", {}], "eigenvalues", {}], "set", {"comment", "spin "<>ToString[is]}], "set", {"comment", "kpoint " <> ToString[ik]}], "r"], 1][[ifat0;;ifat1]]), {ik, NumKpoint}], {is, ISPIN}];

  fat = Complement[ParseXML[ParseXML[xml, "projected", {}], "array", {}], ParseXML[ParseXML[ParseXML[xml, "projected", {}], "eigenvalues", {}], "array", {}]];
  orbits = Flatten[ParseXML[fat, "field", {}]];
  Print[orbits];

  character["orb"] = orbits;
  Do[character[If[is==1, "up", is==2, "dn"]] = Table[SparseArray[Flatten[ParseXMLData[ParseXML[ParseXML[ParseXML[fat, "set", {"comment", "spin"<>ToString[is]}], "set", {"comment", "kpoint "<>ToString[ik]}], "set", {"comment", "band "<>ToString[ib]}], "r"], 1]], {ik, NumKpoint}, {ib, ifat0, ifat1}], {is, ISPIN}];

  Return[{spinupdn, character}]
]

ImportPOSCAR[f_] := Module[{inp, Latt, xyz, xyzType, EleType, EleNum, TotalNum},
  inp = OpenRead[f];
  Read[inp, String];
  Latt = Quiet@ReadList[inp, Number, RecordLists -> True];
  Latt = (Latt[[1, 1]]*Latt[[2 ;; 4]])\[Transpose];
  EleType = ReadList[StringToStream@Read[inp, String], Word];
  EleNum = ReadList[StringToStream@Read[inp, String], Number];
  TotalNum = Total[EleNum];
  xyzType = ToLowerCase[Read[inp, String]];
  xyz = ReadList[inp, String, RecordLists -> True][[1 ;; TotalNum]];
  xyz = {ToExpression[First[StringSplit[#1]][[1 ;; 3]]], #2} & @@@ ({xyz, Flatten[ConstantArray @@@ ({EleType, EleNum}\[Transpose])]}\[Transpose]);
  Close[inp];
  Return[{Latt, xyz}]
]

ExportPOSCAR[dir_, fname_, f_] := Module[{i, pos, EleList, POSCAR},
  EleList = SimplifyElementSymbol[#] & /@ (f[[2]]\[Transpose][[2]]);
  pos = {f[[1]], {f[[2]]\[Transpose][[1]], EleList}\[Transpose]};
  POSCAR = Flatten /@ Join[{{"From Mathematica"}, {1}}, 
                           Table[ToString[NumberForm[N[i], {16, 15}]], {i, #}] & /@ pos[[1]], 
                           Through[{Keys, Values}[Counts[pos[[2]]\[Transpose][[-1]]]]], 
                           {{"Direct"}}, 
                           {Table[ToString[NumberForm[N[i], {17, 16}]], {i, #[[1]]}], #[[2]]} & /@ (pos[[2]])];
  Export[dir <> "/" <> fname, POSCAR, "Table", "FieldSeparators" -> "    "]
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

Kpoints2Kpath[kp_] := Module[{},
  {Accumulate[Join[{0}, Norm[#] & /@ Differences[kp]]], kp}\[Transpose]
]

VaspKLabel[kpt_] := Module[{},
  Which[StringContainsQ[kpt, "_"], 
        ToString[Subscript[VaspKLabel[#1], #2], StandardForm] & @@ StringSplit[kpt, "_"], 
        StringContainsQ[kpt, "\\"], ToString@ToExpression["\\[" <> "Capital" <> StringReplace[kpt, "\\" -> ""] <> "]"], 
        True, kpt]
]

ImportKPATH[kpath_] := Module[{KPATH, klist, kpt, i},
  KPATH = ReadList[kpath, String][[5 ;;]];
  klist = Table[kpt = StringSplit[KPATH[[i]]]; {ToExpression@kpt[[1 ;; 3]], VaspKLabel[kpt[[5]]]}, {i, Length@KPATH}];
  Prepend[Table[If[klist[[i + 1]] === klist[[i]], ## &[], klist[[i + 1]]], {i, Length[klist] - 1}], klist[[1]]]
]

ImportWannierCHK[chk_] := Module[{wchkfile, NumBands, NumExcludeBands, Lattice, NKPT, KGrid, KPOINTS, NumNNkpt, NumWann, chkpos, DisentanglementQ, omega, lwindow, ndimwin, UMatOpt, UMat, Mmat, WannCenters, WannSpreads, i, ik, iw, ib, jw, in},
  wchkfile = OpenRead[chk];
  SetStreamPosition[wchkfile, 0];
  ReadLine[wchkfile];
  NumBands = ToExpression[ReadLine[wchkfile]];
  NumExcludeBands = ToExpression[ReadLine[wchkfile]];
  Do[ReadLine[wchkfile], {i, NumExcludeBands}];
  Lattice = Partition[ToExpression[StringSplit[ReadLine[wchkfile]]], 3];
  ReadLine[wchkfile];
  NKPT = ToExpression[ReadLine[wchkfile]];
  KGrid = ToExpression[StringSplit[ReadLine[wchkfile]]];
  KPOINTS = Table[ToExpression[StringSplit[ReadLine[wchkfile]]], {i, NKPT}];
  NumNNkpt = ToExpression[ReadLine[wchkfile]];
  NumWann = ToExpression[ReadLine[wchkfile]];
  chkpos = ReadLine[wchkfile];
  DisentanglementQ = If[ToExpression[ReadLine[wchkfile]] == 0, False, True];
  If[DisentanglementQ,
     omega = ToExpression[ReadLine[wchkfile]];
     lwindow = Table[ToExpression[ReadLine[wchkfile]], {ik, NKPT}, {ib, NumBands}];
     ndimwin = Table[ToExpression[ReadLine[wchkfile]], {ik, NKPT}];
     UMatOpt = SparseArray[Table[ToExpression[StringSplit[ReadLine[wchkfile]]].{1, I}, {ik, NKPT}, {iw, NumWann}, {ib, NumBands}]]];
  UMat = Table[ToExpression[StringSplit[ReadLine[wchkfile]]].{1, I}, {ik, NKPT}, {iw, NumWann}, {jw, NumWann}];
  Mmat = Table[ToExpression[StringSplit[ReadLine[wchkfile]]].{1, I}, {ik, NKPT}, {in, NumNNkpt}, {iw, NumWann}, {jw, NumWann}];
  WannCenters = Mod[Inverse[Lattice\[Transpose]].#, 1] & /@ Table[ToExpression[StringSplit[ReadLine[wchkfile]]], {iw, NumWann}];
  WannSpreads = Table[ToExpression[ReadLine[wchkfile]], {iw, NumWann}];
  If[DisentanglementQ, Return[{UMatOpt, Umat, lwindow, ndimwin, Mmat, WannCenters}], Return[{UMat, lwindow, ndimwin, Mmat, WannCenters}]]
]

VaspBSPlot[bsxml_, klabels_, OptionsPattern[{"PlotRange" -> {All, All}, "FigRatio"->1/GoldenRatio}]] := Module[{kpath, updata, dndata, upplot, dnplot},
  kpath = {Accumulate[Join[{0}, Norm[#] & /@ Differences[bsxml["k"]]]], bsxml["k"]}\[Transpose];
  updata = {kpath\[Transpose][[1]], #\[Transpose][[1]]}\[Transpose] & /@ (bsxml["up"]\[Transpose]);
  upplot = ListLinePlot[updata, 
                        PlotStyle -> {{Black, Thick}}, 
                        Joined -> True, 
                        PlotRange -> OptionValue["PlotRange"], 
                        AspectRatio -> OptionValue["FigRatio"], 
                        Frame -> True, 
                        GridLines -> {{klabels\[Transpose][[1]], ConstantArray[Thick, Length[klabels]]}\[Transpose], Automatic},
                        FrameTicks -> {{Automatic, None}, {klabels, None}},
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
  NKPT = ToExpression[templine[[4]]];
  NBANDS = ToExpression[templine[[8]]];
  NumAtoms = ToExpression[templine[[12]]];
  bands = Table[templine = ReadNonEmptyLine[procar];
                weight = ToExpression[StringSplit[StringTrim[StringSplit[StringSplit[templine, ":"][[2]], "weight"][[2]]]][[2]]];
                kpt = ToExpression@StringPartition[" " <> StringTrim[StringSplit[StringSplit[templine, ":"][[2]], "weight"][[1]]], 11];
                Table[templine = StringSplit@ReadNonEmptyLine[procar];
                      ene = ToExpression[templine[[5]]];
                      occ = ToExpression[templine[[8]]];
                      templine = StringSplit@ReadNonEmptyLine[procar];
                      orbits = templine[[2 ;; -2]];
                      {kpt, ene, occ, weight, 
                       SparseArray[Table[ToExpression[StringSplit[ReadLine[procar]]][[2 ;; -2]], {i,NumAtoms + 1}][[1 ;; NumAtoms]]]}, {ib, NBANDS}], {ik, NKPT}];
  Return[{orbits, bands}]
]

BerryPol[vasprun_, coord_] := Module[{ZVAL, lattice, sites, pos, units, PolIon, PolEle, Ptot},
  ZVAL = Flatten[ConstantArray[{#2, #4}, #1] & @@@ Partition[ToExpression[#] & /@ Flatten@ParseXML[ParseXML[vasprun, "array", {"name", "atomtypes"}], "c", {}], 5], 1]\[Transpose];
  lattice = First@ParseXMLData[ParseXML[ParseXML[vasprun, "structure", {"name", "finalpos"}], "varray", {"name", "basis"}], "v"];
  sites = First@ParseXMLData[ParseXML[ParseXML[vasprun, "structure", {"name", "finalpos"}], "varray", {"name", "positions"}], "v"];
  pos = {lattice, Join[{sites}, ZVAL]\[Transpose]};
  Print[ZVAL];

  units = 1.6021766 10^-13 10^16/Det[lattice];
  PolIon =Chop@ToExpression@StringSplit[First@First@ParseXML[vasprun, "v", {"name", "PION"}]];
  PolEle =Chop@ToExpression@StringSplit[First@First@ParseXML[vasprun, "v", {"name", "PELC"}]];
  PolIon = Total[#1 #3 & @@@ (Join[{sites}, ZVAL]\[Transpose])];
  Ptot = units lattice\[Transpose].PbcDiff[Mod[PolIon - Inverse[(lattice\[Transpose])].PolEle, 1]];
  Ptot = If[coord == "Cartesian", 
            units lattice\[Transpose].PbcDiff[Mod[PolIon - Inverse[(lattice\[Transpose])].PolEle, 1]], 
            units (Norm[#] &@lattice)*PbcDiff[Mod[PolIon - Inverse[(lattice\[Transpose])].PolEle, 1]]];
  Return[Chop[Ptot, 10^-4]]
]

Rot2FFTGrid[NGrid_, k_] := Module[{},
  If[#2 > #1/2, #2 - #1, #2] & @@@ ({NGrid, k}\[Transpose])
]

GetBlochFunc[file_, is_, ik_, ib_, OptionsPattern[{"refine" -> 1}]] := Module[{Bohr2Ang, Ry2eV, EleKin, wavecar, LineRecord, NumSpin, ComplexTag, WFPrec, info, NumKPTs, NumBands, Ecut, lattice, NumGrid, Bands, PWCoeff, BinaryShift, temp, NumPW, icoeff, unkG, unkr, G, Rmesh1, Rmesh2, Rmesh, eikr, Ekin, kvec, En, fn, Phi, rR},
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
  NumGrid = OptionValue["refine"] (2 Ceiling[Norm[#]/(2 Pi Bohr2Ang) Sqrt[Ecut/Ry2eV]] + 1 & /@ lattice);

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
  
  kvec = Bands[[2]];
  {En, fn} = Bands[[3, ib]];
  Rmesh = Table[lattice\[Transpose].({i, j, k}/NumGrid), {i, 0, NumGrid[[1]] - 1}, {j, 0, NumGrid[[2]] - 1}, {k, 0, NumGrid[[3]] - 1}];
  eikr = Table[Exp[I 2 Pi kvec.({i, j, k}/NumGrid)], {i, 0, NumGrid[[1]] - 1}, {j, 0, NumGrid[[2]] - 1}, {k, 0, NumGrid[[3]] - 1}];
  icoeff = 0;
  unkG = TensorTranspose[Table[G = Rot2FFTGrid[NumGrid, {i, j, k}];
               Ekin = EleKin Norm[2 Pi Inverse[lattice].(kvec + G)]^2;
               If[Ekin < Ecut, icoeff = icoeff + 1; PWCoeff[[icoeff]], 0.0], 
               {k, 0, NumGrid[[3]] - 1}, {j, 0, NumGrid[[2]] - 1}, {i, 0, NumGrid[[1]] - 1}], {3,2,1}];
  unkr = InverseFourier[unkG, FourierParameters -> {0, -1}];
  Phi = eikr unkr;
  Return[{En, fn, Rmesh, Phi}]
]

GetTDM[fwavecar_, poscar_, skb_, LatticePadding_:{0,0,0}, OptionsPattern[{"refine" -> 1}]] := Module[{is, ik, ib, iskb, site, Bohr2Ang, Ry2eV, EleKin, LineRecord, wavecar, pos, NumSpin, ComplexTag, WFPrec, info, NumKPTs, NumBands, Ecut, lattice, NumGrid, Bands, PWCoeff, BinaryShift, temp, NumPW, icoeff, unkG, unkr, G, eikr, Ekin, DeltaK={0,0,0}, kvec, En, fn, Phi, rR, tdm, i, j, k, Ri, Rj, Rk},
  Bohr2Ang = 0.529177249;
  Ry2eV = 13.605826;
  EleKin = 3.8100198740807945;
  wavecar = OpenRead[fwavecar, BinaryFormat -> True];
  pos = ImportPOSCAR[poscar];

  (*Read info*)
  SetStreamPosition[wavecar, 0];
  {LineRecord, NumSpin, ComplexTag} = Rationalize[BinaryReadList[wavecar, "Real64", 3]];
  WFPrec = Which[ComplexTag == 45200, "Complex64", ComplexTag == 45210, "Complex128"];
  SetStreamPosition[wavecar, LineRecord];
  info = BinaryReadList[wavecar, "Real64", 12];
  {NumKPTs, NumBands} = Rationalize[info[[1 ;; 2]]];
  Ecut = info[[3]];
  lattice = Partition[info[[4 ;; 12]], 3];
  NumGrid = OptionValue["refine"] (2 Ceiling[Norm[#]/(2 Pi Bohr2Ang) Sqrt[Ecut/Ry2eV]] + 1 & /@ lattice);

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

    tdm = Sum[Sum[rR = lattice\[Transpose].({i, j, k}/NumGrid - #[[1]] + {Ri, Rj, Rk});
                  Exp[I 2 Pi DeltaK.{Ri,Rj,Rk}] Phi[[1,i+1,j+1,k+1]]\[Conjugate] Phi[[2,i+1,j+1,k+1]] rR Exp[-Norm[rR]](*/Norm[rR]^3*), 
                  {i, 0, NumGrid[[1]] - 1}, {j, 0, NumGrid[[2]] - 1}, {k, 0, NumGrid[[3]] - 1}], 
              {Ri, 0, LatticePadding[[1]]}, {Rj, 0, LatticePadding[[2]]}, {Rk, 0, LatticePadding[[3]]}] &/@ (pos[[2]]);
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

  OverlapMat = Table[Table[Table[ph1 = GetBlochFunc[dir0 <> "WAVECAR", is, ik-1, ib, "refine" -> 1][[4]];
                                 ph2 = GetBlochFunc[dir0 <> "WAVECAR", is, ik, jb, "refine" -> 1][[4]];
                                 Norm[Total[ph1\[Conjugate] ph2, 3]], {jb, ib, NumBands}], {ib, NumBands}], {is, NumSpin}, {ik, 2, NumKPTS}];
  PermChain = Table[Table[First@First@Position[OverlapMat[[is, ik, ib]], Max[OverlapMat[[is, ik, ib]]]]+ib-1, {ib, NumBands}], {is, NumSpin}, {ik, NumKPTS-1}];
  BlochLink = Table[Prepend[PermuteThrough[PermChain], Range[NumBands]], {is, NumSpin}];
  Return[BlochLink]
]
(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
