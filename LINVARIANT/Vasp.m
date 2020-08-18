BeginPackage["LINVARIANT`Vasp`", {"LINVARIANT`INVARIANT`", "LINVARIANT`Structure`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
ImportPOSCAR                 ::usage "ImportPOSCAR[f]"
ExportPOSCAR                 ::usage "ExportPOSCAR[dir, fname, f]"
ParseXML                     ::usage "ParseXML[xml, tag, label, level]"
ParseXMLData                 ::usage "ParseXMLData[xml, DataType]"
ParseVasprunBands            ::usage "ParseVasprunBands[xml]"
Kpoints2Kpath                ::usage "Kpoints2Kpath[kp]"
VaspBSPlot                   ::usage "VaspBSPlot[bsxml, klabels]"
(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
ParseXML[xml_, tag_, label_, level_: Infinity] := Module[{xmldata, DataType},
  xmldata = Cases[xml, XMLElement[tag, #1 -> #2 & @@@ label, x_] :> x, level];
  Return[xmldata]
]

ParseXMLData[xml_, DataType_] := Module[{xmldata},
  Flatten[Cases[#, XMLElement[DataType, {}, x_] :> ToExpression[StringSplit[x]]], 1] & /@ xml
]

ParseVasprunBands[xml_] := Module[{ik, klist, NumKpoint, spinupdn},
  klist = Flatten[ParseXMLData[ParseXML[xml, "varray", {{"name", "kpointlist"}}], "v"], 1];
  NumKpoint = Length[klist];
  spinupdn = <|"k" -> klist, "up" -> Table[Flatten[ParseXMLData[ParseXML[ParseXML[ParseXML[xml, "projected", {}], "set", {{"comment", "spin 1"}}], "set", {{"comment", "kpoint " <> ToString[ik]}}], "r"], 1], {ik, NumKpoint}], "dn" -> Table[Flatten[ParseXMLData[ParseXML[ParseXML[ParseXML[xml, "projected", {}], "set", {{"comment", "spin 2"}}], "set", {{"comment", "kpoint " <> ToString[ik]}}], "r"], 1], {ik, NumKpoint}]|>;
  Return[spinupdn]
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

VaspBSPlot[bsxml_, klabels_, OptionsPattern[{"ERange" -> All, "FigRatio"->1/GoldenRatio}]] := Module[{kpath, updata, dndata, upplot, dnplot},
  kpath = {Accumulate[Join[{0}, Norm[#] & /@ Differences[bsxml["k"]]]], bsxml["k"]}\[Transpose];
  updata = {kpath\[Transpose][[1]], #\[Transpose][[1]]}\[Transpose] & /@ (bsxml["up"]\[Transpose]);
  dndata = {kpath\[Transpose][[1]], #\[Transpose][[1]]}\[Transpose] & /@ (bsxml["dn"]\[Transpose]);
  upplot = ListLinePlot[updata, 
                        PlotStyle -> {{Black, Thick}}, 
                        Joined -> True, 
                        PlotRange -> {All, OptionValue["ERange"]}, 
                        AspectRatio -> OptionValue["FigRatio"], 
                        Frame -> True, 
                        GridLines -> {{klabels\[Transpose][[1]], ConstantArray[Thick, Length[klabels]]}\[Transpose], Automatic},
                        FrameTicks -> {{Automatic, None}, {klabels, None}},
                        ImageSize -> Medium];
  dnplot = ListLinePlot[dndata,  
                        PlotStyle -> {{Black, Dashed, Thin}},  
                        Joined -> True,  
                        PlotRange -> {All, {-3, 2}},  
                        AspectRatio -> OptionValue["FigRatio"],  
                        Frame -> True, 
                        GridLines -> {{klabels\[Transpose][[1]], ConstantArray[Thick, Length[klabels]]}\[Transpose], Automatic},
                        FrameTicks -> {{Automatic, None}, {klabels, None}},
                        ImageSize -> Medium];

  Return[Show[{upplot,dnplot}]]
]
(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
