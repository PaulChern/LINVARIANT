BeginPackage["LINVARIANT`TBHamiltonian`", {"LINVARIANT`Structure`", "LINVARIANT`GroupTheory`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
GetTi0Bonds                  ::usage "GetTi0Bonds[spg0, tij, latt, AllSites]"
ReadWannier90                ::usage "ReadWannier90[filename]"
TBHk                         ::usage "TBHk[Ti0, pos, k]"
TBBandsPlot                  ::usage "TBBandsPlot[Ti0, klist, kintv, OptionsPattern[{"range" -> All}]]"
TBBandsSpinCharacterPlot     ::usage "TBBandsSpinCharacterPlot[Ti0, klist, kintv, si]"
TBPlotSpinTextureSurface     ::usage "TBGetSpinTextureSurface[Ti0, latt, E0, klimit, kintv, tol, si]"
BondsToTi0                   ::usage "BondsToTi0[Bonds, ExtBonds, Atoms]"
(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
GetTi0Bonds[spg0_, tij_, pos_] := Module[{latt, AllSites, c1, c2, c3, sol, site, AllBonds, ReducedBonds}, 
  {latt, AllSites} = pos;
  AllBonds = Flatten[Table[{GrpxV[xyz, AllSites[[#1[[1]], 1]]], GrpxV[xyz, AllSites[[#1[[2]], 1]] + #2], #3}, {xyz, Keys[spg0]}] & @@@ tij, 1];
  ReducedBonds = Flatten[Table[sol = Rationalize@Chop@First@Values@Solve[site[[2]] + IdentityMatrix[3].{c1, c2, c3} == site[[1]]];If[AllTrue[sol, IntegerQ], {#1 - sol, #2 - sol, #3}, ## &[]], {site, Tuples[{{#1, #2}, AllSites\[Transpose][[1]]}]}] & @@@ AllBonds, 1];
  ReducedBonds = DeleteDuplicates[If[Rationalize[Chop[#1 - Mod[#1, 1]]] == {0, 0, 0}, {#1, #2, #3}, {#2, #1, #3\[Transpose]}] & @@@ ReducedBonds, Chop[#1[[1]] - #2[[1]]] == {0, 0, 0} && Chop[#1[[2]] - #2[[2]]] == {0, 0, 0} &];
  ReducedBonds = {pos2index[AllSites, {#1, Mod[#2, 1]}], Rationalize[#2 - Mod[#2, 1]], 1, #3} & @@@ ReducedBonds;
  Return[ReducedBonds]]

ReadWannier90[filename_] := Module[{Wannier90, SystemName, NBands, NumCells, WeightList, Ti0}, 
  Wannier90 = OpenRead[filename];
  SystemName = Read[Wannier90, String];
  NBands = Read[Wannier90, Number];
  NumCells = Read[Wannier90, Number];
  WeightList = Table[Read[Wannier90, Number], {NumCells}];
  Ti0 = {#[[1]][[1 ;; 3]], Normal[SparseArray[{#4, #5} -> #6 + I #7 & @@@ #, {NBands, NBands}]]} & /@ Partition[ReadList[Wannier90, Number, RecordLists -> True], NBands^2];
  Return[{Ti0\[Transpose][[1]], WeightList, Ti0\[Transpose][[2]]}\[Transpose]]
]

TBHk[Ti0_, pos_, k_] := Module[{latt, sites, i, j, TBBlock, NumPos},
  {latt, sites} = pos;
  NumPos = Length[sites];
  TBBlock = Merge[{#1 -> #4/#3 Exp[I 2 Pi k.(sites[[#1[[2]], 1]] + #2 - sites[[#1[[1]], 1]])]} & @@@ Ti0, Total];
  ArrayFlatten[TensorTranspose[Table[If[MemberQ[Keys[TBBlock],{i,j}], TBBlock[{i, j}], 0.0 First@Values[TBBlock]], {i, NumPos}, {j, NumPos}], {3, 4, 1, 2}]]
]

TBBandsPlot[Ti0_, pos_, klist_, kintv_, OptionsPattern[{"range" -> All}]] := Module[{pdata, k, TBHij, sol, kpath, xticks, BandsPlot},
  {kpath, xticks} = GetKpath[klist, kintv];
  pdata = Table[TBHij = TBHk[Ti0, pos, k[[2]]];
                sol = Chop@Eigenvalues[TBHij];
                {k[[1]], #} & /@ sol, {k, kpath}];
  BandsPlot = ListPlot[pdata\[Transpose],
                       PlotStyle -> Black,
                       Joined -> False, 
                       PlotRange -> {All, OptionValue["range"]}, 
                       AspectRatio -> 1/GoldenRatio, 
                       Frame -> True, 
                       ImageSize -> Medium, 
                       GridLines -> {{xticks\[Transpose][[1]], ConstantArray[Thick, Length[xticks]]}\[Transpose], Automatic}, 
                       ImageSize -> Medium, 
                       FrameTicks -> {{Automatic, None}, {xticks, None}}
                 ];
  Return[BandsPlot]
]

TBBandsSpinCharacterPlot[Ti0_, pos_, klist_, kintv_, si_, OptionsPattern[{"range" -> All}]] := Module[{pdata, k, TBHij, sol, kpath, xticks, SpinOp, TBDim, SpinCharacterPlot, BandsPlot},
  {kpath, xticks} = GetKpath[klist, kintv];
  TBDim = Length[TBHk[Ti0, pos, {0, 0, 0}]];
  SpinOp = ArrayFlatten[PauliMatrix[#]\[TensorProduct]IdentityMatrix[TBDim/2]] & /@ Range[3];
  pdata = Table[TBHij = TBHk[Ti0, pos, k[[2]]];
                sol = {ToString[Chop@#1] -> Chop[#2\[Conjugate].SpinOp[[3]].#2]} & @@@ (Eigensystem[TBHij]\[Transpose]);
                {k[[1]], ToExpression[#1], Rescale[Chop[#2], {-1,1}]} & @@@ Normal[Merge[Flatten[sol], Total]], {k, kpath}];
  (*SpinCharacterPlot = BubbleChart[Flatten[pdata, 1], 
                                  Frame -> True,
                                  PlotRange -> {All, OptionValue["range"]},
                                  ColorFunction -> ({Opacity[#3], Blend[{Blue, White, Red}, #3]} &),
                                  BubbleSizes -> {0.03, 0.03},
                                  ChartStyle -> Directive[EdgeForm[None]], 
                                  ImageSize -> Medium,
                                  AspectRatio -> 1/GoldenRatio,
                                  FrameTicks -> {{Automatic, None}, {xticks, None}},
                                  GridLines -> {{xticks\[Transpose][[1]], ConstantArray[Thick, Length[xticks]]}\[Transpose], Automatic}];*)

  SpinCharacterPlot = Graphics[{Blend[{Blue, White, Red}, #3], Opacity[#3], PointSize[0.015], Point[{#1, #2}]} & @@@ Flatten[pdata, 1], Frame -> True, AspectRatio -> 1/GoldenRatio, ImageSize -> Medium];

  BandsPlot = ListPlot[Flatten[pdata, 1]\[Transpose][[1;;2]]\[Transpose],
                       PlotStyle -> Black,
                       Joined -> False,
                       PlotRange -> {All, OptionValue["range"]},
                       AspectRatio -> 1/GoldenRatio,
                       Frame -> True,
                       ImageSize -> Medium,
                       GridLines -> {{xticks\[Transpose][[1]], ConstantArray[Thick, Length[xticks]]}\[Transpose], Automatic},
                       ImageSize -> Medium,
                       FrameTicks -> {{Automatic, None}, {xticks, None}}
                  ];
  Return[Show[BandsPlot, SpinCharacterPlot]]
]

TBPlotSpinTextureSurface[Ti0_, pos_, E0_, klimit_, kintv_: 0.01, tol_: 0.1, si_, OptionsPattern[{"srange"->1}]] := Module[{latt, sites, TBDim, SpinOp, kx, ky, xlimit, ylimit, spindata, TBHij, sol, SpinTexture, TexturePlot, srange},
  {latt, sites} = pos;
  {xlimit, ylimit} = klimit;
  srange = OptionValue["srange"];
  TBDim = Length[TBHk[Ti0, pos, {0, 0, 0}]];
  SpinOp = ArrayFlatten[PauliMatrix[#]\[TensorProduct]IdentityMatrix[TBDim/2]] & /@ Range[3];
  spindata = Table[TBHij = TBHk[Ti0, pos, {kx, ky, 0}];
                   sol = {ToString[Chop@#1] -> Chop[#2\[Conjugate].SpinOp[[si]].#2]} & @@@ (Eigensystem[TBHij]\[Transpose]);
                   {ToExpression[#1], kx, ky, Rescale[Chop[#2], {-srange, srange}]} & @@@ Normal[Merge[Flatten[sol], Total]], {kx, -xlimit, xlimit, kintv}, {ky, -ylimit, ylimit, kintv}];
  SpinTexture = If[Chop[#1 - E0, tol] == 0, Join[2 Pi (Inverse[latt\[Transpose]].{#2, #3, 0})[[1 ;; 2]], {#4}], ## &[]] & @@@ Flatten[spindata, 2];
  TexturePlot = Graphics[{Blend[{Blue, White, Red}, #3], Opacity[#3], PointSize[0.01], Point[{#1, #2}]} & @@@ SpinTexture, Frame -> True, AspectRatio -> 1/(Divide @@ (Inverse[latt\[Transpose]].{xlimit, ylimit, 0})[[1 ;; 2]]), ImageSize -> Medium];
  Return[TexturePlot]
]

BondsToTi0[Bonds_, ExtBonds_, Atoms_] := Module[{NumAtom, Ti0},
  NumAtom = Length[Atoms\[Transpose][[1]]];
  Ti0 = {Rationalize[#2 - Mod[#2, 1]], 1, ArrayFlatten[#3\[TensorProduct]Normal[SparseArray[{pos2index[Atoms, #1], pos2index[Atoms, #2]} -> 1, {NumAtom, NumAtom}]]]} & @@@ Join[Bonds, ExtBonds];
  Return[Ti0]
]

(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
