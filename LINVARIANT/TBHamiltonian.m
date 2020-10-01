BeginPackage["LINVARIANT`TBHamiltonian`", {"LINVARIANT`Structure`", "LINVARIANT`GroupTheory`", "LINVARIANT`INVARIANT`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
GetTi0Bonds                  ::usage "GetTi0Bonds[spg0, tij, latt, AllSites]"
ReadWannier90                ::usage "ReadWannier90[filename]"
TBHk                         ::usage "TBHk[Ti0, pos, k]"
TBBandsPlot                  ::usage "TBBandsPlot[Ti0, klist, kintv, OptionsPattern[{"range" -> All}]]"
TBBandsSpinCharacterPlot     ::usage "TBBandsSpinCharacterPlot[Ti0, klist, kintv, si]"
TBPlotSpinTextureSurface     ::usage "TBGetSpinTextureSurface[Ti0, latt, E0, klimit, kintv, tol, si]"
BondsToTi0                   ::usage "BondsToTi0[Bonds, ExtBonds, Atoms]"
Hso                          ::usage "Hso[l]"
CheckTBInvariants            ::usage "CheckTBInvariants[def, rules, sub]"
PlotTBInvariants             ::usage "PlotTBInvariants[pos, thop]"
GetAtomTBPair                ::usage "GetAtomTBPair[tb, atomlist]"

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
  AllBonds = Flatten[Table[{GrpV[xyz, AllSites[[#1[[1]], 1]]], GrpV[xyz, AllSites[[#1[[2]], 1]] + #2], #3}, {xyz, Keys[spg0]}] & @@@ tij, 1];
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

Hso[l_] := Module[{\[Mu], \[Tau], \[Theta], m, m1, m2, h11, h22, h12, h21},
  \[Tau][m_] := If[m >= 0, 1, 0];
  \[Theta][m_] := (1-KroneckerDelta[m,0])*(\[Tau][m]+I \[Tau][-m])/Sqrt[2]+KroneckerDelta[m,0]/2;
  \[Mu][m1_,m2_] := Conjugate[\[Theta][m1]]*\[Theta][m2](1+\[Tau][-Abs[m1*m2]]);
  h11 = Table[I*m2*KroneckerDelta[m1,-m2]/2, {m1,-l,l}, {m2,-l,l}];
  h22 = -h11;
  h12 = Table[\[Mu][m1,m2](KroneckerDelta[Abs[m1],Abs[m2]+1]
                           -(-1)^(\[Tau][m1]+\[Tau][m2])*KroneckerDelta[Abs[m1],Abs[m2]-1])*
                           Sqrt[(l+m1^2-Abs[m1*m2])(l+m2^2-Abs[m1*m2])]/2, {m1,-l,l}, {m2,-l,l}];
  h21 = Table[\[Mu][m1,m2](KroneckerDelta[Abs[m1],Abs[m2]-1]
                           -(-1)^(\[Tau][m1]+\[Tau][m2])*KroneckerDelta[Abs[m1],Abs[m2]+1])*
                           Sqrt[(l+m1^2-Abs[m1*m2])(l+m2^2-Abs[m1*m2])]/2, {m1,-l,l}, {m2,-l,l}];
 hso = ArrayFlatten[{{h11, h12}, {h21, h22}}];
 Return[hso]
]

CheckTBInvariants[def_, rules_, sub_] := Module[{s, orb, ind, vars, tb, invariants},
  vars = Table[
    s = Which[def[[i, 4]] == "up", "\[UpArrow]", def[[i, 4]] == "dn", "\[DownArrow]"];
    orb = Which[def[[i, 2]] == 1, Subscript[ToExpression["p"], def[[i, 1]], def[[i, 3]], s], 
                def[[i, 2]] == 2, Subscript[ToExpression["d"], def[[i, 1]], def[[i, 3]], s]];
    ind = First@First@Position[Values[sub], orb];
    Subscript[ToExpression["eIso"], ind], {i, Length[def]}];
  invariants = Expand[Total[Times @@ vars /. # & /@ rules] /. sub];
  tb = If[invariants===0, 0, SimplifyCommonFactor[invariants]];
  Return[tb]
]

PlotTBInvariants[pos_, thop_] := Module[{tb, latt, sites, SiteData, BondData, LabelData, tt, HoppingPairs},
  {latt, sites} = pos;
  tb = If[Head[thop]===Plus, Level[thop, 1], {thop}];
  HoppingPairs = Table[{First@First@Position[sites\[Transpose][[2]], #2], #3, #4} & @@@ Flatten[Which[NumberQ[#], ##&[], MatchQ[#, _Power]&&NumberQ[Level[#,1][[1]]], ##&[], MatchQ[#, _Subscript], #, MatchQ[#, _Power], ConstantArray[Level[#, 1][[1]], Level[#, 1][[2]]]] &/@ Which[MatchQ[tt, _Times], Level[tt, 1], MatchQ[tt, _Power], {tt}]], {tt, tb}];
  SiteData = {ElementData[SimplifyElementSymbol[#2], "IconColor"], Specularity[White, 20], Sphere[latt.#1, QuantityMagnitude[ElementData[SimplifyElementSymbol[#2], "AtomicRadius"], "Angstroms"]/2]} & @@@ sites;
  BondData = {Gray, Tube[{latt.sites[[#1[[1]], 1]], latt.sites[[#2[[1]], 1]]}]} & @@@ HoppingPairs;
  LabelData = Text@@# & /@ Flatten[{{ToString[Style[sites[[#1[[1]], 2]] <> ":", Medium, Bold, Black], StandardForm]<>ToString[Style[#1[[2]] <> #1[[3]], Medium, Bold, Red], StandardForm], latt.sites[[#1[[1]], 1]]},
                                     {ToString[Style[sites[[#2[[1]], 2]] <> ":", Medium, Bold, Black], StandardForm]<>ToString[Style[#2[[2]] <> #2[[3]], Medium, Bold, Red], StandardForm], latt.sites[[#2[[1]], 1]]}} & @@@ HoppingPairs, 1];
  Graphics3D[Join[SiteData, BondData, LabelData],
                    ImageSize -> 500, Axes -> True,
                    AxesLabel -> (Style[#, Bold, 64] & /@ {"a", "b", "c"}),
                    ViewPoint -> {0, 0, \[Infinity]}]
]

GetAtomTBPair[tb_, atomlist_] := Module[{SitePairs},
  SitePairs = {#1[[2]], #2[[2]]} & @@@ (Flatten[Which[NumberQ[#], ## &[], MatchQ[#, _Power] && NumberQ[Level[#, 1][[1]]], ## &[], MatchQ[#, _Subscript], #, MatchQ[#, _Power], ConstantArray[Level[#, 1][[1]], Level[#, 1][[2]]]] & /@ Which[MatchQ[#, _Times], Level[#, 1], MatchQ[#, _Power], {#}]] & /@ Level[tb, 1]);
  Table[If[Complement[SitePairs[[i]], atomlist] == {}, {i, Level[tb, 1][[i]]}, ## &[]], {i, Length@SitePairs}]
]

(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
