BeginPackage["LINVARIANT`Perovskite`",{"LINVARIANT`Structure`", "LINVARIANT`ISODISTORT`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
GetPmode           ::usage = "GetPmode[DspMode]"
GetX5mode          ::usage = "GetX5mode[DspMode]"
GetStrain          ::usage = "GetStrain[parent, sub]"
GenTraningSets     ::usage = "GenTraningSets[AllModes, num, PhaseName, type]"
GenStrainSets      ::usage = "GenStrainSets[LatticeVector, type, num]"
GoalFunction       ::usage = "GoalFunctionEnergy[func, volume, trainingset, coeff]"
GetElasticModuli   ::usage = "GetElasticModuli[file, volume]"
ReadInvariant      ::usage = "ReadInvariant[invariants_, \[CapitalGamma]4_]"
GetTensor          ::usage = "GetTensor[IvariantTerm]"
LandauInvariant    ::usage = "LandauInvariant[]"
MeshMask           ::usage = "MeshMask[i, Dim, x, y, z]"
SetInitial         ::usage = "SetInitial[vars, iconfig]"
Plotop             ::usage = "Plotop[Data, op]"
PlotPvector        ::usage = "PlotPvector[min]"
PlotSurface        ::usage = "PlotSurface[Data]"
MinimizeBoracite   ::usage = "MinimizeBoracite[GF, iconfig, Niter1, Etype, Efield, Niter2, OptionsPattern[{NumTimeStep -> 0}]]"
PathPlot           ::usage = "PathPlot[min, frame, y0, z0, param]"
CompareOP          ::usage = "CompareOP[pdata, min, frame, range]"
GetOpTables        ::usage = "GetOpTables[minlist, Dim]"
SimplifyBoracite   ::usage = "SimplifyBoracite[pos]"

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)
InvariantString = Association[
"3rd\[CapitalGamma]4" -> "3 n1 n2 n3",
"3rdX" -> "3   n4 n6 n9 - n5 n7 n8",
"3rd\[CapitalGamma]4X" -> "3   n1 n6 n7 + n2 n4 n5 + n3 n8 n9",
"4thX" -> 
"4   n4^4 + 2 n4^2 n5^2 + 2 n4^2 n6^2 + 2 n4^2 n7^2 + 2 n4^2 n8^2 + 2 n4^2 n9^2 + n5^4 + 2 n5^2 n6^2 + 2 n5^2 n7^2 + 2 n5^2 n8^2 + 2 n5^2 n9^2 + n6^4 + 2 n6^2 n7^2 + 2 n6^2 n8^2 + 2 n6^2 n9^2 + n7^4 + 2 n7^2 n8^2 + 2 n7^2 n9^2 + n8^4 + 2 n8^2 n9^2 + n9^4
4   n4^4 + n5^4 + n6^4 + n7^4 + n8^4 + n9^4
4   n4^2 n5^2 + n6^2 n7^2 + n8^2 n9^2
4   n4^2 n6^2 + n4^2 n9^2 + n5^2 n7^2 + n5^2 n8^2 + n6^2 n9^2 + n7^2 n8^2
4   n4^2 n7^2 + n5^2 n9^2 + n6^2 n8^2",
"4th\[CapitalGamma]4X" -> 
"4   n1^2 n4^2 + n1^2 n5^2 + n1^2 n6^2 + n1^2 n7^2 + n1^2 n8^2 + n1^2 n9^2 + n2^2 n4^2 + n2^2 n5^2 + n2^2 n6^2 + n2^2 n7^2 + n2^2 n8^2 + n2^2 n9^2 + n3^2 n4^2 + n3^2 n5^2 + n3^2 n6^2 + n3^2 n7^2 + n3^2 n8^2 + n3^2 n9^2
4   n1^2 n4^2 + n1^2 n8^2 + n2^2 n7^2 + n2^2 n9^2 + n3^2 n5^2 + n3^2 n6^2
4   n1^2 n5^2 + n1^2 n9^2 + n2^2 n6^2 + n2^2 n8^2 + n3^2 n4^2 + n3^2 n7^2
4   n1 n2 n8 n9 + n1 n3 n4 n5 + n2 n3 n6 n7
4   n1 n4 n7 n9 - n1 n5 n6 n8 - n2 n4 n7 n8 + n2 n5 n6 n9 + n3 n4 n6 n8 - n3 n5 n7 n9",
"6thX" -> 
"6   n4^6 + n5^6 + n6^6 + n7^6 + n8^6 + n9^6
6   n4^4 n5^2 + n4^2 n5^4 + n6^4 n7^2 + n6^2 n7^4 + n8^4 n9^2 + n8^2 n9^4
6   n4^4 n6^2 + n4^2 n9^4 + n5^4 n8^2 + n5^2 n7^4 + n6^4 n9^2 + n7^2 n8^4
6   n4^4 n7^2 + n4^2 n7^4 + n5^4 n9^2 + n5^2 n9^4 + n6^4 n8^2 + n6^2 n8^4
6   n4^4 n8^2 + n4^2 n8^4 + n5^4 n6^2 + n5^2 n6^4 + n7^4 n9^2 + n7^2 n9^4"];

\[CapitalDelta]1; \[CapitalDelta]2; \[CapitalDelta]3;
NumTsteps; ntstep;
x; y; z;
aa; AA; bb; BB; cc; CC;
n1; n2; n3; n4; n5; n6; n7; n8; n9;
Subscript[\[Epsilon], 1, 1]; Subscript[\[Epsilon], 2, 2]; Subscript[\[Epsilon], 3, 3]; 
Subscript[\[Epsilon], 2, 3]; Subscript[\[Epsilon], 1, 3]; Subscript[\[Epsilon], 1, 2];
Subscript[P, 1]; Subscript[P, 2]; Subscript[P, 3];
Subscript[u, 1]; Subscript[u, 2]; Subscript[u, 3];
Subscript[\[Alpha], xx]; 
Subscript[\[Alpha], xxxx]; Subscript[\[Alpha], xxXX]; Subscript[\[Alpha], xxyy]; Subscript[\[Alpha], xxYY]; 
Subscript[\[Alpha], xxxxxx]; Subscript[\[Alpha], xxxxXX]; Subscript[\[Alpha], xxxxyy]; Subscript[\[Alpha], xxxxYY]; Subscript[\[Alpha], XXXXyy];
Subscript[\[Beta], 11];
Subscript[\[Beta], 1111]; Subscript[\[Beta], 1122]; 
Subscript[\[Beta], 111111]; Subscript[\[Beta], 111122]; Subscript[\[Beta], 112233];
Subscript[C, 1111]; Subscript[C, 1122]; Subscript[C, 1212];
Subscript[Q, 1111]; Subscript[Q, 1122]; Subscript[Q, 1212];
Subscript[Gx, 1111]; Subscript[Gx, 1122]; Subscript[Gx, 1212];
Subscript[Gp, 1111]; Subscript[Gp, 1122]; Subscript[Gp, 1212];
Subscript[\[Gamma], x5];
Subscript[\[Gamma], p];
Subscript[\[Gamma], c];
Subscript[\[Gamma], pc];
Subscript[\[Gamma], xc];
Subscript[\[Gamma], xp];
Subscript[\[Eta], c1c2p3]; Subscript[\[Eta], a2b1p3]; Subscript[\[Eta], a1b2p3];

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)

GetPmode[DspMode_] := Module[{GM4NoneZero, Pmode, Pindex, Px, Py, Pz},
  Pindex = Partition[{43, 44, 45, 55, 56, 57, 67, 68, 69, 313, 314, 315, 340, 341, 342, 343, 344, 345, 412, 413, 414, 415, 416, 417}, 3]\[Transpose];
  GM4NoneZero = If[StringCases[#[[2]], RegularExpression["GM4"]] != {}, #, ## &[]] & /@ DspMode;
  Px = Sign[DspMode[[340]][[4]]] Norm[(If[MemberQ[Pindex[[1]], #[[1]]], #, ## &[]] & /@ GM4NoneZero)\[Transpose][[4]]];
  Py = Sign[DspMode[[341]][[4]]] Norm[(If[MemberQ[Pindex[[2]], #[[1]]], #, ## &[]] & /@ GM4NoneZero)\[Transpose][[4]]];
  Pz = Sign[DspMode[[342]][[4]]] Norm[(If[MemberQ[Pindex[[3]], #[[1]]], #, ## &[]] & /@ GM4NoneZero)\[Transpose][[4]]];
  Return[{Px, Py, Pz}]
]

GetX5mode[DspMode_] := Module[{X5string, X5a, X5b, X5c, X5d, X5e, X5f},
  X5string = "X5(\\()[abcdef,]+(\\))(\[)[[:upper:]]([[:lower:]]?)\\d+:\\w:dsp(\])\\w+(\\()";
  X5a = Sign[DspMode[[379]][[4]]] Norm[(If[StringCases[#[[2]], RegularExpression[X5string <> "a" <> "(\\))"]] != {}, #, ## &[]] & /@ DspMode)\[Transpose][[4]]];
  X5b = Sign[DspMode[[380]][[4]]] Norm[(If[StringCases[#[[2]], RegularExpression[X5string <> "b" <> "(\\))"]] != {}, #, ## &[]] & /@ DspMode)\[Transpose][[4]]];
  X5c = Sign[DspMode[[381]][[4]]] Norm[(If[StringCases[#[[2]], RegularExpression[X5string <> "c" <> "(\\))"]] != {}, #, ## &[]] & /@ DspMode)\[Transpose][[4]]];
  X5d = Sign[DspMode[[382]][[4]]] Norm[(If[StringCases[#[[2]], RegularExpression[X5string <> "d" <> "(\\))"]] != {}, #, ## &[]] & /@ DspMode)\[Transpose][[4]]];
  X5e = Sign[DspMode[[383]][[4]]] Norm[(If[StringCases[#[[2]], RegularExpression[X5string <> "e" <> "(\\))"]] != {}, #, ## &[]] & /@ DspMode)\[Transpose][[4]]];
  X5f = Sign[DspMode[[384]][[4]]] Norm[(If[StringCases[#[[2]], RegularExpression[X5string <> "f" <> "(\\))"]] != {}, #, ## &[]] & /@ DspMode)\[Transpose][[4]]];
  Return[{X5a, X5b, X5c, X5d, X5e, X5f}]
]

GetStrain[parent_, sub_] := Module[{e, strain}, 
  e = N[sub[[6]] /. sub[[8]]].Inverse[N[parent[[6]] /. parent[[8]]]] - IdentityMatrix[3];
  strain = 1/2 (e + Transpose[e] + Transpose[e].e);
  Return[{strain[[1, 1]], strain[[2, 2]], strain[[3, 3]], strain[[1, 2]], strain[[2, 3]], strain[[1, 3]]}]
]

GenTraningSets[AllModes_, num_, PhaseName_, type_] := Module[{i, str, N, ModesNum, DirName, ModeTable, X5NoneZero, GM4NoneZero, Pmode, Pxmode, Pymode, Pzmode, Pindex, X5string, X5a1,X5a2, X5b1, X5b2, X5c1, X5c2},
  Pindex = Partition[{43, 44, 45, 55, 56, 57, 67, 68, 69, 313, 314, 315, 340, 341, 342, 343, 344, 345, 412, 413, 414, 415, 416, 417}, 3]\[Transpose];
  X5NoneZero = If[StringCases[#[[2]], RegularExpression["X5"]] != {} && Abs[#[[4]]] > 0.0001, #, ## &[]] & /@ AllModes;
  X5string = "X5(\\()[abcdef,]+(\\))(\[)[[:upper:]]([[:lower:]]?)\\d+:\\w:dsp(\])\\w+(\\()";
  {X5a1, X5a2, X5b1, X5b2, X5c1, X5c2} = Table[If[StringCases[#[[2]], RegularExpression[X5string <> str <> "(\\))"]] != {}, #, ## &[]] & /@ AllModes, {str, {"a", "b", "c", "d", "e", "f"}}];
  GM4NoneZero = If[StringCases[#[[2]], RegularExpression["GM4"]] != {} && Abs[#[[4]]] > 0.000, #, ## &[]] & /@ AllModes;
  Pmode = If[MemberQ[Flatten[Pindex], #[[1]]], #, ## &[]] & /@ GM4NoneZero;
  {Pxmode, Pymode, Pzmode} = Table[If[MemberQ[Pindex[[str]], #[[1]]], #, ## &[]] & /@ AllModes, {str, 3}];
  DirName = PhaseName <> "-" <> type;
  Print[DirName<>" Generated."];
  N = Quotient[num, 2];
  Which[
   type == "xp",
   ModesNum = Join[Pmode\[Transpose][[1]], X5NoneZero\[Transpose][[1]]];
   ModeTable = Table[If[MemberQ[ModesNum, #[[1]]], {#[[1]], #[[2]], #[[3]], (1 + i/100) #[[4]]}, {#[[1]], #[[2]], #[[3]], #[[4]]}] & /@ AllModes, {i, -N, N}],
   type == "x",
   ModesNum = Join[X5NoneZero\[Transpose][[1]]];
   ModeTable = Table[If[MemberQ[ModesNum, #[[1]]], {#[[1]], #[[2]], #[[3]], (1 + i/100) #[[4]]}, {#[[1]], #[[2]], #[[3]], #[[4]]}] & /@ AllModes, {i, -N, N}],
   type == "p",
   ModesNum = Join[Pmode\[Transpose][[1]]];
   ModeTable = Table[If[MemberQ[ModesNum, #[[1]]], {#[[1]], #[[2]], #[[3]], (1 + i/100) #[[4]]}, {#[[1]], #[[2]], #[[3]], #[[4]]}] & /@ AllModes, {i, -N, N}],
   type == "p1p2",
   ModesNum = {Pxmode\[Transpose][[1]], Pymode\[Transpose][[1]]};
   ModeTable = Table[If[Position[ModesNum, #[[1]]] != {}, {#[[1]], #[[2]], #[[3]],Sign[i]^Position[ModesNum, #[[1]]][[1]][[1]] i/500}, {#[[1]], #[[2]], #[[3]], #[[4]]}] & /@ AllModes, {i, -N, N}],
   type == "p2",
   ModesNum = Join[Pymode\[Transpose][[1]]];
   ModeTable = Table[If[MemberQ[ModesNum, #[[1]]], {#[[1]], #[[2]], #[[3]], i/500}, {#[[1]], #[[2]], #[[3]], #[[4]]}] & /@ AllModes, {i, -N, N}],
   type == "p3",
   ModesNum = Join[Pzmode\[Transpose][[1]]];
   ModeTable = Table[If[MemberQ[ModesNum, #[[1]]], {#[[1]], #[[2]], #[[3]], i/500}, {#[[1]], #[[2]], #[[3]], #[[4]]}] & /@ AllModes, {i, -N, N}],
   type == "a2b1",
   ModesNum = {X5a2\[Transpose][[1]], X5b1\[Transpose][[1]]};
   ModeTable = Table[If[Position[ModesNum, #[[1]]] != {}, {#[[1]], #[[2]], #[[3]], Sign[i]^Position[ModesNum, #[[1]]][[1]][[1]] 2*i/5000}, {#[[1]], #[[2]], #[[3]], #[[4]]}] & /@ AllModes, {i, -N, N}],
   type == "a1b2",
   ModesNum = {X5a1\[Transpose][[1]], X5b2\[Transpose][[1]]};
   ModeTable = Table[If[Position[ModesNum, #[[1]]] != {}, {#[[1]], #[[2]], #[[3]], Sign[i]^Position[ModesNum, #[[1]]][[1]][[1]] 2*i/5000}, {#[[1]], #[[2]], #[[3]], #[[4]]}] & /@ AllModes, {i, -N, N}],
   type == "c1c2",
   ModesNum = {X5c1\[Transpose][[1]], X5c2\[Transpose][[1]]};
   ModeTable = Table[If[Position[ModesNum, #[[1]]] != {}, {#[[1]], #[[2]], #[[3]], Sign[i]^Position[ModesNum, #[[1]]][[1]][[1]] 2*i/5000}, {#[[1]], #[[2]], #[[3]], #[[4]]}] & /@ AllModes, {i, -N, N}],
   type == "a1b1c2",
   ModesNum = Join[X5a1\[Transpose][[1]], X5b1\[Transpose][[1]], X5c2\[Transpose][[1]]];
   ModeTable = Table[If[MemberQ[ModesNum, #[[1]]], {#[[1]], #[[2]], #[[3]], 2*i/5000}, {#[[1]], #[[2]], #[[3]], #[[4]]}] & /@ AllModes, {i, -N, N}]
  ]; 
  Return[{ModeTable, DirName}]
]

GenStrainSets[LatticeVector_, PhaseName_, type_, num_] := Module[{i, N, R, DirName, StrainTable},
  DirName = type <> "-" <> PhaseName;
  Print[DirName<>" Generated."];
  N = Quotient[num, 2];
  R = Normal[Symmetrize[LatticeVector]];
  Which[
   type == "strain",
   StrainTable = Table[R + DiagonalMatrix[i/1000 (Diagonal@R)], {i, -N, N}],
   type == "Q1111",
   StrainTable = Table[DiagonalMatrix[{1.0, 1.0, 1.0 + i/1000}].R, {i, -N, N}],
   type == "Q1122",
   StrainTable = Table[DiagonalMatrix[{1.0 + i/1000, 1.0 + i/1000, 1.0}].R, {i, -N, N}],
   type == "Q1212",
   StrainTable = Table[{{1.0, i/1000, 0}, {i/1000, 1.0, 0}, {0, 0, 1.0}}.R, {i, -N, N}],
   type == "gammaxc",
   StrainTable = Table[{{1.0, i/1000, 0}, {i/1000, 1.0, 0}, {0, 0, 1.0}}.R, {i, -N, N}],
   type == "gammac",
   StrainTable = Table[{{1.0, 2 i/10000, 2 i/10000}, {2 i/10000, 1.0, 2 i/10000}, {2 i/10000, 2 i/10000, 1.0}}.R, {i, -N, N}]
   ];
  Return[{StrainTable, DirName}]
]

GoalFunction[func_, volume_, ts_, coe_] := Module[{},
  Total[((func /. coe /. Thread[{aa, AA, bb, BB, CC, cc, Subscript[P, 1], Subscript[P, 2], Subscript[P, 3], Subscript[\[Epsilon], 1, 1], Subscript[\[Epsilon], 2, 2], Subscript[\[Epsilon], 3, 3], Subscript[\[Epsilon], 1, 2], Subscript[\[Epsilon], 2, 3], Subscript[\[Epsilon], 1, 3]} -> #[[3]]]) - #[[2]]/volume)^2 & /@ ts]
]

GetElasticModuli[file_, volume_] := Module[{KBar, KBar2meV, Cijkl},
  KBar2meV = 0.1*1000/160.21766208;
  KBar = Select[Import[file, {"Data"}], UnsameQ[#, {}] &];
  Cijkl = SparseArray[{{i_, i_, i_, i_} -> Subscript[C, 1111],
                       {i_, i_, j_, j_} /; i != j -> Subscript[C, 1122],
                       {i_, j_, i_, j_} /; i != j -> Subscript[C, 1212]}, {3, 3, 3, 3}];
  DeleteDuplicates@Flatten[Table[Which[i == j && k == l, Cijkl[[i, j, k, l]] -> KBar[[i, k]] KBar2meV, i == k && j == l && i != j, Cijkl[[i, j, k, l]] -> KBar[[i + 3, i + 3]] KBar2meV, True, ## &[]], {i, 3}, {j, 3}, {k, 3}, {l, 3}]]
]   

ReadInvariant[invariants_, \[CapitalGamma]4_] := Module[{R, X5, ops, IsoVarsX5, IsoVars\[CapitalGamma]4, IsoInvOrder, InvariantTerms},
  R = {x, y, z};
  X5 = {aa, AA, bb, BB, CC, cc};
  IsoVarsX5 = {n4, n5, n6, n7, n8, n9};
  IsoVars\[CapitalGamma]4 = {n1, n2, n3};
  Which[\[CapitalGamma]4 == "ShearStrains", ops = {Subscript[\[Epsilon], 2, 3], Subscript[\[Epsilon], 1, 3], Subscript[\[Epsilon], 1, 2]}, 
        \[CapitalGamma]4 == "Polarizations", ops = {Subscript[P, 1], Subscript[P, 2], Subscript[P, 3]}
       ];
  {IsoInvOrder, InvariantTerms} = ReadList[StringToStream[invariants], {Number, Expression}]\[Transpose];
  InvariantTerms = # /. Thread[IsoVarsX5 -> Through[X5 @@ R]] /. Thread[IsoVars\[CapitalGamma]4 -> Through[ops @@ R]] & /@ InvariantTerms;
  Return[InvariantTerms]
]

GetTensor[IvariantTerm_] := Module[{CoeffRules, X5, tensor, sa, XX},
  X5 = {aa, AA, bb, BB, CC, cc};
  tensor = Table[CoeffRules = CoefficientRules[ter, X5];
    sa = SparseArray[
      Table[Flatten[Join[ConstantArray[Position[Thread[ru[[1]] == 6], True], 3], 
                         ConstantArray[Position[Thread[ru[[1]] == 4], True], 2], 
                         ConstantArray[Position[Thread[ru[[1]] == 2], True], 1]]] -> ru[[2]], {ru, CoeffRules}], 6]; 
    1/2 (sa + Transpose[sa]), {ter, IvariantTerm /. {XX_[x, y, z] -> XX}}];
  Return[tensor]
]

LandauInvariant[type_, OptionsPattern[{"MeshType"->"square"}]] := Module[{Boracite, R, r, i, j, k, l, m, n, ndim, Pk, Pl, Pi2, Pj2, Pk2, Pij, Pkl, Pmn, \[Epsilon]ij, \[Epsilon]kl, \[Epsilon]mn, epsilon, X, X5, X5k, X5l, X5i2, X5j2, X5k2, X5ij, X5kl, X5mn, \[CapitalGamma]4P, \[CapitalGamma]4\[Epsilon], \[Alpha]i2, TensorX5, TensorX56th, \[Alpha]i2\[Alpha]j2, \[Alpha]i2\[Alpha]j2\[Alpha]k2, \[Beta]i2, \[Beta]iijj, \[Beta]iijjkk, Cijkl, Cijmn, Cklmn, Qmnkl, qijkl, Gijkl, EX5termsX5\[CapitalGamma]4P, fx5, fx56, fg, fp, fc, fpc, fxc, fxp, \[Epsilon]2u, PDisc, X5Disc, zbc, xbc, fDisc0, fEfield},

(**********************************************-----Variables-----************************************************)
  R = {x, y, z};
  ndim = Length@R;
  Subscript[r, i_] := R[[i]];
  Pk = Pl = Array[Subscript[P, #] @@ R &, {3}];
  Pi2 = Pj2 = Pk2 = Array[(Subscript[P, #] @@ R)^2 &, {3}];
  Pij = Pkl = Pmn = Array[D[Subscript[P, #1] @@ R, Subscript[r, #2]] &, {3, ndim}];
  \[Epsilon]ij = \[Epsilon]kl = \[Epsilon]mn = SparseArray[{{i_, j_} /; i == j -> Subscript[\[Epsilon], i, j] @@ R, 
                                                            {i_, j_} /; i < j -> Subscript[\[Epsilon], i, j] @@ R, 
                                                            {i_, j_} /; i > j -> Subscript[\[Epsilon], j, i] @@ R}, {3, 3}];
  epsilon = Table[Subscript[\[Epsilon], i, j] @@ R -> Array[1/2 (D[Subscript[u, #1] @@ R, Subscript[r, #2]] + D[Subscript[u, #2] @@ R, Subscript[r, #1]]) &, {ndim, ndim}][[i, j]], {i, 1, ndim}, {j, 1, ndim}] // Flatten;
  X5 = {aa, AA, bb, BB, CC, cc};
  X5k = X5l = Array[X5[[#]] @@ R &, {6}];
  X5i2 = X5j2 = X5k2 = Array[(X5[[#]] @@ R)^2 &, {6}];
  X5ij = X5kl = X5mn = Array[D[X5[[#1]] @@ R, Subscript[r, #2]] &, {6, ndim}];
  X = Through[X5 @@ R];
  \[CapitalGamma]4P = {Subscript[P, 1], Subscript[P, 2], Subscript[P, 3]};
  \[CapitalGamma]4\[Epsilon] = {Subscript[\[Epsilon], 2, 3], Subscript[\[Epsilon], 1, 3], Subscript[\[Epsilon], 1, 2]};

(**********************************************-----Tensors-----************************************************)
  TensorX5 = GetTensor[ReadInvariant[InvariantString["4thX"], "Polarizations"]];
  TensorX56th = GetTensor[ReadInvariant[InvariantString["6thX"], "Polarizations"]];
  EX5termsX5\[CapitalGamma]4P = ReadInvariant[InvariantString["4th\[CapitalGamma]4X"], "Polarizations"];
  \[Alpha]i2 = SparseArray[{{i_} -> Subscript[\[Alpha], xx]}, {6}];
  \[Alpha]i2\[Alpha]j2 = Total[Normal[#1*#2] & @@ {TensorX5, {2 Subscript[\[Alpha], xxxx], 0, 2 Subscript[\[Alpha], xxXX], 2 Subscript[\[Alpha], xxyy], 2 Subscript[\[Alpha], xxYY]}}];
  \[Alpha]i2\[Alpha]j2\[Alpha]k2 = Total[Normal[#1*#2] & @@ {TensorX56th, {Subscript[\[Alpha], xxxxxx], Subscript[\[Alpha], xxxxXX], Subscript[\[Alpha], xxxxyy], Subscript[\[Alpha], xxxxYY], Subscript[\[Alpha], XXXXyy]}}];
  \[Beta]i2 = SparseArray[{{i_} -> Subscript[\[Beta], 11]}, {3}];
  \[Beta]iijj = SparseArray[{{i_, i_} -> 2 Subscript[\[Beta], 1111], {i_, j_} /; i != j -> Subscript[\[Beta], 1122]}, {3, 3}];
  \[Beta]iijjkk = SparseArray[{{i_, i_, i_} -> 2 Subscript[\[Beta], 111111], 
                               {i_, i_, j_} /; i != j -> 2/3 Subscript[\[Beta], 111122], 
                               {i_, j_, i_} /; i != j -> 2/3 Subscript[\[Beta], 111122], 
                               {j_, i_, i_} /; i != j -> 2/3 Subscript[\[Beta], 111122], 
                               {i_, j_, k_} /; i != j && i != k && j != k -> 1/3 Subscript[\[Beta], 112233]}, {3, 3, 3}];
  Cijkl = Cijmn = Cklmn = SparseArray[{{i_, i_, i_, i_} -> Subscript[C, 1111], 
                                       {i_, i_, j_, j_} /; i != j -> Subscript[C, 1122], 
                                       {i_, j_, i_, j_} /; i != j -> Subscript[C, 1212]}, {3, 3, 3, 3}];
  Qmnkl = SparseArray[{{m_, m_, m_, m_} -> Subscript[Q, 1111], 
                       {m_, m_, n_, n_} /; m != n -> Subscript[Q, 1122], 
                       {m_, n_, m_, n_} /; m != n -> Subscript[Q, 1212]}, {3, 3, 3, 3}];
  qijkl = TensorContract[2 Cijmn\[TensorProduct]Qmnkl, {{3, 5}, {4, 6}}];
  Gxijkl = SparseArray[{{i_, i_, i_, i_} -> Subscript[Gx, 1111], 
                       {i_, i_, j_, j_} /; i != j -> Subscript[Gx, 1122], 
                       {i_, j_, i_, j_} /; i != j -> Subscript[Gx, 1212]}, {6, 3,6, 3}];
  Gpijkl = SparseArray[{{i_, i_, i_, i_} -> Subscript[Gp, 1111],
                       {i_, i_, j_, j_} /; i != j -> Subscript[Gp, 1122],
                       {i_, j_, i_, j_} /; i != j -> Subscript[Gp, 1212]}, {6, 3,6, 3}];

(**********************************************-----EnergyOnSite-----************************************************)
  fx5 = TensorContract[\[Alpha]i2\[TensorProduct]X5i2, {1, 2}] 
        + TensorContract[1/2 \[Alpha]i2\[Alpha]j2\[TensorProduct]X5i2\[TensorProduct]X5j2,{{1, 3}, {2, 4}}] 
        + 0 TensorContract[\[Alpha]i2\[Alpha]j2\[Alpha]k2\[TensorProduct]X5i2\[TensorProduct]X5j2\[TensorProduct]X5k2, {{1, 4}, {2, 5}, {3, 6}}] 
        + Subscript[\[Gamma], x5] (aa[x, y, z] bb[x, y, z] cc[x, y, z] - AA[x, y, z] BB[x, y, z] CC[x, y, z]); 
  fg = TensorContract[1/2 X5ij\[TensorProduct]Gxijkl\[TensorProduct]X5kl, {{1, 3}, {2, 4}, {5, 7}, {6, 8}}]
     + TensorContract[1/2 Pij\[TensorProduct]Gpijkl\[TensorProduct]Pkl, {{1, 3}, {2, 4}, {5, 7}, {6, 8}}];
  fp = TensorContract[1/2 \[Beta]i2\[TensorProduct]Pi2, {{1, 2}}] 
     + TensorContract[1/2 \[Beta]iijj\[TensorProduct]Pi2\[TensorProduct]Pj2, {{1, 3}, {2, 4}}] 
     + 0 TensorContract[1/2 \[Beta]iijjkk\[TensorProduct]Pi2\[TensorProduct]Pj2\[TensorProduct]Pk2, {{1, 4}, {2, 5}, {3, 6}}] 
     + Subscript[\[Gamma], p] Subscript[P, 1][x, y, z] Subscript[P, 2][x, y, z] Subscript[P, 3][x, y, z];
  fx56 = Subscript[\[Alpha], xxxxXX] (aa[x, y, z]^2 AA[x, y, z]^2 (aa[x, y, z]^2 + AA[x, y, z]^2) + bb[x, y, z]^2 BB[x, y, z]^2 (bb[x, y, z]^2 + BB[x, y, z]^2) + CC[x, y, z]^2 cc[x, y, z]^2 (CC[x, y, z]^2 + cc[x, y, z]^2)) + Subscript[\[Alpha], xxxxxx] ((aa[x, y, z]^2 + AA[x, y, z]^2)^3 + (bb[x, y, z]^2 + BB[x, y, z]^2)^3 + (CC[x, y, z]^2 + cc[x, y, z]^2)^3);
  fc = TensorContract[1/2 \[Epsilon]ij\[TensorProduct]Cijkl\[TensorProduct]\[Epsilon]kl, {{1, 3}, {2, 4}, {5, 7}, {6, 8}}] 
     + 0 Subscript[\[Gamma], c] Subscript[\[Epsilon], 1, 2][x, y, z] Subscript[\[Epsilon], 2, 3][x, y, z] Subscript[\[Epsilon], 1, 3][x, y, z];
  fpc = TensorContract[-(1/2) \[Epsilon]ij\[TensorProduct]qijkl\[TensorProduct]Pk\[TensorProduct]Pl, {{1, 3}, {2, 4}, {5, 7}, {6, 8}}] + Subscript[\[Gamma], pc] {Subscript[\[Epsilon], 2, 3][x, y, z], Subscript[\[Epsilon], 1, 3][x, y, z], Subscript[\[Epsilon], 1, 2][x, y, z]}.{Subscript[P, 1][x, y, z], Subscript[P, 2][x, y, z], Subscript[P, 3][x, y, z]};
  fxc = Subscript[\[Gamma], xc] (Subscript[\[Epsilon], 2, 3][x, y, z] bb BB + Subscript[\[Epsilon], 1, 3][x, y, z] aa AA + Subscript[\[Epsilon], 1, 2][x, y, z] cc CC) /. Thread[X5 -> X];
  fxp = (Subscript[\[Gamma], xp] (Subscript[P, 1][x, y, z] bb BB + Subscript[P, 2][x, y, z] aa AA + Subscript[P, 3][x, y, z] cc CC) /. Thread[X5 -> X]) 
      + EX5termsX5\[CapitalGamma]4P.{Subscript[\[Eta], c1c2p3], Subscript[\[Eta], a2b1p3], Subscript[\[Eta], a1b2p3], 0, 0};
  Boracite = fx5 + fp + fc + fpc + fxc + fxp;
(**********************************************-----Discretization-----************************************************)
  \[Epsilon]2u = Dispatch[epsilon /. Subscript[\[Epsilon], i_, j_][x_, y_, z_] -> 
                                     Subscript[\[Epsilon], i, j, x, y, z] 
                                  /. Derivative[dx_, dy_, dz_][Subscript[u, i_]][x_, y_, z_] -> 
                                     (Subscript[u, i, x + dx, y + dy, z + dz] - Subscript[u, i, x, y, z])/(\[CapitalDelta]1 dx + \[CapitalDelta]2 dy + \[CapitalDelta]3 dz)
                                  /. Subscript[\[Epsilon], i_, j_, __] -> 
                                     Subscript[\[Epsilon], i, j, x_, y_, z_]];
  \[Epsilon]Disc = Dispatch[{Subscript[\[Epsilon], i_, j_][x_, y_, z_] -> Subscript[\[Epsilon], i, j, x, y, z]}];
  PDisc = Dispatch[{Subscript[P, i_][x_, y_, z_] -> 
                    Subscript[P, i, x, y, z], 
                    Derivative[dx_, dy_, dz_][Subscript[P, i_]][x_, y_, z_] -> 
                    (Subscript[P, i, x + dx, y + dy, z + dz] - Subscript[P, i, x, y, z])/(\[CapitalDelta]1 dx + \[CapitalDelta]2 dy + \[CapitalDelta]3 dz)}];
  X5Disc = Dispatch[Flatten@{Thread[X -> (Subscript[#, x, y, z] & /@ X5)], 
                    Thread[Derivative[dx_, dy_, dz_][#][x_, y_, z_] -> 
                           (Subscript[#, x + dx, y + dy, z + dz] - Subscript[#, x, y, z])/(\[CapitalDelta]1 dx + \[CapitalDelta]2 dy + \[CapitalDelta]3 dz)] & /@ X5}];
  zbc = Dispatch[{Subscript[u, 3, x_, y_, Lz + 1] :> Subscript[u, 3, x, y, 1]}];
  xbc = Dispatch[{Subscript[u, 3, Lx + 1, y_, z_] :> Subscript[u, 3, 1, y, z]}];

  fDisc0 = Boracite + fg /. \[Epsilon]Disc /. PDisc /. X5Disc /. \[Epsilon]2u;

  Which[type=="OnSite", Return[Boracite], type=="Ginzburg", Return[Boracite + fg], type=="Discretized", Return[fDisc0]]
]

MeshMask[i_, Dim_, x_, y_, z_] := Module[{s, Lx, Ly, Lz},
  {Lx, Ly, Lz} = Dim;
  s = Switch[i,
      "square",          If[0 < x <= Lx && 0 < y <= Ly && 0 < z <= Lz, 1, 0],
      "circle",          If[(x - Lx/2)^2 + (y - Ly/2)^2 < (Min[Lx, Ly]/2)^2, 1, 0],
      "TiltedRec",       If[x - Ly/4 <= y <= x + Ly/4 && -y + 7 Ly/4 >= x >= -y + Ly/4, 1, 0],
      "TiltedRecUp",     If[x - Ly/4 <= y <= x + Ly/4 && -y + 7 Ly/4 >= x >= -y + Ly/4, 1, 0],
      "TiltedRecDn",     If[x - Ly/4 <= y <= x + Ly/4 && -y + 7 Ly/4 >= x >= -y + Ly/4, 1, 0],
      "TiltedRecCenter", If[x - Ly/4 <= y <= x + Ly/4 && -y + 5 Ly/4 >= x >= -y + 3 Ly/4, 1,0],
      "Ehh",             If[x >= -y + 3 Ly/2 && x - Ly/8 <= y <= x + Ly/8, 1, 0],
      "Ett",             If[x <= -y + Ly/2 && x - Ly/8 <= y <= x + Ly/8, 1, 0],
      "Enarrow1",        If[x - Ly/8 <= y <= x + Ly/8, 1, 0],
      "Enarrow2",        If[-x - Ly/8 <= y <= -x + Ly/8, 1, 0],
      "EnarrowVertical", If[Lx/4 <= x <= Lx 3/4, 1, 0]];
  Return[s]
]

SetInitial[vars_, iconfig_] := Module[{init, v, op, x, y, z},
   RandomSeed[12356];
   Which[
    iconfig === "random",
    init = Table[{op = v[[1]], x = v[[2]] - Lx/2, y = v[[3]] - Ly/2, z = v[[4]] - Lz/2}; 
                 Which[op === cc || op === CC || op === u || op === P, 0, True, RandomReal[{-1, 1}]], {v, vars}],
    iconfig === "1",
    init = Table[{op = v[[1]], x = v[[2]] - Lx/2, y = v[[3]] - Ly/2, z = v[[4]] - Lz/2};
                 Which[op === aa && (y - x) (y + x) < 0, 1,
                       op === AA && (y - x) (y + x) < 0 && x > 0, 1,
                       op === AA && (y - x) (y + x) < 0 && x < 0, -1,
                       op === bb && (y - x) (y + x) > 0, 1,
                       op === BB && (y - x) (y + x) > 0 && y < 0, 1,
                       op === BB && (y - x) (y + x) > 0 && y > 0, -1, True, 0], {v, vars}],
    iconfig === "2",
    init = Table[{op = v[[1]], x = v[[2]] - Lx/2, y = v[[3]] - Ly/2, z = v[[4]] - Lz/2};
                 Which[op === aa && (y - x) (y + x) > 0, 1,
                       op === AA && (y - x) (y + x) > 0 && y < 0, 1,
                       op === AA && (y - x) (y + x) > 0 && y > 0, 1,
                       op === bb && (y - x) (y + x) < 0, 1,
                       op === BB && (y - x) (y + x) < 0 && x < 0, -1,
                       op === BB && (y - x) (y + x) < 0 && x > 0, 1, True, 0], {v, vars}],
    iconfig === "3",
    init = Table[{op = v[[1]], x = v[[2]] - Lx/2, y = v[[3]] - Ly/2, z = v[[4]] - Lz/2};
                 Which[op === aa && (y - x) (y + x) < 0, 1,
                       op === AA && (y - x) (y + x) < 0 && x < 0, 1,
                       op === AA && (y - x) (y + x) < 0 && x > 0, 1,
                       op === bb && (y - x) (y + x) > 0, 1,
                       op === BB && (y - x) (y + x) > 0 && y < 0, -1,
                       op === BB && (y - x) (y + x) > 0 && y > 0, 1, True, 0], {v, vars}],
    iconfig === "4",
    init = Table[{op = v[[1]], x = v[[2]] - Lx/2, y = v[[3]] - Ly/2, z = v[[4]] - Lz/2};
                 Which[op === aa && x > 0, 1, op === AA && x > 0 && y - x < 0, 1,
                       op === AA && x > 0 && y - x > 0, -1, op === bb && x < 0, 1,
                       op === BB && x < 0 && y + x > 0, 1,
                       op === BB && x < 0 && y + x < 0, -1, True, 0], {v, vars}],
    iconfig === "5",
    init = Table[{op = v[[1]], x = v[[2]] - Lx/2, y = v[[3]] - Ly/2, z = v[[4]] - Lz/2};
                 Which[op === aa && y + x > 0, 1,
                       op === AA && y + x > 0 && x < 0, 1,
                       op === AA && y + x > 0 && x > 0, -1, op === bb && y + x < 0, 1,
                       op === BB && y + x < 0 && y > 0, 1,
                       op === BB && y + x < 0 && y < 0, -1, True, 0], {v, vars}],
    iconfig === "6",
    init = Table[{op = v[[1]], x = v[[2]] - Lx/2, y = v[[3]] - Ly/2, z = v[[4]] - Lz/2};
                 Which[op === aa && y - x > 0, 1,
                       op === AA && y - x > 0 && x < 0, -1,
                       op === AA && y - x > 0 && x > 0, 1, op === bb && y - x < 0, 1,
                       op === BB && y - x < 0 && y < 0, -1,
                       op === BB && y - x < 0 && y > 0, 1, True, 0], {v, vars}],
    iconfig === "exp",
    init = Table[{op = v[[1]], x = v[[2]], y = v[[3]], z = v[[4]]};
                 Which[op === aa && y - x - 2/4 Ly > 0 && y - x - 3/4 Ly < 0, 1,
                       op === aa && y - x > 0 && y - x - 1/4 Ly < 0, 1,
                       op === aa && y - x + 2/4 Ly > 0 && y - x + 1/4 Ly < 0, 1,
                       op === aa && y - x + 3/4 Ly < 0, 1,
                       op === AA && y - x - 2/4 Ly > 0 && y - x - 3/4 Ly < 0 && x < 1/4 Lx, 1,
                       op === AA && y - x - 2/4 Ly > 0 && y - x - 3/4 Ly < 0 && x > 1/4 Lx, -1,
                       op === AA && y - x > 0 && y - x - 1/4 Ly < 0 && ((1/4 Lx < x < 2/4 Lx) || (x > 3/4 Lx)), 1,
                       op === AA && y - x > 0 && y - x - 1/4 Ly < 0 && ((x < 1/4 Lx) || (2/4 Lx < x < 3/4 Lx)), -1,
                       op === AA && y - x + 2/4 Ly > 0 && y - x + 1/4 Ly < 0 && 2/4 Lx < x < 3/4 Lx, 1,
                       op === AA && y - x + 2/4 Ly > 0 && y - x + 1/4 Ly < 0 && ((x < 2/4 Lx) || (x > 3/4 Lx)), -1,
                       op === AA && y - x + 3/4 Ly < 0, 1,
                       op === bb && y - x - 3/4 Ly > 0, 1,
                       op === bb && y - x - 1/4 Ly > 0 && y - x - 2/4 Ly < 0, 1,
                       op === bb && y - x + 1/4 Ly > 0 && y - x < 0, 1,
                       op === bb && y - x + 3/4 Ly > 0 && y - x + 2/4 Ly < 0, 1,
                       op === BB && y - x - 3/4 Ly > 0, 1,
                       op === BB && y - x - 1/4 Ly > 0 && y - x - 2/4 Ly < 0 && 2/4 Ly < y < 3/4 Ly, 1,
                       op === BB && y - x - 1/4 Ly > 0 && y - x - 2/4 Ly < 0 &&((y > 3/4 Ly) || (y < 2/4 Ly)), -1,
                       op === BB && y - x + 1/4 Ly > 0 && y - x < 0 && ((1/4 Ly < y < 2/4 Ly) || (y > 3/4 Ly)), 1,
                       op === BB && y - x + 1/4 Ly > 0 && y - x < 0 && ((2/4 Ly < y < 3/4 Ly) || (y < 1/4 Ly)), -1,
                       op === BB && y - x + 3/4 Ly > 0 && y - x + 2/4 Ly < 0 && y < 1/4 Ly, 1,
                       op === BB && y - x + 3/4 Ly > 0 && y - x + 2/4 Ly < 0 && y > 1/4 Ly, -1, 
                       op === P  && v[[2]] ==3, -1, True, 0], {v, vars}],
    iconfig === "monopole",
    init = Table[{op = v[[1]], x = v[[2]], y = v[[3]], z = v[[4]]};
                 Which[op === aa && y - x - 1/2 Ly > 0 && y + x - Ly < 0, 1,
                       op === AA && y - x - 1/2 Ly > 0 && y + x - Ly < 0, -1,
                       op === bb && y - x - 1/2 Ly > 0 && y + x - Ly > 0, 1,
                       op === BB && y - x - 1/2 Ly > 0 && y + x - Ly > 0, 1, 
                       op === aa && y - x > 0 && y - x - 1/2 Ly < 0 && y + x - 1/2 Ly < 0, 1,
                       op === AA && y - x > 0 && y - x - 1/2 Ly < 0 && y + x - 1/2 Ly < 0, 1,
                       op === bb && y - x > 0 && y - x - 1/2 Ly < 0 && y + x - 1/2 Ly > 0 && y + x - Ly < 0, 1,
                       op === BB && y - x > 0 && y - x - 1/2 Ly < 0 && y + x - 1/2 Ly > 0 && y + x - Ly < 0, -1,
                       op === aa && y - x > 0 && y - x - 1/2 Ly < 0 && y + x - Ly > 0 && y + x - 3/2 Ly < 0, 1,
                       op === AA && y - x > 0 && y - x - 1/2 Ly < 0 && y + x - Ly > 0 && y + x - 3/2 Ly < 0, 1,
                       op === bb && y - x > 0 && y - x - 1/2 Ly < 0 && y + x - 3/2 Ly > 0, 1,
                       op === BB && y - x > 0 && y - x - 1/2 Ly < 0 && y + x - 3/2 Ly > 0, -1,
                       op === bb && y - x < 0 && y - x + 1/2 Ly > 0 && y + x - 1/2 Ly < 0, 1,
                       op === BB && y - x < 0 && y - x + 1/2 Ly > 0 && y + x - 1/2 Ly < 0, 1,
                       op === aa && y - x < 0 && y - x + 1/2 Ly > 0 && y + x - 1/2 Ly > 0 && y + x - Ly < 0, 1,
                       op === AA && y - x < 0 && y - x + 1/2 Ly > 0 && y + x - 1/2 Ly > 0 && y + x - Ly < 0, -1,
                       op === bb && y - x < 0 && y - x + 1/2 Ly > 0 && y + x - Ly > 0 && y + x - 3/2 Ly < 0, 1,
                       op === BB && y - x < 0 && y - x + 1/2 Ly > 0 && y + x - Ly > 0 && y + x - 3/2 Ly < 0, 1,
                       op === aa && y - x < 0 && y - x + 1/2 Ly > 0 && y + x - 3/2 Ly > 0, 1,
                       op === AA && y - x < 0 && y - x + 1/2 Ly > 0 && y + x - 3/2 Ly > 0, -1,
                       op === bb && y - x + 1/2 Ly < 0 && y + x - Ly < 0, 1,
                       op === BB && y - x + 1/2 Ly < 0 && y + x - Ly < 0, -1,
                       op === aa && y - x + 1/2 Ly < 0 && y + x - Ly > 0, 1,
                       op === AA && y - x + 1/2 Ly < 0 && y + x - Ly > 0, 1, True, 0], {v, vars}],
    iconfig === "toroidal",
    init = Table[{op = v[[1]], x = v[[2]], y = v[[3]], z = v[[4]]};
                 Which[op === bb && y - x - 1/2 Ly > 0 && y + x - Ly < 0, 1,
                       op === BB && y - x - 1/2 Ly > 0 && y + x - Ly < 0, -1,
                       op === aa && y - x - 1/2 Ly > 0 && y + x - Ly > 0, 1,
                       op === AA && y - x - 1/2 Ly > 0 && y + x - Ly > 0, -1,
                       op === bb && y - x > 0 && y - x - 1/2 Ly < 0 && y + x - 1/2 Ly < 0, 1,
                       op === BB && y - x > 0 && y - x - 1/2 Ly < 0 && y + x - 1/2 Ly < 0, 1,
                       op === aa && y - x > 0 && y - x - 1/2 Ly < 0 && y + x - 1/2 Ly > 0 && y + x - Ly < 0, 1,
                       op === AA && y - x > 0 && y - x - 1/2 Ly < 0 && y + x - 1/2 Ly > 0 && y + x - Ly < 0, 1,
                       op === bb && y - x > 0 && y - x - 1/2 Ly < 0 && y + x - Ly > 0 && y + x - 3/2 Ly < 0, 1,
                       op === BB && y - x > 0 && y - x - 1/2 Ly < 0 && y + x - Ly > 0 && y + x - 3/2 Ly < 0, 1,
                       op === aa && y - x > 0 && y - x - 1/2 Ly < 0 && y + x - 3/2 Ly > 0, 1,
                       op === AA && y - x > 0 && y - x - 1/2 Ly < 0 && y + x - 3/2 Ly > 0, 1,
                       op === aa && y - x < 0 && y - x + 1/2 Ly > 0 && y + x - 1/2 Ly < 0, 1,
                       op === AA && y - x < 0 && y - x + 1/2 Ly > 0 && y + x - 1/2 Ly < 0, -1,
                       op === bb && y - x < 0 && y - x + 1/2 Ly > 0 && y + x - 1/2 Ly > 0 && y + x - Ly < 0, 1,
                       op === BB && y - x < 0 && y - x + 1/2 Ly > 0 && y + x - 1/2 Ly > 0 && y + x - Ly < 0, -1,
                       op === aa && y - x < 0 && y - x + 1/2 Ly > 0 && y + x - Ly > 0 && y + x - 3/2 Ly < 0, 1,
                       op === AA && y - x < 0 && y - x + 1/2 Ly > 0 && y + x - Ly > 0 && y + x - 3/2 Ly < 0, -1,
                       op === bb && y - x < 0 && y - x + 1/2 Ly > 0 && y + x - 3/2 Ly > 0, 1,
                       op === BB && y - x < 0 && y - x + 1/2 Ly > 0 && y + x - 3/2 Ly > 0, -1,
                       op === aa && y - x + 1/2 Ly < 0 && y + x - Ly < 0, 1,
                       op === AA && y - x + 1/2 Ly < 0 && y + x - Ly < 0, 1,
                       op === bb && y - x + 1/2 Ly < 0 && y + x - Ly > 0, 1,
                       op === BB && y - x + 1/2 Ly < 0 && y + x - Ly > 0, 1, True, 0], {v, vars}],
    iconfig === "wav",
    init = Table[{op = v[[1]], x = v[[2]], y = v[[3]], z = v[[4]]};
                 Which[op === aa, 1,
                       op === AA && 1/3 Lx < x < 2/3 Lx, 1,
                       op === AA && x < 1/3 Lx, -1,
                       op === AA && x > 2/3 Lx, -1,
                       True, 0], {v, vars}],
    iconfig === "hill",
    init = Table[{op = v[[1]], x = v[[2]], y = v[[3]], z = v[[4]]};
                 Which[op === cc && (x < 1/4 Lx) || (x > 3/4 Lx), 1,
                       op === CC && (x < 1/4 Lx) || (x > 3/4 Lx), 1,
                       op === aa && 1/4 Lx < x < 3/4 Lx, 1,
                       op === AA && 1/4 Lx < x < 2/4 Lx, -1,
                       op === AA && 2/4 Lx < x < 3/4 Lx, 1,
                       True, 0], {v, vars}]                       
    ];
   Return[init]
]

Plotop[Dict_, op_, Dim_, OptionsPattern[{"RangeType" -> "whole", "FrameNum" -> 1}]] := Module[{pl, range, disrange, iop, Data, TablePlot},
  {Lx, Ly, Lz} = Dim;
  Data = Dict[op][[2]];
  disrange = Which[OptionValue["RangeType"] == "whole", ConstantArray[Dict[op][[1]], Length@Data], OptionValue["RangeType"] == "separated", MinMax[#] & /@ Data];
  (*Print[disrange];*)
  Table[
   ListDensityPlot[Data[[iop]],
                   ColorFunctionScaling -> False,
                   ColorFunction -> (ColorData["Rainbow"][Rescale[#, disrange[[iop]]]] &), 
                   PlotLegends -> Placed[BarLegend[{(ColorData["Rainbow"][Rescale[#, range]] &) /. range -> disrange[[iop]], disrange[[iop]]}, 
                                         LabelStyle -> {Black, Bold, 10}, LegendLayout -> "Column", 
                                         LegendMargins -> 0, 
                                         LegendLabel -> Placed[Style[ToString[Subscript[op, iop]] <> "@" <> ToString[OptionValue["FrameNum"] - 1], 20], Top]], Right],
                   PlotRange -> disrange[[iop]],
                   InterpolationOrder -> 1,
                   ClippingStyle -> {Blue, Red},
                   AspectRatio -> Ly/Lx,
                   ImageSize -> Small
                   ], 
   {iop, Length@Data}]
]

PlotPvector[min_, Dim_, density_, OptionsPattern[{"FrameNum" -> 1, "meshtype"->"square"}]] := Module[{lla, x, y, z},
  Show[Graphics3D[
                  Join[
                   Transpose[Table[lla = 1 MeshMask[OptionValue["meshtype"], Dim, x, y, z]; 
                                   {Hue[Sign[lla Subscript[P, 3, x, y, z]]/4 + 0.5], 
                                    Arrowheads[0.02], 
                                    Arrow[Tube[{{x, y, z} - lla {Subscript[P, 1, x, y, z], 
                                                Subscript[P, 2, x, y, z], 
                                                Subscript[P, 3, x, y, z]}, 
                                               {x, y, z} + lla {Subscript[P, 1, x, y, z], 
                                               Subscript[P, 2, x, y, z], 
                                               Subscript[P, 3, x, y, z]}}]]}, 
                                   {x, 1, Lx, density}, {y, 1, Ly, density}, {z, 1, Lz}]] /. min, 
                   {Black, Thick, Dashed, Line[{{1, Ly, 1}, {Lx, 1, 1}}], Line[{{1, 1, 1}, {Lx, Ly, 1}}]}], 
                  ViewPoint -> {0, 0, 10}, 
                  ImageSize -> Large, 
                  PlotLabel -> "@" <> ToString[OptionValue["FrameNum"] - 1]], 
                  FrameTicks -> True]
]


PlotSurface[Data_] := Module[{},
  Show[ListPlot3D[Data["u"][[2]][[3]], 
       ColorFunction -> "Rainbow",
       PlotLegends -> Automatic, 
       PlotRange -> Automatic,
       AspectRatio -> Automatic, 
       ImageSize -> Small,
       ViewPoint -> {10, -10, 100}]]
]

PathPlot[min_, Dim_, frame_, y0_, z0_, param_] := Module[{xx, i, Lx, Ly, Lz},
  {Lx, Ly, Lz} = Dim;
  Grid[{Table[ListPlot[Chop[(Table[Subscript[P, i, xx, y0, z0], {xx, 1, Lx}] /. Dispatch[min[[#]]])&/@frame, 10^-6], 
                       PlotRange -> All, 
                       ImageSize -> Medium, 
                       Joined -> True, 
                       InterpolationOrder -> 1, 
                       PlotMarkers -> Automatic, 
                       PlotLabel -> Subscript[P, i], 
                       PlotLegends -> frame, 
                       Frame -> True], 
                {i, Range[3]}], 
            {ListPlot[Chop[
                (Table[(Subscript[u, 1, xx + 1, y0, z0] - Subscript[u, 1, xx, y0, z0])/\[CapitalDelta]1 /. param, {xx, 1, Lx-1}] /. Dispatch[min[[#]]])&/@frame, 10^-4], 
                       PlotRange -> All, 
                       ImageSize -> Medium, 
                       Joined -> True, 
                       InterpolationOrder -> 1, 
                       PlotMarkers -> Automatic, 
                       PlotLabel -> Subscript[\[Epsilon], 1, 1], 
                       PlotLegends -> frame, 
                       Frame -> True],
             ListPlot[Chop[
                 (Table[(Subscript[u, 2, xx, y0 + 1, z0] - Subscript[u, 2, xx, y0, z0])/\[CapitalDelta]1 /. param, {xx, 1, Lx-1}] /. Dispatch[min[[#]]])&/@frame, 10^-4],
                        PlotRange -> All, 
                        ImageSize -> Medium,
                        Joined -> True, 
                        InterpolationOrder -> 1, 
                        PlotMarkers -> Automatic, 
                        PlotLabel -> Subscript[\[Epsilon], 2, 2],
                        PlotLegends -> frame,
                        Frame -> True],
             ListPlot[Chop[
                 (Table[(Subscript[u, 3, xx, y0, z0 + 1] - Subscript[u, 3, xx, y0, z0])/\[CapitalDelta]1 /. param, {xx, 1, Lx-1}] /. Dispatch[min[[#]]])&/@frame, 10^-4],
                        PlotRange -> All,
                        ImageSize -> Medium,
                        Joined -> True,
                        InterpolationOrder -> 1,
                        PlotMarkers -> Automatic,
                        PlotLabel -> Subscript[\[Epsilon], 3, 3],
                        PlotLegends -> frame,
                        Frame -> True]},
            {ListPlot[Chop[
                (Table[1/2 (Subscript[u, 2, xx, y0, z0 + 1] - Subscript[u, 2, xx, y0, z0])/\[CapitalDelta]3 
                     + 1/2 (Subscript[u, 3, xx, y0 + 1, z0] - Subscript[u, 3, xx, y0, z0])/\[CapitalDelta]2 /. param, {xx, 1, Lx-1}] /. Dispatch[min[[#]]])&/@frame, 10^-4],
                       PlotRange -> All, 
                       ImageSize -> Medium,
                       Joined -> True, 
                       InterpolationOrder -> 1, 
                       PlotMarkers -> Automatic, 
                       PlotLabel -> Subscript[\[Epsilon], 2, 3],
                       PlotLegends -> frame,
                       Frame -> True],
             ListPlot[Chop[
                 (Table[1/2 (Subscript[u, 1, xx, y0, z0 + 1] - Subscript[u, 1, xx, y0, z0])/\[CapitalDelta]3
                      + 1/2 (Subscript[u, 3, xx + 1, y0, z0] - Subscript[u, 3, xx, y0, z0])/\[CapitalDelta]1 /. param, {xx, 1, Lx-1}] /. Dispatch[min[[#]]])&/@frame, 10^-4],
                        PlotRange -> All,
                        ImageSize -> Medium,
                        Joined -> True,
                        InterpolationOrder -> 1,
                        PlotMarkers -> Automatic,
                        PlotLabel -> Subscript[\[Epsilon], 1, 3],
                        PlotLegends -> frame,
                        Frame -> True],
             ListPlot[Chop[
                 (Table[1/2 (Subscript[u, 1, xx, y0 + 1, z0] - Subscript[u, 1, xx, y0, z0])/\[CapitalDelta]2
                      + 1/2 (Subscript[u, 2, xx + 1, y0, z0] - Subscript[u, 2, xx, y0, z0])/\[CapitalDelta]1 /. param, {xx, 1, Lx-1}] /. Dispatch[min[[#]]])&/@frame, 10^-4],
                        PlotRange -> All,
                        ImageSize -> Medium,
                        Joined -> True,
                        InterpolationOrder -> 1,
                        PlotMarkers -> Automatic,
                        PlotLabel -> Subscript[\[Epsilon], 1, 2],
                        PlotLegends -> frame,
                        Frame -> True]}
        }
    ]
]

CompareOP[pdata_, min_, frame_, Dim_, range_, OptionsPattern[{"VectorDensity"->1}]] := Module[{},
  Print[Grid@Flatten[(Partition[#, 6] & /@ Table[Flatten[Plotop[pdata[[n]], #, Dim, "RangeType" -> range] & /@ {"X5", "P", "u"}], {n, frame}])\[Transpose], 1]];
  Print[Grid@{Flatten[Table[{PlotPvector[Dispatch[min[[n]]], Dim, OptionValue["VectorDensity"]], PlotSurface[pdata[[n]]]}, {n, frame}]\[Transpose]]}];
]

MinimizeBoracite[iconfig_, Dim_, GF_, Niter1_, Etype_, Efields_, Niter2_, OptionsPattern[{"TimeSeries"->True, "dir"->""}]] := Module[{vars, fEfields, efield, init, minlist = {}, min, minEfield, XTable, PTable, uTable, PlotData, op, upx, n, t, x, y, z, BoundarySub},
  {Lx, Ly, Lz} = Dim;
  init = If[Length@iconfig==0, vars=Variables@GF;{vars, SetInitial[vars, iconfig]}\[Transpose], iconfig];
  XTable = Table[op, {op, {Subscript[aa, x, y, z], Subscript[AA, x, y, z], Subscript[bb, x, y, z], Subscript[BB, x, y, z], Subscript[CC, x, y, z], Subscript[cc, x, y, z]}}, {x, Lx}, {y, Ly}, {z, Lz}];
  PTable = Table[op, {op, {Subscript[P, 1, x, y, z], Subscript[P, 2, x, y, z], Subscript[P, 3, x, y, z]}}, {x, Lx}, {y, Ly}, {z, Lz}];
  uTable = Table[op, {op, {Subscript[u, 1, x, y, z], Subscript[u, 2, x, y, z], Subscript[u, 3, x, y, z]}}, {x, Lx}, {y, Ly}, {z, Lz}];

  Print["Minimizing without electric field."];
  t = Timing[min = FindMinimum[GF, init, MaxIterations -> Niter1];];
  BoundarySub = If[#1[[1]] === u && MemberQ[Mod @@@ ({{#1[[3]], #1[[4]]}, {Lx, Ly}}\[Transpose]), 1], #1 -> #2, Unevaluated[Sequence[]]] & @@@ min[[2]];
  Print["Initializing minimization took " <> ToString[First@t] <> " seconds."];
  init = {#1, #2} & @@@ (min[[2]]);
  init = If[#1[[1]] === u && MemberQ[Mod @@@ ({{#1[[3]], #1[[4]]}, {Lx, Ly}}\[Transpose]), 1], Unevaluated[Sequence[]], {#1, #2}] & @@@ (min[[2]]);
  AppendTo[minlist, min[[2]]];
  If[OptionValue["dir"]!="", Export[OptionValue["dir"]<>"/"<>OptionValue["dir"]<>".1"<>".txt",{#1,#2}&@@@min[[2]],"Table"], Unevaluated[Sequence[]]];

  If[Efields!={},
     fEfields =  Table[Sum[MeshMask[Etype, {Lx, Ly, Lz}, x, y, z] MeshMask["square", {Lx, Ly, Lz}, x, y, z] 
                        efield.{Subscript[P, 1, x, y, z], Subscript[P, 2, x, y, z], Subscript[P, 3, x, y, z]}, {x, Lx}, {y, Ly}, {z, Lz}], {efield, Efields}];
     Print["Applying electric field "];
     Do[t = Timing[minEfield = FindMinimum[(GF + fEfields[[n]])/.Dispatch[BoundarySub], init, MaxIterations -> Niter2];];
        Print[ToString[n]<>": Electric time step took " <> ToString[First@t] <> " seconds."];
        If[OptionValue["TimeSeries"], init = {#1, #2} & @@@ (minEfield[[2]]), Unevaluated[Sequence[]]];
        AppendTo[minlist, Join[minEfield[[2]], BoundarySub]];
        If[OptionValue["dir"]!="", 
           Export[OptionValue["dir"]<>"/"<>OptionValue["dir"]<>"."<>ToString[n+1]<>".txt",{#1,#2}&@@@minEfield[[2]],"Table"], 
           Unevaluated[Sequence[]]], 
        {n, Length@fEfields}
       ];
     Unevaluated[Sequence[]]
  ];

  PlotData = Table[Association[Thread[{"X5", "P", "u"} -> Table[{{-1, 1}*Max[Abs[upx /. Dispatch[min]]], Transpose[#[[;; , ;; , Lz]]] /. Dispatch[min] & /@ upx}, {upx, {XTable, PTable, uTable}}]]], {min, minlist}];
  EmitSound@Play[Sin[300 t Sin[20 t]], {t, 0, 3}];
  Return[{PlotData, minlist}]
]

GetOpTables[minlist_, Dim_]:=Module[{Lx, Ly, Lz, XTable, PTable, uTable, PlotData, min, op, upx, x, y, z, msub},
  {Lx, Ly, Lz} = Dim;
  XTable = Table[op, {op, {Subscript[aa, x, y, z], Subscript[AA, x, y, z], Subscript[bb, x, y, z], Subscript[BB, x, y, z], Subscript[CC, x, y, z], Subscript[cc, x, y, z]}}, {x, Lx}, {y, Ly}, {z, Lz}];
  PTable = Table[op, {op, {Subscript[P, 1, x, y, z], Subscript[P, 2, x, y, z], Subscript[P, 3, x, y, z]}}, {x, Lx}, {y, Ly}, {z, Lz}];
  uTable = Table[op, {op, {Subscript[u, 1, x, y, z], Subscript[u, 2, x, y, z], Subscript[u, 3, x, y, z]}}, {x, Lx}, {y, Ly}, {z, Lz}];

  PlotData = Table[msub=Flatten[{#[[1]]->#[[2]]}&/@min];Association[Thread[{"X5", "P", "u"} -> Table[{{-1, 1}*Max[Abs[upx /. Dispatch[msub]]], Transpose[#[[;; , ;; , Lz]]] /. Dispatch[msub] & /@ upx}, {upx, {XTable, PTable, uTable}}]]], {min, minlist}];
  Return[PlotData]
]

SimplifyBoracite[pos_] := Module[{SimplifyPos},
  SimplifyPos = If[StringCases[#2, RegularExpression["Cu|Cl"]] != {}, {#1, #2}, ## &[]] & @@@ pos;
  Return[SimplifyPos]
]


(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
