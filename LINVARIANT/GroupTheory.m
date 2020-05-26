BeginPackage["LINVARIANT`GroupTheory`"]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
CifImportSpg              ::usage "CifImportSpg[file]"
GroupQ                    ::usage "GroupQ[grp]"
GetGroupK                 ::usage "GetGroupK[grp0, vec0, lattice0]"
GetStarK                  ::usage "GetStarK[grp0, k]"
GetSiteSymmetry           ::usage "GetSiteSymmetry[grp0, vec0]"
CifImportOperators        ::usage "CifImportOperators[file]"
GetGenerators             ::usage "GetGenerators[grp]"
Generators2Group          ::usage "Generators2Group[gen]"
xyzStr2Expression         ::usage "xyzStr2Expression[xyzStrData]"
xyz2Rot                   ::usage "xyz2Rot[expr]"
xyzStr2M4                 ::usage "xyzStr2M4[xyzStrData]"
RotTran2M4                ::usage "RotTran2M4[rot, tran]"
M42xyzStr                 ::usage "M42xyzStr[m4]"
M42xyzStrEle              ::usage "M42xyzStrEle[m4]"
GetClasses                ::usage "GetClasses[GrpMat]"
GetCharacters             ::usage "GetCharacters[GrpMat]"
ExpandSimplify            ::usage "ExpandSimplify[in]"
Mat2EulerVector           ::usage "Mat2EulerVector[mat]"
EulerVector2Mat           ::usage "EulerVector2Mat[axis, ang, \[Epsilon]]"
Mat2EulerAngles           ::usage "Mat2EulerAngles[latt, mat]"
GetEleLabel               ::usage "GetEleLabel[StrMat]"
GetRegRep                 ::usage "GetRegularRepresentation[grp]"
GetEleOrder               ::usage "GetEleOrder[mat]"
GetIrep                   ::usage "GetIrep[grp, ct, p]"
ModM4                     ::usage "ModM4[m1, m2]"
GrpMultiply               ::usage "GrpMultiply[lgrp, rgrp]"
GrpxV                     ::usage "GrpxV[grp, v]"
SortByOrder               ::usage "SortByOrder[grp]"
GetElePosition            ::usage "GetElePosition[grp, ele]"
GetSubGroups              ::usage "GetSubGroups[grp, ord]"
GetInvSubGroups           ::usage "GTInvSubGroups[grp]"
GetSpgCosetRepr           ::usage "GetSpgCosetRepr[grp0, invsubgrp]"
SolidSphericalHarmonicY   ::usage "SolidSphericalHarmonicY[l, m, x1, x2, x3, coord]"
SolidTesseralHarmonicY    ::usage "SolidTesseralHarmonicY[l, m, x1, x2, x3, coord]"
GetAngularMomentumRep     ::usage "GetAngularMomentumRep[latt, mat, l, Harmonic]"
GetEulerRodrigues         ::usage "GetEulerRodrigues[n, \[Theta], \[Epsilon]]"
GInverse                  ::usage "GInverse[ele]"
GTimes                    ::usage "GTimes[list]"
GetLeftCosets             ::usage "GetLeftCosets[grp0, invsub]"
(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
ModM4[m0_] := Module[{m},
  m = m0;
  If[Length[m0] == 4, m[[1 ;; 3, 4]] = Mod[m0[[1 ;; 3, 4]], 1], Unevaluated[Sequence[]]];
  Return[m]
]

ExpandSimplify[in_] := Expand[Simplify[in]];

GetWyckoffImages[grp_, pos_] := Module[{xyzStrData, xyzExpData, WyckoffImages},
  xyzStrData = Keys[grp];
  xyzExpData = Table[ToExpression["{" <> xyzStrData[[i]] <> "}"], {i, Length[xyzStrData]}];
  WyckoffImages = Flatten[Table[Sort[Union[{Mod[#,1], pos[[i, 2]]} & /@ (xyzExpData /. Thread[ToExpression[{"x", "y", "z"}] -> pos[[i, 1]]])], Total[#1[[1]]] > Total[#2[[1]]] &], {i, Length[pos]}], 1];
  Return[WyckoffImages];
]

GetClasses[grp0_] := Module[{i, j, g, InvGrpMat, bab, classes, Cout, GrpSub, GrpMat},
  GrpMat = Values[grp0];
  GrpSub = Dispatch[Thread[GrpMat -> Keys[grp0]]];
  g = Length[grp0];
  InvGrpMat = GInverse[grp0];
  bab = ParallelTable[ModM4[GrpMat[[i]].GrpMat[[j]].InvGrpMat[[i]]], {j, 1, g}, {i, 1, g}, DistributedContexts -> {"LINVARIANT`GroupTheory`Private`"}] /. GrpSub;
  classes = Union[Map[Union, bab]];
  Cout = RotateLeft[classes, Position[classes, {"x,y,z"}][[1]] - 1];
Return[Cout]
]

GetCharacters[grp0_, OptionsPattern["kpoint"->"\[CapitalGamma]"]] := Module[{m, n, i, j, k, h, g, nc, Classes, H0Mat, H1Mat, HMat, Hi, hh, Vt, V, Vi, cc, vec, IrrepDim, CharacterTable, CharacterTableSorted, ct, IrrepLabel, ClassesLabel, GrpSub},
  GrpSub = Thread[Keys[grp0] -> Values[grp0]];
  ClassesLabel = GetClasses[grp0];
  Classes = ClassesLabel /. GrpSub;
  h = Length[Classes];
  g = Length[grp0];
  nc = Map[Length, Classes];

  If[g > 1,
    (*---calculate the class constant matrices---*)
    H0Mat = Table[Table[Flatten[Table[ExpandSimplify[Classes[[j, m]].Classes[[i, n]]], {m, 1, nc[[j]]}, {n, 1, nc[[i]]}], 1], {i, 1, h}], {j, 1, h}];
    H1Mat = Table[Table[Length[Position[H0Mat[[j, i]], ExpandSimplify[Classes[[k, 1]]]]], {j, 1, h}, {i, 1, h}], {k, 1, h}];
    HMat = Table[H1Mat[[k, i, j]], {i, 1, h}, {j, 1, h}, {k, 1, h}];

    (*---find the matrix V that diagonalizes the class constant matrices---*)
    hh = Transpose[Sum[SetAccuracy[RandomReal[], 50]*HMat[[i]], {i, 2, h}]];
    Vt = Chop[Eigenvectors[hh]];
    V = ToRadicals@RootApproximant[Vt];
    Vi = ToRadicals@RootApproximant[Inverse[Vt]];

    (*---calculate the class constants---*)
    cc = Transpose@Table[Hi = HMat[[k]];
    vec = Diagonal[V.Hi.Vi], {k, 1, h}];

    (*---calculate the dimensions of the irreducible representations---*)
    IrrepDim = Sqrt[Table[g/(cc[[k]]/nc).Conjugate[cc[[k]]], {k, 1, h}]];

    (*---calculate the character table---*)
    CharacterTable = FullSimplify@Table[Table[cc[[i, j]]*IrrepDim[[i]]/nc[[j]], {j, 1, h}], {i, 1, h}];

    (*---sortby---*)
    CharacterTableSorted = SortBy[CharacterTable, First];
    ct = Prepend[Delete[CharacterTableSorted, First@Flatten@Position[CharacterTableSorted, Table[1, {i, 1, Length[CharacterTable]}]]], Table[1, {i, 1, Length[CharacterTable]}]];,
    ct = {{1}};
  ];


  IrrepLabel = Table["\!\(\*SubscriptBox[\("<>OptionValue["kpoint"]<>"\), \("<>ToString[y]<>"\)]\)", {y, 1, h}];
  Print[TableForm[ct, TableHeadings -> {IrrepLabel, ToString[Length[#]]<>GetEleLabel[First[#]] & /@ ClassesLabel}, TableAlignments -> {Right, Center}]];
  Return[{ClassesLabel, ct, IrrepLabel}]
]

CifImportOperators[file_] := Module[{CifData, CifFlags, xyzName, xyzStrData},
  CifData = Import[file, "string"] <> "\n";(*---fix for some strange behaviour---*)
  CifData = StringReplace[CifData, Thread[DeleteDuplicates[StringCases[CifData, RegularExpression[";"]]] -> ","]]; 
  CifData = StringReplace[CifData, Thread[DeleteDuplicates[StringCases[CifData, RegularExpression["[[:upper:]].{3,6}(\[)\\d,\\d,\\d(\])"]]] -> ""]];
  
  CifData = ImportString[CifData, "CIF"];
  CifFlags = Table[Level[CifData[[i]], 1][[1]], {i, Length[CifData]}];
  
  xyzName = Part[CifFlags, Flatten[Position[StringCount[CifFlags, "xyz"], 1]]];
  xyzStrData = StringDelete[#, " "] & /@ Flatten[xyzName /. CifData];
  
  Return[xyzStrData]
]

xyzStr2Expression[xyzStrData_] := Module[{xyzRotTranData, xyzTranslation, xyzRotData},
  xyzRotTranData = ToExpression["{" <> xyzStrData <> "}"]; 
  xyzTranslation = xyzRotTranData /. {ToExpression["x"] -> 0, ToExpression["y"] -> 0, ToExpression["z"] -> 0}; 
  xyzRotData = xyzRotTranData - xyzTranslation;
  Return[{xyzRotData, xyzTranslation}]
]

xyzStr2M4Ele[xyzStrData_] := Module[{RT, rot, tran},
  {rot, tran} = xyzStr2Expression[xyzStrData];
  RT = Join[Join[xyz2Rot[rot]\[Transpose], {tran}]\[Transpose], {{0., 0., 0., 1.}}];
  Return[Rationalize[RT]]
]

xyzStr2M4[xyzStrData_] := Module[{},
  Which[StringQ[xyzStrData], xyzStr2M4Ele[xyzStrData], ListQ[xyzStrData], xyzStr2M4[#] &/@ xyzStrData, True, Print["Wrong Input type!"];Abort[]]
]

RotTran2M4[rot_, tran_] := Module[{RT},
  RT = Join[Join[rot\[Transpose], {tran}]\[Transpose], {{0., 0., 0., 1.}}];
  Return[Rationalize[RT, 10^-12]]
]

M42xyzStrEle[m4_] := Module[{i, rot, tran},
  rot = ToString[#] & /@ (Rationalize[ModM4[m4].{ToExpression["x"], ToExpression["y"], ToExpression["z"], 1}] - (Rationalize[ModM4[m4].{ToExpression["x"], ToExpression["y"], ToExpression["z"], 1}] /. {ToExpression["x"] -> 0, ToExpression["y"] -> 0, ToExpression["z"] -> 0}));
  tran = (Rationalize[ModM4[m4].{ToExpression["x"], ToExpression["y"], ToExpression["z"], 1}] /. {ToExpression["x"] -> 0, ToExpression["y"] -> 0, ToExpression["z"] -> 0});
  tran = Which[#1 == 0, "", 
               Abs[#1] < Abs[#2], If[Sign[#1] > 0, "+", Unevaluated[Sequence[]]] <> ToString[#1] <> "/" <> ToString[#2], 
               Abs[#1] >= Abs[#2], If[Sign[#1] > 0, "+", Unevaluated[Sequence[]]] <> ToString[#1/#2]] 
         & @@@ (Through[{Numerator, Denominator}[#]] & /@ tran)[[1 ;; 3]];
  Return[StringDelete[StringJoin@Riffle[Table[rot[[i]] <> tran[[i]], {i, 3}],","], " "]]
]

M42xyzStr[m4_] := Module[{},
  If[MatrixQ[m4], M42xyzStrEle[m4], M42xyzStr[#] &/@ m4]
]

xyz2Rot[expr_] := Module[{m},
  m = Coefficient[#, {ToExpression["x"], ToExpression["y"], ToExpression["z"]}] &/@ expr;
  If[Abs[Det[m]] != 1, Print["xyz2Rot gives wrong determinant"]; Abort[], Unevaluate[Sequence[]]];
  Return[N[m]]
]

CifImportSpg[file_] := Module[{xyzStrData, ele, RT},
  xyzStrData = CifImportOperators[file];
  ele = xyzStrData;
  RT = xyzStr2M4[#] &/@ xyzStrData;
  Return[SortByOrder[Association[Thread[ele->RT]]]]
]

GetGroupK[grp0_, vec0_] := Module[{grp, go, trvec, eq, sol, add, grpout, c, vec},

  go = Length[grp0];
  grp = #[[1;;3,1;;3]] &/@ Values[grp0];

  (*--- Main algortihm ---*)
  trvec = Table[Inverse[grp[[i]]].vec0, {i, 1, go}];
  eq = Table[Thread[trvec[[i]] == vec0 + Table[c[i],{i,3}]], {i, 1,go}];
  sol = Rationalize[Round[Table[Table[c[j], {j, 3}] /. Flatten[Solve[eq[[i]], Table[c[j], {j, 3}]]], {i, 1, go}], 10^-3]];
  add = Table[IntegerQ[sol[[i, 1]]] && IntegerQ[sol[[i, 2]]] && IntegerQ[sol[[i, 3]]], {i, 1, go}];
  grpout = <||>;
  Do[If[add[[i]], grpout = Append[grpout, Keys[grp0][[i]]->Values[grp0][[i]]],None], {i, 1, go}];

  Return[grpout]
]

GetStarK[grp0_, k_] := Module[{grpk, LeftCoset, LeftRep},
  grpk = GetGroupK[grp0, k];
  LeftCoset = GetLeftCosets[grp0, grpk];
  LeftRep = LeftCoset[[;; , 1]];
  Return[Mod[GrpxV[LeftRep, k, "k"->True],1]]
]

GetSiteSymmetry[grp0_, vec0_] := Module[{grp, go, trvec, eq, sol, add, grpout, c, vec},

  go = Length[grp0];
  grp = #[[1;;3,1;;3]] &/@ Values[grp0];

  (*--- Main algortihm ---*)
  trvec = Table[Inverse[grp[[i]]].vec0, {i, 1, go}];
  eq = Table[Thread[trvec[[i]] == vec0 + Table[c[i],{i,3}]], {i, 1,go}];
  sol = Rationalize[Round[Table[Table[c[j], {j, 3}] /. Flatten[Solve[eq[[i]], Table[c[j], {j, 3}]]], {i, 1, go}], 10^-3]];
  add = Table[IntegerQ[sol[[i, 1]]] && IntegerQ[sol[[i, 2]]] && IntegerQ[sol[[i, 3]]], {i, 1, go}];
  grpout = <||>;
  Do[If[add[[i]], grpout = Append[grpout, Keys[grp0][[i]]->Values[grp0][[i]]],None], {i, 1, go}];

  Return[grpout]
]

EulerVector2Mat[latt_, axis_, ang_, \[Epsilon]_] := Module[{i, j, k, n, R},
  n = Normalize[latt.axis];
  R = Sum[Table[Cos[ang] KroneckerDelta[i, j] + (\[Epsilon] - Cos[ang]) n[[i]] n[[j]] - Sin[ang] LeviCivitaTensor[3][[i, j, k]] n[[k]], {i, 3}, {j, 3}], {k, 3}];
  Return[Rationalize[Inverse[latt].R.latt]]
]

Mat2EulerVector[mat_] := Module[{\[Phi], axes, es, s, \[Epsilon]},
  \[Epsilon] = Rationalize@Det[mat];

  \[Phi] = ArcCos[1/2 (Tr[\[Epsilon] mat] - 1)];
  If[\[Phi] == 0., Return[{{0, 0, 0}, 0}]];
  s = Sign@N[Chop[Im[Det[Eigensystem[\[Epsilon] mat][[2]]]]]];
  \[Phi] = If[s==0., \[Phi], s \[Phi]];
  es = Eigensystem[Rationalize[\[Epsilon] mat]];
  axes = es[[2]][[First@First@Position[es[[1]], 1]]];

Return[{axes, Rationalize[\[Phi]/Pi]*Pi}]
]

GetEulerRodrigues[n_, \[Theta]_, \[Epsilon]_:1] := Module[{K, R, nx, ny, nz},

  If[MatrixQ[n], K = n, 
    {nx, ny, nz} = n;
    K = {{0, -nz, ny}, {nz, 0, -nx}, {-ny, nx, 0}}
    ];

  If[\[Theta] == 0., 
    R = IdentityMatrix[3],
    R = Expand[(IdentityMatrix[3] + Sin[\[Theta]] K + (1 - Cos[\[Theta]]) K.K)]
    ];
  Return[\[Epsilon] IdentityMatrix[3].R]
]

Mat2EulerAngles[latt_, mat_] := Module[{\[Epsilon], axis, \[Phi], R, K, nx, ny, nz, \[Alpha], \[Beta], \[Gamma]},
  \[Epsilon] = Det[mat];
  {axis, \[Phi]} = Mat2EulerVector[mat];
  If[\[Phi] == 0., Return[{{0, 0, 0}, \[Epsilon]}]];
  {nx, ny, nz} = Normalize[latt.axis];
  K = {{0, -nz, ny}, {nz, 0, -nx}, {-ny, nx, 0}};
  (*R = Expand[(\[Epsilon] IdentityMatrix[3] + Sin[\[Phi]] K + \[Epsilon] (1 - \[Epsilon] Cos[\[Phi]]) K.K)];*)
  R = Expand[(IdentityMatrix[3] + Sin[\[Phi]] K + (1 - Cos[\[Phi]]) K.K)];
  {\[Alpha], \[Beta], \[Gamma]} = EulerAngles[R];
  (*\[Beta] = ArcCos[R[[3, 3]]];
  Which[
    \[Beta] == 0.,
      \[Alpha] = 0; \[Gamma] = \[Phi],
    \[Beta] == Pi,
      \[Alpha] = 0; \[Gamma] = Pi - \[Phi],
    True,
       \[Alpha] = ArcTan[nz nx Sin[\[Phi]/2] + ny Cos[\[Phi]/2], nz ny Sin[\[Phi]/2] - nx Cos[\[Phi]/2]];
       \[Gamma] = ArcTan[-nz nx Sin[\[Phi]/2] + ny Cos[\[Phi]/2], nz ny Sin[\[Phi]/2] + nx Cos[\[Phi]/2]];
      ];*)
  Return[{{\[Alpha], \[Beta], \[Gamma]}, \[Epsilon]}]
]

GetEleLabel[StrMat_] := Module[{det, mat, axes, phi, T, CS, label},
  Which[ListQ[StrMat]&&(!MatrixQ[StrMat]),
        GetEleLabel[#] &/@ StrMat,
        AssociationQ[StrMat],
        GetEleLabel[#] &/@ Keys[StrMat],
        StringQ[StrMat]||MatrixQ[StrMat],
        mat = If[StringQ[StrMat], xyzStr2M4[StrMat], StrMat];
        T = Rationalize[mat[[1 ;; 3, 4]]];
        det = Det[mat[[1 ;; 3, 1 ;; 3]]];
        {axes, phi} = Mat2EulerVector[mat[[1;;3, 1;;3]]];
        If[phi != 0,
             CS = If[det == 1, "C", "S"]; 
             label = \!\(\*TagBox[StyleBox[RowBox[{"\"\<\\!\\(\\*SubsuperscriptBox[\\(\>\"", "<>", "CS", "<>", "\"\<\\), \\(\>\"", "<>", RowBox[{"ToString", "[", RowBox[{"If", "[", RowBox[{RowBox[{"phi", "\\[NotEqual]", "0"}], ",", RowBox[{"2", " ", RowBox[{"Pi", "/", "phi"}]}], ",", "0"}], "]"}], "]"}], "<>", "\"\<\\), \\(\>\"", "<>", RowBox[{"StringReplace", "[", RowBox[{RowBox[{"StringJoin", "[", RowBox[{"Map", "[", RowBox[{"ToString", ",", "axes"}], "]"}], "]"}], ",", RowBox[{"\"\<-1\>\"", "\\[Rule]", "\"\<i\>\""}]}], "]"}], "<>", "\"\<\\)]\\)\>\""}], ShowSpecialCharacters->False,ShowStringCharacters->True, NumberMarks->True], FullForm]\) <> "{" <> StringJoin[Map[ToString, StringRiffle[T, ","]]] <> "}",
             label = If[det == 1, "E", "I"] <> "{" <> StringJoin[Map[ToString, StringRiffle[T, ","]]] <> "}"];
        Return[label]]
];

GetRegRep[grp_] := Module[{grpt, gt, el0, el, g, pos, i, j, k},
  grpt = Values[grp];
  g = Length[grpt];
  gt = Expand@Table[ModM4[grpt[[i]].Inverse[grpt[[j]]]], {i, 1, g}, {j, 1, g}];
  Return[Table[Table[If[gt[[i, j]] == grpt[[k]], 1, 0], {i, 1, g}, {j, 1, g}], {k, 1, g}]]
]

GetEleOrder[StrMat_] := Module[{ord, idm, mat, mat1},
  Which[AssociationQ[StrMat],
        GetEleOrder[#] &/@ Keys[StrMat],
        ListQ[StrMat]&&(!MatrixQ[StrMat]),
        GetEleOrder[#] &/@ StrMat,
        StringQ[StrMat]||MatrixQ[StrMat],
        mat = If[StringQ[StrMat], xyzStr2M4[StrMat], StrMat];
        ord = 1;
        idm=N[IdentityMatrix[Length[mat[[1]]]]];
        mat1 = mat;
        While[!idm == mat1, mat1 = Chop@GTimes[{mat, mat1}]; ord = ord + 1];
        Return[ord]]
]

GetIrep[grp_, p_, ct0_:{}] := Module[{i, j, k, g, lp, h, ct, ProjMat, RegRep, es, ev, vv, m, n}, 
  (*---Group order and dimension of irreducible representation---*)
  ct = If[Length[ct0]==0, GetCharacters[grp], ct0];
  g = Length[grp];
  lp = ct[[2]][[p]][[1]];
  h = Length[ct[[1]]];
  (*---regular representation and sort by order of element---*)
  RegRep = GetRegRep[grp];
  (*---Char Projection operator---*)
  ProjMat = Sum[{m, n} = First@Position[ct[[1]], Keys[grp][[i]]]; lp/g Conjugate[ct[[2]][[p, m]]] RegRep[[i]], {i, g}];
  (*---Find an eigenvector with non degenerate eigenvalue and apply other elements until lp linear independent vectors were found---*)
  Do[es = Eigensystem[RegRep[[i]]];
   Do[ev = ProjMat.es[[2, j]];
    If[Norm[ev] == 0, Continue[]];
    vv = Complement[Orthogonalize[Table[RegRep[[k]].ev, {k, 1, g}]], {Table[0, {i, 1, g}]}];
    If[Length[vv] == lp, Break[], None], {j, 1, g}], {i, g, g - 1, -1}];
  (*---calculate representation matrices---*)
  
  Return[Table[vv.RegRep[[i]].ConjugateTranspose[vv], {i, 1, g}]]
]

GetElePosition[grp_, ele_] := Module[{ind, type, axis, ang, shift},
  type = ele[[1]];
  ind = Which[type == "E" || type == "I", shift = "{" <> StringJoin[Map[ToString, StringRiffle[ele[[2]], ","]]] <> "}";
              Position[GetEleLabel[#] & /@ Keys[grp], \!\(\*TagBox[StyleBox[RowBox[{"type", "<>", "shift"}],ShowSpecialCharacters->False,ShowStringCharacters->True,NumberMarks->True],FullForm]\)], 
              type == "C" || type == "S", ang=ToString[ele[[2]]];axis=ele[[3]];shift="{"<>StringJoin[Map[ToString,StringRiffle[ele[[4]],","]]]<>"}";
              Position[GetEleLabel[#] & /@ Keys[grp], \!\(\*TagBox[StyleBox[RowBox[{"\"\<\\!\\(\\*SubsuperscriptBox[\\(\>\"", "<>", "type", "<>", "\"\<\\), \\(\>\"", "<>", "ang", "<>", "\"\<\\), \\(\>\"", "<>","axis", "<>", "\"\<\\)]\\)\>\"", "<>", "shift"}],ShowSpecialCharacters->False,ShowStringCharacters->True,NumberMarks->True],FullForm]\)], 
              True, Print["No such element"]; Abort[]];
  If[ind == {}, Print["No such element"]; Abort[]];
  Return[First@First@ind]
]

GroupQ[grp_] := Module[{},
  Return[Length@Complement[Union[Flatten[GrpMultiply[grp, grp]]], Keys[grp]] == 0]
]

GetSubGroups[grp_, ord_: {}] := Module[{GrpSub, g, div, sub, sg, test, j},
  GrpSub = Dispatch[Thread[Keys[grp] -> Values[grp]]];
  g = Length[grp];
  div = Divisors[g];
  div = Complement[If[ord == {}, Print["no index are specified, this may take very long!"]; Divisors[g], Intersection[Divisors[g], ord]], {1, g}];
  sg = {};
  Do[sub = Subsets[Complement[Keys[grp], {"x,y,z"}], {div[[i]] - 1}];
    Do[test = Association[# -> (# /. GrpSub) & /@ Union[sub[[j]], {"x,y,z"}]];
      AppendTo[sg, If[GroupQ[test], test, {}]], {j, Length@sub}], {i, Length@div}];
  Print[ToString[Length@Flatten[sg]] <> " sub group found."];
  Return[Flatten@sg]
]

SortByOrder[grp_] := SortBy[grp, (GetEleOrder[#]&)]

GetInvSubGroups[grp_, ind_:0] := Module[{g, GrpSub, classes, subs, invsub, test},
  GrpSub = Dispatch[Thread[Keys[grp] -> Values[grp]]];
  g = Length[grp];
  If[ind != 0, 
     If[Divisible[g, ind], Unevaluated[Sequence[]], Print["Index is not compatible!"];Abort[]], 
     Unevaluated[Sequence[]]];
  classes = Complement[GetClasses[grp], {{"x,y,z"}}];
  subs = If[ind ==0, Union[{{"x,y,z"}}, #] &/@ Subsets[classes], If[Length[Flatten[#]] == g/ind-1, Union[{{"x,y,z"}}, #], {}] &/@ Subsets[classes]];
  invsub = {};
  Do[test = Association[# -> (# /. GrpSub) & /@ Flatten[subs[[i]]]];
     If[GroupQ[test], If[Length[test] > 1 && Length[test] <= g, AppendTo[invsub, test], Null], Null], {i, Length@subs}];
  Return[SortBy[#, (GetEleOrder[#1]&)]&/@invsub]
]

Generators2Group[gen_] := Module[{grp0, grp1, GenMat, ln, i, j, ord1, ord0, grpout},
   GenMat = Values[gen];
   grp0 = grp0 = Union[ModM4[#] & /@ Flatten[Table[Fold[Dot, Table[#, i]], {i, GetEleOrder[#]}] & /@ GenMat, 1]];
   ord0 = Length[grp0];
   ord1 = 0;
   While[ord0 > ord1,
    grp1 = Union[Flatten[Table[ModM4[grp0[[i]].grp0[[j]]], {i, 1, Length[grp0]}, {j, 1, Length[grp0]}], 1]];
    ord1 = Length[grp0]; 
    ord0 = Length[grp1];
    grp0 = grp1];
   grpout = SortBy[grp1, GetEleOrder[#] &];
  Return[Association[M42xyzStr[#] -> Rationalize[#] & /@ grpout]]
]


GetGenerators[grp_] := Module[{GrpMat, g, ord, gen, t, co, ig, o1, i, whinp, subs, sl, log, final},
  GrpMat = Values[grp];
  g = Length[grp];
  (*---Search for generators---*)
  gen = {GrpMat[[g]]};
  t = ModM4[#] & /@ FoldList[Dot, Table[GrpMat[[g]], GetEleOrder[GrpMat[[g]]]]];
  co = Complement[GrpMat, t];
  ig = 1;
  While[Length[co] > 0,
   ig = ig + 1;
   o1 = Length[co];
   gen = Append[gen, co[[o1]]];
   t = Union[t, ModM4[gen[[ig]].#] & /@ t];
   co = Complement[GrpMat, t];];
  (*---minimize number of generators---*)
  ord = Length[gen];
  subs = Subsets[gen, {1, ord}];
  sl = Length[subs];
  log = True;
  i = 1;
  While[log,
   t = Generators2Group[Association[M42xyzStr[#] -> # & /@ subs[[i]]]];
   If[Complement[GrpMat, Values@t] == {},
    final = Association[M42xyzStr[#] -> # & /@ subs[[i]]];
    log = False, If[i < sl, i = i + 1, log = False]]];
  Return[final]
]

GetSpgCosetRepr[grp0_, invsubgrp_] := Module[{ind, g, comp, qel, reordgroup, Tm, groupfound},
  ind = Length[grp0]/Length[invsubgrp];
  g = Length[invsubgrp];
  comp = Complement[Keys[grp0], Keys[invsubgrp]];
  reordgroup = {};
  groupfound = False;
  Tm = 0;
  While[Not[groupfound] && Tm < g,
   Tm = Tm + 1;
   qel = comp[[Tm]];
   reordgroup =
     Which[
       ind == 2,
       Flatten[{Keys[invsubgrp], GTimes[{qel, Keys[invsubgrp]}]}, 1],
       ind == 3,
       Flatten[{Keys[invsubgrp], GTimes[{qel, Keys[invsubgrp]}], GTimes[{qel, qel, Keys[invsubgrp]}]}, 1]
       ];
     groupfound = Length[Complement[reordgroup, Keys[grp0]]] === 0;
     ];
  If[Tm > g, Print["Error: Could not identify coset representative!"]; Abort[]];
  Return[qel];
] 

SolidSphericalHarmonicY[l_?IntegerQ, m_?IntegerQ, x1_, x2_, x3_, coord_: "Cartesian"] := Module[{rr, \[Theta]\[Theta], \[Phi]\[Phi], xx, yy, zz},
  Which[coord == "Cartesian", 
                 FullSimplify@Evaluate[TransformedField["Spherical" -> "Cartesian", rr^l SphericalHarmonicY[l, m, \[Theta]\[Theta], \[Phi]\[Phi]], {rr, \[Theta]\[Theta], \[Phi]\[Phi]} -> {xx, yy, zz}]] /. {xx -> x1, yy -> x2, zz -> x3},
        coord == "Polar", 
                 FullSimplify@Evaluate[rr^l SphericalHarmonicY[l, m, \[Theta]\[Theta], \[Phi]\[Phi]]] /. {rr -> x1, \[Theta]\[Theta] -> x2, \[Phi]\[Phi] -> x3}, 
        True, 
        Print["Wrong coordinate"]]
]

SolidTesseralHarmonicY[l_?IntegerQ, m_?IntegerQ, x1_, x2_, x3_, coord_: "Cartesian"] := Module[{}, 
  Simplify@Which[
                 m > 0, 
                 Sqrt[2]/2 (-1)^m (SolidSphericalHarmonicY[l, m, x1, x2, x3, coord] + (-1)^m SolidSphericalHarmonicY[l, -m, x1, x2, x3, coord]), 
                 m < 0, 
                 Sqrt[2]/(2 I ) (-1)^m (SolidSphericalHarmonicY[l, Abs[m], x1, x2, x3, coord] - (-1)^m SolidSphericalHarmonicY[l, -Abs[m], x1, x2, x3, coord]), 
                 m == 0, 
                 SolidSphericalHarmonicY[l, 0, x1, x2, x3, coord]
                ]
]

Dl2El[m_, mm_] := Which[m > 0, (-1)^m 1/Sqrt[2] (KroneckerDelta[m, mm] + KroneckerDelta[-m, mm]), m == 0, KroneckerDelta[m, mm], m < 0, (-1)^m 1/(I Sqrt[2]) (KroneckerDelta[m, -mm] - KroneckerDelta[m, mm])]

GetAngularMomentumRep[latt_, mat_, l_, Harmonic_: "Tesseral"] := Module[{\[Alpha], \[Beta], \[Gamma], \[Epsilon], m1, m2, Dlmn, T},
  {{\[Alpha], \[Beta], \[Gamma]}, \[Epsilon]} = Mat2EulerAngles[latt, mat];
  Dlmn = Simplify[\[Epsilon]^l Table[WignerD[{l, m1, m2}, \[Alpha], \[Beta], \[Gamma]], {m1, -l, l}, {m2, -l, l}]];
  Which[Harmonic == "Tesseral",
   T = Table[Dl2El[m, mm], {mm, -l, l}, {m, -l, l}];
   Dlmn = Inverse[T].Dlmn.T,
   Harmonic == "Spherical",
   Dlmn];
  Return[Dlmn]
]

GrpMultiply[lgrp_, rgrp_] := Module[{},
  Which[
   AssociationQ[lgrp] && AssociationQ[rgrp],
     GrpMultiply[#1, #2] & @@@ Tuples[{Values@lgrp, Values@rgrp}],
   AssociationQ[lgrp] && MatrixQ[rgrp],
      GrpMultiply[#, rgrp] &/@ Values@lgrp,
   AssociationQ[lgrp] && StringQ[rgrp],
      GrpMultiply[#, xyzStr2M4@rgrp] &/@ Values@lgrp,
   MatrixQ[lgrp] && AssociationQ[rgrp],
      GrpMultiply[lgrp, #] &/@ Values@rgrp,
   StringQ[lgrp] && AssociationQ[rgrp],
      GrpMultiply[xyzStr2M4@lgrp, #] &/@ Values@rgrp,
   AssociationQ[lgrp] && Length[Dimensions[rgrp]] == 3,
      GrpMultiply[#1, #2] & @@@ Tuples[{Values@lgrp, rgrp}],
   AssociationQ[lgrp] && Length[Dimensions[rgrp]] == 1,
      GrpMultiply[#1, xyzStr2M4@#2] & @@@ Tuples[{Values@lgrp, rgrp}],
   Length[Dimensions[lgrp]] == 3 && AssociationQ[rgrp],
      GrpMultiply[#1, #2] & @@@ Tuples[{lgrp, Values@rgrp}],
   Length[Dimensions[lgrp]] == 1 && AssociationQ[rgrp],
       GrpMultiply[xyzStr2M4@#1, #2] & @@@ Tuples[{lgrp, Values@rgrp}],
   Length[Dimensions[lgrp]] == 3 && Length[Dimensions[rgrp]] == 3,
     GrpMultiply[#1, #2] & @@@ Tuples[{lgrp, rgrp}],
   Length[Dimensions[lgrp]] == 3 && MatrixQ[rgrp],
     GrpMultiply[#, rgrp] & /@ lgrp,
   MatrixQ[lgrp] && Length[Dimensions[rgrp]] == 3,
     GrpMultiply[lgrp, #] & /@ rgrp,
   Length[Dimensions[lgrp]] == 3 && Length[Dimensions[rgrp]] == 1,
     GrpMultiply[#1, #2] & @@@ Tuples[{lgrp, (xyzStr2M4[#] &/@ rgrp)}],
   Length[Dimensions[lgrp]] == 1 && Length[Dimensions[rgrp]] == 3,
     GrpMultiply[#1, #2] & @@@ Tuples[{(xyzStr2M4[#] &/@ lgrp), rgrp}],
   Length[Dimensions[lgrp]] == 3 && Length[Dimensions[rgrp]] == 0,
     GrpMultiply[#, xyzStr2M4[rgrp]] &/@ lgrp,
   Length[Dimensions[lgrp]] == 0 && Length[Dimensions[rgrp]] == 3,
     GrpMultiply[xyzStr2M4[lgrp], #] &/@ rgrp,
   MatrixQ[lgrp] && MatrixQ[rgrp],
     M42xyzStr[ModM4[lgrp.rgrp]],
   MatrixQ[lgrp] && Length[Dimensions[rgrp]] == 1,
     GrpMultiply[lgrp, xyzStr2M4[#]] &/@ rgrp,
   Length[Dimensions[lgrp]] == 1 && MatrixQ[rgrp],
     GrpMultiply[xyzStr2M4[#], rgrp] &/@ lgrp,
   MatrixQ[lgrp] && Length[Dimensions[rgrp]] == 0,
     GrpMultiply[lgrp, xyzStr2M4[rgrp]],
   Length[Dimensions[lgrp]] == 0 && MatrixQ[rgrp],
     GrpMultiply[xyzStr2M4[lgrp], rgrp],
   Length[Dimensions[lgrp]] == 1 && Length[Dimensions[rgrp]] == 1,
     GrpMultiply[#1, #2] &@@@ Tuples[{lgrp, rgrp}],
   Length[Dimensions[lgrp]] == 1 && Length[Dimensions[rgrp]] == 0,
     GrpMultiply[#, rgrp] &/@ lgrp,
   Length[Dimensions[lgrp]] == 0 && Length[Dimensions[rgrp]] == 1,
     GrpMultiply[lgrp, #] &/@ rgrp,
   Length[Dimensions[lgrp]] == 0 && Length[Dimensions[rgrp]] == 0,
     GrpMultiply[xyzStr2M4[lgrp], xyzStr2M4[rgrp]]
  ]
]

GrpxV[grp_, v_, OptionsPattern[{"k"->False}]] := Module[{},
  Which[AssociationQ[grp]&&Length[Dimensions@v]==1, GrpxV[#,v] &/@ Keys[grp],
        AssociationQ[grp]&&Length[Dimensions@v]==2, GrpxV[#1,#2] &@@@ Tuples[{Keys[grp],v}],
        StringQ[grp]&&(Length[Dimensions@v]==1), GrpxV[xyzStr2M4[grp],v],
        StringQ[grp]&&(Length[Dimensions@v]==2), GrpxV[xyzStr2M4[grp],#] &/@ v,
        (ListQ[grp])&&(!MatrixQ[grp])&&(Length[Dimensions@v]==1), GrpxV[#,v] &/@ grp,
        (ListQ[grp])&&(!MatrixQ[grp])&&(Length[Dimensions@v]==2), GrpxV[#1,#2] &@@@ Tuples[{grp,v}],
        MatrixQ[grp]&&(Length[Dimensions@v]==2), GrpxV[grp,#] &/@ v,
        MatrixQ[grp]&&(Length[Dimensions@v]==1), If[OptionValue["k"], grp[[1;;3,1;;3]].v, (grp.Append[v, 1])[[1;;3]]]
  ]
]
        

GTimes[list_] := Module[{},
  Fold[GrpMultiply, list]
]

GInverse[ele_] := Module[{},
  Which[AssociationQ[ele], GInverse[Values[ele]],
        ListQ[ele]&&!MatrixQ[ele], GInverse[#] &/@ ele,
        StringQ[ele], GInverse[xyzStr2M4[ele]],
        MatrixQ[ele], ModM4[Inverse[ele]]
  ]
]

GetLeftCosets[grp0_, invsub_] := Module[{i, j},
  Union[M42xyzStr[SortByOrder[#]] & /@ Table[GTimes[{grp0[[i]], invsub[[j]]}], {i, Length[grp0]}, {j, Length[invsub]}]]
]
(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
