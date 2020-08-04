BeginPackage["LINVARIANT`GroupTheory`", {"LINVARIANT`SphericalHarmonics`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
CifImportSpg               ::usage "CifImportSpg[file]"
GrpQ                       ::usage "GrpQ[grp]"
GetGrpK                    ::usage "GetGrpK[grp0, vec0, lattice0]"
GetStarK                   ::usage "GetStarK[grp0, k]"
GetSiteSymmetry            ::usage "GetSiteSymmetry[grp0, vec0]"
CifImportOperators         ::usage "CifImportOperators[file]"
GetGenerators              ::usage "GetGenerators[grp]"
Generators2Grp             ::usage "Generators2Grp[gen]"
xyzStr2Grp                 ::usage "xyzStr2Grp[keys]"
xyzStr2Expression          ::usage "xyzStr2Expression[xyzStrData]"
xyz2Rot                    ::usage "xyz2Rot[expr]"
xyzStr2M4                  ::usage "xyzStr2M4[xyzStrData]"
xyzStr2TRot                ::usage "xyzStr2TRot[xyzStrData]"
RotTran2M4                 ::usage "RotTran2M4[rot, tran]"
M42xyzStr                  ::usage "M42xyzStr[m4]"
M42xyzStrEle               ::usage "M42xyzStrEle[m4]"
GetClasses                 ::usage "GetClasses[GrpMat]"
GetPGCharacters            ::usage "GetPGCharacters[GrpMat]"
ExpandSimplify             ::usage "ExpandSimplify[in]"
GetEleLabel                ::usage "GetEleLabel[StrMat]"
GetRegRep                  ::usage "GetRegularRepresentation[grp]"
GetEleOrder                ::usage "GetEleOrder[mat]"
GetPGIreps                 ::usage "GetPGIreps[grp, ct, p]"
ModM4                      ::usage "ModM4[m1, m2]"
GrpMultiply                ::usage "GrpMultiply[lgrp, rgrp]"
GrpV                       ::usage "GrpV[grp, v]"
GrpO                       ::usage "GrpO[grp, O, vars]"
SortByOrder                ::usage "SortByOrder[grp]"
SortGrp                    ::usage "SortGrp[grp]"
GetElePosition             ::usage "GetElePosition[grp, ele]"
GetSubGrp                  ::usage "GetSubGrp[grp, ord]"
GetInvSubGrp               ::usage "GTInvSubGrp[grp]"
GetLeftCosetRepr           ::usage "GetLeftCosetRepr[grp0, invsubgrp]"
GetAngularMomentumRep      ::usage "GetAngularMomentumRep[latt, mat, l, Harmonic]"
GetAngularMomentumChars    ::usage "GetAngularMomentumChars[grp, j]"
GInverse                   ::usage "GInverse[ele]"
GTimes                     ::usage "GTimes[list]"
GPower                     ::usage "GPower[key, n]"
GetLeftCosets              ::usage "GetLeftCosets[grp0, invsub]"
GetGrpKIreps               ::usage "GetGrpKIreps[spg, kvec]"
GetSpgIreps                ::usage "GetSpgIreps[spg, kvec]"
GetProjOperator            ::usage "GetProjOperator[grp, ireps]"
ProjectOnOperator          ::usage "ProjectOnOperator[grp, ireps, O, vars]"
ProjectOnBasis             ::usage "ProjectOnBasis[grp, ireps, OpMat, basis]"
GetOrbitInducedIreps       ::usage "GetOrbitInducedIreps[kvec, grp, grpH, q, OrbitIreps, perm]"
GetUnitaryTransMatrix      ::usage "GetUnitaryTransMatrix[Amat, Bmat, Cmat, n]"
GetSpgCT                   ::usage "GetSpgCT[spg0, kvec, irreps]"
DecomposeIreps             ::usage "DecomposeIreps[grp, ct, character]"
GetCGCoefficients          ::usage "GetCGCoefficients[grp, ireps, p1p2]"
GetCGjmjm2jm               ::usage "GetCGjmjm2jm[jm1, jm2]"
GetSuperCellGrp            ::usage "GetSuperCellGrp[t]"
GetEulerVector             ::usage "GetEulerVector[ele]"
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
  Which[StringQ[m0], ModM4[xyzStr2M4[m0]],
        AssociationQ[m0], ModM4[#] &/@ Values[m0],
        Length[Dimensions[m0]]==1, ModM4[#] &/@ m0,
        Length[Dimensions[m0]]==3, ModM4[#] &/@ m0,
        MatrixQ[m0],
        m = m0;
        If[Length[m0] == 4, m[[1 ;; 3, 4]] = Mod[m0[[1 ;; 3, 4]], 1], Unevaluated[Sequence[]]];
        Return[M42xyzStr@m]
       ]
]

ExpandSimplify[in_] := Expand[Simplify[in]];

GetWyckoffImages[grp_, pos_] := Module[{xyzStrData, xyzExpData, WyckoffImages},
  xyzStrData = Keys[grp];
  xyzExpData = Table[ToExpression["{" <> xyzStrData[[i]] <> "}"], {i, Length[xyzStrData]}];
  WyckoffImages = Flatten[Table[Sort[Union[{Mod[#,1], pos[[i, 2]]} & /@ (xyzExpData /. Thread[ToExpression[{"x", "y", "z"}] -> pos[[i, 1]]])], Total[#1[[1]]] > Total[#2[[1]]] &], {i, Length[pos]}], 1];
  Return[WyckoffImages];
]

GetClasses[grp0_] := Module[{i, j, g, bab, cl},
  g = Length[grp0];
  bab = ParallelTable[M42xyzStr[ModM4@GTimes[{grp0[[i]],grp0[[j]],GInverse[grp0[[i]]]}]], {j, 1, g}, {i, 1, g}, DistributedContexts -> {"LINVARIANT`GroupTheory`Private`"}];

  cl = Flatten[Values[KeySortBy[GroupBy[#, First@Union[GetEulerVector[#]\[Transpose][[2]]] &], Minus]]
               & /@ Values[KeySortBy[GroupBy[#, First@Union[GetEulerVector[#]\[Transpose][[3]]] &], Minus]]
               & /@ Values[KeySort[GroupBy[Union[Map[Union, bab]], Length[#] &]]], 3];
  cl = Flatten[(Values[KeySortBy[GroupBy[#, GetEulerVector[#][[1]] &], Minus]]
               & /@ Values[KeySortBy[GroupBy[#, GetEulerVector[#][[2]] &], Minus]]), 2] & /@ cl;

  cl = Prepend[DeleteCases[cl, {"x,y,z"}], {"x,y,z"}];

  Return[cl]
]

GetPGCharacters[grp0_, OptionsPattern[{"print"->True}]] := Module[{m, n, i, j, k, h, g, nc, Classes, H0Mat, H1Mat, HMat, Hi, hh, Vt, V, Vi, cc, vec, IrepDim, CharacterTable, CharacterTableSorted, ct, IrepLabel, ClassesLabel},
  Classes = GetClasses[grp0];
  h = Length[Classes];
  g = Length[grp0];
  nc = Map[Length, Classes];

  If[g > 1,
    (*---calculate the class constant matrices---*)
    H0Mat = Table[Table[Flatten[Table[GTimes[{Classes[[j, m]],Classes[[i, n]]}], {m, 1, nc[[j]]}, {n, 1, nc[[i]]}], 1], {i, 1, h}], {j, 1, h}];
    H1Mat = Table[Table[Length[Position[H0Mat[[i, j]], Classes[[k, 1]]]], {i, 1, h}, {j, 1, h}], {k, 1, h}];
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
    IrepDim = Sqrt[Table[g/(cc[[k]]/nc).Conjugate[cc[[k]]], {k, 1, h}]];

    (*---calculate the character table---*)
    CharacterTable = FullSimplify@Table[Table[cc[[i, j]]*IrepDim[[i]]/nc[[j]], {j, 1, h}], {i, 1, h}];

    (*---sortby---*)
    CharacterTableSorted = SortBy[CharacterTable, First];
    ct = Prepend[Delete[CharacterTableSorted, First@Flatten@Position[CharacterTableSorted, Table[1, {i, 1, Length[CharacterTable]}]]], Table[1, {i, 1, Length[CharacterTable]}]];,
    ct = {{1}};
  ];


  IrepLabel = Table["\!\(\*SubscriptBox[\(\[CapitalGamma]\), \("<>ToString[y]<>"\)]\)", {y, 1, h}];
  If[OptionValue["print"],
  Print[TableForm[ct, TableHeadings -> {IrepLabel, ToString[Length[#]]<>GetEleLabel[First[#]] & /@ Classes}, TableAlignments -> {Right, Center}]]];
  Return[ct]
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
  RT = Join[Join[xyz2Rot[rot]\[Transpose], {tran}]\[Transpose], {{0, 0, 0, 1}}];
  Return[Rationalize[RT]]
]

xyzStr2M4[xyzStrData_] := Module[{},
  Which[StringQ[xyzStrData], xyzStr2M4Ele[xyzStrData], ListQ[xyzStrData], xyzStr2M4[#] &/@ xyzStrData, True, Print["Wrong Input type!"];Abort[]]
]

xyzStr2TRotEle[xyzStrData_] := Module[{rot, tran},
  {rot, tran} = xyzStr2Expression[xyzStrData];
  Return[{xyz2Rot[rot], tran}]
]

xyzStr2TRot[xyzStrData_] := Module[{},
  Which[StringQ[xyzStrData], xyzStr2TRotEle[xyzStrData], ListQ[xyzStrData], xyzStr2TRot[#] &/@ xyzStrData, True, Print["Wrong Input type!"];Abort[]]
]

RotTran2M4[rot_, tran_] := Module[{RT},
  RT = Join[Join[rot\[Transpose], {tran}]\[Transpose], {{0., 0., 0., 1.}}];
  Return[Rationalize[RT, 10^-12]]
]

M42xyzStrEle[m4_] := Module[{i, rot, tran},
  rot = ToString[#] & /@ (Rationalize[m4.{ToExpression["x"], ToExpression["y"], ToExpression["z"], 1}] - (Rationalize[m4.{ToExpression["x"], ToExpression["y"], ToExpression["z"], 1}] /. {ToExpression["x"] -> 0, ToExpression["y"] -> 0, ToExpression["z"] -> 0}));
  tran = (Rationalize[m4.{ToExpression["x"], ToExpression["y"], ToExpression["z"], 1}] /. {ToExpression["x"] -> 0, ToExpression["y"] -> 0, ToExpression["z"] -> 0});
  tran = Which[#1 == 0, "",
               #1 != 0 && #2 == 1, If[Sign[#1] > 0, "+", Unevaluated[Sequence[]]] <> ToString[#1],
               True, If[Sign[#1] > 0, "+", Unevaluated[Sequence[]]] <> ToString[#1] <> "/" <> ToString[#2]]
         & @@@ (Through[{Numerator, Denominator}[#]] & /@ tran)[[1 ;; 3]];
  Return[StringDelete[StringJoin@Riffle[Table[rot[[i]] <> tran[[i]], {i, 3}],","], " "]]
]

M42xyzStr[m4_] := Module[{},
  If[MatrixQ[m4], M42xyzStrEle[m4], M42xyzStr[#] &/@ m4]
]

xyz2Rot[expr_] := Module[{m},
  m = Coefficient[#, {ToExpression["x"], ToExpression["y"], ToExpression["z"]}] &/@ expr;
  If[Abs[Det[m]] != 1, Print["xyz2Rot gives wrong determinant "<>ToString[Det[m]]]; Abort[], Unevaluate[Sequence[]]];
  Return[m]
]

CifImportSpg[file_] := Module[{xyzStrData, ele, RT},
  xyzStrData = CifImportOperators[file];
  ele = ModM4[xyzStr2M4[#]] &/@ xyzStrData;
  Return[SortGrp[ele]]
]

GetGrpK[grp0_, vec0_] := Module[{grp, g, trvec, eq, sol, posmap, grpk, c, vec},

  trvec = GrpV[grp0, vec0, "k" -> True];
  eq = Thread[# == vec0 + Table[c[i],{i,3}]] &/@ trvec;
  sol = Rationalize[Round[Table[c[j], {j, 3}] /. Flatten[Solve[#, Table[c[j], {j, 3}]]] &/@ eq, 10^-3]];
  posmap = AllTrue[#, IntegerQ] &/@ sol;
  grpk = xyzStr2Grp@Extract[Keys[grp0], Position[posmap,True]];
  Return[grpk]
]

GetStarK[grp0_, k_] := Module[{grpk, LeftCoset, LeftRep},
  grpk = GetGrpK[grp0, k];
  LeftCoset = GetLeftCosets[grp0, grpk];
  LeftRep = LeftCoset[[;; , 1]];
  Return[Mod[GrpV[LeftRep, k, "k"->True],1]]
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

GetEleLabel[StrMat_] := Module[{det, mat, axes, phi, T, CS, label},
  Which[ListQ[StrMat]&&(!MatrixQ[StrMat]),
        GetEleLabel[#] &/@ StrMat,
        AssociationQ[StrMat],
        GetEleLabel[#] &/@ Keys[StrMat],
        StringQ[StrMat]||MatrixQ[StrMat],
        mat = If[StringQ[StrMat], xyzStr2M4[StrMat], StrMat];
        T = Rationalize[mat[[1 ;; 3, 4]]];
        {axes, phi, det} = Mat2EulerVector[mat[[1;;3, 1;;3]]];
        If[phi != 0,
             CS = If[det == 1, "C", "S"]; 
             label = \!\(\*TagBox[StyleBox[RowBox[{"\"\<\\!\\(\\*SubsuperscriptBox[\\(\>\"", "<>", "CS", "<>", "\"\<\\), \\(\>\"", "<>", RowBox[{"ToString", "[", RowBox[{"If", "[", RowBox[{RowBox[{"phi", "\\[NotEqual]", "0"}], ",", RowBox[{"2", " ", RowBox[{"Pi", "/", "phi"}]}], ",", "0"}], "]"}], "]"}], "<>", "\"\<\\), \\(\>\"", "<>", RowBox[{"StringReplace", "[", RowBox[{RowBox[{"StringJoin", "[", RowBox[{"Map", "[", RowBox[{"ToString", ",", "axes"}], "]"}], "]"}], ",", RowBox[{"\"\<-1\>\"", "\\[Rule]", "\"\<i\>\""}]}], "]"}], "<>", "\"\<\\)]\\)\>\""}], ShowSpecialCharacters->False,ShowStringCharacters->True, NumberMarks->True], FullForm]\) <> "{" <> StringJoin[Map[ToString, StringRiffle[T, ","]]] <> "}",
             label = If[det == 1, "E", "I"] <> "{" <> StringJoin[Map[ToString, StringRiffle[T, ","]]] <> "}"];
        Return[label]]
];

GetRegRep[grp_] := Module[{gt, g, i, j, k},
  g = Length[grp];
  gt = Expand@Table[GTimes[{grp[[i]], GInverse[grp[[j]]]}], {i, 1, g}, {j, 1, g}];
  Return[Table[Table[If[gt[[i, j]] == Keys[grp][[k]], 1, 0], {i, 1, g}, {j, 1, g}], {k, 1, g}]]
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
        While[!idm == mat1, mat1 = Chop@ModM4@GTimes[{mat, mat1}]; ord = ord + 1];
        Return[ord]]
]

GetPGIreps[grp_] := Module[{i, j, k, g, p, lp, h, ct, classes, ireps, ProjMat, RegRep, es, ev, vv, m, n}, 
  (*---Group order and dimension of irreducible representation---*)
  classes = GetClasses[grp];
  ct = GetPGCharacters[grp];
  g = Length[grp];
  h = Length[ct];
  RegRep = GetRegRep[grp];
  ireps = Table[lp = ct[[p]][[1]];
                ProjMat = Sum[{m, n} = First@Position[classes, Keys[grp][[i]]]; 
                              lp/g Conjugate[ct[[p, m]]] RegRep[[i]], {i, g}];
                (*Find an eigenvector with non degenerate eigenvalue and apply other elements until lp linear independent vectors were found*)
                Do[es = Eigensystem[RegRep[[i]]];
                 Do[ev = ProjMat.es[[2, j]];
                  If[Norm[ev] == 0, Continue[]];
                  vv = Complement[Orthogonalize[Table[RegRep[[k]].ev, {k, 1, g}]], {Table[0, {i, 1, g}]}];
                  If[Length[vv] == lp, Break[], None], {j, 1, g}], {i, g, g - 1, -1}];
                (*---calculate representation matrices---*)
                Table[vv.RegRep[[i]].ConjugateTranspose[vv], {i, 1, g}], {p, 1, h}];
  
  Return[ireps]
]

GetElePosition[grp_, ele_] := Module[{ind, type, axis, ang, shift},
  Which[StringQ[ele],
        type = StringCases[ele, RegularExpression["\([A-Z]\)"]];
        ind = If[type != {},
                 Position[GetEleLabel[grp], ele], 
                 Position[Keys[grp], ele]],
        Length[Dimensions[ele]] == 2,
        ind = GetElePosition[grp,M42xyzStr[ele]]
        ];
  If[ind == {}, Print["No such element"]; Abort[]];
  Return[First@First@ind]
]

GrpQ[grp_] := Module[{},
  Return[Length@Complement[Union[ModM4@Flatten[GTimes[{grp, grp}]]], Which[AssociationQ[grp], Keys[grp], Length@Dimensions[grp] == 1, grp, Length@Dimensions[grp] > 1, M42xyzStr[#] &/@ grp]] == 0]
]

GetSubGrp[grp_, ord_: {}] := Module[{GrpSub, g, div, sub, sg, test, j},
  g = Length[grp];
  div = Divisors[g];
  div = Complement[If[ord == {}, Print["no index are specified, this may take very long!"]; Divisors[g], Intersection[Divisors[g], ord]], {1, g}];
  sg = {};
  Do[sub = Subsets[Complement[Keys[grp], {"x,y,z"}], {div[[i]] - 1}];
    Do[test = Union[sub[[j]], {"x,y,z"}];
      AppendTo[sg, If[GrpQ[test], xyzStr2Grp@test, {}]], {j, Length@sub}], {i, Length@div}];
  Print[ToString[Length@Flatten[sg]] <> " sub group found."];
  Return[Flatten@sg]
]

xyzStr2Grp[keys_, OptionsPattern[{"fast"->False}]] := Module[{xyz},
  xyz = Prepend[DeleteCases[keys, "x,y,z"], "x,y,z"];
  If[OptionValue["fast"],
     Return[Association[Thread[xyz->xyzStr2M4[xyz]]]],
     If[GrpQ[xyz],
     Return[Association[Thread[xyz->xyzStr2M4[xyz]]]],
     Print["keys don't form a group!"];Abort[]]
  ]
]

SortByOrder[grp_] := SortBy[grp, (GetEleOrder[#]&)]

SortGrp[grp_] := Module[{cl,grpstd},
  cl = GetClasses[grp];
  cl = Flatten[Values[KeySortBy[GroupBy[#, First@Union[GetEulerVector[#]\[Transpose][[2]]] &], Minus]] 
               & /@ Values[KeySortBy[GroupBy[#, First@Union[GetEulerVector[#]\[Transpose][[3]]] &], Minus]] 
               & /@ Values[KeySort[GroupBy[cl, Length[#] &]]], 3];
  cl = Flatten[(Values[KeySortBy[GroupBy[#, GetEulerVector[#][[1]] &], Minus]] 
               & /@ Values[KeySortBy[GroupBy[#, GetEulerVector[#][[2]] &], Minus]]), 2] & /@ cl;
  cl = Flatten[cl];

  grpstd = Association[Thread[cl->xyzStr2M4[cl]]];
  Return[grpstd]
];

GetInvSubGrp[grp_, ind_:0, OptionsPattern[{"all"->True}]] := Module[{g, GrpListxyzStr, classes, subs, invsub, test, invfound, ipos},
  g = Length[grp];
  If[ind != 0, 
     If[Divisible[g, ind], Unevaluated[Sequence[]], Print["Index is not compatible!"];Abort[]], 
     Unevaluated[Sequence[]]];
  classes = DeleteCases[GetClasses[grp], {"x,y,z"}];

  subs = If[ind ==0, Union[{{"x,y,z"}}, #] &/@ Subsets[classes], If[Length[Flatten[#]] == g/ind-1, Union[{{"x,y,z"}}, #], {}] &/@ Subsets[classes]];
  invsub = {};
  ipos = 0;
  invfound = False;
  If[OptionValue["all"],
     Do[test = Flatten[subs[[i]]];
        If[GrpQ[test], 
           If[Length[test] >= 1 && Length[test] <= g, 
              AppendTo[invsub, xyzStr2Grp@test],
              Null], 
           Null], {i, Length@subs}],
     While[Not[invfound]&&ipos<Length[subs]&&Length[subs]>0, 
           ipos = ipos + 1;
           test = Flatten[subs[[ipos]]];
           If[GrpQ[test], If[Length[test] >= 1 && Length[test] <= g, invfound = True;invsub={xyzStr2Grp@test}]];]
  ];

  Return[invsub]
]

Generators2Grp[gen_] := Module[{grp0, grp1, ln, i, j, ord1, ord0, grpout},
   grp0 = Union@ModM4@Flatten[Table[GTimes@Table[#, i], {i, GetEleOrder[#]}] & /@ gen];
   ord0 = Length[grp0];
   ord1 = 0;
   While[ord0 > ord1,
    grp1 = Union@Flatten@Table[ModM4@GTimes[{grp0[[i]],grp0[[j]]}], {i, 1, Length[grp0]}, {j, 1, Length[grp0]}];
    ord1 = Length[grp0]; 
    ord0 = Length[grp1];
    grp0 = grp1];
   grpout = SortGrp@grp1;
  Return[grpout]
]


GetGenerators[grp_] := Module[{g, ord, gen, t, co, ig, o1, i, whinp, subs, sl, log, final},
  g = Length[grp];
  (*---Search for generators---*)
  gen = {Keys[grp][[g]]};
  t = ModM4[#] & /@ FoldList[Dot, Table[grp[[g]], GetEleOrder[grp[[g]]]]];
  co = Complement[Keys[grp], t];
  ig = 1;
  While[Length[co] > 0,
   ig = ig + 1;
   o1 = Length[co];
   gen = Append[gen, co[[o1]]];
   t = Union[t, ModM4[GTimes[{gen[[ig]],#}]] & /@ t];
   co = Complement[Keys[grp], t];];
  (*---minimize number of generators---*)
  ord = Length[gen];
  subs = Subsets[gen, {1, ord}];
  sl = Length[subs];
  log = True;
  i = 1;
  While[log,
   t = Generators2Grp[subs[[i]]];
   If[Complement[Keys@grp, Keys@t] == {},
      final = subs[[i]];
      log = False, 
      If[i < sl, i = i + 1, log = False]]];
  Return[final]
]

GetLeftCosetRepr[grp0_, invsubgrp_] := Module[{ind, g, comp, qel, reordgroup, Tm, groupfound},
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

Dl2El[m_, mm_] := Which[m > 0, (-1)^m 1/Sqrt[2] (KroneckerDelta[m, mm] + KroneckerDelta[-m, mm]), m == 0, KroneckerDelta[m, mm], m < 0, (-1)^m 1/(I Sqrt[2]) (KroneckerDelta[m, -mm] - KroneckerDelta[m, mm])]

GetAngularMomentumRep[latt_, grp_, j_, Harmonic_: "Tesseral"] := Module[{mat, \[Alpha], \[Beta], \[Gamma], \[Epsilon], m1, m2, Dlmn, T},
  Which[AssociationQ[grp],
        GetAngularMomentumRep[latt, #, j, Harmonic] &/@ Keys[grp],
        StringQ[grp],
        mat = First@xyzStr2TRot[grp];
        GetAngularMomentumRep[latt, mat, j, Harmonic],
        Length[Dimensions[grp]]==1,
        GetAngularMomentumRep[latt, #, j, Harmonic] &/@ grp,
        Length[Dimensions[grp]]==3,
        GetAngularMomentumRep[latt, #, j, Harmonic] &/@ grp,
        Length[Dimensions[grp]]==2,
        mat = grp[[1;;3,1;;3]];
        {{\[Alpha], \[Beta], \[Gamma]}, \[Epsilon]} = Mat2EulerAngles[latt, mat];
        Dlmn = Simplify[\[Epsilon]^j Table[WignerD[{j, m1, m2}, \[Alpha], \[Beta], \[Gamma]], {m1, -j, j}, {m2, -j, j}]];
        Which[Harmonic == "Tesseral",
              T = Table[Dl2El[m, mm], {mm, -j, j}, {m, -j, j}];
              Dlmn = Inverse[T].Dlmn.T,
              Harmonic == "Spherical",
              Dlmn
              ];
        Return[Rationalize@Chop@Dlmn]
    ]
]

GetAngularMomentumChars[grp_, j_] := Module[{char, jm, mat, axis, \[Phi], \[Epsilon]},
  Which[AssociationQ[grp],
        GetAngularMomentumChars[#, j] &/@ Keys[grp],
        StringQ[grp],
        mat = First@xyzStr2TRot[grp];
        GetAngularMomentumChars[mat, j],
        Length[Dimensions[grp]]==1,
        GetAngularMomentumChars[#, j] &/@ grp,
        Length[Dimensions[grp]]==3,
        GetAngularMomentumChars[#, j] &/@ grp,
        Length[Dimensions[grp]]==2,
        mat = grp[[1;;3,1;;3]];
        {axis, \[Phi], \[Epsilon]} = Mat2EulerVector[mat];
        (*char = \[Epsilon]^j*Sin[(j+1/2)*\[Phi]]/Sin[\[Phi]/2];*)
        char = Simplify[\[Epsilon]^j*Sum[Exp[-I jm \[Phi]], {jm, -j, j}]];
        Return[char]
    ]

]

GrpMultiply[lgrp_, rgrp_] := Module[{i,j},
  Which[
   AssociationQ[lgrp] && AssociationQ[rgrp],
      Table[GrpMultiply[Values[lgrp][[i]], Values[rgrp][[j]]], {i, Length@Values@lgrp}, {j, Length@Values@rgrp}],
   AssociationQ[lgrp] && MatrixQ[rgrp],
      GrpMultiply[#, rgrp] &/@ Values@lgrp,
   AssociationQ[lgrp] && StringQ[rgrp],
      GrpMultiply[#, xyzStr2M4@rgrp] &/@ Values@lgrp,
   MatrixQ[lgrp] && AssociationQ[rgrp],
      GrpMultiply[lgrp, #] &/@ Values@rgrp,
   StringQ[lgrp] && AssociationQ[rgrp],
      GrpMultiply[xyzStr2M4@lgrp, #] &/@ Values@rgrp,
   AssociationQ[lgrp] && Length[Dimensions[rgrp]] == 3,
      Table[GrpMultiply[Values[lgrp][[i]], rgrp[[j]]], {i, Length@Values@lgrp}, {j, Length@rgrp}],
   AssociationQ[lgrp] && Length[Dimensions[rgrp]] == 1,
      Table[GrpMultiply[Values[lgrp][[i]], xyzStr2M4[rgrp[[j]]]], {i, Length@Values@lgrp}, {j, Length@rgrp}],
   Length[Dimensions[lgrp]] == 3 && AssociationQ[rgrp],
      Table[GrpMultiply[lgrp[[i]], Values[rgrp][[j]]], {i, Length@lgrp}, {j, Length@Values@rgrp}],
   Length[Dimensions[lgrp]] == 1 && AssociationQ[rgrp],
      Table[GrpMultiply[xyzStr2M4[lgrp[[i]]], Values[rgrp][[j]]], {i, Length@lgrp}, {j, Length@Values@rgrp}],
   Length[Dimensions[lgrp]] == 3 && Length[Dimensions[rgrp]] == 3,
      Table[GrpMultiply[lgrp[[i]], rgrp[[j]]], {i, Length@lgrp}, {j, Length@rgrp}],
   Length[Dimensions[lgrp]] == 3 && MatrixQ[rgrp],
      GrpMultiply[#, rgrp] & /@ lgrp,
   MatrixQ[lgrp] && Length[Dimensions[rgrp]] == 3,
      GrpMultiply[lgrp, #] & /@ rgrp,
   Length[Dimensions[lgrp]] == 3 && Length[Dimensions[rgrp]] == 1,
      Table[GrpMultiply[lgrp[[i]], xyzStr2M4[rgrp[[j]]]], {i, Length@lgrp}, {j, Length@rgrp}],
   Length[Dimensions[lgrp]] == 1 && Length[Dimensions[rgrp]] == 3,
      Table[GrpMultiply[xyzStr2M4[lgrp[[i]]], rgrp[[j]]], {i, Length@lgrp}, {j, Length@rgrp}],
   Length[Dimensions[lgrp]] == 3 && Length[Dimensions[rgrp]] == 0,
      GrpMultiply[#, xyzStr2M4[rgrp]] &/@ lgrp,
   Length[Dimensions[lgrp]] == 0 && Length[Dimensions[rgrp]] == 3,
      GrpMultiply[xyzStr2M4[lgrp], #] &/@ rgrp,
   StringQ[lgrp] && MatrixQ[rgrp],
      M42xyzStr[xyzStr2M4[lgrp].rgrp],
   MatrixQ[lgrp] && StringQ[rgrp],
      M42xyzStr[lgrp.xyzStr2M4[rgrp]],
   MatrixQ[lgrp] && MatrixQ[rgrp],
      M42xyzStr[lgrp.rgrp],
   MatrixQ[lgrp] && Length[Dimensions[rgrp]] == 1,
      GrpMultiply[lgrp, xyzStr2M4[#]] &/@ rgrp,
   Length[Dimensions[lgrp]] == 1 && MatrixQ[rgrp],
      GrpMultiply[xyzStr2M4[#], rgrp] &/@ lgrp,
   MatrixQ[lgrp] && Length[Dimensions[rgrp]] == 0,
      GrpMultiply[lgrp, xyzStr2M4[rgrp]],
   Length[Dimensions[lgrp]] == 0 && MatrixQ[rgrp],
      GrpMultiply[xyzStr2M4[lgrp], rgrp],
   Length[Dimensions[lgrp]] == 1 && Length[Dimensions[rgrp]] == 1,
      Table[GrpMultiply[xyzStr2M4[lgrp[[i]]], xyzStr2M4[rgrp[[j]]]], {i, Length@lgrp}, {j, Length@rgrp}],
   Length[Dimensions[lgrp]] == 1 && Length[Dimensions[rgrp]] == 0,
      GrpMultiply[#, rgrp] &/@ lgrp,
   Length[Dimensions[lgrp]] == 0 && Length[Dimensions[rgrp]] == 1,
      GrpMultiply[lgrp, #] &/@ rgrp,
   Length[Dimensions[lgrp]] == 0 && Length[Dimensions[rgrp]] == 0,
      GrpMultiply[xyzStr2M4[lgrp], xyzStr2M4[rgrp]]
  ]
]

GrpV[grp_, v_, OptionsPattern[{"k" -> False}]] := Module[{i,j},
  Which[AssociationQ[grp]&&Length[Dimensions@v]==1, GrpV[#,v,"k"->OptionValue["k"]] &/@ Keys[grp],
        AssociationQ[grp]&&Length[Dimensions@v]==2, Table[GrpV[Values[grp][[i]], v[[j]]], {i, Length@Values@grp}, {j, Length@v}],
        StringQ[grp]&&(Length[Dimensions@v]==1), GrpV[xyzStr2M4[grp],v,"k"->OptionValue["k"]],
        StringQ[grp]&&(Length[Dimensions@v]==2), GrpV[xyzStr2M4[grp],#,"k"->OptionValue["k"]] &/@ v,
        (ListQ[grp])&&(!MatrixQ[grp])&&(Length[Dimensions@v]==1), GrpV[#,v,"k"->OptionValue["k"]] &/@ grp,
        (ListQ[grp])&&(!MatrixQ[grp])&&(Length[Dimensions@v]==2), Table[GrpV[grp[[i]], v[[j]]], {i, Length@grp}, {j, Length@v}],
        MatrixQ[grp]&&(Length[Dimensions@v]==2), GrpV[grp,#,"k"->OptionValue["k"]] &/@ v,
        MatrixQ[grp]&&(Length[Dimensions@v]==1), If[OptionValue["k"], grp[[1;;3,1;;3]].v, (grp.Append[v, 1])[[1;;3]]]
  ]
]

GrpO[grp_, O_, vars_] := Module[{sub, xyz},
  Which[AssociationQ[grp], GrpO[#, O, vars] & /@ Keys[grp],
        StringQ[grp], GrpO[xyzStr2M4[grp], O, vars],
        Length@Dimensions[grp] == 1, GrpO[#, O, vars] & /@ grp,
        Length@Dimensions[grp] == 2, 
        xyz = (grp.Join[vars, {1}])[[1 ;; 3]];
        sub = Thread[vars -> xyz];
        Return[Expand[O /. sub]]
        ]
]

GetProjOperator[ireps_] := Module[{g, proj},
  g = Length[ireps[[1]]];
  proj = Length[#[[1]]]/g # &/@ ireps;
  Return[proj]
]

ProjectOnOperator[grp_, ireps_, O_, vars_] := Module[{i, l, proj},
  proj = GetProjOperator[ireps];
  Flatten[Table[Table[{Simplify[(Conjugate[#[[l, l]]] & /@ proj[[i]]).GrpO[grp, O, vars]], ToString[Subscript["\[Phi]"^ToString[i], l], StandardForm]}, {l, Length@proj[[i, 1]]}], {i, Length@proj}], 1]
]

ProjectOnBasis[ireps_, OpMat_, basis_] := Module[{i, l, proj},
  proj = GetProjOperator[ireps];
  Flatten[Table[Table[{Simplify[(Conjugate[#[[l, l]]] & /@ proj[[i]]).(#.basis &/@ OpMat)], ToString[Subscript["\[Phi]"^ToString[i], l], StandardForm]}, {l, Length@proj[[i, 1]]}], {i, Length@proj}],1]
]

GTimes[list_] := Module[{},
  Fold[GrpMultiply, list]
]

GPower[key_, n_] := Module[{},
  Fold[GrpMultiply, ConstantArray[key, n]]
]

GInverse[ele_] := Module[{},
  Which[AssociationQ[ele], GInverse[Values[ele]],
        ListQ[ele]&&!MatrixQ[ele], GInverse[#] &/@ ele,
        StringQ[ele], GInverse[xyzStr2M4[ele]],
        MatrixQ[ele], M42xyzStr[Inverse[ele]]
  ]
]

GetLeftCosets[grp0_, invsub_] := Module[{i, j, grp},
  grp = DeleteDuplicates[ModM4[#] &/@ GTimes[{grp0, invsub}], Complement[#1, #2]==={}&];
  Return[grp]
]

GetSpgIreps[spg_, kvec_] := Module[{grpk, orbits, ireps, irepk, gk, Tk, Tj, TTT, T, R, tran, p, pos, GammaK},
  grpk = GetGrpK[spg, kvec];
  irepk = GetGrpKIreps[spg, kvec];
  gk = Length[grpk];
  orbits = GetLeftCosets[spg, GetGrpK[spg, kvec]]\[Transpose][[1]];
  Print["orbits: ", orbits];
  Print["k star: ", GetStarK[spg, kvec]];
  tran = {ToExpression["\!\(\*SubscriptBox[\(t\), \(1\)]\)"], ToExpression["\!\(\*SubscriptBox[\(t\),  \(2\)]\)"], ToExpression["\!\(\*SubscriptBox[\(t\), \(3\)]\)"]};
  ireps = Table[ArrayFlatten@Table[TTT = GTimes[{Tk, T, GInverse@Tj}];
                R = First@xyzStr2TRot[Tk];
                pos = Position[Keys[grpk], ModM4@TTT];
                GammaK = If[T == "x,y,z", Exp[-I 2 Pi R.kvec.tran], 1];
                If[pos == {}, 0, GammaK irepk[[p, First@First@pos]]], {Tk, orbits}, {Tj, orbits}], {p, Length@irepk}, {T, Keys[spg]}];
  Return[ireps]
]

GetGrpKIreps[spg_, kvec_] := Module[{i, j, k, n, gk, grpk, gH, grpH, perm, irepsn, ireps, pos, \[CapitalGamma]k, lcosets, q, qHq, orbits, lorbit, OrbitIreps, charsfinal, classes},
  grpk = GetGrpK[spg, kvec];
  gk = Length[grpk];
  If[Length[grpk] == 1,
     ireps = {{{{1}}}};
     Return[ireps]];
  
  n = If[Divisible[Length[grpk], 2], 2, 3];
  grpH = GetInvSubGrp[grpk, n, "all"->False];
  grpH = If[grpH == {} && n == 2, n = 3; GetInvSubGrp[grpk, n, "all"->False][[1]], grpH[[1]]];
  gH = Length[grpH];
  
  irepsn = GetGrpKIreps[grpH, kvec];
  lcosets = GetLeftCosets[grpk, grpH];
  q = Which[n==2, {lcosets[[1,1]], lcosets[[2,1]]}, n==3, {lcosets[[1,1]], lcosets[[2,1]], ModM4@GPower[lcosets[[2,1]], 2]}];  
  perm = FindPermutation[ModM4@Flatten[GTimes[{q, grpH}]], Keys[grpk]];

  orbits = Table[qHq = GTimes[{q[[i]], grpH[[j]], GInverse[q[[i]]]}];
    pos = First@First@Position[Keys[grpH], ModM4@qHq];
    \[CapitalGamma]k = Exp[-I 2 Pi kvec.(xyzStr2TRot[qHq][[2]] - xyzStr2TRot[Keys[grpH][[pos]]][[2]])];
    {pos, \[CapitalGamma]k, Keys[grpH][[pos]]}, {i, Length@q}, {j, gH}];
  
  ireps = Rationalize[Flatten[Table[
      OrbitIreps = Table[Rationalize[irepsn[[i, #1]] #2] & @@@ (orbits[[j]]), {j, Length@q}];
      GetOrbitInducedIreps[kvec, grpk, grpH, q, OrbitIreps, perm], {i, Length@irepsn}], 1]];

  ireps = DeleteDuplicates[ireps, Map[Tr, #1] == Map[Tr, #2] &];

  Return[ireps]
]

GetOrbitInducedIreps[kvec_, grp_, grpH_, q_, OrbitIreps_, perm_] := Module[{i, j, k, n, pos, hpos, hMat, qhqMat, qnMat, UMat, lorbit, ireps, InducedIreps, Irepq, IrepGrpH, IrepDim, gH, \[CapitalGamma]k, qn, q2, q3, qh, qqh},
  lorbit = Length@Union[Map[Tr, #] & /@ OrbitIreps];
  n = Length[q];
  gH = Length[grpH];
  IrepDim = Length[OrbitIreps[[1]][[1]]];
  (*debug
  Print["H: ", Grid[{Keys@grpH,GetEleLabel[grpH],Simplify@Tr[#]&/@OrbitIreps[[1]],MatrixForm[#]&/@OrbitIreps[[1]]}]];
  Print["q: ", Grid[{q, GetEleLabel[q]}]];
  *)
  ireps = If[lorbit == 1,
  (* Lenght of orbit 1 *)
    
    hpos = Table[If[OrbitIreps[[1]][[i]] != OrbitIreps[[2]][[i]], i, ##&[]], {i,gH}];
    hpos = If[hpos == {}, gH, hpos[[1]]];
    (*debug
    Print["h: ", GetEleLabel[Keys[grpH][[hpos]]]];
    *)
    hMat = OrbitIreps[[1]][[hpos]];
    qhqMat = OrbitIreps[[2]][[hpos]];
    qn = Which[n == 2, GPower[q[[2]], 2],
               n == 3, GPower[q[[2]], 3]];
    pos = First@First@Position[Keys[grpH], ModM4@qn];
    \[CapitalGamma]k = Exp[-I 2 Pi kvec.(xyzStr2TRot[qn][[2]] - xyzStr2TRot[ModM4@qn][[2]])];
    qnMat = \[CapitalGamma]k OrbitIreps[[1]][[pos]];
    UMat = GetUnitaryTransMatrix[hMat, qhqMat, qnMat, n];
    
    Which[
      n == 2,
      InducedIreps = Table[Flatten[{OrbitIreps[[1]], Table[qh = GTimes[{q[[2]], grpH[[i]]}];
         Exp[I j Pi] Exp[I 2 Pi (xyzStr2TRot[qh][[2]] - xyzStr2TRot[ModM4@qh][[2]]).kvec] UMat.OrbitIreps[[1]][[i]], {i, 1, gH}]}, 1], {j, 0, 1}],
      n == 3,
      InducedIreps = Table[Flatten[{OrbitIreps[[1]], Table[qh = GTimes[{q[[2]], grpH[[i]]}]; 
         Exp[-I 2 Pi j/3] Exp[I 2 Pi (xyzStr2TRot[qh][[2]] - xyzStr2TRot[ModM4@qh][[2]]).kvec] UMat.OrbitIreps[[1]][[i]], {i, 1, gH}],
           Table[qqh = GTimes[{q[[3]], grpH[[i]]}]; 
         Exp[-I 4 Pi j/3] Exp[I 2 Pi (xyzStr2TRot[qqh][[2]] - xyzStr2TRot[ModM4@qqh][[2]]).kvec] UMat.UMat.OrbitIreps[[1]][[i]], {i, 1, gH}]}, 1], {j, 0, 2}];
      ];
    (*debug
    Print["G: ", Grid[Join[{Flatten@GTimes[{q,grpH}]},{GetEleLabel[Flatten@GTimes[{q,grpH}]]},Simplify@Map[Tr,#]&/@InducedIreps, Map[MatrixForm,#]&/@InducedIreps]]];
    *)
    Table[Permute[InducedIreps[[i]], perm], {i, 1, n}],
    
    (* Lenght of orbit 2 and 3*)
    Which[
      n == 2,
      IrepGrpH = Table[ArrayFlatten[{{OrbitIreps[[1]][[i]], 0}, {0, OrbitIreps[[2]][[i]]}}], {i, gH}];
      q2 = GPower[q[[2]], 2];
      pos = First@First@Position[Keys[grpH], ModM4@q2];
      \[CapitalGamma]k = Exp[-I 2 Pi kvec.(xyzStr2TRot[q2][[2]] - xyzStr2TRot[ModM4@q2][[2]])];
      Irepq = ArrayFlatten[{{0, IdentityMatrix[IrepDim]}, {\[CapitalGamma]k OrbitIreps[[1]][[pos]], 0}}];
      InducedIreps = Flatten[{IrepGrpH, Table[qh = GTimes[{q[[2]], grpH[[i]]}];
          \[CapitalGamma]k = Exp[I 2 Pi kvec.(xyzStr2TRot[qh][[2]] - xyzStr2TRot[ModM4@qh][[2]])];
          \[CapitalGamma]k Irepq.IrepGrpH[[i]], {i, gH}]}, 1],
      n == 3,
      IrepGrpH = Table[ArrayFlatten[{{OrbitIreps[[1]][[i]], 0, 0}, {0, OrbitIreps[[2]][[i]], 0}, {0, 0, OrbitIreps[[3]][[i]]}}], {i, gH}];
      q3 = GPower[q[[2]], 3];
      pos = First@First@Position[Keys[grpH], ModM4@q3];
      \[CapitalGamma]k = Exp[-I 2 Pi kvec.(xyzStr2TRot[q3][[2]] - xyzStr2TRot[ModM4@q3][[2]])];
      Irepq = ArrayFlatten[{{0, IdentityMatrix[IrepDim], 0}, {0, 0, IdentityMatrix[IrepDim]}, {\[CapitalGamma]k OrbitIreps[[1]][[pos]], 0, 0}}];
      InducedIreps = Flatten[{IrepGrpH, Table[qh = GTimes[{q[[2]], grpH[[i]]}];
          \[CapitalGamma]k = Exp[I 2 Pi kvec.(xyzStr2TRot[qh][[2]] - xyzStr2TRot[ModM4@qh][[2]])];
          \[CapitalGamma]k Irepq.IrepGrpH[[i]], {i, gH}], Table[qqh = GTimes[{q[[3]], grpH[[i]]}];
          \[CapitalGamma]k = Exp[I 2 Pi kvec.(xyzStr2TRot[qqh][[2]] - xyzStr2TRot[ModM4@qqh][[2]])];
          \[CapitalGamma]k Irepq.Irepq.IrepGrpH[[i]], {i, gH}]}, 1]
    ];
    (*debug
    Print["G: ", Grid[Join[{Flatten@GTimes[{q,grpH}]}, {GetEleLabel[Flatten@GTimes[{q,grpH}]]}, Simplify@Map[Tr,#]&/@{InducedIreps}, Map[MatrixForm,#]&/@{InducedIreps}]]];
    *)
    {Permute[InducedIreps, perm]}
  ];
  Return[ireps];
]

GetUnitaryTransMatrix[Hmat_, qhqmat_, Unmat_, n_] := Module[{mdim, uu, U, inf1, inf2, eq1, eq2, eq3, sol}, 
  mdim = Length[Hmat];
  uu = Table[U[i, j], {i, 1, mdim}, {j, 1, mdim}];
  inf1 = Simplify[uu.Hmat - qhqmat.uu];
  eq1 = Table[inf1[[i, j]] == 0, {i, 1, mdim}, {j, 1, mdim}];
  inf2 = Simplify[MatrixPower[uu, n]];
  eq2 = Table[inf2[[i, j]] - Unmat[[i, j]] == 0, {i, 1, mdim}, {j, 1, mdim}];
  eq3 = Det[uu] == 1;
  sol = If[mdim > 1, FindInstance[Flatten[{eq1, eq2, eq3}], Flatten[uu]], FindInstance[Flatten[{eq1, eq2}], Flatten[uu]]];
  If[sol == {}, Print["No Umat is found for q!"];Print[Flatten[{eq1,eq2,eq3}/.{U[i_,j_]:>ToExpression["u"<>ToString[i]<>ToString[j]]}]];Abort[]];
  Return[uu /. Flatten[sol]]
]        

GetSpgCT[spg0_, k_, ireps_, OptionsPattern[{"print" -> True}]] := Module[{i, l, pos, spgk, lp, characters, head1, head2, classes, kvec, klabel, tran},
  {kvec, klabel} = k;
  spgk = GetGrpK[spg0, kvec];
  lp = Length[ireps];
  classes = GetClasses[spgk];
  characters = Table[pos=First@First@Position[Keys[spgk], classes[[i,1]]];
                     Simplify@Tr@ireps[[l, pos]], {l, lp}, {i, lp}];
  head1 = Superscript[klabel, ToString[#]] & /@ Range[lp];
  head2 = ToString[Length[#]]<>GetEleLabel[#[[1]]] &/@ classes;
  tran = Thread[{ToExpression["\!\(\*SubscriptBox[\(t\), \(1\)]\)"], ToExpression["\!\(\*SubscriptBox[\(t\),  \(2\)]\)"], ToExpression["\!\(\*SubscriptBox[\(t\), \(3\)]\)"]}->{0,0,0}];
  If[OptionValue["print"], Print[TableForm[characters/.tran, TableHeadings -> {head1, head2}, TableAlignments -> Right]]];
  Return[characters]
]

GetCGCoefficients[grp_, ireps_, pp_] := Module[{i, j, k, t, p, classes, ct, clpos, lc, g, character, plist, rep, dim, Es, npq, A, CG, CGList, DirectProductRep},
  classes = GetClasses[grp];
  lc =Length[classes];
  g = Length[grp];
  clpos = First@First@Position[Keys[grp], #] &/@ (First[#] &/@ classes);
  ct = Table[Simplify@Tr[ireps[[i,j]]], {i, lc}, {j, clpos}];

  character = Fold[Times, ct[[#]] &/@ pp];
  plist = DecomposeIreps[grp, ct, character, "print"->True];
  rep = (ireps[[#]] &/@ pp)\[Transpose];
  DirectProductRep = Fold[Simplify@ArrayFlatten[#1 \[TensorProduct] #2]&, #] &/@ rep;

  CGList = Table[If[plist[[p]] == 0, ##&[], 
                 dim = Length[ireps[[p,1]]];
                 A[i_, j_] := dim/g*Sum[DirectProductRep[[t]]*Conjugate[ireps[[p, t, i, j]]], {t, 1, g}];
                 Es = Eigensystem[Simplify@A[1, 1]];
                 npq = Total[Es[[1]]];
                 CG = Simplify@Table[Table[A[i, 1].Es[[2, k]], {i, 1, dim}], {k, 1, npq}];
                 Simplify@Map[Orthogonalize, CG]], {p, lc}];
  Return[CGList]
]

GetCGjmjm2jm[jm1_, jm2_] := Module[{j, m, j1, j2, m1, m2, c1, c2, c3, c4, c5, t1, t2, zmin, zmax, CGSum, CGmat},
  {j1, m1} = jm1;
  {j2, m2} = jm2;
  CGmat=Table[Table[
            If[m != m1+m2,
               0,
               (*c1=j1+j2-j;
               c2=j1-m1;
               c3=j2+m2;
               c4=j2-j-m1;
               c5=j1-j+m2;
               zmin=Max[{0,c4,c5}];
               zmax=Min[c1,c2,c3];
               CGSum=Sum[(-1)^z/(z! (c1-z)! (c2-z)! (c3-z)! (z-c4)! (z-c5)!), {z,zmin,zmax}];
               t1 = (2 j + 1) (j1+j2-j)! (j2+j-j1)! (j+j1-j2)!;
               t2 = (j1+m1)! (j1-m1)! (j2+m2)! (j2-m2)! (j+m)! (j-m)!;
               Sqrt[t1 t2/(j1+j2+j+1)!] CGSum*)
               ClebschGordan[jm1, jm2, {j, m}]
               ], {m, Range[-j,j]}], {j, Range[Abs[j1-j2], j1+j2]}];
  Return[CGmat]
]

DecomposeIreps[grp_, ct_, character_, OptionsPattern[{"print"->True}]] := Module[{i, p, g, classes, np},
  classes = GetClasses[grp];
  lc = Length[classes];
  If[Length[character] != lc, Print["character list length doesn't match with the class!"];Abort[]];
  g = Length[grp];
  np = Rationalize@Simplify[Table[Sum[Length[classes[[i]]]/g Conjugate[character[[i]]] ct[[p,i]], {i,1,lc}], {p,1,lc}]];
  If[OptionValue["print"],
  Print[StringRiffle[Table[ToString[np[[p]]] <> ToString[Superscript["\[CapitalGamma]", ToString[p]], StandardForm], {p,lc}], "\[CirclePlus]"]]];
  Return[np]
]

GetSuperCellGrp[t_] := Module[{xyz, grp},
  xyz = ToString@StringRiffle[Map[ToString[#, StandardForm] &, {ToExpression["x"], ToExpression["y"],  ToExpression["z"]} + #], ","] & /@ t;
  grp = xyzStr2Grp[M42xyzStr[xyzStr2M4[#]] & /@ xyz];
  Return[grp]
]

GetEulerVector[ele_] := Module[{},
  Which[AssociationQ[ele],
        GetEulerVector[#] &/@ Keys[ele],
        Length@Dimensions[ele] == 1,
        GetEulerVector[#] &/@ ele,
        Length@Dimensions[ele] == 2,
        Mat2EulerVector[ele[[1;;3,1;;3]]],
        Length@Dimensions[ele] == 3,
        GetEulerVector[#] &/@ ele,
        StringQ[ele],
        GetEulerVector[xyzStr2M4[ele]]
       ]
]

(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
