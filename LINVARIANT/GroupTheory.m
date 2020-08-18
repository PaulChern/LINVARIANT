BeginPackage["LINVARIANT`GroupTheory`", {"LINVARIANT`SphericalHarmonics`", "LINVARIANT`MathematicaPlus`", "LINVARIANT`Structure`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
GetWyckoffImages           ::usage = "GetWyckoffImages[grp, pos]"
CifImportSpg               ::usage "CifImportSpg[file]"
GrpQ                       ::usage "GrpQ[grp]"
GetGrpK                    ::usage "GetGrpK[grp0, vec0, lattice0]"
GetStarK                   ::usage "GetStarK[grp0, k]"
GetSiteSymmetry            ::usage "GetSiteSymmetry[grp0, vec0]"
CifImportOperators         ::usage "CifImportOperators[file]"
GetGenerators              ::usage "GetGenerators[grp]"
Generators2Grp             ::usage "Generators2Grp[gen]"
DoubleGroupQ               ::usage "DoubleGroupQ[grp]"
xyz2Grp                    ::usage "xyz2Grp[keys]"
xyz2Expr                   ::usage "xyz2Expr[xyzStrData]"
Expr2Rot                   ::usage "Expr2Rot[expr]"
xyz2Mat                    ::usage "xyz2Mat[xyzStrData]"
xyz2RotTSU2                ::usage "xyz2RotTSU2[xyzStrData]"
RotTran2M4                 ::usage "RotTran2M4[rot, tran]"
Mat2xyz                    ::usage "Mat2xyz[m4]"
Mat2xyzEle                 ::usage "Mat2xyzEle[m4]"
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
xyz2SU2                    ::usage "xyz2SU2[latt, grp]"
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
PrintIreps                 ::usage "PrintIreps[classes, kvec, irreps, character, reality]"
DecomposeIreps             ::usage "DecomposeIreps[grp, ct, character]"
GetCGCoefficients          ::usage "GetCGCoefficients[grp, ireps, p1p2]"
GetCGjmjm2jm               ::usage "GetCGjmjm2jm[jm1, jm2]"
GetSuperCellGrp            ::usage "GetSuperCellGrp[t]"
GetEulerVector             ::usage "GetEulerVector[ele]"
GetCellFromGrp             ::usage "GetCellFromGrp[grp]"
GetDoubleGrp               ::usage "GetDoubleGrp[grp, latt]"
(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
ModM4[ele_] := Module[{m4, su2},
  Which[AssociationQ[ele], ModM4[#] &/@ Keys[ele],
        ListQ[ele], ModM4[#] &/@ ele,
        StringQ[ele],
        {m4, su2} = xyz2Mat[ele];
        m4[[1 ;; 3, 4]] = Mod[m4[[1 ;; 3, 4]], 1];
        Return[Mat2xyz[{m4, su2}]]
       ]
]

ExpandSimplify[in_] := Expand[Simplify[in]];

GetWyckoffImages[grp0_, wyckoff_] := Module[{i, j, k, cell, grpt, grp, xyzStrData, xyzExpData, WyckoffImages0, WyckoffImages1},
  cell = GetCellFromGrp[grp0];
  grp = Which[ListQ[grp0], Flatten@GTimes[grp0], AssociationQ[grp0], Keys[grp0]];
  WyckoffImages0 = Table[Union[{ModCell[#[[1]].wyckoff[[i,1]]+#[[2]], cell], wyckoff[[i,2]]} &/@ xyz2RotTSU2@grp, SameTest -> (Norm[Chop[#1[[1]] - #2[[1]]]] == 0. &)], {i, Length[wyckoff]}];
  WyckoffImages1 = Table[MapIndexed[{#1[[1]], SimplifyElementSymbol[#1[[2]]] <> ToString[i] <> "." <> ToString[First@#2]} &,  WyckoffImages0[[i]]], {i, Length[wyckoff]}];
  Return[WyckoffImages1]
]


(*GetWyckoffImages[grp_, pos_] := Module[{xyzStrData, xyzExpData, WyckoffImages, xyzvar},
  xyzvar = {"x", "y", "z"};
  xyzExpData = ToExpression["{" <> # <> "}"][[1;;3]] &/@ Keys[grp];
  WyckoffImages = Flatten[Table[Sort[Union[{Mod[#,1], pos[[i, 2]]} & /@ (xyzExpData /. Thread[ToExpression[xyzvar] -> pos[[i, 1]]])], Total[#1[[1]]] > Total[#2[[1]]] &], {i, Length[pos]}], 1];
  Return[WyckoffImages];
]
*)

GetClasses[grp0_] := Module[{dgQ, i, j, g, bab, cl},
  g = Length[grp0];
  dgQ = DoubleGroupQ[grp0];
  bab = ParallelTable[ModM4@GTimes[{Keys[grp0][[i]],Keys[grp0][[j]],GInverse[Keys[grp0][[i]]]}], {j, 1, g}, {i, 1, g}, DistributedContexts -> {"LINVARIANT`GroupTheory`Private`"}];
  cl = Flatten[Values[KeySortBy[GroupBy[#, First@Union[GetEulerVector[#]\[Transpose][[2]]] &], Minus]]
               & /@ Values[KeySortBy[GroupBy[#, First@Union[GetEulerVector[#]\[Transpose][[3]]] &], Minus]]
               & /@ Values[KeySort[GroupBy[Union[Map[Union, bab]], Length[#] &]]], 3];
  cl = Flatten[(Values[KeySortBy[GroupBy[#, GetEulerVector[#][[1]] &], Minus]]
               & /@ Values[KeySortBy[GroupBy[#, GetEulerVector[#][[2]] &], Minus]]), 2] & /@ cl;
  cl = If[dgQ, 
          Prepend[DeleteCases[cl, {"x,y,z,up,dn"}], {"x,y,z,up,dn"}],
          Prepend[DeleteCases[cl, {"x,y,z"}], {"x,y,z"}]
       ];
  Return[cl]
]

GetPGCharacters[grp0_, latt_, OptionsPattern[{"print"->True}]] := Module[{m, n, i, j, k, h, g, nc, Classes, H0Mat, H1Mat, HMat, Hi, hh, Vt, V, Vi, cc, vec, IrepDim, CharacterTable, CharacterTableSorted, ct, IrepLabel, ClassesLabel},
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
  Print[TableForm[ct, TableHeadings -> {IrepLabel, ToString[Length[#]]<>GetEleLabel[First[#], latt] & /@ Classes}, TableAlignments -> {Right, Center}]]];
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

xyz2Expr[xyzStrData_] := Module[{xyzRotTranData, xyzTranslation, xyzRotData, SpinRot, dgQ},
  dgQ = DoubleGroupQ[xyzStrData];
  xyzRotTranData = ToExpression["{" <> xyzStrData <> "}"]; 
  xyzTranslation = xyzRotTranData[[1;;3]] /. {ToExpression["x"]->0, ToExpression["y"]->0, ToExpression["z"]->0};
  xyzRotData = xyzRotTranData[[1;;3]] - xyzTranslation;
  SpinRot = If[dgQ, xyzRotTranData[[4;;5]], Null];
  Return[{xyzRotData, xyzTranslation, SpinRot}]
]

xyz2Mat[xyzStrData_] := Module[{RT, rot, tran, su2},
  Which[ListQ[xyzStrData], 
        xyz2Mat[#] &/@ xyzStrData, 
        StringQ[xyzStrData], 
        {rot, tran, su2} = xyz2Expr[xyzStrData];
        RT = Join[Join[Expr2Rot[rot]\[Transpose], {tran}]\[Transpose], {{0, 0, 0, 1}}];
        RT = RotTran2M4[Expr2Rot[rot], tran];
        {RT, Expr2Rot[su2]},
        True, 
        Print["xyz2Mat: Wrong Input type!"];Abort[]
  ]
]

xyz2RotTSU2[xyzStrData_] := Module[{rot, tran, su2},
  Which[StringQ[xyzStrData], 
        {rot, tran, su2} = xyz2Expr[xyzStrData];
        {Expr2Rot[rot], tran, Expr2Rot[su2]},
        ListQ[xyzStrData], 
        xyz2RotTSU2[#] &/@ xyzStrData, 
        True, 
        Print["xyz2RotTSU2: Wrong Input type!"];Abort[]
  ]
]

RotTran2M4[rot_, tran_] := Module[{RT},
  RT = Join[Join[rot\[Transpose], {tran}]\[Transpose], {{0, 0, 0, 1}}];
  Return[Rationalize[RT]]
]

Mat2xyzEle[mat_] := Module[{i, rot, tran, xyzvar, svar, m4, su2, xyz, updn},
  {m4, su2} = mat;
  xyzvar = {ToExpression["x"], ToExpression["y"], ToExpression["z"]};
  svar = {ToExpression["up"], ToExpression["dn"]};
  rot = ToString[#] & /@ (Rationalize[m4.Join[xyzvar, {1}]] - (Rationalize[m4.Join[xyzvar, {1}]] /. Thread[xyzvar->{0,0,0}]));
  tran = Rationalize[m4.Join[xyzvar, {1}]] /. Thread[xyzvar->{0,0,0}];
  tran = Which[#1 == 0, "",
               #1 != 0 && #2 == 1, If[Sign[#1] > 0, "+", Unevaluated[Sequence[]]] <> ToString[#1],
               True, If[Sign[#1] > 0, "+", Unevaluated[Sequence[]]] <> ToString[#1] <> "/" <> ToString[#2]]
         & @@@ (Through[{Numerator, Denominator}[#]] & /@ tran)[[1 ;; 3]];
  xyz = StringDelete[StringJoin@Riffle[Table[rot[[i]] <> tran[[i]], {i, 3}],","], " "];
  If[su2 === Null, 
     Return[xyz],
     su2 = Map[Complex2Exp[ExpandSimplify[#]]&, #] &/@ su2;
     updn = StringJoin[Riffle[ToString[#, StandardForm] & /@ (su2.svar), ","]];
     Return[xyz<>","<>updn]
  ]
]

Mat2xyz[mat_] := Module[{},
  Which[Length@Dimensions[mat]==1, 
        Mat2xyzEle[mat], 
        Length@Dimensions[mat]==2,
        Mat2xyz[#] &/@ mat]
]

DoubleGroupQ[grp_] := Module[{dgQ},
  dgQ= Which[AssociationQ[grp], Length@ToExpression["{" <> Keys[grp][[1]] <> "}"] > 3,
             ListQ[grp], AllTrue[grp, Length[ToExpression["{" <> # <> "}"]] > 3 &],
             StringQ[grp], Length@ToExpression["{" <> grp <> "}"] > 3];
  Return[dgQ]
]

Expr2Rot[expr_] := Module[{m},
  m = Which[Length[expr]==3,
            Coefficient[#, {ToExpression["x"], ToExpression["y"], ToExpression["z"]}] &/@ expr,
            Length[expr]==2,
            Coefficient[#, {ToExpression["up"], ToExpression["dn"]}] &/@ expr,
            Length[expr]==0,
            Null];
  If[(Abs[Det[m]] != 1)&&(m!=Null), Print["Expr2Rot gives wrong determinant "<>ToString[Det[m]]]; Abort[], Unevaluate[Sequence[]]];
  Return[m]
]

CifImportSpg[file_] := Module[{xyzStrData, ele, RT},
  xyzStrData = CifImportOperators[file];
  ele = xyz2Grp[ModM4[#] &/@ xyzStrData, "fast"->True];
  Return[SortGrp[ele]]
]

GetGrpK[grp0_, vec0_] := Module[{g, trvec, eq, sol, posmap, grpk, c, vec},
  trvec = GrpV[grp0, vec0, "k" -> True];
  eq = Thread[# == vec0 + Table[c[i],{i,3}]] &/@ trvec;
  sol = Rationalize[Round[Table[c[j], {j, 3}] /. Flatten[Solve[#, Table[c[j], {j, 3}]]] &/@ eq, 10^-3]];
  posmap = AllTrue[#, IntegerQ] &/@ sol;
  grpk = xyz2Grp@Extract[Keys[grp0], Position[posmap,True]];
  Return[grpk]
]

GetStarK[grp0_, k_] := Module[{grpk, LeftCoset, LeftRep},
  grpk = GetGrpK[grp0, k];
  LeftCoset = GetLeftCosets[grp0, grpk];
  LeftRep = LeftCoset[[;; , 1]];
  Return[Mod[GrpV[LeftRep, k, "k"->True],1]]
]

GetSiteSymmetry[grp0_, vec0_] := Module[{grp, g, trvec, eq, sol, add, grpout, c, vec},

  g = Length[grp0];
  grp = #[[1]][[1;;3,1;;3]] &/@ Values[grp0];

  (*--- Main algortihm ---*)
  trvec = Table[Inverse[grp[[i]]].vec0, {i, 1, g}];
  eq = Table[Thread[trvec[[i]] == vec0 + Table[c[i],{i,3}]], {i, 1,g}];
  sol = Rationalize[Round[Table[Table[c[j], {j, 3}] /. Flatten[Solve[eq[[i]], Table[c[j], {j, 3}]]], {i, 1, g}], 10^-3]];
  add = Table[IntegerQ[sol[[i, 1]]] && IntegerQ[sol[[i, 2]]] && IntegerQ[sol[[i, 3]]], {i, 1, g}];
  grpout = <||>;
  Do[If[add[[i]], grpout = Append[grpout, Keys[grp0][[i]]->Values[grp0][[i]]],None], {i, 1, g}];

  Return[grpout]
]

GetEleLabel[xyz_, latt_] := Module[{det, m4, su2, axes, phi, T, CS, label, dgQ},
  Which[ListQ[xyz],
        GetEleLabel[#, latt] &/@ xyz,
        AssociationQ[xyz],
        GetEleLabel[#, latt] &/@ Keys[xyz],
        StringQ[xyz],
        {m4, su2} = xyz2Mat[xyz];
        T = Rationalize[m4[[1 ;; 3, 4]]];
        dgQ = Which[su2===Null, 1, su2==xyz2SU2[latt, xyz], 1, su2==Simplify[-xyz2SU2[latt, xyz]], -1];
        {axes, phi, det} = Mat2EulerVector[m4[[1;;3, 1;;3]]];
        If[phi != 0,
             CS = If[det == 1, If[dgQ==1, "C", ToString[OverBar["C"], StandardForm]], If[dgQ==1, "S", ToString[OverBar["S"], StandardForm]]]; 
             label = \!\(\*TagBox[StyleBox[RowBox[{"\"\<\\!\\(\\*SubsuperscriptBox[\\(\>\"", "<>", "CS", "<>", "\"\<\\), \\(\>\"", "<>", RowBox[{"ToString", "[", RowBox[{"If", "[", RowBox[{RowBox[{"phi", "\\[NotEqual]", "0"}], ",", RowBox[{"2", " ", RowBox[{"Pi", "/", "phi"}]}], ",", "0"}], "]"}], "]"}], "<>", "\"\<\\), \\(\>\"", "<>", RowBox[{"StringReplace", "[", RowBox[{RowBox[{"StringJoin", "[", RowBox[{"Map", "[", RowBox[{"ToString", ",", "axes"}], "]"}], "]"}], ",", RowBox[{"\"\<-1\>\"", "\\[Rule]", "\"\<i\>\""}]}], "]"}], "<>", "\"\<\\)]\\)\>\""}], ShowSpecialCharacters->False,ShowStringCharacters->True, NumberMarks->True], FullForm]\) <> "{" <> StringJoin[Map[ToString, StringRiffle[T, ","]]] <> "}",
             label = If[det == 1, If[dgQ==1, "E", ToString[OverBar["E"], StandardForm]], If[dgQ==1, "I", ToString[OverBar["I"], StandardForm]]] <> "{" <> StringJoin[Map[ToString, StringRiffle[T, ","]]] <> "}"];
        Return[label]]
];

GetRegRep[grp_] := Module[{gt, g, i, j, k},
  g = Length[grp];
  gt = Expand@Table[GTimes[{Keys[grp][[i]], GInverse[Keys[grp][[j]]]}], {i, 1, g}, {j, 1, g}];
  Return[Table[Table[If[gt[[i, j]] == Keys[grp][[k]], 1, 0], {i, 1, g}, {j, 1, g}], {k, 1, g}]]
]

GetEleOrder[grp_] := Module[{ord, grpn},
  Which[AssociationQ[grp],
        GetEleOrder[#] &/@ Keys[grp],
        ListQ[grp],
        GetEleOrder[#] &/@ grp,
        StringQ[grp],
        ord = 1;
        grpn = grp;
        While[grpn != "x,y,z,up,dn", grpn = ModM4@GTimes[{grp, grpn}]; ord = ord + 1];
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
                              lp/g Conjugate[ct[[p, m]]] RegRep[[i]][[1]], {i, g}];
                (*Find an eigenvector with non degenerate eigenvalue and apply other elements until lp linear independent vectors were found*)
                Do[es = Eigensystem[RegRep[[i]][[1]]];
                 Do[ev = ProjMat.es[[2, j]];
                  If[Norm[ev] == 0, Continue[]];
                  vv = Complement[Orthogonalize[Table[RegRep[[k]][[1]].ev, {k, 1, g}]], {Table[0, {i, 1, g}]}];
                  If[Length[vv] == lp, Break[], None], {j, 1, g}], {i, g, g - 1, -1}];
                (*---calculate representation matrices---*)
                Table[vv.RegRep[[i]][[1]].ConjugateTranspose[vv], {i, 1, g}], {p, 1, h}];
  
  Return[ireps]
]

GetElePosition[grp_, ele_, latt_] := Module[{ind, type, axis, ang, shift, dgQ},
  Which[StringQ[ele],
        type = StringContainsQ[ele, "{"];
        ind = If[type,
                 Position[GetEleLabel[grp, latt], ele],
                 Position[Keys[grp], ele]],
        ListQ[ele],
        GetElePosition[grp, #, latt] &/@ ele,
        AssociationQ[ele],
        GetElePosition[grp, #, latt] &/@ Keys[ele]
       ];
  If[ind == {}, Print["No such element"]; Abort[]];
  Return[First@First@ind]
]

GrpQ[grp_] := Module[{},
  Return[Length@Complement[Union[ModM4@Flatten[GTimes[{grp, grp}]]], Which[AssociationQ[grp], Keys[grp], ListQ[grp], grp]] == 0]
]

GetSubGrp[grp_, ord_: {}] := Module[{GrpSub, dgQ, g, div, sub, sg, test, j, IdentiyEle},
  dgQ = DoubleGroupQ[grp];
  IdentityEle = If[dgQ, "x,y,z,up,dn", "x,y,z"];
  g = Length[grp];
  div = Divisors[g];
  div = Complement[If[ord == {}, Print["no index are specified, this may take very long!"]; Divisors[g], Intersection[Divisors[g], ord]], {1, g}];
  sg = {};
  Do[sub = Subsets[Complement[Keys[grp], {IdentityEle}], {div[[i]] - 1}];
    Do[test = Union[sub[[j]], {IdentityEle}];
      AppendTo[sg, If[GrpQ[test], xyz2Grp@test, {}]], {j, Length@sub}], {i, Length@div}];
  Print[ToString[Length@Flatten[sg]] <> " sub group found."];
  Return[Flatten@sg]
]

xyz2Grp[keys_, OptionsPattern[{"fast"->False}]] := Module[{xyz, dgQ, IdentityEle},
  dgQ = DoubleGroupQ[keys];
  IdentityEle = If[dgQ, "x,y,z,up,dn", "x,y,z"];
  xyz = Prepend[DeleteCases[keys, IdentityEle], IdentityEle];
  If[OptionValue["fast"],
     Return[Association[Thread[xyz->xyz2Mat[xyz]]]],
     If[GrpQ[xyz],
     Return[Association[Thread[xyz->xyz2Mat[xyz]]]],
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

  grpstd = Association[Thread[cl->xyz2Mat[cl]]];
  Return[grpstd]
];

GetInvSubGrp[grp_, ind_:0, OptionsPattern[{"all"->True}]] := Module[{g, GrpListxyzStr, classes, subs, invsub, test, invfound, ipos, dgQ, IdentityEle},
  dgQ = DoubleGroupQ[grp];
  IdentityEle = If[dgQ, "x,y,z,up,dn", "x,y,z"];
  g = Length[grp];
  If[ind != 0, 
     If[Divisible[g, ind], Unevaluated[Sequence[]], Print["Index is not compatible!"];Abort[]], 
     Unevaluated[Sequence[]]];
  classes = DeleteCases[GetClasses[grp], {IdentityEle}];

  subs = If[ind ==0, Union[{{IdentityEle}}, #] &/@ Subsets[classes], If[Length[Flatten[#]] == g/ind-1, Union[{{IdentityEle}}, #], {}] &/@ Subsets[classes]];
  invsub = {};
  ipos = 0;
  invfound = False;
  If[OptionValue["all"],
     Do[test = Flatten[subs[[i]]];
        If[GrpQ[test], 
           If[Length[test] >= 1 && Length[test] <= g, 
              AppendTo[invsub, xyz2Grp@test],
              Null], 
           Null], {i, Length@subs}],
     While[Not[invfound]&&ipos<Length[subs]&&Length[subs]>0, 
           ipos = ipos + 1;
           test = Flatten[subs[[ipos]]];
           If[GrpQ[test], If[Length[test] >= 1 && Length[test] <= g, invfound = True;invsub={xyz2Grp@test}]];]
  ];

  Return[invsub]
]

Generators2Grp[gen_] := Module[{grp0, grp1, ln, i, j, ord1, ord0, grpout},
   grp0 = Union@ModM4@Flatten[Table[GTimes@Table[#, i], {i, GetEleOrder[#]}] & /@ gen];
   ord0 = Length[grp0];
   ord1 = 0;
   While[ord0 > ord1,
    grp1 = Union@Flatten@Table[ModM4@GTimes[{Keys[grp0][[i]],Keys[grp0][[j]]}], {i, 1, Length[grp0]}, {j, 1, Length[grp0]}];
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
  t = ModM4[#] &/@ FoldList[GrpMultiply, Table[Keys[grp][[g]], GetEleOrder[Keys[grp][[g]]]]];
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

xyz2SU2[latt_, grp_] := Module[{},
  GetAngularMomentumRep[latt, grp, 1/2, "Spherical"]
]

GetAngularMomentumRep[latt_, grp_, j_, Harmonic_: "Tesseral"] := Module[{rot, tran, su2, \[Alpha], \[Beta], \[Gamma], \[Epsilon], m1, m2, Dlmn, T},
  Which[AssociationQ[grp],
        GetAngularMomentumRep[latt, #, j, Harmonic] &/@ Keys[grp],
        ListQ[grp],
        GetAngularMomentumRep[latt, #, j, Harmonic] &/@ grp,
        StringQ[grp],
        {rot, tran, su2} = xyz2RotTSU2[grp];
        {{\[Alpha], \[Beta], \[Gamma]}, \[Epsilon]} = Mat2EulerAngles[latt, rot];
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

GetAngularMomentumChars[grp_, j_] := Module[{char, jm, mat, axis, \[Phi], \[Epsilon], dgQ},
  Which[AssociationQ[grp],
        GetAngularMomentumChars[#, j] &/@ Keys[grp],
        ListQ[grp],
        GetAngularMomentumChars[#, j] &/@ grp,
        StringQ[grp],
        mat = First@xyz2RotTSU2[grp][[1;;3,1;;3]];
        {axis, \[Phi], \[Epsilon]} = Mat2EulerVector[mat];
        (*char = \[Epsilon]^j*Sin[(j+1/2)*\[Phi]]/Sin[\[Phi]/2];*)
        char = Simplify[\[Epsilon]^j*Sum[Exp[-I jm \[Phi]], {jm, -j, j}]];
        Return[char]
    ]

]

GrpMultiply[lgrp_, rgrp_] := Module[{i, j, dgQ, m4, su2},
  Which[
   AssociationQ[lgrp] && AssociationQ[rgrp],
      Table[GrpMultiply[Keys[lgrp][[i]], Keys[rgrp][[j]]], {i, Length@lgrp}, {j, Length@rgrp}],
   AssociationQ[lgrp] && StringQ[rgrp],
      GrpMultiply[#, rgrp] &/@ Keys[lgrp],
   StringQ[lgrp] && AssociationQ[rgrp],
      GrpMultiply[lgrp, #] &/@ Keys[rgrp],
   AssociationQ[lgrp] && ListQ[rgrp],
      Table[GrpMultiply[Keys[lgrp][[i]], rgrp[[j]]], {i, Length@lgrp}, {j, Length@rgrp}],
   ListQ[lgrp] && AssociationQ[rgrp],
      Table[GrpMultiply[lgrp[[i]], Keys[rgrp][[j]]], {i, Length@lgrp}, {j, Length@rgrp}],
   ListQ[lgrp] && ListQ[rgrp],
      Table[GrpMultiply[lgrp[[i]], rgrp[[j]]], {i, Length@lgrp}, {j, Length@rgrp}],
   ListQ[lgrp] && StringQ[rgrp],
      GrpMultiply[#, rgrp] &/@ lgrp,
   StringQ[lgrp] && ListQ[rgrp],
      GrpMultiply[lgrp, #] &/@ rgrp,
   StringQ[lgrp] && StringQ[rgrp],
      If[DoubleGroupQ[lgrp] != DoubleGroupQ[rgrp], 
         Print["Error: double group times single group"];Abort[],
         dgQ = DoubleGroupQ[lgrp];
         m4 = xyz2Mat[lgrp][[1]].xyz2Mat[rgrp][[1]];
         su2 = If[dgQ, xyz2Mat[lgrp][[2]].xyz2Mat[rgrp][[2]], Null];
         Mat2xyz[{m4, su2}]
      ]
  ]
]

GrpV[grp_, v_, OptionsPattern[{"k" -> False}]] := Module[{i,j},
  Which[AssociationQ[grp]&&Length[Dimensions@v]==1, GrpV[#,v,"k"->OptionValue["k"]] &/@ Keys[grp],
        AssociationQ[grp]&&Length[Dimensions@v]==2, Table[GrpV[Keys[grp][[i]], v[[j]]], {i, Length@grp}, {j, Length@v}],
        ListQ[grp]&&(Length[Dimensions@v]==1), GrpV[#,v,"k"->OptionValue["k"]] &/@ grp,
        ListQ[grp]&&(Length[Dimensions@v]==2), Table[GrpV[grp[[i]], v[[j]]], {i, Length@grp}, {j, Length@v}],
        StringQ[grp]&&(Length[Dimensions@v]==2), GrpV[grp,#,"k"->OptionValue["k"]] &/@ v,
        StringQ[grp]&&(Length[Dimensions@v]==1), If[OptionValue["k"], 
                                                    xyz2Mat[grp][[1]][[1;;3,1;;3]].v, 
                                                    (xyz2Mat[grp][[1]].Append[v, 1])[[1;;3]]
                                                   ]
  ]
]

GrpO[grp_, O_, vars_] := Module[{sub, xyz, spinorQ},
  spinorQ = ListQ[O]&&(Length[O]==2);
  Which[AssociationQ[grp], GrpO[#, O, vars] & /@ Keys[grp],
        ListQ[grp], GrpO[#, O, vars] & /@ grp,
        StringQ[grp], 
        xyz = (xyz2Mat[grp][[1]].Join[vars, {1}])[[1 ;; 3]];
        sub = Thread[vars -> xyz];
        Simplify@If[spinorQ, xyz2Mat[grp][[2]].(Expand[O /. sub]), Expand[O /. sub]]
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

GPower[grp_, n_] := Module[{},
  Which[AssociationQ[grp], GPower[#, n] &/@ Keys[grp],
        ListQ[grp], GPower[#,n] &/@ grp,
        StringQ[grp], Fold[GrpMultiply, ConstantArray[grp, n]]
  ]
]

GInverse[ele_] := Module[{m4, su2, xyz, dgQ},
  Which[AssociationQ[ele], GInverse[Keys[ele]],
        ListQ[ele], GInverse[#] &/@ ele,
        StringQ[ele], 
        {m4, su2} = xyz2Mat[ele];
        If[su2 === Null, 
           Mat2xyz[{Simplify@Inverse[m4], Null}], 
           Mat2xyz[{Simplify@Inverse[m4], Inverse[su2]}]
        ]
  ]
]

GetLeftCosets[grp0_, invsub_] := Module[{i, j, grp},
  grp = DeleteDuplicates[ModM4[#] &/@ GTimes[{grp0, invsub}], Complement[#1, #2]==={}&];
  Return[grp]
]

GetSpgIreps[spg_, k_, latt_, OptionsPattern[{"print"->True}]] := Module[{grpk, classes, kvec, characters, reality, klabel, orbits, ireps, irepk, gk, Tk, Tj, TTT, T, R, tran, l, lp, lc, p, pos, GammaK, dgQ, IdentityEle},
  dgQ = DoubleGroupQ[spg];
  IdentityEle = If[dgQ, "x,y,z,up,dn", "x,y,z"];
  {kvec, klabel} = k;
  classes = GetClasses[spg];
  grpk = GetGrpK[spg, kvec];
  irepk = GetGrpKIreps[spg, kvec];
  lp = Length[irepk];
  gk = Length[grpk];
  lc = Length[classes];
  reality = Table[1/gk Sum[pos=First@First@Position[Keys[grpk], ModM4@GPower[T, 2]];
                           Simplify[Tr@irepk[[l, pos]]], {T, Keys[grpk]}], {l, lp}];
  orbits = GetLeftCosets[spg, grpk]\[Transpose][[1]];
  Print["orbits: ", orbits];
  Print["k star: ", GetStarK[spg, kvec]];
  tran = {ToExpression["\!\(\*SubscriptBox[\(t\), \(1\)]\)"], ToExpression["\!\(\*SubscriptBox[\(t\),  \(2\)]\)"], ToExpression["\!\(\*SubscriptBox[\(t\), \(3\)]\)"]};
  ireps = Table[ArrayFlatten@Table[TTT = GTimes[{Tk, T, GInverse@Tj}];
                R = First@xyz2RotTSU2[Tk];
                pos = Position[Keys[grpk], ModM4@TTT];
                GammaK = If[T == IdentityEle, Exp[-I 2 Pi R.kvec.tran], 1];
                If[pos == {}, 0, GammaK irepk[[p, First@First@pos]]], {Tk, orbits}, {Tj, orbits}], {p, Length@irepk}, {T, Keys[spg]}];
  characters = Table[pos=First@First@Position[Keys[spg], classes[[i,1]]];
                     Simplify@Tr@ireps[[l, pos]], {l, lp}, {i, lc}];
  If[OptionValue["print"], PrintIreps[latt, classes, k, ireps, characters, reality]];
  Return[{ireps, characters}]
]

GetGrpKIreps[spg_, kvec_] := Module[{i, j, k, n, gk, grpk, gH, grpH, perm, irepsn, ireps, pos, \[CapitalGamma]k, lcosets, q, qHq, orbits, lorbit, OrbitIreps, charsfinal, classes},
  grpk = GetGrpK[spg, kvec];
  gk = Length[grpk];
  If[Length[grpk] == 1,
     ireps = {{{{1}}}};
     Return[ireps]];
  
  n = If[Divisible[Length[grpk], 2], 2, 3];
  grpH = GetInvSubGrp[grpk, n, "all"->True];
  grpH = If[grpH == {} && n == 2, n = 3; GetInvSubGrp[grpk, n, "all"->True][[-1]], grpH[[-1]]];
  gH = Length[grpH];

  irepsn = GetGrpKIreps[grpH, kvec];
  lcosets = GetLeftCosets[grpk, grpH];
  q = Which[n==2, {lcosets[[1,1]], lcosets[[2,1]]}, n==3, {lcosets[[1,1]], lcosets[[2,1]], ModM4@GPower[lcosets[[2,1]], 2]}];  
  perm = FindPermutation[ModM4@Flatten[GTimes[{q, grpH}]], Keys[grpk]];

  orbits = Table[qHq = GTimes[{q[[i]], Keys[grpH][[j]], GInverse[q[[i]]]}];
    pos = First@First@Position[Keys[grpH], ModM4@qHq];
    \[CapitalGamma]k = Exp[-I 2 Pi kvec.(xyz2RotTSU2[qHq][[2]] - xyz2RotTSU2[Keys[grpH][[pos]]][[2]])];
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
    \[CapitalGamma]k = Exp[-I 2 Pi kvec.(xyz2RotTSU2[qn][[2]] - xyz2RotTSU2[ModM4@qn][[2]])];
    qnMat = \[CapitalGamma]k OrbitIreps[[1]][[pos]];
    UMat = GetUnitaryTransMatrix[hMat, qhqMat, qnMat, n];
    
    Which[
      n == 2,
      InducedIreps = Table[Flatten[{OrbitIreps[[1]], Table[qh = GTimes[{q[[2]], Keys[grpH][[i]]}];
         Exp[I j Pi] Exp[I 2 Pi (xyz2RotTSU2[qh][[2]] - xyz2RotTSU2[ModM4@qh][[2]]).kvec] UMat.OrbitIreps[[1]][[i]], {i, 1, gH}]}, 1], {j, 0, 1}],
      n == 3,
      InducedIreps = Table[Flatten[{OrbitIreps[[1]], Table[qh = GTimes[{q[[2]], Keys[grpH][[i]]}]; 
         Exp[-I 2 Pi j/3] Exp[I 2 Pi (xyz2RotTSU2[qh][[2]] - xyz2RotTSU2[ModM4@qh][[2]]).kvec] UMat.OrbitIreps[[1]][[i]], {i, 1, gH}],
           Table[qqh = GTimes[{q[[3]], Keys[grpH][[i]]}]; 
         Exp[-I 4 Pi j/3] Exp[I 2 Pi (xyz2RotTSU2[qqh][[2]] - xyz2RotTSU2[ModM4@qqh][[2]]).kvec] UMat.UMat.OrbitIreps[[1]][[i]], {i, 1, gH}]}, 1], {j, 0, 2}];
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
      \[CapitalGamma]k = Exp[-I 2 Pi kvec.(xyz2RotTSU2[q2][[2]] - xyz2RotTSU2[ModM4@q2][[2]])];
      Irepq = ArrayFlatten[{{0, IdentityMatrix[IrepDim]}, {\[CapitalGamma]k OrbitIreps[[1]][[pos]], 0}}];
      InducedIreps = Flatten[{IrepGrpH, Table[qh = GTimes[{q[[2]], Keys[grpH][[i]]}];
          \[CapitalGamma]k = Exp[I 2 Pi kvec.(xyz2RotTSU2[qh][[2]] - xyz2RotTSU2[ModM4@qh][[2]])];
          \[CapitalGamma]k Irepq.IrepGrpH[[i]], {i, gH}]}, 1],
      n == 3,
      IrepGrpH = Table[ArrayFlatten[{{OrbitIreps[[1]][[i]], 0, 0}, {0, OrbitIreps[[2]][[i]], 0}, {0, 0, OrbitIreps[[3]][[i]]}}], {i, gH}];
      q3 = GPower[q[[2]], 3];
      pos = First@First@Position[Keys[grpH], ModM4@q3];
      \[CapitalGamma]k = Exp[-I 2 Pi kvec.(xyz2RotTSU2[q3][[2]] - xyz2RotTSU2[ModM4@q3][[2]])];
      Irepq = ArrayFlatten[{{0, IdentityMatrix[IrepDim], 0}, {0, 0, IdentityMatrix[IrepDim]}, {\[CapitalGamma]k OrbitIreps[[1]][[pos]], 0, 0}}];
      InducedIreps = Flatten[{IrepGrpH, Table[qh = GTimes[{q[[2]], Keys[grpH][[i]]}];
          \[CapitalGamma]k = Exp[I 2 Pi kvec.(xyz2RotTSU2[qh][[2]] - xyz2RotTSU2[ModM4@qh][[2]])];
          \[CapitalGamma]k Irepq.IrepGrpH[[i]], {i, gH}], Table[qqh = GTimes[{q[[3]], Keys[grpH][[i]]}];
          \[CapitalGamma]k = Exp[I 2 Pi kvec.(xyz2RotTSU2[qqh][[2]] - xyz2RotTSU2[ModM4@qqh][[2]])];
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

PrintIreps[latt_, classes_, k_, ireps_, characters_, reality_] := Module[{i, l, pos, lc, lp, head1, head2, kvec, klabel, tran},
  tran = Thread[{ToExpression["\!\(\*SubscriptBox[\(t\), \(1\)]\)"], ToExpression["\!\(\*SubscriptBox[\(t\),  \(2\)]\)"], ToExpression["\!\(\*SubscriptBox[\(t\), \(3\)]\)"]}->{0,0,0}];
  {kvec, klabel} = k;
  lc = Length[classes];
  lp = Length[ireps];
  head1 = ToString[Superscript[klabel, ToString[#]], StandardForm] <> " (" <> ToString[reality[[#]]] <> ")" & /@ Range[lp];
  head2 = ToString[Length[#]]<>GetEleLabel[#[[1]], latt] &/@ classes;
  head3 = ToString[GetEleLabel[#, latt]<>"\[NewLine]"<>#, StandardForm] &/@ Flatten[classes];
  Print[TableForm[characters\[Transpose]/.tran, TableHeadings -> {head2, head1}, TableAlignments -> Right]];
  Print["---------------------------------------------------------------------------"];
  Print[TableForm[Map[MatrixForm, #] &/@ (ireps\[Transpose]/.tran), TableHeadings -> {head3, head1}, TableAlignments -> Left]];
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

GetSuperCellGrp[t_, OptionsPattern[{"latt"->{}}]] := Module[{xyz, grp, grpd},
  xyz = "x+" <> ToString[#1, StandardForm] <> "," <> "y+" <> ToString[#2, StandardForm] <> "," <> "z+" <> ToString[#3, StandardForm] & @@@ t;
  grp = xyz2Grp[Mat2xyz[xyz2Mat[#]] & /@ xyz];
  grpd = If[OptionValue["latt"]==={}, grp, GetDoubleGrp[grp, OptionValue["latt"]]];
  Return[grpd]
]

GetEulerVector[ele_] := Module[{},
  Which[AssociationQ[ele],
        GetEulerVector[#] &/@ Keys[ele],
        ListQ[ele],
        GetEulerVector[#] &/@ ele,
        MatrixQ[ele],
        Mat2EulerVector[ele[[1;;3,1;;3]]],
        Length@Dimensions[ele] == 3,
        GetEulerVector[#] &/@ ele,
        StringQ[ele],
        GetEulerVector[xyz2Mat[ele][[1]]]
       ]
]

GetCellFromGrp[grp_] := Module[{cell},
  cell = Which[AssociationQ[grp], 
               Length@DeleteDuplicates[#] & /@ ((Rationalize[# - Mod[#, 1]]  & /@ DeleteDuplicates[xyz2RotTSU2[Keys@grp]\[Transpose][[2]]])\[Transpose]),
               ListQ[grp], Length@DeleteDuplicates[#] & /@ ((Rationalize[# - Mod[#, 1]] & /@ DeleteDuplicates[xyz2RotTSU2[Keys@Merge[grp,Identity]]\[Transpose][[2]]])\[Transpose])];
  Return[cell]
]

GetDoubleGrp[grp_, latt_] := Module[{su2, so3, grpd},
  su2 = xyz2SU2[latt, grp];
  so3 = First[Values[grp]\[Transpose]];
  grpd = Association[Join[Mat2xyz[#] -> # & /@ Transpose[{so3, su2}], 
                          Mat2xyz[#] -> # & /@ Transpose[{so3, -su2}]]];
  Return[grpd]
]
(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
