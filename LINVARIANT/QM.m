BeginPackage["LINVARIANT`QM`", {"LINVARIANT`Structure`","LINVARIANT`GroupTheory`", "LINVARIANT`MathematicaPlus`", "LINVARIANT`INVARIANT`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
QMJQ                        ::usage "QMJQ[j]"
QMJ\[Dagger]                ::usage "QMJ\[Dagger][j]"
QMJ                         ::usage "QMJ[j]"
QMJI                        ::usage "QMJI[j]"
QMJz                        ::usage "QMJz[j]"
QMJx                        ::usage "QMJx[j]"
QMJy                        ::usage "QMJy[j]"
SUMatrix                    ::usage "SUMatrix[j, n, t]"
SpinMatrices                ::usage "SpinMatrices[j]"
AddQNumber                  ::usage "AddOperator[a]"
ShowQNumberList             ::usage "ShowQNumberList[]"
CommutatorDefinition        ::usage "CommutatorDefinition"
AntiCommutatorDefinition    ::usage "AntiCommutatorDefinition"
getOccNumRepresentation     ::usage "getOccNumRepresentation[siteNum, particleNum]"
PauliNormDecompose          ::usage "PauliNormDecompose[M]"
PauliDecompose              ::usage "PauliDecompose[M]"
GetHamJ                     ::usage "GetHamJ[J, latt, spg0, "Spherical"]"
GetHamBasis                 ::usage "GetHamBasis[dim]"
GetInvariantHk              ::usage "GetSymmHk[latt, spg0, Jorb, k, nth, htype]"
HkApplySymm0                ::usage "HkApplySymm0[latt, s, k, hk, AngBasis, HarmonicType]"
HkSpinorApplySymm           ::usage "HkSpinorApplySymm0[latt, s, k, hk, AngBasis, HarmonicType]"
HkApplySymm                 ::usage "HkApplySymm[latt, s, k, hk, AngBasis, HarmonicType]"
TimeReversalU               ::usage "TimeReversalU[j]"

(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)
QNumberList = {};
(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
QMJQ[j_] := Module[{}, Return[IntegerQ[2 j] && j >= 0]]

QMJ\[Dagger][j_?QMJQ] := Module[{m},
  Normal[SparseArray[Band[{1, 2}] -> Table[ToExpression["\[HBar]"] Sqrt[(j - m) (j + m + 1)], {m, j - 1, -j, -1}], {2 j + 1, 2 j + 1}]]
]
QMJ[j_?QMJQ] := Module[{m},
  Normal[SparseArray[Band[{2, 1}] -> Table[ToExpression["\[HBar]"] Sqrt[(j + m) (j - m + 1)], {m, j, -j + 1, -1}], {2 j + 1, 2 j + 1}]]
]

QMJx[j_?QMJQ] := Module[{},
  1/2 (QMJ\[Dagger][j] + QMJ[j])
]

QMJy[j_?QMJQ] := Module[{},
  1/(2 I) (QMJ\[Dagger][j] - QMJ[j])
]

QMJz[j_?QMJQ] := Module[{m},
  Normal[SparseArray[Band[{1, 1}] -> Table[ToExpression["\[HBar]"] m, {m, j, -j, -1}], {2 j + 1, 2 j + 1}]]
] 

QMJI[j_?QMJQ] := Module[{}, Return[IdentityMatrix[2 j + 1]]]

HilbertRotMatrix[j_?QMJQ, n_, t_, e_] := Module[{},
  ComplexExpand[e^j MatrixExp[-I t (n/Norm[n]).{QMJx[j], QMJy[j], QMJz[j]}/ToExpression["\[HBar]"]]]
]

TimeReversalU[j_?QMJQ]:=HilbertRotMatrix[j, {0, 1, 0}, Pi, 1]

SpinMatrices[j_?QMJQ] := Module[{out},
  out = {QMJx[j], QMJy[j], QMJz[j], QMJI[j]};
  Return[out]
]

AddQNumber[x_] := Module[{},
  If[ListQ[x], AddQNumber[#] &/@ x,
  QNumberList = DeleteDuplicates[Append[QNumberList, x]];
  CommutatorDefinition[x[i__], x[j__]] := x[i] ** x[j] - x[j] ** x[i];
  AntiCommutatorDefinition[x[i__], x[j__]] := x[i] ** x[j] + x[j] ** x[i];
  PrintTemporary[Labeled[QNumberList, "activated QNumbers:", Top, LabelStyle -> Bold]];
  Pause[2];
  ];
]

ShowQNumberList[] := Module[{},
  Print[Labeled[QNumberList, "activated QNumbers:", Top, LabelStyle -> Bold]];
]

getOccNumRepresentation[siteNum_, particleNum_] := ReverseSort@Join[Permutations[PadRight[#, siteNum]] & /@ IntegerPartitions[particleNum, siteNum]]

PauliDecompose[M_] := Module[{MI, Mx, My, Mz},
  MI = (M[[1 ;; ;; 2, 1 ;; ;; 2]] + M[[2 ;; ;; 2, 2 ;; ;; 2]])/2;
  Mz = (M[[1 ;; ;; 2, 1 ;; ;; 2]] - M[[2 ;; ;; 2, 2 ;; ;; 2]])/2;
  Mx = (M[[1 ;; ;; 2, 2 ;; ;; 2]] + M[[2 ;; ;; 2, 1 ;; ;; 2]])/2;
  My = (M[[1 ;; ;; 2, 2 ;; ;; 2]] - M[[2 ;; ;; 2, 1 ;; ;; 2]])/(-2 I);
  Return[{Mx, My, Mz, MI}]
]

PauliNormDecompose[M_] := Module[{MI, Mx, My, Mz, evec},
  {Mx, My, Mz, MI} = PauliDecompose[M];
  evec = Normalize[Tr[#] & /@ {Mx, My, Mz}];
  Return[Total[evec {Mx, My, Mz}]]
]

GetHamJ[J_?QMJQ, latt_, spg0_, htype_ : "Spherical"] := Module[{SU, dim, Hij, out, s},
  dim = 2 J + 1;
  Hij = Table[ToExpression@StringRiffle[{"R", i, j}, ""] + I ToExpression@StringRiffle[{"I", i, j}, ""], {i, dim}, {j, dim}];
  Hij = SimplifyTensor[Hij, Hij\[ConjugateTranspose] + Hij];
  out = ComplexExpand@Sum[SU = GetAngularMomentumRep[latt, s, J, htype];
        SU.Hij.(SU\[ConjugateTranspose]), {s, Keys@spg0}];
  Return[SimplifyTensor[Hij, out]]
]

GetHamBasis[dim_] := Module[{null, tmp, i, j, out},
  null = ConstantArray[0, {dim, dim}];
  out = DeleteCases[Flatten[Table[Table[tmp = null;
                    tmp[[i, j]] = 1 + I; 
                    (Sign[#[(tmp\[ConjugateTranspose] + tmp)]] & /@ {Re, Im}) {1, I}, {j, i, dim}], {i, dim}], 2], null];
  Return[out]
]

GetInvariantHk[latt_, spg0_, Jorb_, spinorQ_, k_, nth_, htype_ : "Spherical", OptionsPattern[{"wannier90" -> True}]] := Module[{su2, SU, U, dim, out, s, ksub, HamBasis, hf, HamFamily, p, poly, hb, null, i, m, models, dg, order, TRU, TRQ, sspg, UHU},

  dim = Which[AllTrue[Jorb, IntegerQ], If[spinorQ, 2, 1] Total[2 # + 1 & /@ Jorb], 
              NoneTrue[Jorb, IntegerQ], Total[2 # + 1 & /@ Jorb],
              True, Print["Angular momentum basis should be all integers or all half integers"];Abort[]];

  sspg = If[DoubleGrpQ[#], ##&[], #] &/@ Keys[spg0];

  order = If[ListQ[nth], nth, Range[0, nth]];
  HamBasis = GetHamBasis[dim];
  null = ConstantArray[0, {dim, dim}];
  poly = Table[Flatten[GenMonomialList[k, i]], {i, order}];
  HamFamily = Table[Table[{p, hb}, {p, poly[[i]]}, {hb, HamBasis}], {i, Length@order}]; 
  Print["Number of Hamiltonian families: " <> ToString[Length[Flatten[HamFamily, 2]]]];

  UHU = If[spinorQ, HkSpinorApplySymm, HkApplySymm];

  models = Table[Table[SimplifyTensorCommonFactor[Expand@Sum[UHU[latt, s, k, #, Jorb, htype], {s, sspg}]] & /@ hf, {hf, HamFamily[[i]]}], {i, Length@order}];

  out = Table[DeleteDuplicates[DeleteCases[Flatten[m, 1], null], 
                               ComplexExpand[#1/First@DeleteCases[Flatten[GetConstantFactor[Flatten[#1]]], 0] 
                                           - #2/First@DeleteCases[Flatten[GetConstantFactor[Flatten[#2]]], 0]] == null &], {m, models}];
  Return[out]
]

HkApplySymm0[latt_, s_, k_, hk_, AngBasis_, HarmonicType_, OptionsPattern[{"wannier90" -> True}]] := Module[{dg, TRQ, TRU, su2, SU, U, w90Q, ksub, out, spinorQ},

  If[!AllTrue[AngBasis, IntegerQ] && !NoneTrue[AngBasis, IntegerQ], 
     Print["Angular momentum basis should be all integers or all half integers"];Abort[]];

  spinorQ = Length[hk[[2]]] == 2 Total[2 # + 1 &/@ AngBasis];
  w90Q = OptionValue["wannier90"];
  {dg, TRQ} = xyz2Expr[s][[3;;4]];

  ksub = GetReciprocalTRules[latt, s, k];
  su2 = xyz2su2[latt, s];
  TRU = Which[(!spinorQ)&&TRQ, Fold[ArrayFlatten[{{#1, 0}, {0, #2}}] &, (TimeReversalU[#] & /@ AngBasis)], 
              spinorQ&&TRQ, TimeReversalU[1/2],
              (!spinorQ)&&(!TRQ), IdentityMatrix[Total[2 # + 1 &/@ AngBasis]],
              spinorQ&&(!TRQ), IdentityMatrix[2]];
  U = Fold[ArrayFlatten[{{#1, 0}, {0, #2}}] &, (GetAngularMomentumRep[latt, s, #, HarmonicType] & /@ AngBasis)];
  SU = If[spinorQ, ArrayFlatten[If[w90Q, (TRU.su2)\[TensorProduct]U, U\[TensorProduct](TRU.su2)]], dg TRU.U];

  out = {First[hk] /. ksub, (SU.If[TRQ, Conjugate[hk[[2]]], hk[[2]]].(SU\[ConjugateTranspose]))};

  Return[out]
]

HkSpinorApplySymm[latt_, s_, k_, hk_, AngBasis_, HarmonicType_] := Module[{dg, TRQ, TRU, su2, SU, U, w90Q, ksub, out},

  {dg, TRQ} = xyz2Expr[s][[3;;4]];

  ksub = GetReciprocalTRules[latt, s, k];
  su2 = xyz2su2[latt, s];
  TRU = If[TRQ, TimeReversalU[1/2], IdentityMatrix[2]];
  U = Fold[ArrayFlatten[{{#1, 0}, {0, #2}}] &, (GetAngularMomentumRep[latt, s, #, HarmonicType] & /@ AngBasis)];
  SU = ArrayFlatten[(TRU.su2)\[TensorProduct]U];

  out = (First[hk]/.ksub)*(SU.Conjugate[hk[[2]]].(SU\[ConjugateTranspose]));

  Return[out]
]

HkApplySymm[latt_, s_, k_, hk_, AngBasis_, HarmonicType_] := Module[{dg, TRQ, TRU, su2, SU, U, w90Q, ksub, out},

  {dg, TRQ} = xyz2Expr[s][[3;;4]];

  ksub = GetReciprocalTRules[latt, s, k];
  TRU = If[TRQ, Fold[ArrayFlatten[{{#1, 0}, {0, #2}}] &, (TimeReversalU[#] & /@ AngBasis)], IdentityMatrix[Total[2 # + 1 &/@ AngBasis]]];
  U = Fold[ArrayFlatten[{{#1, 0}, {0, #2}}] &, (GetAngularMomentumRep[latt, s, #, HarmonicType] & /@ AngBasis)];
  SU = dg TRU.U;

  out = (First[hk]/.ksub)*(SU.Conjugate[hk[[2]]].(SU\[ConjugateTranspose]));

  Return[out]
]

(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
