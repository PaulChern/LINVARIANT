BeginPackage["LINVARIANT`QM`"]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
qmJQ                        ::usage "qmJQ[j]"
qmJ\[Dagger]                ::usage "qmJ\[Dagger][j]"
qmJ                         ::usage "qmJ[j]"
qmJI                        ::usage "qmJI[j]"
qmJz                        ::usage "qmJz[j]"
qmJx                        ::usage "qmJx[j]"
qmJy                        ::usage "qmJy[j]"
qmSU                        ::usage "qmSU[j, n, t]"
(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
qmJQ[j_] := Module[{}, Return[IntegerQ[2 j] && j >= 0]]

qmJ\[Dagger][j_?qmJQ] := Module[{m},
  Normal[SparseArray[Band[{1, 2}] -> Table[ToExpression["\[HBar]"] Sqrt[(j - m) (j + m + 1)], {m, j - 1, -j, -1}], {2 j + 1, 2 j + 1}]]
]
qmJ[j_?qmJQ] := Module[{m},
  Normal[SparseArray[Band[{2, 1}] -> Table[ToExpression["\[HBar]"] Sqrt[(j + m) (j - m + 1)], {m, j, -j + 1, -1}], {2 j + 1, 2 j + 1}]]
]

qmJx[j_?qmJQ] := Module[{},
  1/2 (qmJ\[Dagger][j] + qmJ[j])
]

qmJy[j_?qmJQ] := Module[{},
  1/(2 I) (qmJ\[Dagger][j] - qmJ[j])
]

qmJz[j_?qmJQ] := Module[{m},
  Normal[SparseArray[Band[{1, 1}] -> Table[ToExpression["\[HBar]"] m, {m, j, -j, -1}], {2 j + 1, 2 j + 1}]]
] 

qmJI[j_?qmJQ] := Module[{}, Return[IdentityMatrix[2 j + 1]]]

qmSU[j_?qmJQ, n_, t_] := Module[{},
  MatrixExp[-I t (n/Norm[n]).{qmJx[j], qmJy[j], qmJz[j]}/ToExpression["\[HBar]"]]
]

(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
