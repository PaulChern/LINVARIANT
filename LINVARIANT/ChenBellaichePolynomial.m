BeginPackage["LINVARIANT`ChenBellaichePolynomial`", {"LINVARIANT`INVARIANT`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
CBInversionQ              ::usage "CBInversionQ[H, vars]"

(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
CBInversionQ[H_, vars_] := Module[{sub, s, Q},
  sub = Thread[Flatten@vars -> # * Flatten[vars]] &/@ Tuples[{-1, 1}, Length[Flatten@vars]];
  Q = Table[H === Total[Times @@ {#/First@Cases[#, Exp[__]], First@Cases[#, Exp[__]] /. s} & /@ Level[H, {1}]], {s, sub}];
  Return[And@@Q]
]
(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
