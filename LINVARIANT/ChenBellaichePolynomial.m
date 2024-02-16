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
CBInversionQ[H_, vars_] := Module[{v},
  Table[{H === Total[Times @@ {#/First@Cases[#, Exp[__]], First@Cases[#, Exp[__]] /. {v -> -v}} & /@ Level[H, {1}]], v}, {v, vars}]
]
(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
