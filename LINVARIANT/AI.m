BeginPackage["LINVARIANT`AI`", {"LINVARIANT`INVARIANT`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
BernsteinLayer              ::usage "BernsteinLayer[degree, n, opts]" 
InvariantLayer              ::usage "InvariantLayer[vars, invariants]"

(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
BernsteinLayer[degree_?IntegerQ, n_ : Automatic, opts___Rule] := Module[{func},
  func = Function[{x}, Evaluate[Refine[BernsteinBasis[degree, i, x] // PiecewiseExpand, 0 < x < 1 && 0 <= i <= degree]]];
  NetGraph[Join[Table[NetChain[{ElementwiseLayer[func], 
                                LinearLayer[n, Sequence["Biases" -> None, opts]]}], {i, 0, degree}], 
                {TotalLayer[]}], Table[i -> degree + 2, {i, 1, degree + 1}]]
]

InvariantLayer[vars_, invariants_] := Module[{out, layer, func, i},
  func = Function[Evaluate[ReplaceAll[Thread[Flatten[vars] -> Table[Inactive[Part][#, i], {i, Length[Flatten[vars]]}]]][{Flatten[invariants]}]]];
  layer = FunctionLayer[Activate[func], "Input" -> Length[Flatten@vars]];
  Return[layer]
]
(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
