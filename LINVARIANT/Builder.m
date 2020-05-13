BeginPackage["LINVARIANT`Builder`"]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
FitTrainingSet      ::usage "FitTrainingSet[goalfunction]"

(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
FitTrainingSet[gf_] := Module[{coefficients, Coeff},
  coefficients = Variables[gf];
  (*Coeff=FindMinimum[GoalFunction[model,volume,ts,ConstrainCoeff],{coefficients,ConstantArray[0,Length@coefficients]}^\[Transpose],MaxIterations\[Rule]10000];*)
  Coeff = Minimize[gf, coefficients];
  Print["Fitting error is " <> ToString[Coeff[[1]]]];
  Return[Coeff[[2]]]
]

(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
