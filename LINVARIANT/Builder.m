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
FitTrainingSet[ham_, ts_, OptionsPattern[{"plot"->True}]] := Module[{gf, t, vars, Coeff, plt},
  gf = Sum[((ham /. t[[2]]) - t[[1]])^2, {t, ts}];
  vars = Variables[gf];
  (*Coeff=FindMinimum[GoalFunction[model,volume,ts,ConstrainCoeff],{coefficients,ConstantArray[0,Length@coefficients]}^\[Transpose],MaxIterations\[Rule]10000];*)
  Coeff = Minimize[gf, vars];
  plt = Show[ListPlot[ts\[Transpose][[1]], 
                      PlotMarkers -> {Automatic, Medium}, 
                      Frame -> True, 
                      GridLines -> Automatic,
                      PlotLabel -> "Fitting error: " <> ToString[Coeff[[1]]]],
         ListLinePlot[Table[ham /. t[[2]] /. Coeff[[2]], {t, ts}], 
                      PlotStyle -> Directive[Red, Dashed], 
                      Frame -> True, 
                      GridLines -> Automatic]];
  If[OptionValue["plot"], Print[plt], Print["Fitting error is " <> ToString[Coeff[[1]]]]];
  Return[Coeff[[2]]]
]

(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
