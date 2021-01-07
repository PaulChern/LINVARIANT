BeginPackage["LINVARIANT`Builder`"]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
FitTrainingSet      ::usage "FitTrainingSet[goalfunction]"
GenRandTrainSet     ::usage "GenRandTrainSet[ref, gs, basis, BasisLabel, modes, NTS, range]"

(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
FitTrainingSet[ham_, ts_, OptionsPattern[{"plot"->True}]] := Module[{gf, t, x, vars, Coeff, plt, ene, nts},
  nts = Length[ts];
  gf = Sum[t[[3]] ((ham /. t[[2]]) - t[[1]])^2, {t, ts}] + Total[(100 D[ham, {#, 1}]^2 & /@ Cases[ts[[1, 2]], (x_ -> v_) -> x]) /. ts[[1, 2]]];
  gf = Sum[Sum[t[[3]] ((ham /. t[[2]]) - t[[1]])^2, {t, ts[[i]]}] + Total[(ts[[i, 1, 3]] D[ham, {#, 1}]^2 & /@ Cases[ts[[i, 1, 2]], (x_ -> v_) -> x]) /. ts[[i, 1, 2]]], {i, nts}];
  vars = Variables[gf];
  (*Coeff=FindMinimum[GoalFunction[model,volume,ts,ConstrainCoeff],{coefficients,ConstantArray[0,Length@coefficients]}^\[Transpose],MaxIterations\[Rule]10000];*)
  Coeff = Minimize[gf, vars];
  ene = ham /. Coeff[[2]] /. #2 & @@@ Flatten[ts, 1];
  plt = Grid[{{Show[ListPlot[Flatten[ts, 1]\[Transpose][[1]], 
                             PlotMarkers -> {Automatic, Medium}, 
                             Frame -> True, 
                             GridLines -> Automatic,
                             PlotLabel -> "Fitting error: " <> ToString[Coeff[[1]]]],
                    ListLinePlot[Table[ham /. t[[2]] /. Coeff[[2]], {t, Flatten[ts, 1]}], 
                                 PlotStyle -> Directive[Red, Dashed], 
                                 Frame -> True, 
                                 GridLines -> Automatic]], 
               ListPlot[{Flatten[ts, 1]\[Transpose][[1]], ene}\[Transpose], 
                        PlotMarkers -> {Automatic, Small},
                        Frame -> True,
                        AspectRatio -> 1,
                        GridLines -> Automatic, 
                        Epilog -> {Red, Line[{Min[Flatten[ts, 1]\[Transpose][[1]]] {1, 1}, Max[Flatten[ts, 1]\[Transpose][[1]]] {1, 1}}], Green, PointSize[Medium], Point[{Flatten[ts, 1]\[Transpose][[1, 1]], ene[[1]]}]}]}}];
  If[OptionValue["plot"], Print[plt], Print["Fitting error is " <> ToString[NumberForm[DecimalForm[N@Coeff[[1]]], {16, 15}]]]];
  Return[Coeff[[2]]]
]

GenRandTrainSet[ref_, gs_, basis_, BasisLabel_, id_, NTS_, range_] := Module[{v, i, j, NModes, dist, phonon, modes, pos, tsrand},
  NModes = Length@id;
  dist = ISODISTORT[gs[[1]], ref[[2]], {PosMatchTo[ref[[2]]\[Transpose][[1]], gs[[2]]\[Transpose][[1]]][[2]] ref[[2]]\[Transpose][[2]]}\[Transpose], basis, BasisLabel\[Transpose][[2]]];
  modes = Table[dist[[#]] & /@ (id[[i]]), {i, NModes}];
  tsrand = Table[v = RandomReal[range[[1]] {-1, 1}, NModes];
                 phonon = Table[{#1, #2, #3, If[#4 == 0, range[[2]] v[[i]], #4  v[[i]] - #4]} & @@@ (modes[[i]]), {i, NModes}];
                 pos = ImposeMode[gs[[1]], gs[[2]], basis, Flatten[phonon, 1], 1];
                 Chop@{gs[[1]], pos}, {NTS-1}];
  Join[{Chop@gs}, tsrand]                 
]
(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
