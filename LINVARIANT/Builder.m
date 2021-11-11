BeginPackage["LINVARIANT`Builder`", {"LINVARIANT`Structure`", "LINVARIANT`Vasp`", "LINVARIANT`INVARIANT`", "LINVARIANT`MathematicaPlus`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
FitTrainingSet      ::usage "FitTrainingSet[goalfunction]"
GenRandTrainSet     ::usage "GenRandTrainSet[ref, gs, basis, BasisLabel, modes, NTS, range]"
GetTrainingSet      ::usage "GetTrainingSet[ref, phases, spg, rho]"
WeightTrainSet      ::usage "WeightTrainSet[data, tau]"
ExportTrainingSet   ::usage "ExportTrainingSet[dir, data, tau]"
BoundHighOrderTerm  ::usage "BoundHighOrderTerm[invariants, expr, vars, cuts]"

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
  dist = ISODISTORT[gs[[1]], ref[[2]], {PosMatchTo[ref[[1]], ref[[2]]\[Transpose][[1]], gs[[2]]\[Transpose][[1]]][[2]], ref[[2]]\[Transpose][[2]]}\[Transpose], basis, BasisLabel\[Transpose][[2]]];
  modes = Table[dist[[#]] & /@ (id[[i]]), {i, NModes}];
  tsrand = Table[v = RandomReal[range[[1]] {-1, 1}, NModes];
                 phonon = Table[{#1, #2, #3, If[#4 == 0, range[[2]] v[[i]], #4  v[[i]] - #4]} & @@@ (modes[[i]]), {i, NModes}];
                 pos = ImposeMode[gs[[1]], gs[[2]], basis, Flatten[phonon, 1], 1];
                 Chop@{gs[[1]], pos}, {NTS-1}];
  Join[{Chop@gs}, tsrand]                 
]

GetTrainingSet[ref_, phases_, spg_, rho_] := Module[{v, i, j, NModes, dist, phonon, modes, pos, tsrand, DomainSet, TrainingSet},
  DomainSet = Flatten[Table[DeleteDuplicates[GetDomains[ref, pos, #] & /@ Keys[spg],StructDist[#1, #2] < 10^-3 &], {pos, phases}], 1];
  TrainingSet = Flatten[StructInterpolation[#1, #2, rho] & @@@ Tuples[DomainSet, {2}],1];
  Return[{DomainSet, TrainingSet}]
]

WeightTrainSet[data_, tau_] := Module[{ts, d},
  ts = Flatten[data, 1];
  Table[Sum[d = StructDist[data[[1, i]], ts[[j]]]; tau/(d^2 + tau^2), {i, Length[data[[1]]]}], {j, Length@ts}]
]

ExportTrainingSet[dir_, data_, tau_] := Module[{ts, weight},
  ts = Flatten[data, 1];
  weight = WeightTrainSet[data, tau];
  Do[ExportPOSCAR[dir <> "/", "POSCAR." <> ToString[i] <> ".vasp", ts[[i]], "head" -> ToString[weight[[i]]]], {i, Length@ts}]
]

BoundHighOrderTerm[invariants_, expr_, vars_, cuts_] := Module[{exprnew, invnew},
  exprnew = SimplifyCommonFactor[NOrderResponse[expr, vars, 1]];
  invnew = SimplifyCommonFactor[NOrderResponse[expr, vars, 1] /. cuts];
  If[! MemberQ[Flatten[invariants], invnew], BoundHighOrderTerm[invariants, exprnew, vars, cuts], exprnew /. cuts]
]

RoundPhase[phase_] := Module[{i},
  Table[If[Chop[phase[[i]]] == 0., 0, 1], {i, Length[phase]}]
]

PlanFitting[invariants_, vars_, spgmat_, ts_, seeds_] := Module[{p, m, t, i, models, character, modes, phases, plan},
  models = Flatten[Values[GroupBy[#, InvariantCharacter[#, vars] &] & /@ Values[GroupBy[invariants, NumPolynomialVar[#] &]]], 1];
  Print["Number of invariants :" <> ToString[Length[Flatten@models]]];
  Print["Groups of models :" <> ToString[Length[models]]];

  plan = Table[character = InvariantCharacter[m[[1]], vars];
               {m, Join[{ts[[1]]}, Table[modes = DeleteDuplicates[RoundPhase[# . (t[[3]] + t[[4]])] & /@ spgmat];
                                         phases = DeleteDuplicates[RoundPhase[# . t[[3]]] & /@ spgmat];
                                         If[MemberQ[modes, character], If[MemberQ[seeds, character], If[MemberQ[phases, character], t, ## &[]], t], ## &[]], {t, ts}]]}, {m, models}];
  Return[plan]
]

LinvariantFit[invariants_, vars_, order_, ts_, spgmat_, ws_ : None, OptionsPattern[{"orderby" -> "AICc", "round" -> 10^-6.0, "offset" -> None}]] := 
 Module[{NumBasis, models, modelorders, OffsetFunc, fits, TopModels, BoundRadius, tsdata, seeds, plan, measure, i, m, p, t, mf, gc, log, out, phaseweight, minimumtest, sol, modelfitted, modeldata, tsinfo},
  NumBasis = Length[Flatten[vars]];
  OffsetFunc = If[OptionValue["offset"] === None, None, (Total[Times @@ OptionValue["offset"]] /. Thread[Flatten[vars] -> {##}] &)];
  measure = Which[OptionValue["orderby"] === "AIcc", 4, OptionValue["orderby"] === "BIC", 5, OptionValue["orderby"] === "AdjustedRSquared", 6, OptionValue["orderby"] === "RSquared", 7, True, 4];
  
  modeldata = If[OptionValue["offset"] === None, Flatten[invariants[[#]] & /@ Range[order]], Delete[Flatten[invariants[[#]] & /@ Range[order]], Position[Flatten[invariants[[#]] & /@ Range[order]], #][[1]] & /@ (OptionValue["offset"][[2]])]];
  seeds = RoundPhase[#] & /@ (ws\[Transpose][[1]]);
  plan = PlanFitting[modeldata, Flatten[vars], spgmat, ts, seeds];
  If[MemberQ[plan\[Transpose][[2]], {}], Print["Error! Traningset is not complete!"]; Abort[]];

  out = If[OptionValue["offset"] === None, {{{}, {}}}, {OptionValue["offset"]}];
  
  TopModels = Table[
    models = p[[1]];
    tsdata = Flatten[{#3 + #4, #2}] & @@@ (p[[2]]);
    tsinfo = GroupBy[{RoundPhase[#3], RoundPhase[#4]} & @@@ (p[[2]]), First -> Last, DeleteDuplicates];
    
    Quiet[phaseweight = If[Length[ws] == 0, Normalize[1.0 & /@ tsdata], Normalize[1.0 + Sum[gc[[3, 2]] Exp[-(Norm[#[[1 ;; -2]] - gc[[1]]]^2/(2. gc[[3, 1]]^2))], {gc, ws}] & /@ tsdata]]];
    
    mf = Quiet[LinearModelFit[tsdata, models, Flatten[vars], IncludeConstantBasis -> False, Weights -> phaseweight, LinearOffsetFunction -> OffsetFunc][{"AICc", "BIC", "AdjustedRSquared", "RSquared", "BestFitParameters", "BestFit"}]];
    
    If[NumPolynomialVar[models[[1]]] == 1, BoundRadius = Max[#] & /@ partitionBy[(Max[#] & /@ Abs[tsdata\[Transpose][[1 ;; NumBasis]]]), Length[#] & /@ vars];
     BoundRadius = Flatten[Thread[ToExpression[GetVariationVar[#1]] -> ConstantArray[#2, Length[#1]]] & @@@ ({vars, BoundRadius}\[Transpose])];
     modelorders = Total@Exponent[#[[1]], Flatten[vars]] & /@ models;
     Do[If[modelorders[[i]] == order && mf[[5, i]] < 0., models[[i]] = BoundHighOrderTerm[invariants, models[[i]], Flatten[vars], BoundRadius]], {i, Length@models}]];
    
    fits = Join[{Length@models + Length[Flatten[out\[Transpose][[1]]]]}, {Join[Flatten[(Times @@ #) & /@ out], models mf[[5]]]}, mf];
    
    AppendTo[out, {Chop[fits[[7]], OptionValue["round"]], models}];
    modelfitted = Total[(Times @@ #) & /@ out, 2];
    OffsetFunc = (modelfitted /. Thread[Flatten[vars] -> {##}] &);
    
    minimumtest = DeleteDuplicates[Flatten[Table[
       sol = Quiet@FindMinimum[modelfitted, {Flatten[vars], m gc[[1]]}\[Transpose]]; 
             If[gc[[3, 2]] != 0, {(Flatten[vars] /. sol[[2]]), {Norm[gc[[1]] - (Flatten[vars] /. sol[[2]])], gc[[2]] - fits[[8]] /. Thread[Flatten[vars] -> gc[[1]]]}}, ## &[]], {gc, ws}, {m, {1., 2.}}], 1], Chop[Norm[#1[[1]] - #2[[1]]], OptionValue["round"]] == 0. &]\[Transpose];
       Join[fits[[;; -3]], {tsinfo}, {MatrixForm[minimumtest[[1]]]}, {MatrixForm[minimumtest[[2]]]}], {p, plan}];
  
  log = Style[# /. x_Real :> Chop[x, OptionValue["round"]]] &@ Grid[{{"Length", "BestFit", "AICc", "BIC", "\!\(\*SuperscriptBox[\(R\), \(2\)]\)", "\!\(\*SuperscriptBox[\(R0\), \(2\)]\)", "trainingset", "Minimums", "distances"}, ## & @@ TopModels}, Dividers -> All];
  Print[log];
  Return[out]
]

(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
