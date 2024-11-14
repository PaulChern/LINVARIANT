BeginPackage["LINVARIANT`Builder`", {"LINVARIANT`Structure`", "LINVARIANT`Vasp`", "LINVARIANT`INVARIANT`", "LINVARIANT`MathematicaPlus`",  "LINVARIANT`Parser`", "LINVARIANT`Fortran`", "LINVARIANT`LatticeHamiltonian`", "LINVARIANT`LatticeHamiltonianFortran`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
GetPhononBasis              ::usage "GetPhononBasis[pos0, Fij]"
ShowPhononBasis             ::usage "ShowPhononBasis[pos0, phonon]"
PhononDecomposition         ::usage "PhononDecomposition[pos0, pos, phonon]"
ExportPhononBasis           ::usage "ExportPhononBasis[dir, cif, pos0, phonon]"
ExportBasisCif              ::usage "ExportBasisCif[dir, cif, pos0, phonon]"
ExportBasis                 ::usage "ExportBasis[dir, phonon]"
PhononLinearRecombination   ::usage "PhononLinearRecombination[PhononBasis, rules]"
FitTrainingSet              ::usage "FitTrainingSet[goalfunction]"
GenRandTrainSet             ::usage "GenRandTrainSet[ref, gs, basis, BasisLabel, modes, NTS, range]"
GetTrainingSet              ::usage "GetTrainingSet[ref, phases, spg, rho]"
WeightTrainSet              ::usage "WeightTrainSet[data, tau]"
ExportTrainingSet           ::usage "ExportTrainingSet[dir, data, tau]"
BoundHighOrderTerm          ::usage "BoundHighOrderTerm[invariants, expr, vars, cuts]"
PlanFitting                 ::usage "PlanFitting[invariants, vars, spgmat, ts, seeds]"
FitGammaModel               ::usage "FitGammaModel[invariants, vars, order, ts, spgmat]"
FitHoppingModel             ::usage "FitHoppingModel[dir, pos0, spg0, tgrp, TRules, Basis, nn, qgrid, vars, svars, varmap]"
DomainGrid                  ::usage "DomainGrid[domains, vars]"
FitSHoppingModel            ::usage "FitSHoppingModel[dir0, pos0, spg0, tgrp, TRules, Basis, nn, qgrid, vars, svars, ie, npt, maxeps]"
FetchSeeds                  ::usage "FetchSeeds[inv, vars, n]"
SeedsByInvariants           ::usage "SeedsByInvariants[invariants, vars]"
BasisComplement             ::usage "BasisComplement[BasisMatrix, FullBasis]"
ImportDMDPath               ::usage "ImportDMDPath[file, pos0, Basis, Basislabel]"
SamplePhaseSpace            ::usage "SamplePhaseSpace[dir0, s0, spgmat, Basis, potim]"
CollectSeeds                ::usage "CollectSeeds[dir0, nphase, Basis, OpMat]"
CollectPhases               ::usage "CollectPhases[dir0, n, vars, Basis]"
ShowSeedPhases              ::usage "ShowSeedPhases[phases, seeds, vars]"
SamplePolynomialSpace       ::usage "SamplePolynomialSpace[polyphase, dir, spgmat0, Basis, npt, smearing, mim0]"
PlanSampling                ::usage "PlanSampling[dir0_, phases, invariants, vars, Basis, spgmat, npt, smearing]"
SamplingByDMDPath           ::usage "SamplingByDMDPath[dir0, DMDPath, vars, svars, BasisMatrix, FullBasis]"
SampleAround                ::usage "SampleAround[pos0, dir, fname, spgmat0, Basis, npt, potim]"
ImposePhase                 ::usage "ImposePhase[pos0, dir, fname, phase, Basis]"
RoundPhase                  ::usage "RoundPhase[phase]"
LoadTrainingset             ::usage "LoadTrainingset[dir0, BasisMatrix, FullBasis]"
PreparePhases               ::usage " PreparePhases[dir0, pname, num]"
GetHam                      ::usage "GetHam[models]"
UpdateModel                 ::usage "UpdateModel[models, ham, vars, svars]"
TrimmedCoefficient          ::usage "TrimmedCoefficient[inv, vars, coeff, trim]"
CheckMinimums               ::usage "CheckMinimums[models, seeds, amp, vars, svars]"
DMDPath2ts                  ::usage "DMDPath2ts[dmd, vars, svars]"
Expression2Model            ::usage "Expression2Model[ham, vars, svars]"
GetActivePolynomial         ::usage "GetActivePolynomial[modeldata, vars, svars, mim]"
SortModelData               ::usage "SortModelData[modeldata, vars, svars]"
JoinModels                  ::usage "JoinModels[m0, m1]"
CheckFitting                ::usage "CheckFitting[invariants, model, ts, vars, svars]"
GetTopOrder                 ::usage "GetTopOrder[polynomials, vars, svars]"
TopOrderQ                   ::usage "TopOrderQ[p, polynomials, vars, svars]"
Model2List                  ::usage "Model2List[model, vars, svars]"
ExplosionQ                  ::usage "ExplosionQ[model, vars, svars, inv]"
EnergyGainQ                 ::usage "EnergyGainQ[model, vars, svars, inv]"
GetPolynomialTrainingSet    ::usage "GetPolynomialTrainingSet[polynomials, tsdata]"
AdjustModel                 ::usage "AdjustModel[model, modeladjust]"
GetBoundingInvariants       ::usage "GetBoundingInvariants[invariant, vars, svars, cut, TRules, lvar, mvar]"
GetOptModel                 ::usage "GetOptModel[model, vars, svars, tsused, cut, OpMat, TRules, lvar, mvar]"
ExtractMinimumTS            ::usage "ExtractMinimumTS[tsdata, OpMat]"
Model2Association           ::usage "Model2Association[model]"
GetFixMimModel              ::usage "GetFixMimModel[vars, order, tsdata, OpMat, TRules, lvar, mvar]"
ModelOptimization           ::usage "ModelOptimization[model, polynomials, ts, vars, svars, OpMat, mim]"
QuickFit                    ::usage "QuickFit[poly, tsdata, vars, svars, prefix]"
FitHessian                  ::usage "FitHessian[Basis, LTB, HJij, pos, vars]"
GetAcousticEnergy           ::usage "GetAcousticEnergy[spg0, tgrp, LTB, pos0, avars]"
GetAcousticEps              ::usage "GetAcousticEps[spg0, tgrp, LTB, pos0, avars]"
BuildHeterostructureStrain  ::usage "BuildHeterostructureStrain[dir]"
ExportHamiltonian           ::usage "ExportHamiltonian[dir, fullmodel, vars, svars]"
ExportIO                    ::usage "ExportIO[dir]"
BuildLINVARIANT             ::usage "BuildLINVARIANT[dir0, fullmodel, vars, svars, EwaldFieldQ]"
ImposeAcousticSumRule       ::usage "ImposeAcousticSumRule[Hacoustic]"
UpdatePT                    ::usage "UpdatePT[dir, file, WriteQ]"
SqueezeCell                 ::usage "SqueezeCell[dir0, pos0, eta, nstrain, qgrid]"
FromJij2Energy              ::usage "FromJij2Energy[model, grid]"
FitSupercellHessian         ::usage "FitSupercellHessian[Basis, Fij, qgrid, HJij, pos0, vars]"
CookModel                   ::usage "CookModel[fullmodel, vars]"
OriginShiftFC               ::usage "OriginShiftFC[pos0, Fij0, NSC, shift]"
TrainingSet2Supercell       ::usage "TrainingSet2Supercell[ts, grid, vars, uvars]"
FoldGamma                   ::usage "FoldGamma[gmodel, ah]"
LMesh                       ::usage "LMesh[dir, grid, mfunc]"
TrainLINVARIANT             ::usage "TrainLINVARIANT[Plan, vars, ts]"
UpdateADAM                  ::usage "UpdateADAM[Param, mParam, vParam, grad, iter, \[Rho], \[Beta], \[Epsilon]]"
GaussianPolynomial          ::usage "GaussianPolynomial[order, var]"
GetLossFunc                 ::usage "GetLossFunc[H, vars, ts]"
BoundData                   ::usage "BoundData[data, npts, dist]"
ChopGaussianPolynomial      ::usage "ChopGaussianPolynomial[expr, vars, ts, IntRange, tol]"

(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
GetPhononBasis[pos0_, Fij_, OptionsPattern[{"table" -> True, "fontsize" -> 12, "toiso"->None, "IRTol"->0.01, "roundup"->10^-6}]] := Module[{sol, PhononByIR, PhononRotated, ind, frequency, IR, i, j, iv, eigenvectors, rot, NumAtom, tab, simplifytab, color, isolabels, Vasp2THz = 15.633302300230191, latt, sites, alatt},
  {latt, sites} = pos0;
  alatt = Inverse[Transpose@latt];
  NumAtom = Length[sites];
  simplifytab = {0 -> "-"};
  sol = Eigensystem[ArrayFlatten[Fij]]\[Transpose];
  PhononByIR = Gather[MapIndexed[Join[#2, {Chop[#1[[1]], OptionValue["IRTol"]], Chop[#1[[2]], OptionValue["roundup"]]}] &, sol], Chop[#1[[2]] - #2[[2]], 10^-3] == 0 &];
  isolabels = If[OptionValue["toiso"]===None, 
                 MapIndexed[ConstantArray["IR"<>ToString[First@#2], #1]&, Length[#]&/@PhononByIR],
                 CloneReshape2D[PhononByIR, OptionValue["toiso"]]];

  PhononRotated = Table[{ind, frequency, IR} = PhononByIR[[i]]\[Transpose];
                        eigenvectors = DeleteCases[Orthogonalize[Flatten[Table[rot = PseudoInverse[IR\[Transpose][[3 (j - 1) + 1 ;; 3 (j - 1) + 3]]\[Transpose]]; Normalize[#] & /@ (rot . IR), {j, NumAtom}], 1], Tolerance -> 10.^-9], v_ /; Norm[v] == 0]; 
                        Table[{ind[[iv]], isolabels[[i,iv]], frequency[[iv]], 
                               Flatten[Chop[alatt.#, OptionValue["roundup"]] &/@ Partition[eigenvectors[[iv]], 3]]}, {iv, Length[ind]}], {i, Length[PhononByIR]}];
  
  tab = Grid[Join[{Join[{"#","IR","frequency (THz)"}, Flatten[Table[Subscript[#, \[Alpha]], {\[Alpha], {"x", "y", "z"}}] & /@ (sites\[Transpose][[2]])]]}, Flatten[{#[[1;;3]], #[[4]]/.simplifytab}] & /@ Flatten[PhononRotated, 1]], 
             Alignment -> Center, 
             Frame -> True, 
             Background -> {Flatten[{Blue, Blue, Blue, Table[color = If[OddQ[i], White, None]; 
                                                       ConstantArray[color, 3], {i, NumAtom}]}], 
                            Flatten[{LightBlue, Table[color = If[OddQ[i], Yellow, Pink]; 
                                                      ConstantArray[color, Length[PhononRotated[[i]]]], {i, Length@PhononRotated}]}]}, 
             ItemSize -> Full, 
             ItemStyle -> Directive[FontSize -> OptionValue["fontsize"], Black]];
  If[OptionValue["table"], Print[tab]];
  Return[PhononRotated]
]

ShowPhononBasis[pos0_, phonon_, OptionsPattern[{"fontsize" -> 12}]] := Module[{sol, PhononByIR, PhononRotated, ind, frequency, IR, i, j, iv, eigenvectors, rot, NumAtom, tab, simplifytab, color},
   NumAtom = Length[pos0[[2]]];
   simplifytab = {0 -> "-"};
   tab = Grid[Join[{Join[{"#","IR","frequency (THz)"}, Flatten[Table[Subscript[#, \[Alpha]], {\[Alpha], {"x", "y", "z"}}] & /@ (pos0[[2]]\[Transpose][[2]])]]}, 
              {Flatten[{{"", "", "sites"}, pos0[[2]]\[Transpose][[1]]}]},
              Flatten[{#[[1;;3]], #[[4]]/.simplifytab}] & /@ Flatten[phonon, 1]],
              Alignment -> Center,
              Frame -> True,
              Background -> {Flatten[{Blue, Blue, Blue, Table[color = If[OddQ[i], White, None];
                                                        ConstantArray[color, 3], {i, NumAtom}]}],
                             Flatten[{LightBlue, LightBlue, Table[color = If[OddQ[i], Yellow, Pink];
                                                       ConstantArray[color, Length[phonon[[i]]]], {i, Length@phonon}]}]},
              ItemSize -> Full,
              ItemStyle -> Directive[FontSize -> OptionValue["fontsize"], Black]];
   Print[NumberForm[tab, {4,4}]];
]


ExportPhononBasis[dir_, cif_, sites0_, phonon_] := Module[{i, phononordered, cifdata, latt, sites},
  cifdata = ImportIsodistortCIF[cif];
  latt = N[cifdata[[6]] /. cifdata[[8]]];
  sites = cifdata[[10]];
  phononordered = SortIR2POSCAR[latt, sites, sites0, phonon\[Transpose]]\[Transpose];

  Table[Export[dir <> "/phonon-" <> ToString[i] <> ".dat", cif2mcif[cif, i, phononordered, sites]];
        Run["mv " <> dir <> "/phonon-" <> ToString[i] <> ".dat " <> dir <> "/phonon-" <> ToString[i] <> ".mcif"], {i, Length@Transpose[phonon]}];
]

ExportBasis[dir_, basis_] := Module[{data},
  data = Flatten[{StringPadLeft[ToString@DecimalForm[#1, {16, 10}, NumberPadding -> {" ", " "}], 16, " "],
                  StringPadLeft[StringRiffle[Level[ToExpression@#2, Infinity], ""], 16, " "], 
                  ToString@DecimalForm[N@#3, {16, 10}, NumberSigns -> {"-", " "}, NumberPadding -> {" ", "0"}], 
                  ToString[DecimalForm[N@#, {16, 10}, NumberSigns -> {"-", " "}, NumberPadding -> {" ", "0"}]] & /@ #4}] & @@@ Flatten[PhononBasis, 1];
  Export[dir <> "/basis.inp", data\[Transpose], "Table"];
]
  
ExportBasisCif[dir_, cif_, sites0_, basis_] := Module[{i, basisordered, cifdata, latt, sites},
  cifdata = ImportIsodistortCIF[cif];
  latt = N[cifdata[[6]] /. cifdata[[8]]];
  sites = cifdata[[10]];
  basisordered = SortIR2POSCAR[latt, sites, sites0, basis\[Transpose]]\[Transpose];

  Table[Export[dir <> "/basis-" <> ToString[i] <> ".dat", cif2mcif[cif, i, N@basisordered, sites]];
        Run["mv " <> dir <> "/basis-" <> ToString[i] <> ".dat " <> dir <> "/basis-" <> ToString[i] <> ".mcif"], {i, Length@Transpose[basis]}];
]

PhononDecomposition[pos0_, pos_, phonon_] := Module[{i, labels, modes, frequencies, color},
  labels = Flatten[phonon, 1]\[Transpose][[2]];
  frequencies = Flatten[phonon, 1]\[Transpose][[3]];
  modes = ISODISTORT[pos0[[1]], pos0[[2]], pos[[2]], Flatten[phonon, 1]\[Transpose][[4]]\[Transpose], labels];
  Grid[Prepend[Insert[modes\[Transpose], frequencies, 3]\[Transpose], {"#", "IR", "frequency", "norm", "amplitude"}], 
       Background -> {None, 
                      Flatten[{LightBlue, Table[color = If[OddQ[i], Yellow, Pink]; 
                              ConstantArray[color, Length[phonon[[i]]]], {i, Length@phonon}]}]}, 
       ItemSize -> Full]
]

PhononLinearRecombination[PhononBasis_, rules_] := Module[{rule, rot, m, s, pos, id, phonon, i, ind, ind0, IR, frequency, eigenvectors},
  ind = DeleteCases[Sort[Flatten[#\[Transpose][[1;;2]]]], 0] & /@ rules;
  Table[ind0 = PhononBasis[[i]]\[Transpose][[1]];
        IR = PhononBasis[[i]]\[Transpose][[2]];
        frequency = PhononBasis[[i]]\[Transpose][[3]];
        pos = Position[ind, ind0];
        rot = If[pos != {},
                 rule = Extract[rules, First[pos]];
                 1/Sqrt[2] Flatten[Table[Which[Length[id] == 1, m = ConstantArray[0, Length[ind0]]; m[[id[[1]]-Min[ind0]+1]] = Sqrt[2]; {m},
                                               Length[id] == 3, RotateLeft[Table[m = ConstantArray[0, Length[ind0]];
                                                                                 m[[id[[1]]-Min[ind0]+1]] = 1; m[[id[[2]]-Min[ind0]+1]] = (-1)^(s - 1);
                                                                                 m, {s, 2}], id[[3]]]], {id, rule}], 1], 
                 IdentityMatrix[Length[ind0]]];
        eigenvectors = rot . (PhononBasis[[i]]\[Transpose][[4]]);
        {ind0, IR, frequency, Chop@eigenvectors}\[Transpose], {i, Length[PhononBasis]}]
]

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
  exprnew = SimplifyCommonFactor[NOrderResponse[expr, vars, 1], 10^-6];
  invnew = SimplifyCommonFactor[NOrderResponse[expr, vars, 1] /. cuts, 10^-6];
  If[DisjointQ[Flatten[invariants], {invnew,-invnew}], BoundHighOrderTerm[invariants, exprnew, vars, cuts], exprnew /. cuts]
]

RoundPhase[phase_] := Module[{i},
  Table[If[Chop[phase[[i]]] == 0., 0, 1], {i, Length[phase]}]
]

SampleAround[dir0_, its_, pos_, TRules_, Basis_, npt_, potim_, OptionsPattern[{"round" -> 10.^-3}]] := Module[{BasisDim, NumBasis, sets, samples, modes, pos0, pos1, phase, t, m, i, spgmat0, spgmat, vars, eij},
  {BasisDim, NumBasis} = Dimensions[Basis];
  vars = #1 & @@@ (TRules[[1]]);
  spgmat0 = Table[Coefficient[# /. t, vars] & /@ vars, {t, TRules}];
  pos0 = ImportPOSCAR[dir0 <> "/POSCAR0"];
  phase = Join[ISODISTORT[pos0[[1]], pos0[[2]], pos[[2]], Basis, Range[NumBasis], "round" -> OptionValue["round"]]\[Transpose][[4]],
               eij2eta[Chop@GetStrainTensor[pos0[[1]], pos[[1]], "iso" -> False]]];
  spgmat = Table[If[Chop[Norm[m . Sign[phase] - Sign[phase]], OptionValue["round"]] == 0, m, ## &[]], {m, spgmat0}];
  
  sets = DeleteDuplicates[Flatten[Table[SortBy[# Sign[Abs[potim]] & /@ Tuples[{0, -i, i}, {NumBasis + 6}], Norm[N@#] &], {i, Range[0, npt]}], 1]];
  
  samples = Reverse@Round[DeleteDuplicates[Table[{i, #.(sets[[i]]) & /@ spgmat}, {i, Length[sets]}], Complement[#1[[2]], #2[[2]]] == {} &]\[Transpose][[2]]\[Transpose][[1]]];
  
  Print["Number of samples: " <> ToString[Length[samples]]];
  
  Do[modes = Table[{m, "modes", 1, potim[[m]] samples[[i, m]]}, {m, NumBasis + 6}];
     eij = eta2eij[potim[[NumBasis + 1 ;;]] samples[[i, NumBasis + 1 ;;]]];
     pos1 = {(eij + IdentityMatrix[3]) . (pos0[[1]]), ImposeMode[pos0[[1]], pos[[2]], Basis, modes[[1 ;; NumBasis]], 1.0]};
     ExportPOSCAR[dir0, "/trainingset/ts-" <> ToString[its] <> "/sample." <> ToString[i] <> ".vasp", pos1], {i, Length@samples}]
]

ImposePhase[pos0_, dir_, fname_, phase_, Basis_] := Module[{m, modes, BasisDim, NumBasis, pos},
  {BasisDim, NumBasis} = Dimensions[Basis];
  modes = Table[{m, "modes", 1, phase[[m]]}, {m, NumBasis}];
  pos = {pos0[[1]], ImposeMode[pos0[[1]], pos0[[2]], Basis, modes, 1.0]}; 
  ExportPOSCAR[dir, fname, pos]
]

FetchSeeds[inv_, vars_, n_, OptionsPattern[{"zero" -> False, "log"->0, "limit"->Infinity}]] := Module[{character, ndigits, shortseed, list, map, seeds0, order}, 
  character = First@InvariantCharacter[inv, vars];
  list = If[OptionValue["log"]==0,
            If[OptionValue["zero"], Range[n,-n,-Sign[n]], DeleteCases[Range[n,-n,-Sign[n]],0]],
            Join[-Reverse[OptionValue["log"]^Range[1,n-1]], OptionValue["log"]^Range[1,n-1]]];
  order = Length[list];
  ndigits = Count[character, 1];
  shortseed = If[order^ndigits < OptionValue["limit"], Tuples[list, ndigits]\[Transpose], DeleteDuplicates@Table[If[# >= 0, 1, -1] & /@ RandomReal[{-1, 1}, ndigits], OptionValue["limit"]]\[Transpose]];
  map = Flatten[Position[character, 1]];
  seeds0 = Table[If[i == 0, ConstantArray[0, Dimensions[shortseed][[2]]], ## &[]], {i, character}];
  Do[seeds0 = Insert[seeds0, shortseed[[i]], map[[i]]], {i, Length@map}];
  Return[SparseArray[seeds0\[Transpose]]]
]

SeedsByInvariants[invariants_, vars_, svars_, OptionsPattern[{"table"->True, "exseeds"->None, "odd"->False, "exclude"->{}}]] := Module[{id, inv, sub, i, seeds, s, dict, sampleseeds, tab, tabdata, color, fullvars, simplifytab, varlabels, exseeds, invdict, models, t, spinv, pinv, sinv},
  fullvars = Append[vars, svars];

  sinv = Select[Flatten[invariants], And @@ XInvQ[#, {svars[[1, 1]]}, ContainsExactly] &];
  pinv = Select[Flatten[invariants], And @@ XInvQ[#, {svars[[1, 1]]}, ContainsNone] &];
  spinv = Complement[Flatten[invariants], Join[sinv, pinv]];
  
  models = Flatten[SortInvariants[#, Flatten@fullvars] & /@ {pinv, spinv, sinv}];

  varlabels = MapIndexed["(" <> ToString[First@#2] <> ")" &, Flatten@fullvars];
  simplifytab = {-1 -> Style["-1", Bold, Blue], 1 -> Style["1", Bold, Red], 0 -> "-"};

  exseeds = If[OptionValue["exseeds"]===None, 
               {},
               {0,"f["<>StringJoin[Riffle[ToString[#, StandardForm] & /@ DeleteCases[Flatten@fullvars Abs[#], 0], ","]] <> "]", #} &/@ OptionValue["exseeds"]];

  seeds = Prepend[Flatten[Table[inv = If[Head[models[[id]]] === Plus, First[models[[id]]], models[[id]]];
                        s = FetchSeeds[models[[id]], Flatten@fullvars, 1, "zero"->False];
                        dict = Values@GroupBy[Table[sub = Thread[Flatten@fullvars -> If[OptionValue["odd"], s[[i]], Abs@s[[i]]]]; 
                                                    {inv /. sub, s[[i]]}, {i, Length[s]}], First];
                        {id, #1 inv, #2} & @@@ First[dict\[Transpose]], {id, Length@models}], 1],
                  {0, "f[0]", ConstantArray[0,Length[Flatten@fullvars]]}];
  sampleseeds = Join[DeleteDuplicates[seeds, #1[[3]]==#2[[3]]&], exseeds];
  sampleseeds = If[OptionValue["exclude"]==={}, 
                   sampleseeds, 
                   If[MemberQ[OptionValue["exclude"]\[Transpose][[3]], #[[3]]], ##&[], #] &/@ sampleseeds];
  tabdata = {#1, #2, #3 /. simplifytab} & @@@ sampleseeds;
  tab = Grid[Prepend[Prepend[MapIndexed[Flatten[{#2, #1}] &, tabdata], Flatten[{{"", "", ""}, Flatten@fullvars}]], Flatten[{{"", "", ""}, varlabels}]],
             Background -> {Flatten[{LightBlue, White,LightBlue, 
                            Table[color = If[OddQ[i], Yellow, Pink];
                                  ConstantArray[color, Length[fullvars[[i]]]], {i, Length@fullvars}]}], None}];
  If[OptionValue["table"], Print[tab]];
  Return[sampleseeds]
]

BasisComplement[BasisMatrix_, FullBasis_] := Module[{ind, FullBasisMatrix, ComplementBasisMatrix, exlabels, IRs},
  FullBasisMatrix = Flatten[FullBasis, 1]\[Transpose][[4]]\[Transpose];
  ComplementBasisMatrix = Complement[FullBasisMatrix\[Transpose], BasisMatrix\[Transpose], SameTest -> (Chop[Norm[#1 - #2]] == 0 &)];
  ind = Sort[First@First[Position[FullBasisMatrix\[Transpose], #]] & /@ ComplementBasisMatrix];
  ComplementBasisMatrix = (FullBasisMatrix\[Transpose][[#]] & /@ ind)\[Transpose];
  IRs = First@First@Position[FullBasis, #] & /@ ind;
  exlabels = "("<>ToString[#]<>") "<>(Flatten[FullBasis, 1]\[Transpose][[2]][[#]]) &/@ ind;
  exlabels = SplitBy[exlabels, Extract[IRs, First@Position[exlabels, #]] &];
  Return[{ComplementBasisMatrix, exlabels}]
]

ImportDMDPath[dir0_, vars_, svars_, uvars_, BasisMatrix_, FullBasis_, seeds_, OptionsPattern[{"fontsize"->12, "ts"->{}, "table"->True, "pdf"->True, "plot"->True, "interpolate"->0, "chop" -> 0.01}]] := Module[{run, vasprun, lattices, latt0, sites0, sites, forces, stress, ene0, ene, coordinates, labels, color, iscf, j, i, path, epsplot, scfplot, roseplot, out, eta, fullvars, OtherBasisMatrix, ind, exlabels0, exlabels, FullBasisMatrix, ExInfo, plots, tsdata, data, ts, ibrion, inv, plt, pos, pos0, nvars, grid, ncells, varsub},
  pos0 = ImportPOSCAR[dir0 <> "/POSCAR0"];
  {latt0, sites0} = pos0;
  FullBasisMatrix = Flatten[FullBasis, 1]\[Transpose][[4]]\[Transpose];
  OtherBasisMatrix = Complement[FullBasisMatrix\[Transpose], BasisMatrix\[Transpose], SameTest -> (Chop[Norm[#1 - #2]] == 0 &)];
  ind = Sort[First@First[Position[FullBasisMatrix\[Transpose], #]] & /@ OtherBasisMatrix];
  ExInfo = Tally[First@First@Position[FullBasis, #] & /@ ind];
  OtherBasisMatrix = (FullBasisMatrix\[Transpose][[#]] & /@ ind)\[Transpose];
  exlabels0 = "("<>ToString[#]<>") "<>(Flatten[FullBasis, 1]\[Transpose][[2]][[#]]) &/@ ind;
  nvars = Length[Flatten[Append[vars, uvars]]];

  data = Import[dir0 <> "/seeds/" <> "ene.dat", "Table"];
  ene0 = data[[1,2]];
  tsdata = If[OptionValue["ts"] === {}, data, If[AllTrue[OptionValue["ts"], MemberQ[data[[;;,1]], #]&], 
                                                 Extract[data, First@Position[data\[Transpose][[1]], #]] &/@ OptionValue["ts"],
                                                 Print["You many need to collect the seeds again before importing the seeds data!"];
                                                 Abort[]]];
  plots = {};

  out = Table[
        ts = tsdata[[i, 1]];

        inv = If[ts<=Length[seeds], ToString[seeds[[ts,2]], StandardForm], "unknown"];
        vasprun = Table[Import[dir0<>"/seeds/seed"<>ToString[ts]<>"/vasprun_ibrion"<>ToString[ibrion]<>".xml"], {ibrion, {3, 1}}];
        sites = Join[ParseXMLData[ParseXML[ParseXML[vasprun[[1]], "structure", {}], "varray", {"name", "positions"}], "v"],
                     {Last@ParseXMLData[ParseXML[ParseXML[vasprun[[2]], "structure", {}], "varray", {"name", "positions"}], "v"]}];
        forces = Join[ParseXMLData[ParseXML[vasprun[[1]], "varray", {"name", "forces"}], "v"],
                     {Last@ParseXMLData[ParseXML[vasprun[[2]], "varray", {"name", "forces"}], "v"]}];
        stress = Join[ParseXMLData[ParseXML[vasprun[[1]], "varray", {"name", "stress"}], "v"],
                      {Last@ParseXMLData[ParseXML[vasprun[[2]], "varray", {"name", "stress"}], "v"]}];
        lattices = Join[ParseXMLData[ParseXML[ParseXML[vasprun[[1]], "structure", {}], "varray", {"name", "basis"}], "v"],
                        {Last@ParseXMLData[ParseXML[ParseXML[vasprun[[2]], "structure", {}], "varray", {"name", "basis"}], "v"]}];
        ncells = Times@@Round[(Norm[#] &/@ lattices[[1]])/(Norm[#] &/@ latt0)];
        eta = Round[eij2eta[GetStrainTensor[latt0, #, "iso" -> False]], 10.0^-4]&/@lattices;
        ene = Chop[ParseFortranNumber@Join[Flatten[ParseXML[Last[ParseXML[#, "energy", {}]], "i", {"name", "e_fr_energy"}] & /@ ParseXML[vasprun[[1]], "calculation", {}]], {Last@Flatten[ParseXML[Last[ParseXML[#, "energy", {}]], "i", {"name", "e_fr_energy"}] & /@ ParseXML[vasprun[[2]], "calculation", {}]]}] - ene0, 10^-6];
        coordinates = Table[pos = {lattices[[iscf]], {sites[[iscf]], Flatten[ConstantArray[#, ncells] &/@ (Transpose[sites0][[2]])]}\[Transpose]};
                            Flatten@Join[ModeDecomposition[pos0,pos,BasisMatrix,vars,svars,uvars, 
                                                           "lattnorm"->False,"table"->False], 
                                         ModeDecomposition[pos0,pos,OtherBasisMatrix,exlabels0,svars,{}, 
                                                           "lattnorm"->False,"table"->False][[1]]], {iscf, Length[sites]}];
        
        varsub = First[coordinates];
        grid = VarSubGrid[varsub];
        ncells = Times@@grid;
        labels = MapIndexed["("<>ToString[First@#2]<>") "<>ToString[VarBare@#1,StandardForm]&, Sub2List[varsub][[1]][[1;;nvars*ncells+6]]];
        exlabels = Flatten@ConstantArray[exlabels0, ncells];
        fullvars = Append[Flatten[ConstantArray[Append[vars, uvars], ncells], 1], svars];

        coordinates = Transpose[Join[{Flatten[{Chop[First[coordinates][[1;;nvars*ncells]], OptionValue["chop"]],
                                               Chop[First[coordinates][[nvars*ncells+1;;nvars*ncells+6]], OptionValue["chop"]/10],
                                               Chop[First[coordinates][[nvars*ncells+7;;]], 10.0^-3]}]},
                                     coordinates[[2;;-2]],
                                     {Flatten[{Chop[Last[coordinates][[1;;nvars*ncells]], OptionValue["chop"]/20],
                                               Chop[Last[coordinates][[nvars*ncells+1;;nvars*ncells+6]], 10.0^-3],
                                               Chop[Last[coordinates][[nvars*ncells+7;;]], 10.0^-3]}]}]];

        path = Grid[Prepend[Prepend[Prepend[Round[Sub2List[#][[2]]&/@coordinates, 10.0^-3], ene], Range[Length[ene]]]\[Transpose], Join[{"#", "energy"}, labels, exlabels]],
                    Background -> {Join[{Gray,LightBlue}, Flatten[Table[color = If[OddQ[j], Yellow, Pink];
                                   ConstantArray[color, Length[fullvars[[j]]]], {j, Length@fullvars}]], 
                                   Flatten[Table[ConstantArray[If[OddQ[j], Gray, LightGray], ExInfo[[Mod[j,Length[ExInfo],1],2]]], {j,ncells*Length[ExInfo]}]]],
                                   None},
                    ItemStyle -> Directive[FontSize -> OptionValue["fontsize"]],
                    ItemSize -> Full
                   ]; 
        epsplot = Partition[Table[If[Norm[Sub2List[coordinates[[j]]][[2]]] != 0, ListPlot[{Sub2List[coordinates[[j]]][[2]], ene}\[Transpose], Frame -> True, GridLines -> Automatic, PlotRange -> All, PlotLabel -> labels[[j]], ImageSize -> 300], ## &[]], {j, Length[labels]}], UpTo[4]];
        scfplot = {ListLinePlot[Sub2List[#][[2]] &/@ (coordinates[[1;;nvars*ncells]]), PlotRange -> All, Frame -> True, PlotMarkers -> Automatic, 
                                                   GridLines -> Automatic, ImageSize -> 300, ImagePadding -> {{50, 1}, {20, 1}}, 
                                                   PlotLegends -> LineLegend[labels[[;;Length[Flatten@vars]]], LegendLayout -> {"Column", 3}]], 
                         ListLinePlot[Sub2List[#][[2]] &/@ (coordinates[[nvars*ncells+1;;nvars*ncells+6]]), PlotRange -> All, Frame -> True, PlotMarkers -> Automatic,
                                                    GridLines -> Automatic, ImageSize -> 300, ImagePadding -> {{50, 1}, {20, 1}},
                                                    PlotLegends -> LineLegend[labels[[Length[Flatten@vars]+1;;]], LegendLayout -> {"Column", 2}]],
                          ListLinePlot[Sub2List[#][[2]] &/@ (coordinates[[nvars*ncells+7;;]]), PlotRange -> All, Frame -> True, PlotMarkers -> Automatic,
                                                     GridLines -> Automatic, ImageSize -> 300, ImagePadding -> {{50, 1}, {20, 1}},
                                                     PlotLegends -> LineLegend[exlabels, LegendLayout -> {"Column", 3}]],
                         ListLinePlot[ene, PlotRange -> All, Frame -> True, PlotMarkers -> Automatic, GridLines -> Automatic, 
                                           ImageSize -> 300, ImagePadding -> {{50, 1}, {20, 1}}]};
        roseplot = Partition[Rasterize[Show[ListPlot3D[{Sub2List[coordinates[[#1]]][[2]], Sub2List[coordinates[[#2]]][[2]], ene}\[Transpose], 
                                                       AxesLabel -> {labels[[#1]], labels[[#2]], ""}, 
                                                       PlotRange -> All, 
                                                       ColorFunction -> "Rainbow", 
                                                       InterpolationOrder -> OptionValue["interpolate"], 
                                                       MeshFunctions -> {#3 &}, 
                                                       Mesh -> Length[ene], 
                                                       ViewPoint -> {0, 0, 100}, 
                                                       ImageSize -> 300], 
                                            ListPointPlot3D[{Sub2List[coordinates[[#1]]][[2]], Sub2List[coordinates[[#2]]][[2]], ene}\[Transpose], 
                                                       PlotStyle -> {White, PointSize[Large]}, 
                                                       AxesLabel -> {labels[[#1]], labels[[#2]],""}, 
                                                       PlotRange ->  All, 
                                                       ViewPoint -> {0, 0, 100}, ImageSize -> 300], 
                                            ParametricPlot3D[BSplineFunction[{Sub2List[coordinates[[#1]]][[2]],Sub2List[coordinates[[#2]]][[2]],ene}\[Transpose]][s], {s, 0, 1}, 
                                                       PlotStyle -> Directive[White, Thick], 
                                                       PlotRange -> All, 
                                                       ViewPoint -> {0, 0, 100}]], 
                                                       ImageSize->300] 
                   & @@@ Subsets[Table[If[Norm[Sub2List[coordinates[[j]]][[2]]] == 0, ## &[], j], {j, nvars*ncells+6}], {2}], UpTo[4]];
        plt = {epsplot, scfplot, roseplot};

        If[OptionValue["table"], Print[path]];
        If[OptionValue["plot"], Print[Grid[{{Grid@epsplot},{Grid@roseplot},{Grid@Partition[scfplot, 2]}}]]];
        AppendTo[plots, {ts, inv, plt}];
        {ts, Prepend[coordinates, ene], {lattices, sites, stress, forces}\[Transpose]}, {i, Length[tsdata]}];

  If[OptionValue["pdf"],
     Export[dir0<>"/seeds_info.pdf", CreateDocument[Grid[Join[{{Graphics@Text[Style["ts-"<>ToString[#[[1]]]<>": "<>#[[2]], Red, Bold, 12]]}},
                                                              {#[[3,2]]}\[Transpose],
                                                              {Flatten[#[[3,1]]]}\[Transpose],
                                                              {Flatten[#[[3,3]]]}\[Transpose]], Frame -> All, FrameStyle -> Red] & /@ plots,
                                          PageBreakBelow -> True, Visible -> False]]];

  Return[out]
]

SamplingByDMDPathOld[dir0_, dmdpath_, invariants_, vars_, svars_, Basis_, boundingfactor_, OptionsPattern[{"fast" -> True}]] := Module[{tsdir, bddir, dist, p, icount, icount0, i, j, numdmd, idmd, jdmd, id, s, im, m, models, character, initcharacter, coordinates, trajectory, shiftphase, shiftmode, shifteij, strainlatt, pos, pos0, pos1, fullvars, NumBasis, invdict, seeds0, disp, info, phase0, phase1, phase, latt, sites, ij, numsamples, slavemode}, 
  NumBasis = Length[Flatten[vars]];
  fullvars = Flatten[{vars, svars}];
  pos0 = ImportPOSCAR[dir0 <> "/POSCAR0"];
  invdict = GroupBy[Flatten[invariants], StrainInvQ[#, svars[[1, 1]]] &];
  models = Flatten[Table[SortInvariants[t, fullvars], {t, {invdict[False], invdict[True]}}], 1];

  If[! DirectoryQ[dir0 <> "/trainingset"], CreateDirectory[dir0 <> "/trainingset"]; Run["chmod o-w " <> dir0 <> "/trainingset"]];
  If[! DirectoryQ[dir0 <> "/boundingset"], CreateDirectory[dir0 <> "/boundingset"]; Run["chmod o-w " <> dir0 <> "/boundingset"]];
  Print["Number of invariants :" <> ToString[Length[Flatten@models]]];
  Print["Groups of models :" <> ToString[Length[models]]];
  
  info = {};
  numsamples = 0;
  Do[m = models[[im]];
     tsdir = dir0 <> "/trainingset/ts-" <> ToString[im];
     bddir = dir0 <> "/boundingset/ts-" <> ToString[im];
     icount = 0;
     icount0 = 0;
     If[! DirectoryQ[tsdir], CreateDirectory[tsdir]; Run["chmod o-w " <> tsdir]];
     If[! DirectoryQ[bddir], CreateDirectory[bddir]; Run["chmod o-w " <> bddir]];
     character = InvariantCharacter[m[[1]], fullvars];
   
     dist = Table[phase = First[First[p]\[Transpose]][[2;;]];
                  initcharacter = RoundPhase[phase]; 
                  Abs[Normalize[v].Normalize[initcharacter]], {v, character}, {p, dmdpath}];
     Do[{i, j} = ij;
        {coordinates, trajectory} = dmdpath[[j]];
        numdmd = Length[trajectory];
 
        phase0 = First[coordinates\[Transpose]][[2;;]];
        phase1 = Last[coordinates\[Transpose]][[2;;]];

        If[RoundPhase[phase0[[1;;NumBasis]]]!=RoundPhase[phase1[[1;;NumBasis]]],
           Do[numsamples = numsamples + 1;
              icount = icount + 1;
              icount0 = icount0 + 1;
              phase = coordinates\[Transpose][[idmd]][[1 ;; -2]];
              shiftphase = (character[[i]] - 1)*phase;
              shiftmode = Table[{id, "modes", 1, shiftphase[[id]]}, {id, NumBasis}];
              strainlatt = First[trajectory[[idmd]]];
              pos = {strainlatt, ImposeMode[pos0[[1]], {trajectory[[idmd]][[2]], pos0[[2]]\[Transpose][[2]]}\[Transpose], Basis, shiftmode, 1.0]};
              ExportPOSCAR[tsdir, "sample." <> ToString[icount] <> ".vasp", pos], {idmd, numdmd}]
          ];
  
        phase1 = Join[phase1[[1;;NumBasis]], Chop[phase1[[NumBasis+1;;]], 10^-3]];
        seeds0 = If[OptionValue["fast"], 
                    #*Sign[phase1] &/@ Range[boundingfactor], 
                    DeleteCases[FetchSeeds[First[m][[i]], fullvars, boundingfactor, "zero" -> True], a_ /; Norm[a] == 0]];
        disp = Table[If[phase1[[id]] == 0., 
                        0, 
                        If[Abs[StandardDeviation[coordinates[[id]]]/phase1[[id]]] < 2, 
                           Abs[phase1[[id]]], 
                           Max@Abs[coordinates[[id]]]]], {id, NumBasis + 6}];
        shiftphase = (character[[i]] - 1)*phase1;
        shiftmode = Table[{id, "modes", 1, shiftphase[[id]]}, {id, NumBasis}];
        latt = First[Last[trajectory]];
        sites = ImposeMode[pos0[[1]], {Last[trajectory][[2]], pos0[[2]]\[Transpose][[2]]}\[Transpose], Basis, shiftmode, 1.0];
        Do[If[AllTrue[s*phase1, TrueQ[# >= 0.] &],
           numsamples = numsamples + 1;
           icount = icount + 1;
           shiftphase = disp s;
           shiftmode = Table[{id, "modes", 1, shiftphase[[id]]}, {id, NumBasis}];
           shifteij = Chop@eta2eij[shiftphase[[NumBasis + 1 ;;]]];
           strainlatt = (shifteij + IdentityMatrix[3]).latt;
           pos = {strainlatt, ImposeMode[pos0[[1]], sites, Basis, shiftmode, 1.0]};
           ExportPOSCAR[bddir, "sample." <> ToString[icount] <> ".vasp", pos]], {s, seeds0}], {ij, Position[dist, Max[dist], {2}]}];

     AppendTo[info, "ts-" <> ToString[im] <> " (" <> ToString[icount0] <> "+" <> ToString[icount - icount0] <> "): " <> ToString[m[[1]], StandardForm]], {im, Length[models]}];

  Print[Grid[Partition[info, UpTo[2]], ItemSize -> Full, Alignment -> ":"]];
  Print["Total number of DMD: " <> ToString[Length[Flatten[#\[Transpose] & /@ (dmdpath\[Transpose][[1]]), 1]]]];
  Print["Total number of samples: " <> ToString[numsamples]];
]

SamplingByDMDPath[dir0_, DMDPath_, vars_, svars_, uvars_, opotim_, spotim_, BasisMatrix_, FullBasis_, seeds_, OptionsPattern[{"SlaveFreeQ"->False, "SwitchOffHiddenQ"->False, "FrozenHiddenQ"->True, "spline" -> 0, "mddensity" -> 1.0, "sigma" -> 1.1, "extension"->0, "table"->False, "fontsize"->12, "write"->True, "plot"->True, "imagesize"->200, "interpolate"->3}]] := Module[{ref, mim0, character, scharacter, supervars, fullvars, excharacter, miminfo, data, i, j, k, samples, DMDSamples, InterpolatedSamples, nvars, dft, sigma, osparse, ssparse, plt, plots, pos0, pos, FullBasisMatrix, ComplementBasisMatrix, tcount, scount, dir, shiftmode, shifteij, latt, sites, PickQ, imd, dmd, inv, exlabels, labels, info, numsamples, data2d, enedata, exi, exene, exvars, exsvars, exovars, extension, exseeds, s, s1, bsplinepath, nsamples, splinedata, tab, color, slavemode, primarymode, hiddenmode, mask, primarymask, slavemask, hiddenmask, tsdata, ts, t, pickQ, ncells, grid, ExInfo, ind},
  If[! DirectoryQ[dir0 <> "/trainingset"], CreateDirectory[dir0 <> "/trainingset"]; Run["chmod o-w " <> dir0 <> "/trainingset"]];
  pos0 = ImportPOSCAR[dir0 <> "/POSCAR0"];
  FullBasisMatrix = Flatten[FullBasis, 1]\[Transpose][[4]]\[Transpose];
  {ComplementBasisMatrix, exlabels} = BasisComplement[BasisMatrix, FullBasis];
  FullBasisMatrix = Join[BasisMatrix\[Transpose], ComplementBasisMatrix\[Transpose]]\[Transpose];

  ind = Sort[First@First[Position[Flatten[FullBasis, 1]\[Transpose][[4]], #]] & /@ (ComplementBasisMatrix\[Transpose])];
  ExInfo = Tally[First@First@Position[FullBasis, #] & /@ ind];

  nvars = Length[Flatten[Append[vars, uvars]]];

  sigma = OptionValue["sigma"];
  osparse = Flatten[ConstantArray[#2, Length[#1]] & @@@ ({Append[vars, uvars], Append[opotim, 1]}\[Transpose])];
  ssparse = Join[ConstantArray[spotim[[1]], 3], ConstantArray[spotim[[2]], 3]];

  tsdata = {};
  plots = {};
  info = {};
  numsamples = 0;
  Do[scount = 0;
     tcount = DMDPath[[imd]][[1]];
     inv = If[tcount<=Length[seeds], ToString[seeds[[tcount,2]], StandardForm], "unknown"];
     dir = dir0 <> "/trainingset/ts-" <> ToString[tcount];
     If[! DirectoryQ[dir], CreateDirectory[dir]; Run["chmod o-w " <> dir]];
     supervars = First@Sub2List[First[DMDPath[[imd]][[2]]\[Transpose]][[2;;]]];
     grid = VarSubGrid[First[DMDPath[[imd]][[2]]\[Transpose]][[2;;]]];
     ncells = Times@@grid;
     fullvars = Append[Flatten[ConstantArray[Append[vars, uvars], ncells], 1], svars];
     labels = MapIndexed["("<>ToString[First@#2]<>") "<>ToString[VarBare@#1,StandardForm]&, supervars];
     dft = Join[{#[[1]]}, Sub2List[#[[2;;]]][[2]]] &/@ (DMDPath[[imd]][[2]]\[Transpose]);
     (*
     character = Sign[Abs@dft[[1]][[2 ;; nvars + 1]]];
     scharacter = Sign[Abs@dft[[1]][[nvars + 2 ;; nvars + 7]]];
     excharacter = Sign[Mean[Abs@dft][[nvars + 8 ;;]]];
     dft = Join[{1}, character, scharacter, excharacter] # & /@ (DMDPath[[imd]][[2]]\[Transpose]);
     *)
     mim0 = Mean[MinimalBy[dft, First]];
     mim0 = Last@dft;

     mask = Sign@Mean@Abs[dft[[;;,2;;]]];
     hiddenmask = Join[ConstantArray[0, ncells*nvars+6], ConstantArray[1, Length[supervars]-ncells*nvars-6]];
     primarymask = Sign@Join[Chop[Abs[First[dft][[2;;1+ncells*nvars]]], 0.01], 
                             Chop[Abs[First[dft][[2+ncells*nvars;;6+ncells*nvars+1]]], 0.001],
                             Chop[Abs[First[dft][[6+ncells*nvars+2;;]]]]]*(1-hiddenmask);

     slavemask = (mask - primarymask)*(1-hiddenmask);
     hiddenmode = hiddenmask*mim0[[2;;]];

     miminfo = Table[{i - 1, mim0[[i]], Sqrt@Mean[(# - mim0[[i]])^2 & /@ (dft\[Transpose][[i]])]}, {i, ncells*nvars + 7}]; 
     samples = Table[PickQ = Fold[And, Table[Abs[dft[[i]][[j + 1]] - miminfo[[j + 1, 2]]] <= sigma*miminfo[[j + 1, 3]], {j, ncells*nvars + 6}]];
                     If[PickQ, 
                        primarymode = primarymask*dft[[i,2;;]];
                        slavemode = If[OptionValue["SlaveFreeQ"], 0, 1] slavemask*mim0[[2;;]];
                        (*slavemode = If[OptionValue["SlaveFreeQ"], 0, 1] 
                                      MapIndexed[If[Chop[#1]==0., 0, Extract[slavemask*dft[[i,2;;]], #2]]&, slavemask*mim0[[2;;]]];*)
                        shiftmode = primarymode + slavemode + If[OptionValue["FrozenHiddenQ"], If[OptionValue["SwitchOffHiddenQ"], 0*hiddenmode, hiddenmode], hiddenmask*dft[[i,2;;]]];
                        {i, 
                         dft[[i]][[1]], 
                         Table[N@Round[shiftmode[[j]], osparse[[Mod[j,nvars,1]]]], {j, ncells*nvars}], 
                         Table[N@Round[shiftmode[[ncells*nvars + j]], ssparse[[j]]], {j, 6}],
                         shiftmode[[ncells*nvars + 7 ;;]]}, 
                        ## &[]], {i, Length[dft]}];

     ts = {};
     ref = mim0;
     Do[i = Length@dft - j + 1;
        pickQ = ! And @@ Thread[Part[Abs[(dft[[i]] - ref)], Span[2, 1 + ncells*nvars + 6]] <= OptionValue["mddensity"]*Join[Flatten@ConstantArray[osparse, ncells], ssparse]];
        If[pickQ, AppendTo[ts, dft[[i]]]; ref = dft[[i]]], {j, Length@dft}];
     AppendTo[tsdata, Table[Flatten[{imd, i, Reverse[ts][[i]]}], {i, Length@ts}]];

     exene   = dft[[1]][[1]] - dft[[2]][[1]];
     exvars  = (Sign[dft[[1]][[2;;ncells*nvars+7]]]*(miminfo\[Transpose][[3]][[2;;]]))[[1;;ncells*nvars]];
     exsvars = (Sign[dft[[1]][[2;;ncells*nvars+7]]]*(miminfo\[Transpose][[3]][[2;;]]))[[ncells*nvars+1;;]];
     exovars = Chop[dft[[-1]][[ncells*nvars + 8 ;;]], 0.01];

     extension = Table[{Length[dft]+s, dft[[1]][[1]] + exene*s, exvars*sigma*s, exsvars*sigma*s, exovars}, {s, 1, OptionValue["extension"]}];

     DMDSamples = Join[extension,
                       DeleteDuplicates[samples, Chop[Norm[Flatten[#1[[3 ;; 4]] - #2[[3 ;; 4]]]]] == 0 &],
                       {{Last[samples][[1]]+1, mim0[[1]], mim0[[2;;ncells*nvars+1]], mim0[[ncells*nvars+2;;ncells*nvars+7]], If[OptionValue["SwitchOffHiddenQ"], 0, 1]*mim0[[ncells*nvars+8;;]]}}];

     bsplinepath = BSplineFunction[Join[{#2}, #3, #4, #5] &@@@ DMDSamples];

     InterpolatedSamples = If[OptionValue["spline"] == 0, 
                              DMDSamples, 
                              nsamples = Ceiling[Length[DMDSamples]*OptionValue["spline"]];
                              Table[splinedata = bsplinepath[N[i/nsamples]];
                                    {i, splinedata[[1]], splinedata[[2;;ncells*nvars+1]], splinedata[[ncells*nvars+2;;ncells*nvars+7]], splinedata[[ncells*nvars+8;;]]}, {i, 0, nsamples}]
                              ];
     PrependTo[InterpolatedSamples, {0, First[dft][[1]], First[dft][[2;;ncells*nvars+1]], First[dft][[ncells*nvars+2;;ncells*nvars+7]], If[OptionValue["SwitchOffHiddenQ"], 0, 1]*mim0[[ncells*nvars + 8 ;;]]}];
     PrependTo[InterpolatedSamples, {0, First[dft][[1]], First[dft][[2;;ncells*nvars+1]], First[dft][[ncells*nvars+2;;ncells*nvars+7]], If[OptionValue["SwitchOffHiddenQ"], 0, 1]*mim0[[ncells*nvars + 8 ;;]]}];
     PrependTo[InterpolatedSamples, {0, First[dft][[1]], First[dft][[2;;ncells*nvars+1]], First[dft][[ncells*nvars+2;;ncells*nvars+7]], 0.0*First[dft][[ncells*nvars + 8 ;;]]}];

     tab = Grid[Prepend[Flatten[#] &/@ InterpolatedSamples, Join[{"#", "energy"}, labels]],
                Background -> {Join[{Gray,LightBlue}, Flatten[Table[color = If[OddQ[i], Yellow, Pink];
                               ConstantArray[color, Length[fullvars[[i]]]], {i, Length@fullvars}]],
                               Flatten[Table[ConstantArray[If[OddQ[i], Gray, LightGray], ExInfo[[Mod[i,Length[ExInfo],1],2]]], {i,       ncells*Length[ExInfo]}]]],
                               None},
                ItemStyle -> Directive[FontSize -> OptionValue["fontsize"]],
                ItemSize -> Full
               ];
 
     If[OptionValue["table"], Print[tab]];

     data2d = (Flatten[#] & /@ (InterpolatedSamples[[3;;]]\[Transpose][[3 ;; 4]]\[Transpose]))\[Transpose];
     enedata = InterpolatedSamples[[3;;]]\[Transpose][[2]];
     plt = Flatten[{ListPlot[InterpolatedSamples\[Transpose][[3]][[3;;]]\[Transpose], 
                             PlotRange -> All,
                             ImageSize -> OptionValue["imagesize"], Frame -> True, GridLines -> Automatic, 
                             PlotMarkers -> Automatic, PlotLabel -> "trajectory", 
                             FrameLabel -> {"Damped MD step (#)", "coordinates (\[Angstrom])"}, 
                             ImagePadding -> {{60, 10}, {40, 10}}],
                    ListPlot[InterpolatedSamples[[2;;]]\[Transpose][[4]][[3;;]]\[Transpose],
                             PlotRange -> All,
                             ImageSize -> OptionValue["imagesize"], Frame -> True, GridLines -> Automatic,
                             PlotMarkers -> Automatic, PlotLabel -> "trajectory",
                             FrameLabel -> {"Damped MD step (#)", "strains"},
                             ImagePadding -> {{60, 10}, {40, 10}}],
                    Table[If[Norm[InterpolatedSamples[[3;;]]\[Transpose][[3]]\[Transpose][[i]]] != 0, 
                             ListPlot[{InterpolatedSamples[[3;;]]\[Transpose][[3]]\[Transpose][[i]], enedata}\[Transpose], 
                                      PlotRange -> All,
                                      ImageSize -> OptionValue["imagesize"], Frame -> True, GridLines -> Automatic, 
                                      PlotMarkers -> Automatic, 
                                      PlotLabel -> "(" <> ToString[i] <> ") " <> ToString[Flatten[fullvars][[i]], StandardForm], 
                                      FrameLabel -> {"coordinates (\[Angstrom])", "energy (eV)"}, 
                                      ImagePadding -> {{60, 10}, {40, 10}}], ## &[]], {i, nvars}],
                    Table[If[Norm[InterpolatedSamples[[3;;]]\[Transpose][[4]]\[Transpose][[i]]] != 0,
                              ListPlot[{InterpolatedSamples[[3;;]]\[Transpose][[4]]\[Transpose][[i]], enedata}\[Transpose],
                                       PlotRange -> All,
                                       ImageSize -> OptionValue["imagesize"], Frame -> True, GridLines -> Automatic,
                                       PlotMarkers -> Automatic,
                                       PlotLabel -> "(" <> ToString[i] <> ") " <> ToString[svars[[i]], StandardForm],
                                       FrameLabel -> {"coordinates (\[Angstrom])", "energy (eV)"},
                                       ImagePadding -> {{60, 10}, {40, 10}}], ## &[]], {i, 6}],
                    Rasterize[Show[ListPlot3D[{data2d[[#1]], data2d[[#2]], enedata}\[Transpose],
                                         AxesLabel -> {ToString[Flatten[fullvars][[#1]],StandardForm], ToString[Flatten[fullvars][[#2]],StandardForm], ""},
                                         PlotRange -> All, ColorFunction -> "Rainbow",
                                         InterpolationOrder -> OptionValue["interpolate"],
                                         MeshFunctions -> {#3 &}, Mesh -> Length[enedata],
                                         ViewPoint -> {0, 0, 100}, ImageSize -> OptionValue["imagesize"]],
                                   ListPointPlot3D[{data2d[[#1]], data2d[[#2]], enedata}\[Transpose],
                                          PlotStyle -> {White, PointSize[Large]},
                                          AxesLabel -> {ToString[Flatten[fullvars][[#1]],StandardForm], ToString[Flatten[fullvars][[#2]],StandardForm], ""},
                                          PlotRange -> All,
                                          ViewPoint -> {0, 0, 100}, ImageSize -> OptionValue["imagesize"]],
                                   ParametricPlot3D[BSplineFunction[{data2d[[#1]],data2d[[#2]],enedata}\[Transpose]][s], {s, 0, 1},
                                          PlotStyle -> Directive[White, Thick],
                                          PlotRange -> All,
                                          ViewPoint -> {0, 0, 100}, ImageSize -> OptionValue["imagesize"]]],
                                    ImageSize -> OptionValue["imagesize"]]
                                         & @@@ Subsets[Table[If[Norm[data2d[[i]]] == 0, ## &[], i], {i, Length[Flatten@fullvars]}], {2}]}];

     If[OptionValue["plot"], Print[Grid[{{Text[Style["ts-"<>ToString[tcount]<>": "<>inv, Red, Bold, 12]]}, {plt[[1]], plt[[2]]}, plt[[3;;]]}]]];
     AppendTo[plots, {tcount, inv, plt}];

     Do[scount = scount + 1;
        numsamples = numsamples + 1;
        shiftmode = Table[i = Mod[j, nvars, 1];
                          If[i <= nvars - 3, {i, supervars[[j]], 1, s[[3, j]]}, ##&[]], {j, ncells*nvars}];
        shiftmode = Join[shiftmode, 
                         Table[i = Mod[j, Length[FullBasisMatrix] - nvars + 3, 1];
                               {i+nvars-3, supervars[[ncells*nvars+6+j]], 1, s[[5, j]]}, {j, Length[supervars]-ncells*nvars-6}]];
        shiftmode = GroupBy[shiftmode, Mod[VarGrid[#[[2]]] + 1, grid, 1]&];
        shifteij = Chop@eta2eij[s[[4]]];
        pos = ImposeDomains[pos0, FullBasisMatrix, shifteij, shiftmode, grid];
        If[OptionValue["write"],  ExportPOSCAR[dir, "sample." <> ToString[scount] <> ".vasp", pos]], {s, InterpolatedSamples}];

     AppendTo[info, "ts-" <> ToString[tcount] <> " (" <> ToString[scount] <> "): " <> inv], {imd, Length@DMDPath}];

     tsdata = MapIndexed[{First@#2, #1[[1]], #1[[2]], #1[[3]], #1[[4;;3+ncells*nvars]], #1[[4+ncells*nvars;;9+ncells*nvars]], #1[[10+ncells*nvars;;]]}&, Flatten[tsdata, 1]];
     Export[dir0 <> "/trainingset/tsmd.dat", GatherBy[tsdata, #[[2]] &], "ExpressionML"];
  
  Print[Grid[Join[{{"Total number of samples: "<>ToString[numsamples]}}, Partition[info, UpTo[5]]], ItemSize -> Full, Alignment -> ":"]];
]

DMDPath2ts[dmd_, vars_, svars_] := Module[{npath, nsample, ncount, ts, tsdata, dft, pos, nvars, out},
  nvars = Length[Flatten[vars]];
  npath = Length[dmd];
  ncount = 0;
  out = Table[{ts, tsdata, pos} = dmd[[i]];
              nsample = Length[tsdata\[Transpose]];
              ncount = ncount + 1;
              Table[dft = tsdata\[Transpose][[j]]; 
                    {ncount, ts, j, dft[[1]], dft[[2 ;; nvars + 1]], dft[[nvars + 2 ;; nvars + 7]], dft[[nvars + 8 ;;]]}, {j, nsample}], {i, npath}];
  Return[out]
]

SamplePhaseSpace[dir0_, s0_, spgmat_, Basis_, vars_, opotim_, spotim_, OptionsPattern[{"seeds" -> {}}]] := Module[{BasisDim, NumBasis, modesets, seeds, modes, m, i, pos, pos0, eij, strainlatt, potim, seedname},
  seedname = StringRiffle[Join[opotim, spotim], "_"];
  {BasisDim, NumBasis} = Dimensions[Basis];
  potim = Flatten[ConstantArray[#2, Length[#1]] & @@@ ({vars, opotim}\[Transpose])];
  If[! DirectoryQ[dir0<>"/seeds"<>"-"<>seedname], 
     CreateDirectory[dir0<>"/seeds"<>"-"<>seedname]; 
     Run["chmod o-w " <> dir0<>"/seeds"<>"-"<>seedname]];
  pos0 = ImportPOSCAR[dir0 <> "/" <> s0];

  If[OptionValue["seeds"] == {}, 
     modesets = SortBy[Tuples[{-1, 0, 1}, {NumBasis}], Norm[N@#] &];
     seeds = Round[DeleteDuplicates[Table[Sort[# . modesets[[i]] & /@ spgmat], {i, Length[modesets]}]]\[Transpose][[1]]],
     seeds = OptionValue["seeds"]
  ];

  Print["Number of seeds: " <> ToString[Length[seeds]]<>" (potim="<>seedname<>")"];
  Do[eij=Chop@eta2eij[Flatten[ConstantArray[#, 3]&/@spotim]*seeds[[i]][[-6;;]]];
     strainlatt=(eij+IdentityMatrix[3]).Transpose[pos0[[1]]];
     modes = Table[{m, "modes", 1, potim[[m]] seeds[[i, m]]}, {m, NumBasis}]; 
     pos = {strainlatt, ImposeMode[pos0[[1]], pos0[[2]], Basis, modes, 1.0]}; 
     ExportPOSCAR[dir0<>"/seeds"<>"-"<>seedname, "seed." <> ToString[i] <> ".vasp", pos];
     ExportOPTCELL[dir0<>"/seeds"<>"-"<>seedname<>"/OPTCELL." <> ToString[i] <> ".inp", Abs[seeds[[i]][[-6;;]]], "constrain"->False], {i, Length@seeds}]
]

CollectSeeds[dir0_, vars_, svars_, BasisMatrix_, FullBasis_, seeds_, OpMat_, OptionsPattern[{"round"->10^-6, "file"->"POSCAR", "all"->False, "table" -> True, "fontsize"->12}]] := Module[{pos, mode, eta, pos0, NumDim, NumBasis, Basislabel, i, file, phases, m, ind, OtherBasisMatrix, ExMode, ExBasislabel, out, FullBasisMatrix, ExInfo, t, ene, ene0, isg0, isg, sg0, sg, sginfo, data, fullvars},
  fullvars = Flatten[{vars, svars}];
  FullBasisMatrix = Flatten[FullBasis, 1]\[Transpose][[4]]\[Transpose];
  OtherBasisMatrix = Complement[FullBasisMatrix\[Transpose], BasisMatrix\[Transpose], SameTest -> (Chop[Norm[#1 - #2]] == 0 &)];
  ind = Sort[First@First[Position[FullBasisMatrix\[Transpose], #]] & /@ OtherBasisMatrix];
  ExInfo = Tally[First@First@Position[FullBasis, #] & /@ ind];
  OtherBasisMatrix = (FullBasisMatrix\[Transpose][[#]] & /@ ind)\[Transpose];

  data = Import[dir0 <> "/seeds/" <> "ene.dat", "Table"];
  ene0 = data[[1,2]];

  {NumDim, NumBasis} = Dimensions[BasisMatrix];
  pos0 = ImportPOSCAR[dir0 <> "/POSCAR0"];
  Basislabel = Table["x" <> ToString[i], {i, NumBasis}];
  ExBasislabel = "("<>ToString[#]<>") "<>(Flatten[FullBasis, 1]\[Transpose][[2]][[#]]) &/@ ind;
  phases = Table[t = data[[i, 1]];
                 ene = data[[i, 2]] - ene0;
                 isg0 = data[[i, 3]];
                 sg0 = data[[i, 4]];
                 isg = data[[i, 5]];
                 sg = data[[i, 6]];
                 sginfo = "#"<>ToString[isg0]<>" "<>sg0<>"\[RightArrow]"<>"#"<>ToString[isg]<>" "<>sg;
                 file = StringReplace[OptionValue["file"], "xxx"->ToString[i]];
                 pos = ImportPOSCAR[dir0 <> "/seeds/seed" <> ToString[t] <> "/" <> file];
                 mode = Chop[ISODISTORT[pos0[[1]], pos0[[2]], pos[[2]], BasisMatrix, Basislabel, "round" -> OptionValue["round"]]\[Transpose][[4]], OptionValue["round"]];
                 ExMode = Chop[ISODISTORT[pos0[[1]], pos0[[2]], pos[[2]], OtherBasisMatrix, ExBasislabel, "round" -> OptionValue["round"]]\[Transpose][[4]], OptionValue["round"]];
                 eta = N@Round[eij2eta[GetStrainTensor[pos0[[1]], pos[[1]], "iso" -> False]], 10^-3];
                 {t, sginfo, ene, mode, eta, ExMode, pos, Table[Sign[#]&/@Chop[m.mode,OptionValue["round"]], {m, OpMat}]}, {i, Length[data]}];

  out = If[OptionValue["all"], phases\[Transpose][[1 ;; 7]]\[Transpose],
                                 DeleteDuplicates[phases, Complement[#1[[8]],#2[[8]]] == {} &]\[Transpose][[1 ;; 7]]\[Transpose]];
  If[OptionValue["table"], Print@ShowSeedPhases[out, seeds, Append[vars, svars], ExBasislabel, ExInfo, "fontsize" -> OptionValue["fontsize"]]];

  Return[out] 
]

PreparePhases[dir0_, pname_, num_] := Module[{pos0, sites, n, pos},
  If[! DirectoryQ[dir0 <> "/phases"], 
     CreateDirectory[dir0 <> "/phases"];
     Run["chmod o-w " <> dir0 <> "/phases"]];
  pos0 = ImportPOSCAR[dir0 <> "/POSCAR0"];
  Do[pos = ImportPOSCAR[dir0 <> "/" <> pname <> ToString[n] <> ".vasp"];
     sites = SortByPOSCAR[pos0, pos[[2]]\[Transpose][[1]]];
     ExportPOSCAR[dir0 <> "/phases", "phase." <> ToString[n] <> ".vasp", {pos0[[1]], sites}], {n, num}];
]

CollectPhases[dir0_, vars_, svars_, BasisMatrix_, FullBasis_, OptionsPattern[{"AUTable"->True, "round" -> 10^-6, "fontsize" -> 12, "table" -> True, "file"->"POSCAR"}]] := Module[{NumDim, NumBasis, pos0, Basislabel, i, mpdata, pos, sites, eta, ExPhase, phase, simplifytab, fullvars, tab, ind, OtherBasisMatrix, ExBasislabel, color, fname, FullBasisMatrix, ExInfo, data, t, ev2har, ang2latt, ene, ene0, isg0, sg0, isg, sg, sginfo},

  FullBasisMatrix = Flatten[FullBasis, 1]\[Transpose][[4]]\[Transpose];
  OtherBasisMatrix = Complement[FullBasisMatrix\[Transpose], BasisMatrix\[Transpose], SameTest -> (Chop[Norm[#1 - #2]] == 0 &)];
  ind = Sort[First@First[Position[FullBasisMatrix\[Transpose], #]] & /@ OtherBasisMatrix];
  ExInfo = Tally[First@First@Position[FullBasis, #] & /@ ind];
  OtherBasisMatrix = (FullBasisMatrix\[Transpose][[#]] & /@ ind)\[Transpose];

  data = Import[dir0 <> "/phases/" <> "ene.dat", "Table"];
  ene0 = data[[1,2]];

  fullvars = Append[vars, svars];
  {NumDim, NumBasis} = Dimensions[BasisMatrix];
  pos0 = ImportPOSCAR[dir0 <> "/POSCAR0"];
  Basislabel = Table["x" <> ToString[i], {i, NumBasis}];
  ExBasislabel = "("<>ToString[#]<>") "<>(Flatten[FullBasis, 1]\[Transpose][[2]][[#]]) &/@ ind;
  simplifytab = {-1 -> Style["-1", Bold, Blue], 1 -> Style["1", Bold, Red], 0 -> "-"};
  If[OptionValue["AUTable"], ev2har=1/27.2114;ang2latt=1/Norm[pos0[[1]]], ev2har=1;ang2latt=1];

  mpdata = Table[t = data[[i, 1]];
                 ene = data[[i, 2]] - ene0;
                 isg0 = data[[i, 3]];
                 sg0 = data[[i, 4]];
                 isg = data[[i, 5]];
                 sg = data[[i, 6]];
                 sginfo = "#"<>ToString[isg0]<>" "<>sg0<>"\[RightArrow]"<>"#"<>ToString[isg]<>" "<>sg;
                 fname = If[MemberQ[{"POSCAR", "CONTCAR"}, OptionValue["file"]], OptionValue["file"], OptionValue["file"]<>"."<>ToString[i]<>".vasp"];
                 pos = ImportPOSCAR[dir0 <> "/phases/phase" <> ToString[t] <> "/" <> fname];
                 sites = {PosMatchTo[pos0[[1]], pos0[[2]]\[Transpose][[1]], pos[[2]]\[Transpose][[1]]][[2]], pos0[[2]]\[Transpose][[2]]}\[Transpose];
                 eta = N@Round[eij2eta[GetStrainTensor[pos0[[1]], pos[[1]], "iso" -> False]], 10^-4];
                 phase = ISODISTORT[pos0[[1]], pos0[[2]], sites, BasisMatrix, Basislabel, "round" -> OptionValue["round"], "match"->False]\[Transpose][[4]];
                 ExPhase = ISODISTORT[pos0[[1]], pos0[[2]], sites, OtherBasisMatrix, ExBasislabel, "round" -> OptionValue["round"], "match"->False]\[Transpose][[4]];
                 {i, sginfo, ene, phase, eta, ExPhase}, {i, 1, Length[data]}];
  tab = Grid[Prepend[Flatten[{Join[{#1, StringSplit[#2, "\[RightArrow]"][[2]], #3*ev2har}, #4*ang2latt, #5, #6*ang2latt], Join[{" ", " ", " "}, Sign[#4], Sign[#5], Sign[#6]] /. simplifytab} & @@@ mpdata, 1], Join[{"#", "sg", "ene"}, Flatten@fullvars, ExBasislabel]],
             Background -> {Flatten[{{LightGray, LightBlue, Green}, Table[color = If[OddQ[i], Yellow, Pink]; ConstantArray[color, Length[fullvars[[i]]]], {i, Length@fullvars}], Table[ConstantArray[If[OddQ[i], Gray, LightGray], ExInfo[[i,2]]], {i,Length[ExInfo]}]}], Prepend[Flatten[Table[If[EvenQ[i], White, None], {i, 2 Length[data]}]], White]}, 
             Dividers -> {False, Table[If[EvenQ[i], i -> Black, ## &[]], {i, 2 Length[data]}]}, 
             ItemStyle -> Directive[FontSize -> OptionValue["fontsize"]], 
             ItemSize -> Full];
  If[OptionValue["table"], Print[tab]];
  Return[mpdata]
]

ShowSeedPhases[phases_, seeds_, vars_, exlabels_, exinfo_, OptionsPattern[{"fontsize" -> 12}]] := Module[{i, color,labels, fullvars, simplifytab},
  fullvars = Flatten[vars];
  labels = MapIndexed["("<>ToString[First@#2]<>") "<>ToString[#1,StandardForm]&, fullvars];
  simplifytab = {-1 -> Style["-1", Bold, Blue], 1 -> Style["1", Bold, Red], 0 -> "-"};
  Grid[Prepend[Flatten[{Join[{ToString[Times@@DeleteCases[seeds[[#[[1]],3]],0] seeds[[#[[1]],2]], StandardForm]}, {StringSplit[#[[2]], "\[RightArrow]"][[1]], " "}, seeds[[#[[1]],3]] /. simplifytab, ConstantArray["-",Length[exlabels]]], 
                             Flatten[{Style[#[[1]],Red,Bold], Style[StringSplit[#[[2]], "\[RightArrow]"][[2]], Bold], #[[3;;6]]}]} 
                        & /@ (phases\[Transpose][[1 ;; 6]]\[Transpose]), 1], Flatten[Join[{"#", "sg", "ene"}, labels, exlabels]]], 
       Background -> {Flatten[{{LightGray, LightBlue, Green}, Table[color = If[OddQ[i], Yellow, Pink]; ConstantArray[color, Length[vars[[i]]]], {i, Length@vars}], Table[ConstantArray[If[OddQ[i], Gray, LightGray], exinfo[[i,2]]], {i,Length[exinfo]}]}],
                      Prepend[Flatten[Table[If[EvenQ[i], White, None], {i, 2 Length[phases]}]], White]},
       Dividers -> {False, Table[If[EvenQ[i], i -> Black, ## &[]], {i, 2 Length[phases]}]},
       ItemStyle -> Directive[FontSize -> OptionValue["fontsize"]],
       ItemSize->Full
      ]
]

SamplePolynomialSpace[polyphase_, dir_, spgmat0_, Basis_, npt_, smearing_, boundingfactor_, mim0_, OptionsPattern[{"round" -> 10^-8}]] := Module[{BasisDim, NumBasis, smodesets, vmodesets, modes, pos, pos1, phase, m, n, i, j, spgmat, character, sdisp, vdisp, c0, SlaveQ, EtaQ, vcharacter, scharacter,sv, sig0, signpt, eij, strainlatt, NumSamples, svartuple, vartuple, samplecount, varseeds},
  If[! DirectoryQ[dir], CreateDirectory[dir]; Run["chmod o-w " <> dir]];
  {BasisDim, NumBasis} = Dimensions[Basis];
  SlaveQ = polyphase[[1]];
  EtaQ = polyphase[[2]];
  character = polyphase[[3]];
  phase = polyphase[[4]];
  pos = polyphase[[5]];
  vcharacter = character[[1;;NumBasis]];
  scharacter = character[[NumBasis+1;;]];
  sig0 = If[SlaveQ, 4, 1.0];
  signpt = If[EtaQ, 2, 1];
  If[EtaQ, ExportOPTCELL[dir<>"/OPTCELL", scharacter]];

  spgmat = Table[If[Chop[Norm[m.Sign[phase] - Sign[phase]], OptionValue["round"]] == 0., m, ## &[]], {m, spgmat0}];

  varseeds = If[EtaQ, {0, -1, 1}, Flatten[{0, Table[{-i, i}, {i, Join[Range[npt], 2^Range[boundingfactor]*npt]}]}]];

  vartuple=SortBy[Select[Tuples[varseeds, {Length[vcharacter]}], Total[Abs[#]^1.0]^(1/1.0) <= 2^boundingfactor*npt &], Norm[N@#] &];

  svartuple=Flatten[Table[SortBy[Tuples[{0, i, -i}, {Length[scharacter]}], Norm[N@#] &], {i, Join[Range[signpt*npt], 2^Range[boundingfactor*signpt]*(signpt*npt)]}], 1];

  smodesets=DeleteDuplicates[If[Sign[Abs[#]] == scharacter, #, ## &[]] &/@ svartuple];
  vmodesets=DeleteDuplicates[If[(Norm[# (1 - vcharacter)] == 0) && (Sign[Abs[#] + Abs[phase]] == vcharacter), #, ## &[]] &/@ vartuple];
  PrependTo[vmodesets, ConstantArray[0, Length[vcharacter]]];
  vmodesets=Round[DeleteDuplicates[Table[{i, #.vmodesets[[i]] & /@ spgmat}, {i, Length[vmodesets]}], Complement[#1[[2]], #2[[2]]] == {} &]\[Transpose][[2]]\[Transpose][[1]]];

  NumSamples=Length[vmodesets]*Length[smodesets];
  samplecount=0;
  Do[modes = Table[c0 = sig0*If[phase[[m]] == 0., 0.5*mim0[[m]], phase[[m]]];
                   vdisp = If[Chop[c0]==0., Abs[0.1 (Sqrt[1. + Sign[0.1 vmodesets[[j, m]]] Sqrt[smearing]] - 1)]/(signpt*npt),
                                            If[Abs[vmodesets[[j, m]]]>(signpt*npt)&&Sign[c0*vmodesets[[j, m]]]<0,-1.0,1.0]*Abs[c0 (Sqrt[1. + Sign[c0 vmodesets[[j, m]]] Sqrt[smearing]] - 1)]/(signpt*npt)];
                   {m, "modes", 1, vdisp*vmodesets[[j, m]]}, {m, NumBasis}];
     sdisp={0.0005,0.0005,0.0005,0.0001,0.0001,0.0001};
     eij=Chop@eta2eij[sdisp*smodesets[[i]]];
     strainlatt=(eij+IdentityMatrix[3]).(pos[[1]]);
     pos1 = {strainlatt, ImposeMode[pos[[1]], pos[[2]], Basis, modes, 1.0]};
     samplecount=samplecount+1;
     ExportPOSCAR[dir, "sample." <> ToString[samplecount] <> ".vasp", pos1], {j,Length@vmodesets}, {i,Length@smodesets}];
  Return[NumSamples]
]

PlanSampling[dir0_, phases_, invariants_, vars_, svars_, Basis_, spgmat_, npt_, smearing_, boundingfactor_, OptionsPattern["strain"->"\[Epsilon]"]] := Module[{dist, SlaveQ, EtaQ, polyphase, poly0phase, cubicphase, p, v, m, t, i, j, id, models, modeldata, character, modes, plan, minimum, minimum0, shiftphase, shiftmode, pos0, ts1, ts0, fullvars, NumBasis, scharacter, vcharacter, invdict, NumSamples, info}, 
  NumBasis = Length[Flatten[vars]];
  fullvars = Flatten[{vars, svars}];
  pos0 = ImportPOSCAR[dir0<>"/POSCAR0"];
  minimum0 = Mean[(Table[Norm[vars[[i]]], {i, Length@vars}] /. (Thread[Flatten[vars] -> #]) & /@ (phases\[Transpose][[2]]))];
  minimum = Flatten[Table[Table[minimum0[[i]], {Length[vars[[i]]]}], {i, Length[vars]}]];

  invdict = GroupBy[Flatten[invariants], StrainInvQ[#, svars[[1, 1]]] &];
  models = Flatten[Table[SortInvariants[t, fullvars], {t, {invdict[False], invdict[True]}}], 1];
  Print["Number of invariants :" <> ToString[Length[Flatten@models]]];
  Print["Groups of models :" <> ToString[Length[models]]];

  info = {};
  ts0 = 1;
  ts1 = 1;
  Do[m=models[[id]];
     character = InvariantCharacter[m[[1]], fullvars];
     vcharacter = #[[1;;NumBasis]] &/@ character;
     scharacter = #[[NumBasis+1;;]] &/@ character;
     SlaveQ = AllTrue[ToString[Level[#, 1][[1]]] & /@ Variables[m[[1]]], LowerCaseQ];
     EtaQ = StrainInvQ[m[[1]], svars[[1, 1]]];
     dist = Table[Abs[Normalize[v].Normalize[RoundPhase[p[[2]]]]], {v, vcharacter}, {p, phases}];
     {i, j} = First@Position[dist, Max[dist]];
     shiftphase = (vcharacter[[i]] - 1)*phases[[j, 2]];
     shiftmode = Table[{id, "modes", 1, shiftphase[[id]]}, {id, Length[shiftphase]}];
     poly0phase = {SlaveQ, EtaQ, character[[i]], phases[[j, 2]]+shiftphase, {pos0[[1]], ImposeMode[pos0[[1]],phases[[j, 5]][[2]],Basis,shiftmode,1.0]}};
     polyphase = {SlaveQ, EtaQ, character[[i]], phases[[j, 2]]+shiftphase, {phases[[j, 5]][[1]], ImposeMode[pos0[[1]],phases[[j, 5]][[2]],Basis,shiftmode,1.0]}};
     cubicphase = {SlaveQ, EtaQ, character[[i]], phases[[1, 2]], phases[[1, 5]]};
     If[EtaQ,
       NumSamples=SamplePolynomialSpace[poly0phase, dir0 <> "/ts1-" <> ToString[ts1], spgmat, Basis, npt, smearing, boundingfactor, minimum];
       AppendTo[info, "ts1-"<>ToString[ts1]<>" ("<>ToString[NumSamples]<>"): "<>ToString[m[[1]], StandardForm]];
       ts1 = ts1 + 1,
       If[NumPolynomialVar[m[[1]]]==1,
          NumSamples=SamplePolynomialSpace[cubicphase, dir0 <> "/ts0-" <> ToString[ts0], spgmat, Basis, npt, smearing, boundingfactor, minimum];
          AppendTo[info, "ts0-"<>ToString[ts0]<>" ("<>ToString[NumSamples]<>"): "<>ToString[m[[1]], StandardForm]];
          ts0=ts0+1];
       NumSamples=SamplePolynomialSpace[polyphase, dir0 <> "/ts0-" <> ToString[ts0], spgmat, Basis, npt, smearing, boundingfactor, minimum];
       AppendTo[info, "ts0-"<>ToString[ts0]<>" ("<>ToString[NumSamples]<>"): "<>ToString[m[[1]], StandardForm]];
       ts0 = ts0 + 1;], {id, Length[models]}];
  Print[Grid[Partition[info, UpTo[3]], ItemSize -> Full, Alignment -> ":"]];
]

LoadTrainingset[dir0_, BasisMatrix_, FullBasis_, vars_, svars_, uvars_, OptionsPattern[{"round" -> 10^-6, "plot"->False, "pdf"->True, "imagesize"->200, "interpolate"->3}]] := Module[{i, t, s, BasisDim, NumBasis, ene, pos, mode0, vasprun, mode, exmode, data, eta, out, ws, seeds, plan, p, E0, pos0, FullBasisMatrix, ComplementBasisMatrix, exlabels, plt, pdf, fullvars, enedata, data2d},
  {BasisDim, NumBasis} = Dimensions[BasisMatrix];
  pos0 = ImportPOSCAR[dir0 <> "/POSCAR0"];
  FullBasisMatrix = Flatten[FullBasis, 1]\[Transpose][[4]]\[Transpose];
  {ComplementBasisMatrix, exlabels} = BasisComplement[BasisMatrix, FullBasis]; 
  FullBasisMatrix = Join[BasisMatrix\[Transpose], ComplementBasisMatrix\[Transpose]]\[Transpose];
  data = Import[dir0 <> "/trainingset/" <> "ene.dat", "Table"];
  E0 = data[[1, 3]];
  out = Table[t = data[[i, 1]]; 
              s = data[[i, 2]]; 
            ene = data[[i, 3]];
            pos = ImportPOSCAR[dir0 <> "/trainingset/ts" <> "-" <> ToString[t] <> "/sample" <> ToString[s] <> "/POSCAR"];
    {mode, eta} = ModeDecomposition[pos0, pos, BasisMatrix, vars, svars, uvars, 
                                            "lattnorm" -> False, "table" -> False, "round" -> OptionValue["round"]];
         exmode = ModeDecomposition[pos0, pos, ComplementBasisMatrix, Flatten@exlabels, svars, {}, 
                                       "lattnorm" -> False, "table" -> False, "round" -> OptionValue["round"]][[1]];
          (*eta = Chop[eij2eta[GetStrainTensor[pos0[[1]], pos[[1]], "iso" -> False]]];
           mode = ISODISTORT[pos0[[1]], pos0[[2]], pos[[2]], BasisMatrix, Range[NumBasis], "round" -> OptionValue["round"]]\[Transpose][[4]];
         exmode = ISODISTORT[pos0[[1]], pos0[[2]], pos[[2]], ComplementBasisMatrix, Flatten@exlabels, "round" -> OptionValue["round"]]\[Transpose][[4]];*)
         {i, t, s, ene - Length[mode]*E0, Flatten[mode], eta, Flatten[exmode]}, {i, Length[data]}];
  out = SortBy[#, First] & /@ GatherBy[out, #[[2]] &];

  pdf = Table[fullvars = Join[T5Vars[First[t]], svars];
              data2d = Transpose[T56Data[#] &/@ (t[[3;;]])]; 
              enedata = t[[3;;]]\[Transpose][[4]];
              plt = Join[{ListLinePlot[enedata, PlotRange -> All, Frame -> True, PlotMarkers -> Automatic, GridLines -> Automatic,
                                            ImageSize -> OptionValue["imagesize"], ImagePadding -> {{60, 10}, {40, 10}}]},
                         Table[If[Norm[T5Data[#][[i]]&/@(t[[3;;]])] != 0, 
                                  ListPlot[{T5Data[#][[i]]&/@(t[[3;;]]), enedata}\[Transpose], 
                                           PlotRange -> All,
                                           ImageSize -> OptionValue["imagesize"], 
                                           Frame -> True, 
                                           GridLines -> Automatic, 
                                           PlotMarkers -> Automatic, 
                                           PlotLabel -> "(" <> ToString[i] <> ") " <> ToString[fullvars[[i]], StandardForm], 
                                           FrameLabel -> {"coordinates (\[Angstrom])", "energy (eV)"}, 
                                           ImagePadding -> {{60, 10}, {40, 10}}], ## &[]], {i, Length[fullvars]-6}], 
                         Table[If[Norm[T6Data[#][[i]]&/@(t[[3;;]])] != 0,
                                  ListPlot[{T6Data[#][[i]]&/@(t[[3;;]]), enedata}\[Transpose],
                                           PlotRange -> All,
                                           ImageSize -> OptionValue["imagesize"],
                                           Frame -> True,
                                           GridLines -> Automatic,
                                           PlotMarkers -> Automatic,
                                           PlotLabel -> "(" <> ToString[i] <> ") " <> ToString[svars[[i]], StandardForm],
                                           FrameLabel -> {"coordinates (\[Angstrom])", "energy (eV)"},
                                           ImagePadding -> {{60, 10}, {40, 10}}], ## &[]], {i, 6}],
                         Rasterize[Show[ListPlot3D[{data2d[[#1]], data2d[[#2]], enedata}\[Transpose], 
                                              AxesLabel -> {ToString[fullvars[[#1]],StandardForm], ToString[fullvars[[#2]],StandardForm], ""}, 
                                              PlotRange -> All, ColorFunction -> "Rainbow", 
                                              InterpolationOrder -> OptionValue["interpolate"], 
                                              MeshFunctions -> {#3 &}, Mesh -> Length[enedata], 
                                              ViewPoint -> {0, 0, 100}, ImageSize -> OptionValue["imagesize"]],
                                        ListPointPlot3D[{data2d[[#1]], data2d[[#2]], enedata}\[Transpose],
                                               PlotStyle -> {White, PointSize[Large]},
                                               AxesLabel -> {ToString[fullvars[[#1]],StandardForm], ToString[fullvars[[#2]],StandardForm], ""},
                                               PlotRange -> All,
                                               ViewPoint -> {0, 0, 100}, ImageSize -> OptionValue["imagesize"]],
                                        ParametricPlot3D[BSplineFunction[{data2d[[#1]],data2d[[#2]],enedata}\[Transpose]][s], {s, 0, 1}, 
                                               PlotStyle -> Directive[White, Thick],
                                               PlotRange -> All,
                                               ViewPoint -> {0, 0, 100}, ImageSize -> OptionValue["imagesize"]]],
                                   ImageSize -> OptionValue["imagesize"]] 
                                              & @@@ Subsets[Table[If[Norm[data2d[[i]]] == 0, ## &[], i], {i, Length[fullvars]}], {2}]];
              If[OptionValue["plot"], Print[Grid@Partition[plt, UpTo[4]]]];
              {t[[1, 2]], plt}, {t, out[[2;;]]}];

  If[OptionValue["pdf"],
     Export[dir0<>"/ts_info.pdf", CreateDocument[Grid[Join[{{Graphics@Text[Style["ts-"<>ToString[#[[1]]], Red, Bold, 12]], First[#[[2]]]}},
                                                               Partition[#[[2]][[2;;]], UpTo[3]]], Frame -> All, FrameStyle -> Red] & /@ pdf,
                                          PageBreakBelow -> True, Visible -> False]]];


  Return[out]
]

PlanFittingOld[invariants_, ts_, vars_, svars_] := Module[{p, m, invdict, t, i, models, tcharacter, characters, modes, plan, fullvars},
  fullvars = Flatten[{vars, svars}];
  models = Flatten[SortInvariants[#, fullvars] &/@ Reverse@GatherBy[Flatten[invariants], StrainInvQ[#, svars[[1, 1]]] &], 1];
  Print["Number of invariants :" <> ToString[Length[Flatten@models]]];
  Print["Groups of models :" <> ToString[Length[models]]];
  
  plan = Table[characters = InvariantCharacter[m[[1]], fullvars];
               {m, Table[tcharacter=Sign@Flatten[Total[Abs[t\[Transpose][[5 ;; 6]]\[Transpose]]]];
                   If[MemberQ[characters, tcharacter], t, ##&[]], {t, ts}]}, {m, models}];
  Return[plan]
]

GetPolynomialTrainingSet[polynomial_, tsdata_, OptionsPattern[{"chop" -> 0.09999}]] := Module[{out, m, t, characters, tcharacter, p, tsp1, dft0},
  tsp1 = Select[tsdata, (Flatten[Sign@Abs@Chop[T5Data[First[#]], OptionValue["chop"]]] T5Vars[First[#]] === T5Vars[First[#]]) &];
  m = If[Head[polynomial]===Plus, polynomial[[1]], polynomial];
  characters = InvariantCharacter[m, Flatten@TSVars[First@tsdata]];
  out = Table[dft0 = Flatten[{Chop[T5Data[First[t]], OptionValue["chop"]], Chop[T6Data[First[t]], OptionValue["chop"]/10]}];
              tcharacter = Sign@Abs[dft0];
              If[MemberQ[characters, tcharacter], t, ## &[]], {t, tsdata}];
  If[out === {}, out = tsp1];
  Return[out]
]

PlanFitting[invariants_, modelfixed_, ts_, vars_, svars_, uvars_, OpMat_, TRules_, lvar_, mvar_, ordercut_, OptionsPattern[{"optimization" -> False, "call" -> 0, "fontsize" -> 12, "table" -> True}]] := Module[{p, ip, invdict, v, v1, v2, t, i, j, m, models, tcharacter, vcharacter, characters, modes, plan0, plan1, plan2, plan3, plan4, fullvars, minimums, color, simplifytab, mim, tmask, dependency, PEM, sol, PrimaryPolynomial, SlavePolynomial, HiddenPolynomial, tmp, sorted, len, foundQ, PickQ, tab, prefix, polynomials, tsdata, tag, slavets, OptPlan, mimlist, tmpmodel, initmodel, optmodel, submodel, sinv, pinv, spinv, inv2befit, EquivalentPath, dft, DependencyGraph, gdata, varsub, nvars0, VarTally}, 

  nvars0 = Length[Flatten[Append[vars, uvars]]];
  
  EquivalentPath[dft_] := DeleteDuplicates@Flatten[Table[
      tmask = Flatten[If[Norm[#] != 0, ConstantArray[1, Length[#]], Chop@#] & /@ Flatten[CloneReshape2D[Append[vars, uvars], #] &/@ Partition[T5Data[Last[t]], Length[Flatten[Append[vars, uvars]]]], 1]];
      mim = Sign@Chop[Mean[T5Data[#] tmask & /@ t], 0.1];
      varsub = Thread[T5Vars[Last[t]] -> Flatten[Round[m].# &/@ Partition[mim, nvars0 ]]];
      Sub2List[#][[2]] &/@ SwapVarSub[varsub], {t, dft}, {m, OpMat}], 2];

  tag = ToString[OptionValue["call"]];
  simplifytab = {-1 -> Style["-1", Bold, Blue], 1 -> Style["1", Bold, Red], 0 -> "-"};
  fullvars = TSVars[First@ts];
  VarTally = Append[Flatten[ConstantArray[Tally[VarHead[#] & /@ Flatten[Append[vars, uvars]]]\[Transpose][[2]], (Length[Flatten@fullvars] - 6)/nvars0]], 6];

  prefix = (Flatten[#] & /@ (modelfixed\[Transpose]))\[Transpose];
  inv2befit = Complement[Flatten[invariants], Transpose[prefix][[2]]];

  sinv = Select[inv2befit, And @@ XInvQ[#, {svars[[1, 1]]}, ContainsExactly] &];
  pinv = Select[inv2befit, And @@ XInvQ[#, {svars[[1, 1]]}, ContainsNone] &];
  spinv = Complement[inv2befit, Join[sinv, pinv]];

  models = Flatten[SortInvariants[#, Flatten@fullvars] & /@ {sinv, pinv, spinv}, 1];

  If[OptionValue["table"], Print["Number of invariants : " <> ToString[Length[Flatten@models]]]];
  
  (* polynomials <-> trainingset by First *)
  plan0 = Table[{m, GetPolynomialTrainingSet[m[[1]], ts]}, {m, models}];
  (* polynomials by minimums *)
  plan1 = Gather[plan0, IntersectingQ[Abs@EquivalentPath[#1[[2]]], Abs@EquivalentPath[#2[[2]]]] &];

  dependency = Table[PEM = Flatten[{Sub2List[#5][[2]], Sub2List[#6][[2]], #4}] & @@@ Flatten[plan1[[mim]]\[Transpose][[2]], 2];
                     (* polynomial to be fitted /  starting seed *)
                     PrimaryPolynomial = Flatten[plan1[[mim]]\[Transpose][[1]]];
                     (* primary polynomials couple to other modes *)
                     HiddenPolynomial = Table[vcharacter = InvariantCharacter[v, Flatten@fullvars];
                                              tcharacter = Sign@Abs@PEM[[;; , 1 ;; -2]];
                                              PickQ = FreeQ[Flatten@Table[v1 . v2 == Total[v1], {v1, vcharacter}, {v2, tcharacter}], True];
                                              If[PickQ, v, ## &[]], {v, Flatten@invariants}];
                     sol = Quiet[LinearModelFit[PEM, Complement[Flatten@invariants, HiddenPolynomial], Flatten@fullvars, 
                                                IncludeConstantBasis -> False][{"BestFitParameters", "BasisFunctions"}]]\[Transpose];
                     (* couplings inside the primary polynomial *)
                     SlavePolynomial = Complement[If[! MemberQ[PrimaryPolynomial, #2], #2, ## &[]] & @@@ sol, Transpose[prefix][[2]]];
    {Complement[DeleteDuplicates[First@First@Position[Table[Flatten[plan1[[i]]\[Transpose][[1]]], {i, Length@plan1}], #, {2}] & /@ SlavePolynomial], {mim}], 
    mim, 
    SlavePolynomial}, {mim, Length@plan1}];

  gdata = GroupBy[Flatten@Table[If[v[[1]] != {}, DirectedEdge[#, v[[2]]] & /@ (v[[1]]), {v[[2]]}], {v, dependency}], Head];
  DependencyGraph = Graph[If[MemberQ[Keys@gdata, Integer], gdata[Integer], {}], 
                          If[MemberQ[Keys@gdata, DirectedEdge], gdata[DirectedEdge], {}], VertexLabels -> "Name"];


  Print[Grid[dependency, Frame -> All]];
  Print[Grid[{Prepend[EquivalentPath[#[[1, 2]]], Flatten[#\[Transpose][[1]], 1]]} & /@ plan1, Frame -> All]];
  Print[Graph[DependencyGraph, GraphLayout -> "HighDimensionalEmbedding",
                               EdgeStyle -> Thickness[0.015],
                               EdgeShapeFunction -> GraphElementData[GraphElementData["EdgeShapeFunction", "FilledArrow"][[4]], "ArrowSize" -> 0.12], 
                               VertexSize -> Small,
                               VertexLabels -> (# -> Placed[Style[FromCharacterCode[9311 + #], 48], Center] & /@ Range[Length@dependency])]];

  sorted = Reverse[ConnectedComponents[DependencyGraph]];

  plan2 = Table[tmp = Transpose[Flatten[Part[plan1, Sort@v], 1]];
                Append[tmp, Complement[Flatten[Part[dependency\[Transpose][[3]], Sort@v]], Flatten[tmp[[1]]]]], {v, sorted}];

  (*tmp = dependency;
   sorted = Cases[tmp, {{}, __, {}}];
   tmp = DeleteCases[tmp, {{}, __, {}}];
   While[tmp != {}, 
         len = Length[tmp];
         foundQ = False;
         i = 1;
         While[! foundQ && i <= len, 
               If[ContainsAll[sorted\[Transpose][[2]], First[tmp[[i]]]], 
                  AppendTo[sorted, tmp[[i]]]; tmp = Drop[tmp, {i}]; foundQ = True, 
                  i = i + 1]];
         If[i > len && ! foundQ, Print["dependency error!"]; Abort[]];];

   plan2 = Append[plan1[[#2]]\[Transpose], #3] &@@@ sorted;*)
  
  plan3 = Table[polynomials = Flatten[plan2[[ip]][[1]]];
                SlavePolynomial = If[NumPolynomialVar[#]>1, #, ##&[]] &/@ Complement[plan2[[ip]][[3]], Keys[Model2Association[modelfixed]]];
                slavets = Flatten[Table[GetPolynomialTrainingSet[m, tsdata], {m, SlavePolynomial}], 1];
                tsdata = Join[Flatten[plan2[[ip]][[2]], 1], slavets];
                submodel = SortModelData[DeleteCases[{Transpose[MapIndexed[{ToExpression["PLAN"<>tag<>"PEM"<>ToString[ip]<>"Primary"<>ToString[First@#2]], #1}&, 
                                                             Complement[polynomials, Flatten@prefix]]]}, {}], Flatten@fullvars];
                {submodel,
                 DeleteCases[{(If[MemberQ[polynomials, #[[2]]], #, ## &[]] & /@ prefix)\[Transpose]}, {}], 
                 {},
                 DeleteCases[{Transpose[MapIndexed[{ToExpression["PLAN"<>tag<>"PEM"<>ToString[ip]<>"Slave"<>ToString[First@#2]], #1}&, SlavePolynomial]]}, {}],
                 tsdata}, {ip, Length@plan2}];

  mimlist = {};
  OptPlan = If[OptionValue["optimization"],
     initmodel = {{ConstantArray[2, Length@Flatten[invariants]], Flatten[invariants]}};
     slavets = ExtractMinimumTS[ts, TRules, OpMat];
     mimlist = Flatten[Sign@T56Data[Last[#]]] &/@ slavets;
     SlavePolynomial = DeleteDuplicates@Flatten[DeleteCases[Table[Flatten[GetActivePolynomial[initmodel, vars, svars, m][[2 ;;]], 1]\[Transpose], {m, mimlist}], {}]\[Transpose][[2]]];
     optmodel = GetOptModel[initmodel, vars, svars, ts, ordercut, OpMat, TRules, lvar, mvar];
     tsdata = Join[slavets, Flatten[GetPolynomialTrainingSet[#, ts] &/@ (Flatten[optmodel, 1][[2]]), 1]];
     {{{}, 
       {},
       optmodel,
       DeleteCases[{Transpose[MapIndexed[{ToExpression["OPT"<>tag<>"Slave"<>ToString[First@#2]], #1}&, SlavePolynomial]]}, {}],
       tsdata}}, {}];

  plan4 = Join[plan3, OptPlan];

  mimlist = Append[Table[Sort[Sign@T56Data[Last@#] &/@ ExtractMinimumTS[Flatten[p[[2]],1], TRules, OpMat, "type"->"StrainFree"]], {p, plan2}], mimlist];

  tab = {};
  If[OptionValue["table"], Print["Number of Potential Energy minimums: " <> ToString[Length[plan4]]]];
  If[OptionValue["table"],
     Do[minimums = mimlist[[i]];
        AppendTo[tab, {Style["Potential Energy minimum " <> ToString[i] <> ":", Bold, Black], Style["Polynomial Group " <> ToString[i] <> ":", Black, Bold]}];
        AppendTo[tab, {DomainGrid[(Prepend[#, "ene"] &/@ minimums) /. simplifytab, fullvars],
                            If[MemberQ[Flatten@prefix,#2],Style[#2,Darker@Gray],#2] &@@@ Model2List@SortModelData[Fold[JoinModels, plan4[[i]][[1;;3]]],Flatten@fullvars]}], {i, Length[plan4]}];
     Print[Grid[tab, ItemSize -> Automatic, Frame -> All, Alignment -> {Center, Center}]]];

  OptPlan = If[OptionValue["optimization"],
               initmodel = {{ConstantArray[2, Length@Flatten[invariants]], Flatten[invariants]}};
               mimlist = Flatten[Sign@Last[#][[5;;6]]] &/@ ExtractMinimumTS[ts, TRules, OpMat, "type"->"all"];
               SlavePolynomial = DeleteDuplicates@Flatten[DeleteCases[Table[Flatten[GetActivePolynomial[initmodel, vars, svars, m][[2 ;;]], 1]\[Transpose], {m, mimlist}], {}]\[Transpose][[2]]];
               optmodel = GetFixMimModel[initmodel, vars, 6, ts, OpMat, TRules, lvar, mvar];
               tsdata = Flatten[GetPolynomialTrainingSet[#, ts] &/@ Join[Keys[Model2Association[optmodel]], SlavePolynomial], 1];
               {{{},
                 {},
                 optmodel,
                 DeleteCases[{Transpose[MapIndexed[{ToExpression["OPT"<>tag<>"Slave"<>ToString[First@#2]], #1}&, SlavePolynomial]]}, {}],
                 tsdata}}, {}];

  Return[Join[plan4, OptPlan]]
]


LinvariantFitOld[invariants_, vars_, svars_, order_, ts_, spgmat_, ws_ : None, OptionsPattern[{"orderby" -> "AICc", "round" -> 10^-6.0, "offset" -> None, "alat" -> 1.0, "FontSize" -> 12}]] := Module[{alat, NumBasis, models, modelorders, OffsetFunc, fits, TopModels, BoundRadius, tsdata, seeds, plan, measure, x, i, m, p, t, mf, gc, log, out, phaseweight, minimumtest, sol, modelfitted, modeldata, tsinfo, fullvars, ibg, checkQ},
  fullvars = Flatten[{vars, svars}];
  alat = OptionValue["alat"];
  NumBasis = Length[vars];
  OffsetFunc = If[OptionValue["offset"] === None, None, (Total[Times @@ OptionValue["offset"]] /. Thread[fullvars -> {##}] &)];
  measure = Which[OptionValue["orderby"] === "AIcc", 4, 
                  OptionValue["orderby"] === "BIC", 5, 
                  OptionValue["orderby"] === "AdjustedRSquared", 6, 
                  OptionValue["orderby"] === "RSquared", 7, 
                  True, 4];
  
  modeldata = If[OptionValue["offset"] === None, 
                 Flatten[invariants[[#]] & /@ Range[order]], 
                 Delete[Flatten[invariants[[#]] & /@ Range[order]], Position[Flatten[invariants[[#]] & /@ Range[order]], #][[1]] & /@ (OptionValue["offset"][[2]])]];
  seeds = If[ws === None, {vars /. {Subscript[x_, i_] -> 0}}, RoundPhase[#] & /@ (ws\[Transpose][[2]])];
  out = If[OptionValue["offset"] === None, {{{}, {}}}, {OptionValue["offset"]}];

  plan = PlanFitting[modeldata, vars, svars, spgmat, ts, seeds];
  Print["Checking training set."];
  checkQ = False;
  Do[If[p[[2]]=={}, checkQ=True;Print["Error! Trainingset is not complete:"];Print[p[[1]]]], {p, plan}];
  If[checkQ, Abort[]];
  
  TopModels = Table[models = p[[1]];
                    tsdata = Flatten[{(#6 + #7)/alat, #8, #5}] & @@@ (p[[2]]);
                    tsinfo = Join[{Grid[{{"Number: ", Length[tsdata]}}]}, 
                                   Grid[Join[{{#1[[1]]<>": "<> ToString[#1[[2]]], SpanFromLeft}}, {#1[[3]]}, #2], 
                                        Background -> {None, Join[{Yellow}, {Pink}, Table[LightBlue, {Length[#2]}]]}, 
                                        Alignment -> {Center, Center}] & @@@ Normal[GroupBy[{{#2, #3, RoundPhase[#6]}, RoundPhase[#7]} & @@@ (p[[2]]), First -> Last, DeleteDuplicates]]];
    
    Quiet[phaseweight = If[ws === None, Normalize[1.0 & /@ tsdata], Expand[Normalize[1.0 + Sum[gc[[3, 2]] Exp[-(Norm[(#[[1 ;; NumBasis]] - gc[[2]])/alat]^2/(2.0*gc[[3, 1]]^2))], {gc, ws}] & /@ tsdata]]]];

    mf = Quiet[LinearModelFit[tsdata, models, fullvars, IncludeConstantBasis -> False, Weights -> phaseweight, LinearOffsetFunction -> OffsetFunc][{"AdjustedRSquared", "RSquared", "BestFitParameters", "BestFit"}]];
    
    (*If[NumPolynomialVar[models[[1]]]\[Equal]1,
    BoundRadius=Max[#]&/@partitionBy[(Max[#]&/@Abs[
    tsdata\[Transpose][[1;;NumBasis]]]),Length[#]&/@vars];
    BoundRadius=Flatten[Thread[ToExpression[
    GetVariationVar[#1]]\[Rule]ConstantArray[#2,
    Length[#1]]]&@@@({vars,BoundRadius}\[Transpose])];
    modelorders=PolynomialOrder[models,Flatten[vars]];
    Do[If[modelorders[[i]]\[Equal]order&&mf[[3,i]]<0.,models[[i]]=- 
    1models[[i]]+BoundHighOrderTerm[invariants,models[[i]],Flatten[
    vars],BoundRadius]],{i,Length@models}]];*)
    
    fits = Join[{Length@models + Length[Flatten[out\[Transpose][[1]]]]}, {Grid[{Join[Flatten[(Times @@ #) & /@ out], models mf[[3]]]}\[Transpose], ItemSize -> Full]}, mf];
    
    AppendTo[out, {Chop[fits[[5]], OptionValue["round"]], models}];
    modelfitted = Total[(Times @@ #) & /@ out, 2];
    OffsetFunc = (modelfitted /. Thread[fullvars -> {##}] &);
    
    minimumtest = DeleteDuplicates[
                    Flatten[Table[sol = Quiet@FindMinimum[modelfitted, {fullvars, m Flatten[{gc[[2]],{0.,0.,0.,0.,0.,0.}}]}\[Transpose]]; 
                      If[gc[[3, 2]] != 0, {{(vars /. sol[[2]]), svars /. sol[[2]]}, {Norm[gc[[2]] - (vars /. sol[[2]])], gc[[1]] - fits[[6]] /. Thread[vars -> gc[[2]]] /. sol[[2]]}}, ## &[]], {gc, ws}, {m, {1., 2.}}], 1], 
                    Chop[Norm[Flatten[#1[[1]] - #2[[1]]]], OptionValue["round"]] == 0. &]\[Transpose];
    Join[fits[[;; -3]], {tsinfo}, {Grid[{NumberForm[Grid[#, ItemSize -> Full], {6,6}]} & /@ (minimumtest[[1]]), Background -> {None, Table[If[OddQ[ibg], LightGray, White], {ibg, Length[minimumtest[[1]]]}]}, ItemSize->Full]}, {MatrixForm[minimumtest[[2]]]}], {p, plan}];
  
  log = Style[# /. x_Real :> Chop[x, OptionValue["round"]], OptionValue["FontSize"]] &@ Grid[{{"Length", "BestFit", "\!\(\*SuperscriptBox[\(R\), \(2\)]\)", "\!\(\*SuperscriptBox[\(R0\), \(2\)]\)", "trainingset", "Minimums", "Distances"}, ## & @@ TopModels}, Dividers -> All];
  Print[log];
  Return[out]
]

LinvariantFitOld[dir0_, invariants_, order_, vars_, svars_, ts_, TRules_, isovars_, OptionsPattern[{"exinvariants" -> {}, "round" -> 10^-6.0, "offset" -> None, "alat" -> 1.0, "FontSize" -> 12, "constraints" -> 10.0, "imagesize" -> 200, "pdf" -> True, "plot" -> False, "interpolate" -> 3, "trim" -> 0.05}]] := Module[{alat, NumBasis, models, modelorders, OffsetFunc, fits, TopModels, BoundRadius, tsdata, seeds, plan, x, i, m, p, t, v, mfdata, gc, log, out, phaseweight, minimumtest, sol, init, modelfitted, modeldata, tsinfo, fullvars, ibg, checkQ, character, data2d, enemodel, enedata, cplt, plt, lvar, mvar, ExplosionList, ExInvariants, AddInvariants, tm, ExplosionQ, AllInvariants, tmpmodel, exorder, FoundQ, convergency, trim, rsquare, coeff},
  {lvar, mvar} = isovars;
  fullvars = Flatten[{vars, svars}];
  alat = OptionValue["alat"];
  trim = OptionValue["trim"];
  NumBasis = Length[Flatten@vars];
  OffsetFunc = If[OptionValue["offset"] === None, 
                  None, 
                  (Total[Times @@ OptionValue["offset"]] /. Thread[fullvars -> {##}] &)];
  out = If[OptionValue["offset"] === None, {{{}, {}}}, {OptionValue["offset"]}];
  
  AllInvariants = DeleteDuplicates[Flatten[{invariants, OptionValue["exinvariants"]}]];
  modeldata = If[OptionValue["offset"] === None, 
                 AllInvariants, 
                 Complement[AllInvariants, OptionValue["offset"][[2]]]];
  
  plan = PlanFitting[modeldata, ts, vars, svars];
  Print["Checking training set."];
  checkQ = False;
  Do[If[p[[2]] == {}, checkQ = True; Print["Error! Trainingset is not complete:"]; Print[p[[1]]]], {p, plan}];
  If[checkQ, Abort[], Print["Training set check passed"]];
  
  TopModels = {};
  AddInvariants = {};
  Do[ExplosionQ = False;
     models = p[[1]];
     character = InvariantCharacter[models[[1]], fullvars];
     seeds = Join[Flatten[FetchSeeds[#, fullvars, 2, "zero" -> False] & /@ models, 1],
                  DeleteDuplicates[Flatten[#[[-1]][[5 ;; 6]]] & /@ ts, Chop[Norm[#1 - #2] Norm[#1 + #2], OptionValue["round"]] == 0. &]];
     tsdata = Select[Reverse[(Flatten[{#5/alat, #6, #4}] & @@@ Flatten[p[[2]], 1])], Sign[Abs[#[[1 ;; Length[fullvars]]]]] == character[[1]] &];
     tsinfo = Grid[{Grid[{{Style["ts-" <> ToString[#[[1, 2]]] <> ":", Red, Bold], SpanFromLeft}, {"Number: " <> ToString[Length[tsdata]]}}, 
                         Alignment -> {Center, Center}] & /@ (p[[2]])}];
   
     {rsquare, convergency} = Table[Quiet[LinearModelFit[tsdata[[1 ;; t]], models, fullvars, 
                                IncludeConstantBasis -> False, 
                                LinearOffsetFunction -> OffsetFunc][{"RSquared", "BestFitParameters"}]], 
                                {t, Length[models], Length@tsdata}]\[Transpose];

     coeff = Table[TrimmedCoefficient[convergency\[Transpose][[t]], trim, 0.1], {t, Length@models}];

     cplt = ListPlot[convergency\[Transpose], 
                     PlotLabel -> NumberForm[models, {3, 3}], PlotRange -> All, 
                     PlotMarkers -> Automatic, Frame -> True, 
                     ImagePadding -> {{40, 1}, {40, 1}}, 
                     Epilog -> Join[{Green, Thick, Dashed, Line[{{1, TrimmedMean[#, trim]}, 
                                                                 {Length[convergency], TrimmedMean[#, trim]}}]} & /@ (convergency\[Transpose]), 
                                    {Red, Thick, Dashed, Line[{{1, Sort[#, Less][[1 + Floor[OptionValue["trim"][[1]] Length[#]]]]}, 
                                                               {Length[convergency], Sort[#, Less][[1 + Floor[trim[[1]] Length[#]]]]}}]} & /@ (convergency\[Transpose]), 
                                    {Red, Thick, Dashed, Line[{{1, Sort[#, Less][[Length[#] - Floor[trim[[2]] Length[#]]]]}, 
                                                               {Length[convergency], Sort[#, Less][[Length[#] - Floor[trim[[2]] Length[#]]]]}}]} & /@ (convergency\[Transpose]), 
                                    Text[ToString[TrimmedMean[#, trim]], {Length[convergency]/2, TrimmedMean[#, trim]}] & /@ (convergency\[Transpose])]];
   
     tmpmodel = {Chop[coeff, OptionValue["round"]], models};
     coeff = If[MemberQ[OptionValue["exinvariants"], #2] && #1 < 0, -#1, #1] & @@@ (tmpmodel\[Transpose]);
     fits = Grid[Join[{Range[Length@models + Length[Flatten[out\[Transpose][[1]]]]]}, 
                       Join[Flatten[#\[Transpose] & /@ out, 1], {coeff, models}\[Transpose]]\[Transpose]]\[Transpose], 
                 ItemSize -> Full, Alignment -> Left, Background -> {{2 -> Lighter[Green, 0.5], 3 -> Lighter[Cyan, 0.5]}, None}];
     AppendTo[out, {Chop[coeff, OptionValue["round"]], models}];
     modelfitted = Total[(Times @@ #) & /@ out, 2];
     OffsetFunc = (modelfitted /. Thread[fullvars -> {##}] &);
   
     data2d = tsdata[[;; , 1 ;; Length[fullvars]]];
     enedata = tsdata[[;; , -1]];
     enemodel = modelfitted /. Thread[fullvars -> #] & /@ data2d;
     plt = Flatten[{Table[If[Norm[tsdata[[;; , i]]] != 0, 
                             ListPlot[{{tsdata[[;; , i]], enedata}\[Transpose], {tsdata[[;; , i]], enemodel}\[Transpose]}, 
                                      PlotRange -> All, 
                                      PlotStyle -> {Black, Red}, 
                                      ImageSize -> OptionValue["imagesize"], Frame -> True, 
                                      GridLines -> Automatic, PlotMarkers -> Automatic, 
                                      PlotLabel -> models, 
                                      FrameLabel -> {"("<>ToString[i]<>") "<>ToString[fullvars[[i]],StandardForm]<>" (\[Angstrom])","energy (eV)"}, 
                                      ImagePadding -> {{60, 10}, {40, 10}}], ## &[]], {i, NumBasis}],
                    Table[If[Norm[tsdata[[;; , NumBasis + i]]] != 0, 
                             ListPlot[{{tsdata[[;; , NumBasis + i]], enedata}\[Transpose], {tsdata[[;; , NumBasis + i]], enedata}\[Transpose]}, 
                                      PlotRange -> All, 
                                      PlotStyle -> {Black, Red}, 
                                      ImageSize -> OptionValue["imagesize"], Frame -> True, 
                                      GridLines -> Automatic, PlotMarkers -> Automatic, 
                                      PlotLabel -> models, 
                                      FrameLabel -> {"("<>ToString[i]<>") "<>ToString[fullvars[[i]],StandardForm]<>" (\[Angstrom])","energy (eV)"}, 
                                      ImagePadding -> {{60, 10}, {40, 10}}], ## &[]], {i, 6}]}];
   
     minimumtest = DeleteDuplicates[Flatten[Table[init = {fullvars, m s, ConstantArray[-OptionValue["constraints"], Length@fullvars], 
                                                                         ConstantArray[OptionValue["constraints"], Length@fullvars]}\[Transpose]; 
                                                  sol = Quiet@FindMinimum[modelfitted, init];
                                                  Prepend[fullvars /. (sol[[2]]), sol[[1]]], {s, seeds}, {m, {0.5, 3.}}], 1], 
                                    Chop[Norm[Flatten[#1[[2;;]] - #2[[2;;]]]], OptionValue["round"]] == 0. &];
   
     tm = Join[{fits}, {Grid[{{cplt}, {NumberForm[Grid[convergency],{4,4}]}}]}, {tsinfo}, {NumberForm[Grid[Prepend[minimumtest, Join[{"ene"}, fullvars]], 
        Background -> {Flatten[{LightBlue, Table[ConstantArray[If[OddQ[v], Yellow, Pink], Length[Append[vars, svars][[v]]]], {v, Length@Append[vars, svars]}]}], 
                       Prepend[Table[If[OddQ[ibg], White, None], {ibg, Length[minimumtest]}], Cyan]}, 
        ItemSize -> Full], {4, 4}]}, {Mean[rsquare]}, {Grid[{plt}]}];
     AppendTo[TopModels, tm];
   
     ExplosionList = DeleteDuplicates[Table[If[MemberQ[Abs[minimumtest[[i]]], 10.], 
                                               character[[1]]*minimumtest[[i]][[2;;]], 
                                               ## &[]], {i, Length[minimumtest]}]];
     If[ExplosionList != {}, 
        ExplosionQ = True;
        exorder = order;
        FoundQ = False;
        While[! FoundQ,
              exorder = exorder + 2;
              ExInvariants = Complement[Flatten[GetInvariants[TRules,Variables[models]/.Var2Var[lvar,mvar,1],{exorder},"round"->10^-6]/.Var2Var[lvar,mvar,2]], OptionValue["exinvariants"]];
              AddInvariants = Select[DeleteDuplicates[Flatten@Table[If[MemberQ[InvariantCharacter[#, fullvars], Sign@Abs@ExplosionList[[i]]], #, ## &[]] & /@ ExInvariants, {i, Length@ExplosionList}]], AllTrue[PolynomialOrder[#, fullvars, "tot" -> False], EvenQ] &];
              FoundQ = (AddInvariants != {})];
              Print[Style["WARNING: explosion detected, please following the suggestions to add bounding sets:", Red, Bold]]; 
              Print[AddInvariants];
              Print[tsinfo];
              Print[convergency\[Transpose]];
              Print[cplt];
            (*Print[Grid@Partition[ListPlot[If[#2\[GreaterEqual]0,{#1,#2},##&[]]&@@@Table[{i,modelfitted/.Thread[fullvars\[Rule]i*Normalize[#]]},{i,0,10,0.01/alat}],PlotRange\[Rule]All,PlotMarkers\[Rule]Automatic,Frame\[Rule]True,PlotLabel\[Rule]NumberForm[DeleteCases[fullvars*Normalize[#],0.],{3,3}],ImagePadding\[Rule]{{40,1}, {40,1}}]&/@ExplosionList,UpTo[5]]];*)
              Break[]], {p, plan}];
  
  log = Style[# /. x_Real :> Chop[x, OptionValue["round"]], OptionValue["FontSize"]] &@ Grid[{{"BestFit", "convergency", "trainingset", "Minimums", "R02","plots"}, ## & @@ TopModels}, Dividers -> All, ItemSize -> Full];
  
  Print[log];
  
  Export[dir0 <> "/fitting_info.pdf", CreateDocument[Style[# /. x_Real :> Chop[x, OptionValue["round"]], OptionValue["FontSize"]] &@ Grid[{{"BestFit", "convergency", "trainingset", "Minimums", "R02", "plots"}, #}\[Transpose], ItemSize -> Full] & /@ TopModels, PageBreakBelow -> True, Visible -> False]];
  
  Return[{out, AddInvariants}]
]

FitGammaModel[dir0_, pos0_, invariants_, grid_, vars_, svars_, uvars_, ts_, OpMat_, TRules_, isovars_, OptionsPattern[{"gseneweight"->1, "mimweight"->1, "ahinvariants" -> {}, "exinvariants" -> {}, "ordercut"->8, "nsw"->10^6, "optimization" -> 0, "adjust"->0.1, "round" -> 10^-6.0, "modelfixed" -> {{{}, {}}}, "hmodel" -> {}, "unit" -> "atomic", "FontSize" -> 12, "constraints" -> 10.0, "imagesize" -> 200, "pdf" -> True, "check"->True, "plot" -> False, "interpolate" -> 3}]] := Module[{latt, sites, alat, ev2hartree, slavemodel, submodel, modelorders, offset, OffsetModel, OffsetHam, OffsetFunc, FitInfo, logdata, tsdata, seeds, plan, x, i, m, ip, p, t, t1, t2, v, log, tsout, modeladjusted, modelout, minimumtest, sol, init, model4fit, ModelFuncFitted, modeldata, tsinfo, fullvars, checkQ, character, data2d, enemodel, enedata, plt, lvar, mvar, ExplosionList, ExInvariants, infodata, ModelExplosionQ, AllInvariants, modeltmp, exorder, FoundQ, rsquare, coeff, tmp, FitData, modelfixed, submodelfixed, prefix, lfitvars, lfitpolynomials, constraints, BoundOrder, polyseeds, submodelbound, boundcount, lfitinit, PolyOrderTable, newplan, addts, tsused, poly, PolyCheckList, FoundExplosionQ, lfitadjust, AdjustRate, GF, m1, m2, m3, m4, optts, mimlist, optmodel, var, Ham, tssub, superts, ahinvariants, hmodel, HoppingFix},

  {lvar, mvar} = isovars;
  {latt, sites} = pos0;
  If[OptionValue["unit"]==="evangstrom", alat = 1.0; ev2hartree=1.0, alat=Mean[Norm[#] & /@ latt]; ev2hartree=1/27.2114];

  ahinvariants = OptionValue["ahinvariants"];

  {OffsetModel, HoppingFix} = If[Times@@grid == 1 || OptionValue["hmodel"] === {},
                                 {{},{}},
                                 offset = GroupBy[Flatten[HModelFolding[OptionValue["hmodel"], grid], 1], 
                                                  OnSiteExprQ[#[[2]]]&&XInvQ[#[[2]], {ToExpression["u"]}, ContainsNone] &];
                                 {If[MemberQ[Keys@offset, False], offset[False], {}], If[MemberQ[Keys@offset, True], offset[True], {}]}];

  OffsetHam = Expand[Total[Times @@ # & /@ OffsetModel]];

  modelfixed = {Transpose[DeleteDuplicates[Join[HoppingFix, Model2List[Expr2Gamma[OptionValue["modelfixed"], grid]]], #1[[2]] === #2[[2]] &]]};

  prefix = Model2List[modelfixed];
  AllInvariants = Expr2Gamma[DeleteDuplicates[Flatten[{invariants, OptionValue["exinvariants"], FoldHopping[ahinvariants, grid]}]], grid];
  
  modeldata = Complement[AllInvariants, Keys[Model2Association[modelfixed]]];

  superts = TrainingSet2Supercell[ts, grid];
  fullvars = TSVars[First@superts];

  plan = PlanFitting[AllInvariants, modelfixed, superts, vars, svars, uvars, OpMat, TRules, lvar, mvar, OptionValue["ordercut"], 
                     "table" -> True, "optimization" -> !(OptionValue["optimization"]==0)];
  
  Print[Style["Checking training set.", Black, Bold]];
  checkQ = False;
  Do[If[p[[5]] == {}, checkQ = True; Print[Style["Error! Trainingset is not complete:", Black, Bold]]; Print[p[[1]]]], {p, plan}];
  If[checkQ, Abort[], Print[Style["Training set check passed", Black, Bold]]];

  logdata = {};
  AdjustRate = ConstantArray[OptionValue["adjust"], Length@plan];
  If[OptionValue["optimization"] != 0, AdjustRate[[-1]] = OptionValue["optimization"]];
  modelout = {SortModelData[DeleteCases[{(If[MemberQ[Keys[Model2Association[Flatten[plan\[Transpose][[2]], 1]]], #2], ## &[], {#1, #2}] & @@@ prefix)\[Transpose]}, {}], fullvars]};
  tsout = {};
  
  Do[(* loop plans *)
     submodelbound = {};
     submodel = plan[[ip]][[1]];
     submodelfixed = plan[[ip]][[2]];
     optmodel = plan[[ip]][[3]];
     slavemodel = plan[[ip]][[4]];
     (*OffsetFunc = (modelfitted /. Thread[Flatten[fullvars] -> {##}] &);*)
  
     tsused = plan[[ip]][[5]];
     tsdata = Flatten[{Thread[Sub2List[#5][[1]] -> Sub2List[#5][[2]]/alat], #6, #4*ev2hartree}] & @@@ Flatten[tsused, 1];

     mimlist = {Last[#][[1]], 
                Last[#][[2]], 
                Last[#][[3]], 
                Last[#][[4]] ev2hartree, 
                Thread[Sub2List[Last[#][[5]]][[1]] -> Sub2List[Last[#][[5]]][[2]]/alat],
                Last[#][[6]]} & /@ ExtractMinimumTS[tsused, TRules, OpMat, "type" -> "all"];

     (*{rsquare, coeff} = If[submodel != {},
                           Quiet[LinearModelFit[tsdata, submodel, Flatten[fullvars], 
                                                IncludeConstantBasis -> False, 
                                                LinearOffsetFunction -> OffsetFunc][{"RSquared", "BestFitParameters"}]], 
                           {0, {}}];*)

     ModelExplosionQ = True;
     modeltmp = submodel;
     boundcount = 1;
     sol = {};
     addts = {};
     While[(* loop bounding explosion*)
       ModelExplosionQ,
       submodel = JoinModels[submodelbound, submodel];

       lfitinit = Join[Table[{poly[[1]], If[Head[poly[[1]]/.sol] === Symbol, 0, poly[[1]]/.sol], -Infinity, Infinity}, {poly, Model2List[submodel]}],
                       Table[coeff = First@First@Select[Model2List[Last@modelout], #[[2]] === poly[[2]]&];
                             {poly[[1]], 0, -Abs[coeff*AdjustRate[[ip]]], Abs[coeff*AdjustRate[[ip]]]}, {poly, Model2List[slavemodel]}],
                       Table[{poly[[1]], 1, 0, Infinity}, {poly, Model2List[optmodel]}]];

       Do[If[!MemberQ[Transpose[Flatten[tsused, 1]][[2]], First[t][[2]]], 
             AppendTo[tsdata, Flatten[{#5/alat, #6, #4*ev2hartree}]] &@@@ t;
             AppendTo[tsused, t]], {t, addts}];
       tsinfo = Grid[Partition[Style["ts-" <> ToString[#[[1, 2]]] <> " (" <> ToString[Length[#]] <> ")", Black, Bold] & /@ tsused, 1]];

       (* fitting by minimizing GF *)
       model4fit = Fold[JoinModels, {Last@modelout, submodelfixed, submodel, optmodel, slavemodel}];
       Ham = GetHam[model4fit] + OffsetHam;

       GF = Expand[Sum[tssub = Flatten[t1[[1;;-2]]];
                       (Last[t1] - Ham)^2 /. tssub, {t1, tsdata[[2;;]]}] 
                  + OptionValue["mimweight"]*Sum[tssub = Flatten[t2[[5;;6]]];
                             (OptionValue["gseneweight"]*(t2[[4]] - Ham)^2 + Sum[Expand[D[Ham, {var, 1}]]^2, {var, Flatten@fullvars}]) /. tssub, {t2, mimlist}]];

       If[submodel === {} && slavemodel === {} && optmodel === {},
          {rsquare, sol} = {0, {}},
          {rsquare, sol} = FindMinimum[GF, lfitinit, MaxIterations -> OptionValue["nsw"]];
          sol = Thread[Transpose[lfitinit][[1]] -> N@Round[Transpose[lfitinit][[1]] /. sol, OptionValue["round"]]]];

       Print[Style["\[FilledSquare] " <> "Fitting Potential Energy Minimum: PEM-"<>ToString[ip] 
                                      <> ", bounding trial step: " 
                                      <> ToString[boundcount-1] 
                                      <> ", add terms: "
                                      <> ToString[GetHam[submodelbound], StandardForm], Bold, Black]];

       (*{rsquare, sol} = If[submodel != {},
                           Quiet[NonlinearModelFit[tsdata, 
                                                   (GetHam[model4fit] + GetHam[slavemodel]), 
                                                   lfitinit,
                                                   Flatten[fullvars]][{"RSquared", "BestFitParameters"}]], 
                           {0, {}}];*)

       modeltmp = SortModelData[Fold[JoinModels, {submodelfixed, submodel, optmodel} /. sol], fullvars];
       modeladjusted = AdjustModel[Last@modelout, slavemodel /. sol];

       (* Interaction terms *)
       PolyCheckList = If[slavemodel === {},
                          SortInvariants[Keys@Model2Association[submodel], Flatten@fullvars],
                          SortInvariants[Join[Keys@Model2Association[submodel], Flatten[slavemodel,1][[2]]], Flatten@fullvars]];
       PolyCheckList = Flatten[Table[MaximalBy[poly, PolynomialOrder[#, Flatten[fullvars]] &], 
                                     {poly, PolyCheckList}]];
       PolyCheckList = If[!StrainInvQ[#, svars[[1, 1]]], #, ##&[]] &/@ PolyCheckList;

       (* take the first (ordered) explosion polynomial *)
       polyseeds = {};
       Do[polyseeds = If[ExplosionQ[Fold[JoinModels, {modeladjusted, modeltmp, {Transpose@CollectPolynomial[OffsetHam]}}], poly], {poly}, {}];
          If[polyseeds != {}, Break[]], {poly, PolyCheckList}];

       (* calculate energy invariants on the fly *)
       ExInvariants = If[polyseeds === {}, {}, Print[polyseeds];Abort[];GetBoundingInvariants[polyseeds[[1]], vars, svars, OptionValue["ordercut"], TRules, lvar, mvar]];

       If[ExInvariants === {} && polyseeds != {},
          {m1, m2, m3, m4} = First@Position[{submodel, slavemodel}, First@polyseeds];
          Which[m1 == 1, submodel[[m2,1,m4]] = 0, m1 == 2, slavemodel[[m2,1,m4]] = 0]];

       If[ExInvariants === {},
          ModelExplosionQ = False;
          submodelbound === {{}};
          addts = {},
          newplan = Transpose[PlanFitting[ExInvariants, modelfixed, ts, vars, svars, OpMat, TRules, lvar, mvar, OptionValue["ordercut"],
                                          "table"->False, "call"->ToString[ip]<>"bound"<>ToString[boundcount]]];
          submodelbound = Flatten[newplan[[1]],1];
          addts = Flatten[newplan[[5]], 1]];

       boundcount = boundcount + 1;
       If[boundcount>10,Break[]]
     ];

     (* loop bounding minimums *)

     (* build models and summarys*)

     AppendTo[modelout, JoinModels[modeladjusted, modeltmp, "simplify"->False]];
   
     FitInfo = Grid[MapIndexed[Which[#1[[1]]!=0&&MemberQ[Flatten[prefix], #1[[2]]], {#2[[1]], Style[#1[[1]], Darker@Gray], #1[[2]]}, 
                                     #1[[1]]!=0&&MemberQ[Flatten[slavemodel], #1[[2]]], {#2[[1]], Style[#1[[1]], Red], #1[[2]]},
                                     #1[[1]]!=0, {#2[[1]], Style[#1[[1]], Black], #1[[2]]},
                                     True, ##&[]] &, Model2List[Last@modelout]],
                    ItemSize -> Full, Alignment -> Left, Frame -> True, 
                    Background -> {{2 -> Lighter[Green, 0.5], 3 -> Lighter[Cyan, 0.5]}, None},
                    Dividers -> {{1 -> True}, Thread[(Accumulate[Length[DeleteCases[Transpose[#], {0., _}]] & /@ Last[modelout]] + 1) -> True]}];
   
     AppendTo[tsout, plan[[ip]][[4]]]; 
   
     ModelFuncFitted = GetHam[Last@modelout] + OffsetHam;
   
     plt = If[OptionValue["check"], CheckFitting[pos0, Keys@Model2Association[submodel], Append[Last@modelout, Transpose@CollectPolynomial[OffsetHam]], tsused, fullvars, "plot" -> False, "unit" -> OptionValue["unit"]], {}];

     (*seeds = Join[DeleteDuplicates@Fold[Join, FetchSeeds[#, Flatten[fullvars], 1, "zero" -> False] & /@ Keys@Model2Association[submodel]],
                  DeleteDuplicates[T56Data[Last@#] & /@ superts, Chop[Norm[#1 - #2] Norm[#1 + #2], OptionValue["round"]] == 0. &]];*)

     seeds = Join[Fold[Join, FetchSeeds[#, Flatten[fullvars], 1, "zero" -> False, "limit"->2^10] & /@ Keys@Model2Association[submodel]],
                  SparseArray@DeleteDuplicates[T56Data[Last@#] & /@ superts, Chop[Norm[#1 - #2] Norm[#1 + #2], OptionValue["round"]] == 0. &]];

     Print[Style["Number of Directions During Explosion Test: " <> ToString@Length@seeds, Black, Bold]]; 

     minimumtest = SortBy[DeleteDuplicates[Flatten[Table[init = {Flatten[fullvars], 
                                                                m s, 
                                                                ConstantArray[-OptionValue["constraints"], Length@Flatten[fullvars]],                       
                                                                ConstantArray[OptionValue["constraints"], Length@Flatten[fullvars]]}\[Transpose];
                                                        sol = Quiet@FindMinimum[ModelFuncFitted, init, MaxIterations->10000];
                                                        Prepend[Chop[Flatten[fullvars] /. (sol[[2]]), OptionValue["round"]], 
                                                                N@Round[sol[[1]], OptionValue["round"]]], {s, seeds}, {m, {0.5, 3.}}], 1], 
                                          (Chop[Norm[Flatten[#1[[2 ;;]] - #2[[2 ;;]]]], OptionValue["round"]] == 0. || Chop[#1[[1]]-#2[[1]]] == 0.)&], First, Greater];

     FitInfo = Grid[{{FitInfo}, {NumberForm[DomainGrid[Join[{First[#]}/(Times@@grid), #[[2;;]]] &/@ minimumtest, fullvars], {4, 4}]}}, Alignment -> Left];
   
     infodata = Join[{FitInfo}, 
                     {tsinfo}, 
                     {rsquare}, 
                     {Grid[Partition[plt, UpTo[10]]]}];
   
     AppendTo[logdata, infodata];
   
     ExplosionList = Append[DeleteDuplicates[Table[If[MemberQ[Abs[minimumtest[[i]]][[2 ;; Length[Flatten@fullvars] + 1]], 10.], 
                                                   If[Chop[Abs@# - 10.] == 0., 1, 0] & /@ (minimumtest[[i]][[2 ;; Length[Flatten@fullvars] + 1]]), 
                                                   ## &[]], {i, Length[minimumtest]}]],
                           ConstantArray[0,Length@Flatten@fullvars]]; 

     Print[Style["Explosion list:", Black, Bold]];
     Print[DomainGrid[(Prepend[#, "\[Infinity]"] &/@ ExplosionList) /. {0->"-"}, fullvars, "framecolor" -> If[Length@ExplosionList > 1, Red, Black]]], {ip, Length@plan}];
  
  log = Style[# /. x_Real :> Chop[x, OptionValue["round"]], OptionValue["FontSize"]] &@ Grid[{{"BestFit", "trainingset", "R02", "plots"}, ## & @@ logdata}, Dividers -> All, ItemSize -> Full];
  Print[log];
  
  If[OptionValue["pdf"],  Export[dir0 <> "/fitting_info.pdf", CreateDocument[Style[# /. x_Real :> Chop[x, OptionValue["round"]], OptionValue["FontSize"]] &@ Grid[{{"BestFit", "trainingset", "R02", "plots"}, #}\[Transpose], ItemSize -> Full] & /@ logdata, PageBreakBelow -> True, Visible -> False]]];

  Return[{modelout, tsout, OffsetHam}]
]

DomainGrid[domains_, vars_, OptionsPattern["framecolor"->Black]] := Module[{tmp, data, d, i, thead, VarTally},
  thead = Join[{"hartree/u.c.", "site"}, VarBare@Flatten[#@vars & /@ {First, Last}]];
  VarTally = Split[Flatten@VarHead[#@vars & /@ {First, Last}]];
  data = Flatten[Table[tmp = CloneReshape2D[vars, d[[2 ;;]]];
                       MapIndexed[Join[If[#2 == {1}, {First@d}, {""}], 
                                       {MatrixForm[{VarGrid[First[vars[[First@#2]]]]}]}, 
                                       #1, 
                                       If[#2 == {1}, Last@tmp, ConstantArray["", 6]]] &, tmp[[;; -2]]], {d, domains}], 1];
  Grid[Prepend[data, thead], 
       Frame -> True,
       FrameStyle -> OptionValue["framecolor"],
       Background -> {Flatten[{LightBlue, LightPurple, Table[ConstantArray[If[OddQ[i], Yellow, Pink], Length[VarTally[[i]]]], {i, Length@VarTally}]}],
                     Prepend[Flatten@Table[If[OddQ[i], ConstantArray[Gray, Length@tmp - 1], ConstantArray[None, Length@tmp - 1]], {i, Length[domains]}], Cyan]}, ItemSize -> Full]
]

GetHamOld[models_, vars_, svars_, OptionsPattern[{"subset"->{}}]] := Module[{out, subvars, fixsub, fullvars, polynomials, ham, s},
  fullvars = Append[vars, svars];
  ham = Total[(Times @@ #) & /@ models, 2];
  polynomials = Flatten[If[Head[#] === Plus, Level[#, 1], {#}] &/@ OptionValue["subset"]];
  subvars = If[polynomials === {}, {Flatten[fullvars]}, Variables[#] &/@ polynomials];
  out = Sum[fixsub = Thread[Complement[Flatten@fullvars, s] -> 0]; 
            ham /. fixsub, {s, subvars}];
  Return[out]
]

GetHam[models_, OptionsPattern[{"subset"->{}}]] := Module[{out, subvars, fixsub, fullvars, polynomials, ham, s},
  fullvars = If[models ==={}, {}, Variables[Model2List[models]\[Transpose][[2]]]];
  ham = Total[(Times @@ #) & /@ models, 2];
  polynomials = Flatten[If[Head[#] === Plus, Level[#, 1], {#}] &/@ OptionValue["subset"]];

  subvars = If[polynomials === {}, {fullvars}, Variables[#] &/@ polynomials];

  out = Sum[fixsub = Thread[Complement[fullvars, s] -> 0];
            ham /. fixsub, {s, subvars}];
  Return[out]
]


Expression2Model[ham_, fullvars_] := Module[{data, a, b},
  data = If[Head[ham] === Plus,
            Cases[ham, Times[a_, b_] /; And[AtomQ[a], ! AtomQ[b]] -> {a, b}],
            Cases[{ham}, Times[a_, b_] /; And[AtomQ[a], ! AtomQ[b]] -> {a, b}]];
  Return[SortModelData[DeleteCases[{data\[Transpose]}, {}], fullvars]]
]

GetActivePolynomial[modeldata_, vars_, svars_, mim_] := Module[{fullvars, PolynomialList, FullPolynomialList, i, out, models, tmpvars, subset},
  fullvars = Append[vars, svars];
  tmpvars = DeleteCases[Sign[Chop@mim[[1;;Length[Flatten@fullvars]]]]*Flatten[fullvars], 0];
  subset = Times @@ # & /@ Subsets[tmpvars, {1, Length@tmpvars}];
  models = Transpose[Flatten[#] &/@ (Expression2Model[GetHam[modeldata, "subset" -> subset], vars, svars]\[Transpose])];
  PolynomialList = If[Head[#2] === Plus, Level[#2, 1], {#2}] & @@@ models;
  FullPolynomialList = If[Head[#] === Plus, Level[#, 1], {#}] & /@ Flatten[modeldata\[Transpose][[2]]];
  out = Flatten[Table[If[IntersectingQ[FullPolynomialList[[i]], #], {i, Flatten[modeldata\[Transpose][[2]]][[i]]}, ## &[]], {i, Length@FullPolynomialList}] & /@ PolynomialList, 1];
  Return[GatherBy[SortBy[out, First], NumPolynomialVar[#[[2]]] &]]
]

UpdateModel[models_, ham_, vars_, svars_] := Module[{data, a, b, i, srt, modelsnew},
  modelsnew = models;
  data = Model2List@Expression2Model[ham, vars, svars];

  Do[srt = Position[models, data[[i, 2]]];
     If[srt==={},
        AppendTo[modelsnew, {{data[[i, 1]]}, {data[[i, 2]]}}],
        modelsnew[[First[srt][[1]], 1, First[srt][[3]]]] = data[[i, 1]]], {i, Length[data]}];
  Return[modelsnew]
]

TrimmedCoefficient[coeff_, trim_, tol_] := Module[{LinearInvQ, trimr, triml, trimmed, converged, convergeQ, n, coeffavg},
  trimr = Sort[coeff, Less][[Length[coeff] - Floor[trim[[2]] Length[coeff]]]];
  triml = Sort[coeff, Less][[1 + Floor[trim[[2]] Length[coeff]]]];
  trimmed = (If[triml <= # <= trimr, #, ## &[]] & /@ coeff);
  converged = Reverse[{Length[#], Mean[#]} & /@ Split[trimmed, Abs[#2/#1 - 1] < tol &]];
  convergeQ = False;
  n = 1;
  While[! convergeQ && n <= Length[converged], 
        coeffavg = If[converged[[n]][[1]]>=3, 
                      convergeQ=True;converged[[n]][[2]], 
                      convergeQ=False; 0]; 
        n++];
  coeffavg = If[convergeQ, coeffavg, Mean[converged\[Transpose][[2]]]];
  Return[coeffavg]
]

CheckMinimums[models_, seeds_, amp_, fullvars_, OptionsPattern[{"round" -> 10^-4, "nsw"->10000, "all"->True}]] := Module[{sol, init, s, minimums, simplifytab, a, tabdata, grid},
  
  grid = Length[DeleteDuplicates[#]] & /@ Transpose[VarGrid[#] & /@ Flatten[fullvars[[1 ;; 2]]]];

  simplifytab = DeleteDuplicates@Flatten[Table[{-1*a -> Style["-1", Bold, Blue], 1*a -> Style["1", Bold, Red], 0 -> "-"}, {a, amp}]];
  minimums = Table[init = {Flatten@fullvars, Round[a*s,1], -(Abs[amp]+1)*ConstantArray[1,Length[s]], (Abs[amp]+1)*ConstantArray[1,Length[s]]}\[Transpose];
                   sol = Quiet@FindMinimum[GetHam[models], init, MaxIterations -> OptionValue["nsw"]];
                   Prepend[# & /@ Chop@N@Round[Flatten@fullvars /. sol[[2]], OptionValue["round"]], Chop[sol[[1]]]/(Times@@grid)] /. simplifytab, {s, seeds}, {a, amp}];

  tabdata = Reverse@SortBy[DeleteDuplicates[Flatten[minimums, 1], (Chop[Norm[#1[[1]] - #2[[1]]], 10^-6] == 0) &], #[[1]] &];
  Print[DomainGrid[tabdata, fullvars]];

]

CheckFitting[pos0_, invariants_, model_, ts_, fullvars_, OptionsPattern[{"imagesize" -> 200, "col" -> 4, "plot" -> True, "plotrange" -> All, "unit" -> "atomic"}]] := Module[{inv, character, data2d, enemodel, s, m, i, t, plt, alat, ev2hartree, latt, sites},
  {latt, sites} = pos0;
  If[OptionValue["unit"]==="evangstrom", alat = 1.0; ev2hartree=1.0, alat=Mean[Norm[#] & /@ latt]; ev2hartree=1/27.2114];

  plt = Flatten@Table[inv = If[MatchQ[Head@m, Plus], m[[1]], m];
                      character = If[MemberQ[Variables@inv, #], 1, 0] & /@ Flatten[fullvars];
                      character = First@InvariantCharacter[inv, fullvars];
                      Table[If[Sign@Abs@Join[T5Data[First[ts[[t]]]], T6Data[First[ts[[t]]]]] == character,
                               data2d = Flatten[{Sub2List[#5][[2]]/alat, Sub2List[#6][[2]], #4*ev2hartree}] & @@@ (ts[[t]]);
                               enemodel = (GetHam[model] /. Thread[Flatten[fullvars] -> #[[1 ;; -2]]]) & /@ data2d;
                               Prepend[Table[If[Norm[data2d[[;; , i]]] != 0, 
                                                ListPlot[{{data2d[[;; , i]], data2d[[;; , -1]]}\[Transpose], {data2d[[;; , i]], enemodel}\[Transpose]}, 
                                                         PlotRange -> OptionValue["plotrange"], 
                                                         PlotStyle -> {Black, Red}, 
                                                         ImageSize -> OptionValue["imagesize"], Frame -> True, 
                                                         GridLines -> Automatic, PlotMarkers -> Automatic, 
                                                         PlotLabel -> Style[ToString[inv, StandardForm] <> " (ts-" <> ToString[First[ts[[t]]][[2]]] <> ")", Blue, 8], 
                                                         FrameLabel -> {"(" <> ToString[i] <> ") " <> ToString[Flatten[fullvars][[i]], StandardForm] <> " (\[Angstrom])", "energy (eV)"}, 
                                                         ImagePadding -> {{60, 10}, {40, 10}}], ## &[]], {i, Length[Flatten@fullvars]}], 
                                       If[Chop@Norm[Last[data2d][[1 ;; -2]]] != 0., 
                                          Plot[GetHam[model] /. Thread[Flatten[fullvars] -> s Normalize[Last[data2d][[1 ;; -2]]]], {s, 0, 1.5 Norm[Last[data2d][[1 ;; -2]]]}, 
                                               PlotStyle -> Red, 
                                               PlotRange -> OptionValue["plotrange"], 
                                               ImageSize -> OptionValue["imagesize"], Frame -> True, 
                                               GridLines -> Automatic, 
                                               PlotLabel -> Style[ToString[inv, StandardForm] <> " (ts-" <> ToString[t] <> ")", Blue, 8], 
                                               FrameLabel -> {ToString[inv, StandardForm] <> " (\[Angstrom])", "energy (eV)"}, 
                                               ImagePadding -> {{60, 10}, {40, 10}}], {}]], ## &[]], {t, Length@ts}], {m, If[ListQ[invariants], invariants, {invariants}]}];
                                               
  If[OptionValue["plot"], Print[Grid@Partition[plt, UpTo[OptionValue["col"]]]]];

  Return[plt]
]

SortModelData[modeldata_, vars_] := Module[{sortedinv, out, model},
  out = Table[sortedinv = Flatten[SortInvariants[model[[2]], Flatten@vars]];
              SortBy[model\[Transpose], Position[sortedinv, #[[2]]] &]\[Transpose], {model, modeldata}];
  Return[out]
]

JoinModels[m0_, m1_, OptionsPattern[{"simplify"->True}]] := Module[{m00, data0, data1, tmp, out, m, p1, p2, p3, p},
  m00 = AdjustModel[m0, m1];
  data0 = Flatten[#] & /@ (m00\[Transpose]);
  data1 = Flatten[#] & /@ (m1\[Transpose]);
  tmp = If[MemberQ[Flatten[data0], #2], ##&[], {#1, #2}] &@@@ (data1\[Transpose]);
  tmp = DeleteCases[Append[m00, Transpose[tmp]], {}];
  out = DeleteCases[If[OptionValue["simplify"], 
                       {Flatten[#] &/@ Transpose[tmp]}, 
                       tmp], {}];
  Return[out]
]

GetTopOrder[polynomials_, vars_, svars_] := Module[{fullvars, tab, out},
  fullvars = Append[vars, svars];
  tab = PolynomialOrder[#, Flatten@fullvars, "tot" -> False] & /@ polynomials;
  out = MapIndexed[If[(Length[Flatten@fullvars] - Count[#1, 0] == 1) && Max[tab[[;; , First@First@Position[#1, Max[#1]]]]] == Max[#1], {First@#2, Extract[polynomials, #2]}, ## &[]] &, tab];
  Return[out]
]

TopOrderQ[p_, polynomials_, vars_, svars_] := Module[{fullvars, tab, orders},
  fullvars = Append[vars, svars];
  tab = PolynomialOrder[#, Flatten@fullvars, "tot" -> False] & /@ polynomials;
  orders = PolynomialOrder[p, Flatten@fullvars, "tot" -> False];
  Max[orders]==Total[orders] && Max[orders]==Max[tab[[;; , First@First@Position[orders, Max@orders]]]]
]

Model2List[model_] := Module[{out},
  out = Transpose[Flatten[#] &/@ Transpose[model]];
  Return[out]
]

EnergyGainQ[model_, vars_, svars_, invariant_] := Module[{tmpvars, tmpmodel, saddle, xx, fullvars, sol, s, outQ, inv, DropDirections}, 
  fullvars = Append[vars, svars];
  inv = If[Head[invariant] === Plus, invariant[[1]], invariant];
  tmpvars = Variables[inv];
  tmpmodel = Chop[GetHam[model, "subset" -> {inv}]];
  saddle = NSolve[Table[D[tmpmodel, {xx, 1}] == 0, {xx, tmpvars}], tmpvars, Reals];
  DropDirections = If[(tmpvars Abs@Sign[tmpvars /. #] === tmpvars), #, ## &[]] & /@ saddle;
  outQ = ! DropDirections === {};
  Return[{outQ, DropDirections}]
]

ExplosionQ[model_, invariant_] := Module[{tmpvars, tmpmodel, saddle, xx, sol, s, outQ, inv, DropDirections},
  inv = If[Head[invariant] === Plus, invariant[[1]], invariant];
  tmpvars = Variables[inv];
  tmpmodel = Chop[GetHam[model, "subset" -> {inv}]];
  saddle = NSolve[Table[D[tmpmodel, {xx, 1}] == 0, {xx, tmpvars}], tmpvars, Reals];
  DropDirections = If[(tmpvars Abs@Sign[tmpvars /. #] === tmpvars), #, ## &[]] & /@ saddle;
  outQ = Or @@ Table[sol = Chop@Quiet@FindMinimum[tmpmodel, {#, 2 # /. s, -10000, 10000} & /@ tmpvars]; 
                     MemberQ[Abs[tmpvars /. sol[[2]]], 10000.], {s, DropDirections}];
  Return[outQ]
]

AdjustModel[model_, modeladjust_] := Module[{ModelAdjustList, out},
  ModelAdjustList = Model2List[modeladjust];
  out = If[ModelAdjustList === {},
           model,
           Table[Table[If[MemberQ[ModelAdjustList\[Transpose][[2]], m\[Transpose][[i, 2]]], 
                          {First@First@Select[ModelAdjustList, #[[2]] === m\[Transpose][[i, 2]] &] + m\[Transpose][[i, 1]], m\[Transpose][[i, 2]]}, 
                       m\[Transpose][[i]]], {i, Length[m\[Transpose]]}]\[Transpose], {m, model}]];
  Return[out]
]

GetBoundingInvariants[invariant_, vars_, svars_, cut_, TRules_, lvar_, mvar_, OptionsPattern[{"ExplosionQ"->True}]] := Module[{ExInvariantCandidates, ExInvCandidate, order, exorder,  fullvars, tmpvars},
  fullvars = Append[vars, svars];
  tmpvars = Variables@If[Head[invariant] === Plus, First[invariant], invariant];
  ExInvCandidate = {};
  order = PolynomialOrder[invariant, Flatten@fullvars];
  exorder = If[EvenQ[order], order + 2, order + 1];
  While[ExInvCandidate === {} && exorder <= cut,
        ExInvCandidate = If[OptionValue["ExplosionQ"], 
                            DeleteDuplicates@Flatten[GetInvariants[TRules, 
                                                                   {#} /. Var2Var[lvar, mvar, 1], 
                                                                   {exorder}, 
                                                                   "round" -> 10^-6,
                                                                   "eventerms" -> True,
                                                                   "check" -> False] /. Var2Var[lvar, mvar, 2] &/@ tmpvars],
                            Select[Flatten[GetInvariants[TRules,
                                                         tmpvars /. Var2Var[lvar, mvar, 1],
                                                         {exorder},
                                                         "round" -> 10^-6,
                                                         "eventerms" -> True,
                                                         "check" -> False] /. Var2Var[lvar, mvar, 2]],
                                  (Complement[InvariantCharacter[invariant, Flatten@fullvars], 
                                              InvariantCharacter[#, Flatten@fullvars]] === {})&]];
        exorder = exorder + 2];
  If[ExInvCandidate === {},
     Print[Style["  "<> ToString[invariant, StandardForm]
                     <> " unbounded "
                     <> "until order "
                     <> ToString[cut]
                     <> ", will be dropped!", Red]]];
  Return[ExInvCandidate]
]

GetOptModelOld[model_, vars_, svars_, tsused_, cut_, OpMat_, TRules_, lvar_, mvar_] := Module[{optts, m, mimlist, PolyCheckList, polyseeds, ExInvariants, submodel},
  optts = ExtractMinimumTS[tsused, TRules, OpMat];
  mimlist = Flatten[Sign@Last[#][[5 ;; 6]]] & /@ optts;

  ExInvariants = DeleteDuplicates@Flatten@Table[PolyCheckList = GetActivePolynomial[model, vars, svars, m];
                                                If[Length@PolyCheckList > 1, 
                                                   polyseeds = (Flatten[PolyCheckList[[2 ;;]], 1]\[Transpose][[2]]);
                                                   If[EnergyGainQ[model, vars, svars, #][[1]], 
                                                      GetBoundingInvariants[#,vars,svars,cut,TRules,lvar,mvar,"ExplosionQ"->False], ## &[]] & /@ polyseeds, ## &[]], {m, mimlist}];
  submodel = DeleteCases[{Transpose[MapIndexed[{ToExpression["OptBoundPEM" <> ToString[First@#2]], #1} &, ExInvariants]]}, {}];
  Return[submodel]
]

GetOptModel[model_, vars_, svars_, tsused_, cut_, OpMat_, TRules_, lvar_, mvar_] := Module[{optts, m, mimlist, mimjoint, PolyCheckList, polyseeds, ExInvariants, submodel, vcharacters, ActivePolyElementary, ActiveCharacter, tmpvars, ExInvCandidate, exorder, invariants},
  vcharacters = Flatten[InvariantCharacter[#, Flatten@vars] & /@ (Model2List[model]\[Transpose][[2]]), 1];
  mimlist = Sign@Last[#][[5]] &/@ ExtractMinimumTS[tsused, TRules, OpMat];
  mimjoint = Select[mimlist, ! MemberQ[vcharacters, Abs@#]&];

  ExInvariants = Table[tmpvars = DeleteCases[m*Flatten[vars], 0];
                       exorder = 4;
                       invariants = {};
                       While[invariants === {} && exorder <= cut,
                             invariants = Flatten@Select[Flatten[GetInvariants[TRules,
                                                                               tmpvars /. Var2Var[lvar, mvar, 1],
                                                                               {exorder},
                                                                               "round" -> 10^-6,
                                                                               "eventerms" -> True,
                                                                               "check" -> False] /. Var2Var[lvar, mvar, 2]],
                                                         !IntersectingQ[vcharacters, InvariantCharacter[#, Flatten@vars]]&];
                             exorder = exorder + 2];
                       invariants, {m, mimjoint}];

  submodel = DeleteCases[{Transpose[MapIndexed[{ToExpression["OptBoundPEM" <> ToString[First@#2]], #1} &, Flatten@ExInvariants]]}, {}];
  Return[submodel]
]

GetFixMimModel[model_, vars_, order_, tsused_, OpMat_, TRules_, lvar_, mvar_, OptionsPattern[{"all"->True}]] := Module[{mimlist, pool, m, ExInvariants, submodel, vcharacters, ch},
  vcharacters = Flatten[InvariantCharacter[#, Flatten@vars] & /@ (Model2List[model]\[Transpose][[2]]), 1];
  mimlist = Sign@Abs@Last[#][[5]] &/@ ExtractMinimumTS[tsused, TRules, OpMat];

  pool = Flatten@GetInvariants[TRules,
                               vars /. Var2Var[lvar, mvar, 1],
                               {order},
                               "round" -> 10^-6,
                               "eventerms" -> True,
                               "check" -> False] /. Var2Var[lvar, mvar, 2];

  ExInvariants = Intersection@@Table[Select[pool, And@@Table[Norm[(1-m)*ch]!=0, {ch, InvariantCharacter[#, Flatten@vars]}]&], {m, mimlist}];
  If[!OptionValue["all"], ExInvariants = Select[ExInvariants, !IntersectingQ[vcharacters, InvariantCharacter[#, Flatten@vars]]&]];

  submodel = DeleteCases[{Transpose[MapIndexed[{ToExpression["OptBoundSaddle" <> ToString[First@#2]], #1} &, Flatten@ExInvariants]]}, {}];
  Return[submodel]
]

ExtractMinimumTS[tsdata_, TRules_, OpMat_, OptionsPattern[{"type"->"Strained"}]] := Module[{m, ts, n, ssub, StrainMat, SymMat, tsfull, tscubic, i, j, xx},
  ssub = Cases[#, (Subscript[ToExpression["Epsilon"], i_, j_] -> xx_) :> xx] & /@ TRules;
  StrainMat = Table[Coefficient[#, ssub[[1]]] &/@ m, {m, ssub}];
  SymMat = Round@Transpose[{OpMat, StrainMat}];
  
  tsfull= Select[DeleteDuplicates[tsdata, 
                 Complement[VarSubStar[Last[#1][[5]], Last[#1][[6]], SymMat],
                            VarSubStar[Last[#2][[5]], Last[#2][[6]], SymMat]] === {} &], 
                 Norm[T5Data[Last@#]] != 0. && Norm[T6Data[Last@#]] != 0. &];
  tscubic = Select[DeleteDuplicates[tsdata, 
                   Complement[VarSubStar[Last[#1][[5]], Last[#1][[6]], SymMat],
                              VarSubStar[Last[#2][[5]], Last[#2][[6]], SymMat]] === {} &], 
                   Norm[T6Data[Last@#]] == 0. &];

  ts = Which[OptionValue["type"]==="StrainFree", tscubic,
             OptionValue["type"]==="Strained", tsfull,
             OptionValue["type"]==="all", Join[tscubic, tsfull]];
  Return[ts]
]

Model2Association[model_] := Module[{list, dict},
  list = Model2List[model];
  dict = Association[#2 -> #1 &@@@ list];
  Return[dict]
]

ModelOptimization[model_, adjustlist_, ts_, vars_, svars_, OpMat_, mim_, OptionsPattern[{"round" -> 10^-6, "weight"->1}]] := Module[{m, t, optmodel, Ham, tssub, fullvars, GF, sol, character, tsdata, out, init, mdict},
  fullvars = Append[vars, svars];
  character = Flatten[Table[Round[m].#[[1 ;; Length[Flatten@vars]]] & /@ mim, {m, OpMat}], 1];
  tsdata = Select[ts, Norm[Last[#][[6]]] == 0. && MemberQ[character, Sign@Last[#][[5]]] &];
  Print["training set used:"];
  Print[Last[#] & /@ tsdata];
  init = MapIndexed[Which[Length[#1] == 1, {ToExpression["Adjust" <> ToString[First@#2]], 0, -Infinity, Infinity, Last@#1},
                          Length[#1] == 2, {ToExpression["Adjust" <> ToString[First@#2]], #1[[1]], #1[[1]], #1[[1]], Last@#1},
                          Length[#1] == 4, {ToExpression["Adjust" <> ToString[First@#2]], #1[[1]], #1[[2]], #1[[3]], Last@#1}] &, adjustlist];
  optmodel = {Transpose[{#1, #5} & @@@ init]};
  Ham = GetHam[JoinModels[model, optmodel]];
  GF = Sum[tssub = Thread[Flatten[fullvars] -> Flatten[Last[t][[5 ;; 6]]]];
           (OptionValue["weight"]*(Last[t][[4]] - Ham)^2 + Sum[Expand[D[Ham, {var, 1}]]^2, {var, Flatten@fullvars}]) /. tssub, {t, tsdata}];
  sol = FindMinimum[GF, {#1, #2, #3, #4} & @@@ init];
  Print["Residule: " <> ToString[sol[[1]], StandardForm]];
  out = Chop[optmodel /. sol[[2]], OptionValue["round"]];
  Return[out]
]

QuickFit[poly_, tsdata_, vars_, svars_, prefix_] := Module[{data, model, OffsetFunc, fullvars},
  fullvars = Append[vars, svars];
  OffsetFunc = (prefix /. Thread[Flatten@fullvars -> {##}] &);
  data = Flatten[{#5, #6, #4}] & @@@ Flatten[GetPolynomialTrainingSet[#, tsdata] & /@ poly, 2];
  model = LinearModelFit[data, poly, Flatten@fullvars, IncludeConstantBasis -> False, LinearOffsetFunction -> OffsetFunc]["BestFit"];
  Return[model]
]

FitHessian[Basis_, LTB_, HJij_, pos_, vars_, OptionsPattern[{"unit"->"atomic"}]] := Module[{TB, HR, latt, sites, EPhonon, FCC, FCCC, shift, t, d, v, i, j, h, vl, vr, var2, out1, out2, out, alat, ev2hartree, na, nb, cbasis, param, terms},
  {latt, sites} = pos;
  {na, nb} = Dimensions[Basis];
  cbasis = Table[Flatten[latt\[Transpose].# & /@ Partition[Basis\[Transpose][[i]], 3]], {i, nb}]\[Transpose];
  If[OptionValue["unit"]==="evangstrom", alat = 1.0; ev2hartree=1.0, alat=Mean[Norm[#] & /@ latt]; ev2hartree=1/27.2114];

  TB = GroupBy[LTB, #[[2]] &];

  EPhonon = Expand@Sum[FCC = ArrayFlatten@Table[HR = GroupBy[TB[shift], First][{i, j}];
                                                Total[#4/#3 & @@@ HR], {i, Length[sites]}, {j, Length[sites]}];
                       FCCC = Chop[cbasis\[Transpose].FCC.cbasis];
                       vl = (Flatten[vars] /. {Subscript[v_, d_] :> Subscript[v, d, 0, 0, 0]});
                       vr = (Flatten[vars] /. {Subscript[v_, d_] :> Subscript[v, d, shift[[1]], shift[[2]], shift[[3]]]});
                       If[shift == {0, 0, 0}, 0.5, 1] ev2hartree alat^2 vl.FCCC.vr, {shift, Keys@TB}];

  out1 = If[HJij==={},
            CollectPolynomial[EPhonon, {}, "round" -> 1.0*10^-6], 
            Table[terms = If[Head[h] === Plus, Level[h, {1}], {h}];
                  param = Table[v = Variables[t];
                                var2 = If[Length@v == 1, ConstantArray[v[[1]], 2], v];
                                Length[v]/2 D[EPhonon, {#1, 1}, {#2, 1}] & @@ var2, {t, terms}];
                                {Mean[param], h}, {h, HJij}]];

  out2 = CollectPolynomial[Chop[EPhonon - Expand[Dot@@Transpose[out1]]], {}, "round" -> 1.0*10^-6];

  out = out1;

  (*out = CollectPolynomial[EPhonon];*)

  Return[{out, out2}]
]

FitSupercellHessian[Basis_, Fij_, qgrid_, TRules4Jij_, pos0_, vars_, OptionsPattern[{"sort"->True, "rotate"->True, "SymmetricBasisQ" -> True, "unit"->"atomic"}]]:=Module[{TB, HR, latt, sites, EPhonon, ECCC, FCC, FCCC, shift, t, d, v, i, j, h, vl, vr, var2, out1, out2, out, alat, ev2hartree, na, nb, cbasis, param, terms, supervars, FC, SFC},
  {latt, sites} = pos0;
  {na, nb} = Dimensions[Basis];
  {cbasis, supervars} = BasisInSupercell[pos0, Basis, vars, qgrid, OptionValue["SymmetricBasisQ"], "rotate" -> OptionValue["rotate"]];
  cbasis = Transpose[Flatten@Table[Chop[latt.(Partition[#, 3][[i]])], {i, Length@Fij}] &/@ Transpose[cbasis]];
 If[OptionValue["unit"]==="evangstrom", alat = 1.0; ev2hartree=1.0, alat=Mean[Norm[#] & /@ latt]; ev2hartree=1/27.2114];

  (*SFC = 0.5*cbasis\[Transpose].(ArrayFlatten[Fij]+ArrayFlatten[Fij]\[Transpose]).cbasis;*)
  SFC = cbasis\[Transpose].ArrayFlatten[Fij].cbasis;
  SFC = Chop[0.5 ev2hartree alat^2 (SFC\[Transpose] + SFC)];
  
  EPhonon = Expand[supervars[[1;;nb]].SFC[[1;;nb,;;]].supervars];
  EPhonon = Total[If[OnSiteExprQ[#], 0.5, 1] # & /@ Level[EPhonon, {1}]];

  out = CollectPolynomial[EPhonon, TRules4Jij, "round" -> 1.0*10^-9, "sort"->OptionValue["sort"]];

  Return[out]
]

GetAcousticEnergy[spg0_, tgrp_, nn_, LTB_, pos0_, OptionsPattern[{"unit"->"atomic", "SymmetricBasisQ" -> False}]] := Module[{AcousticEnergy, AcousticOpMat, AcousticTRules, latt, sites, TB, ECCC, EPhonon, FCC, FCCC, shift, HR, i, j, v, d, t, h, vl, var2, vr, ix, iy, iz, isovars, varsub, terms, param, out, out1, out2, alat, ev2hartree, cbasis, avars},
  {latt, sites} = pos0;
  avars = Subscript[ToExpression["u"], #] & /@ Range[3];
  cbasis = (Normalize[#] & /@ (Flatten[Table[1.0 IdentityMatrix[3], {Length[sites]}], 1]\[Transpose]))\[Transpose];
  If[OptionValue["unit"]==="evangstrom", alat = 1.0; ev2hartree=1.0, alat=Mean[Norm[#] & /@ latt]; ev2hartree=1/27.2114];

  TB = GroupBy[LTB, #[[2]] &];
  ECCC = Expand@Sum[FCC = ArrayFlatten@Table[HR = GroupBy[TB[shift], First][{i, j}];
                                             Total[#4/#3 & @@@ HR], {i, Length[sites]}, {j, Length[sites]}];
                    FCCC = Chop[cbasis\[Transpose].FCC.cbasis];
                    vl = (avars /. {Subscript[v_, d_] :> Subscript[v, d, 0, 0, 0]});
                    vr = (avars /. {Subscript[v_, d_] :> Subscript[v, d, shift[[1]], shift[[2]], shift[[3]]]});
                    If[shift == {0, 0, 0}, 0.5, 0.5] vl.FCCC.vr, {shift, Keys@TB}];

  EPhonon = Expand@Sum[FCC = ArrayFlatten@Table[HR = GroupBy[TB[shift], First][{i, j}];
                                                Total[#4/#3 & @@@ HR], {i, Length[sites]}, {j, Length[sites]}];
                       FCCC = Chop[cbasis\[Transpose] . FCC . cbasis];
                       vl = (avars /. {Subscript[v_, d_] :> Subscript[v, d, 0, 0, 0]});
                       vr = (avars /. {Subscript[v_, d_] :> Subscript[v, d, shift[[1]], shift[[2]], shift[[3]]]});
                       If[shift == {0, 0, 0}, 0.5, 0.5] ev2hartree*alat^2*vl.FCCC.vr, {shift, Keys@TB}];
  
  isovars = ToExpression[#] & /@ {"\!\(\*SubscriptBox[\(Iso\), \(1\)]\)", 
                                  "\!\(\*SubscriptBox[\(Iso\), \(2\)]\)", 
                                  "\!\(\*SubscriptBox[\(Iso\), \(3\)]\)"};
  varsub = Table[Subscript[isovars[[1, 1]], i, ix_, iy_, iz_] -> Subscript[avars[[1, 1]], i, ix, iy, iz], {i, 3}];
  
  AcousticOpMat = Chop[GetMatrixRep[spg0, tgrp, pos0, cbasis, Range[3], "disp"], 10^-4];
  AcousticTRules = Join[#1, #2] & @@@ ({GetIsoTransformRules[AcousticOpMat, "disp"], GetIsoStrainTransformRules[pos0[[1]], spg0]}\[Transpose]);
  AcousticEnergy = Jij[spg0, nn, AcousticTRules, {isovars, isovars}, "OnSite" -> True] /. varsub;
 
  out1 = If[OptionValue["SymmetricBasisQ"],
            Table[terms = If[Head[h] === Plus, Level[h, {1}], {h}];
                  param = Table[v = Variables[t];
                                var2 = If[Length@v == 1, ConstantArray[v[[1]], 2], v];
                                Length[v]/2 ev2hartree*alat^2*D[ECCC, {#1, 1}, {#2, 1}] & @@ var2, {t, terms}];
                                {Mean[param], h}, {h, AcousticEnergy}]
            CollectPolynomial[EPhonon, {}, "round" -> 1.0*10^-6]];

  out2 = CollectPolynomial[Chop[EPhonon - Expand[Dot @@ Transpose[out1]]], {}];

  Return[{out1, out2}]
]

GetAcousticEps[spg0_, tgrp_, nn_, pos0_, avars_, OptionsPattern[{"unit"->"atomic"}]] := Module[{AcousticEnergy, AcousticOpMat, AcousticTRules, AcousticBasis, latt, sites, TB, ECCC, FCC, FCCC, shift, HR, i, j, v, d, h, vl, var2, vr, ix, iy, iz, isovars, varsub, param, out, alat, ev2hartree, seeds, svars},
  {latt, sites} = pos0;
  If[OptionValue["unit"]==="evangstrom", alat = 1.0; ev2hartree=1.0, alat=Mean[Norm[#] & /@ latt]; ev2hartree=1/27.2114];

  AcousticBasis = (Normalize[#] & /@ (Flatten[Table[1.0 IdentityMatrix[3], {Length[sites]}], 1]\[Transpose]))\[Transpose];

  isovars = ToExpression[#] & /@ {"\!\(\*SubscriptBox[\(Iso\), \(1\)]\)",
                                  "\!\(\*SubscriptBox[\(Iso\), \(2\)]\)",
                                  "\!\(\*SubscriptBox[\(Iso\), \(3\)]\)"};

  svars = {Subscript[ToExpression["Epsilon"], 1, 1], 
           Subscript[ToExpression["Epsilon"], 2, 2],
           Subscript[ToExpression["Epsilon"], 3, 3],
           Subscript[ToExpression["Epsilon"], 2, 3],
           Subscript[ToExpression["Epsilon"], 1, 3],
           Subscript[ToExpression["Epsilon"], 1, 2]};

  varsub = Join[Table[Subscript[isovars[[1, 1]], i, ix_, iy_, iz_] -> Subscript[avars[[1, 1]], i, ix, iy, iz], {i, 3}],
                (Subscript[ToExpression["Epsilon"], #2, #3, ix_, iy_, iz_] -> Subscript[ToExpression["\[Epsilon]"], #2, #3, ix, iy, iz]) &@@@ svars];

  AcousticOpMat = Chop[GetMatrixRep[spg0, tgrp, pos0, AcousticBasis, Range[3], "disp"], 10^-4];
  AcousticTRules = Join[#1, #2] & @@@ ({GetIsoTransformRules[AcousticOpMat, "disp"], GetIsoStrainTransformRules[pos0[[1]], spg0]}\[Transpose]);
  seeds = {MonomialList[Expand[Total[isovars] Total[svars]]], isovars};
  AcousticEnergy = Jij[spg0, nn, AcousticTRules, seeds, "OnSite" -> True] /. varsub;

  Return[AcousticEnergy]
]

BuildHeterostructureStrain[dir_, vars_, svars_, nn_, ie_] := Module[{i, j, ix, iy, iz, v, d, eps, StrainFromuFortran, dfields, xyzsub, f90, body, eu, eusub, SymmetrizeStrain},

  SymmetrizeStrain = {{ToString@FortranForm[ToExpression["euij"][2, 1, ToExpression["x0"], ToExpression["y0"], ToExpression["z0"]]], 
                       {ToExpression["euij"][1, 2, ToExpression["x0"], ToExpression["y0"], ToExpression["z0"]]}},
                      {ToString@FortranForm[ToExpression["euij"][3, 1, ToExpression["x0"], ToExpression["y0"], ToExpression["z0"]]], 
                       {ToExpression["euij"][1, 3, ToExpression["x0"], ToExpression["y0"], ToExpression["z0"]]}},
                      {ToString@FortranForm[ToExpression["euij"][3, 2, ToExpression["x0"], ToExpression["y0"], ToExpression["z0"]]], 
                       {ToExpression["euij"][2, 3, ToExpression["x0"], ToExpression["y0"], ToExpression["z0"]]}}};

  xyzsub = {ToExpression["x"] -> 1, ToExpression["y"] -> 2, ToExpression["z"] -> 3};

  dfields = Join[vars /. {Subscript[v_, d_] -> Subscript[v, d, 0, 0, 0]} /. xyzsub, 
                 {ToExpression[#] & /@ {"\!\(\*SubscriptBox[\(u\), \(1, 0, 0, 0\)]\)", 
                                        "\!\(\*SubscriptBox[\(u\), \(2, 0, 0, 0\)]\)", 
                                        "\!\(\*SubscriptBox[\(u\), \(3, 0, 0, 0\)]\)"}}, 
                 {ToExpression[#] & /@ {"\!\(\*SubscriptBox[\(\[Epsilon]0\), \(1, 1\)]\\)", 
                                        "\!\(\*SubscriptBox[\(\[Epsilon]0\), \(2, 2\)]\)", 
                                        "\!\(\*SubscriptBox[\(\[Epsilon]0\), \(3, 3\)]\)", 
                                        "\!\(\*SubscriptBox[\(\[Epsilon]0\), \(2, 3\)]\)", 
                                        "\!\(\*SubscriptBox[\(\[Epsilon]0\), \(1, 3\)]\)", 
                                        "\!\(\*SubscriptBox[\(\[Epsilon]0\), \(1, 2\)]\)"}}];
  
  eusub = StrainFromu[eu, ToExpression["u"], ie];
  StrainFromuFortran = {Keys[eusub] /. {Subscript[eu, i_, j_, x0_, y0_, z0_] :> ToString@FortranForm[ToExpression["euij"][i, j, ToExpression["x0"], ToExpression["y0"], ToExpression["z0"]]]}, 
                        Table[TimesFactor2Real[#] & /@ If[MatchQ[eps, Plus[_, __]], Level[eps, {1}], Level[eps, {0}]] /. FortranVarSub[dfields], {eps, Values[eusub]}]}\[Transpose];
  
  f90 = "GetHeterostructureStrain";
  body = Flatten[FortranExprBlock[#1, #2, nn, 1] & @@@ Join[StrainFromuFortran, SymmetrizeStrain], 1];
  MMA2FORTRAN[dir <> "/src/hamiltonian/" <> f90, StrainFromuF90[body, nn]];
]

ExportHamiltonianOld[dir_, fullmodel_, vars_, svars_, ie_, EwaldField_:True] := Module[{gmodelextend, hmodel, gmodel, x0, y0, z0, ix, iy, iz, dx, dy, dz, g1, g2, g3, g3u, h1, h11, h12, h2, h3, h3u, d, v, var0sub, varxsub, xyzsub, dfields, totstrainsub, eu, eusub, e0mute, varsmute, site0sub, uxsub, coefficients, eps0sub, nn, ExprFunc},

  totstrainsub = {Subscript[ToExpression["\[Epsilon]"], i_, j_] :> Subscript[ToExpression["\[Epsilon]0"], i, j] + Subscript[ToExpression["\[Epsilon]u"], i, j, 0, 0, 0], Subscript[ToExpression["\[Epsilon]"], i_, j_, x0_, y0_, z0_] :> Subscript[ToExpression["\[Epsilon]0"], i, j] + Subscript[ToExpression["\[Epsilon]u"], i, j, x0, y0, z0]};
  eps0sub = {Subscript[ToExpression["\[Epsilon]"], i_, j_] :> Subscript[ToExpression["\[Epsilon]0"], i, j], Subscript[ToExpression["\[Epsilon]"], i_, j_, x0_, y0_, z0_] :> Subscript[ToExpression["\[Epsilon]0"], i, j]};
  eusub = Normal@StrainFromu[ToExpression["\[Epsilon]u"], ToExpression["u"], 1];
  xyzsub = {ToExpression["x"] -> 1, ToExpression["y"] -> 2, ToExpression["z"] -> 3};
  e0mute = {Subscript[ToExpression["\[Epsilon]0"], i_, j_] :> Subscript[ToExpression["\[Epsilon]0"], i, j, 0]};
  varsmute = Subscript[#1, #2 /. xyzsub, ix_, iy_, iz_] :> Subscript[#1, #2 /. xyzsub, ix, iy, iz, 0] & @@@ Flatten[vars];
  site0sub = {Subscript[v_, d_] -> Subscript[v, d, 0, 0, 0], Subscript[v_, i_, j_] -> Subscript[v, i, j, 0, 0, 0]};
  var0sub = Thread[(Flatten[vars] /. xyzsub) -> Flatten[vars /. xyzsub /. site0sub]];
  varxsub = Thread[(Flatten[vars] /. xyzsub) -> Flatten[vars /. xyzsub /. {Subscript[v_, d_] -> Subscript[v, d, dx, dy, dz, 0]}]];
  uxsub = {Subscript[ToExpression["u"], d_, x_, y_, z_] :> Subscript[ToExpression["u"], d, x + dx, y + dy, z + dz]};
  
  dfields = Join[vars /. site0sub /. xyzsub, {ToExpression[#] & /@ {"\!\(\*SubscriptBox[\(u\), \(1, 0, 0, 0\)]\)", 
                                                                    "\!\(\*SubscriptBox[\(u\), \(2, 0, 0, 0\)]\)", 
                                                                    "\!\(\*SubscriptBox[\(u\), \(3, 0, 0, 0\)]\)"}},
                                              {ToExpression[#] & /@ {"\!\(\*SubscriptBox[\(\[Epsilon]0\), \(1, 1\)]\)", 
                                                                     "\!\(\*SubscriptBox[\(\[Epsilon]0\), \(2, 2\)]\)", 
                                                                     "\!\(\*SubscriptBox[\(\[Epsilon]0\), \(3, 3\)]\)", 
                                                                     "\!\(\*SubscriptBox[\(\[Epsilon]0\), \(2, 3\)]\)", 
                                                                     "\!\(\*SubscriptBox[\(\[Epsilon]0\), \(1, 3\)]\)", 
                                                                     "\!\(\*SubscriptBox[\(\[Epsilon]0\), \(1, 2\)]\)"}}];
  
  hmodel = Select[Model2List[fullmodel /. xyzsub], JijInvQ[#[[2]]] &];
  gmodel = Select[Model2List[fullmodel /. xyzsub], ! JijInvQ[#[[2]]] &];
  
  g1 = Expand[Select[gmodel, And@@StrainInvQ[Variables[#[[2]]], ToExpression["\[Epsilon]"]] &]];
  g2 = Expand[Select[gmodel, Nor@@StrainInvQ[Variables[#[[2]]], ToExpression["\[Epsilon]"]] &]];
  g3 = Expand[Complement[gmodel, Join[g1, g2]] /. var0sub /. totstrainsub];

  (*g1 = Expand[g1 /. totstrainsub /. eps0sub];*)
  g1 = Expand[g1 /. eps0sub];
  g2 = g2 /. var0sub;
  g3u = Expand[g3 /. eusub /. e0mute /. varsmute];

  h11 = Select[hmodel, And@@XInvQ[Variables[#[[2]]], {ToExpression["u"]}] &];
  h12 = Select[hmodel, (Or@@XInvQ[Variables[#[[2]]], {ToExpression["u"]}]) && (Or@@XInvQ[Variables[#[[2]]], {ToExpression["\[Epsilon]"]}]) &];
  h2 = Select[hmodel, Nor@@XInvQ[Variables[#[[2]]], {ToExpression["\[Epsilon]"], ToExpression["\[Epsilon]0"], ToExpression["u"]}] &];
  h3 = Expand[Complement[hmodel, Join[h11, h12, h2]] /. totstrainsub];
  h12 = {#1, SimplifyCommonFactor[#2, 10^-6]} & @@@ (h12 /. eps0sub);
  h1 = Join[h11, h12];
  h3u = Expand[h3 /. eusub /. e0mute /. varsmute];

  nn = Max[Max[#] & /@ (Abs@Cases[Variables@h3, Subscript[ToExpression["\[Epsilon]u"], __, __, ix_, iy_, iz_] -> {ix, iy, iz}]\[Transpose])];        
  BuildHeterostructureStrain[dir, vars, svars, nn, ie];

  coefficients= {{"CoeffGstrain",     g1\[Transpose][[1]]}, 
                 {"CoeffGprimary",    g2\[Transpose][[1]]}, 
                 {"CoeffGpscoupling", g3\[Transpose][[1]]}, 
                 {"CoeffHstrain",     h1\[Transpose][[1]]}, 
                 {"CoeffHprimary",    h2\[Transpose][[1]]}, 
                 {"CoeffHpscoupling", h3\[Transpose][[1]]}};

  WriteLatticeModelF90[dir<>"/src/hamiltonian", "GetVariationGstrain",   "Variation", {{"Gstrain",      g1\[Transpose][[2]]}}, dfields, "LineLimit"->OptionValue["LineLimit"]];
  WriteLatticeModelF90[dir<>"/src/hamiltonian", "GetVariationGprimary",  "Variation", {{"Gprimary",     g2\[Transpose][[2]]}}, dfields];
  WriteLatticeModelF90[dir<>"/src/hamiltonian", "GetVariationGpscoupling","Variation", {{"Gpscoupling", g3\[Transpose][[2]] 
                                                                                                       +g3u\[Transpose][[2]]}}, dfields];

  WriteLatticeModelF90[dir<>"/src/hamiltonian", "GetVariationG",         "Variation", {{"Gstrain",     g1\[Transpose][[2]]},
                                                                                       {"Gprimary",    g2\[Transpose][[2]]},
                                                                                       {"Gpscoupling", g3\[Transpose][[2]]
                                                                                                      +g3u\[Transpose][[2]]}}, dfields];

  WriteLatticeModelF90[dir<>"/src/hamiltonian", "GetVariationHstrain",   "Variation", {{"Hstrain",      h1\[Transpose][[2]]}}, dfields];
  WriteLatticeModelF90[dir<>"/src/hamiltonian", "GetVariationHprimary",  "Variation", {{"Hprimary",     h2\[Transpose][[2]]}}, dfields];
  WriteLatticeModelF90[dir<>"/src/hamiltonian", "GetVariationHpscoupling","Variation", {{"Hpscoupling", h3\[Transpose][[2]] 
                                                                                                       +h3u\[Transpose][[2]]}}, dfields];

  WriteLatticeModelF90[dir<>"/src/hamiltonian", "GetVariationH",         "Variation", {{"Hstrain",     h1\[Transpose][[2]]}, 
                                                                                       {"Hprimary",    h2\[Transpose][[2]]}, 
                                                                                       {"Hpscoupling", h3\[Transpose][[2]] 
                                                                                                      +h3u\[Transpose][[2]]}}, dfields];

  ExprFunc = WriteLatticeModelF90[dir<>"/src/hamiltonian", "GetVariation",          "Variation", {{"Gstrain",     g1\[Transpose][[2]]},
                                                                                                  {"Gprimary",    g2\[Transpose][[2]]},
                                                                                                  {"Gpscoupling", g3\[Transpose][[2]]
                                                                                                                 +g3u\[Transpose][[2]]},
                                                                                                  {"Hstrain",     h1\[Transpose][[2]]}, 
                                                                                                  {"Hprimary",    h2\[Transpose][[2]]}, 
                                                                                                  {"Hpscoupling", h3\[Transpose][[2]] 
                                                                                                                 +h3u\[Transpose][[2]]}}, dfields];

  MMA2FORTRAN[dir<>"/src/hamiltonian/" <> "VariationExpr", ExprFunc];
  
  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetForcesGstrain",    "Forces", {{"Gstrain",     g1\[Transpose][[2]]}}, dfields, "AllSites" -> True];
  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetForcesGprimary",   "Forces", {{"Gprimary",    g2\[Transpose][[2]]}}, dfields, "AllSites" -> True];
  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetForcesGpscoupling","Forces", {{"Gpscoupling", g3\[Transpose][[2]]
                                                                                                +g3u\[Transpose][[2]]}}, dfields, "AllSites" -> True];

  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetForcesG",         "Forces", {{"Gstrain",      g1\[Transpose][[2]]},
                                                                                {"Gprimary",     g2\[Transpose][[2]]},
                                                                                {"Gpscoupling",  g3\[Transpose][[2]]
                                                                                                +g3u\[Transpose][[2]]}}, dfields, "AllSites" -> True];

  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetForcesHstrain",    "Forces", {{"Hstrain",     h1\[Transpose][[2]]}}, dfields, "AllSites" -> True];
  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetForcesHprimary",   "Forces", {{"Hprimary",    h2\[Transpose][[2]]}}, dfields, "AllSites" -> True];
  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetForcesHpscoupling","Forces", {{"Hpscoupling", h3\[Transpose][[2]]
                                                                                                +h3u\[Transpose][[2]]}}, dfields, "AllSites" -> True];

  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetForcesH",         "Forces", {{"Hstrain",      h1\[Transpose][[2]]}, 
                                                                                {"Hprimary",     h2\[Transpose][[2]]}, 
                                                                                {"Hpscoupling",  h3\[Transpose][[2]]
                                                                                                +h3u\[Transpose][[2]]}}, dfields, "AllSites" -> True];

  ExprFunc = WriteLatticeModelF90[dir <> "/src/hamiltonian", "GetForces",      "Forces", {{"Gstrain",     g1\[Transpose][[2]]},
                                                                                          {"Gprimary",    g2\[Transpose][[2]]},
                                                                                          {"Gpscoupling", g3\[Transpose][[2]]
                                                                                                         +g3u\[Transpose][[2]]},
                                                                                          {"Hstrain",     h1\[Transpose][[2]]},
                                                                                          {"Hprimary",    h2\[Transpose][[2]]},
                                                                                          {"Hpscoupling", h3\[Transpose][[2]]
                                                                                                         +h3u\[Transpose][[2]]}}, dfields, "AllSites" -> True];
  MMA2FORTRAN[dir<>"/src/hamiltonian/" <> "ForcesExpr", ExprFunc];

  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetSiteEnergyGstrain",   "SiteEnergy",{{"Gstrain",      g1\[Transpose][[2]]}}, dfields];
  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetSiteEnergyGprimary",  "SiteEnergy",{{"Gprimary",     g2\[Transpose][[2]]}}, dfields];
  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetSiteEnergyGpscoupling","SiteEnergy",{{"Gpscoupling", g3\[Transpose][[2]]}}, dfields];

  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetSiteEnergyG",         "SiteEnergy",{{"Gstrain",      g1\[Transpose][[2]]},
                                                                                       {"Gprimary",     g2\[Transpose][[2]]},
                                                                                       {"Gpscoupling",  g3\[Transpose][[2]]}}, dfields];

  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetSiteEnergyHstrain",    "SiteEnergy",{{"Hstrain",     FixDoubleCounting[h1\[Transpose][[2]]]}}, dfields];
  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetSiteEnergyHprimary",   "SiteEnergy",{{"Hprimary",    FixDoubleCounting[h2\[Transpose][[2]]]}}, dfields];
  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetSiteEnergyHpscoupling","SiteEnergy",{{"Hpscoupling", FixDoubleCounting[h3\[Transpose][[2]]]}}, dfields];

  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetSiteEnergyH",         "SiteEnergy",{{"Hstrain",    FixDoubleCounting[h1\[Transpose][[2]]]},
                                                                                       {"Hprimary",   FixDoubleCounting[h2\[Transpose][[2]]]},
                                                                                       {"Hpscoupling",FixDoubleCounting[h3\[Transpose][[2]]]}}, dfields];

  ExprFunc = WriteLatticeModelF90[dir<>"/src/hamiltonian","GetSiteEnergy",          "SiteEnergy",{{"Gstrain",     g1\[Transpose][[2]]},
                                                                                                  {"Gprimary",    g2\[Transpose][[2]]},
                                                                                                  {"Gpscoupling", g3\[Transpose][[2]]},
                                                                                                  {"Hstrain",     FixDoubleCounting[h1\[Transpose][[2]]]},
                                                                                                  {"Hprimary",    FixDoubleCounting[h2\[Transpose][[2]]]},
                                                                                                  {"Hpscoupling", FixDoubleCounting[h3\[Transpose][[2]]]}}, dfields];
  MMA2FORTRAN[dir<>"/src/hamiltonian/" <> "SiteEnergyExpr", ExprFunc];
  
  Return[{{g1, g2, g3, h1, h2, h3}, coefficients}]
]

CookModel[fullmodel_, vars_] := Module[{gmodel, hmodel, v, d, i, j, xyzsub, site0sub, var0sub, eps0sub, Gs, Gp, Gsp, Hu, Hp, Hup, Hsu, Hsp, Hsup},

   xyzsub = {ToExpression["x"] -> 1, ToExpression["y"] -> 2, ToExpression["z"] -> 3};
   site0sub = {Subscript[v_, d_] -> Subscript[v, d, 0, 0, 0], Subscript[v_, i_, j_] -> Subscript[v, i, j, 0, 0, 0]};
   var0sub = Thread[(Flatten[vars] /. xyzsub) -> Flatten[vars /. xyzsub /. site0sub]];   
   eps0sub = {Subscript[ToExpression["\[Epsilon]"], i_, j_] :> Subscript[ToExpression["\[Epsilon]0"], i, j]};

   hmodel = Select[Model2List[fullmodel /. xyzsub], JijInvQ[#[[2]]] &];
   gmodel = Select[Model2List[fullmodel /. xyzsub], ! JijInvQ[#[[2]]] &];
 
   Gs   = Expand[Select[gmodel, And@@XInvQ[#[[2]], {ToExpression["\[Epsilon]"]}, ContainsExactly] &]];
   Gp   = Expand[Select[gmodel, And@@XInvQ[#[[2]], {ToExpression["\[Epsilon]"]}, ContainsNone] &]];
   Gsp  = Expand[Complement[gmodel, Join[Gs, Gp]]];
 
   Gs   = Gs  /. eps0sub;
   Gp   = Gp  /. var0sub;
   Gsp  = Gsp /. var0sub /. eps0sub;
 
   Hu   = Select[hmodel, And@@XInvQ[#[[2]], {ToExpression["u"]}, ContainsExactly] &];
   Hp   = Select[hmodel, And@@XInvQ[#[[2]], {ToExpression["\[Epsilon]"], ToExpression["u"]}, ContainsNone] &];
 
   Hsu  = Select[hmodel, And@@XInvQ[#[[2]], {ToExpression["u"], ToExpression["\[Epsilon]"]}, ContainsExactly] &];
   Hup  = Select[Expand[Complement[hmodel, Join[Hu, Hp, Hsu]]], And@@XInvQ[#[[2]], {ToExpression["\[Epsilon]"]}, ContainsNone] &];
   Hsp  = Select[Expand[Complement[hmodel, Join[Hu, Hp, Hsu, Hup]]], And@@XInvQ[#[[2]], {ToExpression["u"]}, ContainsNone] &];
   Hsup = Expand[Complement[hmodel,        Join[Hu, Hp, Hsu, Hup, Hsp]]];
 
   Hsu  = Hsu /. eps0sub;
   Hsp  = Hsp /. eps0sub;
   Hsup = Hsup /. eps0sub;

   Gs   = ReverseSortBy[Gs,   Abs[#[[1]]]&];
   Gp   = ReverseSortBy[Gp,   Abs[#[[1]]]&];
   Gsp  = ReverseSortBy[Gsp,  Abs[#[[1]]]&];
 
   Hu   = ReverseSortBy[Hu,   Abs[#[[1]]]&];
   Hp   = ReverseSortBy[Hp,   Abs[#[[1]]]&];
   Hup  = ReverseSortBy[Hup,  Abs[#[[1]]]&];
   Hsu  = ReverseSortBy[Hsu,  Abs[#[[1]]]&];
   Hsp  = ReverseSortBy[Hsp,  Abs[#[[1]]]&];
   Hsup = ReverseSortBy[Hsup, Abs[#[[1]]]&];
 
   If[Hsu  === {}, Hsu  = {{0, 1}}];
   If[Hsp  === {}, Hsp  = {{0, 1}}];
   If[Hsup === {}, Hsup = {{0, 1}}];

   Return[{Gs, Gp, Gsp, Hu, Hp, Hup, Hsu, Hsp, Hsup}]
]

ExportHamiltonian[dir_, fullmodel_, vars_, svars_, ie_, EwaldField_:True, OptionsPattern[{"LineLimit"->500}]] := Module[{gmodelextend, hmodel, gmodel, x0, y0, z0, ix, iy, iz, dx, dy, dz, Gs, Gp, Gsp, Hu, Hp, Hup, Hsu, Hsp, Hsup, d, v, varxsub, dfields, eu, eusub, e0mute, varsmute, site0sub, uxsub, coefficients, nn, ExprFunc},

  site0sub = {Subscript[v_, d_] -> Subscript[v, d, 0, 0, 0]};
  dfields = Join[vars /. site0sub, {ToExpression[#] & /@ {"\!\(\*SubscriptBox[\(u\), \(1, 0, 0, 0\)]\)",
                                                          "\!\(\*SubscriptBox[\(u\), \(2, 0, 0, 0\)]\)",
                                                          "\!\(\*SubscriptBox[\(u\), \(3, 0, 0, 0\)]\)"}},
                                   {ToExpression[#] & /@ {"\!\(\*SubscriptBox[\(\[Epsilon]0\), \(1, 1\)]\)",
                                                          "\!\(\*SubscriptBox[\(\[Epsilon]0\), \(2, 2\)]\)",
                                                          "\!\(\*SubscriptBox[\(\[Epsilon]0\), \(3, 3\)]\)",
                                                          "\!\(\*SubscriptBox[\(\[Epsilon]0\), \(2, 3\)]\)",
                                                          "\!\(\*SubscriptBox[\(\[Epsilon]0\), \(1, 3\)]\)",
                                                          "\!\(\*SubscriptBox[\(\[Epsilon]0\), \(1, 2\)]\)"}}];

  {Gs, Gp, Gsp, Hu, Hp, Hup, Hsu, Hsp, Hsup} = CookModel[fullmodel, vars];

  nn = Max[Max[#] & /@ (Abs@Cases[Variables@Hu, Subscript[ToExpression["\[Epsilon]u"], __, __, ix_, iy_, iz_] -> {ix, iy, iz}]\[Transpose])];
  BuildHeterostructureStrain[dir, vars, svars, nn, ie];

  coefficients= {{"CoeffGs",     Gs\[Transpose][[1]]},
                 {"CoeffGp",     Gp\[Transpose][[1]]},
                 {"CoeffGsp",   Gsp\[Transpose][[1]]},
                 {"CoeffHu",     Hu\[Transpose][[1]]},
                 {"CoeffHp",     Hp\[Transpose][[1]]},
                 {"CoeffHup",   Hup\[Transpose][[1]]},
                 {"CoeffHsu",   Hsu\[Transpose][[1]]},
                 {"CoeffHsp",   Hsp\[Transpose][[1]]},
                 {"CoeffHsup", Hsup\[Transpose][[1]]}};

  WriteLatticeModelF90[dir<>"/src/hamiltonian", "GetVariationGs",  "Variation",         {{"Gs",   Gs\[Transpose][[2]]}}, dfields, 
                       "LineLimit"->OptionValue["LineLimit"]];
  WriteLatticeModelF90[dir<>"/src/hamiltonian", "GetVariationGp",  "Variation",         {{"Gp",   Gp\[Transpose][[2]]}}, dfields, 
                       "LineLimit"->OptionValue["LineLimit"]];
  WriteLatticeModelF90[dir<>"/src/hamiltonian", "GetVariationGsp", "Variation",         {{"Gsp", Gsp\[Transpose][[2]]}}, dfields, 
                       "LineLimit"->OptionValue["LineLimit"]];

  WriteLatticeModelF90[dir<>"/src/hamiltonian", "GetVariationG",   "Variation",         {{"Gs",   Gs\[Transpose][[2]]},
                                                                                         {"Gp",   Gp\[Transpose][[2]]},
                                                                                         {"Gsp",  Gsp\[Transpose][[2]]}}, dfields, 
                       "LineLimit"->OptionValue["LineLimit"]];

  WriteLatticeModelF90[dir<>"/src/hamiltonian", "GetVariationHu",  "Variation",         {{"Hu",   Hu\[Transpose][[2]]}}, dfields, 
                       "LineLimit"->OptionValue["LineLimit"]];
  WriteLatticeModelF90[dir<>"/src/hamiltonian", "GetVariationHp",  "Variation",         {{"Hp",   Hp\[Transpose][[2]]}}, dfields, 
                       "LineLimit"->OptionValue["LineLimit"]];
  WriteLatticeModelF90[dir<>"/src/hamiltonian", "GetVariationHup", "Variation",         {{"Hup", Hup\[Transpose][[2]]}}, dfields, 
                       "LineLimit"->OptionValue["LineLimit"]];
  WriteLatticeModelF90[dir<>"/src/hamiltonian", "GetVariationHsu", "Variation",         {{"Hsu", Hsu\[Transpose][[2]]}}, dfields, 
                       "LineLimit"->OptionValue["LineLimit"]];
  WriteLatticeModelF90[dir<>"/src/hamiltonian", "GetVariationHsp", "Variation",         {{"Hsp", Hsp\[Transpose][[2]]}}, dfields, 
                       "LineLimit"->OptionValue["LineLimit"]];
  WriteLatticeModelF90[dir<>"/src/hamiltonian", "GetVariationHsup", "Variation",        {{"Hsup", Hsup\[Transpose][[2]]}}, dfields,
                       "LineLimit"->OptionValue["LineLimit"]];

  WriteLatticeModelF90[dir<>"/src/hamiltonian", "GetVariationH",   "Variation",         {{"Hu",     Hu\[Transpose][[2]]},
                                                                                         {"Hp",     Hp\[Transpose][[2]]},
                                                                                         {"Hup",   Hup\[Transpose][[2]]},
                                                                                         {"Hsu",   Hsu\[Transpose][[2]]},
                                                                                         {"Hsp",   Hsp\[Transpose][[2]]},
                                                                                         {"Hsup", Hsup\[Transpose][[2]]}}, dfields, 
                       "LineLimit"->OptionValue["LineLimit"]];

  WriteLatticeModelF90[dir<>"/src/hamiltonian", "GetVariationEps", "Variation",         {{"Gs",     Gs\[Transpose][[2]]},
                                                                                         {"Gsp",   Gsp\[Transpose][[2]]},
                                                                                         {"Hsu",   Hsu\[Transpose][[2]]},
                                                                                         {"Hsp",   Hsp\[Transpose][[2]]},
                                                                                         {"Hsup", Hsup\[Transpose][[2]]}}, dfields,
                       "LineLimit"->OptionValue["LineLimit"]];

  ExprFunc = WriteLatticeModelF90[dir<>"/src/hamiltonian", "GetVariation", "Variation", {{"Gs",     Gs\[Transpose][[2]]},
                                                                                         {"Gp",     Gp\[Transpose][[2]]},
                                                                                         {"Gsp",   Gsp\[Transpose][[2]]},
                                                                                         {"Hu",     Hu\[Transpose][[2]]},
                                                                                         {"Hp",     Hp\[Transpose][[2]]},
                                                                                         {"Hup",   Hup\[Transpose][[2]]},
                                                                                         {"Hsu",   Hsu\[Transpose][[2]]},
                                                                                         {"Hsp",   Hsp\[Transpose][[2]]}, 
                                                                                         {"Hsup", Hsup\[Transpose][[2]]}}, dfields, 
                                  "LineLimit"->OptionValue["LineLimit"]];

  MMA2FORTRAN[dir<>"/src/hamiltonian/" <> "VariationExpr", ExprFunc];

  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetForcesGs",  "Forces",            {{"Gs",   Gs\[Transpose][[2]]}}, dfields, 
                       "AllSites" -> True, "LineLimit"->OptionValue["LineLimit"]];
  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetForcesGp",  "Forces",            {{"Gp",   Gp\[Transpose][[2]]}}, dfields, 
                       "AllSites" -> True, "LineLimit"->OptionValue["LineLimit"]];
  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetForcesGsp", "Forces",            {{"Gsp", Gsp\[Transpose][[2]]}}, dfields, 
                       "AllSites" -> True, "LineLimit"->OptionValue["LineLimit"]];

  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetForcesG",   "Forces",            {{"Gs",   Gs\[Transpose][[2]]},
                                                                                     {"Gp",   Gp\[Transpose][[2]]},
                                                                                     {"Gsp", Gsp\[Transpose][[2]]}}, dfields, 
                       "AllSites" -> True, "LineLimit"->OptionValue["LineLimit"]];

  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetForcesHu",  "Forces",            {{"Hu",   Hu\[Transpose][[2]]}}, dfields, 
                       "AllSites" -> True, "LineLimit"->OptionValue["LineLimit"]];
  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetForcesHp",  "Forces",            {{"Hp",   Hp\[Transpose][[2]]}}, dfields, 
                       "AllSites" -> True, "LineLimit"->OptionValue["LineLimit"]];
  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetForcesHup", "Forces",            {{"Hup", Hup\[Transpose][[2]]}}, dfields, 
                       "AllSites" -> True, "LineLimit"->OptionValue["LineLimit"]];
  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetForcesHsu", "Forces",            {{"Hsu", Hsu\[Transpose][[2]]}}, dfields, 
                       "AllSites" -> True, "LineLimit"->OptionValue["LineLimit"]];
  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetForcesHsp", "Forces",            {{"Hsp", Hsp\[Transpose][[2]]}}, dfields, 
                       "AllSites" -> True, "LineLimit"->OptionValue["LineLimit"]];
  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetForcesHsup", "Forces",           {{"Hsup", Hsup\[Transpose][[2]]}}, dfields,
                       "AllSites" -> True, "LineLimit"->OptionValue["LineLimit"]];

  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetForcesH",   "Forces",            {{"Hu",     Hu\[Transpose][[2]]},
                                                                                     {"Hp",     Hp\[Transpose][[2]]},
                                                                                     {"Hup",   Hup\[Transpose][[2]]},
                                                                                     {"Hsu",   Hsu\[Transpose][[2]]},
                                                                                     {"Hsp",   Hsp\[Transpose][[2]]}, 
                                                                                     {"Hsup", Hsup\[Transpose][[2]]}}, dfields, 
                       "AllSites" -> True, "LineLimit"->OptionValue["LineLimit"]];

  WriteLatticeModelF90[dir<>"/src/hamiltonian", "GetForcesEps", "Forces",           {{"Gs",     Gs\[Transpose][[2]]},
                                                                                     {"Gsp",   Gsp\[Transpose][[2]]},
                                                                                     {"Hsu",   Hsu\[Transpose][[2]]},
                                                                                     {"Hsp",   Hsp\[Transpose][[2]]},
                                                                                     {"Hsup", Hsup\[Transpose][[2]]}}, dfields,
                       "AllSites" -> True, "LineLimit"->OptionValue["LineLimit"]];

  ExprFunc = WriteLatticeModelF90[dir <> "/src/hamiltonian", "GetForces", "Forces", {{"Gs",     Gs\[Transpose][[2]]},
                                                                                     {"Gp",     Gp\[Transpose][[2]]},
                                                                                     {"Gsp",   Gsp\[Transpose][[2]]},
                                                                                     {"Hu",     Hu\[Transpose][[2]]},
                                                                                     {"Hp",     Hp\[Transpose][[2]]},
                                                                                     {"Hup",   Hup\[Transpose][[2]]},
                                                                                     {"Hsu",   Hsu\[Transpose][[2]]},
                                                                                     {"Hsp",   Hsp\[Transpose][[2]]},
                                                                                     {"Hsup", Hsup\[Transpose][[2]]}}, dfields, 
                                  "AllSites" -> True, "LineLimit"->OptionValue["LineLimit"]];

  MMA2FORTRAN[dir<>"/src/hamiltonian/" <> "ForcesExpr", ExprFunc];

  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetSiteEnergyGs",  "SiteEnergy", {{"Gs",   Gs\[Transpose][[2]]}}, dfields, 
                       "LineLimit"->OptionValue["LineLimit"]];
  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetSiteEnergyGp",  "SiteEnergy", {{"Gp",   Gp\[Transpose][[2]]}}, dfields, 
                       "LineLimit"->OptionValue["LineLimit"]];
  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetSiteEnergyGsp", "SiteEnergy", {{"Gsp", Gsp\[Transpose][[2]]}}, dfields];

  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetSiteEnergyG",   "SiteEnergy", {{"Gs",   Gs\[Transpose][[2]]},
                                                                                  {"Gp",   Gp\[Transpose][[2]]},
                                                                                  {"Gsp", Gsp\[Transpose][[2]]}}, dfields, 
                       "LineLimit"->OptionValue["LineLimit"]];

  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetSiteEnergyHu",  "SiteEnergy", {{"Hu",  FixDoubleCounting[Hu\[Transpose][[2]]]}}, dfields, 
                       "LineLimit"->OptionValue["LineLimit"]];
  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetSiteEnergyHp",  "SiteEnergy", {{"Hp",  FixDoubleCounting[Hp\[Transpose][[2]]]}}, dfields, 
                       "LineLimit"->OptionValue["LineLimit"]];
  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetSiteEnergyHup", "SiteEnergy", {{"Hup", FixDoubleCounting[Hup\[Transpose][[2]]]}}, dfields, 
                       "LineLimit"->OptionValue["LineLimit"]];
  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetSiteEnergyHsu", "SiteEnergy", {{"Hsu", FixDoubleCounting[Hsu\[Transpose][[2]]]}}, dfields, 
                       "LineLimit"->OptionValue["LineLimit"]];
  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetSiteEnergyHsp", "SiteEnergy", {{"Hsp", FixDoubleCounting[Hsp\[Transpose][[2]]]}}, dfields, 
                       "LineLimit"->OptionValue["LineLimit"]];
  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetSiteEnergyHsup", "SiteEnergy", {{"Hsup", FixDoubleCounting[Hsup\[Transpose][[2]]]}}, dfields,
                       "LineLimit"->OptionValue["LineLimit"]];


  WriteLatticeModelF90[dir<>"/src/hamiltonian","GetSiteEnergyH",   "SiteEnergy", {{"Hu",   FixDoubleCounting[Hu\[Transpose][[2]]]},
                                                                                  {"Hp",   FixDoubleCounting[Hp\[Transpose][[2]]]},
                                                                                  {"Hup",  FixDoubleCounting[Hup\[Transpose][[2]]]},
                                                                                  {"Hsu",  FixDoubleCounting[Hsu\[Transpose][[2]]]},
                                                                                  {"Hsp",  FixDoubleCounting[Hsp\[Transpose][[2]]]},
                                                                                  {"Hsup", FixDoubleCounting[Hsup\[Transpose][[2]]]}}, dfields, 
                       "LineLimit"->OptionValue["LineLimit"]];

  ExprFunc = WriteLatticeModelF90[dir<>"/src/hamiltonian","GetSiteEnergy", "SiteEnergy",{{"Gs",   Gs\[Transpose][[2]]},
                                                                                         {"Gp",   Gp\[Transpose][[2]]},
                                                                                         {"Gsp", Gsp\[Transpose][[2]]},
                                                                                         {"Hu",   FixDoubleCounting[Hu\[Transpose][[2]]]},
                                                                                         {"Hp",   FixDoubleCounting[Hp\[Transpose][[2]]]},
                                                                                         {"Hup",  FixDoubleCounting[Hup\[Transpose][[2]]]},
                                                                                         {"Hsu",  FixDoubleCounting[Hsu\[Transpose][[2]]]},
                                                                                         {"Hsp",  FixDoubleCounting[Hsp\[Transpose][[2]]]},
                                                                                         {"Hsup", FixDoubleCounting[Hsup\[Transpose][[2]]]}}, dfields, 
                                  "LineLimit"->OptionValue["LineLimit"]];
  MMA2FORTRAN[dir<>"/src/hamiltonian/" <> "SiteEnergyExpr", ExprFunc];

  Return[{{Gs, Gp, Gsp, Hu, Hp, Hup, Hsu, Hsp, Hsup}, coefficients}]
]

ExportIO[dir_] := Module[{},
  MMA2FORTRAN[dir <> "/src/io/Constants", FortranConstants[]]
  MMA2FORTRAN[dir <> "/src/io/FileParser", FortranParsers[]];
  MMA2FORTRAN[dir <> "/src/io/Inputs",     FortranInputs[]];
  MMA2FORTRAN[dir <> "/src/io/Outputs",    FortranOutputs[]];
]

BuildLINVARIANT[dir0_, pos0_, fullmodel_, vars_, svars_, ie_, param_, ndipoles_, EwaldFieldQ_:True, OptionsPattern[{"LineLimit"->500}]] := Module[{dir, d, latt, sites, modelextend, coefficients, f90dir, angstrom2bohr, nn, ix, iy, iz},
  angstrom2bohr = 1.88973;
  {latt, sites} = pos0;
  dir = dir0 <> "/LINVARIANT";

  nn = Max[Max[#] & /@ Abs@Cases[Variables[fullmodel\[Transpose][[2]]], Subscript[__, __, ix_, iy_, iz_] -> {ix, iy, iz}]];

  f90dir = FileNameJoin[{$UserBaseDirectory, "Applications", "LINVARIANT", "f90"}];
  If[! DirectoryQ[dir], CreateDirectory[dir]];
  Run["cp " <> f90dir <> "/makefile " <> dir];
  Run["cp " <> f90dir <> "/makefile.include " <> dir];
  
  If[! DirectoryQ[dir <> "/src"], CreateDirectory[dir <> "/src"]; Run["chmod o-w " <> dir <> "/src"]]; 
  Run["cp " <> f90dir <> "/src/main.f90 " <> dir <> "/src"];
  Run["cp " <> f90dir <> "/src/symbol.inc " <> dir <> "/src"];
  Run["cp " <> f90dir <> "/src/makefile " <> dir <> "/src"];
  Run["cp " <> f90dir <> "/src/.objects " <> dir <> "/src"];

  Do[If[! DirectoryQ[dir <> "/" <> d], CreateDirectory[dir <> "/" <> d]; 
     Run["chmod o-w " <> dir <> "/" <> d]], {d, {"build", "bin", "example", "ext", "database"}}];

  Do[If[! DirectoryQ[dir <> "/src/" <> d], CreateDirectory[dir <> "/src/" <> d]; 
     Run["chmod o-w " <> dir <> "/src/" <> d]], {d, {"core", "io", "solvers", "hamiltonian", "fft", "common", "parallel", "xc"}}];

  Do[Run["cp "<> f90dir <> "/src/" <> d <> "/* " <> dir <> "/src/" <> d], {d, {"core", "io", "solvers", "hamiltonian", "fft", "common", "parallel", "xc"}}];

  Run["sed -i 's/NumIRFields/" <> ToString[ndipoles] <> "/g' "   <> dir <> "/src/solvers/mc.f90"];
  Run["sed -i 's/NumIRFields/" <> ToString[ndipoles] <> "/g' "   <> dir <> "/src/core/Ewald.f90"];
  Run["sed -i 's/neighbourcut/" <> ToString[nn] <> "/g' "        <> dir <> "/src/core/Ewald.f90"];
  Run["sed -i 's/NumIRFields/" <> ToString[ndipoles] <> "/g' "   <> dir <> "/src/core/LINVARIANT.f90"];
  Run["sed -i 's/FrozenQ           =.*/FrozenQ           = "<> "(\/" <> StringRiffle[ConstantArray[".false.", Length[vars]+1], ", "] <> "\/)" <>"/g' " <> dir <> "/src/io/" <> "/Inputs.f90"];

  {modelextend, coefficients} = ExportHamiltonian[dir, fullmodel, vars, svars, ie, EwaldFieldQ, "LineLimit"->OptionValue["LineLimit"]];
  MMA2FORTRAN[dir <> "/src/core/ReadCoefficients", FortranReadCoeff[Join[coefficients, param]]];
  MMA2FORTRAN[dir <> "/src/core/Parameters", FortranParamModule[Join[coefficients, param], vars]];
  MMA2FORTRAN[dir <> "/src/io/ExportModelData", FortranWriteMLModel[Join[coefficients, param], vars]];
  GenerateCoefficientsFile[dir, angstrom2bohr*latt, Join[coefficients, param]];

  Makefile[dir <> "/build/Makefile", dir];
  Return[{modelextend, Join[coefficients, param]}]
]

FitHoppingModelOld[dir_, pos0_, spg0_, tgrp_, TRules_, Basis_, nn_, qgrid_, vars_, svars_, varmap_, OptionsPattern[{"neps" -> 6, "epsmax" -> 0.03}]] := 
 Module[{Fij, neps, epsmax, LTB, lvar, mvar, Hdisp, HdispCoeff, AcousticBasis, Hacoustic, HacousticCoeff, uvar, Hdispeps, HdispepsCoeff, \[Sigma], v, i, j, ii, jj, CoeffFunc, dispeps11, dispeps12, Hepsu, ueps11, ueps12, HepsuCoeff, out, ux, uy, uz, farray},
  
  {lvar, mvar} = varmap;
  neps = OptionValue["neps"];
  epsmax = OptionValue["epsmax"];
  uvar = Subscript[ToExpression["u"], #] & /@ Range[3];
  AcousticBasis = (Normalize[#] & /@ (Flatten[Table[1.0 IdentityMatrix[3], {Length[pos0[[2]]]}], 1]\[Transpose]))\[Transpose];
  
  Fij = ReadForceConstants[dir <> "/supercell_" <> StringRiffle[qgrid, "x"] <> "/e11/e11.0/phonon_dfpt" <> "/FORCE_CONSTANTS"]; 
  LTB = FC2Fi0Bonds[Fij, pos0, qgrid][[1]];
  Hdisp = Jij[spg0, nn, TRules, {vars[[1]], vars[[1]]} /. Var2Var[lvar, mvar, 1], "OnSite" -> True] /. Var2Var[lvar, mvar, 2, "site" -> True];
  HdispCoeff = FitHessian[Basis, LTB, Hdisp, pos0, vars, "unit" -> "atomic"];
 
  {HacousticCoeff, Hacoustic} = Transpose@GetAcousticEnergy[spg0, tgrp, nn, LTB, pos0, uvar, "unit" -> "atomic"];
  farray = Prepend[ConstantArray[0.5, Length[Hacoustic] - 1], 1];
  HacousticCoeff[[1]] = First@HacousticCoeff - ImposeAcousticSumRule[{HacousticCoeff, Hacoustic}];
  
  Hdispeps = Jij[spg0, nn, TRules, {MonomialList[Expand[Total[vars[[1]]] Total[svars]]], vars[[1]]} /. Var2Var[lvar, mvar, 1], "OnSite" -> True] /. Var2Var[lvar, mvar, 2, "site" -> True];
  dispeps11 = Table[Fij = ReadForceConstants[dir <> "/supercell_" <> StringRiffle[qgrid, "x"] <> "/e11/e11." <> ToString[\[Sigma]] <> "/phonon_dfpt" <> "/FORCE_CONSTANTS"]; 
                    LTB = FC2Fi0Bonds[Fij, pos0, qgrid][[1]]; 
                    FitHessian[Basis, LTB, (First[#] /. {Subscript[ToExpression["\[Epsilon]"], __] -> 1}) & /@ Hdispeps, pos0, vars, "unit" -> "atomic"], {\[Sigma], -neps, neps}];
  dispeps12 = Table[Fij = ReadForceConstants[dir <> "/supercell_" <> StringRiffle[qgrid, "x"] <> "/e12/e12." <> ToString[\[Sigma]] <> "/phonon_dfpt" <> "/FORCE_CONSTANTS"]; 
                    LTB = FC2Fi0Bonds[Fij, pos0, qgrid][[1]]; 
                    FitHessian[Basis, LTB, (First[#] /. {Subscript[ToExpression["\[Epsilon]"], __] -> 1}) & /@ Hdispeps, pos0, vars, "unit" -> "atomic"], {\[Sigma], -neps, neps}];
  
  HdispepsCoeff = Table[CoeffFunc = If[D[Hdispeps[[i, 1]] /. {Subscript[ToExpression["P"], __] -> 1, Subscript[ToExpression["\[Epsilon]"], ii_, jj_, __, __, __] -> Subscript[ToExpression["\[Epsilon]"], ii, jj]}, Subscript[ToExpression["\[Epsilon]"], 1, 1]] != 0, LinearModelFit[{epsmax/neps Range[-neps, neps], dispeps11\[Transpose][[i]]}\[Transpose], v, v], LinearModelFit[{epsmax/neps Range[-neps, neps], dispeps12\[Transpose][[i]]}\[Transpose], v, v]]; 
                        CoeffFunc["BestFitParameters"][[2]], {i, Length[Hdispeps]}];
  
  Hepsu = GetAcousticEps[spg0, tgrp, nn, pos0, uvar, "unit" -> "atomic"]; 
  ueps11 = Table[Fij = ReadForceConstants[dir <> "/supercell_" <> StringRiffle[qgrid, "x"] <> "/e11/e11." <> ToString[\[Sigma]] <> "/phonon_dfpt" <> "/FORCE_CONSTANTS"]; 
                 LTB = FC2Fi0Bonds[Fij, pos0, qgrid][[1]]; 
                 FitHessian[AcousticBasis, LTB, (First[#] /. {Subscript[ToExpression["\[Epsilon]"], __] -> 1}) & /@ Hepsu, pos0, {uvar}, "unit" -> "atomic"], {\[Sigma], -neps, neps}];
  ueps12 = Table[Fij = ReadForceConstants[dir <> "/supercell_" <> StringRiffle[qgrid, "x"] <> "/e12/e12." <> ToString[\[Sigma]] <> "/phonon_dfpt" <> "/FORCE_CONSTANTS"]; 
                 LTB = FC2Fi0Bonds[Fij, pos0, qgrid][[1]]; 
                 FitHessian[AcousticBasis, LTB, (First[#] /. {Subscript[ToExpression["\[Epsilon]"], __] -> 1}) & /@ Hepsu, pos0, {uvar}, "unit" -> "atomic"], {\[Sigma], -neps, neps}]; 
  
  HepsuCoeff = Table[CoeffFunc = If[D[Hepsu[[i, 1]] /. {Subscript[ToExpression["u"], __] -> 1, Subscript[ToExpression["\[Epsilon]"], 1, 1, __, __, __] -> Subscript[ToExpression["\[Epsilon]"], 1, 1]}, Subscript[ToExpression["\[Epsilon]"], 1, 1]] != 0, LinearModelFit[{epsmax/neps Range[-neps, neps], ueps11\[Transpose][[i]]}\[Transpose], v, v], LinearModelFit[{epsmax/neps Range[-neps, neps], ueps12\[Transpose][[i]]}\[Transpose], v, v]]; 
                     CoeffFunc["BestFitParameters"][[2]], {i, Length[Hepsu]}];
  
  out = {{HdispCoeff, Hdisp}, {HacousticCoeff, Hacoustic}, {HdispepsCoeff, Hdispeps}, {HepsuCoeff, Hepsu}};
  
  Return[out]
]

FitHoppingModel[dir_, Fij_, pos0_, spg0_, tgrp_, TRules4Jij_, Basis_, nn_, qgrid_, vars_, svars_, ie_, OptionsPattern[{"sort"->True, "rotate"->True, "SymmetricBasisQ" -> True, "Hu" -> {}, "Hup" -> {}}]] := Module[{neps, epsmax, HJij, Hp, Hu, Hup, Hsu, Hsp, uvar, eu, x, dd, \[Sigma], v, \[Alpha], \[Beta], i, j, k, ii, jj, t, CoeffFunc, dispeps11, dispeps12, Hepsu, ueps11, ueps12, HepsuCoeff, out, ux, uy, uz, inv, farray, FullBasis, terms, AcousticBasis, Heff, BoolQ, param, expr},
  
  uvar = Subscript[ToExpression["u"], #] & /@ Range[3];
 
  AcousticBasis = (Normalize[#] & /@ (Flatten[Table[1.0 IdentityMatrix[3], {Length[pos0[[2]]]}], 1]\[Transpose]))\[Transpose];
  FullBasis = Transpose@Join[Basis\[Transpose], AcousticBasis\[Transpose]];

  Heff = FitSupercellHessian[FullBasis, Fij, qgrid, TRules4Jij, pos0, Append[vars, uvar], 
                             "sort"->OptionValue["sort"], "SymmetricBasisQ" -> OptionValue["SymmetricBasisQ"], "rotate" -> OptionValue["rotate"]];

  Hu   = Select[Heff, And@@XInvQ[#[[2]], DeleteDuplicates[First[#] & /@ uvar], ContainsExactly] &];
  Hp   = Select[Heff, And@@XInvQ[#[[2]], DeleteDuplicates[First[#] & /@ uvar], ContainsNone] &];
  Hup  = Complement[Heff, Join[Hu, Hp]];
 
  Heff = {If[OptionValue["sort"], ReverseSortBy[Hu,  Abs[#[[1]]]&], SortBy[Hu,  #[[2]]&]],
          If[OptionValue["sort"], ReverseSortBy[Hp,  Abs[#[[1]]]&], SortBy[Hp,  #[[2]]&]],
          If[OptionValue["sort"], ReverseSortBy[Hup, Abs[#[[1]]]&], SortBy[Hup, #[[2]]&]]};
 
  Hp = Heff[[2]];

  Hu = If[!(OptionValue["Hu"] === {}),
          {param, expr} = OptionValue["Hu"];
          HJij = {param, 
                  Table[terms = Sum[inv /. {Subscript[ToExpression["\[Epsilon]"], i_, j_] -> Subscript[ToExpression["\[Epsilon]"], i, j, dd[[1]], dd[[2]], dd[[3]]]} 
                                        /. Normal[StrainFromu[ToExpression["\[Epsilon]"],ToExpression["u"],ie]], {dd, Tuples[Range[-nn, 0], {3}]}];
                        GetJijStar[Expand@terms, 0], {inv, expr}]}\[Transpose];
          CollectPolynomial[Expand@Total[Times@@@HJij], TRules4Jij, "round"->1.0*10^-6, "sort"->OptionValue["sort"]],
          Heff[[1]]];
  
  Hup = If[!(OptionValue["Hup"] === {}),
           {param, expr} = OptionValue["Hup"];
           HJij = {param/Length[Tuples[Range[-nn, 0], {3}]], 
                   Table[terms = Sum[inv /. {Subscript[ToExpression["\[Epsilon]"], i_, j_] -> Subscript[ToExpression["\[Epsilon]"], i, j, dd[[1]], dd[[2]], dd[[3]]], 
                                             Subscript[v_, i_] -> Subscript[v, i, dd[[1]], dd[[2]], dd[[3]]]} 
                                         /. Normal[StrainFromu[ToExpression["\[Epsilon]"], ToExpression["u"], ie]], {dd, Tuples[Range[-nn, 0], {3}]}];
                         GetJijStar[Expand@terms, 0], {inv, expr}]}\[Transpose];
            CollectPolynomial[Expand@Total[Times@@@HJij], TRules4Jij, "round"->1.0*10^-6, "sort"->OptionValue["sort"]],
            Heff[[3]]];
  
  out = {Hp, Hu, Hup};

  Return[out]
]

FitSHoppingModel[dir0_, pos0_, spg0_, tgrp_, TRules4Jij_, Basis_, nn_, qgrid_, vars_, svars_, ie_, npt_, maxeps_, OptionsPattern[{"chop"->1.0*10^-6, "originshift"->{0,0,0}, "rotate" -> True, "SymmetricBasisQ" -> True}]] := Module[{x, t, i, j, d, eps, Hepp, eijsub, Fij, Hij, ts, param, model, expr3, out, fcfile, eta}, 

  eijsub = Flatten@Table["e" <> ToString[i] <> ToString[j] -> Subscript[ToExpression["\[Epsilon]"], i, j], {i, 3}, {j, 3}];
  eta = {"e11", "e22", "e33", "e23", "e13", "e12"};

  expr3 = Sum[Hepp = Table[fcfile = dir0 <> "/supercell_" <> StringRiffle[qgrid, "x"] <> "/FORCE_CONSTANTS_" <> eps <> "." <> ToString[d];
                           Fij = OriginShiftFC[pos0, ReadForceConstants[fcfile], qgrid, OptionValue["originshift"], "rotate" -> OptionValue["rotate"]];
                           Hij = FitHoppingModel[dir0, Fij, pos0, spg0, tgrp, {}, Basis, nn, qgrid, vars, svars, ie, 
                                                 "sort" -> False, 
                                                 "rotate" -> OptionValue["rotate"], 
                                                 "SymmetricBasisQ" -> OptionValue["SymmetricBasisQ"]];
                           Table[{d maxeps/npt, #1, #2} & @@@ Model2P1[Hij[[i]]], {i, 3}], {d, -npt, npt}]\[Transpose];
              
              Table[Chop@Expand@Total@Table[ts = If[Length[t] == 2 npt + 1, t, Append[t, {0, 0, t[[1, 3]]}]];
                                            param = LinearModelFit[ts[[;; , 1 ;; 2]], x, x, IncludeConstantBasis -> False]["BestFitParameters"][[1]];
                                            model = {param, ts[[1, 3]] eps /. eijsub};
                                            Expand[Chop[Times @@ model]], {t, GatherBy[Flatten[Hepp[[i]], 1], Last]}], {i, 3}], {eps, eta}];
 
  out = CollectPolynomial[#, TRules4Jij, "sort" -> True, "chop"->OptionValue["chop"], "SymmetricQ"->True] & /@ expr3;
  Return[out]
]

ImposeAcousticSumRule[model0_] := Module[{model1, v, a, ix, iy, iz, out, sub, Etot},
  sub = {Subscript[v_, a_, ix_, iy_, iz_] -> Subscript[v, a, 0, 0, 0]};
  Etot = Total[Times @@@ ({model0\[Transpose][[1]], FixDoubleCounting[model0\[Transpose][[2]]]}\[Transpose])];
  model1 = {-#1, #2} & @@@ CollectPolynomial[Chop@Expand[Etot] /. sub];
  out = Model2List[AdjustModel[{model0\[Transpose]}, {model1\[Transpose]}]];
  Return[out]
]

ImposeAcousticSumRuleOld[model_] := Module[{HacousticCoeff, Hacoustic, terms, out, residule, ind, t, i, j, k, v, x},
  {HacousticCoeff, Hacoustic} = model\[Transpose];
  terms = {(# /. Subscript[__] -> 1), SimplifyCommonFactor[#, 10^-6]} & /@ Level[Expand[HacousticCoeff . Hacoustic], {1}];
  residule = Merge[Table[ind = (Variables[t[[2]]] /. Subscript[_, i_, __] -> i);
                         ind = If[Length[ind] == 1, {First@ind, First@ind}, ind];
                         ind -> t[[1]], {t, terms}], Mean];

  out = SortBy[CollectPolynomial[Sum[ind = (Variables[t] /. {Subscript[v_, i_, __] -> i});
                                     Chop[t[[1]] - residule[If[Length[ind] == 1, {First@ind, First@ind}, ind]]] t[[2]], {t, terms}], {}], 
               (Variables[#[[2]]] /. {Subscript[_, x_, __] -> {x}}) &];

  Return[out]
]

UpdatePT[dir_, file_, WriteQ_, OptionsPattern[{"round" -> 1}]] := Module[{ptdata, ene, func, \[CapitalDelta]x, Tlist, npt, i, out, plt},
  Run["grep Epot " <> dir <> "/LINVARIANT.PTMC.out | tail -n 120 | awk '{print $6, $9}' > " <> dir <> "/" <> file];
  ptdata = ReadList[dir <> "/" <> file, Number, RecordLists -> True];
  npt = Length[ptdata];
  Tlist = Sort[ptdata\[Transpose][[1]]];
  ene = Sort[10^-6 RandomReal[{-1, 1}] + # &/@ (ptdata\[Transpose][[2]])];
  func = Interpolation[MapIndexed[{#1, Tlist[[First@#2]]} &, ene], InterpolationOrder -> 1];
  plt = Grid@{{Plot[func[x], {x, Min[ene], Max[ene]}, 
                    Frame -> True, PlotRange -> All, GridLines -> Automatic, AspectRatio->1, ImageSize -> 300,
                    FrameTicksStyle -> Directive[Black, 12], 
                    FrameLabel -> (Style[#, Bold, Black, 12] & /@ {"Energy (Hartree)", "Temperature (K)"})], 
           ListPlot[MapIndexed[{Tlist[[First@#2]], #1} &, ene], 
                    Frame -> True, PlotRange -> All, GridLines -> Automatic, AspectRatio->1, ImageSize -> 300,
                    FrameTicksStyle -> Directive[Black, 12], 
                    FrameLabel -> (Style[#, Bold, Black, 12] & /@ {"Temperature (K)", "Energy (Hartree)"})]}};

  \[CapitalDelta]x = First@Differences@MinMax[ene]/(npt - 1);
  out = Table[{N@Round[func[Min[ene] + (i - 1) \[CapitalDelta]x], OptionValue["round"]]}, {i, npt}];
  If[WriteQ, Export[dir <> "/REPLICAS.dat", out]];
  Print[plt];
  Return[{Flatten@out, Tlist}]
]

SqueezeCell[dir0_, pos0_, nstrain_, ds_, qgrid_] := Module[{\[Delta], s, i, latt, latt0, sites0, eij, exx, eta, svar},
  {latt0, sites0} = pos0;
  svar = ToExpression["Epsilon"];
  eij = Normal[SparseArray[{{i_, j_} /; i == j -> Subscript[svar, i, j], {i_, j_} /; i < j -> Subscript[svar, i, j], {i_, j_} /; i > j -> Subscript[svar, j, i]}, {3, 3}]];
  eta = eij2eta[eij];

  Table[\[Delta] = ds i;
        latt = (IdentityMatrix[3] + eij /. {s -> \[Delta]} /. {Subscript[First[s], __] -> 0}) . latt0;
        exx = "e" <> ToString[s[[2]]] <> ToString[s[[3]]];
        ExportPOSCAR[dir0 <> "/supercell_" <> StringRiffle[qgrid, "x"] <> "/" <> exx, "POSCAR." <> exx <> "." <> ToString[i] <> ".vasp", {latt, sites0}], {s, eta}, {i, -nstrain, nstrain}];
]

FromJij2Energy[model_, grid_, PBCS_] := Module[{Nx, Ny, Nz, Heff0, ene, v, a, x0, y0, z0, dx, dy, dz, Etot, sub},
  {Nx, Ny, Nz} = grid;
  sub = {Subscript[v_, a_, x0_, y0_, z0_] -> Subscript[v, a, If[PBCS[[1]], Mod[x0 + dx, Nx, 1], x0 + dx], 
                                                             If[PBCS[[2]], Mod[y0 + dy, Ny, 1], y0 + dy], 
                                                             If[PBCS[[3]], Mod[z0 + dz, Nz, 1], z0 + dz]]};
  Heff0 = Expand[model\[Transpose][[1]] . FixDoubleCounting[model\[Transpose][[2]]]];
  Etot = Sum[Heff0 /. sub, {dx, Nx}, {dy, Ny}, {dz, Nz}];
  Return[Etot]
]

OriginShiftFC[pos0_, Fij0_, NSC_, shift_, OptionsPattern[{"rotate" -> True}]] := Module[{Nij, latt, sites, Nx, Ny, Nz, NeighborList, i, j, k, s, ppos, spos, NumPpos, NumSpos, p2s, Fij, ind},
  Nij = DiagonalMatrix[NSC];
  {latt, sites} = pos0;
  {Nx, Ny, Nz} = NSC; 
  NeighborList = If[OptionValue["rotate"], 
                    Flatten[Table[{k, j, i} - {1, 1, 1}, {i, Nz}, {j, Ny}, {k, Nx}], 2], 
                    Flatten[Table[{i, j, k} - {1, 1, 1}, {i, Nx}, {j, Ny}, {k, Nz}], 2]];
  
  ppos = {pos0[[1]], {Mod[#1 + shift, {1, 1, 1}], #2} & @@@ (pos0[[2]])};
  spos = If[OptionValue["rotate"], 
            {Nij . latt, Flatten[Table[{#[[1]] + {k - 1, j - 1, i - 1}, #[[2]]}, {i, Nz}, {j, Ny}, {k, Nx}] & /@ (ppos[[2]]), 3]}, 
            {Nij . latt, Flatten[Table[{#[[1]] + {i - 1, j - 1, k - 1}, #[[2]]}, {i, Nx}, {j, Ny}, {k, Nz}] & /@ (ppos[[2]]), 3]}];
  
  NumPpos = Length[ppos[[2]]];
  NumSpos = Length[spos[[2]]];
  
  p2s = Association[Table[{i,s}->First@First@Position[Rationalize[Chop@Norm[Mod[sites[[i]][[1]]+s+shift,NSC]-#1]] &@@@ (spos[[2]]), 0], {i, NumPpos}, {s, NeighborList}]];
  
  ind = Flatten[Table[p2s[{i, s}], {i, NumPpos}, {s, NeighborList}]];
  Fij = Chop@Table[Fij0[[i, j]], {i, ind}, {j, ind}];
  
  Return[Fij]
]

TrainingSet2Supercell[ts_, grid_] := Module[{tgrid, nts, ndft, enescale, i, j, s, shiftlist, supervarsub, varsub, dgx, dgy, dgz, out},
  nts = Length[ts];
  out = Table[tgrid = TSGrid[ts[[i]]];
              If[! And @@ (IntegerQ[#] & /@ (grid/tgrid)), Print["ERROR: ts is not a divider of supercell!"]; Abort[]];
              {dgx, dgy, dgz} = grid/tgrid;
              enescale = (Times @@ grid)/(Times @@ tgrid);
              shiftlist = Flatten[Table[{ix, iy, iz} - {1, 1, 1}, {iz, dgz}, {iy, dgy}, {ix, dgx}], 2];
              ndft = Length[ts[[i]]];
              Table[varsub = ts[[i, j]][[5]];
                    supervarsub = Flatten@Table[(VarGridShift[#1, s, grid, True] -> #2) & @@@ varsub, {s, shiftlist}];
                    supervarsub = SortVarSub[supervarsub];
                    Join[ts[[i, j]][[1 ;; 3]], enescale ts[[i, j]][[4 ;; 4]], {supervarsub}, ts[[i, j]][[6 ;; 6]]], {j, ndft}], {i, nts}];
  Return[out]
]

FoldGamma[gmodel_, ah_] := Module[{ahsub, grid, out, v, a, ix, iy, iz},
  grid = Length@DeleteDuplicates[#] & /@ Transpose[DeleteDuplicates[VarGrid[#] & /@ Variables[Transpose[Model2List[gmodel]][[2]]]]];
  ahsub = Thread[FoldHopping[ah, grid] -> ah];
  out = {#1, If[OnSiteExprQ[#2], SimplifyCommonFactor[#2 /. {Subscript[v_, a_, ix_, iy_, iz_] -> Subscript[v, a]}, 10^-9], SimplifyCommonFactor[#2, 10^- 9]]} & @@@ Model2List[gmodel /. ahsub];
  Return[{out\[Transpose]}]
]

LMesh[dir_, grid_, mfunc_, OptionsPattern[{"write"->True, "plot"->True}]] := Module[{NGx, NGy, NGz, box, boxdata, pltdata, ix, iy, iz, occ, bcs, boxdict},
  {NGx, NGy, NGz} = grid;
  pltdata = Table[{occ, bcs} = mfunc[ix, iy, iz];
                  If[occ == 1, 1, If[First@bcs == 1, -1, 0]], {iz, 1, NGz}, {iy, 1, NGy}, {ix, NGx}];
  If[OptionValue["plot"], Print[ArrayPlot3D[pltdata, Mesh -> True, PlotRange -> {All, All, All}, ColorRules -> {1->Gray, -1->Blue}]]];
  
  box = Table[Flatten@{ix, iy, iz, mfunc[ix, iy, iz]}, {iz, NGz}, {iy, NGy}, {ix, NGx}];
  boxdict = GroupBy[Flatten[box, 2], #[[4]]&];
  boxdata = Prepend[Join[boxdict[1], boxdict[0]], grid];
  If[OptionValue["write"], Export[dir <> "/Supercell.inp", boxdata, "Table", "FieldSeparators" -> " "]];
  Return[pltdata]
]

TrainLINVARIANT[Plan_, vars_, ts_, OptionsPattern[{"InitRange" -> {-1, 1}, "NonlinearSteps" -> 10000, "MaxIterationSteps" -> 10^4, "BoundDist" -> 1, "TapeRate" -> 100, "BatchSize" -> 100, "Beta" -> {0.9, 0.999}, "Epsilon" -> 10.^-8, "XCut" -> 1.5, "tol" -> 10^-6, "FuncPatch"->0}]] := Module[{H, param, order, basis, bound, TrainSet, LossHistory, coeff, func, FuncPatch, mcoeff, vcoeff, GradLoss, batch, nsteps, grad, loss, batchsize, \[Rho], \[Beta], \[Epsilon], i, monitor, LossPlot, snapshot, FittingPlot, Iteration, TrainLoss, ValLoss, sol, path, plotticks, t, v, time0, time1, tmp, tol, XCut},
  
  FuncPatch  = OptionValue["FuncPatch"];
  batchsize  = OptionValue["BatchSize"];
  XCut       = OptionValue["XCut"];
  tol        = OptionValue["tol"];
  \[Beta]    = OptionValue["Beta"];
  \[Epsilon] = OptionValue["Epsilon"];
 
  path = Subsequences[Prepend[Accumulate[Plan\[Transpose][[4]]], 0], {2}];
  func        = FuncPatch;
  TrainSet    = ts;
  LossHistory = {};
  Iteration   = 0; 
  TrainLoss   = 0; 
  ValLoss     = 0;
  time0       = 0;
  time1       = 0;
  
  (*Define the dynamic plot*)
  LossPlot = Dynamic[ListLogPlot[LossHistory, 
                                 PlotRange -> All, 
                                 PlotLabel -> "Training Loss", 
                                 Frame -> True, 
                                 Joined -> True, 
                                 ColorFunction -> Function[{xx, yy}, Hue@Rescale[First@Flatten@Position[If[Between[xx, #], 1, 0] & /@ path, 1], {1, Length[Plan]+0.1}]], 
                                 ColorFunctionScaling -> False, 
                                 FrameLabel -> {"Iterations", "Loss"}, 
                                 ImageSize -> 300, 
                                 ImagePadding -> {{60, 10}, {30, 10}}]];
  
  FittingPlot = Dynamic[ListPlot[{Table[{t, Quiet[func /. Thread[vars -> {t}]]}, {t, TrainSet\[Transpose][[1]]}], TrainSet}, 
                                 PlotStyle -> {Red, Blue}, 
                                 Joined -> {True, False}, 
                                 PlotRange -> All, 
                                 Frame -> True, 
                                 ImageSize -> 300, 
                                 ImagePadding -> {{60, 10}, {30, 10}},
                                 Filling -> {1 -> {{2}, {LightBlue, LightRed}}}]];

  monitor = DynamicModule[{iteration = 0, trainLoss = 0, valLoss = 0},
      Monitor[Do[{order, basis, bound, nsteps, \[Rho]} = Plan[[i]];
                 {param, H} = GaussianPolynomial[order, "x", "basis" -> basis];
                 TrainSet = Join[ts, BoundData[ts, bound, 2]];
                 mcoeff = ConstantArray[0, Length@param];
                 vcoeff = ConstantArray[0, Length@param];
                 coeff = RandomReal[OptionValue["InitRange"], Length@param];
                 sol = Dispatch[Thread[param -> coeff]];

                 {time0, tmp}= AbsoluteTiming[GradLoss = Table[batch = RandomSample[TrainSet, batchsize];
                                                               loss  = GetLossFunc[FuncPatch + H, vars, batch];
                                                               grad  = Grad[loss, param];
                                                               {grad, loss}, {j, nsteps}];];
      
                 Do[{time1, tmp} = AbsoluteTiming[iteration = iteration + 1;
                                                  {grad, loss} = GradLoss[[j]];
       
                                                  {coeff, mcoeff, vcoeff} = UpdateADAM[coeff, mcoeff, vcoeff, Quiet[grad/.sol], j, \[Rho], \[Beta], \[Epsilon]];
                                                  sol = Dispatch[Thread[param -> coeff]];
                                                 ];

                    If[Mod[iteration, OptionValue["TapeRate"]] == 0, Iteration = iteration;
                                                                     trainLoss = Quiet[loss /. sol];
                                                                     valLoss = Quiet[loss /. sol];
                                                                     TrainLoss = trainLoss;
                                                                     ValLoss = valLoss;
                                                                     func = Quiet[FuncPatch + (H /. sol)];
                                                                     AppendTo[LossHistory, {iteration, trainLoss}]], {j, nsteps}];

                 sol = Quiet[NonlinearModelFit[TrainSet, FuncPatch + H, Transpose[{param, coeff}], vars,
                                               MaxIterations -> OptionValue["NonlinearSteps"]]]["BestFitParameters"];
                 sol = Dispatch[#1 -> Round[#2, 10.^-6] & @@@ sol];
                 FuncPatch = ChopGaussianPolynomial[Quiet[FuncPatch + (H /. sol)], vars, TrainSet, XCut, tol], {i, Length@Plan}],

              Panel[Column[{Row[{Style["SGD - with Adaptive Moment Estimation", Black, Bold, 12]}],
                            Row[{Grid[{{"Iteration: ",     Dynamic[Iteration]}, 
                                       {"Training Loss: ", Dynamic[TrainLoss]}, 
                                       {"Validation Loss: ", Dynamic[ValLoss]},
                                       {"Timing Gradients:", Dynamic[time0]},
                                       {"Seconds/Batch:", Dynamic[time1]}}]}],
                            Row[{Dynamic[func]}],
                            Row[{LossPlot, FittingPlot}]}], Style["Training Progress", Red, 14]]
     ]];
  
  monitor;
 
  FittingPlot = ListPlot[{Table[{t, Quiet[FuncPatch /. Thread[vars -> {t}]]}, {t, TrainSet\[Transpose][[1]]}], TrainSet}, 
                         PlotStyle -> {Red, Blue}, 
                         Joined -> {True, False}, 
                         PlotRange -> All, 
                         Frame -> True, 
                         ImageSize -> 300, 
                         ImagePadding -> {{60, 10}, {30, 10}}, 
                         Filling -> {1 -> {{2}, {LightBlue, LightRed}}}];

  snapshot = Panel[Column[{Row[{Style["SGD - with Adaptive Moment Estimation", Black, Bold, 12]}],
                           Row[{Grid[{{"Iteration: ", Iteration}, 
                                      {"Training Loss: ", TrainLoss}, 
                                      {"Validation Loss: ", ValLoss}}]}],
                           Row[{FuncPatch}],
                           Row[{LossPlot, FittingPlot}]}], Style["Training Progress", Red, 14]];
  Print[snapshot];
  
  Return[{FuncPatch, snapshot}]
]

GaussianPolynomial[order_, var_, OptionsPattern[{"basis" -> {1, 1, 1}}]] := Module[{A, B, C, a, b, c, v, n, param, func, s1, s2, s3, OrderList},
  OrderList = If[ListQ[order], order, Range[2, order, 2]];
  v = ToExpression[var];
  {s1, s2, s3} = OptionValue["basis"];
  {param, func} = Table[A = ToExpression["A" <> ToString[n]];
                        B = ToExpression["B" <> ToString[n]];
                        C = ToExpression["C" <> ToString[n]]; 
                        a = ToExpression["a" <> ToString[n]]; 
                        b = ToExpression["b" <> ToString[n]]; 
                        c = ToExpression["c" <> ToString[n]];
                        {{s1 A, s2 B, s2 C, s3 a, s3 b, s3 c}, 
                        s1 A v^n + s2 B v^n Exp[- C^2 v^2/1] + s3 a v^n (Exp[-b^2 ((v - c)^2/1)] + Exp[-b^2 ((v + c)^2/1)])}, {n, OrderList}]\[Transpose];

  Return[{DeleteCases[Flatten[param], 0], Total@func}]
]

GetLossFunc[H_, vars_, ts_, OptionsPattern[{"tol"->10^-6}]] := Module[{loss},
  loss = Quiet@Expand[Mean[(Chop[H /. Thread[vars -> #[[1 ;; -2]]]] - Last[#])^2 & /@ ts]];
  Return[Chop[loss, OptionValue["tol"]]]
]

UpdateADAM[Param_, mParam_, vParam_, grad_, iter_, \[Rho]_, \[Beta]_, \[Epsilon]_] := Module[{\[Beta]1, \[Beta]2, mParamNew, vParamNew, mCorrected, vCorrected, ParamNew, i},
  {\[Beta]1, \[Beta]2} = \[Beta];
  mParamNew = \[Beta]1 mParam + (1 - \[Beta]1) grad;
  vParamNew = \[Beta]2 vParam + (1 - \[Beta]2) grad^2;
  Quiet[mCorrected = mParamNew/(1 - \[Beta]1^iter)];
  Quiet[vCorrected = vParamNew/(1 - \[Beta]2^iter)];
  ParamNew = Table[Param[[i]] - \[Rho] mCorrected[[i]]/(Sqrt[vCorrected[[i]]] + \[Epsilon]), {i, Length@Param}];
  Return[{ParamNew, mParamNew, vParamNew}]
  
]

BoundData[data_, npts_, dist_ : 1] := Module[{i, x, func, FuncMax, FuncMin, f, x0, x1, bound},
  bound = Table[func = Fit[f[data, First, 2], {1, x}, x];
                {x1, x0} = Transpose[f[data, First, 2]][[1]]; 
                Table[{x, func} /. {x -> x1 + i dist (x1 - x0)}, {i, npts}], {f, {MaximalBy, MinimalBy}}];
  Return[Flatten[bound, 1]]
]

ChopGaussianPolynomial[expr_, vars_, ts_, IntRange_, tol_ : 10^-6] := Module[{terms, TermsChopped, residue, range, t, plots},
  range = Flatten[{#1, #2}] & @@@ Thread[vars -> (IntRange MinMax[#] & /@ Transpose[#[[1 ;; -2]] & /@ ts])];
  terms = Level[Expand@expr, {1}];
  TermsChopped = Table[residue = Chop[Quiet[NIntegrate @@ Join[{Abs[t]}, range]], tol]; 
                       If[residue == 0, ## &[], t], {t, terms}];
  Return[Total[TermsChopped]]
]

(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
