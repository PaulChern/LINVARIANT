BeginPackage["LINVARIANT`Builder`", {"LINVARIANT`Structure`", "LINVARIANT`Vasp`", "LINVARIANT`INVARIANT`", "LINVARIANT`MathematicaPlus`",  "LINVARIANT`Parser`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
GetPhononBasis            ::usage "GetPhononBasis[pos0, Fij]"
ShowPhononBasis           ::usage "ShowPhononBasis[pos0, phonon]"
PhononDecomposition       ::usage "PhononDecomposition[pos0, pos, phonon]"
ExportPhononBasis         ::usage "ExportPhononBasis[dir, cif, pos0, phonon]"
PhononLinearRecombination ::usage "PhononLinearRecombination[PhononBasis, rules]"
FitTrainingSet            ::usage "FitTrainingSet[goalfunction]"
GenRandTrainSet           ::usage "GenRandTrainSet[ref, gs, basis, BasisLabel, modes, NTS, range]"
GetTrainingSet            ::usage "GetTrainingSet[ref, phases, spg, rho]"
WeightTrainSet            ::usage "WeightTrainSet[data, tau]"
ExportTrainingSet         ::usage "ExportTrainingSet[dir, data, tau]"
BoundHighOrderTerm        ::usage "BoundHighOrderTerm[invariants, expr, vars, cuts]"
PlanFitting               ::usage "PlanFitting[invariants, vars, spgmat, ts, seeds]"
LinvariantFit             ::usage "LinvariantFit[invariants, vars, order, ts, spgmat]"
FetchSeeds                ::usage "FetchSeeds[inv, vars, n]"
SeedsByInvariants         ::usage "SeedsByInvariants[invariants, vars]"
BasisComplement           ::usage "BasisComplement[BasisMatrix, FullBasis]"
ImportDMDPath             ::usage "ImportDMDPath[file, pos0, Basis, Basislabel]"
SamplePhaseSpace          ::usage "SamplePhaseSpace[dir0, s0, spgmat, Basis, potim]"
CollectSeeds              ::usage "CollectSeeds[dir0, nphase, Basis, OpMat]"
CollectPhases             ::usage "CollectPhases[dir0, n, vars, Basis]"
ShowSeedPhases            ::usage "ShowSeedPhases[phases, seeds, vars]"
SamplePolynomialSpace     ::usage "SamplePolynomialSpace[polyphase, dir, spgmat0, Basis, npt, smearing, mim0]"
PlanSampling              ::usage "PlanSampling[dir0_, phases, invariants, vars, Basis, spgmat, npt, smearing]"
SamplingByDMDPath         ::usage "SamplingByDMDPath[dir0, DMDPath, vars, svars, BasisMatrix, FullBasis]"
SampleAround              ::usage "SampleAround[pos0, dir, fname, spgmat0, Basis, npt, potim]"
ImposePhase               ::usage "ImposePhase[pos0, dir, fname, phase, Basis]"
RoundPhase                ::usage "RoundPhase[phase]"
LoadTrainingset           ::usage "LoadTrainingset[dir0, BasisMatrix, FullBasis]"
PreparePhases             ::usage " PreparePhases[dir0, pname, num]"
GetHam                    ::usage "GetHam[models]"
AdjustHam                 ::usage "AdjustHam[models, ham, vars, svars]"
TrimmedCoefficient        ::usage "TrimmedCoefficient[inv, vars, coeff, trim]"
CheckMinimums             ::usage "CheckMinimums[models, seeds, amp, vars, svars]"
DMDPath2ts                ::usage "DMDPath2ts[dmd, vars, svars]"
Expression2ModelData      ::usage "Expression2ModelData[ham, vars, svars]"
GetActivePolynomial       ::usage "GetActivePolynomial[modeldata, vars, svars, mim]"
SortModelData             ::usage "SortModelData[modeldata, vars, svars]"
JoinModels                ::usage "JoinModels[m0, m1]"
CheckFitting              ::usage "CheckFitting[invariants, model, ts, vars, svars]"
GetTopOrder               ::usage "GetTopOrder[polynomials, vars, svars]"
TopOrderQ                 ::usage "TopOrderQ[p, polynomials, vars, svars]"

(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
GetPhononBasis[pos0_, Fij_, OptionsPattern[{"table" -> True, "fontsize" -> 12, "toiso"->None}]] := Module[{sol, PhononByIR, PhononRotated, ind, frequency, IR, i, j, iv, eigenvectors, rot, NumAtom, tab, simplifytab, color, isolabels, Vasp2THz = 15.633302300230191},
  NumAtom = Length[pos0[[2]]];
  simplifytab = {0 -> "-"};
  sol = Chop[Eigensystem[ArrayFlatten[Fij]]]\[Transpose];
  PhononByIR = Gather[MapIndexed[Join[#2, #1] &, sol], Chop[#1[[2]] - #2[[2]]] == 0 &];
  isolabels = If[OptionValue["toiso"]===None, 
                 MapIndexed[ConstantArray["IR"<>ToString[First@#2], #1]&, Length[#]&/@PhononByIR],
                 CloneReshape2D[PhononByIR, OptionValue["toiso"]]];
  PhononRotated = Table[{ind, frequency, IR} = PhononByIR[[i]]\[Transpose];
                        eigenvectors = DeleteCases[Chop@Orthogonalize@Flatten[Table[rot = PseudoInverse[IR\[Transpose][[3 (j - 1) + 1 ;; 3 (j - 1) + 3]]\[Transpose]]; Normalize[#] & /@ (rot . IR), {j, Length[pos0[[2]]]}], 1], v_ /; Norm[v] == 0]; 
                        Table[{ind[[iv]], isolabels[[i,iv]], frequency[[iv]], eigenvectors[[iv]]}, {iv, Length[ind]}], {i, Length[PhononByIR]}];
  tab = Grid[Join[{Join[{"#","IR","frequency (THz)"}, Flatten[Table[Subscript[#, \[Alpha]], {\[Alpha], {"x", "y", "z"}}] & /@ (pos0[[2]]\[Transpose][[2]])]]}, Flatten[#] /. simplifytab & /@ Flatten[PhononRotated, 1]], 
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
   tab = Grid[Join[{Join[{"#","IR","frequency (THz)"}, Flatten[Table[Subscript[#, \[Alpha]], {\[Alpha], {"x", "y", "z"}}] & /@ (pos0[[2]]\[Transpose][[2]])]]}, Flatten[#] /. simplifytab & /@ Flatten[phonon, 1]],
              Alignment -> Center,
              Frame -> True,
              Background -> {Flatten[{Blue, Blue, Blue, Table[color = If[OddQ[i], White, None];
                                                        ConstantArray[color, 3], {i, NumAtom}]}],
                             Flatten[{LightBlue, Table[color = If[OddQ[i], Yellow, Pink];
                                                       ConstantArray[color, Length[phonon[[i]]]], {i, Length@phonon}]}]},
              ItemSize -> Full,
              ItemStyle -> Directive[FontSize -> OptionValue["fontsize"], Black]];
   Print[tab];
]


ExportPhononBasis[dir_, cif_, sites0_, phonon_] := Module[{i, phononordered, cifdata, latt, sites},
  cifdata = ImportIsodistortCIF[cif];
  latt = N[cifdata[[6]] /. cifdata[[8]]];
  sites = cifdata[[10]];
  phononordered = SortIR2POSCAR[latt, sites, sites0, phonon\[Transpose]]\[Transpose];

  Table[Export[dir <> "/phonon-" <> ToString[i] <> ".dat", cif2mcif[cif, i, phononordered, sites]];
        Run["mv " <> dir <> "/phonon-" <> ToString[i] <> ".dat " <> dir <> "/phonon-" <> ToString[i] <> ".mcif"], {i, Length@phonon}];
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
  ind = Sort[Flatten[#\[Transpose][[1;;2]]]] & /@ rules;
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
  exprnew = SimplifyCommonFactor[NOrderResponse[expr, vars, 1]];
  invnew = SimplifyCommonFactor[NOrderResponse[expr, vars, 1] /. cuts];
  If[DisjointQ[Flatten[invariants], {invnew,-invnew}], BoundHighOrderTerm[invariants, exprnew, vars, cuts], exprnew /. cuts]
]

RoundPhase[phase_] := Module[{i},
  Table[If[Chop[phase[[i]]] == 0., 0, 1], {i, Length[phase]}]
]

SampleAround[pos0_, dir_, fname_, spgmat0_, Basis_, npt_, potim_, OptionsPattern[{"round" -> 10^-8}]] := Module[{BasisDim, NumBasis, modesets, samples, modes, pos, pos1, phase, m, i, spgmat},
  {BasisDim, NumBasis} = Dimensions[Basis];
  pos = ImportPOSCAR[dir <> "/" <> fname];
  phase = ISODISTORT[pos0[[1]], pos0[[2]], pos[[2]], Basis, Range[NumBasis], "round" -> OptionValue["round"]]\[Transpose][[4]];
  spgmat = Table[If[Chop[Norm[m . Sign[phase] - Sign[phase]], OptionValue["round"]] == 0, m, ## &[]], {m, spgmat0}];
  
  modesets = DeleteDuplicates[Flatten[Table[SortBy[Tuples[{0, -i, i}, {NumBasis}], Norm[N@#] &], {i, Range[0, npt]}], 1]];
  samples = Round[DeleteDuplicates[Table[{i, # . modesets[[i]] & /@ spgmat}, {i, Length[modesets]}], Complement[#1[[2]], #2[[2]]] == {} &]\[Transpose][[2]]\[Transpose][[1]]];
  Print["Number of samples: " <> ToString[Length[samples]]];
  Do[modes = Table[{m, "modes", 1, potim[[m]] samples[[i, m]]}, {m, NumBasis}];
     pos1 = {pos[[1]], ImposeMode[pos[[1]], pos[[2]], Basis, modes, 1.0]}; 
     ExportPOSCAR[dir, "sample" <> ToString[i] <> ".vasp", pos1], {i, Length@samples}]
]

ImposePhase[pos0_, dir_, fname_, phase_, Basis_] := Module[{m, modes, BasisDim, NumBasis, pos},
  {BasisDim, NumBasis} = Dimensions[Basis];
  modes = Table[{m, "modes", 1, phase[[m]]}, {m, NumBasis}];
  pos = {pos0[[1]], ImposeMode[pos0[[1]], pos0[[2]], Basis, modes, 1.0]}; 
  ExportPOSCAR[dir, fname, pos]
]

FetchSeeds[inv_, vars_, n_, OptionsPattern[{"zero" -> False, "log"->0}]] := Module[{character, ndigits, shortseed, list, map, seeds0, order}, 
  character = First@InvariantCharacter[inv, vars];
  list = If[OptionValue["log"]==0,
            If[OptionValue["zero"], Range[n,-n,-Sign[n]], DeleteCases[Range[n,-n,-Sign[n]],0]],
            Join[-Reverse[OptionValue["log"]^Range[1,n-1]], OptionValue["log"]^Range[1,n-1]]];
  order = Length[list];
  ndigits = Count[character, 1];
  shortseed = Tuples[list, ndigits]\[Transpose];
  map = Flatten[Position[character, 1]];
  seeds0 = Table[If[i == 0, ConstantArray[0, order^ndigits], ## &[]], {i, character}];
  Do[seeds0 = Insert[seeds0, shortseed[[i]], map[[i]]], {i, Length@map}];
  Return[seeds0\[Transpose]]
]

SeedsByInvariants[invariants_, vars_, OptionsPattern[{"table"->True, "exseeds"->None, "odd"->False}]] := Module[{id, inv, sub, i, seeds, s, dict, sampleseeds, tab, tabdata, color, fullvars, simplifytab, varlabels, exseeds, invdict, models, t},
  fullvars = Flatten[vars];
  invdict = GroupBy[Flatten[invariants], StrainInvQ[#, vars[[-1]][[1, 1]]] &];
  models = Flatten[Table[SortInvariants[t, fullvars], {t, {invdict[False], invdict[True]}}]];

  varlabels = MapIndexed["(" <> ToString[First@#2] <> ")" &, fullvars];
  simplifytab = {-1 -> Style["-1", Bold, Blue], 1 -> Style["1", Bold, Red], 0 -> "-"};

  exseeds = If[OptionValue["exseeds"]===None, 
               {},
               {0,"f["<>StringJoin[Riffle[ToString[#, StandardForm] & /@ DeleteCases[fullvars Abs[#], 0], ","]] <> "]", #} &/@ OptionValue["exseeds"]];

  seeds = Prepend[Flatten[Table[inv = If[Head[models[[id]]] === Plus, First[models[[id]]], models[[id]]];
                        s = FetchSeeds[models[[id]], fullvars, 1, "zero"->False];
                        dict = Values@GroupBy[Table[sub = Thread[fullvars -> If[OptionValue["odd"], s[[i]], Abs@s[[i]]]]; 
                                                    {inv /. sub, s[[i]]}, {i, Length[s]}], First];
                        {id, #1 inv, #2} & @@@ First[dict\[Transpose]], {id, Length@models}], 1],
                  {0, "f[0]", ConstantArray[0,Length[fullvars]]}];
  sampleseeds = Join[DeleteDuplicates[seeds, #1[[3]]==#2[[3]]&], exseeds];
  tabdata = {#1, #2, #3 /. simplifytab} & @@@ sampleseeds;
  tab = Grid[Prepend[Prepend[MapIndexed[Flatten[{#2, #1}] &, tabdata], Flatten[{{"", "", ""}, fullvars}]], Flatten[{{"", "", ""}, varlabels}]],
             Background -> {Flatten[{LightBlue, White,LightBlue, 
                            Table[color = If[OddQ[i], Yellow, Pink];
                                  ConstantArray[color, Length[vars[[i]]]], {i, Length@vars}]}], None}];
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

ImportDMDPath[dir0_, vars_, svars_, BasisMatrix_, FullBasis_, seeds_, OptionsPattern[{"fontsize"->12, "ts"->{}, "table"->True, "pdf"->True, "plot"->True, "interpolate"->0}]] := Module[{run, vasprun, lattices, sites, forces, stress, ene0, ene, coordinates, labels, color, i, path, epsplot, scfplot, roseplot, out, eta, fullvars, OtherBasisMatrix, ind, exlabels, FullBasisMatrix, ExInfo, plots, tsdata, data, ts, ibrion, inv, plt, pos0},
  pos0 = ImportPOSCAR[dir0 <> "/POSCAR0"];
  FullBasisMatrix = Flatten[FullBasis, 1]\[Transpose][[4]]\[Transpose];
  OtherBasisMatrix = Complement[FullBasisMatrix\[Transpose], BasisMatrix\[Transpose], SameTest -> (Chop[Norm[#1 - #2]] == 0 &)];
  ind = Sort[First@First[Position[FullBasisMatrix\[Transpose], #]] & /@ OtherBasisMatrix];
  ExInfo = Tally[First@First@Position[FullBasis, #] & /@ ind];
  OtherBasisMatrix = (FullBasisMatrix\[Transpose][[#]] & /@ ind)\[Transpose];
  fullvars = Append[vars, svars];
  labels = MapIndexed["("<>ToString[First@#2]<>") "<>ToString[#1,StandardForm]&, Flatten@fullvars];
  exlabels = "("<>ToString[#]<>") "<>(Flatten[FullBasis, 1]\[Transpose][[2]][[#]]) &/@ ind;

  data = Import[dir0 <> "/seeds/" <> "ene.dat", "Table"];
  ene0 = data[[1,2]];
  tsdata = If[OptionValue["ts"] === {}, data, If[AllTrue[OptionValue["ts"], MemberQ[data[[;;,1]], #]&], 
                                                 Extract[data, First@Position[data\[Transpose][[1]], #]] &/@ OptionValue["ts"],
                                                 Print["You many need to collect the seeds again before importing the seeds data!"];
                                                 Abort[]]];
  plots = {};

  out = Table[
        ts = tsdata[[i, 1]];

        inv = ToString[seeds[[ts,2]], StandardForm];
        vasprun = Table[Import[dir0<>"/seeds/seed"<>ToString[ts]<>"/vasprun_ibrion"<>ToString[ibrion]<>".xml"], {ibrion, {3, 1}}];
        sites = Flatten[Table[ParseXMLData[ParseXML[ParseXML[run, "structure", {}], "varray", {"name", "positions"}], "v"], {run, vasprun}], 1];
        forces = Flatten[Table[ParseXMLData[ParseXML[run, "varray", {"name", "forces"}], "v"], {run, vasprun}], 1];
        stress = Flatten[Table[ParseXMLData[ParseXML[run, "varray", {"name", "stress"}], "v"], {run, vasprun}], 1];
        lattices = Flatten[Table[ParseXMLData[ParseXML[ParseXML[run, "structure", {}], "varray", {"name", "basis"}], "v"], {run, vasprun}], 1];
        eta = N@Round[eij2eta[GetStrainTensor[pos0[[1]], #, "iso" -> False]], 10^-4]&/@lattices;
        ene = Flatten[Table[ParseFortranNumber@Flatten[ParseXML[Last[ParseXML[#, "energy", {}]], "i", {"name", "e_fr_energy"}] & /@ ParseXML[run, "calculation", {}]], {run, vasprun}], 1] - ene0;
        coordinates = Table[Join[Flatten[{ISODISTORT[pos0[[1]], pos0[[2]], {sites[[i]], pos0[[2]]\[Transpose][[2]]}\[Transpose], BasisMatrix, Range[Length[BasisMatrix\[Transpose]]], "round"->10.0^-3]\[Transpose][[4]], eta[[i]]}], ISODISTORT[pos0[[1]], pos0[[2]], {sites[[i]], pos0[[2]]\[Transpose][[2]]}\[Transpose], OtherBasisMatrix, exlabels, "round"->10.0^-3]\[Transpose][[4]]], {i, Length[sites]}]\[Transpose];
        path = Grid[Prepend[Prepend[Prepend[coordinates, ene], Range[Length[ene]]]\[Transpose], Join[{"#", "energy"}, labels, exlabels]],
                    Background -> {Join[{Gray,LightBlue}, Flatten[Table[color = If[OddQ[i], Yellow, Pink];
                                   ConstantArray[color, Length[fullvars[[i]]]], {i, Length@fullvars}]], 
                                   Flatten[Table[ConstantArray[If[OddQ[i], Gray, LightGray], ExInfo[[i,2]]], {i,Length[ExInfo]}]]],
                                   None},
                    ItemStyle -> Directive[FontSize -> OptionValue["fontsize"]],
                    ItemSize -> Full
                   ]; 
        epsplot = Partition[Table[If[Norm[coordinates[[i]]] != 0, ListPlot[{coordinates[[i]], ene}\[Transpose], Frame -> True, GridLines -> Automatic, PlotRange -> All, PlotLabel -> labels[[i]], ImageSize -> 300], ## &[]], {i, Length[Flatten@fullvars]}], UpTo[4]];
        scfplot = {ListLinePlot[coordinates[[1;;Length[Flatten@vars]]], PlotRange -> All, Frame -> True, PlotMarkers -> Automatic, 
                                                   GridLines -> Automatic, ImageSize -> 300, ImagePadding -> {{50, 1}, {20, 1}}, 
                                                   PlotLegends -> LineLegend[labels[[;;Length[Flatten@vars]]], LegendLayout -> {"Column", 3}]], 
                         ListLinePlot[coordinates[[Length[Flatten@vars]+1;;Length[Flatten@vars]+6]], PlotRange -> All, Frame -> True, PlotMarkers -> Automatic,
                                                    GridLines -> Automatic, ImageSize -> 300, ImagePadding -> {{50, 1}, {20, 1}},
                                                    PlotLegends -> LineLegend[labels[[Length[Flatten@vars]+1;;]], LegendLayout -> {"Column", 2}]],
                          ListLinePlot[coordinates[[Length[Flatten@vars]+7;;]], PlotRange -> All, Frame -> True, PlotMarkers -> Automatic,
                                                     GridLines -> Automatic, ImageSize -> 300, ImagePadding -> {{50, 1}, {20, 1}},
                                                     PlotLegends -> LineLegend[exlabels, LegendLayout -> {"Column", 3}]],
                         ListLinePlot[ene, PlotRange -> All, Frame -> True, PlotMarkers -> Automatic, GridLines -> Automatic, 
                                           ImageSize -> 300, ImagePadding -> {{50, 1}, {20, 1}}]};
        roseplot = Partition[Rasterize[Show[ListPlot3D[{coordinates[[#1]], coordinates[[#2]], ene}\[Transpose], 
                                                       AxesLabel -> {labels[[#1]], labels[[#2]], ""}, 
                                                       PlotRange -> All, 
                                                       ColorFunction -> "Rainbow", 
                                                       InterpolationOrder -> OptionValue["interpolate"], 
                                                       MeshFunctions -> {#3 &}, 
                                                       Mesh -> Length[ene], 
                                                       ViewPoint -> {0, 0, 100}, 
                                                       ImageSize -> 300], 
                                            ListPointPlot3D[{coordinates[[#1]], coordinates[[#2]], ene}\[Transpose], 
                                                       PlotStyle -> {White, PointSize[Large]}, 
                                                       AxesLabel -> {labels[[#1]], labels[[#2]],""}, 
                                                       PlotRange ->  All, 
                                                       ViewPoint -> {0, 0, 100}, ImageSize -> 300], 
                                            ParametricPlot3D[BSplineFunction[{coordinates[[#1]],coordinates[[#2]],ene}\[Transpose]][s], {s, 0, 1}, 
                                                       PlotStyle -> Directive[White, Thick], 
                                                       PlotRange -> All, 
                                                       ViewPoint -> {0, 0, 100}]], 
                                                       ImageSize->300] 
                   & @@@ Subsets[Table[If[Norm[coordinates[[i]]] == 0, ## &[], i], {i, Length[Flatten@fullvars]}], {2}], UpTo[4]];
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

SamplingByDMDPath[dir0_, DMDPath_, vars_, svars_, BasisMatrix_, FullBasis_, seeds_, OptionsPattern[{"spline" -> 0, "sparse" -> 10^-2, "sigma" -> 1.1, "extension"->0, "table"->False, "fontsize"->12, "write"->True, "plot"->True, "imagesize"->200, "interpolate"->3}]] := Module[{mim0, character, scharacter, fullvars, excharacter, miminfo, data, i, j, k, samples, DMDSamples, InterpolatedSamples, nvars, dft, sigma, sparse, plt, plots, pos0, pos, FullBasisMatrix, ComplementBasisMatrix, tcount, scount, dir, shiftmode, shifteij, latt, sites, PickQ, imd, dmd, inv, exlabels, labels, info, numsamples, data2d, enedata, exi, exene, exvars, exsvars, exovars, extension, exseeds, s, s1, bsplinepath, nsamples, splinedata, tab, color, slavemode, primarymode, hiddenmode, mask, primarymask, slavemask, hiddenmask},
  If[! DirectoryQ[dir0 <> "/trainingset"], CreateDirectory[dir0 <> "/trainingset"]; Run["chmod o-w " <> dir0 <> "/trainingset"]];
  fullvars = Append[vars, svars];
  pos0 = ImportPOSCAR[dir0 <> "/POSCAR0"];
  FullBasisMatrix = Flatten[FullBasis, 1]\[Transpose][[4]]\[Transpose];
  {ComplementBasisMatrix, exlabels} = BasisComplement[BasisMatrix, FullBasis];
  labels = MapIndexed["("<>ToString[First@#2]<>") "<>ToString[#1,StandardForm]&, Flatten@fullvars];
  FullBasisMatrix = Join[BasisMatrix\[Transpose], ComplementBasisMatrix\[Transpose]]\[Transpose];
  nvars = Length[Flatten[vars]];

  sigma = OptionValue["sigma"];
  sparse = OptionValue["sparse"];

  plots = {};
  info = {};
  numsamples = 0;
  Do[scount = 0;
     tcount = DMDPath[[imd]][[1]];
     inv = ToString[seeds[[tcount,2]], StandardForm];
     dir = dir0 <> "/trainingset/ts-" <> ToString[tcount];
     If[! DirectoryQ[dir], CreateDirectory[dir]; Run["chmod o-w " <> dir]];
     dft = DMDPath[[imd]][[2]]\[Transpose];
     (*
     character = Sign[Abs@dft[[1]][[2 ;; nvars + 1]]];
     scharacter = Sign[Abs@dft[[1]][[nvars + 2 ;; nvars + 7]]];
     excharacter = Sign[Mean[Abs@dft][[nvars + 8 ;;]]];
     dft = Join[{1}, character, scharacter, excharacter] # & /@ (DMDPath[[imd]][[2]]\[Transpose]);
     *)
     mim0 = Mean[MinimalBy[dft, First]];
     mim0 = Last@dft;

     mask = Sign@Mean@Abs[dft[[;;,2;;]]];
     primarymask = Sign@Abs[First[dft][[2;;]]];
     hiddenmask = Join[ConstantArray[0, nvars+6], ConstantArray[1, Length[Flatten@exlabels]]];
     slavemask = (mask - primarymask)*(1-hiddenmask);
     hiddenmode = hiddenmask*mim0[[2;;]];

     miminfo = Table[{i - 1, mim0[[i]], Sqrt@Mean[(# - mim0[[i]])^2 & /@ (dft\[Transpose][[i]])]}, {i, nvars + 7}]; 
     samples = Table[PickQ = Fold[And, Table[Abs[dft[[i]][[j + 1]] - miminfo[[j + 1, 2]]] <= sigma*miminfo[[j + 1, 3]], {j, nvars + 6}]];
                     If[PickQ, 
                        primarymode = primarymask*dft[[i,2;;]];
                        slavemode = MapIndexed[If[Chop[#1]==0., 0, Extract[slavemask*dft[[i,2;;]], #2]]&, slavemask*mim0[[2;;]]];
                        shiftmode = primarymode + slavemode + hiddenmode;
                        {i, 
                         dft[[i]][[1]], 
                         Table[N@Round[shiftmode[[j]], sparse], {j, nvars}], 
                         Table[N@Round[shiftmode[[nvars + j]], 0.001], {j, 6}],
                         shiftmode[[nvars + 7 ;;]]}, 
                        ## &[]], {i, Length[dft]}];


     exene = dft[[1]][[1]] - dft[[2]][[1]];
     exvars = (Sign[dft[[1]][[2;;nvars+7]]]*(miminfo\[Transpose][[3]][[2;;]]))[[1;;nvars]];
     exsvars = (Sign[dft[[1]][[2;;nvars+7]]]*(miminfo\[Transpose][[3]][[2;;]]))[[nvars+1;;]];
     exovars = Chop[dft[[-1]][[nvars + 8 ;;]], sparse];

     extension = Table[{Length[dft]+s, dft[[1]][[1]] + exene*s, exvars*sigma*s, exsvars*sigma*s, exovars}, {s, 1, OptionValue["extension"]}];

     DMDSamples = Join[extension,
                       DeleteDuplicates[samples, Chop[Norm[Flatten[#1[[3 ;; 4]] - #2[[3 ;; 4]]]]] == 0 &],
                       {{Last[samples][[1]]+1, Last[dft][[1]], Last[dft][[2;;nvars+1]], Last[dft][[nvars+2;;nvars+7]], Last[dft][[nvars+8;;]]}}];

     bsplinepath = BSplineFunction[Join[{#2}, #3, #4, #5] &@@@ DMDSamples];

     InterpolatedSamples = If[OptionValue["spline"] == 0, 
                              DMDSamples, 
                              nsamples = Ceiling[Length[DMDSamples]*OptionValue["spline"]];
                              Table[splinedata = bsplinepath[N[i/nsamples]];
                                    {i, splinedata[[1]], splinedata[[2;;nvars+1]], splinedata[[nvars+2;;nvars+7]], splinedata[[nvars+8;;]]}, {i, 0, nsamples}]
                              ];
     PrependTo[InterpolatedSamples, {0, First[dft][[1]], First[dft][[2;;nvars+1]], First[dft][[nvars+2;;nvars+7]], mim0[[nvars + 8 ;;]]}];
     PrependTo[InterpolatedSamples, {0, First[dft][[1]], First[dft][[2;;nvars+1]], First[dft][[nvars+2;;nvars+7]], First[dft][[nvars + 8 ;;]]}];

     tab = Grid[Prepend[Flatten[#] &/@ InterpolatedSamples, Join[{"#", "energy"}, labels, Flatten@exlabels]],
                Background -> {Join[{Gray,LightBlue}, Flatten[Table[color = If[OddQ[i], Yellow, Pink];
                               ConstantArray[color, Length[fullvars[[i]]]], {i, Length@fullvars}]],
                               Flatten[Table[ConstantArray[If[OddQ[i], Gray, LightGray], Length[exlabels[[i]]]], {i,Length[exlabels]}]]],
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
        shiftmode = Table[{i, "modes", 1, s[[3, i]]}, {i, nvars}];
        shiftmode = Join[shiftmode, Table[{nvars + i, "modes", 1, s[[5, i]]}, {i, Length[Flatten@exlabels]}]];
        shifteij = Chop@eta2eij[s[[4]]];
        latt = (shifteij + IdentityMatrix[3]).(pos0[[1]]);
        sites = ImposeMode[pos0[[1]], pos0[[2]], FullBasisMatrix, shiftmode, 1.0];
        pos = {latt, sites};
        If[OptionValue["write"],  ExportPOSCAR[dir, "sample." <> ToString[scount] <> ".vasp", pos]], {s, InterpolatedSamples}];

     AppendTo[info, "ts-" <> ToString[tcount] <> " (" <> ToString[scount] <> "): " <> inv], {imd, Length@DMDPath}];
  
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

SamplePhaseSpace[dir0_, s0_, spgmat_, Basis_, potim_, OptionsPattern[{"seeds" -> {}}]] := Module[{BasisDim, NumBasis, modesets, seeds, modes, m, i, pos, pos0, eij, strainlatt},
  If[! DirectoryQ[dir0<>"/seeds"<>"-"<>ToString[potim]], 
     CreateDirectory[dir0<>"/seeds"<>"-"<>ToString[potim]]; 
     Run["chmod o-w " <> dir0<>"/seeds"<>"-"<>ToString[potim]]];
  pos0 = ImportPOSCAR[dir0 <> "/" <> s0];
  {BasisDim, NumBasis} = Dimensions[Basis];

  If[OptionValue["seeds"] == {}, 
     modesets = SortBy[Tuples[{-1, 0, 1}, {NumBasis}], Norm[N@#] &];
     seeds = Round[DeleteDuplicates[Table[Sort[# . modesets[[i]] & /@ spgmat], {i, Length[modesets]}]]\[Transpose][[1]]],
     seeds = OptionValue["seeds"]
  ];

  Print["Number of seeds: " <> ToString[Length[seeds]]<>" (potim="<>ToString[potim]<>")"];
  Do[eij=potim*Chop@eta2eij[seeds[[i]][[-6;;]]]/Norm[pos0[[1]]];
     strainlatt=(eij+IdentityMatrix[3]).(pos0[[1]]);
     modes = Table[{m, "modes", 1, potim seeds[[i, m]]}, {m, NumBasis}]; 
     pos = {strainlatt, ImposeMode[pos0[[1]], pos0[[2]], Basis, modes, 1.0]}; 
     ExportPOSCAR[dir0<>"/seeds"<>"-"<>ToString[potim], "seed." <> ToString[i] <> ".vasp", pos];
     ExportOPTCELL[dir0<>"/seeds"<>"-"<>ToString[potim]<>"/OPTCELL." <> ToString[i] <> ".inp", Abs[seeds[[i]][[-6;;]]], "constrain"->False], {i, Length@seeds}]
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
                 eta = N@Round[eij2eta[GetStrainTensor[pos0[[1]], pos[[1]], "iso" -> False]], 10^-4];
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

CollectPhases[dir0_, vars_, svars_, BasisMatrix_, FullBasis_, OptionsPattern[{"round" -> 10^-6, "fontsize" -> 12, "table" -> True, "file"->"POSCAR"}]] := Module[{NumDim, NumBasis, pos0, Basislabel, i, mpdata, pos, sites, eta, ExPhase, phase, simplifytab, fullvars, tab, ind, OtherBasisMatrix, ExBasislabel, color, fname, FullBasisMatrix, ExInfo, data, t, ene, ene0, isg0, sg0, isg, sg, sginfo},
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
  tab = Grid[Prepend[Flatten[{Join[{#1, StringSplit[#2, "\[RightArrow]"][[2]], #3}, #4, #5, #6], Join[{" ", " ", " "}, Sign[#4], Sign[#5], Sign[#6]] /. simplifytab} & @@@ mpdata, 1], Join[{"#", "sg", "ene"}, Flatten@fullvars, ExBasislabel]],
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
  Grid[Prepend[Flatten[{Flatten[{#[[1]], StringSplit[#[[2]], "\[RightArrow]"][[2]], #[[3;;6]]}], 
                        Join[{"("<>ToString@seeds[[#[[1]],1]]<>") "<>ToString[Times@@DeleteCases[seeds[[#[[1]],3]],0] seeds[[#[[1]],2]], StandardForm]}, {StringSplit[#[[2]], "\[RightArrow]"][[1]], " "}, seeds[[#[[1]],3]] /. simplifytab, ConstantArray["-",Length[exlabels]]]} & /@ (phases\[Transpose][[1 ;; 6]]\[Transpose]), 1], Flatten[Join[{"#", "sg", "ene"}, labels, exlabels]]], 
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

LoadTrainingset[dir0_, BasisMatrix_, FullBasis_, vars_, svars_, OptionsPattern[{"round" -> 10^-6, "plot"->False, "pdf"->True, "imagesize"->200, "interpolate"->3}]] := Module[{i, t, s, BasisDim, NumBasis, ene, pos, mode0, vasprun, mode, exmode, data, eta, out, ws, seeds, plan, p, E0, pos0, FullBasisMatrix, ComplementBasisMatrix, exlabels, plt, pdf, fullvars, enedata, data2d},
  {BasisDim, NumBasis} = Dimensions[BasisMatrix];
  pos0 = ImportPOSCAR[dir0 <> "/POSCAR0"];
  fullvars = Flatten[{vars, svars}];
  FullBasisMatrix = Flatten[FullBasis, 1]\[Transpose][[4]]\[Transpose];
  {ComplementBasisMatrix, exlabels} = BasisComplement[BasisMatrix, FullBasis]; 
  FullBasisMatrix = Join[BasisMatrix\[Transpose], ComplementBasisMatrix\[Transpose]]\[Transpose];
  data = Import[dir0 <> "/trainingset/" <> "ene.dat", "Table"];
  E0 = data[[1, 3]];
  out = Table[t = data[[i, 1]]; 
              s = data[[i, 2]]; 
            ene = data[[i, 3]];
            pos = ImportPOSCAR[dir0 <> "/trainingset/ts" <> "-" <> ToString[t] <> "/sample" <> ToString[s] <> "/POSCAR"];
            eta = Chop[eij2eta[GetStrainTensor[pos0[[1]], pos[[1]], "iso" -> False]]];
           mode = ISODISTORT[pos0[[1]], pos0[[2]], pos[[2]], BasisMatrix, Range[NumBasis], "round" -> OptionValue["round"]]\[Transpose][[4]];
         exmode = ISODISTORT[pos0[[1]], pos0[[2]], pos[[2]], ComplementBasisMatrix, Flatten@exlabels, "round" -> OptionValue["round"]]\[Transpose][[4]];
         {i, t, s, ene - E0, mode, eta, exmode}, {i, Length[data]}];
  out = SortBy[#, First] & /@ GatherBy[out, #[[2]] &];

  pdf = Table[data2d = (Flatten[#] & /@ (t[[3;;]]\[Transpose][[5 ;; 6]]\[Transpose]))\[Transpose];
              enedata = t[[3;;]]\[Transpose][[4]];
              plt = Join[Table[If[Norm[t[[3;;]]\[Transpose][[5]]\[Transpose][[i]]] != 0, 
                                  ListPlot[{t[[3;;]]\[Transpose][[5]]\[Transpose][[i]], enedata}\[Transpose], 
                                           PlotRange -> All,
                                           ImageSize -> OptionValue["imagesize"], 
                                           Frame -> True, 
                                           GridLines -> Automatic, 
                                           PlotMarkers -> Automatic, 
                                           PlotLabel -> "(" <> ToString[i] <> ") " <> ToString[fullvars[[i]], StandardForm], 
                                           FrameLabel -> {"coordinates (\[Angstrom])", "energy (eV)"}, 
                                           ImagePadding -> {{60, 10}, {40, 10}}], ## &[]], {i, Length[Flatten@vars]}], 
                         Table[If[Norm[t[[3;;]]\[Transpose][[6]]\[Transpose][[i]]] != 0,
                                  ListPlot[{t[[3;;]]\[Transpose][[6]]\[Transpose][[i]], enedata}\[Transpose],
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
              {t[[1, 2]], plt}, {t, out}];

  If[OptionValue["pdf"],
     Export[dir0<>"/ts_info.pdf", CreateDocument[Grid[Join[{{Graphics@Text[Style["ts-"<>ToString[#[[1]]], Red, Bold, 12]]}},
                                                               Partition[#[[2]], UpTo[3]]], Frame -> All, FrameStyle -> Red] & /@ pdf,
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

PlanFitting[invariants_, modelfixed_, ts_, vars_, svars_, OpMat_, OptionsPattern[{"fontsize" -> 12, "table" -> True}]] := Module[{p, ip, invdict, v, v1, v2, t, i, j, m, models, tcharacter, vcharacter, characters, modes, plan0, plan1, plan2, plan3, fullvars, minimums, color, simplifytab, mim, tmask, dependency, PES, sol, PrimaryPolynomial, SlavePolynomial, HiddenPolynomial, tmp, sorted, len, foundQ, SilentQ, tab, prefix, polynomials, tsdata}, 
  simplifytab = {-1 -> Style["-1", Bold, Blue], 1 -> Style["1", Bold, Red], 0 -> "-"};
  fullvars = Append[vars, svars];
  models = Flatten[SortInvariants[#, Flatten@fullvars] & /@ Reverse@GatherBy[Flatten[invariants], StrainInvQ[#, svars[[1, 1]]] &], 1];
  prefix = (Flatten[#] &/@ (modelfixed\[Transpose]))\[Transpose];
  Print["Number of invariants : " <> ToString[Length[Flatten@models]]];
  plan0 = Table[characters = InvariantCharacter[m[[1]], Flatten@fullvars];
                {m, Table[tcharacter = Sign@Flatten[First[Abs[t\[Transpose][[5 ;; 6]]\[Transpose]]]];
                          If[MemberQ[characters, tcharacter], t, ## &[]], {t, ts}]}, {m, models}];
  plan1 = Gather[plan0, IntersectingQ[Flatten[Table[tmask = Flatten[If[Norm[#] != 0, ConstantArray[1, Length[#]], #] & /@ CloneReshape2D[vars, Last[t][[5]]]];
                                                    mim = Sign@Mean[#[[5]] tmask & /@ t];
                                                    Round[m] . mim, {t, #1[[2]]}, {m, OpMat}], 1], 
                                      Flatten[Table[tmask = Flatten[If[Norm[#] != 0, ConstantArray[1, Length[#]], #] & /@ CloneReshape2D[vars, Last[t][[5]]]];
                                                    mim = Sign@Mean[#[[5]] tmask & /@ t];
                                                    Round[m] . mim, {t, #2[[2]]}, {m, OpMat}], 1]] &];
  dependency = Table[PES = Flatten[{#5, #6, #4}] & @@@ Flatten[plan1[[mim]]\[Transpose][[2]], 2];
                     HiddenPolynomial = Table[vcharacter = InvariantCharacter[v, Flatten@fullvars];
                                              tcharacter = Sign@Abs@PES[[;; , 1 ;; -2]];
                                              SilentQ = FreeQ[Flatten@Table[v1 . v2 == Total[v1], {v1, vcharacter}, {v2, tcharacter}], True];
                                              If[SilentQ, v, ## &[]], {v, Flatten@invariants}];
                     sol = Quiet[LinearModelFit[PES, Complement[Flatten@invariants, HiddenPolynomial], Flatten@fullvars, IncludeConstantBasis -> False][{"BestFitParameters", "BasisFunctions"}]]\[Transpose];
                     PrimaryPolynomial = Flatten[plan1[[mim]]\[Transpose][[1]]];
                     SlavePolynomial = If[! MemberQ[PrimaryPolynomial, #2], #2, ## &[]] & @@@ sol;
                     {DeleteDuplicates[First@First@Position[Table[Flatten[plan1[[i]]\[Transpose][[1]]], {i, Length@plan1}], #] & /@ SlavePolynomial], mim}, {mim, Length@plan1}];
  
  tmp = dependency;
  sorted = Cases[tmp, {{}, __}];
  tmp = DeleteCases[tmp, {{}, __}];
  While[tmp != {}, 
        len = Length[tmp];
        foundQ = False;
        i = 1;
        While[! foundQ && i <= len, 
              If[ContainsAll[sorted\[Transpose][[2]], First[tmp[[i]]]], 
                 AppendTo[sorted, tmp[[i]]]; tmp = Drop[tmp, {i}]; foundQ = True, 
                 i = i + 1]];
        If[i > len && ! foundQ, Print["dependency error!"]; Abort[]];];

  plan2 = plan1[[#2]]\[Transpose] & @@@ sorted;
  
  tab = {};
  Print["Number of Potential Energy minimums: " <> ToString[Length[plan2]]];
  If[OptionValue["table"], 
     Do[minimums = DeleteDuplicates[Sign[Flatten[Last[#][[5 ;; 6]]]] & /@ Flatten[plan2[[i]][[2]], 1]];
        AppendTo[tab, {Style["Potential Energy minimum " <> ToString[i] <> ":", Bold, Black], Style["Polynomial Group " <> ToString[i] <> ":", Black, Bold]}];
        AppendTo[tab, {Grid[Prepend[minimums /. simplifytab, Flatten@fullvars], 
                            Background -> {Flatten[Table[color = If[OddQ[i], Yellow, LightBlue];
                                                         ConstantArray[color, Length[fullvars[[i]]]], {i, Length@fullvars}]], 
                                           Flatten[{Gray, ConstantArray[None, Length[minimums]]}]}, 
                            ItemStyle -> Directive[FontSize -> OptionValue["fontsize"]], 
                            ItemSize -> Full], 
                            Table[If[MemberQ[Flatten@prefix,j],Style[j,Darker@Gray],Style[j,Black]], {j, #}]&/@ (plan2[[i]][[1]])}], {i, Length[plan2]}]; 
     Print[Grid[tab, ItemSize -> Automatic, Frame -> All, Alignment -> {Center, Center}]]];
  
  Which[Length@Flatten[plan2\[Transpose][[1]]] < Length[Flatten[invariants]], 
        Print["Warning: not all polynomials are planned!"], 
        Length@Flatten[plan2\[Transpose][[1]]] > Length[Flatten[invariants]], 
        Print["Warning: some polynomials belong to different mimimums!"]];

  plan3 = Table[polynomials = Flatten[plan2[[ip]][[1]]];
                tsdata = Flatten[plan2[[ip]][[2]], 1];
                {{Transpose[MapIndexed[{ToExpression["LF"<>ToString[ip]<>"C"<>ToString[First@#2]], #1}&, Complement[polynomials, Flatten@prefix]]]},
                 DeleteCases[{(If[MemberQ[polynomials, #[[2]]], #, ## &[]] & /@ prefix)\[Transpose]}, {}], 
                 tsdata}, {ip, Length@plan2}];

  Return[plan3]
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

LinvariantFit[dir0_, invariants_, order_, vars_, svars_, ts_, OpMat_, TRules_, isovars_, OptionsPattern[{"exinvariants" -> {}, "round" -> 10^-6.0, "modelfixed" -> {{{}, {}}}, "alat" -> 1.0, "FontSize" -> 12, "constraints" -> 10.0, "imagesize" -> 200, "pdf" -> True, "plot" -> False, "interpolate" -> 3}]] := Module[{alat, NumBasis, submodel, modelorders, OffsetFunc, FitInfo, logdata, tsdata, seeds, plan, x, i, m, ip, p, t, v, log, tsout, modelout, minimumtest, mimene, sol, init, model4fit, ModelFuncFitted, modeldata, tsinfo, fullvars, checkQ, character, data2d, enemodel, enedata, plt, lvar, mvar, ExplosionList, ExInvariants, AddInvariants, infodata, ExplosionQ, AllInvariants, modeltmp, exorder, FoundQ, rsquare, coeff, tmp, FitData, modelfixed, submodelfixed, prefix, lfitvars, lfitpolynomials, constraints, BoundOrder, polyseeds, submodelbound, boundcount, lfitinit, PolyOrderTable},
  {lvar, mvar} = isovars;
  fullvars = Append[vars, svars];
  alat = OptionValue["alat"];
  NumBasis = Length[Flatten@vars];

  modelfixed = OptionValue["modelfixed"];
  prefix = (Flatten[#] & /@ (modelfixed\[Transpose]))\[Transpose];
  AllInvariants = DeleteDuplicates[Flatten[{invariants, OptionValue["exinvariants"]}]];
  
  modeldata = If[Flatten[modelfixed\[Transpose][[2]]] === {}, 
                 AllInvariants, 
                 Complement[AllInvariants, Flatten[modelfixed\[Transpose][[2]]]]];
  plan = PlanFitting[AllInvariants, modelfixed, ts, vars, svars, OpMat, "table" -> True];
  
  Print["Checking training set."];
  checkQ = False;
  Do[If[p[[3]] == {}, checkQ = True; Print["Error! Trainingset is not complete:"]; Print[p[[1]]]], {p, plan}];
  If[checkQ, Abort[], Print["Training set check passed"]];

  logdata = {};
  AddInvariants = {};
  modelout = {SortModelData[DeleteCases[{(If[MemberQ[Flatten[Flatten[plan\[Transpose][[2]],1]\[Transpose][[2]]], #2], ## &[], {#1, #2}] & @@@ prefix)\[Transpose]},{}], vars, svars]};
  tsout = {};
  
  Do[submodelbound = {};
     submodel = plan[[ip]][[1]];
     submodelfixed = plan[[ip]][[2]];
     (*OffsetFunc = (modelfitted /. Thread[Flatten[fullvars] -> {##}] &);*)
   
     tsdata = Flatten[{#5/alat, #6, #4}] & @@@ Flatten[plan[[ip]][[3]], 1];
     tsinfo = Grid[Partition[Style["ts-" <> ToString[#[[1, 2]]] <> " (" <> ToString[Length[#]] <> ")", Black, Bold] & /@ (plan[[ip]][[3]]), 1]];
   
     (*{rsquare, coeff} = If[submodel != {},
                           Quiet[LinearModelFit[tsdata, submodel, Flatten[fullvars], 
                                                IncludeConstantBasis -> False, 
                                                LinearOffsetFunction -> OffsetFunc][{"RSquared", "BestFitParameters"}]], 
                           {0, {}}];*)

     ExplosionQ = True;
     modeltmp = submodel;
     boundcount = 1;
     sol = {};
     While[ExplosionQ,
       submodel =JoinModels[submodelbound, submodel];
       {lfitvars, lfitpolynomials} = Flatten[#] &/@ Transpose[submodel];
       lfitinit = Transpose[{lfitvars, If[MatchQ[Head@#, Symbol], 0, #] &/@ (lfitvars /. sol)}];

       constraints = lfitvars[[#1]] >=0 &@@@ GetTopOrder[lfitpolynomials, vars, svars];
       constraints = If[TopOrderQ[#2, lfitpolynomials, vars, svars]||PolynomialOrder[#2, Flatten@fullvars]>order, 
                        #1 >= 0, 
                        ##&[]] &@@@ Transpose[{lfitvars, lfitpolynomials}];
       constraints = {};

       model4fit = Fold[JoinModels, Join[modelout, {submodelfixed, submodel}]];

       Print[{ip, boundcount, GetHam[submodelbound, vars, svars]}];
       {rsquare, sol} = If[submodel != {},
                           Quiet[NonlinearModelFit[tsdata, 
                                                   Flatten@{GetHam[model4fit, vars, svars], constraints}, 
                                                   lfitinit,
                                                   Flatten[fullvars]][{"RSquared", "BestFitParameters"}]], 
                           {0, {}}];

       sol = Thread[lfitvars -> N@Round[lfitvars /. sol, OptionValue["round"]]];
       modeltmp = SortModelData[JoinModels[submodelfixed, submodel] /. sol, vars, svars];

       coeff = lfitvars /. sol;
 
       polyseeds = DeleteDuplicates@Flatten[If[TopOrderQ[#2, lfitpolynomials, vars, svars] && #1 <= 0., 
                                               #2, 
                                               ##&[]] &@@@ Transpose[{coeff, lfitpolynomials}]];
 
       ExInvariants = If[polyseeds === {}, 
                         {},
                         Flatten[GetInvariants[TRules, 
                                               Variables[#]/.Var2Var[lvar,mvar,1], 
                                               {PolynomialOrder[#, Flatten@fullvars]+2}, "round" -> 10^-6] /. Var2Var[lvar,mvar,2] &/@ polyseeds]];

       submodelbound = {MapIndexed[If[NumPolynomialVar[#1]==1, 
                                      {ToExpression["LF"<>ToString[ip]<>"B"<>ToString[boundcount]<>"C"<>ToString[First@#2]], #1}, 
                                      ##&[]]&, ExInvariants]\[Transpose]};
       If[submodelbound === {{}}, ExplosionQ = False];
       boundcount = boundcount + 1;
       If[boundcount>10,Break[]]];

     AppendTo[modelout, modeltmp];
   
     FitInfo = Grid[Join[{Range[Length[Flatten[Flatten[modelout,1]\[Transpose][[2]]]]]},
                         {Flatten[Flatten[modelout,1]\[Transpose][[1]]]},
                         {If[MemberQ[Flatten[prefix], #], Style[#, Darker@Gray], Style[#, Black]] &/@ Flatten[Flatten[modelout,1]\[Transpose][[2]]]}]\[Transpose],
                   ItemSize -> Full, Alignment -> Left, Frame -> True, 
                   Background -> {{2 -> Lighter[Green, 0.5], 3 -> Lighter[Cyan, 0.5]}, None},
                   Dividers -> {{1 -> True}, Thread[(Accumulate[Length[Flatten[Transpose[#][[2]]]] & /@ modelout] + 1) -> True]}];
   
     AppendTo[tsout, plan[[ip]][[3]]]; 
   
     ModelFuncFitted = GetHam[Flatten[modelout,1], vars, svars];
   
     plt = CheckFitting[lfitpolynomials, Flatten[modelout,1], plan[[ip]][[3]], vars, svars, "plot" -> False, "alat" -> alat];

     seeds = Join[DeleteDuplicates@Flatten[FetchSeeds[#, Flatten[fullvars], 2, "zero" -> False] & /@ lfitpolynomials, 1],
                  DeleteDuplicates[Flatten[Last[#][[5 ;; 6]]] & /@ ts, Chop[Norm[#1 - #2] Norm[#1 + #2], OptionValue["round"]] == 0. &]];

     
     minimumtest = SortBy[DeleteDuplicates[Flatten[Table[init = {Flatten[fullvars], 
                                                                m s, 
                                                                ConstantArray[-OptionValue["constraints"], Length@Flatten[fullvars]],                       
                                                                ConstantArray[OptionValue["constraints"], Length@Flatten[fullvars]]}\[Transpose];
                                                        sol = Quiet@FindMinimum[ModelFuncFitted, init];
                                                        Prepend[Chop[Flatten[fullvars] /. (sol[[2]]), OptionValue["round"]], 
                                                                N@Round[sol[[1]], OptionValue["round"]]], {s, seeds}, {m, {0.5, 3.}}], 1], 
                                          Chop[Norm[Flatten[#1[[2 ;;]] - #2[[2 ;;]]]], OptionValue["round"]] == 0. &], First];
     mimene = Split[minimumtest\[Transpose][[1]]];
   
     infodata = Join[{FitInfo}, 
                     {tsinfo}, 
                     {NumberForm[Grid[Prepend[minimumtest, Join[{"ene"}, Flatten[fullvars]]], 
                                      Background -> {Flatten[{LightBlue, Table[ConstantArray[If[OddQ[v], Yellow, Pink], Length[Append[vars, svars][[v]]]], {v, Length@Append[vars, svars]}]}], 
                                                             Prepend[Flatten@Table[If[OddQ[i], ConstantArray[Gray, Length[mimene[[i]]]], ConstantArray[None, Length[mimene[[i]]]]], {i, Length[mimene]}], Cyan]}, 
                                      ItemSize -> Full], {4, 4}]}, 
                     {rsquare}, 
                     {Grid[Partition[plt, UpTo[10]]]}];
   
     AppendTo[logdata, infodata];
   
     ExplosionList = DeleteDuplicates[Table[If[MemberQ[Abs[minimumtest[[i]]][[2 ;; NumBasis + 1]], 10.], 
                                              If[Chop[Abs@# - 10.] == 0., 1, 0] & /@ (minimumtest[[i]][[2 ;; NumBasis + 1]]), 
                                              ## &[]], {i, Length[minimumtest]}]]; 
     Print[ExplosionList // MatrixForm], {ip, Length@plan}];
  
  log = Style[# /. x_Real :> Chop[x, OptionValue["round"]], OptionValue["FontSize"]] &@ Grid[{{"BestFit", "trainingset", "Minimums", "R02", "plots"}, ## & @@ logdata}, Dividers -> All, ItemSize -> Full];
  Print[log];

  Export[dir0 <> "/fitting_info.pdf", CreateDocument[Style[# /. x_Real :> Chop[x, OptionValue["round"]], OptionValue["FontSize"]] &@ Grid[{{"BestFit", "trainingset", "Minimums", "R02", "plots"}, #}\[Transpose], ItemSize -> Full] & /@ logdata, PageBreakBelow -> True, Visible -> False]];

  Return[{modelout, AddInvariants, tsout}]
]

GetHam[models_, vars_, svars_, OptionsPattern[{"subset"->{}}]] := Module[{out, subvars, fixsub, fullvars, polynomials, ham, s},
  fullvars = Append[vars, svars];
  ham = Total[(Times @@ #) & /@ models, 2];
  polynomials = Flatten[If[Head[#] === Plus, Level[#, 1], {#}] &/@ OptionValue["subset"]];
  subvars = If[polynomials === {}, {Flatten[fullvars]}, Variables[#] &/@ polynomials];
  out = Sum[fixsub = Thread[Complement[Flatten@fullvars, s] -> 0]; 
            ham /. fixsub, {s, subvars}];
  Return[out]
]

Expression2ModelData[ham_, vars_, svars_] := Module[{data, a, b},
  data = If[Head[ham] === Plus,
            Cases[ham, Times[a_, b_] /; And[AtomQ[a], ! AtomQ[b]] -> {a, b}],
            Cases[{ham}, Times[a_, b_] /; And[AtomQ[a], ! AtomQ[b]] -> {a, b}]];
  Return[SortModelData[DeleteCases[{data\[Transpose]}, {}], vars, svars]]
]

GetActivePolynomial[modeldata_, vars_, svars_, mim_] := Module[{fullvars, PolynomialList, FullPolynomialList, i, out, models},
  fullvars = Append[vars, svars];
  models = Transpose[Flatten[#] &/@ (Expression2ModelData[GetHam[modeldata, vars, svars, "subset" -> {Chop[mim]*Flatten[fullvars]}], vars, svars]\[Transpose])];
  PolynomialList = If[Head[#2] === Plus, Level[#2, 1], {#2}] & @@@ models;
  FullPolynomialList = If[Head[#] === Plus, Level[#, 1], {#}] & /@ Flatten[modeldata\[Transpose][[2]]];
  out = Flatten[Table[If[IntersectingQ[FullPolynomialList[[i]], #], {i, Flatten[modeldata\[Transpose][[2]]][[i]]}, ## &[]], {i, Length@FullPolynomialList}] & /@ PolynomialList, 1];
  Return[GatherBy[SortBy[out, First], NumPolynomialVar[#[[2]]] &]]
]

AdjustHam[models_, ham_, vars_, svars_] := Module[{data, a, b, i, srt, modelsnew},
  modelsnew = models;
  data = Expression2ModelData[ham, vars, svars];

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

CheckMinimums[models_, seeds_, amp_, vars_, svars_, OptionsPattern[{"round" -> 10^-4, "nsw"->10000}]] := Module[{sol, fullvars, init, s, minimums, simplifytab, a, tabdata},
  simplifytab = DeleteDuplicates@Flatten[Table[{-1*a -> Style["-1", Bold, Blue], 1*a -> Style["1", Bold, Red], 0 -> "-"}, {a, amp}]];
  fullvars = Flatten[{vars, svars}];
  minimums = Table[init = {fullvars, Round[a*s,1], -(Abs[amp]+1)*ConstantArray[1,Length[s]], (Abs[amp]+1)*ConstantArray[1,Length[s]]}\[Transpose];
                   sol = Quiet@FindMinimum[GetHam[models, vars, svars], init, MaxIterations -> OptionValue["nsw"]];
                   {Prepend[# & /@ (init\[Transpose][[2]]), " "] /. simplifytab,
                   Prepend[# & /@ Chop@N@Round[fullvars /. sol[[2]], OptionValue["round"]], Chop[sol[[1]]]]}, {s, seeds}, {a, amp}];
  tabdata = Reverse@SortBy[DeleteDuplicates[Flatten[minimums, 1], (Chop[Norm[#1[[2,1]]-#2[[2,1]]], 10^-6]==0)&], #[[2,1]]&];
  Print[Grid[Prepend[Flatten[tabdata, 1], Prepend[fullvars, "ene"]], 
             Background -> {Flatten[{LightBlue, Table[ConstantArray[If[OddQ[v], Yellow, Pink], Length[Append[vars, svars][[v]]]], {v, Length@Append[vars, svars]}]}], Prepend[Table[If[OddQ[ibg], White, None], {ibg, 2 Length[minimums]}], Cyan]}, ItemSize -> Full]];
]

CheckFitting[invariants_, model_, ts_, vars_, svars_, OptionsPattern[{"imagesize" -> 200, "col" -> 4, "plot" -> True, "plotrange" -> All, "alat" -> 1}]] := Module[{inv, fullvars, character, data2d, enemodel, s, m, i, t, plt},
  fullvars = Append[vars, svars];
  plt = Flatten@Table[inv = If[MatchQ[Head@m, Plus], m[[1]], m];
                      character = If[MemberQ[Variables@inv, #], 1, 0] & /@ Flatten[fullvars];
                      Table[If[Sign@Abs@Join[First[ts[[t]]][[5]], First[ts[[t]]][[6]]] == character,
                               data2d = Flatten[{#5/OptionValue["alat"], #6, #4}] & @@@ (ts[[t]]);
                               enemodel = (GetHam[model, vars, svars] /. Thread[Flatten[fullvars] -> #[[1 ;; -2]]]) & /@ data2d;
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
                                          Plot[GetHam[model, vars, svars] /. Thread[Flatten[fullvars] -> s Normalize[Last[data2d][[1 ;; -2]]]], {s, 0, 1.5 Norm[Last[data2d][[1 ;; -2]]]}, 
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

SortModelData[modeldata_, vars_, svars_] := Module[{fullvars, sortedinv, out, model},
  fullvars = Append[vars, svars];
  out = Table[sortedinv = Flatten[SortInvariants[model[[2]], Flatten@fullvars]];
              SortBy[model\[Transpose], Position[sortedinv, #[[2]]] &]\[Transpose], {model, modeldata}];
  Return[out]
]

JoinModels[m0_, m1_] := Module[{data0, data1, tmp, out},
  data0 = Flatten[#] & /@ (m0\[Transpose]);
  data1 = Flatten[#] & /@ (m1\[Transpose]);
  tmp = If[MemberQ[Flatten[data0], #2], ## &[], {#1, #2}] & @@@ (data1\[Transpose]);
  DeleteCases[Append[m0, tmp\[Transpose]], {}]
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

(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
