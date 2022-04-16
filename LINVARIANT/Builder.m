BeginPackage["LINVARIANT`Builder`", {"LINVARIANT`Structure`", "LINVARIANT`Vasp`", "LINVARIANT`INVARIANT`", "LINVARIANT`MathematicaPlus`",  "LINVARIANT`Parser`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
GetPhononBasis            ::usage "GetPhononBasis[pos0, Fij]"
ShowPhononBasis           ::usage "ShowPhononBasis[pos0, phonon]"
PononDecomposition        ::usage "PononDecomposition[pos0, pos, phonon]"
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
ImportDMDPath             ::usage "ImportDMDPath[file, pos0, Basis, Basislabel]"
SamplePhaseSpace          ::usage "SamplePhaseSpace[dir0, s0, spgmat, Basis, potim]"
CollectSeeds              ::usage "CollectSeeds[dir0, nphase, Basis, OpMat]"
CollectPhases             ::usage "CollectPhases[dir0, n, vars, Basis]"
ShowSeedPhases            ::usage "ShowSeedPhases[phases, seeds, vars]"
SamplePolynomialSpace     ::usage "SamplePolynomialSpace[polyphase, dir, spgmat0, Basis, npt, smearing, mim0]"
PlanSampling              ::usage "PlanSampling[dir0_, phases, invariants, vars, Basis, spgmat, npt, smearing]"
SamplingByDMDPath         ::usage "SamplingByDMDPath[dir0, dmdpath, invariants, vars, svars, Basis, boundingfactor]"
SampleAround              ::usage "SampleAround[pos0, dir, fname, spgmat0, Basis, npt, potim]"
ImposePhase               ::usage "ImposePhase[pos0, dir, fname, phase, Basis]"
RoundPhase                ::usage "RoundPhase[phase]"
LoadTrainingset           ::usage "LoadTrainingset[pos0, dir, Basis, spgmat, numts, Ecubic]"

(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
GetPhononBasis[pos0_, Fij_, OptionsPattern[{"table" -> True, "fontsize" -> 12}]] := Module[{sol, PhononByIR, PhononRotated, ind, frequency, IR, i, j, iv, eigenvectors, rot, NumAtom, tab, simplifytab},
  NumAtom = Length[pos0[[2]]];
  simplifytab = {0 -> "-"};
  sol = Chop[Eigensystem[ArrayFlatten[Fij]]]\[Transpose];
  PhononByIR = Gather[MapIndexed[Join[#2, #1] &, sol], Chop[#1[[2]] - #2[[2]]] == 0 &];
  PhononRotated = Table[{ind, frequency, IR} = PhononByIR[[i]]\[Transpose];
                        eigenvectors = DeleteCases[Chop@Orthogonalize@Flatten[Table[rot = PseudoInverse[IR\[Transpose][[3 (j - 1) + 1 ;; 3 (j - 1) + 3]]\[Transpose]]; Normalize[#] & /@ (rot . IR), {j, Length[pos0[[2]]]}], 1], v_ /; Norm[v] == 0]; 
                        Table[{ind[[iv]], frequency[[iv]], eigenvectors[[iv]]}, {iv, Length[ind]}], {i, Length[PhononByIR]}];
  tab = Grid[Join[{Join[{"#","frequency (THz)"}, Flatten[Table[Subscript[#, \[Alpha]], {\[Alpha], {"x", "y", "z"}}] & /@ (pos0[[2]]\[Transpose][[2]])]]}, Flatten[#] /. simplifytab & /@ Flatten[PhononRotated, 1]], 
             Alignment -> Center, 
             Frame -> True, 
             Background -> {Flatten[{Blue, Blue, Table[color = If[OddQ[i], White, None]; 
                                                       ConstantArray[color, 3], {i, NumAtom}]}], 
                            Flatten[{LightBlue, Table[color = If[OddQ[i], Yellow, Pink]; 
                                                      ConstantArray[color, Length[PhononRotated[[i]]]], {i, Length@PhononRotated}]}]}, 
             ItemSize -> Full, 
             ItemStyle -> Directive[FontSize -> OptionValue["fontsize"], Black]];
  If[OptionValue["table"], Print[tab]];
  Return[PhononRotated]
]

ShowPhononBasis[pos0_, phonon_, OptionsPattern[{"fontsize" -> 12}]] := Module[{sol, PhononByIR, PhononRotated, ind, frequency, IR, i, j, iv, eigenvectors, rot, NumAtom, tab, simplifytab},
   NumAtom = Length[pos0[[2]]];
   simplifytab = {0 -> "-"};
   tab = Grid[Join[{Join[{"#","frequency (THz)"}, Flatten[Table[Subscript[#, \[Alpha]], {\[Alpha], {"x", "y", "z"}}] & /@ (pos0[[2]]\[Transpose][[2]])]]}, Flatten[#] /. simplifytab & /@ Flatten[phonon, 1]],
              Alignment -> Center,
              Frame -> True,
              Background -> {Flatten[{Blue, Blue, Table[color = If[OddQ[i], White, None];
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

PononDecomposition[pos0_, pos_, phonon_] := Module[{i},
  Grid[Prepend[ISODISTORT[pos0[[1]], pos0[[2]], pos[[2]], Flatten[phonon, 1]\[Transpose][[3]]\[Transpose], Flatten[phonon, 1]\[Transpose][[2]]], {"#", "frequency", "norm", "amplitude"}], 
       Background -> {None, 
                      Flatten[{LightBlue, Table[color = If[OddQ[i], Yellow, Pink]; 
                              ConstantArray[color, Length[phonon[[i]]]], {i, Length@phonon}]}]}, 
       ItemSize -> Full]
]

PhononLinearRecombination[PhononBasis_, rules_] := Module[{rule, rot, m, s, pos, id, phonon, i, ind, IR, frequency, eigenvectors},
  ind = Sort[Flatten[#\[Transpose][[1;;2]]]] & /@ rules;
  Table[IR = PhononBasis[[i]]\[Transpose][[1]];
        frequency = PhononBasis[[i]]\[Transpose][[2]];
        pos = Position[ind, IR];
        rot = If[pos != {},
                 rule = Extract[rules, First[pos]];
                 1/Sqrt[2] Flatten[Table[Which[Length[id] == 1, m = ConstantArray[0, Length[IR]]; m[[id[[1]]-Min[IR]+1]] = Sqrt[2]; {m},
                                               Length[id] == 3, RotateLeft[Table[m = ConstantArray[0, Length[IR]];
                                                                                 m[[id[[1]]-Min[IR]+1]] = 1; m[[id[[2]]-Min[IR]+1]] = (-1)^(s - 1);
                                                                                 m, {s, 2}], id[[3]]]], {id, rule}], 1], 
                 IdentityMatrix[Length[IR]]];
        eigenvectors = rot . (PhononBasis[[i]]\[Transpose][[3]]);
        {IR, frequency, Chop@eigenvectors}\[Transpose], {i, Length[PhononBasis]}]
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

FetchSeeds[inv_, vars_, n_, OptionsPattern[{"zero" -> False}]] := Module[{character, ndigits, shortseed, list, map, seeds0, order}, 
  character = First@InvariantCharacter[inv, vars];
  list = If[OptionValue["zero"], Range[n,-n,-Sign[n]], DeleteCases[Range[n,-n,-Sign[n]],0]];
  order = Length[list];
  ndigits = Count[character, 1];
  shortseed = Tuples[list, ndigits]\[Transpose];
  map = Flatten[Position[character, 1]];
  seeds0 = Table[If[i == 0, ConstantArray[0, order^ndigits], ## &[]], {i, character}];
  Do[seeds0 = Insert[seeds0, shortseed[[i]], map[[i]]], {i, Length@map}];
  Return[seeds0\[Transpose]]
]

SeedsByInvariants[invariants_, vars_, OptionsPattern[{"table"->True, "exseeds"->None}]] := Module[{id, inv, sub, i, seeds, s, dict, sampleseeds, tab, tabdata, color, fullvars, simplifytab, varlabels, exseeds, invdict, models, t},
  fullvars = Flatten[vars];
  invdict = GroupBy[Flatten[invariants], StrainInvQ[#, vars[[-1]][[1, 1]]] &];
  models = Flatten[Table[SortInvariants[t, fullvars], {t, {invdict[False], invdict[True]}}]];

  varlabels = MapIndexed["(" <> ToString[First@#2] <> ")" &, fullvars];
  simplifytab = {-1 -> Style["-1", Bold, Blue], 1 -> Style["1", Bold, Red], 0 -> "-"};

  exseeds = If[OptionValue["exseeds"]===None, 
               {{0,"f[0]",ConstantArray[0,Length[fullvars]]}},
               {0,"f["<>StringJoin[Riffle[ToString[#, StandardForm] & /@ DeleteCases[fullvars Abs[#], 0], ","]] <> "]", #} &/@ OptionValue["exseeds"]];

  seeds = Flatten[Table[inv = First[models[[id]]];
                        s = FetchSeeds[models[[id]], fullvars, 1, "zero"->False];
                        dict = Values@GroupBy[Table[sub = Thread[fullvars -> s[[i]]]; {inv /. sub, s[[i]]}, {i, Length[s]}], First];
                        {id, #1 inv, #2} & @@@ First[dict\[Transpose]], {id, Length@models}], 1];
  sampleseeds = Join[exseeds, DeleteDuplicates[seeds, #1[[3]]==#2[[3]]&]];
  tabdata = {#1, #2, #3 /. simplifytab} & @@@ sampleseeds;
  tab = Grid[Prepend[Prepend[MapIndexed[Flatten[{#2, #1}] &, tabdata], Flatten[{{"", "", ""}, fullvars}]], Flatten[{{"", "", ""}, varlabels}]],
             Background -> {Flatten[{LightBlue, White,LightBlue, 
                            Table[color = If[OddQ[i], Yellow, Pink];
                                  ConstantArray[color, Length[vars[[i]]]], {i, Length@vars}]}], None}];
  If[OptionValue["table"], Print[tab]];
  Return[sampleseeds]
]

ImportDMDPath[file_, pos0_, Basis_, vars_, OptionsPattern[{"fontsize"->12, "table"->True, "plot"->True, "interpolate"->0}]] := Module[{vasprun, lattices, sites, ene, coordinates, labels, color, i, path, epsplot, scfplot, roseplot, out, eta, fullvars, rescale},
  fullvars = Flatten[vars];
  rescale = Flatten[{ConstantArray[1, Length[fullvars]-6], ConstantArray[20, 6]}];
  vasprun = Import[file];
  labels = MapIndexed["("<>ToString[First@#2]<>") "<>ToString[#1,StandardForm]&, fullvars];
  sites = ParseXMLData[ParseXML[ParseXML[vasprun, "structure", {}], "varray", {"name", "positions"}], "v"];
  lattices = ParseXMLData[ParseXML[ParseXML[vasprun, "structure", {}], "varray", {"name", "basis"}], "v"];
  eta = Chop[eij2eta[GetStrainTensor[pos0[[1]], #, "iso" -> False]], 10^-6]&/@lattices;
  ene = ParseFortranNumber@Flatten[ParseXML[Last[ParseXML[#, "energy", {}]], "i", {"name", "e_fr_energy"}] & /@ ParseXML[vasprun, "calculation", {}]];
  coordinates = Table[Flatten[{ISODISTORT[pos0[[1]], pos0[[2]], {sites[[i]], pos0[[2]]\[Transpose][[2]]}\[Transpose], Basis, Range[Length[Basis\[Transpose]]], "round"->10.0^-12]\[Transpose][[4]], eta[[i]]}], {i, Length[sites]}]\[Transpose];
  path = Grid[Prepend[Prepend[coordinates, ene]\[Transpose], Prepend[labels, "energy"]],
              Background -> {Prepend[Flatten[Table[color = If[OddQ[i], Yellow, Pink];
                             ConstantArray[color, Length[vars[[i]]]], {i, Length@vars}]], LightBlue],
                             None},
              ItemStyle -> Directive[FontSize -> OptionValue["fontsize"]],
              ItemSize->Full
             ]; 
  epsplot = Grid@Partition[Table[If[Norm[coordinates[[i]]] != 0, ListPlot[{coordinates[[i]], ene}\[Transpose], Frame -> True, GridLines -> Automatic, PlotRange -> All, PlotLabel -> labels[[i]], ImageSize -> 300], ## &[]], {i, Length[fullvars]}], UpTo[4]];
  scfplot = Grid[{{ListLinePlot[coordinates*rescale, PlotRange -> All, Frame -> True, PlotMarkers -> Automatic, 
                                             GridLines -> Automatic, ImageSize -> 300, ImagePadding -> {{50, 1}, {20, 1}}, 
                                             PlotLegends -> LineLegend[labels, LegendLayout -> {"Column", 3}]], 
                   ListLinePlot[ene, PlotRange -> All, Frame -> True, PlotMarkers -> Automatic, GridLines -> Automatic, 
                                     ImageSize -> 300, ImagePadding -> {{50, 1}, {20, 1}}]}}];
  roseplot = Grid@Partition[Rasterize[ListPlot3D[{coordinates[[#1]], coordinates[[#2]], ene}\[Transpose], AxesLabel -> {labels[[#1]], labels[[#2]], ""}, PlotRange -> All, ColorFunction -> "Rainbow", InterpolationOrder -> OptionValue["interpolate"], MeshFunctions -> {#3 &}, Mesh -> Length[ene], ViewPoint -> {0, 0, 100}, ImageSize -> 300], ImageSize->300] & @@@ Subsets[Table[If[Norm[coordinates[[i]]] == 0, ## &[], i], {i, Length[coordinates]}], {2}], UpTo[4]];
  If[OptionValue["table"], Print[path]];
  If[OptionValue["plot"], Print[Grid[{{epsplot},{roseplot},{scfplot}}]]];
  out = {Append[coordinates, ene], {lattices, sites}\[Transpose]};
  Return[out]
]

SamplingByDMDPath[dir0_, dmdpath_, invariants_, vars_, svars_, Basis_, boundingfactor_, OptionsPattern[{"fast" -> True}]] := Module[{tsdir, bddir, dist, p, icount, icount0, i, j, numdmd, idmd, jdmd, id, s, im, m, models, character, initcharacter, coordinates, trajectory, shiftphase, shiftmode, shifteij, strainlatt, pos, pos0, pos1, fullvars, NumBasis, invdict, seeds0, disp, info, phase0, phase1, phase, latt, sites, ij, numsamples}, 
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
   
     dist = Table[phase = First[First[p]\[Transpose]][[1 ;; -2]];
                  initcharacter = RoundPhase[phase]; 
                  Abs[Normalize[v].Normalize[initcharacter]], {v, character}, {p, dmdpath}];
     Do[{i, j} = ij;
        {coordinates, trajectory} = dmdpath[[j]];
        numdmd = Length[trajectory];
 
        phase0 = First[coordinates\[Transpose]][[1;;-2]];
        phase1 = Last[coordinates\[Transpose]][[1;;-2]];

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

SamplePhaseSpace[dir0_, s0_, spgmat_, Basis_, potim_, OptionsPattern[{"seeds" -> {}}]] := Module[{BasisDim, NumBasis, modesets, seeds, modes, m, i, pos, pos0, eij, strainlatt},
  If[! DirectoryQ[dir0<>"/seeds"], CreateDirectory[dir0<>"/seeds"]; Run["chmod o-w " <> dir0<>"/seeds"]];
  pos0 = ImportPOSCAR[dir0 <> "/" <> s0];
  {BasisDim, NumBasis} = Dimensions[Basis];

  If[OptionValue["seeds"] == {}, 
     modesets = SortBy[Tuples[{-1, 0, 1}, {NumBasis}], Norm[N@#] &];
     seeds = Round[DeleteDuplicates[Table[Sort[# . modesets[[i]] & /@ spgmat], {i, Length[modesets]}]]\[Transpose][[1]]],
     seeds = OptionValue["seeds"]
  ];

  Print["Number of seeds: " <> ToString[Length[seeds]]];
  Do[eij=potim*Chop@eta2eij[seeds[[i]][[-6;;]]]/40;
     strainlatt=(eij+IdentityMatrix[3]).(pos0[[1]]);
     modes = Table[{m, "modes", 1, potim seeds[[i, m]]}, {m, NumBasis}]; 
     pos = {strainlatt, ImposeMode[pos0[[1]], pos0[[2]], Basis, modes, 1.0]}; 
     ExportPOSCAR[dir0<>"/seeds", "phase." <> ToString[i] <> ".vasp", pos];
     ExportOPTCELL[dir0<>"/seeds/OPTCELL." <> ToString[i] <> ".inp", Abs[seeds[[i]][[-6;;]]]], {i, Length@seeds}]
]

CollectSeeds[dir0_, nphase_, Basis_, OpMat_, OptionsPattern[{"round"->10^-6, "file"->"POSCAR", "all"->False}]] := Module[{pos, mode, eta, pos0, NumDim, NumBasis, Basislabel, i, file, phases, m},
  {NumDim, NumBasis} = Dimensions[Basis];
  pos0 = ImportPOSCAR[dir0 <> "/POSCAR0"];
  Basislabel = Table["x" <> ToString[i], {i, NumBasis}];
  phases = Table[file = StringReplace[OptionValue["file"], "xxx"->ToString[i]];
                 pos = ImportPOSCAR[dir0 <> "/seeds/phase" <> ToString[i] <> "/" <> file];
                 mode = Chop[ISODISTORT[pos0[[1]], pos0[[2]], pos[[2]], Basis, Basislabel, "round" -> OptionValue["round"]]\[Transpose][[4]], OptionValue["round"]];
                 eta = Chop[eij2eta[GetStrainTensor[pos0[[1]], pos[[1]], "iso" -> False]]];
                 {i, mode, eta, pos, Table[Sign[#]&/@Chop[m.mode,OptionValue["round"]], {m, OpMat}]}, {i, Range[nphase]}];
  Return@If[OptionValue["all"], phases\[Transpose][[1 ;; 4]]\[Transpose], 
                                DeleteDuplicates[phases, Complement[#1[[5]],#2[[5]]] == {} &]\[Transpose][[1 ;; 4]]\[Transpose]]
]

CollectPhases[dir0_, n_, vars_, Basis_, OptionsPattern[{"round" -> 10^-6, "fontsize" -> 12, "table" -> True}]] := Module[{NumDim, NumBasis, pos0, Basislabel, i, mpdata, pos, sites, eta, phase, simplifytab, fullvars, tab},
  fullvars = Flatten[vars];
  {NumDim, NumBasis} = Dimensions[Basis];
  pos0 = ImportPOSCAR[dir0 <> "/POSCAR0"];
  Basislabel = Table["x" <> ToString[i], {i, NumBasis}];
  simplifytab = {-1 -> Style["-1", Bold, Blue], 1 -> Style["1", Bold, Red], 0 -> "-"};
  mpdata = Table[pos = ImportPOSCAR[dir0 <> "/phase" <> ToString[i] <> "/POSCAR"];
                 sites = {PosMatchTo[pos0[[1]], pos0[[2]]\[Transpose][[1]], pos[[2]]\[Transpose][[1]]][[2]], pos0[[2]]\[Transpose][[2]]}\[Transpose];
                 eta = eij2eta[GetStrainTensor[pos0[[1]], pos[[1]], "iso" -> False]];
                 phase = ISODISTORT[pos0[[1]], pos0[[2]], sites, Basis, Basislabel, "round" -> OptionValue["round"]]\[Transpose][[4]];
                 {i, phase, eta}, {i, 1, n}];
  tab = Grid[Prepend[Flatten[{Join[{#1}, #2, #3], Join[{" "}, Sign[#2], Sign[#3]] /. simplifytab} & @@@ mpdata, 1], Prepend[fullvars, "#"]],
             Background -> {Prepend[Flatten[Table[color = If[OddQ[i], Yellow, Pink]; ConstantArray[color, Length[vars[[i]]]], {i, Length@vars}]], LightBlue], Prepend[Flatten[Table[If[EvenQ[i], White, None], {i, 2 n}]], White]}, 
             Dividers -> {False, Table[If[EvenQ[i], i -> Black, ## &[]], {i, 2 n}]}, 
             ItemStyle -> Directive[FontSize -> OptionValue["fontsize"]], 
             ItemSize -> Full];
  If[OptionValue["table"], Print[tab]];
  Return[mpdata]
]

ShowSeedPhases[phases_, seeds_, vars_, OptionsPattern[{"fontsize" -> 12}]] := Module[{i, color,labels, fullvars, simplifytab},
  fullvars = Flatten[vars];
  labels = MapIndexed["("<>ToString[First@#2]<>") "<>ToString[#1,StandardForm]&, fullvars];
  simplifytab = {-1 -> Style["-1", Bold, Blue], 1 -> Style["1", Bold, Red], 0 -> "-"};
  Grid[Prepend[Flatten[{Flatten[#], Prepend[seeds[[#[[1]],3]] /. simplifytab, "("<>ToString@seeds[[#[[1]],1]]<>") "<>ToString[Times@@DeleteCases[seeds[[#[[1]],3]],0] seeds[[#[[1]],2]], StandardForm]]} & /@ (phases\[Transpose][[1 ;; 3]]\[Transpose]), 1], Flatten[Join[{"#"}, labels]]], 
       Background -> {Prepend[Flatten[Table[color = If[OddQ[i], Yellow, Pink]; 
                              ConstantArray[color, Length[vars[[i]]]], {i, Length@vars}]], LightBlue],
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

LoadTrainingset[pos0_, dir_, invariants_, vars_, svars_, Basis_, spgmat_, OptionsPattern[{"round" -> 10^-6}]] := Module[{i, ts, t, s, BasisDim, NumBasis, ene, pos, mode0, vasprun, mode, data, eta, out, ws, seeds, plan, p, Ecubic},
  {BasisDim, NumBasis} = Dimensions[Basis];
  data = Import[dir <> "/" <> "ene.dat", "Table"];
  Ecubic = data[[1,4]];
  out=Table[ts = data[[i,1]]; t = data[[i, 2]]; s = data[[i, 3]]; ene = data[[i, 4]];
          pos = ImportPOSCAR[dir <> "/"<>ts<>"-"<> ToString[t] <> "/sample" <> ToString[s] <> "/POSCAR"];
          eta = Chop[eij2eta[GetStrainTensor[pos0[[1]], pos[[1]], "iso" -> False]]];
          If[s == 1, mode0 = ISODISTORT[pos0[[1]], pos0[[2]], pos[[2]], Basis, Range[NumBasis], "round" -> OptionValue["round"]]\[Transpose][[4]]];
          mode = ISODISTORT[pos0[[1]], pos0[[2]], pos[[2]], Basis, Range[NumBasis], "round" -> OptionValue["round"]]\[Transpose][[4]];
          {i, ts, t, s, ene - Ecubic, mode0, mode - mode0, eta, Norm[mode], # . mode & /@ spgmat}, {i, Length[data]}]\[Transpose][[1 ;; 8]]\[Transpose];

  ws=DeleteDuplicates[If[#4 == 1 && Norm[#8] == 0., {#5, #6, {0.00001, 1}}, ## &[]] & @@@ out, Chop[Norm[#1[[1]] - #2[[1]]]] == 0. &];
  seeds=If[ws === None, {vars /. {Subscript[x_, i_] -> 0}}, RoundPhase[#] & /@ (ws\[Transpose][[2]])];
  plan = PlanFitting[invariants, vars, svars, spgmat, out, seeds];
  Print["Checking training set."];
  Do[If[p[[2]] == {}, Print["Warning! Training set is not complete:"]; Print[p[[1]]]], {p, plan}];

  Return[out]
]

PlanFitting[invariants_, vars_, svars_, spgmat_, ts_, seeds_] := Module[{p, m, invdict, t, i, models, character, modes, phases, plan, fullvars},
  fullvars = Flatten[{vars, svars}];
  invdict = GroupBy[Flatten[invariants], StrainInvQ[#, svars[[1, 1]]] &];
  models = Flatten[Table[SortInvariants[t, fullvars], {t, {invdict[False], invdict[True]}}], 1];
  Print["Number of invariants :" <> ToString[Length[Flatten@models]]];
  Print["Groups of models :" <> ToString[Length[models]]];
  
  plan = Table[EtaQ = StrainInvQ[m[[1]], svars[[1, 1]]];
               character = First@InvariantCharacter[m[[1]], Flatten@vars];
               {m, Table[If[Xor[Chop[Norm[t[[8]]]] == 0., EtaQ],
                            modes = DeleteDuplicates[RoundPhase[Normal[#].(t[[6]] + t[[7]])] & /@ spgmat];
                            phases = DeleteDuplicates[RoundPhase[Normal[#].t[[6]]] & /@ spgmat];
                            If[MemberQ[modes, character], If[MemberQ[seeds, character], If[MemberQ[phases, character], t, ## &[]], t], ## &[]], ##&[]], {t, ts}]}, {m, models}];
  Return[plan]
]

LinvariantFit[invariants_, vars_, svars_, order_, ts_, spgmat_, ws_ : None, OptionsPattern[{"orderby" -> "AICc", "round" -> 10^-6.0, "offset" -> None, "alat" -> 1.0, "FontSize" -> 12}]] := Module[{alat, NumBasis, models, modelorders, OffsetFunc, fits, TopModels, BoundRadius, tsdata, seeds, plan, measure, x, i, m, p, t, mf, gc, log, out, phaseweight, minimumtest, sol, modelfitted, modeldata, tsinfo, fullvars, ibg, checkQ},
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


(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
