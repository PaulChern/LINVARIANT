BeginPackage["LINVARIANT`SpinHamiltonian`", {"LINVARIANT`Structure`", "LINVARIANT`GroupTheory`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
GetJi0Bonds          ::usage "GetJi0Bonds[spg0, tij, latt, AllSites]"
Jij2JKD              ::usage "Jij2JKD[J]"
PlotSpinField        ::usage "PlotSpinField[latt, dim]"
SpecificHeat         ::usage "SpecificHeat[T, EOnSiteList]"
MagnonHk             ::usage "SpinHk[Ji0, k]"
MagnonBandsPlot      ::usage "MagnonBandsPlot[Ji0, klist, kintv]"
(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
GetJi0Bonds[spg0_, tij_, pos_] := Module[{latt, AllSites, c1, c2, c3, sol, site, AllBonds, ReducedBonds},
  {latt, AllSites} = pos;
  AllBonds = Flatten[Table[{(xyzStr2M4[xyz].Join[AllSites[[#1[[1]], 1]], {1}])[[1 ;; 3]], (xyzStr2M4[xyz].Join[AllSites[[#1[[2]], 1]] + #2, {1}])[[1 ;; 3]], (latt.Inverse[xyzStr2M4[xyz]][[1 ;; 3, 1 ;; 3]].Inverse[latt])\[Transpose].#3.(latt.Inverse[xyzStr2M4[xyz]][[1 ;; 3, 1 ;; 3]].Inverse[latt])}, {xyz, Keys[spg0]}] & @@@ tij, 1];
  ReducedBonds = Flatten[Table[sol = Rationalize@Chop@First@Values@Solve[site[[2]] + IdentityMatrix[3].{c1, c2, c3} == site[[1]]]; If[AllTrue[sol, IntegerQ], {#1 - sol, #2 - sol, #3}, ## &[]], {site, Tuples[{{#1, #2}, AllSites\[Transpose][[1]]}]}] & @@@ AllBonds, 1];
  ReducedBonds = DeleteDuplicates[If[Rationalize[Chop[#1 - Mod[#1, 1]]] == {0, 0, 0}, {#1, #2, #3}, {#2, #1, #3\[Transpose]}] &@@@ ReducedBonds, Chop[#1[[1]] - #2[[1]]] == {0, 0, 0} && Chop[#1[[2]] - #2[[2]]] == {0,  0, 0} &];
  ReducedBonds = {pos2index[AllSites, {#1, Mod[#2, 1]}], Rationalize[#2 - Mod[#2, 1]], 1, #3} & @@@ ReducedBonds;
  Return[ReducedBonds]
]

Jij2JKD[J_] := Module[{JJ, Dij, Kij},
  Dij = Chop@Expand[1/2 (J - J\[Transpose])];
  Kij = Expand[Chop[Expand[1/2 (J + J\[Transpose])] - 1/3 Tr[Expand[1/2 (J + J\[Transpose])]] IdentityMatrix[3]]];
  JJ = Expand[1/3 Tr[Expand[1/2 (J + J\[Transpose])]] IdentityMatrix[3]];
  Return[Rationalize@{JJ, Kij, {Dij[[2, 3]], Dij[[3, 1]], Dij[[1, 2]]}}]
]

PlotSpinField[latt_, dim_, OptionsPattern[{"vec" -> {}, "ImageSize" -> 400}]] := Module[{Nx, Ny, Nz, ExtVec}, 
  {Nx, Ny, Nz} = dim;
  ExtVec = If[OptionValue["vec"] != {}, {ColorData["TemperatureMap"][Rescale[(latt.#2)[[3]], {-1, 1}]], Arrowheads[0.015], Arrow[Tube[{latt.#1, latt.(#1 + #2)}, 0.3]]} & @@@ OptionValue["vec"], {{}}];
  Graphics3D[Join[ExtVec, 
                  {{Gray, PointSize[0.005], Point[Flatten[Table[latt.{ix, iy, iz}, {ix, Nx}, {iy, Ny}, {iz, Nz}], 2]]}, 
                  {Red, PointSize[0.005], Point[Flatten[Table[latt.{ix + 1/3, iy + 2/3, iz}, {ix, Nx}, {iy, Ny}, {iz, Nz}], 2]]}, 
                  {Pink, PointSize[0.005], Point[Flatten[Table[latt.{ix + 2/3, iy + 1/3, iz}, {ix, Nx}, {iy, Ny}, {iz, Nz}], 2]]}}], 
                  ImageSize -> OptionValue["ImageSize"], ViewPoint -> {0, 0, 100}
            ]
]

MCStep[field_, T_, initdamp_] := Module[{AcceptRatio, diceX, diceY, diceZ, ix, iy, iz, imc, idum, jdum, rdum, accepted, damp, dspin, FlippedSpin, FlippedField, DicedField, \[CapitalDelta]E, AcceptanceProbability, KB = 8.617 10^-5},
  imc = 0;
  accepted = 0;
  AcceptRatio = 0.4;
  damp = initdamp;

  Do[FlippedField = field;
   FlippedSpin = damp RandomReal[{-1, 1}, 3] + field[[is, ix, iy, iz]];
   dspin = FlippedSpin/Norm[FlippedSpin] - field[[is, ix, iy, iz]];
   FlippedField[[is, ix, iy, iz]] = dspin + field[[is, ix, iy, iz]];
   \[CapitalDelta]E = DeltaE[{ix, iy, iz}, {Lx, Ly, Lz}, field, dspin, is];
   AcceptanceProbability = Min[1, Quiet@Exp[-(\[CapitalDelta]E/(KB T)) ]];
   FlippedField = If[RandomReal[] <= AcceptanceProbability, accepted += 1; FlippedField, field];imc += 1, 
   {ix, Lx}, {iy, Ly}, {iz, Lz}, {is, 2}];

  rdum = accepted/(2 Lx Ly Lz);

  If[rdum > AcceptRatio, damp = 1.01 damp, damp/1.01];
  damp = Min[damp, 10.0];
  
  Do[jdum = Mod[idum, 3, 1];
   diceX = RandomInteger[{1, Lx}];
   diceY = RandomInteger[{1, Ly}];
   diceZ = RandomInteger[{1, Lz}];
   DicedField = FlippedField;
   DicedField[[is, diceX, diceY, diceZ]][[jdum]] *= -1;
   dspin = DicedField[[is, diceX, diceY, diceZ]] - FlippedField[[is, diceX, diceY, diceZ]];

   \[CapitalDelta]E = DeltaE[{diceX, diceY, diceZ}, {Lx, Ly, Lz}, FlippedField, dspin, is];
   
   AcceptanceProbability = Min[1, Quiet@Exp[-(\[CapitalDelta]E/(KB T)) ]];
   FlippedField = If[RandomReal[] <= AcceptanceProbability, DicedField, FlippedField];
   , {idum, 1, 3 3 IntegerPart[Sqrt[2 Lx Ly Lz]]}, {is, 2}];
  Return[{damp, FlippedField}]
]

SpecificHeat[T_, EOnSiteList_] := Module[{NumMc, i, KB = 8.617 10^-5, \[Beta]},
  \[Beta] = 1/(KB T);
  NumMc = Length[EOnSiteList];
  (Total[EOnSiteList^2]/NumMc - (Total[EOnSiteList]/NumMc)^2) KB \[Beta]^2
]    

MagnonHk[Ji0_, pos_, k_] := Module[{latt, sites, i, j, HijBlock, H, NumPos},
  {latt, sites} = pos;
  NumPos = Length[sites];
  HijBlock = Merge[{#1 -> #4/#3 Exp[I 2 Pi k.(sites[[#1[[2]], 1]] + #2 - sites[[#1[[1]], 1]])]} & @@@ Ji0, Total];
  H = ArrayFlatten[Table[HijBlock[{i, j}], {i, NumPos}, {j, NumPos}]];
  Return[0.5 (H + H\[ConjugateTranspose])]
]


MagnonBandsPlot[Ji0_, pos_, klist_, kintv_, OptionsPattern[{"range" -> All}]] := Module[{latt, sites, pdata, k, TBHij, sol, kpath, xticks, BandsPlot},
  {latt, sites} = pos;
  {kpath, xticks} = GetKpath[klist, kintv];
  pdata = Table[MagnonHij = MagnonHk[Ji0, pos, k[[2]]];
                sol = Eigenvalues[MagnonHij];
                {k[[1]], ReIm[Sqrt[#]].{1,-1}} & /@ sol, {k, kpath}];
  BandsPlot = ListPlot[pdata\[Transpose],
                       PlotStyle -> Black,
                       Joined -> False,
                       PlotRange -> {All, OptionValue["range"]},
                       AspectRatio -> 1/GoldenRatio,
                       Frame -> True,
                       ImageSize -> Medium,
                       GridLines -> {{xticks\[Transpose][[1]], ConstantArray[Thick, Length[xticks]]}\[Transpose], Automatic},
                       ImageSize -> Medium,
                       FrameTicks -> {{Automatic, None}, {xticks, None}}
                 ];
  Return[BandsPlot]
]

GreenIntersites[spg0_, Ti0_, pos_, \[Epsilon]_, RShift_, Ef_: 0, OptionsPattern[{"kinv" -> 0.1}]] := Module[{Greens1s2ijmn, sites, latt, NumSites, NumOrbitals, GreenBlock, klist, k, i0, j0, Ri, Rj, ii, jj, m, n, \[Sigma]1, \[Sigma]2},
  {latt, sites} = pos;
  NumOrbitals = Length[First[Ti0][[4]]]/2;
  NumSites = Length[sites];
  {ii, jj} = RShift;
  klist = GetBZKList[spg0, OptionValue["kinv"]];
  (* GreenBlock[[\[Sigma]1,\[Sigma]2,i0,j0,m,n]] *)
  Greens1s2ijmn = 1/Length[klist] Sum[GreenBlock = Partition[Partition[k[[2]] Inverse[(\[Epsilon] + Ef) IdentityMatrix[2 NumSites] - TBHk[Ti0, pos, k[[1]]]], {NumOrbitals, NumOrbitals}], {NumSites, NumSites}]; GreenBlock Exp[I 2 Pi k[[1]].(ii - jj)], {k, klist}];
  Return[Chop[Greens1s2ijmn]]
]

DeltaOnSiteSpin[spg0_, Ti0_, pos_, OptionsPattern[{"kinv" -> 0.1}]] := Module[{latt, sites, NumOrbitals, NumSites, klist, Deltaijmn, HBlock},
  {latt, sites} = pos;
  NumOrbitals = Length[First[Ti0][[4]]]/2;
  NumSites = Length[sites];
  klist = GetBZKList[spg0, OptionValue["kinv"]];
  (* Deltaijmn[[i0,j0,m,n]] *)
  Deltaijmn = Sum[HBlock = Partition[Partition[k[[2]] TBHk[Ti0, pos, k[[1]]], {NumOrbitals, NumOrbitals}], {NumSites, NumSites}]; HBlock[[1, 1]] - HBlock[[2, 2]], {k, klist}];
  Return[Deltaijmn]
]

GreenJij[spg0_, Ti0_, pos_, RShift_, Ef_: 0, OptionsPattern[{"kinv" -> 0.1}]] := Module[{latt, sites, NumOrbitals, NumSites, i, j, m, n, p, q, Delta, G, Jij, \[Epsilon]},
  {latt, sites} = pos;
  NumOrbitals = Length[First[Ti0][[4]]]/2;
  NumSites = Length[sites];
  Delta = DeltaOnSiteSpin[spg0, Ti0, pos, "kinv" -> OptionValue["kinv"]];
  Jij = -(1/(2 Pi)) Sum[G = GreenIntersites[spg0, Ti0, pos, \[Epsilon], RShift, Ef, "kinv" -> OptionValue["kinv"]]; Im@Table[Sum[Delta[[i, i, m, n]] G[[2, 2, i, j, n, p]] Delta[[j, j, p, q]] G[[1, 1, j, i, q, m]], {m, NumOrbitals}, {n, NumOrbitals}, {p, NumOrbitals}, {q, NumOrbitals}], {i, NumSites}, {j, NumSites}], {\[Epsilon], -10, Ef}];
  Return[Chop[Jij]]
]

(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
