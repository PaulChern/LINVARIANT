BeginPackage["LINVARIANT`SpinHamiltonian`", {"LINVARIANT`Structure`", "LINVARIANT`GroupTheory`", "LINVARIANT`GreenFunctions`", "LINVARIANT`TBHamiltonian`", "LINVARIANT`MathematicaPlus`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
GetJi0Bonds          ::usage "GetJi0Bonds[spg0, tij, latt, AllSites]"
Jij2JKD              ::usage "Jij2JKD[J]"
PlotSpinField        ::usage "PlotSpinField[latt, dim]"
SpecificHeat         ::usage "SpecificHeat[T, EOnSiteList]"
MagnonHk             ::usage "SpinHk[Ji0, k]"
MagnonBandsPlot      ::usage "MagnonBandsPlot[Ji0, klist, kintv]"
GreenJij             ::usage "GreenJij[Ti0, spg, pos, orb, path, efermi, R, kmesh]"

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
  AllBonds = Flatten[Table[{(xyz2m4[xyz].Join[AllSites[[#1[[1]], 1]], {1}])[[1 ;; 3]], (xyz2m4[xyz].Join[AllSites[[#1[[2]], 1]] + #2, {1}])[[1 ;; 3]], (latt.Inverse[xyz2m4[xyz]][[1 ;; 3, 1 ;; 3]].Inverse[latt])\[Transpose].#3.(latt.Inverse[xyz2m4[xyz]][[1 ;; 3, 1 ;; 3]].Inverse[latt])}, {xyz, Keys[spg0]}] & @@@ tij, 1];
  ReducedBonds = Flatten[Table[sol = Rationalize@Chop@First@Values@Solve[site[[2]] + IdentityMatrix[3].{c1, c2, c3} == site[[1]]]; If[AllTrue[sol, IntegerQ], {#1 - sol, #2 - sol, #3}, ## &[]], {site, Tuples[{{#1, #2}, AllSites\[Transpose][[1]]}]}] & @@@ AllBonds, 1];
  ReducedBonds = DeleteDuplicates[If[Rationalize[Chop[#1 - Mod[#1, 1]]] == {0, 0, 0}, {#1, #2, #3}, {#2, #1, #3\[Transpose]}] &@@@ ReducedBonds, Chop[#1[[1]] - #2[[1]]] == {0, 0, 0} && Chop[#1[[2]] - #2[[2]]] == {0,  0, 0} &];
  ReducedBonds = {pos2index[latt, AllSites, {#1, Mod[#2, 1]}], Rationalize[#2 - Mod[#2, 1]], 1, #3} & @@@ ReducedBonds;
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

GreenJij[Ti0_, spg_, pos_, orb_, path_, efermi_, R_, kmesh_] := Module[{Aij, i, j, J0, a, b},
  Aij=SimpsonIntegrate[VGVG[Ti0, pos, orb, #, efermi, R, kmesh] &, path];
  Print[Aij];
  Jij = If[Length@Ti0==2,
           Table[Im@Aij[[i,j]] IdentityMatrix[3] ,{i, Length@orb}, {j, Length@orb}],
           Table[J0=Im@Total[(Aij[[4,4]] - Aij[[1,1]] - Aij[[2,2]] - Aij[[3,3]])[[i,j]],2];
              Ja=Table[Im@Total[(Aij[[a,b]]+Aij[[b,a]])[[i,j]],2], {a,3}, {b,3}];
              D3=Table[Re@Total[(Aij[[4,a]]-Aij[[a,4]])[[i,j]],2], {a,3}];
              J0 IdentityMatrix[3] + Ja + {{0, D3[[3]], -D3[[2]]}, 
                                           {-D3[[3]], 0, D3[[1]]}, 
                                           {D3[[2]], -D3[[1]], 0}}, {i, Length@orb}, {j, Length@orb}]
  ];
  Ji0 = Table[GetJi0Bonds[spg, {{{orb[[i, 1]], orb[[j, 1]]}, R, Jij[[i, j]]}}, pos], {i, Length@orb}, {j, Length@orb}];
  Return[Ji0]
]

(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
