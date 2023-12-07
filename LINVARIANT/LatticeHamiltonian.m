BeginPackage["LINVARIANT`LatticeHamiltonian`",{"LINVARIANT`Structure`","LINVARIANT`GroupTheory`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
ReadQpoint                   ::usage "ReadQpoint[filename]"
ReadBands                    ::usage "ReadBands[filename]"
ModesSum                     ::usage "ModesSum[list, Modes]"
ImposeTranslationalSymmetry  ::usage "ImposeTranslationalSymmetry[FC]"
FCSymmetrization             ::usage "FCSymmetrization[pos0, spg0, tgrp, FC]"
FC2Fi0Bonds                  ::usage "FC2Fi0Bonds[FC, pos, NSC]"
GetHessianOnSite             ::usage "GetHessianOnSite[InvList, vars, CoeffStr]"
Invariant2Fi0Bonds           ::usage "Invariant2Fi0Bonds[InvList, vars, CoeffStr]"
GetDynamicMatrix             ::usage "GetDynamicMatrix[Fi0, pos, q]"
PhononK                      ::usage "PhononK[Fi0, pos, k]"
PhononBandsPlot              ::usage "PhononBandsPlot[Fi0, pos, klist, kintv]"
Unfolding                    ::usage "Unfolding[PhononNK, G, DMPos, SPCPos, NSC]"
UnfoldingPW                  ::usage "UnfoldingPW[PhononNK, Gsc, DMPos, SPCPos, NSC]"
PhononUnfolding              ::usage "PhononUnfolding[Fi0, DMPos, StdPos, q, Qbz, unfoldDim]"
GetEwaldMatrix               ::usage "GetEwaldMatrix[NGrid, tol]"
ArrayFold2D                  ::usage "ArrayFold2D[list, n]"
LoadEwaldMat                 ::usage "LoadEwaldMat[dir, file, grid]"

(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
ReadQpoint[filename_] := Module[{QpointPhonon, EigenValuesVectors},
    QpointPhonon = Import[filename, "YAML"];
    EigenValuesVectors = QpointPhonon[[4]][[2]][[1]][[2]][[2]];
    Return[{"reciprocal_lattice"/.QpointPhonon,Table[{"frequency", "eigenvector", "group_velocity"} /. EigenValuesVectors[[i]], {i, Length@EigenValuesVectors}]}];
]

ReadBands[filename_] := Module[{PhonopyBands},
  PhonopyBands = Import[filename, "YAML"];
  Return[Table[{"distance" /. ph, "q-position" /. ph, #2 & @@@ Flatten["band" /. ph]}, {ph, PhonopyBands[[10]][[2]]}]];
]

ModesSum[list_, Modes_] := Module[{modes, i},
  Return[Sum[i[[2]] #\[Transpose][[1]] & /@ Modes[[i[[1]]]][[2]], {i, list}]];
]

FCSymmetrization[pos0_, spg0_, tgrp_, FC_, OptionsPattern[{"round" -> 10^-4, "iter" -> 3}]] := Module[{latt, sites, atomics, label, op, ni, nj, na, nb, dim, SymFC},
  {latt, sites} = pos0;
  atomics = IdentityMatrix[3 Length[sites]];
  label = Flatten[Table[# <> \[Alpha], {\[Alpha], {"x", "y", "z"}}] & /@ (pos0[[2]]\[Transpose][[2]])];
  op = Chop[GetMatrixRep[spg0, tgrp, pos0, atomics, label, "disp"], OptionValue["round"]];
  {ni, nj, na, nb} = Dimensions[FC];
  dim = ni na;
  SymFC = ArrayFlatten[FC];
  Do[SymFC = (# - Total[#]/dim) & /@ SymFC;
     SymFC = ((# - Total[#]/dim) & /@ (SymFC\[Transpose]))\[Transpose];
     SymFC = Chop[Total[# . SymFC . (#\[Transpose]) & /@ op]/Length[op]], {OptionValue["iter"]}];
  SymFC = (Partition[#, 3] & /@ ((Partition[#, 3] & /@ SymFC)\[Transpose]))\[Transpose];
  Return[SymFC]
]

ImposeTranslationalSymmetry[FC_] := Module[{ni, nj, na, nb, dim, SymFC},
  {ni, nj, na, nb} = Dimensions[FC];
  dim = ni na;
  Do[SymFC = (# - Total[#]/dim) & /@ ArrayFlatten[FC];
     SymFC = ((# - Total[#]/dim) & /@ (SymFC\[Transpose]))\[Transpose], {3}];
  SymFC = (Partition[#,3] & /@ ((Partition[#, 3] & /@ SymFC)\[Transpose]))\[Transpose];
  Return[SymFC]
]

ArrayFold2D[list_, n_] := Module[{},
    (Partition[#,n] & /@ ((Partition[#, n] & /@ list)\[Transpose]))\[Transpose]
]

GetDynamicMatrixOld[FC_, pos_, q_, NSC_, OptionsPattern["tol" -> 0.1]] := Module[{latt, sites, p2s, Nij, HermiteDM, dm, ppos, spos, NumPpos, NumSpos, mass, NeighborList, ImageVectorList, vectors, multi, Nx, Ny, Nz, i, j, \[Alpha], \[Beta], s1, s2, s3, posdiff, Vasp2THz = 15.633302300230191},
  Nij = DiagonalMatrix[NSC];
  {latt, sites} = pos;
  {Nx, Ny, Nz} = NSC;
  NeighborList = Flatten[Table[{i, j, k} - {1, 1, 1}, {i, Nx}, {j, Ny}, {k, Nz}], 2];
  spos = Join[{latt}, {{Rationalize@Chop[Inverse[latt].Nij.latt.#1], #2} & @@@ sites}];
  ppos = Join[{Inverse[Nij].spos[[1]]}, {DeleteDuplicates[{Rationalize@Mod[Chop[Inverse[latt].Nij.latt.#1], 1], #2} & @@@ sites, Norm[#1[[1]] - #2[[1]]] < OptionValue["tol"] &]}];
  NumPpos = Length[ppos[[2]]];
  NumSpos = Length[spos[[2]]];
  p2s = Association[Thread[NeighborList -> Table[Association[Table[i -> First@First@Position[Rationalize[Norm[ppos[[2]][[i]][[1]] + shift - #1]] & @@@spos[[2]], 0], {i, NumPpos}]], {shift, NeighborList}]]];
  ImageVectorList = Table[vectors = {Norm[ppos[[1]]\[Transpose].#1], #1, #2} & @@@ Flatten[Table[{spos[[2]][[j, 1]] + NSC {s1, s2, s3} - spos[[2]][[i, 1]], {s1, s2, s3}}, {s1, {-1, 0, 1}}, {s2, {-1, 0, 1}}, {s3, {-1, 0, 1}}], 2]; Association[GroupBy[vectors, First]][Min[vectors\[Transpose][[1]]]], {i, NumSpos}, {j, NumSpos}];
  mass = Table[QuantityMagnitude[ElementData[i, "AtomicMass"]], {i, ppos[[2]]\[Transpose][[2]]}];
  dm = ArrayFlatten[Total[Table[multi = 1/Length[ImageVectorList[[p2s[{0, 0, 0}][i], p2s[shift][j]]]];
    Sum[1/Sqrt[mass[[i]] mass[[j]]] multi FC[[p2s[{0, 0, 0}][i], p2s[shift][j]]][[\[Alpha], \[Beta]]] Exp[I 2 Pi q.posdiff[[2]]], {posdiff, ImageVectorList[[p2s[{0, 0, 0}][i], p2s[shift][j]]]}], {shift,Keys@p2s}, {i, NumPpos}, {j, NumPpos}, {\[Alpha], 3}, {\[Beta], 3}]]];
  HermiteDM = 1/2 (dm + dm\[ConjugateTranspose]) Vasp2THz^2;
  Return[{HermiteDM, ppos}]
]

FC2Fi0Bonds[FC_, pos_, NSC_, OptionsPattern[{"tol" -> 0.1, "rotate"->False}]] := Module[{latt, sites, p2s, Fi0, Nij, ppos, spos, poscar, NumPpos, NumSpos, NeighborList, shift, imagevec, ImageVectorList, vectors, multi, Nx, Ny, Nz, i, j, s1, s2, s3, bonds},
  Nij = DiagonalMatrix[NSC];
  {latt, sites} = pos;
  {Nx, Ny, Nz} = NSC;
  NeighborList = If[OptionValue["rotate"],
                    Flatten[Table[{k, j, i} - {1, 1, 1}, {i, Nz}, {j, Ny}, {k, Nx}], 2],
                    Flatten[Table[{i, j, k} - {1, 1, 1}, {i, Nx}, {j, Ny}, {k, Nz}], 2]];

  If[Det[Nij] < 0,
     spos = Join[{latt}, {{Rationalize@Chop[Inverse[latt].Inverse[Nij].latt.#1], #2} & @@@ sites}];
     ppos = Join[{Nij.latt}, {DeleteDuplicates[{Rationalize@Mod[Chop[Inverse[latt].Inverse[Nij].latt.#1], 1], #2} & @@@ sites, ((Norm[#1[[1]] - #2[[1]]] < OptionValue["tol"])&&(#1[[2]]===#2[[2]])) &]}];
     poscar = ppos,
     ppos=pos;
     spos = If[OptionValue["rotate"],
               {Nij.latt, Flatten[Table[{#[[1]] + {k - 1, j - 1, i - 1}, #[[2]]}, {i, Nz}, {j, Ny}, {k, Nx}] &/@ sites, 3]},
               {Nij.latt, Flatten[Table[{#[[1]] + {i - 1, j - 1, k - 1}, #[[2]]}, {i, Nx}, {j, Ny}, {k, Nz}] &/@ sites, 3]}];
     poscar = spos
  ];

  NumPpos = Length[ppos[[2]]];
  NumSpos = Length[spos[[2]]];
  p2s = Association[Thread[NeighborList -> Table[Association[Table[i -> First@First@Position[Rationalize[Norm[ppos[[2]][[i]][[1]] + shift - #1]] & @@@ (spos[[2]]), 0], {i, NumPpos}]], {shift, NeighborList}]]];

  ImageVectorList = Table[bonds = PbcDiff[spos[[2]][[j, 1]] - ppos[[2]][[i, 1]], NSC];
                          Flatten[Table[If[Chop[Norm[bonds] - Norm[bonds + NSC {s1, s2, s3}]] == 0., {Norm[bonds], bonds, NSC {s1, s2, s3}}, ## &[]], {s1, {-1, 0, 1}}, {s2, {-1, 0, 1}}, {s3, {-1, 0, 1}}], 2], {i, NumPpos}, {j, NumSpos}];

  Fi0 = Flatten[SortBy[#, #[[1]] &] & /@ Values@GroupBy[Flatten[Table[multi = Length[ImageVectorList[[i, p2s[shift][j]]]];
        Table[{{i, j}, PbcDiff[shift + imagevec[[3]], NSC], multi, FC[[p2s[{0,0,0}][i], p2s[shift][j]]]}, {imagevec, ImageVectorList[[i, p2s[shift][j]]]}], {shift, Keys@p2s}, {i, NumPpos}, {j, NumPpos}], 3], #[[2]] &], 1];

  Return[{Fi0, poscar}]
]

GetDynamicMatrix[Fi0_, pos_, q_, OptionsPattern[{"mass"->{}}]] := Module[{latt, sites, HermiteDM, dm, dmBlock, NumPos, mass, ele, amu, s, mfactor, Vasp2THz = 15.633302300230191, i, j, TB},
  {latt, sites} = pos;
  NumPos = Length[sites];
  mass = If[OptionValue["mass"]==={},
            Table[amu=ToLowerCase[StringJoin[StringCases[i, RegularExpression["[a,m,u,A,M,U]"]]]]==="amu";
               ele=If[amu, "C", First@StringCases[i, RegularExpression["[[:upper:]]*[[:lower:]]*"]]];
               s = StringCases[i, RegularExpression["\\d+"]];
               mfactor=If[s==={}, 1, ToExpression[First@s]] If[amu, 1/12, 1];
               QuantityMagnitude[ElementData[ele, "AtomicMass"]], {i, sites\[Transpose][[2]]}],
            OptionValue["mass"]];
  dmBlock = Merge[{#1 -> 1/Sqrt[mass[[#1[[1]]]] mass[[#1[[2]]]]] #4/#3 Exp[I 2 Pi q.(sites[[#1[[2]], 1]] + #2 - sites[[#1[[1]], 1]])]} & @@@ Fi0, Total];
  dm = ArrayFlatten[Table[dmBlock[{i, j}], {i, NumPos}, {j, NumPos}]];
  HermiteDM = 1/2 (dm + dm\[ConjugateTranspose]) Vasp2THz^2;
  Return[HermiteDM]
]

PhononBandsPlot[Fi0_, pos_, klist_, kintv_, OptionsPattern[{"range" -> All, "mass1"->False, "AspectRatio" -> 1/GoldenRatio, "imagesize" -> Medium, "plotstyle" -> Automatic}]] := Module[{latt, sites, pdata, q, DM, sol, kpath, xticks, BandsPlot},
  {latt, sites} = pos;
  {kpath, xticks} = GetKpath[latt, klist, kintv];
  pdata = Table[DM = GetDynamicMatrix[Fi0, pos, q[[2]], If[OptionValue["mass1"], "mass" -> ConstantArray[1, Length@sites], {}]];
                sol = Eigenvalues[DM];
                Sort[{q[[1]], ReIm[Sqrt[#]].{1,-1}} &/@ sol], {q, kpath}]; 
  BandsPlot = ListPlot[pdata\[Transpose], 
                       PlotStyle -> OptionValue["plotstyle"],
                       Joined -> True, 
                       PlotRange -> {All,OptionValue["range"]}, 
                       AspectRatio -> OptionValue["AspectRatio"], 
                       Frame -> True, 
                       ImageSize -> OptionValue["imagesize"],
                       GridLines -> {{xticks\[Transpose][[1]], ConstantArray[Thick, Length[xticks]]}\[Transpose], Automatic}, 
                       ImageSize -> Medium, 
                       FrameTicks -> {{Automatic, None}, {xticks, None}}];
  Return[BandsPlot]
]

PhononK[Fi0_, pos_, k_?ListQ] := Module[{latt, sites, DM, sol, w2, phononfield},
  If[Depth[k] > 2, 
     PhononK[Fi0, pos, #] &/@ k,
     {latt, sites} = pos;
     DM = GetDynamicMatrix[Fi0, pos, k];
     sol = Chop[Eigensystem[DM]];
     w2 = sol[[1]];
     phononfield = {sites, Partition[#1, 3]}\[Transpose] & /@ sol[[2]]];
     Return[{w2, phononfield}]
]

Unfolding[PhononNK_, G_, DMPos_, SPCPos_, NSC_] := Module[{w2, K, phonon, PCTrans, PC2SCList, ReorderPosMap, PosStd, i, t, ElementMatched, PCTranMap, Tphonon, WNJG},
  {w2, K, phonon} = PhononNK;
  PCTrans = Join[{{0, 0, 0}}, Inverse[DiagonalMatrix[NSC]]];
  PC2SCList = DeleteDuplicates[Total[#] & /@ DeleteCases[Subsets[PCTrans, 3], {}]];
  ReorderPosMap = Association[#1 -> #2 & @@@ (PosMatchTo[IdentityMatrix[3], DMPos\[Transpose][[1]], SPCPos[[2]]\[Transpose][[1]]][[1]])];
  PosStd = Table[SPCPos[[2]][[ReorderPosMap[i]]], {i, Length[SPCPos[[2]]]}];
  ElementMatched = AllTrue[Table[DMPos[[i]][[2]] === SPCPos[[2]][[ReorderPosMap[i]]][[2]], {i, Length[SPCPos[[2]]]}], # === True &];
  If[Not@ElementMatched, Print["The selected quasi primary cell doesn't have a matched atom list!"]; Abort[]];
  PCTranMap = Association[Table[t -> Association[#1 -> #2 & @@@ (PosMatchTo[IdentityMatrix, PosStd\[Transpose][[1]], Mod[# + t, 1] & /@ (PosStd\[Transpose][[1]])][[1]])], {t, PC2SCList}]];
  WNJG = Norm[1/Length[PC2SCList] Sum[Tphonon =Table[{phonon\[Transpose][[1]][[i]], phonon\[Transpose][[2]][[PCTranMap[t][i]]]}, {i, Length[phonon\[Transpose][[2]]]}];
      Flatten[Tphonon\[Transpose][[2]]] Exp[-I 2 Pi (K + G).t], {t, PC2SCList}]]^2;
  Return[{Inverse[DiagonalMatrix[NSC]].(K+G), w2, Chop[WNJG]}]
]

UnfoldingPW[PhononNK_, Gsc_, DMPos_, SPCPos_, NSC_] := Module[{w2, q, Gpc, phonon, PCTrans, GpcList, BGpc\[Alpha], WqGsc}, 
  {w2, q, phonon} = PhononNK;
  PCTrans = Join[{{0, 0, 0}}, DiagonalMatrix[NSC]];
  GpcList = DeleteDuplicates[Total[#] & /@ DeleteCases[Subsets[PCTrans, 3], {}]];
  BGpc\[Alpha] = Flatten[Table[Total[1/Sqrt[Length[phonon]] #2 Exp[-I 2 Pi (Gsc + Gpc).#1[[1]]] & @@@ phonon], {Gpc, GpcList}]];
  WqGsc = Norm[BGpc\[Alpha]]^2;
  Return[{Inverse[DiagonalMatrix[NSC]].(q + Gsc), w2, Chop[WqGsc]}]
]

PhononUnfolding[Fi0_, DMPos_, StdPos_, q_, Qbz_, unfoldDim_, OptionsPattern[{"pw" -> True}]] := Module[{w2, phonon},
  If[Depth[q] > 2, 
  PhononUnfolding[Fi0, DMPos, StdPos, #, Qbz, unfoldDim, "pw" -> OptionValue["pw"]] &/@ q,
  {w2, phonon} = PhononK[Fi0, DMPos, q];
  ParallelTable[
    If[OptionValue["pw"], UnfoldingPW[{w2[[i]], q, phonon[[i]]}, Qbz, DMPos, StdPos, unfoldDim], 
                            Unfolding[{w2[[i]], q, phonon[[i]]}, Qbz, DMPos, StdPos, unfoldDim]], 
    {i, Length[w2]}, DistributedContexts -> {"LINVARIANT`LatticeHamiltonian`Private`"}]]
]

GetEwaldMatrix[NGrid_, epinf_, Q_, alat_, tol_] := Module[{gcut, ix, iy, iz, NGridx, NGridy, NGridz, SuperLattice, mg, dpij, residue, \[Eta], cc, OriginQ, GList0, GList1, G, RecLatt, ig1, ig2, ig3, fact, angstrom2bohr}, 
  angstrom2bohr = 1.88973;
  {NGridx, NGridy, NGridz} = NGrid;
  SuperLattice = N[IdentityMatrix[3] {NGridx, NGridy, NGridz}];
  RecLatt = 2 Pi Inverse[SuperLattice];
  \[Eta] = N[Sqrt[-Log[tol]]];
  gcut = 2 \[Eta]^2;
  mg = Ceiling[gcut Norm[#]/(2 Pi)] & /@ SuperLattice;
  residue = (4 \[Eta]^3)/(3 Sqrt[Pi]); 
  cc = (8 Pi)/Det[SuperLattice];
  fact = Q^2/(epinf*angstrom2bohr*alat*2);
  GList0 = Flatten[Table[G = RecLatt.{0, ig2, ig3}; 
                         If[10^-8 < Norm[G]^2 < gcut^2, G, ## &[]], {ig2, -mg[[2]], mg[[2]]}, {ig3, -mg[[3]], mg[[3]]}], 1];
                         GList1 = Flatten[Table[G = RecLatt.{ig1, ig2, ig3}; 
                         If[10^-8 < Norm[G]^2 < gcut^2, G, ## &[]], 
                         {ig1, 1, mg[[1]]}, {ig2, -mg[[2]], mg[[2]]}, {ig3, -mg[[3]], mg[[3]]}], 2];
  dpij = ParallelTable[Sum[0.5 cc Cos[G.{ix-1, iy-1, iz-1}] Exp[-(G.G/(4 \[Eta]^2))] G\[TensorProduct]G/G.G, 
                           {G, GList0}] + Sum[cc Cos[G.{ix-1, iy-1, iz-1}] Exp[-(G.G/(4 \[Eta]^2))] G\[TensorProduct]G/G.G, 
                           {G, GList1}], {ix, NGridx}, {iy, NGridy}, {iz, NGridz}, DistributedContexts -> {"LINVARIANT`LatticeHamiltonian`Private`"}];
  dpij[[1, 1, 1]] = dpij[[1, 1, 1]] - residue*fact*IdentityMatrix[3];
  Return[fact*dpij]
]

Invariant2Fi0Bonds[InvList_, vars0_, CoeffStr_] := Module[{x, i, NeighbourList, varsijk, Fi0, vars, nn}, 
  NeighbourList = DeleteDuplicates[Level[#, {1}][[3 ;; 5]] & /@ Variables[If[ListQ[InvList], InvList /. Subscript[ToExpression@"\[Epsilon]0", ___] -> 1, InvList /. Subscript[ToExpression@"\[Epsilon]0", ___] -> 1 /. ToExpression[CoeffStr][___] -> 1]]];
  vars = vars0 /. {Subscript[x_, i_, ix_, iy_, iz_] -> Subscript[x, i]};
  Fi0 = Table[varsijk = vars /. {Subscript[x_, i_] -> (Subscript[x, i, #1, #2, #3] & @@ nn)};
              {nn, Table[D[If[ListQ[InvList], InvList.(ToExpression[CoeffStr][#] & /@ Range[Length@InvList]), InvList], v1, v2], {v1, vars0}, {v2, varsijk}]}, 
              {nn, NeighbourList}];
Return[Fi0]
]

GetHessianOnSite[InvList_, vars0_, CoeffStr_] := Module[{x, i, Hi0i0},
  Hi0i0 = Table[D[If[ListQ[InvList], InvList.(ToExpression[CoeffStr][#] & /@ Range[Length@InvList]), InvList], v1, v2], {v1, vars0}, {v2, vars0}];
  Return[Hi0i0]
]

LoadEwaldMat[dir_, file_, grid_, OptionsPattern[{"xyz"->True}]] := Module[{ix, iy, iz, Nx, Ny, Nz, Ewald, EwaldMat},
  {Nx, Ny, Nz} = grid;
  Ewald = Partition[Chop@Import[dir <> "/" <> file, "Table", RecordLists -> True], 3];
  EwaldMat = If[OptionValue["xyz"], 
                Table[Ewald[[Ny Nx (iz - 1) + Nx (iy - 1) + ix]], {ix, Nx}, {iy, Ny}, {iz, Nz}],
                Table[Ewald[[Ny Nx (iz - 1) + Nx (iy - 1) + ix]], {iz, Nz}, {iy, Ny}, {ix, Nx}]];
  Return[EwaldMat]
]

(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
