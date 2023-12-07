BeginPackage["LINVARIANT`Salmon`", {"LINVARIANT`INVARIANT`", "LINVARIANT`Structure`", "LINVARIANT`MathematicaPlus`", "LINVARIANT`Parser`", "LINVARIANT`Vasp`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
ImportSalmonTrj            ::usage "ImportSalmonTrj[file]"
Salmon2POSCAR              ::usage "Salmon2POSCAR[file]"
SalmonMonitor              ::usage "SalmonMonitor[file, txt, col]"
Export2Salmon              ::usage "Export2Salmon[dir, inpf, pos]"
RunSalmonHessian           ::usage "RunSalmonHessian[dir, inpf]"
SalmonVaspKPT              ::usage "SalmonVaspKPT[dir, IBZKPT]"

(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
ImportSalmonTrj[inp_, xyz_, OptionsPattern[{"round" -> 10.0^-16, "svf" -> False}]] := Module[{latt, xyzdata, natom, sites},
  latt = DiagonalMatrix@Flatten[If[Quiet[StringTake[StringTrim[#], 7]] === "al(1:3)", ToExpression@StringSplit[StringSplit[#, "="][[2]], ","], ## &[]] & /@ ReadList[inp, "String"]];

  xyzdata = ReadList[xyz, "String", RecordLists -> True];

  natom = First@ToExpression[xyzdata[[1]]];

  sites = Partition[If[Length[Flatten[StringSplit[#]]] != 1 && StringTake[StringTrim@First[#], 1] != "#", {Round[Mod[Inverse[latt].ToExpression@Flatten[StringSplit[#]][[2;;4]], 1], OptionValue["round"]], Flatten[StringSplit[#]][[1]], If[OptionValue["svf"], Quiet@ToExpression@Flatten[StringSplit[#]][[6;;8]], ##&[]], If[OptionValue["svf"], Quiet@ToExpression@Flatten[StringSplit[#]][[10;;12]], ##&[]]}, ## &[]] & /@ xyzdata, natom];

  Return[{latt, #} & /@ sites]
]

Salmon2POSCAR[inpf_] := Module[{inp, x, latt, sites, natom, i, is, il},
  inp = ReadList[inpf, "String"];
  natom = ToExpression@First@Flatten[If[StringContainsQ[#, "natom"], StringCases[StringSplit[#], x__ /; PrimeQ[Quiet@ToExpression[x]]], ## &[]] & /@ inp];
  is = First@Flatten@MapIndexed[If[#1 != {} && First[First@#1] <= 2, #2, ## &[]] &, StringPosition[inp, "atomic_red_coor"]];
  il = First@MapIndexed[If[#1, First@#2, ## &[]] &, StringContainsQ[inp, "al(1:3)"]];

  sites = Table[line = StringSplit[inp[[is+i]]]; {ToExpression[line[[2;;4]]], line[[1]]}, {i, natom}];
  latt = DiagonalMatrix@ToExpression["{" <> StringSplit[inp[[il]], "="][[2]] <> "}"];
  Return[{latt, sites}]
]

Export2Salmon[dir_, inpf_, pos_, OptionsPattern[{"fname" -> "mma.inp"}]] := Module[{inp, natom, x, is, il, pseudo, i, latt, sites, out, sub, line},
  inp = ReadList[dir <> "/" <> inpf, "String"];
  {latt, sites} = pos;

  natom = ToExpression@First@Flatten[If[StringContainsQ[#, "natom"], StringCases[StringSplit[#], x__ /; PrimeQ[Quiet@ToExpression[x]]], ## &[]] & /@ inp];
  is = First@Flatten@MapIndexed[If[#1 != {} && First[First@#1] <= 2, #2, ## &[]] &, StringPosition[inp, "atomic_red_coor"]];
  il = First@MapIndexed[If[#1, First@#2, ## &[]] &, StringContainsQ[inp, "al(1:3)"]];

  pseudo = Association[MapIndexed[If[#1, First@StringTake[StringCases[inp[[First@#2]], RegularExpression["\\/\\w[^_/]*\\_"]], {2, -2}] -> First@ToExpression@StringCases[inp[[First@#2]], RegularExpression["\\(\\d\\)"]], ## &[]] &, StringContainsQ[inp, "file_pseudo"]]];
  sub = Table[line = StringRiffle[Flatten[{"  " <> sites[[i, 2]], ToString@NumberForm[#, {20, 16}] & /@ (sites[[i, 1]]), ToString@pseudo[sites[[i, 2]]], "y"}]];
              inp[[is + i]] -> line, {i, natom}];

  AppendTo[sub, inp[[il]] -> "  al(1:3) = " <> StringRiffle[ToString@NumberForm[Norm[#], {20, 16}] & /@ latt, ", "]];
  out = Flatten[If[StringTake[#, 1] === "&", {"  ", #}, {#}] & /@ (inp /. sub)];
  
  Export[dir <> "/" <> OptionValue["fname"], out, "Table"]
]

SalmonMonitor[file_, txt_, col_, OptionsPattern[{"imagesize"->300, "log" -> True, "range"->All}]] := Module[{plt, data},
  Run["grep '" <> txt <> "' " <> file <> " | awk '{print " <> StringRiffle["$" <> ToString[#] & /@ col, ", "] <> "}' > " <> file <> ".mma.dat"];
  plt = If[OptionValue["log"], ListLogPlot, ListPlot];
  data = Transpose[Import[file <> ".mma.dat"]];
  plt[Join[{ConstantArray[Min@data, Length[First@data]]}, data], ImageSize->OptionValue["imagesize"], Frame -> True, PlotRange -> OptionValue["range"], GridLines -> Automatic, Joined -> True, PlotMarkers -> {Automatic, 5}, PlotStyle -> {Black, Blue, Red}]
]

RunSalmonHessian[dir_, inpf_, potim_, npt_, OptionsPattern[{"mod"->False}]] := Module[{natom, latt, sites, sub, s, i, j, k},
  {latt, sites} = Salmon2POSCAR[dir <> "/" <> inpf];
  natom = Length[sites];
  Do[s = sites[[i]];
     sub = s -> If[OptionValue["mod"], {Mod[s[[1]] + (k - npt - 1) potim Inverse[latt\[Transpose]] . (IdentityMatrix[3][[j]]), 1], s[[2]]},
                                       {s[[1]] + (k - npt - 1) potim Inverse[latt\[Transpose]] . (IdentityMatrix[3][[j]]), s[[2]]}];
     Export2Salmon[dir,inpf,{latt,sites/.sub},"fname"->"H"<>ToString[i]<>"_"<>ToString[j]<>"-"<>ToString[k]<>".inp"];
     ExportPOSCAR[dir,"H"<>ToString[i]<>"_"<>ToString[j]<>"-"<>ToString[k]<>".vasp",{latt,sites/.sub}], {i,natom}, {j,3}, {k,2*npt+1}];
]

SalmonVaspKPT[dir_, IBZKPT_, OptionsPattern[{"fname" -> "kpts.dat"}]] := Module[{klist, kw, i, w},
  klist = ReadList[dir <> "/" <> IBZKPT, "Number", RecordLists -> True];
  w = klist\[Transpose][[4]]/Total[1.0 klist\[Transpose][[4]]];
  kw = Prepend[Table[ToString[#] & /@ Flatten[{i, NumberForm[#, {16, 15}] & /@ (klist[[i]][[1 ;; 3]]), NumberForm[w[[i]], {16, 15}]}], {i, Length@klist}], {Length@klist}];
  Export[dir <> "/" <> OptionValue["fname"], kw]
]
(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
