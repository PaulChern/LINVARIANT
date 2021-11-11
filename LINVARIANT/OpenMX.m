BeginPackage["LINVARIANT`OpenMX`", {"LINVARIANT`INVARIANT`", "LINVARIANT`Structure`", "LINVARIANT`MathematicaPlus`", "LINVARIANT`Parser`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
ParseOpenMXBand              ::usage "ParseOpenMXBand[file]"
(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
ParseOpenMXBand[file_] := Module[{band, templine, ip, tmp, nbands, Efermi, Blatt, npath, kpts, ndiv, qpt1, qpt2, hkpts, xticks, dk, kpath, klist, distance, data},
  hkpts = {};
  band = OpenRead[file];
  templine = StringSplit@ReadNonEmptyLine[band];
  {nbands, tmp, Efermi} = ParseFortranNumber[templine];
  templine = StringSplit@ReadNonEmptyLine[band];
  Blatt = Partition[ParseFortranNumber[templine], 3];
  data = <|"latt" -> 2 Pi Inverse[Blatt]|>;
  npath = ParseFortranNumber@ReadNonEmptyLine[band];
  klist = Flatten[Table[templine = StringSplit[ReadNonEmptyLine[band]];
                        ndiv = ParseFortranNumber[templine[[1]]];
                        {qpt1, qpt2} = Partition[ParseFortranNumber[templine[[2 ;; 7]]], 3];
                        AppendTo[hkpts, {{qpt1, qpt2}, templine[[8 ;; 9]]}\[Transpose]];
                        Join[{qpt1}, Table[N[qpt1 + i (qpt2 - qpt1)/(ndiv - 1)], {i, ndiv - 1}]], {ip, npath}], 1];
  data["k"] = klist;
  data["hkpts"] = Flatten[hkpts, 1];
  data["up"] = Table[templine = ReadNonEmptyLine[band];
                     templine = StringSplit@ReadNonEmptyLine[band];
                     {27.2114 #, If[# <= Efermi, 1, 0]} & /@ ParseFortranNumber[templine], {i, Length[data[["k"]]]}];
  Close[band];
  Return[data]
]

(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
