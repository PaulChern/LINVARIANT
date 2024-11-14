BeginPackage["LINVARIANT`QE`", {"LINVARIANT`INVARIANT`", "LINVARIANT`Structure`", "LINVARIANT`MathematicaPlus`", "LINVARIANT`Parser`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
ParseQEStruct                ::usage "ParseQEStruct[xml]"
ParseQEBands                 ::usage "ParseQEBands[xml]"
ParseQEFC                    ::usage "ParseQEFC[file]"
(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
ParseQEStruct[xml_] := Module[{posxml, cellxml, latt, pos, bohr2angstrom},
  bohr2angstrom = 0.529177;
  cellxml = Cases[Cases[xml, XMLElement["output", {}, x_] :> x, Infinity], XMLElement["cell", {}, x_] :> x, Infinity];
  posxml = Cases[Cases[xml, XMLElement["output", {}, x_] :> x, Infinity], XMLElement["atomic_positions", {}, x_] :> x, Infinity];
  latt = Flatten[Table[Cases[cellxml, XMLElement["a" <> ToString[i], {}, x_] :> ParseFortranNumber[x], Infinity], {i, 3}], 2]; 
  pos = Cases[posxml, XMLElement["atom", label_, pos_] :> {Flatten[ParseFortranNumber[pos]], "name" /. label}, Infinity];
  Return[{bohr2angstrom latt, {Mod[Inverse[latt\[Transpose]] . #1, 1], #2} & @@@ pos}]
]

ParseQEBands[xml_] :=Module[{eigen, bohr2angstrom, cellxml, latt, kpath, fermi, Ry2eV},
  bohr2angstrom = 0.529177;
  Ry2eV = 27.2114;
  fermi=First@First@Cases[xml,XMLElement["fermi_energy", {}, x_] :> ParseFortranNumber[x],Infinity];
  cellxml=Cases[Cases[xml, XMLElement["output", {}, x_] :> x, Infinity],XMLElement["cell", {}, x_] :> x, Infinity];
  latt=Flatten[Table[Cases[cellxml,XMLElement["a" <> ToString[i], {}, x_] :> ParseFortranNumber[x],Infinity], {i, 3}], 2];
  kpath = Cases[Cases[Cases[xml, XMLElement["output", {}, x_] :> x, Infinity],XMLElement["starting_k_points", {}, x_] :> x, Infinity],XMLElement["k_point", __, x_] :> Flatten@ParseFortranNumber[x],Infinity];
  eigen = #\[Transpose] & /@ ({Cases[Cases[xml, XMLElement["output", {}, x_] :> x, Infinity],XMLElement["eigenvalues", __, x_] :> Ry2eV Flatten@ParseFortranNumber[x], Infinity], Cases[Cases[xml, XMLElement["output", {}, x_] :> x, Infinity], XMLElement["occupations", __, x_] :> Flatten@ParseFortranNumber[x], Infinity]}\[Transpose]);
  Return[<|"latt" -> latt, "k" -> kpath, "fermi" -> fermi, "up" -> eigen|>]
]

ParseQEFC[file_, OptionsPattern[{"rotate"->False}]] := Module[{fc, NumType, NumAtoms, fcdata, H0, Hij, nq1, nq2, nq3, i, j, k, NL, iq1, iq2, iq3, Q, Rybb2eVaa},
  Rybb2eVaa=13.6056980659/0.529177^2;
  fc = OpenRead[file];
  {NumType, NumAtoms} = ToExpression[StringSplit[ReadLine[fc]]][[1 ;; 2]];
  Do[ReadLine[fc], {3}];
  Do[ReadLine[fc], {NumType + NumAtoms}];
  Q = ReadLine[fc];
  If[Q==="T", Do[ReadLine[fc], {3 + 4 NumAtoms}]];
  {nq1, nq2, nq3} = ToExpression[StringSplit[ReadLine[fc]]];
  fcdata = Partition[ReadList[fc, Number, RecordLists -> True], nq1 nq2 nq3 + 1];
  Close[fc];

  H0 = Association@Table[{{iq1-1, iq2-1, iq3-1}->Table[Table[Chop[Rybb2eVaa fcdata[[(i - 1) 3 NumAtoms NumAtoms + (j - 1) NumAtoms NumAtoms + (ia - 1) NumAtoms + jb, (iq1 - 1) nq2 nq3 + (iq2 - 1) nq3 + iq3 + 1]][[4]]], {i, 3}, {j, 3}], {ia, NumAtoms}, {jb, NumAtoms}]}, {iq3, nq1}, {iq2, nq2}, {iq1, nq3}];

  NL = If[OptionValue["rotate"], Flatten[Table[{k, j, i} - {1, 1, 1}, {i, nq3}, {j, nq2}, {k, nq1}], 2], 
                                 Flatten[Table[{i, j, k} - {1, 1, 1}, {i, nq1}, {j, nq2}, {k, nq3}], 2]];
  Hij = ArrayFlatten[Table[H0[Mod[j - i, 3]], {i, NL}, {j, NL}]];

  Return[Hij]
]
(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
