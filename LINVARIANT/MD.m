BeginPackage["LINVARIANT`MD`", {"LINVARIANT`MathematicaPlus`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
ReadFieldTrajectory  ::usage "ReadFieldTrajectory[file, dof, ifield, grid, time]"

(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
ReadFieldTrajectory[file_, dof_, ifield_, grid_, time_ : {1}] := Module[{f, t, i, nlines, nb, tmp, out, jump},
  f = OpenRead[file];
  nlines = (dof + 1) Times @@ grid + 3;
  Do[ReadLine[f], {nlines}];
  jump = StreamPosition[f] + 1;
  out = Table[SetStreamPosition[f, (t - 1) jump];
              tmp = Table[ReadLine[f], {nlines}];
              Which[ifield == 0, tmp[[1]], 
                    ifield > dof, tmp[[2 ;; 3]], 
                    True, ParseFortranNumber@StringSplit[#[[ifield + 1]]] & /@ Partition[tmp[[4 ;;]], dof + 1]], {t, time}];
  Close[f];
  Return[out]
]
(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
