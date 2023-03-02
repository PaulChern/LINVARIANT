BeginPackage["LINVARIANT`MultipoleExpansion`", {"LINVARIANT`GroupTheory`", "LINVARIANT`SphericalHarmonics`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
GetMultipoleExpansionCoeff ::usage "GetMultipoleExpansionCoeff[latt, grp, fieldseed, Jcut, Harmonic]"
GetCrystalFieldHam         ::usage "GetCrystalFieldHam[Coeff, l]"

(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
GetMultipoleExpansionCoeff[latt_, grp_, fieldseed_, Jcut_, Harmonic_:"Tesseral"] := Module[{i, field, g, J, m, eqns, sol, coeff, Alm, vars, A, jj, mm},
  field = Union@Flatten[Table[{latt\[Transpose].First[#].fieldseed[[i,1]], fieldseed[[i,2]]} &/@ xyz2RotT[Keys[grp]], {i, Length[fieldseed]}], 1];
  Which[Harmonic === "Tesseral",
        Alm[J_, m_] := -1/(2 J +1) Total[#2 SolidTesseralHarmonicY[J, m, #1]/Norm[#1]^(J+1) &@@@ field],
        Harmonic === "Spherical",
        Alm[J_, m_] := -1/(2 J +1) Total[#2 SolidSphericalHarmonicY[J, m, #1]/Norm[#1]^(J+1) &@@@ field]];

  g = Length[grp];
  coeff = Table[
    vars = Table[Subscript[A, J, m], {m, -J, J}];
    vars = Table[A[J, m], {m, -J, J}];
    eqns = Join[Thread[vars == Chop[Simplify[1/g Total[#.vars & /@ GetAngularMomentumRep[latt, grp, J, Harmonic]]]]], {Total[vars^2] != 0}];
    sol = Flatten@Quiet@Solve[eqns, Reverse@vars];
    If[sol == {}, sol = Thread[vars -> 0]];
    vars /. sol /. A[jj_,mm_]->Alm[jj,mm], {J, 1, Jcut}];
  Return[Chop@Join[{{Alm[0,0]}}, coeff]]
]

GetCrystalFieldHam[Coeff_, l_] := Module[{ham, jcut, m1, m2, j, m},
  jcut = Length[Coeff]-1;
  ham = Sum[Table[(-1)^m1 ToExpression["r"]^j Coeff[[j+1]][[j+m+1]] ThreeYIntegral[{l, -m1}, {j, m}, {l, m2}], {m1, -l, l}, {m2, -l, l}], {j, 0, jcut}, {m, -j, j}];
  Return[ham]
]
(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
