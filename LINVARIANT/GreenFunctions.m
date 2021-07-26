BeginPackage["LINVARIANT`GreenFunctions`", {"LINVARIANT`TBHamiltonian`", "LINVARIANT`Structure`", "LINVARIANT`MathematicaPlus`", "LINVARIANT`QM`"}]

H2G                ::usage "H2G[H, e, f]"
TB2G               ::usage "TB2G[Ti0, pos, k, e, f]"
TB2Delta           ::usage "TB2Delta[Ti0up, Ti0dn, pos, kbz]"
TB2GR              ::usage "TB2GR[Ti0, pos, e, fermi, R, kbz]"
VGVG               ::usage "VGVG[Ti0, spg, pos, orb, e, efermi, R, kmesh]"

Begin["`Private`"]

H2G[H_, e_, fermi_] := Module[{dim},
  dim = Length[H];
  Chop@Inverse[(e + fermi) IdentityMatrix[dim] - H]
]

TB2G[Ti0_, pos_, k_, e_, fermi_] := Module[{dim, Hk},
  Hk = TBHk[Ti0, pos, k];
  dim = Length[Hk];
  Chop@Inverse[(e + fermi) IdentityMatrix[dim] - Hk]
]

TB2GR[Ti0_, pos_, e_, fermi_, R_, kbz_] := Module[{k, Hk, Gk},
  Sum[Hk = TBHk[Ti0, pos, k[[1]]];
      Gk = H2G[Hk, e, fermi];
      k[[2]] Gk Exp[-I 2 Pi k[[1]].R], {k, kbz}]
]

TB2Delta[Ti0up_, Ti0dn_, pos_, kbz_] := Module[{k, Hkup, Hkdn},
  Chop@Sum[Hkup = TBHk[Ti0up, pos, k[[1]]];
           Hkdn = TBHk[Ti0dn, pos, k[[1]]];
           k[[2]] (Hkup - Hkdn), {k, kbz}]
]

VGVG[Ti0_, pos_, orb_, e_, efermi_, R_, kmesh_] := Module[{Ti0up, Ti0dn, T00, orbsoc, GRij, GRji, GRup, GRdn, Delta, Aij, kbz, i, j, \[Alpha], \[Beta]},
  kbz = GetBZKList[None, kmesh];
  Aij = If[Length[Ti0] == 2,
           (* Collinear *)
           {Ti0up, Ti0dn} = Ti0;
           Delta = Mat2SiteBlock[TB2Delta[Ti0up, Ti0dn, pos, kbz], orb];
           GRup = Mat2SiteBlock[TB2GR[Ti0up, pos, e, efermi, R, kbz], orb];
           GRdn = Mat2SiteBlock[TB2GR[Ti0dn, pos, e, efermi, -R, kbz], orb];
           Table[1/(4 Pi) Tr[(Delta[[i, i]].GRup[[i, j]]).(Delta[[j, j]].GRdn[[j, i]])],
                 {i, Length@orb}, {j, Length@orb}],
           (* Noncollinear *)
           T00 = If[#[[2]] == {0, 0, 0}, #, ## &[]] & /@ Ti0;
           orbsoc = {#1, 2 #2} & @@@ orb;
           Delta = Mat2SiteBlock[PauliNormDecompose[TBHk[T00, pos, {0, 0, 0}]], orb];
           GRij = Mat2SiteBlock[TB2GR[Ti0, pos, e, efermi, R, kbz], orbsoc];
           GRji = Mat2SiteBlock[TB2GR[Ti0, pos, e, efermi, -R, kbz], orbsoc];
           Table[Table[1/Pi Tr[(Delta[[i, i]].PauliDecompose[GRij[[i, j]]][[\[Alpha]]]).(Delta[[j, j]].PauliDecompose[GRji[[j, i]]][[\[Beta]]])],
                 {i, Length@orbsoc}, {j, Length@orbsoc}], {\[Alpha], 4}, {\[Beta], 4}]];
  Return[Aij]
]

End[]

EndPackage[]
