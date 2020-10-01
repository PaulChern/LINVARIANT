BeginPackage["LINVARIANT`SphericalHarmonics`"]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
Mat2EulerVector            ::usage "Mat2EulerVector[mat]"
EulerVector2Mat            ::usage "EulerVector2Mat[axis, ang, \[Epsilon]]"
Mat2EulerAngles            ::usage "Mat2EulerAngles[latt, mat]"
SolidSphericalHarmonicY    ::usage "SolidSphericalHarmonicY[l, m, x1, x2, x3, coord]"
SolidTesseralHarmonicY     ::usage "SolidTesseralHarmonicY[l, m, x1, x2, x3, coord]"
TesseralLabel              ::usage "TesseralLabel[l, i]"
GetEulerRodrigues          ::usage "GetEulerRodrigues[n, \[Theta], \[Epsilon]]"
ThreeYIntegral             ::usage "ThreeYIntegral[jm1, jm2, jm3]"

(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
EulerVector2Mat[latt_, axis_, ang_, \[Epsilon]_] := Module[{i, j, k, n, R},
  n = Normalize[latt.axis];
  R = Sum[Table[Cos[ang] KroneckerDelta[i, j] + (\[Epsilon] - Cos[ang]) n[[i]] n[[j]] - Sin[ang] LeviCivitaTensor[3][[i, j, k]] n[[k]], {i, 3}, {j, 3}], {k, 3}];
  Return[Rationalize[Inverse[latt].R.latt]]
]

Mat2EulerVector[mat_] := Module[{\[Phi], axes, es, s, \[Epsilon]},
  \[Epsilon] = Rationalize@Det[mat];

  \[Phi] = ArcCos[1/2 (Tr[\[Epsilon] mat] - 1)];
  If[\[Phi] == 0., Return[{{0, 0, 0}, 0, \[Epsilon]}]];
  s = -Sign@N[Chop[Im[Det[Eigensystem[\[Epsilon] mat][[2]]]]]];
  \[Phi] = If[s==0., \[Phi], s \[Phi]];
  es = Eigensystem[Rationalize[\[Epsilon] mat]];
  axes = es[[2]][[First@First@Position[es[[1]], 1]]];

  Return[{axes, Rationalize[\[Phi]/Pi]*Pi, \[Epsilon]}]
]

GetEulerRodrigues[n_, \[Theta]_, \[Epsilon]_:1] := Module[{K, R, nx, ny, nz},

  If[MatrixQ[n], K = n, 
    {nx, ny, nz} = n;
    K = {{0, -nz, ny}, {nz, 0, -nx}, {-ny, nx, 0}}
    ];

  If[\[Theta] == 0., 
    R = IdentityMatrix[3],
    R = Expand[(IdentityMatrix[3] + Sin[\[Theta]] K + (1 - Cos[\[Theta]]) K.K)]
    ];
  Return[\[Epsilon] IdentityMatrix[3].R]
]

Mat2EulerAngles[latt_, mat_] := Module[{\[Epsilon], axis, \[Phi], R, K, nx, ny, nz, \[Alpha], \[Beta], \[Gamma]},
  {axis, \[Phi], \[Epsilon]} = Mat2EulerVector[mat];
  If[\[Phi] == 0., Return[{{0, 0, 0}, \[Epsilon]}]];
  {nx, ny, nz} = Normalize[latt\[Transpose].axis];
  K = {{0, -nz, ny}, {nz, 0, -nx}, {-ny, nx, 0}};
  (*R = Expand[(\[Epsilon] IdentityMatrix[3] + Sin[\[Phi]] K + \[Epsilon] (1 - \[Epsilon] Cos[\[Phi]]) K.K)];*)
  R = Expand[(IdentityMatrix[3] + Sin[-\[Phi]] K + (1 - Cos[-\[Phi]]) K.K)];
  {\[Alpha], \[Beta], \[Gamma]} = EulerAngles[R];
  (*\[Beta] = ArcCos[R[[3, 3]]];
  Which[
    \[Beta] == 0.,
      \[Alpha] = 0; \[Gamma] = \[Phi],
    \[Beta] == Pi,
      \[Alpha] = 0; \[Gamma] = Pi - \[Phi],
    True,
       \[Alpha] = ArcTan[nz nx Sin[\[Phi]/2] + ny Cos[\[Phi]/2], nz ny Sin[\[Phi]/2] - nx Cos[\[Phi]/2]];
       \[Gamma] = ArcTan[-nz nx Sin[\[Phi]/2] + ny Cos[\[Phi]/2], nz ny Sin[\[Phi]/2] + nx Cos[\[Phi]/2]];
      ];*)
  Return[{Chop[{Rationalize[\[Alpha]/Pi]*Pi, Rationalize[\[Beta]/Pi]*Pi, Rationalize[\[Gamma]/Pi]*Pi}], \[Epsilon]}]
]

SolidSphericalHarmonicY[l_?IntegerQ, m_?IntegerQ, xyz_, coord_: "Cartesian"] := Module[{rr, \[Theta]\[Theta], \[Phi]\[Phi], xx, yy, zz, x1, x2, x3},
  {x1, x2, x3} = xyz;
  Which[coord == "Cartesian", 
                 FullSimplify@Evaluate[TransformedField["Spherical" -> "Cartesian", rr^l SphericalHarmonicY[l, m, \[Theta]\[Theta], \[Phi]\[Phi]], {rr, \[Theta]\[Theta], \[Phi]\[Phi]} -> {xx, yy, zz}]] /. {xx -> x1, yy -> x2, zz -> x3},
        coord == "Polar", 
                 FullSimplify@Evaluate[rr^l SphericalHarmonicY[l, m, \[Theta]\[Theta], \[Phi]\[Phi]]] /. {rr -> x1, \[Theta]\[Theta] -> x2, \[Phi]\[Phi] -> x3}, 
        True, 
        Print["Wrong coordinate"]]
]

SolidTesseralHarmonicY[l_?IntegerQ, m_?IntegerQ, xyz_, coord_: "Cartesian"] := Module[{x1, x2, x3}, 
  {x1, x2, x3} = xyz;
  Simplify@Which[
                 m > 0, 
                 Sqrt[2]/2 (-1)^m (SolidSphericalHarmonicY[l, m, xyz, coord] + (-1)^m SolidSphericalHarmonicY[l, -m, xyz, coord]), 
                 m < 0, 
                 Sqrt[2]/(2 I ) (-1)^m (SolidSphericalHarmonicY[l, Abs[m], xyz, coord] - (-1)^m SolidSphericalHarmonicY[l, -Abs[m], xyz, coord]), 
                 m == 0, 
                 SolidSphericalHarmonicY[l, 0, xyz, coord]
                ]
]

TesseralLabel[l_, i_] := Module[{dict},
  dict = <|1 -> {"py", "pz", "px"},
           2 -> {"dxy", "dyz", "dz2", "dxz", "dx2-y2"}|>;
  Return[dict[l][[i]]]
]

ThreeYIntegral[jm1_, jm2_, jm3_] := Module[{j1, j2, j3, m1, m2, m3},
  {j1, m1} = jm1;
  {j2, m2} = jm2;
  {j3, m3} = jm3;
  Quiet[Sqrt[(2 j1 + 1)*(2 j2 + 1)*(2 j3 + 1)/(4 \[Pi])] ThreeJSymbol[{j1, 0}, {j2, 0}, {j3,0}] ThreeJSymbol[{j1, m1}, {j2, m2}, {j3, m3}]]
]
(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
