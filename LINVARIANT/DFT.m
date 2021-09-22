BeginPackage["LINVARIANT`DFT`"]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
DFTscf              ::usage "DFTscf[potF, domain, MaxIterations]"
SchrodingerSolver   ::usage "SchrodingerSolver[var, Vpot, NBands]"
PoissonSolver       ::usage "PoissonSolver[var, \[Rho]]"
BandsPlot           ::usage "BandsPlot[bsxml, klabels]"
Kpoints2Kpath       ::usage "Kpoints2Kpath[kp, latt]"
(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
SchrodingerSolver[var_, Vpot_Function, NBands_] := Module[{eqn, bcs, En, shift, Psi, \[Psi], r, rmin, rmax},
  {r, rmin, rmax} = var;
  shift = 1000000;
  eqn = {shift*\[Psi][r] - 1/2 Laplacian[\[Psi][r], {r}] + Vpot[r] \[Psi][r]};
  bcs = {DirichletCondition[\[Psi][r] == 0, True]};
  {En, Psi} = NDEigensystem[Join[eqn, bcs], \[Psi][r], {r, rmin, rmax}, NBands];
  Return[{# - shift & /@ En, Head[#] & /@ Psi}\[Transpose]]
]

PoissonSolver[var_, \[Rho]_Function] := Module[{eqn, bcs, Phi, \[Phi], r, rmin, rmax},
  {r, rmax} = var;
  eqn = {D[\[Phi][r], {r, 2}] == -\[Rho][r]/r};
  bcs = {DirichletCondition[\[Phi][r] == 0, r == 0],
         DirichletCondition[\[Phi][r] == 1, r == rmax]};
  Phi = NDSolveValue[Join[eqn, bcs], \[Phi][r], {r, 0, rmax}];
  Return[Head@Phi]
]

DFTscf[potF_Function, domain_, MaxIterations_] := Module[{r, \[Psi], \[Phi], scf, En, rmin, rmax, VHartree, Vxc, Vtot, its = 0, energy, xcEnergy, hartreeEnergy},
  {rmin, rmax} = domain;
  VHartree = Function[{U, r}, 2. U[r]/r];
  Vxc = Function[{u, r}, -(3./Pi*2.*u[r]^2/((r^2)*4*Pi))^(1./3)];
  {En, \[Psi]} = SchrodingerSolver[{r, rmin, rmax}, potF, 1][[1]];
  \[Phi] = Quiet@PoissonSolver[{r, rmax}, Function[r, \[Psi][r]\[Conjugate] \[Psi][r]]];

  FixedPoint[Function[its++;
               Vtot = Function[r, Function[{r, \[Psi], \[Phi]}, potF[r] + VHartree[\[Phi], r] + Vxc[\[Psi], r]][r, \[Psi], \[Phi]]];
               {En, \[Psi]} = SchrodingerSolver[{r, rmin, rmax}, Vtot, 1][[1]];
               \[Phi] = Quiet@PoissonSolver[{r, rmax}, Function[r, \[Psi][r]\[Conjugate] \[Psi][r]]];
               Echo[En, its]], 
             0, MaxIterations, SameTest -> (Abs[# - #2] < .001 &)];

  hartreeEnergy = NIntegrate[VHartree[\[Phi], r]*Function[r, \[Psi][r]\[Conjugate] \[Psi][r]][r], {r, rmin, rmax}];
  xcEnergy = NIntegrate[Vxc[\[Psi], r]*Function[r, \[Psi][r]\[Conjugate] \[Psi][r]][r], {r, rmin, rmax}];
  energy = 2 En - (hartreeEnergy + 1/2 xcEnergy);
  scf = <|"Energy" -> energy, "Energies" -> <|"RadialEnergy" -> En, "ExchangeCorrelation" -> xcEnergy, "Hartree" -> hartreeEnergy|>, "Wavefunction" -> \[Psi], "Potential" -> \[Phi], "Iterations" -> its|>;
  Return[scf]
]

BandsPlot[bsxml_, hkpts_, OptionsPattern[{"PlotRange" -> {All, All}, "AspectRatio"->1/GoldenRatio}]] := Module[{latt, blatt, kpath, updata, dndata, upplot, dnplot, klabels},
  latt = bsxml["latt"];
  blatt = 1/(2*Pi)*Inverse[latt]\[Transpose];
  klabels = Kpoints2Kpath[hkpts,blatt];
  kpath = {Accumulate[Join[{0}, Norm[blatt.#] & /@ Differences[bsxml["k"]]]], bsxml["k"]}\[Transpose];
  updata = {kpath\[Transpose][[1]], #\[Transpose][[1]]}\[Transpose] & /@ (bsxml["up"]\[Transpose]);
  upplot = ListLinePlot[updata,
                        PlotStyle -> {{Black, Thick}},
                        Joined -> True,
                        PlotRange -> OptionValue["PlotRange"],
                        AspectRatio -> OptionValue["AspectRatio"],
                        Frame -> True,
                        GridLines -> {{klabels\[Transpose][[2]], ConstantArray[Thick, Length[klabels]]}\[Transpose], Automatic},
                        FrameTicks -> {{Automatic, None}, {klabels[[;;,2;;3]], None}},
                        ImageSize -> Medium];
  If[MemberQ[Keys[bsxml], "dn"],
     dndata = {kpath\[Transpose][[1]], #\[Transpose][[1]]}\[Transpose] & /@ (bsxml["dn"]\[Transpose]);
     dnplot = ListLinePlot[dndata,
                           PlotStyle -> {{Black, Dashed, Thin}},
                           Joined -> True,
                           PlotRange -> OptionValue["PlotRange"],
                           AspectRatio -> OptionValue["FigRatio"],
                           Frame -> True,
                           GridLines -> {{klabels\[Transpose][[1]], ConstantArray[Thick, Length[klabels]]}\[Transpose], Automatic},
                           FrameTicks -> {{Automatic, None}, {klabels, None}},
                           ImageSize -> Medium];
     Return[Show[{upplot,dnplot}]],
     Return[upplot]]
]

Kpoints2Kpath[kp_, latt_] := Module[{},
  Which[Depth[kp]==3,
        {kp, Accumulate[Join[{0}, Norm[latt\[Transpose].#] & /@ Differences[kp]]]}\[Transpose],
        Depth[kp]==4,
        Join[Kpoints2Kpath[kp\[Transpose][[1]], latt]\[Transpose],{kp\[Transpose][[2]]}]\[Transpose]
        ]
]
(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
