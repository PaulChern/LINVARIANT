BeginPackage["LINVARIANT`Mc`"]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
MCStep               ::usage "MCStep[initconfig, T, initdamp, MeshDim]"
(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
MMAEvalOnSiteEnergy[SiteIndex_, Config_] := Module[{sub, ii, x, y, z, xx, yy, zz},
   {x, y, z} = SiteIndex;
   sub = Dispatch[Flatten[Table[{Subscript[sm, ii, x, y, z] -> Config[[1, ii, z, y, x]], Subscript[hf, ii, x, y, z] -> Config[[2, ii, z, y, x]]}, {ii, 3}]]];
   Return[ESiteTable[[z, y, x]] /. sub]
];

MCStep[initconfig_, T_, initdamp_, MeshDim_] := Module[{Lx, Ly, Lz, AcceptRatio, config, diceX, diceY, diceZ, ii, ix, iy, iz, imc, idum, jdum, rdum, accepted, damp, dmode, FlippedConfig, \[CapitalDelta]E, AcceptanceProbability, KB = 8.617 10^-5, Hartree = 27.211386245, debug},
  FlippedConfig = initconfig;
  imc = 0;
  accepted = 0;
  AcceptRatio = 0.2;
  {Lx, Ly, Lz} = MeshDim;
  damp = initdamp;
  Do[config = FlippedConfig;
   
     dmode = damp RandomReal[{-0.5, 0.5}, 3];
     FlippedConfig[[ii, All, iz, iy, ix]] += dmode;
   
     \[CapitalDelta]E = MMAEvalOnSiteEnergy[{ix, iy, iz}, FlippedConfig] - MMAEvalOnSiteEnergy[{ix, iy, iz}, config];
   
     Which[\[CapitalDelta]E < 45.0, 
           AcceptanceProbability = Min[1, Exp[-((\[CapitalDelta]E Hartree)/(KB T)) ]];
           If[RandomReal[] < AcceptanceProbability, accepted += 1, FlippedConfig = config],
           True,
           FlippedConfig = config];
   
     imc += 1, {iz, Lz}, {iy, Ly}, {ix, Lx}, {ii, 2}
    ];
  
  rdum = accepted/(2 Lx Ly Lz);
  If[rdum > AcceptRatio, damp = 1.01 damp, damp/1.01];
  damp = Min[damp, 0.3];
  
  Return[{damp, FlippedConfig, rdum}]
]

(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
