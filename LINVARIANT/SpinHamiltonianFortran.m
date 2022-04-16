BeginPackage["LINVARIANT`SpinHamiltonianFortran`", {"LINVARIANT`Structure`", "LINVARIANT`GroupTheory`"}]

(*--------- Load, Save and Modify Crystal Structure Libraries ------------*)
FortranSpinMCStep    ::usage "FortranSpinMCStep[FieldDim]"
(*--------- Plot and Manipulate Crystal Structures -------------------- ----------------*)

(*--------- Point and Space Group Information ---------------------------*)

(*--------------------------------------------------*)
(*-------------------------- Internal --------------*)
(*--------------------------------------------------*)

(*--------------------------- Options ----------------------------*)

(*--------------------------- external Modules -------------------*)

Begin["`Private`"]

(*--------------------------- Modules ----------------------------*)
FortranSpinMCStep[FieldDim_?ListQ] := Module[{Fint, n, numfield, numcomponent, Fout},
  {numfield, numcomponent} = FieldDim;
  Fint[n_] := StringRepeat["  ", n];
  Fout = {{"Subroutine MCStep(Fields, T, initdamp)"},
    {Fint[1]},
    {Fint[1]<>"Use LINVARIANT"},
    {Fint[1]<>"Use Constants"},
    {Fint[1]<>"Use Inputs"},
    {Fint[1]},
    {Fint[1] <> "Real*8,  Intent(in)    :: T"},
    {Fint[1] <> "Real*8,  Intent(inout) :: initdamp"}, 
    {Fint[1] <> "Real*8,  Intent(inout) :: Fields(NumField, cgrid%n1, cgrid%n2, cgrid%n3, FieldDim)"},
    {Fint[1] <> "integer                :: diceX, diceY, diceZ, i, idelta, ix, iy, iz"},
    {Fint[1] <> "integer                :: seed=0, imc, idum, jdum, accepted"},
    {Fint[1] <> "real*8                 :: FlippedSpin(3),dspin(3), DeltaE"},
    {Fint[1] <> "real*8                 :: rdum, AccProb"},
    {Fint[1]},
    {Fint[1] <> "imc = 0"},
    {Fint[1] <> "accepted = 0"},
    {Fint[1]},
    {Fint[1] <> "do ix = 1, cgrid%n1"},
    {Fint[2] <> "do iy = 1, cgrid%n2"},
    {Fint[3] <> "do iz = 1, cgrid%n3"},
    {Fint[4] <> "do idelta = 1, NumField"},
    {Fint[1]},
    {Fint[5] <> "do i = 1, FieldDim"},
    {Fint[6] <> "FlippedSpin(i) = initdamp*(ran1(seed)-0.5d0) + Fields(idelta,ix,iy,iz,i)"},
    {Fint[5] <> "end do"},
    {Fint[5] <> "do i = 1, FieldDim"},
    {Fint[6] <> "dspin(i) = FlippedSpin(i)/norm2(FlippedSpin) - Fields(idelta,ix,iy,iz,i)"},
    {Fint[5] <> "end do"},
    {Fint[1]},
    {Fint[5] <> "DeltaE = DeltaHJij(ix, iy, iz, Fields, dspin, idelta)"},
    {Fint[1]},
    {Fint[5] <> "AccProb = Min(1.0d0, Exp(-1.0d0*DeltaE/(k_bolt_ev*T)))"},
    {Fint[5] <> "if(ran1(seed) .lt. AccProb) then"},
    {Fint[6] <> "accepted = accepted + 1"},
    {Fint[6] <> "do i = 1, FieldDim"},
    {Fint[7] <> "Fields(idelta,ix,iy,iz,i) = Fields(idelta,ix,iy,iz,i) + dspin(i)"},
    {Fint[6] <> "end do"},
    {Fint[5] <> "end if"},
    {Fint[1]},
    {Fint[5] <> "imc = imc + 1"},
    {Fint[1]},
    {Fint[4] <> "end do"},
    {Fint[3] <> "end do"},
    {Fint[2] <> "end do"},
    {Fint[1] <> "end do"},
    {Fint[1]},
    {Fint[1] <> "rdum = real(accepted)/real(NumField*cgrid%n1*cgrid%n2*cgrid%n3)"},
    {Fint[1] <> "if(rdum .gt. AcceptRatio) then"},
    {Fint[2] <> "initdamp = DampRatio*initdamp"},
    {Fint[1] <> "else"},
    {Fint[2] <> "initdamp = initdamp/DampRatio"},
    {Fint[1] <> "end if"},
    {Fint[1] <> "initdamp = Min(initdamp, 10.0d0)"},
    {Fint[1]},
    {Fint[1] <> "do idum = 1, FieldDim*int((cgrid%n1*cgrid%n2*cgrid%n3)**0.5d0)"},
    {Fint[2] <> "do idelta = 1, NumField"},
    {Fint[3] <> "jdum = mod(idum,3)+1"},
    {Fint[3] <> "diceX = int(ran1(seed)*cgrid%n1) + 1"},
    {Fint[3] <> "diceY = int(ran1(seed)*cgrid%n2) + 1"},
    {Fint[3] <> "diceZ = int(ran1(seed)*cgrid%n3) + 1"},
    {Fint[3] <> "do i = 1, FieldDim"},
    {Fint[4] <> "dspin(i) = 0.0d0"},
    {Fint[3] <> "end do"},
    {Fint[3] <> "dspin(jdum) = -2.0d0*Fields(idelta,diceX,diceY,diceZ,jdum)"},
    {Fint[1]},
    {Fint[3] <> "DeltaE = DeltaHJij(diceX, diceY, diceZ, Fields, dspin, idelta)"},
    {Fint[1]},
    {Fint[3] <> "AccProb = Min(1.0d0, Exp(-1.0d0*DeltaE/(k_bolt_ev*T)))"},
    {Fint[3] <> "if(ran1(seed) .lt. AccProb) then"},
    {Fint[4] <> "Fields(idelta,diceX,diceY,diceZ,jdum) = Fields(idelta,diceX,diceY,diceZ,jdum) + dspin(jdum)"},
    {Fint[3] <> "end if"},
    {Fint[2] <> "end do"},
    {Fint[1] <> "end do"},
    {"End Subroutine MCStep"}};
  Return[Fout]
]

(*-------------------------- Attributes ------------------------------*)

(*Attributes[]={Protected, ReadProtected}*)

End[]

EndPackage[]
