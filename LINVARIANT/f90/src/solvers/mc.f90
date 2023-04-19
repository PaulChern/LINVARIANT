Subroutine MCMCStep(imc, Fields, EwaldField, e0ij, T, udamp, etadamp, dene)
  
  Use LINVARIANT
  Use Constants
  Use Parameters
  Use Variation
  
  Implicit none
  Real*8,  Intent(in)    :: T
  Integer, Intent(inout) :: imc
  Real*8,  Intent(inout) :: udamp(NumField), etadamp, dene
  Real*8,  Intent(inout) :: Fields(FieldDim, NumField, cgrid_a%n1+cgrid_b%n1, cgrid_a%n2+cgrid_b%n2, cgrid_a%n3+cgrid_b%n3)
  Real*8,  Intent(inout) :: EwaldField(3, NumField, cgrid_a%n1+cgrid_b%n1, cgrid_a%n2+cgrid_b%n2, cgrid_a%n3+cgrid_b%n3)
  Real*8,  Intent(inout) :: e0ij(3,3)
  Integer                :: i, j, idelta, ix, iy, iz
  Integer                :: idum, acceptedu(NumField), acceptedeta
  Real*8                 :: dfield(FieldDim), deta(6), DeltaE
  Real*8                 :: ProbTest
  Real*8                 :: rdum, AccProb
  
  acceptedu = 0
  acceptedeta = 0
  DeltaE = 0.0D0
  dene = 0.0D0
  dfield = 0.0D0
 
  do iz = 1, cgrid%n3
    do iy = 1, cgrid%n2
      do ix = 1, cgrid%n1
        do idelta = 1, NumField
          If (.not.FrozenQ(idelta)) then
          
            Call random_number(dfield)
            dfield = udamp(idelta)*(dfield - 0.5D0)
            
            DeltaE = GetVariation(ix, iy, iz, Fields, e0ij, idelta, dfield)
            if (DipoleQ.and.(idelta.lt.NumIRFields)) then
              DeltaE = DeltaE &
                + GetVariationEwald(ix, iy, iz, EwaldField, idelta, dfield)
            end if
  
            AccProb = Min(1.0d0, Exp(-1.0d0*DeltaE/T/k_bolt_ev*Hartree))
            Call random_number(ProbTest)
            if(ProbTest .lt. AccProb) then
              acceptedu(idelta) = acceptedu(idelta) + 1
              dene = dene + DeltaE
              do i = 1, FieldDim
                Fields(i,idelta,ix,iy,iz) = Fields(i,idelta,ix,iy,iz) + dfield(i)*FieldsBinary(i,idelta)
              end do
              if (DipoleQ.and.(idelta.lt.NumIRFields)) then
                Call UpdateEwaldField(ix, iy, iz, EwaldField, idelta, dfield)
              end if
            end if
          end if 
        end do
      end do
    end do
  end do
  
  Call RemoveGridDrifts(Fields)
 
  if(.Not.All(CLAMPQ)) then
    do idum = 1, 2*cgrid%n3+1
      Call random_number(deta)
      deta = etadamp*(deta - 0.5D0)
      do i = 1, 6
        If (CLAMPQ(i)) deta(i) = 0.0D0
      end do

      DeltaE = 0.0D0
      do iz = 1, cgrid%n3
        do iy = 1, cgrid%n2
          do ix = 1, cgrid%n1
            DeltaE = DeltaE + GetVariation(ix, iy, iz, Fields, e0ij, NumField+1, deta)
!            DeltaE = DeltaE &
!              + GetVariationGpscoupling(ix, iy, iz, Fields, e0ij, NumField+1, deta) &
!              + GetVariationHstrain(ix, iy, iz, Fields, e0ij, NumField+1, deta)
          end do
        end do
      end do
!      DeltaE = DeltaE + cgrid%npts*GetVariationGstrain(1, 1, 1, Fields, e0ij, NumField+1, deta)
      
      AccProb = Min(1.0d0, Exp(-1.0d0*DeltaE/T/k_bolt_ev*Hartree))
      Call random_number(ProbTest)
      if(ProbTest .lt. AccProb) then
        acceptedeta = acceptedeta + 1
        dene = dene + DeltaE
        e0ij(1,1) = e0ij(1,1) + deta(1)
        e0ij(2,2) = e0ij(2,2) + deta(2)
        e0ij(3,3) = e0ij(3,3) + deta(3)
        e0ij(2,3) = e0ij(2,3) + deta(4)
        e0ij(1,3) = e0ij(1,3) + deta(5)
        e0ij(1,2) = e0ij(1,2) + deta(6)
      end if
      
    end do
  end if
  
  do idelta = 1, NumField
    If (.not.FrozenQ(idelta)) then
      rdum = real(acceptedu(idelta))/real(cgrid%npts)
      if(rdum .gt. AcceptRatio) then
        udamp(idelta) = DampRatio*udamp(idelta)
      else
        udamp(idelta) = udamp(idelta)/DampRatio
      end if
      udamp(idelta) = Min(udamp(idelta), 3.0d0)
    end if
  end do
  
  if(.Not.All(CLAMPQ)) then
    if(mod(imc,1).eq.0) then
      rdum = real(acceptedeta)/real(2*cgrid%n3+1)
      if(rdum .gt. AcceptRatio) then
        etadamp = DampRatio*etadamp
      else
        etadamp = etadamp/DampRatio
      end if
      etadamp = Min(etadamp, 3.0d0)
    end if
  end if
End Subroutine MCMCStep

Subroutine WLMCStep(imc, Fields, e0ij, udamp, etadamp, wl_s, wl_h, wl_f)
  Use LINVARIANT
  Use Constants
  Use Parameters
  Use Variation

  Implicit none
  Integer, Intent(in)    :: imc
  Real*8,  Intent(inout) :: udamp(NumField), etadamp
  Real*8,  Intent(inout) :: Fields(FieldDim, NumField, cgrid_a%n1+cgrid_b%n1, cgrid_a%n2+cgrid_b%n2, cgrid_a%n3+cgrid_b%n3)
  Real*8,  Intent(inout) :: e0ij(3,3), wl_f, wl_s(2,ContourNPoints(2)), wl_h(2,ContourNPoints(2))
  Real*8                 :: Fields_tmp(FieldDim, NumField, cgrid_a%n1+cgrid_b%n1, cgrid_a%n2+cgrid_b%n2, cgrid_a%n3+cgrid_b%n3)
  Real*8                 :: e0ij_tmp(3,3)
  Integer                :: i, j, idelta, ix, iy, iz, Ei_bin, Ef_bin
  Integer                :: idum, acceptedu(NumField), acceptedeta
  Real*8                 :: dfield(FieldDim), deta(6), Ei, Ef
  Real*8                 :: ProbTest
  Real*8                 :: rdum, AccProb

  acceptedu = 0
  acceptedeta = 0
  dfield = 0.0D0

  Ei = GetEtot(Fields, e0ij)
  Ei_bin = Energy2bin(Ei,wl_s)

  do iz = 1, cgrid%n3
    do iy = 1, cgrid%n2
      do ix = 1, cgrid%n1
        do idelta = 1, NumField

          Call random_number(dfield)
          dfield = udamp(idelta)*(dfield - 0.5D0)
          do i = 1, FieldDim
            Fields_tmp(i,idelta,ix,iy,iz) = Fields(i,idelta,ix,iy,iz) + FieldsBinary(i,idelta)*dfield(i)
          end do
          Ef = GetEtot(Fields_tmp, e0ij)
          Ef_bin = Energy2bin(Ef,wl_s)
          AccProb = Min(1.0d0, Exp(wl_s(2,Ei_bin)-wl_s(2,Ef_bin)))
          if(ProbTest .lt. AccProb) then
            acceptedu(idelta) = acceptedu(idelta) + 1
            wl_h(2,Ef_bin) = wl_h(2,Ef_bin) + 1
            wl_s(2,Ef_bin) = wl_s(2,Ef_bin) + Exp(wl_f)
            Ei = Ef
            Ei_bin = Ef_bin
            Fields = Fields_tmp
          else 
            wl_h(2,Ei_bin) = wl_h(2,Ei_bin) + 1
            wl_s(2,Ei_bin) = wl_s(2,Ei_bin) + Exp(wl_f)
          end if
           
        end do
      end do
    end do
  end do

  do idum = 1, 2*cgrid%n3+1
    Call random_number(deta)
    deta = etadamp*(deta - 0.5D0)
    e0ij_tmp(1,1) = e0ij(1,1) + deta(1)
    e0ij_tmp(2,2) = e0ij(2,2) + deta(2)
    e0ij_tmp(3,3) = e0ij(3,3) + deta(3)
    e0ij_tmp(2,3) = e0ij(2,3) + deta(4)
    e0ij_tmp(1,3) = e0ij(1,3) + deta(5)
    e0ij_tmp(1,2) = e0ij(1,2) + deta(6)
    Ef = GetEtot(Fields, e0ij_tmp)
    Ef_bin = Energy2bin(Ef,wl_s)
    AccProb = Min(1.0d0, Exp(wl_s(2,Ei_bin)-wl_s(2,Ef_bin)))

    Call random_number(ProbTest)
    if(ProbTest .lt. AccProb) then
      acceptedeta = acceptedeta + 1
      wl_h(2,Ef_bin) = wl_h(2,Ef_bin) + 1
      wl_s(2,Ef_bin) = wl_s(2,Ef_bin) + Exp(wl_f)
      Ei = Ef
      Ei_bin = Ef_bin
      e0ij = e0ij_tmp
    else
      wl_h(2,Ei_bin) = wl_h(2,Ei_bin) + 1
      wl_s(2,Ei_bin) = wl_s(2,Ei_bin) + Exp(wl_f)
    end if

  end do

  do idelta = 1, NumField
    rdum = real(acceptedu(idelta))/real(cgrid%npts)
    if(rdum .gt. AcceptRatio) then
      udamp(idelta) = DampRatio*udamp(idelta)
    else
      udamp(idelta) = udamp(idelta)/DampRatio
    end if
    udamp(idelta) = Min(udamp(idelta), 3.0d0)
  end do

  if(mod(imc,1).eq.0) then
    rdum = real(acceptedeta)/real(2*cgrid%n3+1)
    if(rdum .gt. AcceptRatio) then
      etadamp = DampRatio*etadamp
    else
      etadamp = etadamp/DampRatio
    end if
    etadamp = Min(etadamp, 3.0d0)
  end if

  return
  contains
  function Energy2bin(ene,wl_s) result(bin)
    implicit none
    Real(dp),intent(inout)  :: ene
    Real(dp),intent(in)     :: wl_s(2,ContourNPoints(2))
    Integer                 :: i, bin

    if(ene.lt.ContourMin) then
      bin = 1
    else if(ene.lt.ContourMax) then
      bin = ContourNPoints(2)
    else
      do i = 1, ContourNPoints(2)
        if(abs(wl_s(1,i) - ene).le.(ContourMax-ContourMin)/(ContourNPoints(2)-1)/2) bin = i
      end do
    end if
    
  End function

End Subroutine WLMCStep
