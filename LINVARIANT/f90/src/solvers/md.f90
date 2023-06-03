Subroutine TikTok(Fields, dFieldsdt, e0ij, de0ijdt, T0, gm, TimeNow)
  
  Use omp_lib
  Use LINVARIANT
  Use Force
  Use Parameters
  Use Constants
  Use Aux
  
  Implicit none
  Real*8,  Intent(in)    :: T0
  Real*8,  Intent(inout) :: gm(NumField+1)
  Real*8,  Intent(inout) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(inout) :: dFieldsdt(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(inout) :: de0ijdt(3,3)
  Real*8,  Intent(inout) :: e0ij(3,3), TimeNow
  Real*8                 :: Fields0(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8                 :: e0ij0(3,3), mu, sigma, tau(3), t
  Real*8                 :: Forces(Max(FieldDim, 6),NumField+1,cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8                 :: Friction(FieldDim,NumField,cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8                 :: ForcesEta(6), FrictionEta(6),VolumeForce(3,3)
  Integer                :: i, ifield, fi, ix, iy, iz

  Fields0 = Fields
  e0ij0 = e0ij

  t = (TimeNow-ThermoSteps*DeltaT)*time_fs
  sigma = 3*EPhi(1)/sqrt(2*DLog(2.0_dp))/2
  mu = 6*sigma
  tau = 1.0D0
!  tau(1) = 1.0-(1.0-1.0D1)*Exp(-(t-mu)**2/(2*sigma**2))

  Forces = 0.0D0
  If(DipoleQ) then
    Forces(1:3,:,:,:,:) = GetForcesEwald(Fields0)
  End If
  If(EfieldQ) then
    Call ApplyEfield(TimeNow, Forces)
  End If
  Forces = Forces + GetForces(Fields0,e0ij0) 
  Do i = 1, NumField
    Forces(:,i,:,:,:) = Forces(:,i,:,:,:)/mass(i)/tau(i)/alat**2
    If(FrozenQ(i).and.(TimeNow.gt.ThermoSteps*DeltaT)) Forces(:,i,:,:,:) = 0.0D0
  End do
  Do i = 1, 6
    If(CLAMPQ(i)) then
      ForcesEta(i) = 0.0D0
    else
      ForcesEta(i) = Sum(Forces(i,NumField+1,:,:,:))/cgrid%npts/mass(NumField+1)
    End If
  End do
!  Call GetFriction(gm, dFieldsdt, de0ijdt, Friction)
  Call BaroState(e0ij, VolumeForce)

  ! v(t+dt/2)
  Do i = 1, NumField
    dFieldsdt(:,i,:,:,:) = dFieldsdt(:,i,:,:,:) &
                         + 0.5d0*DeltaT*(Forces(:,i,:,:,:)-gm(i)*dFieldsdt(:,i,:,:,:))
  End do
  de0ijdt = de0ijdt + 0.5d0*DeltaT*(eta2eij(ForcesEta)+VolumeForce-gm(NumField+1)*de0ijdt)
  ! r(t+dt)
  Fields = Fields0+dFieldsdt*DeltaT
  e0ij = e0ij + de0ijdt*DeltaT

  Call RemoveGridDrifts(Fields)
  Call RemoveGridDrifts(dFieldsdt)
  Call ReservoirUpdate(gm, T0, dFieldsdt, de0ijdt, TimeNow)

  TimeNow = TimeNow + DeltaT/2

  Forces = 0.0D0
  If(DipoleQ) then
    Forces(1:3,:,:,:,:) = GetForcesEwald(Fields)
  End If
  If(EfieldQ) then
    Call ApplyEfield(TimeNow, Forces)
  End If
  Forces = Forces + GetForces(Fields,e0ij)
  Do i = 1, NumField
    Forces(:,i,:,:,:) = Forces(:,i,:,:,:)/mass(i)/tau(i)/alat**2
    If(FrozenQ(i).and.(TimeNow.gt.ThermoSteps*DeltaT)) Forces(:,i,:,:,:) = 0.0D0
  End do

  Do i = 1, 6
    If(CLAMPQ(i)) then
      ForcesEta(i) = 0.0D0
    else
      ForcesEta(i) = Sum(Forces(i,NumField+1,:,:,:))/cgrid%npts/mass(NumField+1)
    End If
  End do
  Call BaroState(e0ij, VolumeForce)

  ! v(t+dt/2+dt/2)
  Do i = 1, NumField
!    dFieldsdt(:,i,:,:,:) = (dFieldsdt(:,i,:,:,:)+0.5d0*DeltaT*Forces(:,i,:,:,:))/(1.0D0 + 0.5D0*gm(i)*DeltaT)
    dFieldsdt(:,i,:,:,:) = dFieldsdt(:,i,:,:,:)+0.5d0*DeltaT*(Forces(:,i,:,:,:)-gm(i)*dFieldsdt(:,i,:,:,:))
  End do
  de0ijdt = de0ijdt + 0.5d0*DeltaT*(eta2eij(ForcesEta)+VolumeForce-gm(NumField+1)*de0ijdt)

  TimeNow = TimeNow + DeltaT/2

  Call ImposeBCS(Fields)

End Subroutine TikTok

Subroutine GetFriction(gm, dFieldsdt, Friction)
  Use Parameters
  Implicit none
  Real*8,  Intent(in)     :: dFieldsdt(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(in)     :: gm(NumField+1)
  Real*8,  Intent(out)    :: Friction(FieldDim,NumField,cgrid%n1, cgrid%n2, cgrid%n3)
  Integer                 :: i
  
  Do i = 1, NumField
    Friction(:,i,:,:,:) = -gm(i)*dFieldsdt(:,i,:,:,:)
  End do

End Subroutine GetFriction

Subroutine ReservoirUpdate(gm, T0, dFieldsdt, de0ijdt, TimeNow)
  Use Energy
  Use Parameters
  Use Fileparser

  Implicit none
  Real*8,  Intent(in)     :: dFieldsdt(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(in)     :: de0ijdt(3,3), T0, TimeNow
  Real*8,  Intent(inout)  :: gm(NumField+1)
  Real*8                  :: Ekin(NumField+1), mu, sigma, tau, t, tt
  Real*8                  :: collision
  Integer                 :: TotalDim, i

  t = (TimeNow-ThermoSteps*DeltaT)*time_fs
  tt = t
  if(TrainQ) tt = Mod(TimeNow/DeltaT-ThermoSteps,Float(TrainRate))*DeltaT*time_fs
  sigma = ReservoirTau*EPhi(1)/sqrt(2*DLog(2.0_dp))/2
  mu    = 6*sigma
  tau   = 1.0-(1.0-ReservoirMass)*Exp(-(tt-mu)**2/(2*sigma**2)) 

  Ekin = ResolveEkin(dFieldsdt, de0ijdt)

  call caps2small(ThermoState)
  If(.not.ReservoirQ(1).and.TimeNow.gt.ThermoSteps*DeltaT.or.(mod(nint(TimeNow/DeltaT),ReservoirRate).ne.0)) then
    gm = 0.0D0
  else if(ThermoState.eq."nose-hoover") then
    TotalDim = cgrid%npts*3.0D0
    Do i = 1, NumField
      Call random_number(collision)
      gm(i) = gm(i) + (Ekin(i)-0.5D0*TotalDim*k_bolt_ev*T0/Hartree)*DeltaT/NoseMass(i)/tau
      If(FrozenQ(i).and.(TimeNow.gt.ThermoSteps*DeltaT)) gm(i) = 0.0D0
      If(collision.lt.(1-ReservoirRatio)) gm(i) = 0.0D0
    End do
    If(All(CLAMPQ)) then
      gm(NumField+1) = 0.0D0
    else if(ReservoirQ(2)) then
      gm(NumField+1) = gm(NumField+1) + (Ekin(NumField+1) - 0.5D0*6*k_bolt_ev*T0/Hartree)*DeltaT/NoseMass(NumField+1)/tau
    else
      gm(NumField+1) = 0.0D0
    end if
  else 
    If((TimeNow.gt.ThermoSteps*DeltaT)) then
      TotalDim = 0.0D0
      do i = 1, NumField 
        If(.Not.FrozenQ(i)) TotalDim = TotalDim+cgrid%npts*3
      End do
    else
      TotalDim = cgrid%npts*9.0D0
    end if
    Do i = 1, NumField
      Call random_number(collision)
      gm(i) = gm(i) + (Ekin(1)+Ekin(2)+Ekin(3)-0.5D0*TotalDim*k_bolt_ev*T0/Hartree)*DeltaT/NoseMass(i)/tau
      If(FrozenQ(i).and.(TimeNow.gt.ThermoSteps*DeltaT)) gm(i) = 0.0D0
      If(collision.lt.(1-ReservoirRatio)) gm(i) = 0.0D0
    End do
    i = NumField+1
    If(All(CLAMPQ)) then
      gm(NumField+1) = 0.0D0
    else if(ReservoirQ(2)) then
      gm(NumField+1) = gm(NumField+1) + (Ekin(NumField+1) - 0.5D0*6*k_bolt_ev*T0/Hartree)*DeltaT/NoseMass(NumField+1)/tau
    else
      gm(NumField+1) = 0.0D0
    end if
  End if

End Subroutine ReservoirUpdate

Subroutine BaroState(e0ij, VolumeForce)
  Use Parameters
  Implicit none
  Real*8,  Intent(in)     :: e0ij(3,3) 
  Real*8,  Intent(out)    :: VolumeForce(3,3)
  Real*8                  :: delta, L
  Integer                 :: i

  If(All(CLAMPQ)) then
    VolumeForce = 0.0D0
  Else
    VolumeForce(1,1) = -1.0D0 - e0ij(2,2) - e0ij(3,3) - e0ij(2,2)*e0ij(3,3) + e0ij(2,3)**2
    VolumeForce(2,2) = -1.0D0 - e0ij(1,1) - e0ij(3,3) - e0ij(1,1)*e0ij(3,3) + e0ij(1,3)**2
    VolumeForce(3,3) = -1.0D0 - e0ij(1,1) - e0ij(2,2) - e0ij(1,1)*e0ij(2,2) + e0ij(1,2)**2
    VolumeForce(2,3) = 2.0D0*(e0ij(2,3) + e0ij(1,1)*e0ij(2,3) - e0ij(1,2)*e0ij(1,3))
    VolumeForce(1,3) = 2.0D0*(e0ij(1,3) + e0ij(2,2)*e0ij(1,3) - e0ij(1,2)*e0ij(2,3))
    VolumeForce(1,2) = 2.0D0*(e0ij(1,2) + e0ij(3,3)*e0ij(1,2) - e0ij(1,3)*e0ij(2,3))
    VolumeForce(3,2) = VolumeForce(2,3)
    VolumeForce(2,1) = VolumeForce(1,2)
    VolumeForce(3,1) = VolumeForce(1,3)
    
    VolumeForce = VolumeForce*Pressure/(mass(NumField+1)*cgrid%npts)
  End If

End Subroutine BaroState

Subroutine ApplyEfield(TimeNow, Forces)
  Use Parameters
  Use Constants
  Use Fileparser
  Use Aux

  Implicit none
  character(20)          :: keyword
  Real*8,  Intent(in)    :: TimeNow
  Real*8,  Intent(inout) :: Forces(Max(FieldDim, 6),NumField+1,cgrid%n1, cgrid%n2, cgrid%n3)
  Real(dp)               :: Efield(3,cgrid%n1, cgrid%n2, cgrid%n3), Etmp(3)
  Real*8                 :: omega, mu, sigmar, sigmat, beta, t, tt
  Real*8                 :: rr, r0(3), theta, GaussianT, GaussianR
  Integer                :: i, ix, iy, iz

  Efield = 0.0_dp

  if(TimeNow.ge.ThermoSteps*DeltaT) then
    t = (TimeNow-ThermoSteps*DeltaT)*time_fs
    tt = t
    if(TrainQ) tt = Mod(TimeNow/DeltaT-ThermoSteps, Float(TrainRate))*DeltaT*time_fs
    keyword=trim(EfieldType)
    call caps2small(keyword)
    Select case(keyword)
      case('oam')
        omega  = EAmp(1)/1.0D3
        sigmat = 600/sqrt(2*DLog(2.0_dp))/2
        sigmar = EPhi(1)/sqrt(2*DLog(2.0_dp))/2
        mu    = 6*sigmat
        r0 = [cgrid%n1/2-0.5, cgrid%n2/2-0.5, cgrid%n3/2-0.5]
        Do ix = 1, cgrid%n1
        Do iy = 1, cgrid%n2
        Do iz = 1, cgrid%n3
          rr = sqrt((ix-r0(1))**2+(iy-r0(2))**2)
          theta = ATAN2(iy-r0(2), ix-r0(1)) + pi
          GaussianT = Exp(-(tt-mu)**2/(2*sigmat**2))
          GaussianT = 1.0_dp
          GaussianR = Exp(-4*(rr-sigmar)**2/sigmar/Ephi(4))
          Efield(1,ix,iy,iz) = EAmp(2)*GaussianT*GaussianR*Cos(2*pi*omega*tt+EPhi(2)*theta) + GateField(2)
          Efield(2,ix,iy,iz) = EAmp(3)*GaussianT*GaussianR*Sin(2*pi*omega*tt+EPhi(3)*theta) + GateField(3)
          Efield(3,ix,iy,iz) = GateField(4)
        End do
        End do
        End do
      case('oampulse')
        omega  = EAmp(1)/1.0D3
        sigmat = 600/sqrt(2*DLog(2.0_dp))/2
        sigmar = EPhi(1)/sqrt(2*DLog(2.0_dp))/2
        mu    = 6*sigmat
        r0 = [cgrid%n1/2-0.5, cgrid%n2/2-0.5, cgrid%n3/2-0.5]
        Do ix = 1, cgrid%n1
        Do iy = 1, cgrid%n2
        Do iz = 1, cgrid%n3
          rr = sqrt((ix-r0(1))**2+(iy-r0(2))**2)
          theta = ATAN2(iy-r0(2), ix-r0(1)) + pi
          GaussianT = Exp(-(tt-mu)**2/(2*sigmat**2))
          GaussianR = Exp(-4*(rr-sigmar)**2/sigmar/Ephi(4))
          Efield(1,ix,iy,iz) = EAmp(2)*GaussianT*GaussianR*Cos(2*pi*omega*tt+EPhi(2)*theta) + GateField(2)
          Efield(2,ix,iy,iz) = EAmp(3)*GaussianT*GaussianR*Sin(2*pi*omega*tt+EPhi(3)*theta) + GateField(3)
          Efield(3,ix,iy,iz) = GateField(4) 
        End do
        End do
        End do
      case('gaussian')
        omega  = EAmp(1)/1.0D3
        sigmat = EPhi(1)/sqrt(2*DLog(2.0_dp))/2
        mu    = 6*sigmat
        GaussianT = Exp(-(tt-mu)**2/(2*sigmat**2))
        Do i = 1, 3
          Etmp(i) = EAmp(i+1)*GaussianT*Sin(2*pi*omega*tt+EPhi(i+1)*pi)
        End do
      case('chirped')
        omega  = EAmp(1)/1.0D3
        sigmat = EPhi(1)/sqrt(2*DLog(2.0_dp))/2
        mu    = 6*sigmat
        beta  = -0.0002
        GaussianT = Exp(-(tt-mu)**2/(2*sigmat**2))
        Do i = 1, 3
          Etmp(i) = EAmp(i+1)*GaussianT*Sin(2*pi*omega*tt+EPhi(i+1)*pi+beta*(tt-mu)**2)
        End do
      case('neuromorphic')
        mu = 50
        Do i = 1, 3
          Etmp(i) = EAmp(i+1)*Exp(-0.05D0*(tt-mu))*(-0.1D0+1.0D0/Cosh(0.5D0*(tt-mu))**2-0.1D0*Tanh(0.5D0*(tt-mu)))
        End do
      case('sin')
        omega = EAmp(1)/1.0D3
        Do i = 1, 3
          Etmp(i) = EAmp(i+1)*Sin(2*pi*omega*tt+EPhi(i+1)*pi)
        End do
      case('step')
        omega = EAmp(1)/1.0D3
        Do i = 1, 3
          Etmp(i) = 2.0D0*EAmp(i+1)*(Ceiling(Sin(2*pi*omega*tt+EPhi(i+1)*pi))-0.5D0)
        End do
      case('dc')
        Do i = 1, 3
          Etmp(i) = EAmp(i+1)
        End do
      case default
        write(*,*) "Unknown type of electric field!"
        call abort

    End Select

    if(keyword.ne."oam") then
      Do iz = 1, cgrid%n3
      Do iy = 1, cgrid%n2
      Do ix = 1, cgrid%n1
      Do i = 1, 3
        Efield(i,ix,iy,iz) = chop(Etmp(i) + GateField(i+1), 1.0D-20)
      End do
      End do
      End do
      End do
    end if

    Do ix = 1, cgrid%n1
    Do iy = 1, cgrid%n2
    Do iz = 1, cgrid%n3
      if (ix > nint(GateField(1)) .and. ix <= cgrid%n1-nint(GateField(1)) .and. &
          iy > nint(GateField(1)) .and. iy <= cgrid%n2-nint(GateField(1)) .and. &
          iz > nint(GateField(1)) .and. iz <= cgrid%n3-nint(GateField(1))) then
        Do i = 1, NumField
          Forces(1:3,i,ix,iy,iz) = Forces(1:3,i,ix,iy,iz) + Efield(:,ix,iy,iz)*FieldCharge(i)*alat
        End do
      end if
    End do
    End do
    End do
  End If

End Subroutine ApplyEfield
