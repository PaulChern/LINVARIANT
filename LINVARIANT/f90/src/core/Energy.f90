Module Energy
  Use Constants
  Use Aux
  Use Parameters
  Use Lattice
  Use Ewald

  Implicit none
  Contains

  Include "SiteEnergyExpr.f90"
 
  Include "GetSiteEnergy.f90"
  Include "GetSiteEnergyG.f90"
  Include "GetSiteEnergyH.f90"
 
  Include "GetSiteEnergyGp.f90"
  Include "GetSiteEnergyGs.f90"
  Include "GetSiteEnergyGsp.f90"
 
  Include "GetSiteEnergyHp.f90"
  Include "GetSiteEnergyHu.f90"
  Include "GetSiteEnergyHup.f90"
  Include "GetSiteEnergyHsu.f90"
  Include "GetSiteEnergyHsp.f90"
  Include "GetSiteEnergyHsup.f90"

  Function GetEkin(dFieldsdt, de0ijdt) Result(Ekin)
    Implicit none
    Real*8,  Intent(in)  :: dFieldsdt(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8,  Intent(in)  :: de0ijdt(3,3)
    Real*8               :: Ekin(2), detadt(6), etaV
    Integer              :: i, ifield, ix, iy, iz, icell

    Ekin = 0.0D0
    detadt = eij2eta(de0ijdt)

    do icell = 1, cgrid%ncells
      ix = GridFold(1, icell)
      iy = GridFold(2, icell)
      iz = GridFold(3, icell)
      do ifield = 1, NumField 
        do i = 1, FieldDim
          Ekin(1) = Ekin(1) + 0.5d0*mass(ifield)*FieldsBinary(i,ifield)*dFieldsdt(i, ifield, ix, iy, iz)**2*alat**2
        end do
      end do
    end do

    do i = 1, 6
      If(.Not.CLAMPQ(i)) Ekin(2) = Ekin(2) + 0.5d0*cgrid%npts*mass(NumField+1)*(detadt(i)*alat)**2
    end do

  End Function GetEkin

  Function ResolveEkin(dFieldsdt, de0ijdt) Result(Ekin)
    Implicit none
    Real*8,  Intent(in)  :: dFieldsdt(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8,  Intent(in)  :: de0ijdt(3,3)
    Real*8               :: Ekin(NumField+1), detadt(6), etaV
    Integer              :: i, ifield, ix, iy, iz, icell

    Ekin = 0.0D0
    detadt = eij2eta(de0ijdt)

    do icell = 1, cgrid%ncells
      ix = GridFold(1, icell)
      iy = GridFold(2, icell)
      iz = GridFold(3, icell)
      do ifield = 1, NumField
        do i = 1, FieldDim
          Ekin(ifield) = Ekin(ifield) + 0.5d0*mass(ifield)*FieldsBinary(i,ifield)*dFieldsdt(i, ifield, ix, iy, iz)**2*alat**2
        end do
      end do
    end do

    do i = 1, 6
      If(.Not.CLAMPQ(i)) Ekin(NumField+1) = Ekin(NumField+1) + 0.5d0*cgrid%npts*mass(NumField+1)*(detadt(i)*alat)**2
    end do

  End Function ResolveEkin

  Function GetEtot(Fields, e0ij) Result(Etot)
    Implicit none
    Real*8,  Intent(in)    :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8,  Intent(in)    :: e0ij(3,3)
    Real*8                 :: Etot, omptemp
    Integer                :: icell, ix, iy, iz

    Etot = 0.0D0

    if (DipoleQ) then
!      !$OMP    PARALLEL DEFAULT(SHARED) PRIVATE(ix,iy,iz,omptemp)
        omptemp = 0.0_dp
!      !$OMP    DO
        do icell = 1, cgrid%ncells
          ix = GridFold(1, icell)
          iy = GridFold(2, icell)
          iz = GridFold(3, icell)
          omptemp = omptemp + GetSiteEnergy(ix, iy, iz, Fields, e0ij)
        end do
!      !$OMP    END DO
!      !$OMP    CRITICAL
        Etot = Etot + omptemp
!      !$OMP    END CRITICAL
!      !$OMP    END PARALLEL
      do icell = 1, cgrid%ncells
        ix = GridFold(1, icell)
        iy = GridFold(2, icell)
        iz = GridFold(3, icell)
        Etot = Etot + GetSiteEnergyEwald(ix, iy, iz, Fields)
      end do
    else
      !$OMP    PARALLEL DEFAULT(SHARED) PRIVATE(ix,iy,iz,omptemp)
        omptemp = 0.0_dp
      !$OMP    DO
        do icell = 1, cgrid%ncells
          ix = GridFold(1, icell)
          iy = GridFold(2, icell)
          iz = GridFold(3, icell)
          omptemp = omptemp + GetSiteEnergy(ix, iy, iz, Fields, e0ij) 
        end do
      !$OMP    END DO
      
      !$OMP    CRITICAL
        Etot = Etot + omptemp
      !$OMP    END CRITICAL

      !$OMP    END PARALLEL
    end if

  End Function GetEtot

End Module Energy
