Module Energy
  Use Constants
  Use Aux
  Use Parameters
  Use Lattice
  Use Ewald

  Implicit none
  Contains

  Include "GetEOnSiteDisp.f90" 
  Include "GetEOnSiteEps.f90"
  Include "GetEOnSite.f90"
  Include "GetEOnSiteJijsm.f90"
  Include "GetEOnSiteDispu.f90"
  Include "GetEOnSiteEpsu.f90"
  Include "GetEOnSiteJij.f90"
  Include "GetEOnSiteu.f90"
  Include "GetEOnSiteEpsDisp.f90"
  Include "GetEOnSiteJijhf.f90"
  Include "GetEOnSiteJijsmhf.f90"

  Function GetEkin(dFieldsdt, de0ijdt) Result(Ekin)
    Implicit none
    Real*8,  Intent(in)  :: dFieldsdt(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8,  Intent(in)  :: de0ijdt(3,3)
    Real*8               :: Ekin(2), detadt(6), etaV
    Integer              :: i, ifield, ix, iy, iz

    Ekin = 0.0D0
    detadt = eij2eta(de0ijdt)

    do iz = 1, cgrid%n3
      do iy = 1, cgrid%n2
        do ix = 1, cgrid%n1
          do ifield = 1, NumField 
            do i = 1, FieldDimList(ifield)
              Ekin(1) = Ekin(1) + 0.5d0*mass(ifield)*dFieldsdt(i, ifield, ix, iy, iz)**2*(7.464D0)**2
            end do
          end do
        end do
      end do
    end do

!    do iz = 1, cgrid%n3
!      do iy = 1, cgrid%n2
!        do ix = 1, cgrid%n1
!          do i = 1, 3
!            etaV = de0ijdt(i,1)*(ix-1)+de0ijdt(i,2)*(iy-1)+de0ijdt(i,3)*(iz-1)
!            Ekin = Ekin + 0.5d0*mass(3)*etaV**2*(7.464D0)**2
!          end do
!        end do
!      end do
!    end do

    do i = 1, 6
      Ekin(2) = Ekin(2) + 0.5d0*cgrid%npts*mass(NumField+1)*(detadt(i)*7.464D0)**2
    end do

  End Function GetEkin

  Function ResolveEkin(dFieldsdt, de0ijdt) Result(Ekin)
    Implicit none
    Real*8,  Intent(in)  :: dFieldsdt(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8,  Intent(in)  :: de0ijdt(3,3)
    Real*8               :: Ekin(NumField+1), detadt(6), etaV
    Integer              :: i, ifield, ix, iy, iz

    Ekin = 0.0D0
    detadt = eij2eta(de0ijdt)

    do iz = 1, cgrid%n3
      do iy = 1, cgrid%n2
        do ix = 1, cgrid%n1
          do ifield = 1, NumField
            do i = 1, FieldDimList(ifield)
              Ekin(ifield) = Ekin(ifield) + 0.5d0*mass(ifield)*dFieldsdt(i, ifield, ix, iy, iz)**2*(7.464D0)**2
            end do
          end do
        end do
      end do
    end do

    do i = 1, 6
      Ekin(NumField+1) = Ekin(NumField+1) + 0.5d0*cgrid%npts*mass(NumField+1)*(detadt(i)*7.464D0)**2
    end do

  End Function ResolveEkin

  Function GetEtot(Fields, e0ij) Result(Etot)
    Implicit none
    Real*8,  Intent(in)    :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8,  Intent(in)    :: e0ij(3,3)
    Real*8                 :: Etot
    Integer                :: ix, iy, iz

    Etot = 0.0D0

    if (DipoleQ) then
      do iz = 1, cgrid%n3
        do iy = 1, cgrid%n2
          do ix = 1, cgrid%n1
            Etot = Etot + GetEOnSite(ix, iy, iz, Fields, e0ij) + GetEOnSiteEwald(ix, iy, iz, Fields) + GetEOnSiteJijsmhf(ix,iy,iz,Fields,e0ij)
          end do
        end do
      end do
    else
      do iz = 1, cgrid%n3
        do iy = 1, cgrid%n2
          do ix = 1, cgrid%n1
            Etot = Etot + GetEOnSite(ix, iy, iz, Fields, e0ij) + GetEOnSiteJijsmhf(ix,iy,iz,Fields,e0ij)
          end do
        end do
      end do
    end if

  End Function GetEtot

End Module Energy
