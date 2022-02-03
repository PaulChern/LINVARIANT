Module Wannier

  Use Parameters
  Use Aux
  Use Inputs
  Use Lattice
  Use fft

  Implicit none
  Contains

  Function ResolveWannSite(iwann) result(isite)
    Implicit none
    Integer,  Intent(in)  :: iwann
    Integer               :: isite
    Integer               :: i, tmp, lmin, rmin
    integer               :: ind0(NumWannSites), ind(NumWannSites,NumWannSites)

    tmp = 0
    Do i = 1, NumWannSites
      ind0(i) = SiteOrbInfo(2,i)*OrbMul + tmp
      tmp = ind0(i)
    End do
    Do i = 1, NumWannSites
      ind(1,i) = ind0(i) - SiteOrbInfo(2,i)*OrbMul + 1
      ind(2,i) = ind0(i)
    End do

    isite = 0
    Do i = 1, NumWannSites
      rmin = ind(2,i) - iwann
      lmin = iwann - ind(1,i) 
      if(lmin.ge.0.and.rmin.ge.0) then
        isite = SiteOrbInfo(1,i)
        return
      end if
    End do

  End function ResolveWannSite

  Subroutine ResolveWannSiteInd(mul)

    Implicit None
    Integer,intent(in) :: mul
    Integer            :: i, iwann,  tmp, lmin, rmin
    Integer            :: ind0(NumWannSites),ind(2,NumWannSites)

    tmp = 0
    Do i = 1, NumWannSites
      ind0(i) = mul*SiteOrbInfo(2,i) + tmp
      tmp = ind0(i)
    End do
    Do i = 1, NumWannSites
      ind(1,i) = ind0(i) - mul*SiteOrbInfo(2,i) + 1
      ind(2,i) = ind0(i)
    End do

    Allocate(WannSiteInd(NumWann))
    WannSiteInd = 0
    Do iwann = 1, NumWann
      Do i = 1, NumWannSites
        rmin = ind(2,i) - iwann
        lmin = iwann - ind(1,i)
        if(lmin.ge.0.and.rmin.ge.0) WannSiteInd(iwann) = SiteOrbInfo(1,i)
      End do
    End do

  End subroutine ResolveWannSiteInd

  Subroutine GetWannDist
    Implicit none
    Integer            :: m, n

    Allocate(WannDist(3,NumWann,NumWann))
    do m = 1, NumWann
      do n = 1, NumWann
        WannDist(:,m,n) = unitcell%site(:,WannSiteInd(n)) - unitcell%site(:,WannSiteInd(m))
      end do
    end do

  End subroutine

  Subroutine ResolveWannBlockInd(mul)

    Implicit None
    Integer,intent(in)          :: mul
    Integer                     :: tmp
    Integer                     :: i,ind0(NumWannSites)

    Allocate(WannBlockInd(2,NumWannSites))
    tmp = 0
    Do i = 1, NumWannSites
      ind0(i) = mul*SiteOrbInfo(2,i) + tmp
      tmp = ind0(i)
    End do
    Do i = 1, NumWannSites
      WannBlockInd(1,i) = ind0(i) - mul*SiteOrbInfo(2,i) + 1
      WannBlockInd(2,i) = ind0(i)
    End do

  End subroutine ResolveWannBlockInd

  Subroutine ReadWannier90(filename,SpinDim)

    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in)          :: SpinDim
    integer                      :: iow90, iow90_up=2001, iow90_dn=2002
    integer                      :: iow90_socx=1901, iow90_socy=1902, iow90_socz=1903
    integer                      :: iow90_up0=1801, iow90_upx=1802, iow90_upy=1803, iow90_upz=1804
    integer                      :: iow90_upix=1805, iow90_upiy=1806, iow90_upiz=1807
    integer                      :: iow90_dn0=1808, iow90_dnx=1809, iow90_dny=1810, iow90_dnz=1811
    integer                      :: iow90_dnix=1812, iow90_dniy=1813, iow90_dniz=1814
    integer                      :: ioorb=1993, ikbz=1994
    integer                      :: i, j, is, m, n, NumFiles 
    integer                      :: Rx, Ry, Rz, iw1, iw2 
    Real(dp)                     :: Re, Im
    integer, allocatable         :: ndgen(:,:)
    Integer                      :: FileHandle = 1919
    logical                      :: file_exists

    Call GetBZKList
    open(FileHandle,file=trim(Solver)//'.out/kbz.dat',form='formatted',status='unknown')
    write(FileHandle, "(I10)") kgrid%npts
    Do i = 1, kgrid%npts
      write(FileHandle, "(4F25.15)") kbz(:,i)
    End do

    if(trim(Solver).eq."JGreen") then
      If(SpinDim.eq.3) then
        NumFiles = 3
        iow90 = 1900
        call open_input_file(iow90_socx,'wannier/'//trim(filename)//'_socx.dat')
        call open_input_file(iow90_socy,'wannier/'//trim(filename)//'_socy.dat')
        call open_input_file(iow90_socz,'wannier/'//trim(filename)//'_socz.dat')
      else if(SpinDim.eq.2) then
        NumFiles = 2
        iow90 = 2000
        call open_input_file(iow90_up,'wannier/'//trim(filename)//'_up.dat')
        call open_input_file(iow90_dn,'wannier/'//trim(filename)//'_dn.dat')
      else
        write(*,*) "SpinDim need to be set 2 or 3 in LINVARIANT.inp"
        Call Exit
      end if

      Allocate(NumRpts(NumFiles))
      do is = 1, NumFiles
        Read(iow90+is, *)
        Read(iow90+is, *) NumWann
        Read(iow90+is, *) NumRpts(is)
      end do
      Allocate(ndgen(Maxval(NumRpts),3))
      do is = 1, NumFiles
        Read(iow90+is, *) (ndgen(i,is), i=1,NumRpts(is))
      end do
      Allocate(HRmn(NumWann,NumWann,Maxval(NumRpts),NumFiles))
      Allocate(Ti0(4,Maxval(NumRpts),NumFiles))

      do is = 1, NumFiles
        Do i = 1, NumRpts(is)
          Do m = 1, NumWann
            Do n = 1, NumWann
              Read(iow90+is,*) Rx, Ry, Rz, iw1, iw2, Re, Im
              HRmn(iw1,iw2,i,is) = dcmplx(Re, Im)
              Ti0(1,i,is) = Rx
              Ti0(2,i,is) = Ry
              Ti0(3,i,is) = Rz
              Ti0(4,i,is) = ndgen(i,is)
            End do
          End do
        End do
      End do
    else if(trim(Solver).eq."FGreen") then
      iow90 = 1800
      if(SpinDim.eq.2) then
        NumFiles = 2
        call open_input_file(iow90_up0,'wannier/'//trim(filename)//'_up.dat')
!        call open_input_file(iow90_upx,'wannier/'//trim(filename)//'_up_disp+x.dat')
!        call open_input_file(iow90_upy,'wannier/'//trim(filename)//'_up_disp+y.dat')
!        call open_input_file(iow90_upz,'wannier/'//trim(filename)//'_up_disp+z.dat')
!        call open_input_file(iow90_upix,'wannier/'//trim(filename)//'_up_disp-x.dat')
!        call open_input_file(iow90_upiy,'wannier/'//trim(filename)//'_up_disp-y.dat')
!        call open_input_file(iow90_upiz,'wannier/'//trim(filename)//'_up_disp-z.dat')
        call open_input_file(iow90_dn0,'wannier/'//trim(filename)//'_dn.dat')
!        call open_input_file(iow90_dnx,'wannier/'//trim(filename)//'_dn_disp+x.dat')
!        call open_input_file(iow90_dny,'wannier/'//trim(filename)//'_dn_disp+y.dat')
!        call open_input_file(iow90_dnz,'wannier/'//trim(filename)//'_dn_disp+z.dat')
!        call open_input_file(iow90_dnix,'wannier/'//trim(filename)//'_dn_disp-x.dat')
!        call open_input_file(iow90_dniy,'wannier/'//trim(filename)//'_dn_disp-y.dat')
!        call open_input_file(iow90_dniz,'wannier/'//trim(filename)//'_dn_disp-z.dat')
      else
        NumFiles = 1
        call open_input_file(iow90_up0,'wannier/'//trim(filename)//'.dat')
!        call open_input_file(iow90_upx,'wannier/'//trim(filename)//'_disp+x.dat')
!        call open_input_file(iow90_upy,'wannier/'//trim(filename)//'_disp+y.dat')
!        call open_input_file(iow90_upz,'wannier/'//trim(filename)//'_disp+z.dat')
!        call open_input_file(iow90_upix,'wannier/'//trim(filename)//'_disp-x.dat')
!        call open_input_file(iow90_upiy,'wannier/'//trim(filename)//'_disp-y.dat')
!        call open_input_file(iow90_upiz,'wannier/'//trim(filename)//'_disp-z.dat')
      end if
      Allocate(NumRpts(NumFiles))
      do is = 1, NumFiles
        Read(iow90+is, *)
        Read(iow90+is, *) NumWann
        Read(iow90+is, *) NumRpts(is)
      end do
      Allocate(ndgen(Maxval(NumRpts),NumFiles))
      do is = 1, NumFiles
        Read(iow90+is, *) (ndgen(i,is), i=1,NumRpts(is))
      end do
      Allocate(HRmn(NumWann,NumWann,Maxval(NumRpts),NumFiles))
      Allocate(Ti0(4,Maxval(NumRpts),NumFiles))

      Do is = 1, NumFiles
        Do i = 1, NumRpts(is)
          Do m = 1, NumWann
            Do n = 1, NumWann
              Read(iow90+is,*) Rx, Ry, Rz, iw1, iw2, Re, Im
              HRmn(iw1,iw2,i,is) = dcmplx(Re, Im)
              Ti0(1,i,is) = Rx
              Ti0(2,i,is) = Ry
              Ti0(3,i,is) = Rz
              Ti0(4,i,is) = ndgen(i,is)
            End do
          End do
        End do
      End do
    end if

    write(*,*) "TB hamiltonian loaded."

  End subroutine ReadWannier90

  subroutine LoadWannFunc
    Use omp_lib
    implicit none
    character(72)           :: filename
    Integer                 :: iw, jw, is, idisp
    Integer                 :: ix, iy, iz
    Complex(dp)             :: WannCopy(3*rgrid%n1,3*rgrid%n2,3*rgrid%n3)
    Complex(dp)             :: WannFFTQx(3*rgrid%n1*potim,3*rgrid%n2,3*rgrid%n3)
    Complex(dp)             :: WannFFTQy(3*rgrid%n1,3*rgrid%n2*potim,3*rgrid%n3)
    Complex(dp)             :: WannFFTQz(3*rgrid%n1,3*rgrid%n2,3*rgrid%n3*potim)
    Complex(dp)             :: WannFFTRx(3*rgrid%n1*potim,3*rgrid%n2,3*rgrid%n3)
    Complex(dp)             :: WannFFTRy(3*rgrid%n1,3*rgrid%n2*potim,3*rgrid%n3)
    Complex(dp)             :: WannFFTRz(3*rgrid%n1,3*rgrid%n2,3*rgrid%n3*potim)
    Complex(dp)             :: WannShiftx(3*rgrid%n1*potim,3*rgrid%n2,3*rgrid%n3)
    Complex(dp)             :: WannShifty(3*rgrid%n1,3*rgrid%n2*potim,3*rgrid%n3)
    Complex(dp)             :: WannShiftz(3*rgrid%n1,3*rgrid%n2,3*rgrid%n3*potim)
    Complex(dp)             :: charge

!    prefix=['disp0', 'disp+x', 'disp+y', 'disp+z', 'disp-x', 'disp-y', 'disp-z']

    Allocate(WannFunc(3*rgrid%n1,3*rgrid%n2,3*rgrid%n3,NumWann,7,SpinDim))
    WannFunc = cmplx(0.0_dp,0.0_dp)

    Do is = 1, SpinDim
      Do iw = 1, NumWann
        filename='wannier/wannier90_'//trim(int2str5(iw))//'.'//trim(int2str(is))//'.xsf'
        WannCopy = ReadXSF(filename)/sqrt(Real(rgrid%npts,dp))
        charge = Sum(WannCopy*Conjg(WannCopy))
        WannFunc(:,:,:,iw,1,is) = WannCopy
        if(iw.eq.11) WannFunc(:,:,:,iw,1,is) = -WannFunc(:,:,:,iw,1,is)
        if(iw.eq.12) WannFunc(:,:,:,iw,1,is) = -WannFunc(:,:,:,iw,1,is)
        if(iw.eq.14) WannFunc(:,:,:,iw,1,is) = -WannFunc(:,:,:,iw,1,is)
        if(iw.eq.15) WannFunc(:,:,:,iw,1,is) = -WannFunc(:,:,:,iw,1,is)
        if(iw.eq.17) WannFunc(:,:,:,iw,1,is) = -WannFunc(:,:,:,iw,1,is)
        write(*,'(A20,A40,A10,2F10.4)') "wannier basis from", "./"//filename, "norm: ", charge
      end do
    end do

    Do is = 1, SpinDim 
      Do iw = 1, NumWann
        ! shift [100]
        write(*,'(A20,I5,A3,F10.4,A17)') "shift wannier basis ", iw, "by ", 1.0/potim, "along x direction"
        Call fft3d_c2c(1,-1,WannFunc(:,:,:,iw,1,is),3*rgrid%n1,3*rgrid%n2,3*rgrid%n3,WannFFTQx,3*rgrid%n1*potim,3*rgrid%n2,3*rgrid%n3)
        Call fft3d_c2c(-1,1,WannFFTQx,3*rgrid%n1*potim,3*rgrid%n2,3*rgrid%n3,WannFFTRx,3*rgrid%n1*potim,3*rgrid%n2,3*rgrid%n3)
        WannShiftx = CSHIFT(WannFFTRx,-1,DIM=1)
        do iz = 1, 3*rgrid%n3
          do iy = 1, 3*rgrid%n2
            do ix = 1, 3*rgrid%n1
              WannFunc(ix,iy,iz,iw,2,is) = WannShiftx(1+(ix-1)*potim,iy,iz)
            end do
          end do
        end do
        WannShiftx = CSHIFT(WannFFTRx,1,DIM=1)
        do iz = 1, 3*rgrid%n3
          do iy = 1, 3*rgrid%n2
            do ix = 1, 3*rgrid%n1
              WannFunc(ix,iy,iz,iw,5,is) = WannShiftx(1+(ix-1)*potim,iy,iz)
            end do
          end do
        end do

        ! shift [010]
        write(*,'(A20,I5,A3,F10.4,A17)') "shift wannier basis ", iw, "by ", 1.0/potim, "along y direction"
        Call fft3d_c2c(1,-1,WannFunc(:,:,:,iw,1,is),3*rgrid%n1,3*rgrid%n2,3*rgrid%n3,WannFFTQy,3*rgrid%n1,3*rgrid%n2*potim,3*rgrid%n3)
        Call fft3d_c2c(-1,1,WannFFTQy,3*rgrid%n1,3*rgrid%n2*potim,3*rgrid%n3,WannFFTRy,3*rgrid%n1,3*rgrid%n2*potim,3*rgrid%n3)
        WannShifty = CSHIFT(WannFFTRy,-1,DIM=1)
        do iz = 1, 3*rgrid%n3
          do iy = 1, 3*rgrid%n2
            do ix = 1, 3*rgrid%n1
              WannFunc(ix,iy,iz,iw,3,is) = WannShifty(ix,1+(iy-1)*potim,iz)
            end do
          end do
        end do
        WannShifty = CSHIFT(WannFFTRy,1,DIM=1)
        do iz = 1, 3*rgrid%n3
          do iy = 1, 3*rgrid%n2
            do ix = 1, 3*rgrid%n1
              WannFunc(ix,iy,iz,iw,6,is) = WannShifty(ix,1+(iy-1)*potim,iz)
            end do
          end do
        end do

        ! shift [001]
        write(*,'(A20,I5,A3,F10.4,A17)') "shift wannier basis ", iw, "by ", 1.0/potim, "along z direction"
        Call fft3d_c2c(1,-1,WannFunc(:,:,:,iw,1,is),3*rgrid%n1,3*rgrid%n2,3*rgrid%n3,WannFFTQz,3*rgrid%n1,3*rgrid%n2,3*rgrid%n3*potim)
        Call fft3d_c2c(-1,1,WannFFTQz,3*rgrid%n1,3*rgrid%n2,3*rgrid%n3,WannFFTRz,3*rgrid%n1,3*rgrid%n2,3*rgrid%n3*potim)
        WannShiftz = CSHIFT(WannFFTRz,-1,DIM=1)
        do iz = 1, 3*rgrid%n3
          do iy = 1, 3*rgrid%n2
            do ix = 1, 3*rgrid%n1
              WannFunc(ix,iy,iz,iw,4,is) = WannShiftz(ix,iy,1+(iz-1)*potim)
            end do
          end do
        end do
        WannShiftz = CSHIFT(WannFFTRz,1,DIM=1)
        do iz = 1, 3*rgrid%n3
          do iy = 1, 3*rgrid%n2
            do ix = 1, 3*rgrid%n1
              WannFunc(ix,iy,iz,iw,7,is) = WannShiftz(ix,iy,1+(iz-1)*potim)
            end do
          end do
        end do
      End do
    End do

    do is = 1, SpinDim
      do idisp = 1, 7
        do iw = 1, NumWann
          WannCopy = WannFunc(:,:,:,iw,idisp,is)
          do iz = 1, 3*rgrid%n3
            do iy = 1, 3*rgrid%n2
              do ix = 1, 3*rgrid%n1
                WannFunc(ix,iy,iz,iw,idisp,is) = WannCopy(Mod(ix,3*rgrid%n1)+1,Mod(iy,3*rgrid%n2)+1,Mod(iz,3*rgrid%n3)+1)
              end do
            end do
          end do
        end do
      end do
    end do
   
!    call unk2wann
    write(*,*) "All wannier functions are loaded and their displacements are performed."
    
  End subroutine LoadWannFunc

  subroutine LoadWannBasis(prefix,shift)
    Use omp_lib
    implicit none
    Real(dp), intent(in)         :: shift(3)
    character(len=*), intent(in) :: prefix
    character(72)                :: filename
    Integer                      :: iw, jw, is
    Integer                      :: ix, iy, iz
    Complex(dp)                  :: WannCopy(3*rgrid%n1,3*rgrid%n2,3*rgrid%n3)
    Complex(dp)                  :: WannG(3*rgrid%n1,3*rgrid%n2,3*rgrid%n3)
    Complex(dp)                  :: charge


    Allocate(WannBasis(3*rgrid%n1,3*rgrid%n2,3*rgrid%n3,NumWann,SpinDim))
    WannBasis = cmplx(0.0_dp,0.0_dp)

    Do is = 1, SpinDim
      Do iw = 1, NumWann
        filename='wannier/wannier90_'//trim(prefix)//'_'//trim(int2str5(iw))//'.'//trim(int2str(is))//'.xsf'
        WannCopy = ReadXSF(filename)/sqrt(Real(rgrid%npts,dp))
        charge = Sum(WannCopy*Conjg(WannCopy))
        WannBasis(:,:,:,iw,is) = WannCopy
        write(*,'(A20,A40,A10,2F10.4)') "wannier basis from", "./"//filename, "norm: ", charge
      end do
    end do

    Call WannOriginShift(shift)

    do is = 1, SpinDim
      do iw = 1, NumWann
        WannCopy = WannBasis(:,:,:,iw,is)
        do iz = 1, 3*rgrid%n3
          do iy = 1, 3*rgrid%n2
            do ix = 1, 3*rgrid%n1
              WannBasis(ix,iy,iz,iw,is) = WannCopy(Mod(ix,3*rgrid%n1)+1,Mod(iy,3*rgrid%n2)+1,Mod(iz,3*rgrid%n3)+1)
            end do
          end do
        end do
      end do
    end do

!    call unk2wann
    write(*,*) "All wannier basis loaded."

  End subroutine LoadWannBasis

  Subroutine WannOriginShift(shift)
    implicit none
    Real(dp), intent(in)         :: shift(3)
    Integer                      :: iw, jw, is
    Integer                      :: ix, iy, iz
    Integer                      :: igx, igy, igz
    Real(dp)                     :: s1, s2, s3
    Complex(dp)                  :: WannTmp
    Complex(dp)                  :: WannG(3*rgrid%n1,3*rgrid%n2,3*rgrid%n3)

    s1 = shift(1)
    s2 = shift(2)
    s3 = shift(3)

    Do is = 1, SpinDim
      Do iw = 1, NumWann
        Call fft3d_c2c(1,-1,WannBasis(:,:,:,iw,is),3*rgrid%n1,3*rgrid%n2,3*rgrid%n3,WannG,3*rgrid%n1,3*rgrid%n2,3*rgrid%n3)
        WannTmp = cmplx(0.0_dp,0.0_dp)
        do iz = 1, 3*rgrid%n3 
        do iy = 1, 3*rgrid%n2 
        do ix = 1, 3*rgrid%n1
          WannTmp = cmplx(0.0_dp,0.0_dp)
          do igz = 1, 3*rgrid%n3
          do igy = 1, 3*rgrid%n2
          do igx = 1, 3*rgrid%n1
            WannTmp = WannTmp + &
            1/(rgrid%npts)*WannG(igx,igy,igz)*Exp(2*pi*ii*((igx-1)*(ix-1-s1)+(igy-1)*(iy-1-s2)+(igz-1)*(iz-1-s3))/rgrid%npts)
          end do
          end do
          end do
          WannBasis(ix,iy,iz,iw,is) = WannTmp
        end do
        end do
        end do
      End do
    End do

  End Subroutine WannOriginShift


  Function NablaW(WF) result(DeltaW)
    implicit none
    Complex(dp),intent(in)  :: WF(3*rgrid%n1,3*rgrid%n2,3*rgrid%n3)
    Complex(dp)             :: DeltaW(3,3*rgrid%n1,3*rgrid%n2,3*rgrid%n3)
    Real(dp)                :: da, db, dc
    integer                 :: ix, iy, iz
    integer                 :: x1, y1, z1, xi1, yi1, zi1
    integer                 :: x2, y2, z2, xi2, yi2, zi2

    da = sqrt(sum(unitcell%latt_a(1,:)**2))/rgrid%n1
    db = sqrt(sum(unitcell%latt_a(2,:)**2))/rgrid%n2
    dc = sqrt(sum(unitcell%latt_a(3,:)**2))/rgrid%n3

    Do iz = 1, 3*rgrid%n3
      Do iy = 1, 3*rgrid%n2
        Do ix = 1, 3*rgrid%n1
          x1 = GridPbc(ix+1, 3*rgrid%n1)
          y1 = GridPbc(iy+1, 3*rgrid%n2)
          z1 = GridPbc(iz+1, 3*rgrid%n3)
          xi1 = GridPbc(ix-1, 3*rgrid%n1)
          yi1 = GridPbc(iy-1, 3*rgrid%n2)
          zi1 = GridPbc(iz-1, 3*rgrid%n3)

          x2 = GridPbc(ix+2, 3*rgrid%n1)
          y2 = GridPbc(iy+2, 3*rgrid%n2)
          z2 = GridPbc(iz+2, 3*rgrid%n3)
          xi2 = GridPbc(ix-2, 3*rgrid%n1)
          yi2 = GridPbc(iy-2, 3*rgrid%n2)
          zi2 = GridPbc(iz-2, 3*rgrid%n3)

          DeltaW(1,ix,iy,iz) = (-WF(x2,iy,iz) + 8.0*WF(x1,iy,iz) - 8.0*WF(xi1,iy,iz) + WF(xi2,iy,iz))/real(12*da,dp)
          DeltaW(2,ix,iy,iz) = (-WF(ix,y2,iz) + 8.0*WF(ix,y1,iz) - 8.0*WF(ix,yi1,iz) + WF(ix,yi2,iz))/real(12*db,dp)
          DeltaW(3,ix,iy,iz) = (-WF(ix,iy,z2) + 8.0*WF(ix,iy,z1) - 8.0*WF(ix,iy,zi1) + WF(ix,iy,zi2))/real(12*dc,dp)
        End do
      End do
    End do
    
  End function

  subroutine unk2wann
    implicit none
    integer                 :: num_kpts, num_wann, num_bands, num_inc
    real(dp),allocatable    :: kpts(:,:)
    complex(dp),allocatable :: u_matrix(:,:,:), u_matrix_opt(:,:,:), m_matrix(:,:,:,:)
    logical,allocatable     :: lwindow(:, :), inc_band(:)
    integer,allocatable     :: ndimwin(:)
    complex(dp),allocatable :: Unkr(:,:), r_wvfn_tmp(:,:,:), r_wvfn(:,:,:)
    complex(dp)             :: catmp, wmod
    real(dp)                :: scalfac, tmax, tmaxx, ratio, ratmax
    logical                 :: Qdisentangled
    character(72)           :: unkfile
    character(10)           :: prefix(7)
    Integer                 :: idisp, is, ik, ib, iw, counter, npoint, ierr
    Integer                 :: ix, iy, iz, nx, ny, nz, nxx, nyy, nzz

    prefix=['disp0', 'disp+x', 'disp+y', 'disp+z', 'disp-x', 'disp-y', 'disp-z']

    Call ReadCHK(num_kpts,num_bands,num_wann,Qdisentangled,kpts, &
                 u_matrix,u_matrix_opt,m_matrix,lwindow,ndimwin)

    allocate(WannFunc(-rgrid%n1 : 2*rgrid%n1-1, &
                      -rgrid%n2 : 2*rgrid%n2-1, &
                      -rgrid%n3 : 2*rgrid%n3-1, num_wann, 7, SpinDim), stat=ierr)

    allocate(Unkr(rgrid%npts,num_bands))

    if(Qdisentangled) then
      allocate(r_wvfn_tmp(rgrid%npts, maxval(ndimwin), SpinDim), stat=ierr)
    end if
    allocate(r_wvfn(rgrid%npts, num_wann, SpinDim), stat=ierr)
    allocate(inc_band(num_bands))

    do is = 1, SpinDim
    do idisp = 1, 7
    do ik = 1, num_kpts
      unkfile='UNK_'//trim(prefix(idisp))//'_'//trim(int2str5(ik))//'.'//trim(int2str(is))
      Unkr=ReadUnk(unkfile,num_bands,rgrid%npts)
      write(*,*) is, ik, Sum(Unkr(:,1)*Conjg(Unkr(:,1)))/rgrid%npts
      inc_band = .true.
      if(Qdisentangled) then
        inc_band(:) = lwindow(:, ik)
        num_inc = ndimwin(ik)
      end if

      if (Qdisentangled) then
        counter = 1
        do ib = 1, num_bands
          if(counter.gt.num_inc) exit
          do nx = 1, rgrid%npts
            r_wvfn_tmp(nx, counter, is) = Unkr(nz, ib)
          end do
          if(inc_band(ib)) counter = counter + 1
        end do
      else
        do ib = 1, num_bands
          do nx = 1, rgrid%npts
            r_wvfn(nx, ib, is) = Unkr(nx, ib)
          end do
        end do
      end if

      if (Qdisentangled) then
        r_wvfn = cmplx(0.0_dp,0.0_dp)
        do iw = 1, num_wann
          do ib = 1, num_inc
            r_wvfn(:,iw,is)=r_wvfn(:,iw,is)+u_matrix_opt(ib,iw,ik)*r_wvfn_tmp(:,ib,is)
          end do
        end do
      end if

      ! nxx, nyy, nzz span a parallelogram in the real space mesh, of side
      ! 2*nphir, and centered around the maximum of phi_i, nphimx(i, 1 2 3)
      !
      ! nx ny nz are the nxx nyy nzz brought back to the unit cell in
      ! which u_nk(r)=cptwrb(r,n)  is represented
      !
      ! There is a big performance improvement in looping over num_wann
      ! in the inner loop. This is poor memory access for wann_func and
      ! but the reduced number of operations wins out.
      do nzz = -rgrid%n3, 2*rgrid%n3 - 1
        nz = mod(nzz, rgrid%n3)
        if(nz .lt. 1) nz = nz + rgrid%n3
        do nyy = -rgrid%n2, 2*rgrid%n2 - 1
          ny = mod(nyy, rgrid%n2)
          if(ny .lt. 1) ny = ny + rgrid%n2
          do nxx = -rgrid%n1, 2*rgrid%n1 - 1
            nx = mod(nxx, rgrid%n1)
            if(nx .lt. 1) nx = nx + rgrid%n1

            scalfac = kpts(1, ik)*real(nxx - 1, dp)/real(rgrid%n1, dp) + &
                      kpts(2, ik)*real(nyy - 1, dp)/real(rgrid%n2, dp) + &
                      kpts(3, ik)*real(nzz - 1, dp)/real(rgrid%n3, dp)
            npoint = nx + (ny - 1)*rgrid%n1 + (nz - 1)*rgrid%n2*rgrid%n1
            catmp = exp(2.0*pi*ii*scalfac)
            do ib = 1, num_wann
              do iw = 1, num_wann
                WannFunc(nxx, nyy, nzz, iw, idisp, is) = &
                  WannFunc(nxx, nyy, nzz, iw, idisp, is) + &
                  u_matrix(ib, iw, ik)*r_wvfn(npoint, ib, is)*catmp
              end do
            end do
          end do
        end do

      end do

    end do
    end do
    end do

    ! fix the global phase by setting the wannier to
    ! be real at the point where it has max. modulus
    do is = 1, SpinDim
    do idisp = 1, 7
    do iw = 1, num_wann
      tmaxx = 0.0
      wmod = cmplx(1.0_dp, 0.0_dp)
      do nzz = -rgrid%n3, 2*rgrid%n3 - 1
        do nyy = -rgrid%n2, 2*rgrid%n2 - 1
          do nxx = -rgrid%n1, 2*rgrid%n1 - 1
            WannFunc(nxx, nyy, nzz, iw, idisp, is) = WannFunc(nxx, nyy, nzz, iw, idisp, is)/real(num_kpts, dp)
            tmax = real(WannFunc(nxx, nyy, nzz, iw, idisp, is)* &
                        Conjg(WannFunc(nxx, nyy, nzz, iw, idisp, is)), dp)
            if (tmax > tmaxx) then
              tmaxx = tmax
              wmod = WannFunc(nxx, nyy, nzz, iw, idisp, is)
            end if
          end do
        end do
      end do
      wmod = wmod/sqrt(real(wmod)**2 + aimag(wmod)**2)
      WannFunc(:, :, :, iw, idisp, is) = WannFunc(:, :, :, iw, idisp, is)/wmod
    end do
    end do
    end do

    ! Check the 'reality' of the WF
    do is = 1, SpinDim
    do idisp = 1, 7
    do iw = 1, num_wann
      ratmax = 0.0_dp
      do nzz = -rgrid%n3, 2*rgrid%n3 - 1
        do nyy = -rgrid%n2, 2*rgrid%n2 - 1
          do nxx = -rgrid%n1, 2*rgrid%n1 - 1
            if (abs(real(WannFunc(nxx, nyy, nzz, iw, idisp, is), dp)) >= 0.01_dp) then
              ratio = abs(aimag(WannFunc(nxx, nyy, nzz, iw, idisp, is)))/ &
                      abs(real(WannFunc(nxx, nyy, nzz, iw, idisp, is), dp))
              ratmax = max(ratmax, ratio)
            end if
          end do
        end do
      end do
      write (*, '(6x,a,i4,7x,a,f11.6)') 'Wannier Function Num: ', iw, &
        'Maximum Im/Re Ratio = ', ratmax
    end do
    end do
    end do
    return

  contains
    Subroutine ReadCHK(num_kpts,num_bands,num_wann,Qdisentangled,kpts, &
                       u_matrix,u_matrix_opt,m_matrix,lwindow,ndimwin)
      implicit none
      integer,intent(out)                 :: num_kpts, num_wann, num_bands
      complex(dp),allocatable,intent(out) :: u_matrix(:,:,:), u_matrix_opt(:,:,:), m_matrix(:,:,:,:)
      logical,allocatable,intent(out)     :: lwindow(:, :)
      integer,allocatable,intent(out)     :: ndimwin(:)
      real(dp),allocatable,intent(out)    :: kpts(:,:)
      real(dp),allocatable                :: exclude_bands(:)
      real(dp),allocatable                :: wannier_centres(:,:), wannier_spreads(:)
      real(dp)                            :: latt_a(3,3), latt_b(3,3)
      real(dp)                            :: omega_invariant
      integer                             :: chk=2090, i, j, k, l, nkp, ierr, idisp
      integer                             :: nntot, mp_grid(3), num_exclude_bands
      logical                             :: Qdisentangled
      character(20)                       :: checkpoint
      character(33)                       :: header
      character(10)                       :: prefix(7)


      prefix=['disp0', 'disp+x', 'disp+y', 'disp+z', 'disp-x', 'disp-y', 'disp-z']
 
!      Do idisp = 1, 7 
        open(unit=chk, file='wannier90.chk', status='old',form='unformatted')
        ! Read comment line
        read(chk) header
  
        ! Consistency checks
        read(chk) num_bands                                  ! Number of bands
        read(chk) num_exclude_bands                     ! Number of excluded bands
        allocate(exclude_bands(num_exclude_bands), stat=ierr)
        read(chk)(exclude_bands(i), i=1, num_exclude_bands)  ! Excluded bands
        read(chk) ((latt_a(i, j), i=1, 3), j=1, 3)     ! Real lattice
        read(chk) ((latt_b(i, j), i=1, 3), j=1, 3)    ! Reciprocal lattice
        read(chk) num_kpts                                   ! K-points
        read(chk) (mp_grid(i), i=1, 3)                       ! M-P kgrid
        allocate(kpts(3, num_kpts), stat=ierr)
        read(chk) ((kpts(i, nkp), i=1, 3), nkp=1, num_kpts)
        read(chk) nntot                                      ! nntot
        read(chk) num_wann                                   ! num_wann
        read(chk) checkpoint                                 ! checkpoint
        checkpoint = adjustl(trim(checkpoint))
  
        read(chk) Qdisentangled                          ! whether a disentanglement has been performed
        if (Qdisentangled) then
          read(chk) omega_invariant                          ! omega invariant
  
          ! lwindow
          allocate(lwindow(num_bands, num_kpts), stat=ierr)
          read(chk) ((lwindow(i, nkp), i=1, num_bands), nkp=1, num_kpts)
  
          ! ndimwin
          allocate(ndimwin(num_kpts), stat=ierr)
          read(chk) (ndimwin(nkp), nkp=1, num_kpts)
  
          ! U_matrix_opt
          allocate(u_matrix_opt(num_bands, num_wann, num_kpts), stat=ierr)
          read(chk) (((u_matrix_opt(i, j, nkp), i=1, num_bands), j=1, num_wann), nkp=1, num_kpts)
        endif
  
        ! U_matrix
        allocate(u_matrix(num_wann, num_wann, num_kpts), stat=ierr)
        read(chk) (((u_matrix(i, j, k), i=1, num_wann), j=1, num_wann), k=1, num_kpts)
  
        ! M_matrix
        allocate (m_matrix(num_wann, num_wann, nntot, num_kpts), stat=ierr)
        read(chk) ((((m_matrix(i, j, k, l), i=1, num_wann), j=1, num_wann), k=1, nntot), l=1, num_kpts)
  
        ! wannier_centres
        allocate(wannier_centres(3, num_wann), stat=ierr)
        read(chk) ((wannier_centres(i, j), i=1, 3), j=1, num_wann)
  
        ! wannier spreads
        allocate(wannier_spreads(num_wann), stat=ierr)
        read(chk) (wannier_spreads(i), i=1, num_wann)
  
        close (chk)
!      end do
  
    end subroutine ReadCHK

  end subroutine

  Function ReadUnk(filename,num_bands,ngtot) Result(Unkr)
    implicit none
    character(72),intent(in)            :: filename
    Integer,intent(in)                  :: ngtot, num_bands
    Complex(dp)                         :: Unkr(ngtot,num_bands)
    Logical                             :: file_exists
    Integer                             :: ngx, ngy, ngz, nbnd
    Integer                             :: nx, nk, ib
    Integer                             :: unk=2020

    INQUIRE(FILE=trim(filename), EXIST=file_exists)
    if(file_exists) then
      open(unit=unk, file=trim(filename), form='unformatted')
      read(unk) ngx, ngy, ngz, nk, nbnd
      read(unk) ((Unkr(nx, ib), nx=1, ngx*ngy*ngz), ib=1,nbnd)
      close(unk)
    else
      write(*,*) "need "//filename
      call abort
    end if

  end Function ReadUnk

End Module Wannier
