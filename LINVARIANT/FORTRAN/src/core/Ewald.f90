Module Ewald

  Use omp_lib
  Use Parameters
  Use Constants
  Use Lattice
  Use Aux

  Implicit none
  Contains

  Include "GetEOnSiteEwald.f90"
  Include "GetEwaldField.f90"

  Subroutine GetEwaldHessian(EwaldMatrix, EwaldHessian)

    implicit none
    real*8,  intent(in)           :: EwaldMatrix(3, 3, cgrid%n1, cgrid%n2, cgrid%n3)
    real*8,  intent(inout)        :: EwaldHessian(3*cgrid%npts,3*cgrid%npts,NumField)
    integer                       :: i, j, k, ifield, ix, iy, iz, jx, jy, jz, x, y, z

    Do ifield = 1, NumField
    Do i = 1, cgrid%npts
    Do j = 1, cgrid%npts
      Call LinearIndexing(i, ix, iy, iz, cgrid%n1, cgrid%n2, cgrid%n3)
      Call LinearIndexing(j, jx, jy, jz, cgrid%n1, cgrid%n2, cgrid%n3)
      Call Pbc(jx-ix, jy-iy, jz-iz, x, y, z)
      EwaldHessian((i-1)*3+1:(i-1)*3+3, (j-1)*3+1:(j-1)*3+3, ifield) = &
             EwaldMatrix(:,:,x,y,z)*FieldCharge(ifield)**2
      If(i.eq.j) then
        Do k = 1, 3 
        EwaldHessian((i-1)*3+k, (j-1)*3+k, ifield) = &
           2.0D0*EwaldHessian((i-1)*3+k, (j-1)*3+k, ifield)
        End do
      End if
    End Do
    End Do 
    End Do

  End Subroutine GetEwaldHessian

  Function GetEwaldForces(Fields) Result(EwaldForce)
    
    implicit none
    real*8,  intent(in)           :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    integer                       :: i, ifield
    integer                       :: j, jfield, jx, jy, jz
    integer                       :: ix, iy, iz
    integer                       :: x, y, z
    real*8                        :: EwaldForce(Max(FieldDim, 6),NumField+1, cgrid%n1, cgrid%n2, cgrid%n3)
    real*8                        :: force(3*cgrid%npts,NumField - 1)
    real*8                        :: vector(3*cgrid%npts,NumField - 1)
    
    EwaldForce = 0.0D0
    force = 0.0D0
  
    do jfield = 1, NumField - 1
      vector(:,jfield) = FieldTo1D(Fields, jfield)
      call dgemv('N', 3*cgrid%npts, 3*cgrid%npts, &
                 -1.0D0, EwaldHessian(:,:,jfield), 3*cgrid%npts, &
                 vector(:,jfield), 1, 0.0D0, force(:,jfield), 1)
      EwaldForce(1:3,jfield,:,:,:) = Reshape(force(:,jfield), (/3, cgrid%n1, cgrid%n2, cgrid%n3/))
    end do
  
  End Function GetEwaldForces

  Function GetDeltaHEwald(x0, y0, z0, EwaldField, idelta, delta) Result(ene)
    
    Implicit none
    Integer, Intent(in) :: x0, y0, z0, idelta
    Real*8,  Intent(in) :: EwaldField(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8,  Intent(in) :: delta(FieldDimList(idelta))
    Integer             :: i, j, ifield, jfield
    Real*8              :: ene, ThreadEne
    
    ene = 0.0d0
    
    do i = 1, FieldDimList(idelta)
      ene = ene + delta(i)*EwaldField(i,idelta,x0,y0,z0)
    end do ! i
    
    do i = 1, FieldDimList(idelta)
      do j = 1, FieldDimList(idelta)
        ene = ene + 0.5*FieldCharge(idelta)**2*delta(i)*EwaldMat(j,i,1,1,1)*delta(j)
      end do ! j
    end do ! i
    
  End Function GetDeltaHEwald

  Subroutine UpdateEwaldField(ix, iy, iz, EwaldField, idelta, delta)

    Implicit none
    integer, intent(in)     :: ix, iy, iz, idelta
    real*8,  intent(inout)  :: EwaldField(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    real*8,  intent(in)     :: delta(FieldDimList(idelta))
    real*8                  :: dfield(cgrid%npts*3)
    integer                 :: id

    dfield = 0.0D0
    id = (iz-1)*cgrid%n2*cgrid%n1*FieldDimList(idelta) &
              + (iy-1)*cgrid%n1*FieldDimList(idelta) &
              + (ix-1)*FieldDimList(idelta)

    call dgemv('N', 3*cgrid%npts, 3, &
               1.0D0, EwaldHessian(:,id+1:id+3,idelta), 3*cgrid%npts,&
               delta, 1, 0.0D0, dfield, 1)
    EwaldField(:,idelta,:,:,:) = EwaldField(:,idelta,:,:,:) + Reshape(dfield, (/3, cgrid%n1, cgrid%n2, cgrid%n3/))

  End Subroutine UpdateEwaldField

  Subroutine EwaldMatrix
    
    Use Parameters
    Use Constants
    Use Lattice
    Use Aux
    
    Implicit none
    Real*8             :: dpij(3,3,cgrid%n1,cgrid%n2,cgrid%n3), dum(3,3), superlatt(3,3), recalat(3,3), am(3)
    Real*8             :: pi2, celvol, eta, eta4, gcut, gcut2, dum0
    Real*8             :: gx, gy, gz, g2, fact, tol, c, residue
    Character(len=10)  :: NameEwald
    Integer            :: N1, N2, N3
    Integer            :: i, j, mg1, mg2, mg3, ix, iy, iz, ig1, ig2, ig3, rx, ry, rz
    
    Logical            :: origin, file_exists, file_match
  
    file_exists = .false.
    file_match = .false.
  
    INQUIRE(FILE="EwaldMat.dat", EXIST=file_exists)
    if (file_exists) then
      write(*,*) 'EwaldMat.dat exist, checking ...'
      open(ifileno,file='EwaldMat.dat',form='unformatted',status='old')
        read(ifileno) NameEwald, N1, N2, N3
      close(ifileno)
      file_match = N1.eq.cgrid%n1 .and. N2.eq.cgrid%n2 .and. N3.eq.cgrid%n3 .and. NameEwald .eq. NameSim
    end if
      
    if (file_exists .and. file_match) then
      write(*,*) 'Dipole file matches and will be used.'
      open(ifileno,file='EwaldMat.dat',form='unformatted',status='old')
        read(ifileno) NameEwald, N1, N2, N3
        read(ifileno) (((((EwaldMat(j,i,ix,iy,iz),j=1,3),i=1,3),ix=1,cgrid%n1),iy=1,cgrid%n2),iz=1,cgrid%n3)
      close(ifileno)
    else
      if (file_exists) then
        write(*,*) 'Dipole file mismatch, current simulation: ', NameSim, cgrid%n1, cgrid%n2, cgrid%n3
        write(*,*) 'Ewald file gives: ', NameEwald, N1, N2, N3
        write(*,*) 'Now, regenerate new Ewald Matrix!!!'
      else
        write(*,*) 'No EwaldMat.dat is found, generating ...'
      end if
  
      !!!!!!!!!!!!!!!!!!
      ! some constants !
      !!!!!!!!!!!!!!!!!!
      
      pi2 = pi*2.0d0
      
      tol = 1.0d-12
      eta = sqrt(-log(tol))
      gcut = 2.0*eta**2
      gcut2 = gcut**2
      eta4 = 1.0d0/(4*eta**2)
      
      superlatt(:,1) = alat(:,1)*cgrid%n1
      superlatt(:,2) = alat(:,2)*cgrid%n2
      superlatt(:,3) = alat(:,3)*cgrid%n3
      
      call reclat(superlatt,recalat,1)
      celvol = Cell2Volume(superlatt)
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Estimate number of reciprocal alatice vectors to use !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      do i = 1, 3
        am(i) = 0.d0
        do j = 1, 3
          am(i) = am(i) + superlatt(j,i)**2
        enddo
        am(i) = sqrt(am(i))
      enddo
      
      mg1 = int(gcut*am(1)/pi2) + 1
      mg2 = int(gcut*am(2)/pi2) + 1
      mg3 = int(gcut*am(3)/pi2) + 1
      
      write(6,*) 'Simulation name: ', NameSim
      write(6,*) 'Simulation cell volume:', celvol
      write(6,*) 'Simulation Grid:', cgrid%n1, cgrid%n2, cgrid%n3
      write(6,*) 'Gcut: ', gcut, 'mg1,mg2,mg3: ', mg1,mg2,mg3
      write(6,*) 'Tol: ', tol
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Begin the calculation of dpij matrix !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      do ix = 1, cgrid%n1
        do iy = 1, cgrid%n2
          do iz = 1, cgrid%n3
          
            do i = 1, 3
              do j = 1, 3
                dum(i,j) = 0.0d0
              end do
            end do
            
            rx = ix - 1
            ry = iy - 1
            rz = iz - 1
            
            origin = (rx.eq.0) .and. (ry.eq.0) .and. (rz.eq.0)
            
            !
            !The term for -G gives the same contribution as the one for G,
            !so we only sum over half the ig1 range and multiply by 2 (exclude
            !ig1=0)
            !
            !Since we are really computing the "dipole field", 
            !not the total dipole energy, we need another factor of 2:
            !
            !Energy = Sum_ij{ Q_ij u_i u_j}, 
            !Field_k= Partial Energy/Partial u_k = 2*Sum_i{Q_ik u_i} 
            !
            
            c = 8.0d0*pi/celvol
            residue = eta**3*(4.0d0/(3.0d0*sqrt(pi)))
            
            do ig1 = 0, mg1
              do ig2 = -mg2, mg2
                gx = ig1*recalat(1,1) + ig2*recalat(1,2) 
                gy = ig1*recalat(2,1) + ig2*recalat(2,2) 
                do ig3 = -mg3, mg3
                
                  gz = ig3*recalat(3,3)
                  g2 = gx**2 + gy**2 + gz**2
                  
                  if(g2.lt.gcut2.and.g2.gt.1.0d-8) then
                    fact = 1.0d0
                    if(ig1.eq.0) fact = 0.50d0
                    dum0 = fact*cos(gx*rx+gy*ry+gz*rz)*exp(-g2*eta4)/g2
                    dum(1,1) = dum(1,1) + dum0*gx**2
                    dum(2,2) = dum(2,2) + dum0*gy**2
                    dum(3,3) = dum(3,3) + dum0*gz**2
                    dum(2,3) = dum(2,3) + dum0*gz*gy
                    dum(1,3) = dum(1,3) + dum0*gx*gz
                    dum(1,2) = dum(1,2) + dum0*gx*gy
                    
                    dum(3,2) = dum(2,3)
                    dum(3,1) = dum(1,3)
                    dum(2,1) = dum(1,2)
                  end if
                end do ! ig3=-mg3,mg3
              end do ! ig2=-mg2,mg2
            end do ! ig1=0,mg1
                
            do i = 1, 3
              do j = 1, 3
                dpij(i, j, ix, iy, iz) = dum(i,j)*c
              end do
            end do
            if (origin) then
              dpij(1, 1, ix, iy, iz) = dpij(1, 1, ix, iy, iz) - residue
              dpij(2, 2, ix, iy, iz) = dpij(2, 2, ix, iy, iz) - residue
              dpij(3, 3, ix, iy, iz) = dpij(3, 3, ix, iy, iz) - residue
            endif
          end do ! ix
        end do ! iy
      end do ! iz
          
      open(1, file='EwaldMat.dat', form='unformatted', status='unknown')
      write(1) NameSim, cgrid%n1, cgrid%n2, cgrid%n3
      write(1) (((((dpij(j, i, ix,iy,iz),j=1,3),i=1,3),ix=1,cgrid%n1),iy=1,cgrid%n2),iz=1,cgrid%n3)
      close(1)
      
      open(2, file='EwaldMat.mma.dat', form='formatted', status='unknown')
      write(2, '(3e25.14)') (((((dpij(j, i, ix,iy,iz),j=1,3),i=1,3),ix=1,cgrid%n1),iy=1,cgrid%n2),iz=1,cgrid%n3)
      close(2)
    end if
  
  End Subroutine EwaldMatrix
End Module Ewald
