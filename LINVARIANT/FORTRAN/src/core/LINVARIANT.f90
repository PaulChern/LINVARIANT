Module LINVARIANT
  
  Use Parameters
  Use Constants
  Use Energy
  Use Lattice
  Use Aux
  
  Implicit none
  Contains
  
  Subroutine RemoveGridDrifts(Fields)
    Implicit none
    Real*8,  Intent(inout) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8                 :: GridDriftx, GridDrifty, GridDriftz
    Integer                :: ix, iy, iz
    
    GridDriftx = Sum(Fields(1,NumField,:,:,:))/(cgrid%npts)
    GridDrifty = Sum(Fields(2,NumField,:,:,:))/(cgrid%npts)
    GridDriftz = Sum(Fields(3,NumField,:,:,:))/(cgrid%npts)
    
    do iz = 1, cgrid%n3
      do iy = 1, cgrid%n2
        do ix = 1, cgrid%n1
          Fields(1,NumField,ix,iy,iz) = Fields(1,NumField,ix,iy,iz) - GridDriftx
          Fields(2,NumField,ix,iy,iz) = Fields(2,NumField,ix,iy,iz) - GridDrifty
          Fields(3,NumField,ix,iy,iz) = Fields(3,NumField,ix,iy,iz) - GridDriftz
        end do
      end do
    end do
    
  End Subroutine RemoveGridDrifts
  
  Function AllFields1D(Fields) Result(Fields1D)
    Implicit none
    Real*8,  Intent(in)    :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8                 :: Fields1D(FieldDim*NumField*cgrid%npts)
    Integer                :: ix, iy, iz, ifield, i, id
    
    do iz = 1, cgrid%n3
      do iy = 1, cgrid%n2
        do ix = 1, cgrid%n1
          do ifield = 1, NumField
            do i = 1, FieldDim
              id = (iz-1)*cgrid%n2*cgrid%n1*NumField*FieldDim &
              + (iy-1)*cgrid%n1*NumField*FieldDim &
              + (ix-1)*NumField*FieldDim &
              + (ifield-1)*FieldDim &
              + i
              Fields1D(id) = Fields(i,ifield,ix,iy,iz)
            end do
          end do
        end do
      end do
    end do
    
  End Function AllFields1D

  Function AllFieldsFrom1D(Fields1D) Result(Fields)
    Implicit none
    Real*8,  Intent(in)    :: Fields1D(FieldDim*NumField*cgrid%npts)
    Real*8                 :: Fields(FieldDim,NumField,cgrid%n1,cgrid%n2,cgrid%n3)
    
    Fields = Reshape(Fields1D, (/FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3/))
    
  End Function AllFieldsFrom1D

  Function Thermometer(dFieldsdt, de0ijdt) Result(TK)
    Implicit none
    Real*8,  Intent(in) :: dFieldsdt(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8,  Intent(in) :: de0ijdt(3,3)
    Real*8              :: TK(2), Ekin(2)
    Integer             :: i, ifield, ix, iy, iz, TotalDim

    TotalDim = 0.0D0
    do i = 1, NumField
      If(.Not.CLAMPQ(i)) TotalDim = TotalDim+cgrid%npts*3
    end do
    
    Ekin = GetEkin(dFieldsdt, de0ijdt)
    TK(1) = 2*Ekin(1)*Hartree/k_bolt_ev/TotalDim
    TK(2) = 2*Ekin(2)*Hartree/k_bolt_ev/6
    
  End Function Thermometer

  Subroutine AdaptTk(dFieldsdt, de0ijdt)
    Implicit none
    Real*8,  Intent(inout) :: dFieldsdt(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8,  Intent(inout) :: de0ijdt(3,3)
    Real*8                 :: dFieldsdt_tmp(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8                 :: de0ijdt_tmp(3,3), detadt(6)
    Real*8                 :: TK, Ekin(NumField+1)
    Real*8                 :: rand
    Integer                :: i, ifield, ix, iy, iz, NumGrid

    NumGrid = cgrid%npts

    Ekin = ResolveEkin(dFieldsdt, de0ijdt)

    Do i = 1, NumField
      TK = 2*Ekin(i)*Hartree/k_bolt_ev/(FieldDimList(i)*NumGrid)
      dFieldsdt(:,i,:,:,:) = Sqrt(Temp/Tk)*dFieldsdt(:,i,:,:,:)
    End do
    
    TK = 2*Ekin(NumField+1)*Hartree/k_bolt_ev/6
    If(TK .gt. 0.0D0) then
      de0ijdt = Sqrt(Temp/Tk)*de0ijdt
    Else
      do i = 1, 6
        Call random_number(rand)
        detadt(i) = rand-0.50d0
      end do
      de0ijdt = eta2eij(detadt)
      Ekin = ResolveEkin(dFieldsdt, de0ijdt)
      TK = 2*Ekin(NumField+1)*Hartree/k_bolt_ev/6
      de0ijdt = Sqrt(Temp/Tk)*de0ijdt
    End If
    Write(*,'(A25,2F10.4)') "Initial Temperature Tk:", Thermometer(dFieldsdt, de0ijdt)

  End Subroutine AdaptTk

!  Subroutine GetACF(trajectory)
!    Implicit none
!    character(*), Intent(in) :: trajectory
!    Real*8,  Intent(out)     :: ACF(3,3,NumField,(NumSteps-ThermoSteps)/TapeRate-1)
!    Integer                  :: FileHandle = 1111
!    Integer                  :: ix, iy, iz, i, j, n, NN
!    Integer                  :: ifield, jfield, t, h
!    Real*8                   :: f1(cgrid%n1, cgrid%n2, cgrid%n3), f2
!    Real*8                   :: Fields(FieldDim,NumField,cgrid%n1,cgrid%n2,cgrid%n3,(NumSteps-ThermoSteps)/TapeRate)
!    Real*8                   :: e0ij(3,3,(NumSteps-ThermoSteps)/TapeRate)
!    Real*8                   :: ene(3,(NumSteps-ThermoSteps)/TapeRate)
!
!    NN = (NumSteps - ThermoSteps)/TapeRate
!    ACF = 0.0D0
!
!    open(FileHandle,file=trajectory,form='formatted',status='old')
!    do n = 1, NN
!      Read(FileHandle) ene(1,n), ene(2,n), ene(3,n)
!      Read(FileHandle) e0ij(1,1,n), e0ij(2,2,n), e0ij(3,3,n)
!      Read(FileHandle) e0ij(2,3,n), e0ij(1,3,n), e0ij(1,2,n)
!      do iz = 1, cgrid%n3
!      do iy = 1, cgrid%n2
!      do ix = 1, cgrid%n1
!        Read(FileHandle)
!        do ifield = 1, NumField
!          Read(FileHandle) (Fields(i, ifield, ix, iy, iz, n), i=1,FieldDimList(ifield))
!        end do
!      end do
!      end do
!      end do
!    end do
!    close(FileHandle)
!   
!    do ifield = 1, NumField
!      do i = 1, 3
!        do j = 1, 3
!          mu1 = Sum(Fields(i,ifield,:,:,:,:))/NN/(cgrid%npts)
!          mu2 = Sum(Fields(j,ifield,:,:,:,:))/NN/(cgrid%npts)
!          do h = 1, NN-1 
!            do t = 1, NN-h
!              ACF(i,j,ifield,h) = ACF(i,j,ifield,h) + Sum((Fields(i,ifield,:,:,:,t)-mu1)*(Fields(j,ifield,:,:,:,t+h)-mu2))
!            end do
!          end do
!        end do
!      end do
!    end do
!
!  End Subroutine

End Module LINVARIANT
