Subroutine GetInitConfig(Fields, e0ij, inittype, FileCode)
  
  Implicit none
  character(20)                 :: keyword
  character(50)                 :: filename
  character(40), intent(inout)  :: inittype
  character(10), intent(in)     :: FileCode
  real*8,        intent(inout)  :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  real*8,        intent(inout)  :: e0ij(3,3)
  real*8                        :: eta(6)
  real*8                        :: field(FieldDim,NumField)
  real*8                        :: rand
  real(dp)                      :: qdotr
  integer                       :: i, j, ix, iy, iz, icell, ifield
  Integer                       :: FileHandle = 1111
  logical                       :: file_exists
   
  
  filename = trim(Solver)//".out/"//trim(inittype)//"-"//trim(FileCode)//".dat"
  call caps2small(inittype)
  keyword=trim(inittype)
       
  Fields = 0.0_dp

  SELECT CASE (keyword)
    CASE ("random")
      do icell = 1, cgrid%ncells
        ix = GridFold(1, icell)
        iy = GridFold(2, icell)
        iz = GridFold(3, icell)
        do i = 1, NumField
          do j = 1, FieldDim
            Call random_number(rand)
            Fields(j,i,ix,iy,iz) = FieldsBinary(j,i)*(rand-0.5d0)/1.0d1
          end do
        end do
      end do
      Fields(:,NumField,:,:,:) =0.0D0

      do i = 1, 6
      Call random_number(rand)
        eta(i) = (rand-0.50d0)/1.0d3
      end do
      e0ij = eta2eij(eta)
    CASE ("zero")
      Fields = 0.0d0
      e0ij = 0.0d0
    
    CASE ("phase")
      INQUIRE(FILE=trim(Solver)//".out/"//"phase-"//trim(FileCode)//".dat", EXIST=file_exists)
      if (file_exists .eqv. .False.) then
        write(*,*) "need phase-xxxxx.dat file"
        call abort
      else
        open(FileHandle,file=trim(Solver)//".out/"//"phase-"//trim(FileCode)//".dat",form='formatted',status='old')
        Read(FileHandle, "(3E25.15)") e0ij(1,1), e0ij(2,2), e0ij(3,3)
        Read(FileHandle, "(3E25.15)") e0ij(2,3), e0ij(1,3), e0ij(1,2)
        e0ij(3,2) = e0ij(2,3)
        e0ij(3,1) = e0ij(1,3)
        e0ij(2,1) = e0ij(1,2)
    
        do ifield = 1, NumField
          Read(FileHandle, "(3E25.15)") (field(i, ifield), i=1,FieldDim)
          do i = 1, FieldDim
            do iz = 1, cgrid%n3
              do iy = 1, cgrid%n2
                do ix =1, cgrid%n1
                  qdotr=Real(Exp(2*pi*cmplx(0.0_dp,1.0_dp)*(DWq(1,i)*ix+DWq(2,i)*iy+DWq(3,i)*iz)))
                  Fields(i,ifield,ix,iy,iz) = FieldsBinary(i,ifield)*field(i,ifield)*tanh(100.0*qdotr)
                end do
              end do
            end do
          end do
        end do
        Close(FileHandle)
      end if
    CASE DEFAULT
      INQUIRE(FILE=filename, EXIST=file_exists)
      if (file_exists .eqv. .False.) then
        write(*,*) "Init not implemented yet!"
        call abort
      else
        write(*,*) "Init from file: ", filename
        call InitFromFile(filename, Fields, e0ij)
      end if
    END SELECT
  
End Subroutine GetInitConfig

