Module Aux

  Use Parameters
  Use Lattice

  Implicit none
  Contains

  Subroutine LinearIndexing(id, ix, iy, iz, ngx, ngy, ngz)
    Implicit none
    Integer, Intent(in)  :: id, ngx, ngy, ngz
    Integer, Intent(out) :: ix, iy, iz
    
    ix = Mod(id-1, ngx) + 1
    iy = Ceiling(Real(Mod(id-1, ngx*ngy)+1)/Real(ngx))
    iz = Ceiling(Real(id)/Real(ngx*ngy))

  End Subroutine LinearIndexing

  Function int2str5(i) Result(str)
    integer, intent(in) :: i
    character(len=20)   :: str

    write (str, '(I5.5)') i
    str = trim(str)
  end function int2str5

  function int2str(i) result(str)
    implicit none
    character(20)       :: str
    integer, intent(in) :: i

    write (str, *) i
    str = trim(adjustl(str))
  end function int2str

  Function FieldTo1D(Fields, ifield) Result(Field1D)
    Implicit none
    Integer, Intent(in)    :: ifield
    Real*8,  Intent(in)    :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8                 :: Field1D(FieldDim*NumField*(cgrid%npts))
    Integer                :: ix, iy, iz, i, id

    do iz = 1, cgrid%n3
      do iy = 1, cgrid%n2
        do ix = 1, cgrid%n1
          do i = 1, FieldDim
            id = (iz-1)*cgrid%n2*cgrid%n1*FieldDim &
            + (iy-1)*cgrid%n1*FieldDim &
            + (ix-1)*FieldDim &
            + i
            Field1D(id) = Fields(i,ifield,ix,iy,iz)
          end do
        end do
      end do
    end do

  End Function FieldTo1D

  Function GetSysVector(Fields,e0ij) Result(vector)
    Implicit none
    Real*8,  Intent(in)    :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8,  Intent(in)    :: e0ij(3,3)
    Real*8                 :: vector(optdim)
    Real*8                 :: eta(6)
    Integer                :: ix, iy, iz, ifield, i, id

    id = 0

    do iz = 1, cgrid%n3
      do iy = 1, cgrid%n2
        do ix = 1, cgrid%n1
          do ifield = 1, NumField
            do i = 1, FieldDim
              If (ix.le.cgrid%n1.and.iy.le.cgrid%n2.and.iz.le.cgrid%n3) then
                id = id + 1
                vector(id) = Fields(i,ifield,ix,iy,iz)
              end if
            end do
          end do
        end do
      end do
    end do

    eta = eij2eta(e0ij)
    do i = 1, 6
      vector(FieldDim*NumField*cgrid%npts+i) = eta(i)
    end do

  End Function GetSysVector

  Subroutine SysVector2Fields(sysvector, Fields, e0ij)
    Implicit none
    Real*8,  Intent(in)    :: sysvector(optdim)
    Real*8,  Intent(out)   :: e0ij(3,3)
    Real*8,  Intent(out)   :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8                 :: tmp(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Integer                :: ix, iy, iz, ifield, i, id

    tmp = Reshape(sysvector(1:FieldDim*NumField*cgrid%npts), &
                            (/FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3/))

    do iz = 1, cgrid%n3
      do iy = 1, cgrid%n2
        do ix = 1, cgrid%n1
          do ifield = 1, NumField
            do i = 1, FieldDim
              If (ix.le.cgrid%n1.and.iy.le.cgrid%n2.and.iz.le.cgrid%n3) then
                Fields(i, ifield, ix, iy, iz) = tmp(i, ifield, ix, iy, iz)
              else
                Fields(i, ifield, ix, iy, iz) = 0.0_dp
              end if
            end do
          end do
        end do
      end do
    end do
    e0ij = eta2eij(sysvector(optdim-5:optdim))

  End Subroutine SysVector2Fields

  Subroutine CMatSymmetrize(Mat1, Mat2, ndim, prec)
    Implicit none
    Integer,intent(in)        :: ndim
    Real(dp),intent(in)       :: prec
    Complex(dp),intent(inout) :: Mat1(ndim,ndim), Mat2(ndim,ndim)
    Complex(dp)               :: M1(ndim,ndim), M2(ndim,ndim)

      M1 = 0.5*(Mat1+Mat2)
      M2 = 0.5*(Mat1-Mat2)
      Call CMatChop(M1,NumWann,prec)
      Call CMatChop(M2,NumWann,prec)
      Mat1 = M1+M2
      Mat2 = M1-M2

  End Subroutine

  function roundup(x,ndigit) Result(rounded)
    Implicit none
    Real(dp),intent(in) :: x
    Integer,intent(in)  :: ndigit
    Real(dp)            :: rounded

    rounded = real(nint(x * 10**ndigit),dp) / 10**ndigit

  end function roundup

  subroutine MatRoundup(Mat,ndim,ndigit)
    Implicit none
    Integer,intent(in)        :: ndim,ndigit
    Real(dp),intent(inout) :: Mat(ndim,ndim)
    Integer                   :: i,j

    Do i = 1, ndim
      Do j = 1, ndim
        Mat(i,j) = roundup(Mat(i,j),ndigit)
      End do
    End do

  end subroutine MatRoundup

  subroutine CMatRoundup(Mat,ndim,ndigit)
    Implicit none
    Integer,intent(in)        :: ndim,ndigit
    Complex(dp),intent(inout) :: Mat(ndim,ndim)
    Real(dp)                  :: x_real, x_imag
    Integer                   :: i,j

    Do i = 1, ndim
      Do j = 1, ndim
        x_real = roundup(Real(Mat(i,j)),ndigit)
        x_imag = roundup(Aimag(Mat(i,j)),ndigit)
        Mat(i,j) = cmplx(x_real,x_imag)
      End do
    End do

  end subroutine CMatRoundup

  function chop(x,prec) Result(choped)
    Implicit none
    Real(dp),intent(in) :: x
    Real(dp),intent(in) :: prec
    Real(dp)            :: choped

    choped = x
    if(Abs(x).lt.prec) choped = 0.0_dp

  end function chop

  subroutine MatChop(Mat,ndim,prec)
    Implicit none
    Integer,intent(in)        :: ndim
    Real(dp),intent(in)       :: prec
    Real(dp),intent(inout)    :: Mat(ndim,ndim)
    Integer                   :: i,j

    Do i = 1, ndim
      Do j = 1, ndim
        Mat(i,j) = chop(Mat(i,j),prec)
      End do
    End do

  end subroutine MatChop

  subroutine CMatChop(Mat,ndim,prec)
    Implicit none
    Integer,intent(in)        :: ndim
    Real(dp),intent(in)       :: prec
    Complex(dp),intent(inout) :: Mat(ndim,ndim)
    Integer                   :: i,j
    Real(dp)                  :: x_real, x_imag

    Do i = 1, ndim
      Do j = 1, ndim
        x_real = Chop(Real(Mat(i,j)),prec)
        x_imag = Chop(Aimag(Mat(i,j)),prec)
        Mat(i,j) = cmplx(x_real,x_imag)
      End do
    End do

  end subroutine CMatChop

  function mesh_mirror_x(V,xdim,ydim,zdim,ind) result(V_out)
    Implicit none
    Integer,intent(in)      :: xdim, ydim, zdim, ind
    Real(dp),intent(in)     :: V(xdim,ydim,zdim)
    Real(dp)                :: V_out(xdim,ydim,zdim)
    Integer                 :: ix, iy, iz

    if(ind.eq.1) then
      do iz = 1, zdim
        do iy = 1, ydim
          do ix = 2, xdim
            V_out(ix,iy,iz) = V(xdim-ix+2, iy, iz)
          end do
          V_out(1,iy,iz) = V(1, iy, iz)
        end do
      end do
    else if(ind.eq.2) then
      do iz = 1, zdim
        do iy = 2, ydim
          do ix = 1, xdim
            V_out(ix,iy,iz) = V(ix, ydim-iy+2, iz)
          end do
        end do
        V_out(:,1,iz) = V(:,1,iz)
      end do
    else
      do iz = 2, zdim
        do iy = 1, ydim
          do ix = 1, xdim
            V_out(ix,iy,iz) = V(ix, iy, zdim-iz+2)
          end do
        end do
      end do
      V_out(:,:,1) = V(:,:,1)
    end if

  End function mesh_mirror_x

  Function mesh_mirror_xy(V,xdim,ydim,zdim,ind) result(V_out)
    Implicit none
    Integer,intent(in)      :: xdim, ydim, zdim, ind
    Real(dp),intent(in)     :: V(xdim,ydim,zdim)
    Real(dp)                :: V_out(xdim,ydim,zdim)
    Integer                 :: ix, iy, iz

    if(ind.eq.3) then
      do iz = 1, zdim
        do iy = 1, ydim
          do ix = 1, xdim
            V_out(ix,iy,iz) = V(iy, ix, iz)
          end do
        end do
      end do
    else if(ind.eq.2) then
      do iz = 1, zdim
        do iy = 1, ydim
          do ix = 1, xdim
            V_out(ix,iy,iz) = V(iz, iy, ix)
          end do
        end do
      end do
    else
      do iz = 1, zdim
        do iy = 1, ydim
          do ix = 1, xdim
            V_out(ix,iy,iz) = V(ix, iz, iy)
          end do
        end do
      end do
    end if

  End function mesh_mirror_xy

  Function lower(s) result(lower_s)
    implicit none
    character(*), intent(in) :: s
    character(len(s)) :: lower_s
    integer :: i

    do i = 1, len(s)
      if (('A' <= s(i:i)) .and. (s(i:i) <= 'Z')) then
        lower_s(i:i) = char(ichar(s(i:i)) + 32)
      else
        lower_s(i:i) = s(i:i)
      end if
    end do

    return
  End function lower

  function ArrayFlatten2D(A) result(C)
    real*8              ::  A(:,:)
    real*8,allocatable  ::  C(:)

    C = Reshape(A, [Size(A)])

  end function ArrayFlatten2D

 
End Module Aux
