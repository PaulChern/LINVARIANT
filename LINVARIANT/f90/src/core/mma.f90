Module mma
  Use Constants
  Contains
  Function rot3(mat,v0) Result(v)
    Implicit None
    Real*8  :: mat(3,3), v0(3)
    Real*8  :: v(3)
    Integer :: i, j

    v=0.0D0
    Do i = 1, 3
      Do j =1,3
        v(i)=v(i)+v0(j)*mat(i,j)
      End do
    End do
  End Function rot3

  Subroutine LUInverse(mat,ndim)

    implicit none
    integer      :: ndim
    Complex(dp)  :: mat(ndim,ndim)
    Complex(dp)  :: work(ndim)
    integer      :: ipiv(ndim), info

    call ZGETRF(ndim,ndim,mat,ndim,ipiv,info)
    if(info.ne.0) stop "inverse_complex_matrix Error in ZGETRF"
    call ZGETRI(ndim,mat,ndim,ipiv,work,ndim,info)
    if(info.ne.0) stop "inverse_complex_matrix Error in ZGETRI"
 
  End subroutine LUInverse

  Subroutine Inverse(Amat,ndim)

    implicit none
    integer,Intent(in)         :: ndim
    Complex(dp),Intent(inout)  :: Amat(ndim,ndim)
    Complex(dp)                :: Bmat(ndim,ndim)
    integer                    :: i, ipiv(ndim), info=0

    ipiv = 0
    Bmat = cmplx(0.0_dp,0.0_dp)
    Do i = 1, ndim
      Bmat(i,i) = cmplx(1.0_dp,0.0_dp)
    End do

    call zgesv(ndim,ndim,Amat,ndim,ipiv,Bmat,ndim,info)
    if(info.ne.0) print *,'something wrong with linear solve zgesv'
    Amat = Bmat

  End subroutine Inverse

  Subroutine EigenSystem(JOBZ,UPLO,Hmat,ndim,EigenValues)
     implicit none

     !  JOBZ    (input) CHARACTER*1
     !          = 'N':  Compute eigenvalues only;
     !          = 'V':  Compute eigenvalues and eigenvectors.
     !
     !  UPLO    (input) CHARACTER*1
     !          = 'U':  Upper triangle of A is stored;
     !          = 'L':  Lower triangle of A is stored.
     Character*1,intent(in)    :: JOBZ
     Character*1,intent(in)    :: UPLO
     Integer,intent(in)        :: ndim 
     Complex(dp),intent(inout) :: Hmat(ndim,ndim)
     Real(dp),intent(inout)    :: EigenValues(ndim)
     integer                   :: info
     integer                   :: lwork
     Real(dp),allocatable      :: rwork(:)
     Complex(dp),allocatable   :: work(:)

     lwork=16*ndim
     allocate(work(lwork))
     allocate(rwork(lwork))
     rwork= 0.0_dp
     work= cmplx(0.0_dp, 0.0_dp)

     info=1
     EigenValues=0.0_dp

     !> if N==1, you don't have to do the diagonalization
     if (ndim.eq.1) then 
        EigenValues=Hmat(1, 1)
        Hmat(1, 1)= 1.0_dp
        return
     endif

     call zheev(JOBZ,UPLO,ndim,Hmat,ndim,EigenValues,work,lwork,rwork,info)

     if (info.ne.0) then
        write(*, *) 'ERROR : something wrong with zheev', info
        stop
     endif

     deallocate(rwork, work)
     return
  End subroutine EigenSystem

  Function D2Jij(vec) result(Jij)
    Implicit None
    Real(dp),intent(in)  :: vec(3)
    Real(dp)             :: Jij(3,3)

    Jij = cmplx(0.0_dp,0.0_dp)
    Jij(1,2) =  vec(3)
    Jij(2,1) = -vec(3)
    Jij(1,3) = -vec(2)
    Jij(3,1) =  vec(2)
    Jij(2,3) =  vec(1)
    Jij(3,2) = -vec(1)

  End function

  Function PauliDecompose(Mat,mdim,ndim) result(M4Array)

    Implicit None
    Integer,intent(in)       :: mdim,ndim
    Complex(dp),intent(in)   :: Mat(mdim,ndim)
    Complex(dp)              :: M4Array(mdim/2,ndim/2,4)

    M4Array(:,:,1)=(Mat(1::2,1::2)+Mat(2::2,2::2))/2.0
    M4Array(:,:,4)=(Mat(1::2,1::2)-Mat(2::2, 2::2))/2.0
    M4Array(:,:,2)=(Mat(1::2,2::2)+Mat(2::2, 1::2))/2.0
    M4Array(:,:,3)=(Mat(1::2,2::2)-Mat(2::2, 1::2))/(-2.0*ii)

  End Function PauliDecompose

  Function PauliNormDecompose(Mat,ndim) result(MatNorm)

    Implicit None
    Integer,intent(in)       :: ndim
    Complex(dp),intent(in)   :: Mat(ndim,ndim)
    Complex(dp)              :: M4Array(ndim/2,ndim/2,4)
    Complex(dp)              :: MatNorm(ndim/2,ndim/2),evec(3)
    Real(dp)                 :: cnorm
    Integer                  :: i

    M4Array = PauliDecompose(Mat,ndim,ndim)
    evec(1) = sum( (/ (M4Array(i,i,2), i=1,ndim/2) /) )
    evec(2) = sum( (/ (M4Array(i,i,3), i=1,ndim/2) /) )
    evec(3) = sum( (/ (M4Array(i,i,4), i=1,ndim/2) /) )
    cnorm = Norm2([Norm2(Real(evec)),Norm2(Aimag(evec))])
    
    MatNorm=(evec(1)*M4Array(:,:,2)+evec(2)*M4Array(:,:,3)+evec(3)*M4Array(:,:,4))/cnorm
  End Function PauliNormDecompose

  Function Trace(M,ndim) result(tr)
    Implicit none
    Complex(dp),intent(in) :: M(ndim,ndim)
    Integer,intent(in)     :: ndim
    Complex(dp)            :: tr
    Integer                :: i
   
    tr = sum( (/ (M(i,i), i=1, ndim) /) )

  End function Trace

  subroutine RectangleContour
    Use Parameters
    Implicit none
    Integer            :: i, nz, nz1, nz2, nz3
    Complex(dp)        :: zstep

    nz1=ContourNPoints(1)
    nz2=ContourNPoints(2)
    nz3=ContourNPoints(3)

    nz = nz1 + nz2 + nz3
    allocate(ContourPath(nz+1))
    ContourPath = cmplx(0.0,0.0,dp)
 
    ContourPath(1) = cmplx(ContourMin,0.d0,dp)
 
    zstep = cmplx(0.d0,(ContourHeight/nz1),dp)
    do i = 2, nz1+1
      ContourPath(i) = ContourPath(i-1) + zstep
    end do
 
    zstep = cmplx((ContourMax-ContourMin)/nz2,0.d0,dp)
    do i = nz1+2, nz1+nz2
      ContourPath(i) = ContourPath(i-1) + zstep
    end do
 
    zstep = -1.d0*cmplx(0.d0,(ContourHeight/nz3),dp)
    do i = nz1+nz2+1, nz+1
      ContourPath(i) = ContourPath(i-1) + zstep
    end do

  End subroutine RectangleContour

  Subroutine SemiCircleContour
    Use Parameters
    Implicit none
    Integer            :: i, nz, nz1, nz2, nz3
    Real(dp)           :: R0, R, phi
    Complex(dp)        :: zstep

    nz1=ContourNPoints(1)
    nz2=ContourNPoints(2)
    nz3=ContourNPoints(3)

    nz = nz1 + nz2 + nz3
    allocate(ContourPath(nz+1))

    R0= (ContourMin+ContourMax)/2.0
    R= (ContourMax-ContourMin)/2.0
    do i = 1, nz+1
      phi = pi - pi*(i-1)/(nz)
      ContourPath(i) = R0+R*exp(ii*phi)
    end do

  End subroutine SemiCircleContour

  
  function  GCD(a0,b0) result(b)
  ! ---------------------------------------------------------
  ! This program computes the GCD of two positive integers
  ! using the Euclid method.  Given a and b, a >= b, the
  ! Euclid method goes as follows:  (1) dividing a by b yields
  ! a reminder c; (2) if c is zero, b is the GCD; (3) if c is
  ! no zero, b becomes a and c becomes c and go back to
  ! Step (1).  This process will continue until c is zero.
  ! ---------------------------------------------------------

     IMPLICIT  NONE
  
     INTEGER, intent(in)   :: a0, b0
     INTEGER               :: a, b, c
     a = a0
     b = b0
  
     IF (a < b) THEN       ! since a >= b must be true, they
        c = a              ! are swapped if a < b
        a = b
        b = c
     END IF
  
     DO                    ! now we have a <= b
        c = MOD(a, b)      ! compute c, the reminder
        IF (c == 0) EXIT   ! if c is zero, we are done.  GCD = b
        a = b              ! otherwise, b becomes a
        b = c              ! and c becomes b
     END DO                ! go back
  
  END function  GCD

End Module mma

