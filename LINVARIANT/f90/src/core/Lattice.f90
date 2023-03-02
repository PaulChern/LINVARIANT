Module Lattice

  Use Parameters
  Use Constants

  Implicit none
  Contains

  Include "GetHeterostructureStrain.f90"

  Subroutine Pbc(x,y,z,x0,y0,z0)
    Implicit none
    integer, Intent(in)     :: x, y, z
    integer, Intent(out)    :: x0, y0, z0

    x0 = (x+1)-floor(real(x)/real(cgrid%n1))*cgrid%n1
    y0 = (y+1)-floor(real(y)/real(cgrid%n2))*cgrid%n2
    z0 = (z+1)-floor(real(z)/real(cgrid%n3))*cgrid%n3

  End Subroutine Pbc

  Function GridPbc(a,p) Result(b)
    Implicit none
    Integer::a,p,b

    b=a-floor((real(a)-1.0)/real(p))*p
  End Function GridPbc

  Subroutine PbcDiff(diff, numdiff)
    Implicit none
    Integer,intent(in)     :: numdiff
    Real(dp),intent(inout) :: diff(numdiff)
    Integer                :: i

    Do i = 1, numdiff
      if(diff(i).gt.-0.5.and.diff(i).le.0.5) diff(i) = diff(i)
      if(diff(i).gt.0.5) diff(i) = diff(i) - 1.0
      if(diff(i).le.-0.5) diff(i) = diff(i) + 1.0
    End do

  End subroutine PbcDiff

  Function eta2eij(eta) Result(eij)
    Implicit none
    Real*8,  Intent(in)    :: eta(6)
    Real*8                 :: eij(3,3)

    eij(1,1) = eta(1)
    eij(2,2) = eta(2)
    eij(3,3) = eta(3)
    eij(2,3) = eta(4)
    eij(1,3) = eta(5)
    eij(1,2) = eta(6)

    eij(3,2) = eta(4)
    eij(3,1) = eta(5)
    eij(2,1) = eta(6)

  End Function eta2eij

  Function eij2eta(eij) Result(eta)
    Implicit none
    Real*8,  Intent(in)    :: eij(3,3)
    Real*8                 :: eta(6)

    eta(1) = eij(1,1)
    eta(2) = eij(2,2)
    eta(3) = eij(3,3)
    eta(4) = eij(2,3)
    eta(5) = eij(1,3)
    eta(6) = eij(1,2)

  End Function eij2eta

  Function GetStrainVolume(e0ij) Result(Volume)

    Implicit none
    Real*8,  Intent(in)    :: e0ij(3,3)
    Real*8                 :: Volume

    Volume = 1.0D0 + e0ij(1,1) + e0ij(2,2) + e0ij(3,3) &
           - e0ij(1,2)**2 - e0ij(1,3)**2 - e0ij(2,3)**2 &
           + e0ij(1,1)*e0ij(2,2) + e0ij(1,1)*e0ij(3,3) + e0ij(2,2)*e0ij(3,3) &
           - e0ij(1,1)*e0ij(2,3)**2 - e0ij(2,2)*e0ij(1,3)**2 - e0ij(3,3)*e0ij(1,2)**2 &
           + e0ij(1,1)*e0ij(2,2)*e0ij(3,3) + 2*e0ij(1,2)*e0ij(1,3)*e0ij(2,3)

  End Function GetStrainVolume

  Function Cell2Volume(A) Result(VOLCELL)

    ! CALCULATES THE VOLUME OF THE UNIT CELL,POSSIBLY WITH A MINUS SIGN

    Implicit none
    Real*8::A(3,3)
    Real*8::VOLCELL

    VOLCELL=(A(2,1)*A(3,2)-A(3,1)*A(2,2))*A(1,3)+&
    (A(3,1)*A(1,2)-A(1,1)*A(3,2))*A(2,3)+&
    (A(1,1)*A(2,2)-A(2,1)*A(1,2))*A(3,3)

  End Function Cell2Volume

  Subroutine GetBZKList
    Implicit none
    Integer    :: ikpt, ix, iy, iz
    Real(dp)   :: kpt(3)
   
    Allocate(kbz(4,kgrid%npts)) 
    ikpt = 0
    Do ix = 1, kgrid%n1
      Do iy = 1, kgrid%n2
        Do iz = 1, kgrid%n3
          ikpt = ikpt + 1
          kpt = [0.0_dp + 1.0_dp*(ix-1)/kgrid%n1, &
                 0.0_dp + 1.0_dp*(iy-1)/kgrid%n2, &
                 0.0_dp + 1.0_dp*(iz-1)/kgrid%n3]
          Call PbcDiff(kpt, 3)
          kbz(1,ikpt) = kpt(1)
          kbz(2,ikpt) = kpt(2)
          kbz(3,ikpt) = kpt(3)
          kbz(4,ikpt) = 1.0_dp/kgrid%npts
        End do
      End do
    End do

  End subroutine GetBZKList

  Subroutine RECLAT(A,B,IOPT)

    ! CALCULATES RECIPROCAL LATTICE VECTORS.THEIR PRODUCT WITH DIRECT
    ! LATTICE VECTORS IS 1 IF IOPT=0 OR 2*PI IF IOPT=1

    IMPLICIT NONE
    Real*8::A(3,3),B(3,3)
    integer::iopt,i
    Real*8::pi,c,ci

    PI=ACOS(-1.D0)
    B(1,1)=A(2,2)*A(3,3)-A(3,2)*A(2,3)
    B(2,1)=A(3,2)*A(1,3)-A(1,2)*A(3,3)
    B(3,1)=A(1,2)*A(2,3)-A(2,2)*A(1,3)
    B(1,2)=A(2,3)*A(3,1)-A(3,3)*A(2,1)
    B(2,2)=A(3,3)*A(1,1)-A(1,3)*A(3,1)
    B(3,2)=A(1,3)*A(2,1)-A(2,3)*A(1,1)
    B(1,3)=A(2,1)*A(3,2)-A(3,1)*A(2,2)
    B(2,3)=A(3,1)*A(1,2)-A(1,1)*A(3,2)
    B(3,3)=A(1,1)*A(2,2)-A(2,1)*A(1,2)
    C=1.D0

    if(IOPT.eq.1) then
      C=2.D0*PI
    endif

    do i=1,3
      CI=C/(A(1,i)*B(1,i)+A(2,i)*B(2,i)+A(3,i)*B(3,i))
      B(1,i)=B(1,i)*CI
      B(2,i)=B(2,i)*CI
      B(3,i)=B(3,i)*CI
    end do

  End Subroutine

End Module Lattice
