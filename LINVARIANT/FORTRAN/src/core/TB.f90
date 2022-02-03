Module TB
  use Constants
  use Inputs
  use Aux
  use Wannier
  
  Contains

  Function TBHk(k,ispin) Result(Hk)

    implicit none
    Real(dp),Intent(in) :: k(3)
    Integer,Intent(in)  :: ispin
    Real(dp)            :: dRij(3), kdotR
    Integer             :: m, n, i, WannInd(NumWann)
    Complex(dp)         :: Hk(NumWann, NumWann), factor
    Complex(dp)         :: HRmnExp(NumWann, NumWann)

    Hk = cmplx(0.0_dp, 0.0_dp)
    Do i = 1, NumRpts(ispin)
      Do m = 1, NumWann
        Do n = 1, NumWann
!         convention 1
          dRij = WannDist(:,m,n) + Ti0(1:3,i,ispin)
!         convention 2
!          dRij = Ti0(1:3,i,ispin)
          kdotR = k(1)*dRij(1) + k(2)*dRij(2) + k(3)*dRij(3)
          factor = Exp(2.0*pi*ii*kdotR)
          HRmnExp(m,n)=HRmn(m,n,i,ispin)*factor/Ti0(4,i,ispin)
!          HRmnExp(m,n)=HRmn(m,n,i,ispin)*factor
        End do 
      End do
      
      Hk = Hk + 0.5*(HRmnExp + Transpose(Conjg(HRmnExp)))
!      Hk = 0.5*(HRmnExp + Transpose(Conjg(HRmnExp)))
!      write(*,*) sum(Hk)
!      write(*,*) k
!      write(*,*) Ti0(1:3,i,ispin)
!      Pause

    End do

  end function TBHk

  Function TBH00(ispin) Result(H00)

    implicit none
    Integer,Intent(in)  :: ispin
    Integer             :: m, n, i
    Complex(dp)         :: H00(NumWann, NumWann)

    H00 = cmplx(0.0_dp, 0.0_dp)
    Do i = 1, NumRpts(ispin)
      If(Sum(Ti0(1:3,i,ispin)**2).eq.0) then
        Do m = 1, NumWann
          Do n = 1, NumWann
            H00(m,n)=HRmn(m,n,i,ispin)
          End do
        End do
      End if
    End do

  End function TBH00

End module TB
