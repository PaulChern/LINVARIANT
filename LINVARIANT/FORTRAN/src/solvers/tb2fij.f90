  Subroutine GreenFij
    Use omp_lib
    Use GreenFunction

    Implicit None
    Integer        :: i, j, ix, iy, NumEnergy
    Integer        :: omp_done, omp_count
    Integer        :: FileHandle = 1111
    Complex(dp)    :: a, b, m
    Complex(dp)    :: W(NumWann,NumWann,6*SpinDim)
    Complex(dp)    :: W_bare(NumWann,NumWann,6*SpinDim)
    Complex(dp)    :: W_tmp(NumWann,NumWann,6*SpinDim)
    Complex(dp)    :: Fij(NumWann,NumWann,3,3)
    Complex(dp)    :: Fij_local(NumWann,NumWann,3,3)
    Complex(dp)    :: fa(NumWann,NumWann,3,3)
    Complex(dp)    :: fb(NumWann,NumWann,3,3)
    Complex(dp)    :: fm(NumWann,NumWann,3,3)

    Fij = cmplx(0.0_dp,0.0_dp)
    Fij_local = cmplx(0.0_dp,0.0_dp)
    NumEnergy = Size(ContourPath)

    Call ResolveWannSiteInd(OrbMul)
    Call ResolveWannBlockInd(OrbMul)
    Call GetWannDist
    Call GDOS

!    Call DFTpot2W('Ti','tot',W)
    Call DFTpot2W('Ti','eff',W)
!    Call TB2W(W)
!    Call LoadW('W_eff_w200.dat',W)
!    Call LoadW('W_bare_w200.dat',W_bare)
!    Call DFTpot2W('Ti','kinetic',W_tmp)
!    Call DFTpot2W('Ti','bare',W)
!    Call DFTpot2W('Ti','bare',W)
    Call DFTpot2W('Ti','bare',W_bare)
!    Call DFTpot2W('Ti','hartree',W_tmp)
!    Call DFTpot2W('Ti','xc',W_tmp)
!    Call DFTpot2W('Ti','hb',W_tmp)
!    Call DFTpot2W('Ti','hxc',W_tmp)
!    W_bare = W - W_tmp

    omp_count = 0
    omp_done = 0
    !$OMP    PARALLEL  DEFAULT(SHARED) PRIVATE(i,a,b,m,fa,fb,fm,Fij_local)
    !$OMP    DO     Reduction(+:Fij)
    Do i = 1, NumEnergy-1
      a = ContourPath(NumEnergy-i)
      b = ContourPath(NumEnergy-i+1)
      m = (a+b)/2
      fa = WGWG(a,Jij_R,W,W_bare)
      fb = WGWG(b,Jij_R,W,W_bare)
      fm = WGWG(m,Jij_R,W,W_bare)
      Fij_local = (b-a)*(fa+4*fm+fb)/6.0
      Fij = Fij + Fij_local
      !$omp atomic update
      omp_count = omp_count + 1
      !$omp end atomic
      if (omp_count.gt.omp_done) then
        omp_done = omp_count
        write(*,'(A24,6F10.4)') "Integration: "//trim(int2str(omp_count))//"/"//trim(int2str(NumEnergy-1)), a, b, b-a
        write(*, "(3F16.10)") ((Aimag(Sum(Fij_local(:,:,ix,iy)))/ABS(b-a), iy=1,3), ix=1,3)
        Call WriteFijOrbital(a,b,Jij_R,Aimag(Fij_local)/Abs(b-a))
      end if
    End do
    !$OMP    END DO
    !$OMP    END PARALLEL

    Open(FileHandle,file=trim(Solver)//'.out/'//'Fij.dat',form='formatted',status='unknown')
    write(FileHandle, "(2I4)") JijSites(1), JijSites(2)
    write(FileHandle, "(3F16.10)") ((Aimag(Sum(Fij(:,:,ix,iy))), iy=1,3), ix=1,3)
    write(FileHandle, *)
    Close(FileHandle)

  End Subroutine GreenFij

