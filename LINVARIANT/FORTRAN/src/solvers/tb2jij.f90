  Subroutine GreenJij
    Use omp_lib
    Use GreenFunction
    Use Wannier

    Implicit None
    Integer        :: i, NumEnergy
    Integer        :: omp_done, omp_count
    Integer        :: FileHandle = 1111
    Complex(dp)    :: a, b, m
    Complex(dp)    :: V(NumWann/OrbMul,NumWann/OrbMul,2*SpinDim-3)
    Complex(dp)    :: Aij(3*SpinDim-5,3*SpinDim-5,NumWannSites,NumWannSites,2*SpinDim-3)
    Complex(dp)    :: fa(3*SpinDim-5,3*SpinDim-5,NumWannSites,NumWannSites,2*SpinDim-3)
    Complex(dp)    :: fb(3*SpinDim-5,3*SpinDim-5,NumWannSites,NumWannSites,2*SpinDim-3)
    Complex(dp)    :: fm(3*SpinDim-5,3*SpinDim-5,NumWannSites,NumWannSites,2*SpinDim-3)

    Aij = cmplx(0.0_dp,0.0_dp)
    NumEnergy = Size(ContourPath)
    Call ResolveWannSiteInd(OrbMul)
    Call ResolveWannBlockInd(1)
    Call GetWannDist

    Call TB2Delta(V)

    omp_count = 0
    omp_done = 0
    !$OMP    PARALLEL  DEFAULT(SHARED) PRIVATE(i,a,b,m,fa,fb,fm)
    !Reduction(+:Aij)
    !$OMP    DO 
    Do i = 1, NumEnergy-1
      a = ContourPath(i)
      b = ContourPath(i+1)
      m = (a+b)/2
      fa = VGVG(a,Jij_R,V)
      fb = VGVG(b,Jij_R,V)
      fm = VGVG(m,Jij_R,V)
      Aij = Aij + (b-a)*(fa+4*fm+fb)/6.0
      !$omp atomic update
      omp_count = omp_count + 1
      !$omp end atomic
      if (omp_count.gt.omp_done) then
        omp_done = omp_count
        write(*,*) "Integration: "//trim(int2str(omp_count))//"/"//trim(int2str(NumEnergy-1))
      end if
    End do
    !$OMP    END DO
    !$OMP    END PARALLEL

    Call Aij2Jij(Aij)

  End Subroutine GreenJij

