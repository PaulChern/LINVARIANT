Module GreenFunction
  Use aux
  Use mma
  Use Constants
  Use TB
  Use Wannier

  Contains

  Function H2G(Hmat, ene, fermi) result(G)
    Implicit None
    Complex(dp), Intent(in)    :: Hmat(NumWann,NumWann)
    Complex(dp), Intent(in)    :: ene
    Real(dp), Intent(in)       :: fermi
    Complex(dp)                :: G(NumWann,NumWann)
    Integer                    :: i

    G = -1.0D0*Hmat
    Do i = 1, NumWann
       G(i,i)= ene + fermi + G(i,i)
    End Do
    Call Inverse(G, NumWann)

  End Function H2G

  Function Eigen2G(Hmat, ene, fermi) result(G)
    Implicit None
    Complex(dp)                :: Hmat(NumWann,NumWann)
    Complex(dp),Intent(in)     :: ene
    Real(dp), Intent(in)       :: fermi
    Real(dp)                   :: EigenValues(NumWann)
    Complex(dp)                :: MatDiagE(NumWann,NumWann)
    Complex(dp)                :: G(NumWann,NumWann)
    Integer                    :: i

    Call EigenSystem('V','U',Hmat,NumWann,EigenValues)
    
    MatDiagE = cmplx(0.0_dp,0.0_dp)
    Do i = 1, NumWann
      MatDiagE(i,i) = 1.0_dp/(ene + fermi - EigenValues(i))
    End do

    G = MatMul(Hmat,MatMul(MatDiagE,Conjg(Transpose(Hmat))))
!    G = MatDiagE
  
  End Function Eigen2G

  Function TB2GR(ene, fermi, R, ispin) result(GR)
    Implicit None
    
    Complex(dp), Intent(in)   :: ene
    Real(dp), Intent(in)      :: fermi
    Integer,  Intent(in)      :: R(3), ispin
    Complex(dp)               :: Gk(NumWann,NumWann), Hk(NumWann,NumWann)
    Complex(dp)               :: GR(NumWann,NumWann)
    Real(dp)                  :: kdotR
    Integer                   :: i,j

    GR = cmplx(0.0_dp, 0.0_dp)
    Do i = 1, kgrid%npts
      Hk = TBHk(kbz(1:3,i), ispin)
      Gk = H2G(Hk, ene, fermi)
!      Gk = Eigen2G(Hk, ene, fermi)
!      write(*,*) i, sum(H2G(Hk, ene, fermi)-Eigen2G(Hk, ene, fermi))
      kdotR = kbz(1,i)*R(1)+kbz(2,i)*R(2)+kbz(3,i)*R(3)
      GR = GR + kbz(4,i)*Gk*Exp(-2*pi*ii*kdotR)
    End do

  End Function TB2GR

  Subroutine TB2Delta(Delta)

    Implicit None
    Complex(dp),intent(inout) :: Delta(NumWann/OrbMul,NumWann/OrbMul,2*SpinDim-3)
    Complex(dp)               :: H00(NumWann,NumWann)
    Integer                   :: i, i1, i2, is
    
    Delta = cmplx(0.0_dp, 0.0_dp)
    If(SpinDim.eq.2) then
      Do i = 1, kgrid%npts
        Delta(:,:,1) = Delta(:,:,1) + (TBHk(kbz(1:3,i),1) - TBHk(kbz(1:3,i),2))*kbz(4,i)
      End do
    else if(SpinDim.eq.3) then
      Do is = 1, 3
        H00 = TBH00(is)
        Do i = 1, NumWannSites
          i1=WannBlockInd(1,i)
          i2=WannBlockInd(2,i)
          Delta(i1:i2,i1:i2,is) = PauliNormDecompose(H00(2*i1-1:2*i2,2*i1-1:2*i2), (i2-i1+1)*2)
        End do
      End do
    End if

  End Subroutine TB2Delta

  Subroutine TB2W(W)

    Implicit None
    Complex(dp),Intent(inout) :: W(NumWann,NumWann,6*SpinDim)
    Complex(dp)               :: Hmn(NumWann,NumWann,7*SpinDim)
    Complex(dp)               :: H00(NumWann,NumWann)
    Complex(dp)               :: H00_disp(NumWann,NumWann)
    integer                   :: iow=200
    Integer                   :: idisp,ik,is, iw, jw, m, n, i, j, k

    Hmn = cmplx(0.0_dp, 0.0_dp)
    W = cmplx(0.0_dp, 0.0_dp)
    H00 = TBH00(1)
    Hmn(:,:,1) = H00

!    Do ik = 1, kgrid%npts
      Do is = 1, SpinDim
        Do idisp = 1, 6
!          W(:,:,(is-1)*6+ihr) = W(:,:,(is-1)*6+ihr) &
!                + (TBHk(kbz(1:3,ik),(is-1)*7+ihr+1) - TBHk(kbz(1:3,ik),(is-1)*7+1))*kbz(4,ik)
          H00_disp = TBH00((is-1)*7+idisp+1)
          Hmn(:,:,idisp+1) = H00_disp
          W(:,:,(is-1)*6+idisp) = H00_disp - H00
!          Call CMatRoundup(W(:,:,idisp),NumWann,3)
        End do
      End do
!    End do

!    Do i = 1, 2
!      Call CMatSymmetrize(W(:,:,i),W(:,:,i+3),NumWann,2.0E-2_dp)
!    End do

    Call WriteWHmn('tb',W,Hmn)

    write(*,*) "screened potential calculated."

  End Subroutine TB2W

  Subroutine DFTpot2W(site,v,W)
    Implicit none
    complex(dp),intent(out)   :: W(NumWann,NumWann,6*SpinDim)
    complex(dp)               :: T(NumWann,NumWann,6*SpinDim)
    complex(dp)               :: Hmn(NumWann,NumWann,7*SpinDim)
    complex(dp)               :: Tmn(NumWann,NumWann,7*SpinDim) 
    Real(dp)                  :: dftpot(3*rgrid%n1,3*rgrid%n2,3*rgrid%n3,7,3,SpinDim)
    character(*)              :: site, v
    integer                   :: iow=200, ipot
    integer                   :: i, j, k, is, iw, jw, ix, iy, iz
    Real(dp)                  :: Vpot(3*rgrid%n1,3*rgrid%n2,3*rgrid%n3,7*SpinDim)

    W = cmplx(0.0_dp,0.0_dp)
    T = cmplx(0.0_dp,0.0_dp)
    Hmn = cmplx(0.0_dp,0.0_dp)
    Tmn = cmplx(0.0_dp,0.0_dp)

    call LoadDFTpot(dftpot,site)

    SELECT CASE (trim(v))
      CASE ('hxc')
        Do i = 1, 7
          Vpot(:,:,:,i) = dftpot(:,:,:,i,1,1) - dftpot(:,:,:,i,2,1)
        End do
      CASE ('hartree')
        Do i = 1, 7
          Vpot(:,:,:,i) = dftpot(:,:,:,i,3,1) - dftpot(:,:,:,i,2,1)
        End do
      CASE ('xc') 
        Do i = 1, 7
          Vpot(:,:,:,i) = dftpot(:,:,:,i,1,1) - dftpot(:,:,:,i,3,1)
        End do
      CASE ('eff')
        Vpot(:,:,:,:) = dftpot(:,:,:,:,1,1)
      CASE ('tot')
        Vpot(:,:,:,:) = dftpot(:,:,:,:,1,1)
      CASE ('bare')
        Vpot(:,:,:,:) = dftpot(:,:,:,:,2,1)
      CASE ('bh')
        Vpot(:,:,:,:) = dftpot(:,:,:,:,3,1)
      CASE DEFAULT
        Vpot(:,:,:,:) = dftpot(:,:,:,:,2,1)
    END SELECT

    if(trim(v).ne.'kinetic') then
      Do is = 1, SpinDim
        Do i = 1, 7
          Do jw = 1, NumWann
            Do iw = 1, NumWann
              Hmn(iw,jw,i) = Sum(WannFunc(:,:,:,iw,1,is)*Vpot(:,:,:,i)*WannFunc(:,:,:,jw,1,is))
            End do
          End do
        End do
      End do
      Do i = 1, 6
        W(:,:,i) = Hmn(:,:,i+1) - Hmn(:,:,1)
        Call CMatRoundup(W(:,:,i),NumWann,3)
      End do
    end if


    if(trim(v).eq.'tot'.or.trim(v).eq.'kinetic') then
      Do is = 1, SpinDim
        Do i = 1, 7
          Do jw = 1, NumWann
            Do iw = 1, NumWann
              Tmn(iw,jw,i) = hbar2_2me*Sum(NablaW(WannFunc(:,:,:,iw,i,is))*Conjg(NablaW(WannFunc(:,:,:,jw,i,is))))
            End do
          End do
          write(*,'(A25,I)') 'kinetic energy for disp: ', i
        End do
      End do
      Do i = 1, 6
        T(:,:,i) = Tmn(:,:,i+1) - Tmn(:,:,1)
        Call CMatRoundup(T(:,:,i),NumWann,3)
      End do
    end if

    W = T + W

    Call WriteWHmn(v,W,Hmn+Tmn)

    write(*,*) trim(v)//" potential calculated."

  End subroutine DFTpot2W

  Subroutine LoadW(filename,W)
    Implicit none
    character(*),intent(in)   :: filename
    complex(dp),intent(inout) :: W(NumWann,NumWann,6*SpinDim)
    real(dp)                  :: W_real(NumWann,NumWann,6*SpinDim)
    logical                   :: file_exists
    integer                   :: iow=200
    integer                   :: i, j, k, iw, jw

    W = cmplx(0.0_dp,0.0_dp)
    INQUIRE(FILE='dft/'//trim(filename), EXIST=file_exists)
    if (file_exists) then
      call open_input_file(iow,'dft/'//trim(filename))
      Do i = 1, 6*SpinDim
        Read(iow,*) ((W_real(9+iw,9+jw,i), jw=1,9), iw=1,9)
        Read(iow,*)
      End do
      W = W_real*cmplx(1.0_dp,0.0_dp)
      write(*,*) "potentail file ", filename, "loaded."
    else
      write(*,*) "potentail file not exist!"
      call abort
    end if

  End subroutine LoadW

  Subroutine LoadDFTpot(dftpot,site)
    Implicit none
    Real(dp),intent(inout)    :: dftpot(3*rgrid%n1,3*rgrid%n2,3*rgrid%n3,7,3,SpinDim)
    character(*)              :: site
    character(10)             :: pot(3),disp(7)
    character(72)             :: filename
    integer                   :: i, is, idisp

    pot=['tot', 'bare', 'hb']
    disp=['0', '_disp+x', '_disp+y', '_disp+z', '_disp-x', '_disp-y', '_disp-z']

    do is = 1, SpinDim
      do i = 1, 3
        do idisp = 1, 7
          filename='dft/'//trim(site)//'_'//trim(pot(i))//trim(disp(idisp))//'.xsf'
          dftpot(:,:,:,idisp,i,is) = ry_ev*ReadXSF(filename)
        end do
      end do
    end do
    write(*,*) "dft potential on "//trim(site)//" loaded."

  End Subroutine

  Function VGVG(ene,R,V) result(Aij)

    Implicit None
    Complex(dp),intent(in) :: ene
    Integer,intent(in)     :: R(3)
    Complex(dp),intent(in) :: V(NumWann/OrbMul,NumWann/OrbMul,2*SpinDim-3)
    Integer                :: is,i,j,i1,i2,j1,j2,s1,s2
    Integer                :: BlockInd(2,NumWannSites)
    Complex(dp)            :: GRij(NumWann/OrbMul,NumWann/OrbMul,3*SpinDim-5,2*SpinDim-3)
    Complex(dp)            :: GRji(NumWann/OrbMul,NumWann/OrbMul,3*SpinDim-5,2*SpinDim-3)
    Complex(dp)            :: Aij(3*SpinDim-5,3*SpinDim-5,NumWannSites,NumWannSites,2*SpinDim-3)

    If(SpinDim.eq.3) then
      Do is = 1, 3
        GRij(:,:,:,is)  = PauliDecompose(TB2GR(ene,efermi(is),R,is),NumWann,NumWann)
        GRji(:,:,:,is) = PauliDecompose(TB2GR(ene,efermi(is),-1*R,is),NumWann,NumWann)
      End do
    else if(SpinDim.eq.2) then
      GRij(:,:,1,1)  = TB2GR(ene,efermi(1),R,1)
      GRji(:,:,1,1) = TB2GR(ene,efermi(2),-1*R,2)
    End if

    Do is = 1, 2*SpinDim-3
      Do i = 1, NumWannSites
        Do j = 1, NumWannSites
          i1=WannBlockInd(1,i)
          i2=WannBlockInd(2,i)
          j1=WannBlockInd(1,j)
          j2=WannBlockInd(2,j)
          Do s1=1,3*SpinDim-5
            Do s2=1,3*SpinDim-5
              Aij(s1,s2,i,j,is) = &
              Sum( &
              MatMul(V(i1:i2,i1:i2,is),GRij(i1:i2,j1:j2,s1,is))* &
              Transpose(MatMul(V(j1:j2,j1:j2,is),GRji(j1:j2,i1:i2,s2,is))))/pi
            End do
          End do
        End do
      End do
    End do

  End function VGVG

!  Subroutine GreenJij
!
!    Implicit None
!    Integer        :: i, NumEnergy
!    Integer        :: omp_done, omp_count
!    Integer        :: FileHandle = 1111
!    Complex(dp)    :: a, b, m
!    Complex(dp)    :: V(NumWann/OrbMul,NumWann/OrbMul,2*SpinDim-3)
!    Complex(dp)    :: Aij(3*SpinDim-5,3*SpinDim-5,NumWannSites,NumWannSites,2*SpinDim-3)
!    Complex(dp)    :: fa(3*SpinDim-5,3*SpinDim-5,NumWannSites,NumWannSites,2*SpinDim-3)
!    Complex(dp)    :: fb(3*SpinDim-5,3*SpinDim-5,NumWannSites,NumWannSites,2*SpinDim-3)
!    Complex(dp)    :: fm(3*SpinDim-5,3*SpinDim-5,NumWannSites,NumWannSites,2*SpinDim-3)
!
!    Aij = cmplx(0.0_dp,0.0_dp)
!    NumEnergy = Size(ContourPath)
!    Call ResolveWannSiteInd(OrbMul)
!    Call ResolveWannBlockInd(1)
!    Call GetWannDist
!   
!    Call TB2Delta(V) 
! 
!    omp_count = 0
!    omp_done = 0
!    !$OMP    PARALLEL  DEFAULT(SHARED) PRIVATE(i,a,b,m,fa,fb,fm) Reduction(+:Aij)
!    !$OMP    DO 
!    Do i = 1, NumEnergy-1
!      a = ContourPath(i)
!      b = ContourPath(i+1)
!      m = (a+b)/2
!      fa = VGVG(a,Jij_R,V)
!      fb = VGVG(b,Jij_R,V)
!      fm = VGVG(m,Jij_R,V)
!      Aij = Aij + (b-a)*(fa+4*fm+fb)/6.0
!      !$omp atomic update
!      omp_count = omp_count + 1
!      !$omp end atomic
!      if (omp_count.gt.omp_done) then
!        omp_done = omp_count
!        write(*,*) "Integration: "//trim(int2str(omp_count))//"/"//trim(int2str(NumEnergy-1))
!      end if
!    End do
!    !$OMP    END DO
!    !$OMP    END PARALLEL
!
!    Call Aij2Jij(Aij)
!
!  End Subroutine GreenJij

  Subroutine Aij2Jij(Aij)
    Implicit None
    Complex(dp),intent(in) :: Aij(3*SpinDim-5,3*SpinDim-5,NumWannSites,NumWannSites,2*SpinDim-3)
    Real(dp)               :: J0(3,3,NumWannSites,NumWannSites)
    Real(dp)               :: Jani(3,3,NumWannSites,NumWannSites)
    Real(dp)               :: Dvec(3,NumWannSites,NumWannSItes)
    Real(dp)               :: MatTmp(3,3,3), VecTmp(3,3), Rot(3,3)
    Real(dp)               :: IdenMat(3,3), Rot120(3,3), Rot240(3,3)
    Integer                :: is, i, j, s1, s2, m, n
   
    IdenMat = 0.0_dp
    Do i = 1, 3
      IdenMat(i,i) = 1.0_dp
    End do

    Rot120 = 0.0_dp
    Rot120(1,3) = 1.0_dp
    Rot120(2,1) = 1.0_dp
    Rot120(3,2) = 1.0_dp
    Rot240 = MatMul(Rot120,Rot120)

    J0 = 0.0_dp
    Jani = 0.0_dp
    Dvec = 0.0_dp
    MatTmp = 0.0_dp
    VecTmp = 0.0_dp

    Do i = 1, NumWannSites
      Do j = 1, NumWannSites
        If(SpinDim.eq.3) then
          do is = 1, 3
            J0(:,:,i,j) = J0(:,:,i,j) + &
              1.0D3*Aimag(Aij(1,1,i,j,is) - Aij(2,2,i,j,is) - Aij(3,3,i,j,is) - Aij(4,4,i,j,is))*IdenMat/3
          end do
!          J0(:,:,i,j)=1.0D3*Aimag(Aij(1,1,i,j,1)-Aij(2,2,i,j,1)-Aij(3,3,i,j,1)-Aij(4,4,i,j,1))*IdenMat
          do is = 1, 3
            Do s1 = 2, 4
              Do s2 = 2, 4
                MatTmp(s1-1,s2-1,is) = 1.0D3*Aimag(Aij(s1,s2,i,j,is)+Aij(s2,s1,i,j,is))
              End do
            End do
            MatTmp(is,:,is) = 0.0_dp
            MatTmp(:,is,is) = 0.0_dp
          end do
          Jani(:,:,i,j) = (MatTmp(:,:,1) + MatTmp(:,:,2) + MatTmp(:,:,3))/2
!            + MatMul(MatMul(Transpose(Rot120),MatTmp(:,:,1)),Rot120)/2 &
!            + MatMul(MatMul(Transpose(Rot240),MatTmp(:,:,2)),Rot240)/2 &
!          Jani(:,:,i,j) = MatTmp(:,:,1)
          Do is = 1, 3
            Do s1 = 2, 4
              VecTmp(s1-1,is) = 1.0D3*Real(Aij(1,s1,i,j,is)-Aij(s1,1,i,j,is))
            End do
            VecTmp(is,is) = 0.0_dp
          End do
          Dvec(:,i,j) = (VecTmp(:,1) + VecTmp(:,2) + VecTmp(:,3))/2
!          Dvec(:,i,j) = (Rot3(Rot240,VecTmp(:,1)) + Rot3(Rot120,VecTmp(:,2)) + VecTmp(:,3))/2
!          Dvec(:,i,j) = VecTmp(:,1)
        else if(SpinDim.eq.2) then
          J0(:,:,i,j) = 1.0D3*Aimag(Aij(1,1,i,j,1))*IdenMat/4
        End if
      end do
    end do

    Call WriteJij(J0,Jani,Dvec)
 
  End Subroutine Aij2Jij

  Subroutine WriteJij(J0,Jani,Dvec)
    Implicit None
    Integer             :: FileHandle = 1999
    Integer             :: i, j, s1, s2, m, n
    Real(dp)            :: Jij(3,3)
    Real(dp),intent(in) :: J0(3,3,NumWannSites,NumWannSites)
    Real(dp),intent(in) :: Jani(3,3,NumWannSites,NumWannSites)
    Real(dp),intent(in) :: Dvec(3,NumWannSites,NumWannSItes)

    Open(FileHandle,file=trim(Solver)//'.out/'//'Jij.dat',form='formatted',status='unknown')
    Do i = 1, NumWannSites
      Do j = 1, NumWannSites
        Jij(:,:) = J0(:,:,i,j) + Jani(:,:,i,j) + D2Jij(Dvec(:,i,j))
        write(FileHandle, "(2I4)") SiteOrbInfo(1,i), SiteOrbInfo(1,j)
        write(FileHandle, "(A6)") "   Jij" 
        write(FileHandle, "(3F10.4)") (Jij(m,:), m=1,3)
        write(FileHandle, "(A6)") "   J0 " 
        write(FileHandle, "(3F10.4)") (J0(m,:,i,j), m=1,3)
        write(FileHandle, "(A6)") "   Ja " 
        write(FileHandle, "(3F10.4)") (Jani(m,:,i,j), m=1,3)
        write(FileHandle, "(A6)") "   D  " 
        write(FileHandle, "(3F10.4)") Dvec(:,i,j)
        write(FileHandle, *)
      End do
    End do
    Close(FileHandle) 
  End Subroutine WriteJij

  Function WGWG(ene,R,W,W_bare) result(Fij)

    Implicit None
    Complex(dp),intent(in) :: ene
    Integer,intent(in)     :: R(3)
    Complex(dp),intent(in) :: W(NumWann,NumWann,3)
    Complex(dp),intent(in) :: W_bare(NumWann,NumWann,3)
    Integer                :: ifij,ia,ib,i,j,i1,i2,j1,j2,s1,s2,m,n
    Integer                :: NumFij
    Complex(dp)            :: G12(NumWann,NumWann), G21(NumWann,NumWann)
    Complex(dp)            :: Fij(NumWann,NumWann,3,3), F22(3,3), F22_tmp
    Real(dp)               :: lattpotim(3)

    lattpotim = norm2(unitcell%latt_a(1,:))/potim/rgrid%n1
    lattpotim = norm2(unitcell%latt_a(2,:))/potim/rgrid%n2
    lattpotim = norm2(unitcell%latt_a(3,:))/potim/rgrid%n3

    NumFij = 1
    if(SpinDim.eq.2) NumFij = 2

    Fij = cmplx(0.0_dp,0.0_dp)
    F22 = cmplx(0.0_dp,0.0_dp)
    
    Do ifij = 1, NumFij
      G12 = TB2GR(ene,efermi(4*ifij-3),R,7*ifij-6)
      G21 = TB2GR(ene,efermi(4*ifij-3),-1*R,7*ifij-6)
      Call CMatSymmetrize(G12,G21,NumWann,1.0E-6_dp)
   
      i1=WannBlockInd(1,JijSites(1))
      i2=WannBlockInd(2,JijSites(1))
      j1=WannBlockInd(1,JijSites(2))
      j2=WannBlockInd(2,JijSites(2))
      Do ia = 1, 3
        Do ib = 1, 3
          Do s1 = 1, 2
            Do s2 = 1, 2
              Fij(i1:i2,j1:j2,ia,ib) = Fij(i1:i2,j1:j2,ia,ib) + &
              (-2*s1+3)*(-2*s2+3)/(2*pi)/lattpotim(ia)/lattpotim(ib)/4* &
              (MatMul(W_bare(i1:i2,i1:i2,(ifij-1)*6+(s1-1)*3+ia),G12(i1:i2,j1:j2))* &
               Transpose(MatMul(W(j1:j2,j1:j2,(ifij-1)*6+(s2-1)*3+ib),G21(j1:j2,i1:i2))) + &
               MatMul(W(i1:i2,i1:i2,(ifij-1)*6+(s1-1)*3+ia),G12(i1:i2,j1:j2))* &
               Transpose(MatMul(W_bare(j1:j2,j1:j2,(ifij-1)*6+(s2-1)*3+ib),G21(j1:j2,i1:i2))))
            End do
          End do
          Do n = 10, 18
            Do m = 10, 18
              Do i = 10, 18
                Do j = 10, 18
                  F22_tmp = &
                  ((W_bare(n,m,ia)-W_bare(n,m,ia+3))*G12(m,i)*(W(i,j,ib)-W(i,j,ib+3))*G21(j,n) + &
                   (W(n,m,ia)-W(n,m,ia+3))*G12(m,i)*(W_bare(i,j,ib)-W_bare(i,j,ib+3))*G21(j,n))/pi/lattpotim(ia)/lattpotim(ib)/4
                  F22(ia,ib) = F22(ia,ib) + F22_tmp
                  If(Abs(Aimag(F22_tmp)).gt.0.000001) then
                    Write(*,'(A6,I2,I2,I3,I3,I3,I3,2F12.6,F12.6)') "wgwg:", ia, ib, n-9, m-9, i-9, j-9, ene, Aimag(F22_tmp)
                  End if
                  If(Abs(G12(m,i)).gt.0.000001.and.n.eq.10.and.j.eq.10) then
                    Write(*,'(A6,I3,I3,2F12.6,2F12.6)') "Gij:", m-9, i-9, ene, G12(m,i)
                  end if
                  If(Abs(G21(j,n)).gt.0.000001.and.i.eq.10.and.m.eq.10) then
                    Write(*,'(A6,I3,I3,2F12.6,2F12.6)') "Gji:", j-9, n-9, ene, G21(j,n)
                  end if
                End do
              End do
            End do
          End do
          Do n = 10, 18
            Do j = 10, 18
              F22_tmp = cmplx(0.0_dp,0.0_dp)
              Do m = 10, 18
                Do i = 10, 18
                  F22_tmp = F22_tmp + &
                  ((W_bare(n,m,ia)-W_bare(n,m,ia+3))*G12(m,i)*(W(i,j,ib)-W(i,j,ib+3))*G21(j,n)+ &
                   (W(n,m,ia)-W(n,m,ia+3))*G12(m,i)*(W_bare(i,j,ib)-W_bare(i,j,ib+3))*G21(j,n))/pi/lattpotim(ia)/lattpotim(ib)/4
                End do
              End do
              Write(*,'(A6,I2,I2,I3,I3,2F12.6,F12.6)') "mg:", ia, ib, n-9, j-9, ene, Aimag(F22_tmp)
            End do
          End do
        End do
      End do
      if(WriteG) Call WriteGRij(ene,R,G12,G21)
    End do

  End function WGWG

  Function WGWGWG(ene,R,W1,W2,W3) result(Fij)

    Implicit None
    Complex(dp),intent(in) :: ene
    Integer,intent(in)     :: R(3)
    Complex(dp),intent(in) :: W1(NumWann,NumWann,3)
    Complex(dp),intent(in) :: W2(NumWann,NumWann,3)
    Complex(dp),intent(in) :: W3(NumWann,NumWann,3)
    Integer                :: ifij,ia,ib,ic,i,j,i1,i2,j1,j2,s1,s2,s3,m,n
    Integer                :: NumFij, R1(3), R2(3)
    Complex(dp)            :: G12(NumWann,NumWann), G21(NumWann,NumWann)
    Complex(dp)            :: G22(NumWann,NumWann)
    Complex(dp)            :: G23(NumWann,NumWann), G32(NumWann,NumWann)
    Complex(dp)            :: G13(NumWann,NumWann), G31(NumWann,NumWann)
    Complex(dp)            :: Fij(NumWann,NumWann,3,3)
    Real(dp)               :: lattpotim(3)

    lattpotim = norm2(unitcell%latt_a(1,:))/potim/rgrid%n1
    lattpotim = norm2(unitcell%latt_a(2,:))/potim/rgrid%n2
    lattpotim = norm2(unitcell%latt_a(3,:))/potim/rgrid%n3

    R1 = 0
    R2 = 0
    R1(1) = R(1)
    R2(2) = R(2)

    NumFij = 1
    if(SpinDim.eq.2) NumFij = 2

    Fij = cmplx(0.0_dp,0.0_dp)

    Do ifij = 1, NumFij
      G12 = TB2GR(ene,efermi(4*ifij-3),R,7*ifij-6)
      G22 = TB2GR(ene,efermi(4*ifij-3),[0,0,0],7*ifij-6)
      G21 = TB2GR(ene,efermi(4*ifij-3),-1*R,7*ifij-6)

      i1=WannBlockInd(1,JijSites(1))
      i2=WannBlockInd(2,JijSites(1))
      j1=WannBlockInd(1,JijSites(2))
      j2=WannBlockInd(2,JijSites(2))
      Do ia = 1, 3
        Do ib = 1, 3
          Do ic = 1, 3
            Do s1 = 1, 2
              Do s2 = 1, 2
                Do s3 = 1, 2
                  Fij(i1:i2,j1:j2,ia,ic) = Fij(i1:i2,j1:j2,ia,ic) + &
                  (-2*s1+3)*(-2*s2+3)*(-2*s3+3)*pi/lattpotim(ia)/lattpotim(ic)/4* &
                  (MatMul(MatMul(W1(i1:i2,i1:i2,(ifij-1)*6+(s1-1)*3+ia),G12(i1:i2,j1:j2)), &
                          MatMul(W2(i1:i2,i1:i2,(ifij-1)*6+(s2-1)*3+ib),G22(i1:i2,j1:j2)))* &
                  Transpose(MatMul(W3(j1:j2,j1:j2,(ifij-1)*6+(s3-1)*3+ic),G21(j1:j2,i1:i2)))) 
                End do
              End do
            End do
          End do
        End do
      End do
    End do

  End function WGWGWG

!  Subroutine GreenFij
!
!    Implicit None
!    Integer        :: i, j, ix, iy, NumEnergy
!    Integer        :: omp_done, omp_count
!    Integer        :: FileHandle = 1111
!    Complex(dp)    :: a, b, m
!    Complex(dp)    :: W(NumWann,NumWann,6*SpinDim)
!    Complex(dp)    :: W_bare(NumWann,NumWann,6*SpinDim)
!    Complex(dp)    :: W_tmp(NumWann,NumWann,6*SpinDim)
!    Complex(dp)    :: Fij(NumWann,NumWann,3,3)
!    Complex(dp)    :: Fij_local(NumWann,NumWann,3,3)
!    Complex(dp)    :: fa(NumWann,NumWann,3,3)
!    Complex(dp)    :: fb(NumWann,NumWann,3,3)
!    Complex(dp)    :: fm(NumWann,NumWann,3,3)
!
!    Fij = cmplx(0.0_dp,0.0_dp)
!    Fij_local = cmplx(0.0_dp,0.0_dp)
!    NumEnergy = Size(ContourPath)
!
!    Call ResolveWannSiteInd(OrbMul)
!    Call ResolveWannBlockInd(OrbMul)
!    Call GetWannDist
!    Call GDOS
!
!!    Call DFTpot2W('Ti','tot',W)
!    Call DFTpot2W('Ti','eff',W)
!!    Call TB2W(W)
!!    Call TB2W(W_bare)
!!    Call LoadW('W_eff_w200.dat',W)
!!    Call LoadW('W_bare_w200.dat',W_bare)
!!    Call DFTpot2W('Ti','kinetic',W_tmp)
!!    Call DFTpot2W('Ti','bare',W)
!!    Call DFTpot2W('Ti','bare',W)
!    Call DFTpot2W('Ti','bare',W_bare)
!!    Call DFTpot2W('Ti','hartree',W_tmp)
!!    Call DFTpot2W('Ti','xc',W_tmp)
!!    Call DFTpot2W('Ti','hb',W_tmp)
!!    Call DFTpot2W('Ti','hxc',W_tmp)
!!    W_bare = W - W_tmp
!
!    omp_count = 0
!    omp_done = 0
!    !$OMP    PARALLEL  DEFAULT(SHARED) PRIVATE(i,a,b,m,fa,fb,fm,Fij_local)
!    !$OMP    DO     Reduction(+:Fij)
!    Do i = 1, NumEnergy-1
!      a = ContourPath(NumEnergy-i)
!      b = ContourPath(NumEnergy-i+1)
!      m = (a+b)/2
!      fa = WGWG(a,Jij_R,W,W_bare)
!      fb = WGWG(b,Jij_R,W,W_bare)
!      fm = WGWG(m,Jij_R,W,W_bare)
!      Fij_local = (b-a)*(fa+4*fm+fb)/6.0
!      Fij = Fij + Fij_local
!      !$omp atomic update
!      omp_count = omp_count + 1
!      !$omp end atomic
!      if (omp_count.gt.omp_done) then
!        omp_done = omp_count
!        write(*,'(A24,6F10.4)') "Integration: "//trim(int2str(omp_count))//"/"//trim(int2str(NumEnergy-1)), a, b, b-a
!        write(*, "(3F16.10)") ((Aimag(Sum(Fij_local(:,:,ix,iy))), iy=1,3), ix=1,3)
!        Call WriteFijOrbital(a,b,Jij_R,Fij_local)
!      end if
!    End do
!    !$OMP    END DO
!    !$OMP    END PARALLEL
!
!    Open(FileHandle,file=trim(Solver)//'.out/'//'Fij.dat',form='formatted',status='unknown')
!    write(FileHandle, "(2I4)") JijSites(1), JijSites(2)
!    write(FileHandle, "(3F16.10)") ((Aimag(Sum(Fij(:,:,ix,iy))), iy=1,3), ix=1,3)
!    write(FileHandle, *)
!    Close(FileHandle)
!
!  End Subroutine GreenFij

  Subroutine GDOS
    Use omp_lib
    Implicit none
    Integer        :: i, iw, NumEnergy, FileHandle=100
    real(dp)       :: pocc(NumWann), occ(NumWann)
    Complex(dp)    :: a, b, m
    Complex(dp)    :: G00a(NumWann,NumWann), G00b(NumWann,NumWann)
    Complex(dp)    :: G00(NumWann,NumWann), G00m(NumWann,NumWann)

    NumEnergy = Size(ContourPath)
    occ = 0.0_dp
    call io_file_unit(FileHandle)
    Open(FileHandle,file=trim(Solver)//'.out/'//'dos.dat',form='formatted',status='unknown')

    !$OMP    PARALLEL  DEFAULT(SHARED)    PRIVATE(i,iw,a,b,m,G00a,G00b,G00m,G00,pocc)
    !$OMP    DO     Reduction(+:occ)
    Do i = 1, NumEnergy-1
      a = ContourPath(NumEnergy-i)
      b = ContourPath(NumEnergy-i+1)
      m = (a+b)/2
      G00a = TB2GR(a,efermi(1),[0,0,0],1)
      G00b = TB2GR(b,efermi(1),[0,0,0],1)
      G00m = TB2GR(m,efermi(1),[0,0,0],1)
      G00  = (G00a+4*G00m+G00b)/6.0
      do iw = 1, NumWann
        pocc(iw) = -Aimag((b-a)*G00(iw,iw))/pi/ABS(b-a)
      end do
      occ = occ + ABS(b-a)*pocc
      Write(FileHandle,'(2F10.4,31F10.4)') m, sum(pocc), pocc
    End do
    !$OMP    END DO
    !$OMP    END PARALLEL
    
    Write(FileHandle,'(A20,30F10.4)') " ", occ

    Write(*,*) "DOS calculated"
    Close(FileHandle)

  End Subroutine

  Subroutine WriteWHmn(v,W,Hmn)
    Implicit none
    complex(dp),intent(in)    :: W(NumWann,NumWann,6*SpinDim)
    complex(dp),intent(in)    :: Hmn(NumWann,NumWann,7*SpinDim)
    character(*),intent(in)   :: v
    character(10)             :: disp(7)
    integer                   :: iow=200
    integer                   :: i, j, k

    disp=['_disp0', '_disp+x', '_disp+y', '_disp+z', '_disp-x', '_disp-y', '_disp-z']

    Open(iow,file=trim(Solver)//'.out/'//'W_'//trim(v)//'.dat',form='formatted',status='unknown')
    Do i = 1, 6
      Write(iow,'(A15,I10)') "W Ti:", i
      Write(iow,'(A10,A4,3A10,A4,5A10)') "s", "|", "pz", "px", "py", "|", "dz2", "dxz", "dyz", "dx2-y2", "dxy"
      write(iow,'(A98)') '   ------------------------------------------------------------------------------------------------'
      write(iow,'(F10.4,A4,3F10.4,A4,5F10.4)')Real(W(10,10,i)),'|',Real(W(10,11:13,i)), '|', Real(W(10,14:18,i))
      write(iow,'(A98)') '   ------------------------------------------------------------------------------------------------'
      write(iow,'(F10.4,A4,3F10.4,A4,5F10.4)')(Real(W(10+j,10,i)),'|',Real(W(10+j,11:13,i)),'|',Real(W(10+j,14:18,i)),j=1,3)
      write(iow,'(A98)') '   ------------------------------------------------------------------------------------------------'
      write(iow,'(F10.4,A4,3F10.4,A4,5F10.4)')(Real(W(13+j,10,i)),'|',Real(W(13+j,11:13,i)),'|',Real(W(13+j,14:18,i)),j=1,5)
      write(iow,'(A98)') '   ================================================================================================'
      write(iow,'(9F16.10)') ((Real(W(k+9,j+9,i)), j=1,9), k=1,9)
      write(iow,'(A98)') '   ================================================================================================'
    End do
    Close(iow)

    Open(iow,file=trim(Solver)//'.out/'//'Hmn_'//trim(v)//'.dat',form='formatted',status='unknown')
    Do i = 1, 7
      Write(iow,'(A30)') "<iw|H(Ti0)|jw>"//trim(disp(i))//" real:"
      Write(iow,'(A10,A4,3A10,A4,5A10)') "s", "|", "pz", "px", "py", "|", "dz2", "dxz", "dyz", "dx2-y2", "dxy"
      write(iow,'(A98)') '   ------------------------------------------------------------------------------------------------'
      write(iow,'(F10.4,A4,3F10.4,A4,5F10.4)')Real(Hmn(10,10,i)),'|',Real(Hmn(10,11:13,i)), '|', Real(Hmn(10,14:18,i))
      write(iow,'(A98)') '   ------------------------------------------------------------------------------------------------'
      write(iow,'(F10.4,A4,3F10.4,A4,5F10.4)')(Real(Hmn(10+j,10,i)),'|',Real(Hmn(10+j,11:13,i)),'|',Real(Hmn(10+j,14:18,i)),j=1,3)
      write(iow,'(A98)') '   ------------------------------------------------------------------------------------------------'
      write(iow,'(F10.4,A4,3F10.4,A4,5F10.4)')(Real(Hmn(13+j,10,i)),'|',Real(Hmn(13+j,11:13,i)),'|',Real(Hmn(13+j,14:18,i)),j=1,5)
      write(iow,'(A98)') '   ================================================================================================'
      write(iow,'(9F16.10)') ((Real(Hmn(k+9,j+9,i)), j=1,9), k=1,9)
      write(iow,'(A98)') '   ================================================================================================'

      Write(iow,'(A30)') "<iw|H(Ti0)|jw>"//trim(disp(i))//" imag:"
      Write(iow,'(A10,A4,3A10,A4,5A10)') "s", "|", "pz", "px", "py", "|", "dz2", "dxz", "dyz", "dx2-y2", "dxy"
      write(iow,'(A98)') '   ------------------------------------------------------------------------------------------------'
      write(iow,'(F10.4,A4,3F10.4,A4,5F10.4)')Aimag(Hmn(10,10,i)),'|',Aimag(Hmn(10,11:13,i)), '|', Aimag(Hmn(10,14:18,i))
      write(iow,'(A98)') '   ------------------------------------------------------------------------------------------------'
      write(iow,'(F10.4,A4,3F10.4,A4,5F10.4)')(Aimag(Hmn(10+j,10,i)),'|',Aimag(Hmn(10+j,11:13,i)),'|',Aimag(Hmn(10+j,14:18,i)),j=1,3)
      write(iow,'(A98)') '   ------------------------------------------------------------------------------------------------'
      write(iow,'(F10.4,A4,3F10.4,A4,5F10.4)')(Aimag(Hmn(13+j,10,i)),'|',Aimag(Hmn(13+j,11:13,i)),'|',Aimag(Hmn(13+j,14:18,i)),j=1,5)
      write(iow,'(A98)') '   ================================================================================================'
      write(iow,'(9F16.10)') ((Aimag(Hmn(k+9,j+9,i)), j=1,9), k=1,9)
      write(iow,'(A98)') '   ================================================================================================'
    End do
    Close(iow)

  End Subroutine WriteWHmn

  Subroutine WriteGRij(ene,R,GRij,GRji)
    Implicit None
    Integer,intent(in)     :: R(3)
    Complex(dp),intent(in) :: ene
    Complex(dp),intent(in) :: GRij(NumWann,NumWann)
    Complex(dp),intent(in) :: GRji(NumWann,NumWann)
    Integer                :: iogw=1111, i, j

!    open(iogw,file=trim(Solver)//'.out/GWGW.dat',form='formatted',status='unknown')
    Write(*,'(2F10.4,3I10)') ene, R
    write(*,'(A98)') '   ------------------------------------------------------------------------------------------------'
    Write(*,'(A15,2F10.4,3I10,F16.10)') "G R real Ti", ene, R, Sum(Abs(Real(GRij(10:18,10:18))))
    Write(*,'(A10,A4,3A10,A4,5A10)') "s", "|", "pz", "px", "py", "|", "dz2", "dxz", "dyz", "dx2-y2", "dxy"
    write(*,'(A98)') '   ------------------------------------------------------------------------------------------------'
    write(*,'(F10.4,A4,3F10.4,A4,5F10.4)') Real(GRij(10,10)), '|',  Real(GRij(10,11:13)), '|', Real(GRij(10,14:18))
    write(*,'(A98)') '   ------------------------------------------------------------------------------------------------'
    write(*,'(F10.4,A4,3F10.4,A4,5F10.4)')(Real(GRij(10+i,10)),'|',Real(GRij(10+i,11:13)),'|',Real(GRij(10+i,14:18)),i=1,3)
    write(*,'(A98)') '   ------------------------------------------------------------------------------------------------'
    write(*,'(F10.4,A4,3F10.4,A4,5F10.4)')(Real(GRij(13+i,10)),'|',Real(GRij(13+i,11:13)),'|',Real(GRij(13+i,14:18)),i=1,5)
    write(*,'(A98)') '   ------------------------------------------------------------------------------------------------'

    Write(*,'(A15,2F10.4,3I10,F16.10)') "G R imag Ti", ene, R, Sum(Abs(Aimag(GRij(10:18,10:18))))
    Write(*,'(A10,A4,3A10,A4,5A10)') "s", "|", "pz", "px", "py", "|", "dz2", "dxz", "dyz", "dx2-y2", "dxy"
    write(*,'(A98)') '   ------------------------------------------------------------------------------------------------'
    write(*,'(F10.4,A4,3F10.4,A4,5F10.4)') Aimag(GRij(10,10)), '|', Aimag(GRij(10,11:13)), '|', Aimag(GRij(10,14:18))
    write(*,'(A98)') '   ------------------------------------------------------------------------------------------------'
    write(*,'(F10.4,A4,3F10.4,A4,5F10.4)')(Aimag(GRij(10+i,10)),'|',Aimag(GRij(10+i,11:13)),'|',Aimag(GRij(10+i,14:18)),i=1,3)
    write(*,'(A98)') '   ------------------------------------------------------------------------------------------------'
    write(*,'(F10.4,A4,3F10.4,A4,5F10.4)')(Aimag(GRij(13+i,10)),'|',Aimag(GRij(13+i,11:13)),'|',Aimag(GRij(13+i,14:18)),i=1,5)
    write(*,'(A98)') '   ------------------------------------------------------------------------------------------------'

    Write(*,'(A15,2F10.4,3I10,F16.10)') "G -R real Ti", ene, -R, Sum(Abs(Real(GRji(10:18,10:18))))
    Write(*,'(A10,A4,3A10,A4,5A10)') "s", "|", "pz", "px", "py", "|", "dz2", "dxz", "dyz", "dx2-y2", "dxy"
    write(*,'(A98)') '   ------------------------------------------------------------------------------------------------'
    write(*,'(F10.4,A4,3F10.4,A4,5F10.4)') Real(GRji(10,10)), '|',Real(GRji(10,11:13)), '|', Real(GRji(10,14:18))
    write(*,'(A98)') '   ------------------------------------------------------------------------------------------------'
    write(*,'(F10.4,A4,3F10.4,A4,5F10.4)')(Real(GRji(10+i,10)),'|',Real(GRji(10+i,11:13)),'|',Real(GRji(10+i,14:18)),i=1,3)
    write(*,'(A98)') '   ------------------------------------------------------------------------------------------------'
    write(*,'(F10.4,A4,3F10.4,A4,5F10.4)')(Real(GRji(13+i,10)),'|',Real(GRji(13+i,11:13)),'|',Real(GRji(13+i,14:18)),i=1,5)
    write(*,'(A98)') '   ------------------------------------------------------------------------------------------------'

    Write(*,'(A15,2F10.4,3I10,F16.10)') "G -R imag Ti", ene, -R, Sum(Abs(Aimag(GRji(10:18,10:18))))
    Write(*,'(A10,A4,3A10,A4,5A10)') "s", "|", "pz", "px", "py", "|", "dz2", "dxz", "dyz", "dx2-y2", "dxy"
    write(*,'(A98)') '   ------------------------------------------------------------------------------------------------'
    write(*,'(F10.4,A4,3F10.4,A4,5F10.4)') Aimag(GRji(10,10)), '|',Aimag(GRji(10,11:13)), '|', Aimag(GRji(10,14:18))
    write(*,'(A98)') '   ------------------------------------------------------------------------------------------------'
    write(*,'(F10.4,A4,3F10.4,A4,5F10.4)')(Aimag(GRji(10+i,10)),'|',Aimag(GRji(10+i,11:13)),'|',Aimag(GRji(10+i,14:18)),i=1,3)
    write(*,'(A98)') '   ------------------------------------------------------------------------------------------------'
    write(*,'(F10.4,A4,3F10.4,A4,5F10.4)')(Aimag(GRji(13+i,10)),'|',Aimag(GRji(13+i,11:13)),'|',Aimag(GRji(13+i,14:18)),i=1,5)
    write(*,'(A98)') '   ------------------------------------------------------------------------------------------------'

    write(*,'(A98)') '   ================================================================================================'
    write(*,'(9F16.10)') ((Real(GRij(i+9,j+9)), j=1,9), i=1,9)
    write(*,'(A98)') '   ------------------------------------------------------------------------------------------------'
    write(*,'(9F16.10)') ((Aimag(GRij(i+9,j+9)), j=1,9), i=1,9)
    write(*,'(A98)') '   ------------------------------------------------------------------------------------------------'
    write(*,'(9F16.10)') ((Real(GRji(i+9,j+9)), j=1,9), i=1,9)
    write(*,'(A98)') '   ------------------------------------------------------------------------------------------------'
    write(*,'(9F16.10)') ((Aimag(GRji(i+9,j+9)), j=1,9), i=1,9)
    write(*,'(A98)') '   ================================================================================================'

    Write(*,*) "                         "

!    close(iogw)

  End subroutine WriteGRij

  Subroutine WriteFijOrbital(a,b,R,Forb)
    Implicit None
    Integer,intent(in)     :: R(3)
    Complex(dp),intent(in) :: a, b
    Real(dp),intent(in)    :: Forb(NumWann,NumWann,3,3)
    Integer                :: iorbital=1111, i, j, ia, ib

    Write(*,'(A30)') "Fij orbital contributions:"
    Write(*,'(4F10.4,3I10)') a, b, R
    Do ia = 1, 3
      Do ib = 1, 3
        write(*,'(A98)') '   ------------------------------------------------------------------------------------------------'
        Write(*,'(A10,2F10.4,2I3,A10,F10.4)') 'energy:', (a+b)/2, ia, ib, 'F22: ', Sum(Forb(10:18,10:18,ia,ib))
        Write(*,'(A10,A4,3A10,A4,5A10)') "s", "|", "pz", "px", "py", "|", "dz2", "dxz", "dyz", "dx2-y2", "dxy"
        write(*,'(A98)') '   ------------------------------------------------------------------------------------------------'
        write(*,'(F10.4,A4,3F10.4,A4,5F10.4)') Forb(10,10,ia,ib), '|', Forb(10,11:13,ia,ib), '|', Forb(10,14:18,ia,ib)
        write(*,'(A98)') '   ------------------------------------------------------------------------------------------------'
        write(*,'(F10.4,A4,3F10.4,A4,5F10.4)')(Forb(10+i,10,ia,ib),'|',Forb(10+i,11:13,ia,ib),'|',Forb(10+i,14:18,ia,ib),i=1,3)
        write(*,'(A98)') '   ------------------------------------------------------------------------------------------------'
        write(*,'(F10.4,A4,3F10.4,A4,5F10.4)')(Forb(13+i,10,ia,ib),'|',Forb(13+i,11:13,ia,ib),'|',Forb(13+i,14:18,ia,ib),i=1,5)
        write(*,'(A98)') '   ------------------------------------------------------------------------------------------------'

        write(*,'(A98)') '   ================================================================================================'
        write(*,'(9F16.10)') ((Forb(i+9,j+9,ia,ib), j=1,9), i=1,9)
        write(*,'(A98)') '   ================================================================================================'

        Write(*,*) "                         "
      End do
    End do

  End subroutine WriteFijOrbital


End Module GreenFunction
