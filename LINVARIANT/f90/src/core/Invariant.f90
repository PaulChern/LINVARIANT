Module invariant
  Use Parameters
  Use Inputs
  Use Aux
  Use mma

  Contains
  Subroutine BasisTransform(basis,basisnew,sites,rot,tran)

    Implicit None
    Real*8, Intent(out)    ::  basisnew(6,unitcell%nion,3*unitcell%nion)
    Real*8, Intent(in)     ::  basis(6,unitcell%nion,3*unitcell%nion)
    Real*8, Intent(in)     ::  rot(3,3), tran(3), sites(3,unitcell%nion)
    Real*8                 ::  basistemp(6,unitcell%nion), dist
    Integer                ::  id, ia, ia1, ia2, i

    basisnew = 0.0D0

    Do id = 1, 3*unitcell%nion
      basistemp = 0.0D0
      Do ia = 1, unitcell%nion
        basistemp(1:3,ia)=Modulo(rot3(rot, basis(1:3,ia,id))+tran(:),1.0D0)
        basistemp(4:6,ia)=rot3(rot, basis(4:6,ia,id))
      End do
      Do ia1 = 1, unitcell%nion
        Do ia2 = 1, unitcell%nion
          dist = 0.0D0
          Do i = 1, 3
            dist = dist + (basistemp(i,ia2)-sites(i,ia1))**2
          End do
          If(dist.lt.1.0D-12) then
            basisnew(1:3,ia1,id) = basistemp(1:3,ia2)
!            basisnew(1:3,ia1,id) = sites(1:3,ia1)
            basisnew(4:6,ia1,id) = basistemp(4:6,ia2)
            Cycle
          End if
        End do ! ia2
      End do ! ia1
    End do ! id
  End Subroutine BasisTransform

  Subroutine MatrixRep
    IMPLICIT NONE

    CHARACTER(len=200), DIMENSION(:), allocatable :: basislabel
    Real*8, DIMENSION(:),             allocatable :: normfactor, vecleft, vecright
    Real*8, DIMENSION(:,:),           allocatable :: OpMat
    Real*8, DIMENSION(:,:,:),         allocatable :: basis, basisnew
    Real*8, DIMENSION(:,:,:),         allocatable :: spg0, tgrp
    Real*8, DIMENSION(:,:),           allocatable :: sites
    Real*8                    :: dist, rot(3,3), tran(3), cell(3,3)
    INTEGER                   :: ig, id, id1, id2, ia, ia1, ia2, i, j
    INTEGER                   :: NumSpg, NumTgrp
    CHARACTER(len=200)        :: line
  
    allocate(basislabel(3*unitcell%nion))
    allocate(normfactor(3*unitcell%nion))
    allocate(vecleft(3*unitcell%nion))
    allocate(vecright(3*unitcell%nion))
    allocate(OpMat(3*unitcell%nion,3*unitcell%nion))
    allocate(basis(6,unitcell%nion,3*unitcell%nion))
    allocate(basisnew(6,unitcell%nion,3*unitcell%nion))
    allocate(sites(3,unitcell%nion))

!   read in basis
    Open(UNIT=11, FILE='basis.dat')
    Do id = 1,3*unitcell%nion
      Read(11,'(A)') line
      basislabel(id) = trim(line)
      Do ia = 1,unitcell%nion
        Read(11,*) (basis(i,ia,id), i=1,6)
      End Do
    End Do
    Close(11)
!

!   Read in symmetry operations
    Open(UNIT=22, FILE='spg0.dat')
    Read(22,*) NumSpg
    allocate(spg0(4,3,NumSpg))
    Do ig = 1, NumSpg
      Read(22,'(A)') line
!      Write(*,*) trim(line)
      Do i = 1, 4
        Read(22,*) (spg0(i,j,ig), j=1,3)
      End do
    End Do
    Close(22)

    Open(UNIT=33, FILE='tgrp.dat')
    Read(33,*) NumTgrp
    allocate(tgrp(4,3,NumTgrp))
    Do ig = 1, NumTgrp
      Read(33,'(A)') line
!      Write(*,*) trim(line)
      Do i = 1, 4
        Read(33,*) (tgrp(i,j,ig), j=1,3)
      End do
    End Do
    Close(33)
!

    sites = basis(1:3,:,1)

    normfactor = 0.0D0
    Do id = 1, 3*unitcell%nion
      Do ia = 1, unitcell%nion
        normfactor(id) = normfactor(id) + Norm2(rot3(unitcell%latt_a,basis(4:6,ia,id)))**2
      End do
      normfactor(id)=1.0D0/Sqrt(normfactor(id))
    End do

    Do ig = 1, NumSpg
      rot = spg0(1:3,1:3,ig)
      tran = spg0(4,:,ig)
      Open(44,file=trim(Solver)//'.out/basisnew-s'//trim(int2str5(ig))//'.dat',status='unknown')
      Call BasisTransform(basis,basisnew,sites,rot,tran)
      Do id = 1, 3*unitcell%nion
        write(44,*) trim(basislabel(id))
        write(44,'(6F15.8)') (basisnew(:,ia,id), ia=1,unitcell%nion)
      End do
      Close(44)

      Do id1 = 1, 3*unitcell%nion
        Do id2 = 1, 3*unitcell%nion
          Do ia = 1, unitcell%nion
            vecleft((ia-1)*3+1:3*ia) = rot3(unitcell%latt_a,basisnew(4:6,ia,id1))
            vecright((ia-1)*3+1:3*ia) = rot3(unitcell%latt_a,basis(4:6,ia,id2))
          End do
          OpMat(id1, id2) = sum(vecleft(:)*vecright(:))*normfactor(id1)*normfactor(id2)
        End do
      End do

      Open(55,file=trim(Solver)//'.out/OpDispSpgMat-'//trim(int2str5(ig))//'.dat',status='unknown')
      Do id1 = 1, 3*unitcell%nion
        Do id2 = 1, 3*unitcell%nion
          If(Abs(OpMat(id1,id2)).gt.10D-12) write(55,'(I10,I10,F15.8)') id1, id2, OpMat(id1,id2)
        End do
      End do
      Close(55)

    End do ! ig

    Do ig = 1, NumTgrp
      rot = tgrp(1:3,1:3,ig)
      tran = tgrp(4,:,ig)
      Open(44,file=trim(Solver)//'.out/basisnew-t'//trim(int2str5(ig))//'.dat',status='unknown')
      Call BasisTransform(basis,basisnew,sites,rot,tran)
      Do id = 1, 3*unitcell%nion
        write(44,*) trim(basislabel(id))
        write(44,'(6F15.8)') (basisnew(:,ia,id), ia=1,unitcell%nion)
      End do
      Close(44)
 
      Do id1 = 1, 3*unitcell%nion
        Do id2 = 1, 3*unitcell%nion
          Do ia = 1, unitcell%nion
            vecleft((ia-1)*3+1:3*ia) = rot3(unitcell%latt_a,basisnew(4:6,ia,id1))
            vecright((ia-1)*3+1:3*ia) = rot3(unitcell%latt_a,basis(4:6,ia,id2))
          End do
          OpMat(id1, id2) = sum(vecleft(:)*vecright(:))*normfactor(id1)*normfactor(id2)
        End do
      End do
 
      Open(55,file=trim(Solver)//'.out/OpDispTgrpMat-'//trim(int2str5(ig))//'.dat',status='unknown')
      Do id1 = 1, 3*unitcell%nion
        Do id2 = 1, 3*unitcell%nion
          If(Abs(OpMat(id1,id2)).gt.10D-12) write(55,'(I10,I10,F15.8)') id1, id2, OpMat(id1,id2)
        End do
      End do
      Close(55)
 
    End do ! ig

    deallocate(basislabel)
    deallocate(normfactor)
    deallocate(vecleft)
    deallocate(vecright)
    deallocate(OpMat)
    deallocate(basis)
    deallocate(basisnew)
    deallocate(sites)
    deallocate(spg0)
    deallocate(tgrp)


  End Subroutine MatrixRep
End Module invariant
