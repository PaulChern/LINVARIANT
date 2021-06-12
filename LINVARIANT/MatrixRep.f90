Module param
  Implicit none
  INTEGER :: NumAtoms=192, NumSpg=96, NumTgrp=1
End Module param

Module Aux
  Contains
  Function int2str(i) Result(str)
    integer, intent(in) :: i
    character(len=20)   :: str
  
    write (str, '(I5.5)') i
    str = trim(str)
  End function int2str

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

End Module Aux

Module invariant
  Contains
  Subroutine BasisTransform(basis,basisnew,latt,pos,rot,tran)
    Use param
    Use Aux
    Implicit None
    Real*8, Intent(out)    ::  basisnew(6,NumAtoms,NumAtoms*3)
    Real*8, Intent(in)     ::  basis(6,NumAtoms,NumAtoms*3)
    Real*8, Intent(in)     ::  rot(3,3), tran(3), latt(3,3), pos(3,NumAtoms)
    Real*8                 ::  basistemp(6,NumAtoms), dist
    Integer                ::  id, ia, ia1, ia2, i

    basisnew = 0.0D0

    Do id = 1, NumAtoms*3
      basistemp = 0.0D0
      Do ia = 1, NumAtoms
        basistemp(1:3,ia)=Modulo(rot3(rot, basis(1:3,ia,id))+tran(:),1.0D0)
        basistemp(4:6,ia)=rot3(rot, basis(4:6,ia,id))
      End do
      Do ia1 = 1, NumAtoms
        Do ia2 = 1, NumAtoms
          dist = 0.0D0
          Do i = 1, 3
            dist = dist + (basistemp(i,ia2)-pos(i,ia1))**2
          End do
          If(dist.lt.1.0D-12) then
            basisnew(1:3,ia1,id) = basistemp(1:3,ia2)
!            basisnew(1:3,ia1,id) = pos(1:3,ia1)
            basisnew(4:6,ia1,id) = basistemp(4:6,ia2)
            Cycle
          End if
        End do ! ia2
      End do ! ia1
    End do ! id
  End Subroutine BasisTransform

End Module invariant

Program MatrixRep
  Use param
  Use invariant
  Use Aux
  IMPLICIT NONE

  CHARACTER(len=200), DIMENSION(:), allocatable :: basislabel
  Real*8, DIMENSION(:),             allocatable :: normfactor, vecleft, vecright
  Real*8, DIMENSION(:,:),           allocatable :: OpMat
  Real*8, DIMENSION(:,:,:),         allocatable :: basis, basisnew
  Real*8, DIMENSION(:,:,:),         allocatable :: spg0, tgrp
  Real*8, DIMENSION(:,:),           allocatable :: pos
  Real*8                    :: dist, rot(3,3), tran(3), cell(3,3), latt(3,3)
  INTEGER                   :: ig, id, id1, id2, ia, ia1, ia2, i, j
  CHARACTER(len=200)        :: line

  allocate(basislabel(NumAtoms*3)) 
  allocate(normfactor(NumAtoms*3)) 
  allocate(vecleft(NumAtoms*3)) 
  allocate(vecright(NumAtoms*3)) 
  allocate(OpMat(NumAtoms*3,NumAtoms*3)) 
  allocate(basis(6,NumAtoms,NumAtoms*3)) 
  allocate(basisnew(6,NumAtoms,NumAtoms*3)) 
  allocate(pos(3,NumAtoms)) 
  allocate(spg0(4,3,NumSpg))
  allocate(tgrp(4,3,NumTgrp))

  latt = 0.0D0
  latt(1,1) = 24.9525D0
  latt(2,2) = 24.9525D0
  latt(3,3) = 24.9525D0

! read in basis
  Open(UNIT=11, FILE="/home/paulchern/Documents/Workshop/Working/Boracites/LINVARIANT/basis.dat")
  Do id = 1,NumAtoms*3
    Read(11,'(A)') line
    basislabel(id) = trim(line)
    Do ia = 1,NumAtoms
      Read(11,*) (basis(i,ia,id), i=1,6)
    End Do
  End Do
  Close(11)
!

! Read in symmetry operations
  Open(UNIT=22, FILE="/home/paulchern/Documents/Workshop/Working/Boracites/LINVARIANT/grp0.dat")
  Do ig = 1, NumSpg
    Read(22,'(A)') line
!    Write(*,*) trim(line)
    Do i = 1, 4
      Read(22,*) (spg0(i,j,ig), j=1,3)
    End do
  End Do
  Close(22)

  Open(UNIT=33, FILE="/home/paulchern/Documents/Workshop/Working/Boracites/LINVARIANT/tgrp.dat")
  Do ig = 1, NumTgrp
    Read(33,'(A)') line
!    Write(*,*) trim(line)
    Do i = 1, 4
      Read(33,*) (tgrp(i,j,ig), j=1,3)
    End do
  End Do
  Close(33)
!

  pos = basis(1:3,:,1)

  normfactor = 0.0D0
  Do id = 1, NumAtoms*3
    Do ia = 1, NumAtoms
      normfactor(id) = normfactor(id) + Norm2(rot3(latt,basis(4:6,ia,id)))**2
    End do
    normfactor(id)=1.0D0/Sqrt(normfactor(id))
  End do

  Do ig = 1, NumSpg
    rot = spg0(1:3,1:3,ig)
    tran = spg0(4,:,ig)
    Open(44,file="/home/paulchern/Documents/Workshop/Working/Boracites/LINVARIANT/basisnew-s"&
        //trim(int2str(ig))//".dat",status='unknown')
    Call BasisTransform(basis,basisnew,latt,pos,rot,tran)
    Do id = 1, NumAtoms*3
      write(44,*) trim(basislabel(id))
      write(44,'(6F15.8)') (basisnew(:,ia,id), ia=1,NumAtoms)
    End do
    Close(44)

    Do id1 = 1, NumAtoms*3
      Do id2 = 1, NumAtoms*3
        Do ia = 1, NumAtoms
          vecleft((ia-1)*3+1:3*ia) = rot3(latt,basisnew(4:6,ia,id1))
          vecright((ia-1)*3+1:3*ia) = rot3(latt,basis(4:6,ia,id2))
        End do
        OpMat(id1, id2) = sum(vecleft(:)*vecright(:))*normfactor(id1)*normfactor(id2)
      End do
    End do

    Open(55,file="/home/paulchern/Documents/Workshop/Working/Boracites/LINVARIANT/OpDispSpgMat-"&
        //trim(int2str(ig))//".dat",status='unknown')
    Do id1 = 1, NumAtoms*3
      Do id2 = 1, NumAtoms*3
        If(Abs(OpMat(id1,id2)).gt.10D-12) write(55,'(I10,I10,F15.8)') id1, id2, OpMat(id1,id2)
      End do
    End do
    Close(55)

  End do ! ig

  Do ig = 1, NumTgrp
    rot = tgrp(1:3,1:3,ig)
    tran = tgrp(4,:,ig)
    Open(44,file="/home/paulchern/Documents/Workshop/Working/Boracites/LINVARIANT/basisnew-t"&
        //trim(int2str(ig))//".dat",status='unknown')
    Call BasisTransform(basis,basisnew,latt,pos,rot,tran)
    Do id = 1, NumAtoms*3
      write(44,*) trim(basislabel(id))
      write(44,'(6F15.8)') (basisnew(:,ia,id), ia=1,NumAtoms)
    End do
    Close(44)
 
    Do id1 = 1, NumAtoms*3
      Do id2 = 1, NumAtoms*3
        Do ia = 1, NumAtoms
          vecleft((ia-1)*3+1:3*ia) = rot3(latt,basisnew(4:6,ia,id1))
          vecright((ia-1)*3+1:3*ia) = rot3(latt,basis(4:6,ia,id2))
        End do
        OpMat(id1, id2) = sum(vecleft(:)*vecright(:))*normfactor(id1)*normfactor(id2)
      End do
    End do
 
    Open(55,file="/home/paulchern/Documents/Workshop/Working/Boracites/LINVARIANT/OpDispTgrpMat-"&
        //trim(int2str(ig))//".dat",status='unknown')
    Do id1 = 1, NumAtoms*3
      Do id2 = 1, NumAtoms*3
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
  deallocate(pos)
  deallocate(spg0)
  deallocate(tgrp)

End Program MatrixRep
