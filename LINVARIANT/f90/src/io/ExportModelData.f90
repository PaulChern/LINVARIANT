Subroutine ExportMLModel(filename, Epot, avg, eta)
  use Constants
  use Parameters
  use Aux
  implicit none
  Real*8, Intent(in)           :: avg(FieldDim, NumField), eta(6)
  Real*8, Intent(in)           :: Epot
  character(len=*), intent(in) :: filename
  Real*8                       :: Coeff1D(470)
  Integer                      :: FileHandle = 3131
  Integer                      :: ifield, i
  
  
  Coeff1D = (/ CoeffGs,CoeffGp,CoeffGsp,CoeffHu,CoeffHp,CoeffHup,CoeffHsu,CoeffHsp,CoeffHsup,FieldCharge,epinf,mass,NoseMass /)
  open(FileHandle,file=trim(Solver)//'.out/'//filename,form='formatted',action='write',access='append',status='unknown')
  write(FileHandle,'(470F20.10)') Coeff1D
  write(FileHandle,'(F10.6)') Epot
  write(FileHandle,'(3F10.6)') eta(1), eta(2), eta(3)
  write(FileHandle,'(3F10.6)') eta(4), eta(5), eta(6)
  do ifield = 1, NumField
    write(FileHandle,'('//trim(int2str(FieldDim))//'F10.6)') avg(:,ifield)
  end do
  close(FileHandle)
  
End Subroutine ExportMLModel
