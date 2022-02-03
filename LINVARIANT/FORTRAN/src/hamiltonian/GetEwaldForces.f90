Function GetEwaldForces(Fields) Result(EwaldForce)
  
  implicit none
  real*8,  intent(in)           :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  integer                       :: i, ifield
  integer                       :: j, jfield, jx, jy, jz
  integer                       :: ix, iy, iz
  integer                       :: x, y, z
  real*8                        :: EwaldForce(Max(FieldDim, 6),NumField+1, cgrid%n1, cgrid%n2, cgrid%n3)
  real*8                        :: force(3*cgrid%npts,NumField - 1)
  real*8                        :: vector(3*cgrid%npts,NumField - 1)
  
  EwaldForce = 0.0D0
  force = 0.0D0

  do jfield = 1, NumField - 1
    vector(:,jfield) = FieldTo1D(Fields, jfield)
    call dgemv('N', 3*cgrid%npts, 3*cgrid%npts, &
               -1.0D0, EwaldHessian(:,:,jfield), 3*cgrid%npts, &
               vector(:,jfield), 1, 0.0D0, force(:,jfield), 1)
    EwaldForce(1:3,jfield,:,:,:) = Reshape(force(:,jfield), (/3, cgrid%n1, cgrid%n2, cgrid%n3/))
  end do

End Function GetEwaldForces
