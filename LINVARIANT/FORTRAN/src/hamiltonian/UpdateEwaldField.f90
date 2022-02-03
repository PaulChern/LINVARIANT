Subroutine UpdateEwaldField(ix, iy, iz, EwaldField, idelta, delta)
  Use omp_lib
  use Parameters
  
  Implicit none
  integer, intent(in)     :: ix, iy, iz, idelta
  real*8,  intent(inout)  :: EwaldField(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  real*8,  intent(in)     :: delta(FieldDimList(idelta))
  real*8                  :: dfield(cgrid%npts*3)
  integer                 :: id
  
  dfield = 0.0D0
  id = (iz-1)*cgrid%n2*cgrid%n1*FieldDimList(idelta) &
            + (iy-1)*cgrid%n1*FieldDimList(idelta) &
            + (ix-1)*FieldDimList(idelta)

  call dgemv('N', 3*cgrid%npts, 3, &
             1.0D0, EwaldHessian(:,id+1:id+3,idelta), 3*cgrid%npts, &
             delta, 1, 0.0D0, dfield, 1)
  EwaldField(:,idelta,:,:,:) = EwaldField(:,idelta,:,:,:) + Reshape(dfield, (/3, cgrid%n1, cgrid%n2, cgrid%n3/))
  
End Subroutine UpdateEwaldField
