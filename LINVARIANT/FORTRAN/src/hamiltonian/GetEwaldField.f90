Subroutine GetEwaldField(Fields, EwaldField)
  
  implicit none
  real*8,  intent(in)           :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  real*8,  intent(inout)        :: EwaldField(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  integer                       :: i, j, ifield, jfield, ix, iy, iz, jx, jy, jz, x, y, z
  
  EwaldField = 0.0D0
  
  !$OMP    PARALLEL  DEFAULT(SHARED) PRIVATE(x,y,z)
    do jz = 1, cgrid%n3
    do jy = 1, cgrid%n2
    do jx = 1, cgrid%n1
    do jfield = 1, NumField - 1
    do j = 1, FieldDimList(jfield)
  !$OMP    DO COLLAPSE(1)
    do iz = 1, cgrid%n3
    do iy = 1, cgrid%n2
    do ix = 1, cgrid%n1
    do i = 1, FieldDimList(jfield) ! ifield -> jfield
      Call Pbc(jx-ix, jy-iy, jz-iz, x, y, z)
      ! ifield -> jfield
      EwaldField(i,jfield,ix,iy,iz) = EwaldField(i,jfield,ix,iy,iz) &
        + FieldCharge(jfield)*FieldCharge(jfield)*EwaldMat(j,i,x,y,z)*Fields(j,jfield,jx,jy,jz)
    end do
    end do
    end do
    end do
  !$OMP    END DO
    end do
    end do
    end do
    end do
    end do
  !$OMP    END PARALLEL
  
  
End Subroutine GetEwaldField
