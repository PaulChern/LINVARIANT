Function GetEOnSiteEwald(ix, iy, iz, Fields) Result(ene)
  
  Implicit none
  Integer, Intent(in) :: ix, iy, iz
  Real*8,  Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Integer             :: i, j, jx, jy, jz, x, y, z, jfield, ifield
  Real*8              :: ene, ThreadEne
  
  ene = 0.0d0
  !$OMP    PARALLEL DEFAULT(SHARED) PRIVATE(jfield,x,y,z,ThreadEne)
    ThreadEne = 0.0d0
    
    !$OMP    DO COLLAPSE(1)
      do jz = 1, cgrid%n3
      do jy = 1, cgrid%n2
      do jx = 1, cgrid%n1
      do ifield = 1, NumField
      do j = 1, FieldDimList(ifield) ! jfield -> ifield
      do i = 1, FieldDimList(ifield)
        x = (jx-ix+1)-floor(real(jx-ix)/real(cgrid%n1))*cgrid%n1
        y = (jy-iy+1)-floor(real(jy-iy)/real(cgrid%n2))*cgrid%n2
        z = (jz-iz+1)-floor(real(jz-iz)/real(cgrid%n3))*cgrid%n3
        ! jfield -> ifield
        ThreadEne = ThreadEne &
          + 0.5*FieldCharge(ifield)*FieldCharge(ifield)*Fields(i,ifield,ix,iy,iz)*EwaldMat(i,j,x,y,z)*Fields(j,ifield,jx,jy,jz)
      end do ! j
      end do ! j
      end do ! i
      end do ! jz
      end do ! jy
      end do ! jx
    !$OMP    END DO
    
    !$OMP CRITICAL
      ene = ene + ThreadEne
    !$OMP END CRITICAL
    
  !$OMP    END PARALLEL
  
    do ifield = 1, NumField - 1
      do j = 1, FieldDimList(ifield) ! jfield -> ifield
        do i = 1, FieldDimList(ifield)
        ! jfield -> ifield
          ene = ene &
            + 0.5*FieldCharge(ifield)*FieldCharge(ifield)*Fields(i,ifield,ix,iy,iz)*EwaldMat(i,j,1,1,1)*Fields(j,ifield,ix,iy,iz)
        end do ! j
      end do ! i
    end do ! 
  
End Function GetEOnSiteEwald
