Function GetDeltaHEwald(x0, y0, z0, EwaldField, idelta, delta) Result(ene)
  
  Implicit none
  Integer, Intent(in) :: x0, y0, z0, idelta
  Real*8,  Intent(in) :: EwaldField(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(in) :: delta(FieldDimList(idelta))
  Integer             :: i, j, ifield, jfield
  Real*8              :: ene, ThreadEne
  
  ene = 0.0d0
  
  do i = 1, FieldDimList(idelta)
    ene = ene + delta(i)*EwaldField(i,idelta,x0,y0,z0)
  end do ! i
  
  do i = 1, FieldDimList(idelta)
    do j = 1, FieldDimList(idelta)
      ene = ene + 0.5*FieldCharge(idelta)**2*delta(i)*EwaldMat(j,i,1,1,1)*delta(j)
    end do ! j
  end do ! i
  
End Function GetDeltaHEwald
