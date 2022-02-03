Function GetDeltaHEpsu(x0, y0, z0, Fields, e0ij, idelta, delta) Result(DeltaH)
  
  Implicit none
  Integer, Intent(in) :: x0, y0, z0, idelta
  Real*8,  Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(in) :: delta(FieldDimList(idelta))
  Real*8,  Intent(in) :: e0ij(3,3)
  Real*8              :: eij(3,3), euij(3,3)
  Real*8              :: DeltaH
  Real*8              :: DeltaHEpsu
  
  Integer             :: x1,y1,z1,xi1,yi1,zi1
  
  x1 = (x0+1)-floor(real(x0+1-1)/real(cgrid%n1))*cgrid%n1
  y1 = (y0+1)-floor(real(y0+1-1)/real(cgrid%n2))*cgrid%n2
  z1 = (z0+1)-floor(real(z0+1-1)/real(cgrid%n3))*cgrid%n3
  xi1 = (x0-1)-floor(real(x0-1-1)/real(cgrid%n1))*cgrid%n1
  yi1 = (y0-1)-floor(real(y0-1-1)/real(cgrid%n2))*cgrid%n2
  zi1 = (z0-1)-floor(real(z0-1-1)/real(cgrid%n3))*cgrid%n3
  
  DeltaH = 0.0D0
  DeltaHEpsu = 0.0D0
  euij = StrainFromu(x0, y0, z0, Fields)
  eij = e0ij + euij
  
  SELECT CASE (idelta)
  CASE (1)
    DeltaHEpsu = &
      DeltaHEpsu &
      +0.
    
    DeltaH=DeltaHEpsu
  
  CASE (2)
    DeltaHEpsu = &
      DeltaHEpsu &
      +0.
    
    DeltaH=DeltaHEpsu
  
  CASE (3)
    DeltaHEpsu = &
      DeltaHEpsu &
      +0. &
      -0.25*CoeffEpsu(4)*delta(1)*e0ij(1,1) &
      -0.25*CoeffEpsu(3)*delta(2)*e0ij(1,1) &
      -0.25*CoeffEpsu(3)*delta(3)*e0ij(1,1) &
      -0.25*CoeffEpsu(2)*delta(1)*e0ij(1,2) &
      -0.25*CoeffEpsu(2)*delta(2)*e0ij(1,2) &
      -0.25*CoeffEpsu(2)*delta(1)*e0ij(1,3) &
      -0.25*CoeffEpsu(2)*delta(3)*e0ij(1,3) &
      -0.25*CoeffEpsu(3)*delta(1)*e0ij(2,2) &
      -0.25*CoeffEpsu(4)*delta(2)*e0ij(2,2) &
      -0.25*CoeffEpsu(3)*delta(3)*e0ij(2,2) &
      -0.25*CoeffEpsu(2)*delta(2)*e0ij(2,3) &
      -0.25*CoeffEpsu(2)*delta(3)*e0ij(2,3) &
      -0.25*CoeffEpsu(3)*delta(1)*e0ij(3,3) &
      -0.25*CoeffEpsu(3)*delta(2)*e0ij(3,3) &
      -0.25*CoeffEpsu(4)*delta(3)*e0ij(3,3)
    
    DeltaH=DeltaHEpsu
  
  CASE (4)
    DeltaHEpsu = &
      DeltaHEpsu &
      +0.
    
    DeltaH=DeltaHEpsu
  
  CASE DEFAULT
      write(*,*) "mode out of range!"
      call abort
  END SELECT
  
End Function GetDeltaHEpsu
