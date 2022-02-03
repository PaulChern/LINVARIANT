Function GetDeltaHEps(x0, y0, z0, Fields, e0ij, idelta, delta) Result(DeltaH)
  
  Implicit none
  Integer, Intent(in) :: x0, y0, z0, idelta
  Real*8,  Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(in) :: delta(FieldDimList(idelta))
  Real*8,  Intent(in) :: e0ij(3,3)
  Real*8              :: eij(3,3), euij(3,3)
  Real*8              :: DeltaH
  Real*8              :: DeltaHEps
  
  Integer             :: x1,y1,z1,xi1,yi1,zi1
  
  x1 = (x0+1)-floor(real(x0+1-1)/real(cgrid%n1))*cgrid%n1
  y1 = (y0+1)-floor(real(y0+1-1)/real(cgrid%n2))*cgrid%n2
  z1 = (z0+1)-floor(real(z0+1-1)/real(cgrid%n3))*cgrid%n3
  xi1 = (x0-1)-floor(real(x0-1-1)/real(cgrid%n1))*cgrid%n1
  yi1 = (y0-1)-floor(real(y0-1-1)/real(cgrid%n2))*cgrid%n2
  zi1 = (z0-1)-floor(real(z0-1-1)/real(cgrid%n3))*cgrid%n3
  
  DeltaH = 0.0D0
  DeltaHEps = 0.0D0
  euij = StrainFromu(x0, y0, z0, Fields)
  eij = e0ij + euij
  
  SELECT CASE (idelta)
  CASE (1)
    DeltaHEps = &
      DeltaHEps &
      +0.
    
    DeltaH=DeltaHEps
  
  CASE (2)
    DeltaHEps = &
      DeltaHEps &
      +0.
    
    DeltaH=DeltaHEps
  
  CASE (3)
    DeltaHEps = &
      DeltaHEps &
      +0.
    
    DeltaH=DeltaHEps
  
  CASE (4)
    DeltaHEps = &
      DeltaHEps &
      +CoeffEps(1)*delta(1) &
      +0.5*CoeffEps(4)*delta(1)**2 &
      +CoeffEps(1)*delta(2) &
      +CoeffEps(3)*delta(1)*delta(2) &
      +0.5*CoeffEps(4)*delta(2)**2 &
      +CoeffEps(1)*delta(3) &
      +CoeffEps(3)*delta(1)*delta(3) &
      +CoeffEps(3)*delta(2)*delta(3) &
      +0.5*CoeffEps(4)*delta(3)**2 &
      +0.5*CoeffEps(2)*delta(4)**2 &
      +0.5*CoeffEps(2)*delta(5)**2 &
      +0.5*CoeffEps(2)*delta(6)**2 &
      +CoeffEps(4)*delta(1)*e0ij(1,1) &
      +CoeffEps(3)*delta(2)*e0ij(1,1) &
      +CoeffEps(3)*delta(3)*e0ij(1,1) &
      +CoeffEps(2)*delta(6)*e0ij(1,2) &
      +CoeffEps(2)*delta(5)*e0ij(1,3) &
      +CoeffEps(3)*delta(1)*e0ij(2,2) &
      +CoeffEps(4)*delta(2)*e0ij(2,2) &
      +CoeffEps(3)*delta(3)*e0ij(2,2) &
      +CoeffEps(2)*delta(4)*e0ij(2,3) &
      +CoeffEps(3)*delta(1)*e0ij(3,3) &
      +CoeffEps(3)*delta(2)*e0ij(3,3) &
      +CoeffEps(4)*delta(3)*e0ij(3,3)
    
    DeltaH=DeltaHEps
  
  CASE DEFAULT
      write(*,*) "mode out of range!"
      call abort
  END SELECT
  
End Function GetDeltaHEps
