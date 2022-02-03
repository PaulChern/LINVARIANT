Function GetDeltaHEpsDisp(x0, y0, z0, Fields, e0ij, idelta, delta) Result(DeltaH)
  
  Implicit none
  Integer, Intent(in) :: x0, y0, z0, idelta
  Real*8,  Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(in) :: delta(FieldDimList(idelta))
  Real*8,  Intent(in) :: e0ij(3,3)
  Real*8              :: eij(3,3), euij(3,3)
  Real*8              :: DeltaH
  Real*8              :: DeltaHEpsDisp
  
  Integer             :: x1,y1,z1,xi1,yi1,zi1
  
  x1 = (x0+1)-floor(real(x0+1-1)/real(cgrid%n1))*cgrid%n1
  y1 = (y0+1)-floor(real(y0+1-1)/real(cgrid%n2))*cgrid%n2
  z1 = (z0+1)-floor(real(z0+1-1)/real(cgrid%n3))*cgrid%n3
  xi1 = (x0-1)-floor(real(x0-1-1)/real(cgrid%n1))*cgrid%n1
  yi1 = (y0-1)-floor(real(y0-1-1)/real(cgrid%n2))*cgrid%n2
  zi1 = (z0-1)-floor(real(z0-1-1)/real(cgrid%n3))*cgrid%n3
  
  DeltaH = 0.0D0
  DeltaHEpsDisp = 0.0D0
  euij = StrainFromu(x0, y0, z0, Fields)
  eij = e0ij + euij
  
  SELECT CASE (idelta)
  CASE (1)
    DeltaHEpsDisp = &
      DeltaHEpsDisp &
      +CoeffEpsDisp(5)*delta(1)**2*e0ij(1,1) &
      +CoeffEpsDisp(4)*delta(2)**2*e0ij(1,1) &
      +CoeffEpsDisp(4)*delta(3)**2*e0ij(1,1) &
      +CoeffEpsDisp(1)*delta(1)*delta(2)*e0ij(1,2) &
      +CoeffEpsDisp(1)*delta(1)*delta(3)*e0ij(1,3) &
      +CoeffEpsDisp(4)*delta(1)**2*e0ij(2,2) &
      +CoeffEpsDisp(5)*delta(2)**2*e0ij(2,2) &
      +CoeffEpsDisp(4)*delta(3)**2*e0ij(2,2) &
      +CoeffEpsDisp(1)*delta(2)*delta(3)*e0ij(2,3) &
      +CoeffEpsDisp(4)*delta(1)**2*e0ij(3,3) &
      +CoeffEpsDisp(4)*delta(2)**2*e0ij(3,3) &
      +CoeffEpsDisp(5)*delta(3)**2*e0ij(3,3) &
      +CoeffEpsDisp(5)*delta(1)**2*euij(1,1) &
      +CoeffEpsDisp(4)*delta(2)**2*euij(1,1) &
      +CoeffEpsDisp(4)*delta(3)**2*euij(1,1) &
      +CoeffEpsDisp(1)*delta(1)*delta(2)*euij(1,2) &
      +CoeffEpsDisp(1)*delta(1)*delta(3)*euij(1,3) &
      +CoeffEpsDisp(4)*delta(1)**2*euij(2,2) &
      +CoeffEpsDisp(5)*delta(2)**2*euij(2,2) &
      +CoeffEpsDisp(4)*delta(3)**2*euij(2,2) &
      +CoeffEpsDisp(1)*delta(2)*delta(3)*euij(2,3) &
      +CoeffEpsDisp(4)*delta(1)**2*euij(3,3) &
      +CoeffEpsDisp(4)*delta(2)**2*euij(3,3) &
      +CoeffEpsDisp(5)*delta(3)**2*euij(3,3) &
      +2.*CoeffEpsDisp(5)*delta(1)*e0ij(1,1)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(1)*delta(2)*e0ij(1,2)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(1)*delta(3)*e0ij(1,3)*Fields(1,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(4)*delta(1)*e0ij(2,2)*Fields(1,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(4)*delta(1)*e0ij(3,3)*Fields(1,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(5)*delta(1)*euij(1,1)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(1)*delta(2)*euij(1,2)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(1)*delta(3)*euij(1,3)*Fields(1,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(4)*delta(1)*euij(2,2)*Fields(1,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(4)*delta(1)*euij(3,3)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(8)*delta(1)*e0ij(1,1)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(2)*e0ij(1,2)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(3)*e0ij(1,3)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(1)*e0ij(2,2)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(1)*e0ij(3,3)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(8)*delta(1)*euij(1,1)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(2)*euij(1,2)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(3)*euij(1,3)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(1)*euij(2,2)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(1)*euij(3,3)*Fields(1,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(4)*delta(2)*e0ij(1,1)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(1)*delta(1)*e0ij(1,2)*Fields(2,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(5)*delta(2)*e0ij(2,2)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(1)*delta(3)*e0ij(2,3)*Fields(2,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(4)*delta(2)*e0ij(3,3)*Fields(2,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(4)*delta(2)*euij(1,1)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(1)*delta(1)*euij(1,2)*Fields(2,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(5)*delta(2)*euij(2,2)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(1)*delta(3)*euij(2,3)*Fields(2,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(4)*delta(2)*euij(3,3)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(2)*e0ij(1,1)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(1)*e0ij(1,2)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(8)*delta(2)*e0ij(2,2)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(3)*e0ij(2,3)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(2)*e0ij(3,3)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(2)*euij(1,1)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(1)*euij(1,2)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(8)*delta(2)*euij(2,2)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(3)*euij(2,3)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(2)*euij(3,3)*Fields(2,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(4)*delta(3)*e0ij(1,1)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(1)*delta(1)*e0ij(1,3)*Fields(3,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(4)*delta(3)*e0ij(2,2)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(1)*delta(2)*e0ij(2,3)*Fields(3,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(5)*delta(3)*e0ij(3,3)*Fields(3,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(4)*delta(3)*euij(1,1)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(1)*delta(1)*euij(1,3)*Fields(3,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(4)*delta(3)*euij(2,2)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(1)*delta(2)*euij(2,3)*Fields(3,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(5)*delta(3)*euij(3,3)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(3)*e0ij(1,1)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(1)*e0ij(1,3)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(3)*e0ij(2,2)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(2)*e0ij(2,3)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(8)*delta(3)*e0ij(3,3)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(3)*euij(1,1)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(1)*euij(1,3)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(3)*euij(2,2)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(2)*euij(2,3)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(8)*delta(3)*euij(3,3)*Fields(3,2,x0,y0,z0)
    
    DeltaH=DeltaHEpsDisp
  
  CASE (2)
    DeltaHEpsDisp = &
      DeltaHEpsDisp &
      +CoeffEpsDisp(9)*delta(1)**2*e0ij(1,1) &
      +CoeffEpsDisp(7)*delta(2)**2*e0ij(1,1) &
      +CoeffEpsDisp(7)*delta(3)**2*e0ij(1,1) &
      +CoeffEpsDisp(3)*delta(1)*delta(2)*e0ij(1,2) &
      +CoeffEpsDisp(3)*delta(1)*delta(3)*e0ij(1,3) &
      +CoeffEpsDisp(7)*delta(1)**2*e0ij(2,2) &
      +CoeffEpsDisp(9)*delta(2)**2*e0ij(2,2) &
      +CoeffEpsDisp(7)*delta(3)**2*e0ij(2,2) &
      +CoeffEpsDisp(3)*delta(2)*delta(3)*e0ij(2,3) &
      +CoeffEpsDisp(7)*delta(1)**2*e0ij(3,3) &
      +CoeffEpsDisp(7)*delta(2)**2*e0ij(3,3) &
      +CoeffEpsDisp(9)*delta(3)**2*e0ij(3,3) &
      +CoeffEpsDisp(9)*delta(1)**2*euij(1,1) &
      +CoeffEpsDisp(7)*delta(2)**2*euij(1,1) &
      +CoeffEpsDisp(7)*delta(3)**2*euij(1,1) &
      +CoeffEpsDisp(3)*delta(1)*delta(2)*euij(1,2) &
      +CoeffEpsDisp(3)*delta(1)*delta(3)*euij(1,3) &
      +CoeffEpsDisp(7)*delta(1)**2*euij(2,2) &
      +CoeffEpsDisp(9)*delta(2)**2*euij(2,2) &
      +CoeffEpsDisp(7)*delta(3)**2*euij(2,2) &
      +CoeffEpsDisp(3)*delta(2)*delta(3)*euij(2,3) &
      +CoeffEpsDisp(7)*delta(1)**2*euij(3,3) &
      +CoeffEpsDisp(7)*delta(2)**2*euij(3,3) &
      +CoeffEpsDisp(9)*delta(3)**2*euij(3,3) &
      +CoeffEpsDisp(8)*delta(1)*e0ij(1,1)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(2)*e0ij(1,2)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(3)*e0ij(1,3)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(1)*e0ij(2,2)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(1)*e0ij(3,3)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(8)*delta(1)*euij(1,1)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(2)*euij(1,2)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(3)*euij(1,3)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(1)*euij(2,2)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(1)*euij(3,3)*Fields(1,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(9)*delta(1)*e0ij(1,1)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(2)*e0ij(1,2)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(3)*e0ij(1,3)*Fields(1,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(7)*delta(1)*e0ij(2,2)*Fields(1,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(7)*delta(1)*e0ij(3,3)*Fields(1,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(9)*delta(1)*euij(1,1)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(2)*euij(1,2)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(3)*euij(1,3)*Fields(1,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(7)*delta(1)*euij(2,2)*Fields(1,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(7)*delta(1)*euij(3,3)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(2)*e0ij(1,1)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(1)*e0ij(1,2)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(8)*delta(2)*e0ij(2,2)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(3)*e0ij(2,3)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(2)*e0ij(3,3)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(2)*euij(1,1)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(1)*euij(1,2)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(8)*delta(2)*euij(2,2)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(3)*euij(2,3)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(2)*euij(3,3)*Fields(2,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(7)*delta(2)*e0ij(1,1)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(1)*e0ij(1,2)*Fields(2,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(9)*delta(2)*e0ij(2,2)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(3)*e0ij(2,3)*Fields(2,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(7)*delta(2)*e0ij(3,3)*Fields(2,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(7)*delta(2)*euij(1,1)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(1)*euij(1,2)*Fields(2,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(9)*delta(2)*euij(2,2)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(3)*euij(2,3)*Fields(2,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(7)*delta(2)*euij(3,3)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(3)*e0ij(1,1)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(1)*e0ij(1,3)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(3)*e0ij(2,2)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(2)*e0ij(2,3)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(8)*delta(3)*e0ij(3,3)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(3)*euij(1,1)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(1)*euij(1,3)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(3)*euij(2,2)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(2)*euij(2,3)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(8)*delta(3)*euij(3,3)*Fields(3,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(7)*delta(3)*e0ij(1,1)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(1)*e0ij(1,3)*Fields(3,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(7)*delta(3)*e0ij(2,2)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(2)*e0ij(2,3)*Fields(3,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(9)*delta(3)*e0ij(3,3)*Fields(3,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(7)*delta(3)*euij(1,1)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(1)*euij(1,3)*Fields(3,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(7)*delta(3)*euij(2,2)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(2)*euij(2,3)*Fields(3,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(9)*delta(3)*euij(3,3)*Fields(3,2,x0,y0,z0)
    
    DeltaH=DeltaHEpsDisp
  
  CASE (3)
    DeltaHEpsDisp = &
      DeltaHEpsDisp &
      +0.
    
    DeltaH=DeltaHEpsDisp
  
  CASE (4)
    DeltaHEpsDisp = &
      DeltaHEpsDisp &
      +CoeffEpsDisp(5)*delta(1)*Fields(1,1,x0,y0,z0)**2 &
      +CoeffEpsDisp(4)*delta(2)*Fields(1,1,x0,y0,z0)**2 &
      +CoeffEpsDisp(4)*delta(3)*Fields(1,1,x0,y0,z0)**2 &
      +CoeffEpsDisp(8)*delta(1)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(2)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(3)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(9)*delta(1)*Fields(1,2,x0,y0,z0)**2 &
      +CoeffEpsDisp(7)*delta(2)*Fields(1,2,x0,y0,z0)**2 &
      +CoeffEpsDisp(7)*delta(3)*Fields(1,2,x0,y0,z0)**2 &
      +CoeffEpsDisp(1)*delta(6)*Fields(1,1,x0,y0,z0)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(6)*Fields(1,2,x0,y0,z0)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(4)*delta(1)*Fields(2,1,x0,y0,z0)**2 &
      +CoeffEpsDisp(5)*delta(2)*Fields(2,1,x0,y0,z0)**2 &
      +CoeffEpsDisp(4)*delta(3)*Fields(2,1,x0,y0,z0)**2 &
      +CoeffEpsDisp(2)*delta(6)*Fields(1,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(6)*Fields(1,2,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(1)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(8)*delta(2)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(3)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(7)*delta(1)*Fields(2,2,x0,y0,z0)**2 &
      +CoeffEpsDisp(9)*delta(2)*Fields(2,2,x0,y0,z0)**2 &
      +CoeffEpsDisp(7)*delta(3)*Fields(2,2,x0,y0,z0)**2 &
      +CoeffEpsDisp(1)*delta(5)*Fields(1,1,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(5)*Fields(1,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(1)*delta(4)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(4)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(4)*delta(1)*Fields(3,1,x0,y0,z0)**2 &
      +CoeffEpsDisp(4)*delta(2)*Fields(3,1,x0,y0,z0)**2 &
      +CoeffEpsDisp(5)*delta(3)*Fields(3,1,x0,y0,z0)**2 &
      +CoeffEpsDisp(2)*delta(5)*Fields(1,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(5)*Fields(1,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(4)*Fields(2,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(4)*Fields(2,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(1)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(2)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(8)*delta(3)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(7)*delta(1)*Fields(3,2,x0,y0,z0)**2 &
      +CoeffEpsDisp(7)*delta(2)*Fields(3,2,x0,y0,z0)**2 &
      +CoeffEpsDisp(9)*delta(3)*Fields(3,2,x0,y0,z0)**2
    
    DeltaH=DeltaHEpsDisp
  
  CASE DEFAULT
      write(*,*) "mode out of range!"
      call abort
  END SELECT
  
End Function GetDeltaHEpsDisp
