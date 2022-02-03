Function GetForcesJijsm(x0, y0, z0, Fields, e0ij) Result(ForcesJijsm)
  
  Implicit none
  Real*8,  Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(in) :: e0ij(3,3)
  Integer, Intent(in) :: x0, y0, z0
  Real*8              :: eij(3,3), euij(3,3)
  Real*8              :: ForcesJijsm(Max(FieldDim, 6), NumField+1)
  
  Integer             :: x1,y1,z1,xi1,yi1,zi1
  
    
  ForcesJijsm = 0.0D0
    euij = StrainFromu(x0, y0, z0, Fields)
    eij = e0ij + euij
    
  x1 = (x0+1)-floor(real(x0+1-1)/real(cgrid%n1))*cgrid%n1
  y1 = (y0+1)-floor(real(y0+1-1)/real(cgrid%n2))*cgrid%n2
  z1 = (z0+1)-floor(real(z0+1-1)/real(cgrid%n3))*cgrid%n3
  xi1 = (x0-1)-floor(real(x0-1-1)/real(cgrid%n1))*cgrid%n1
  yi1 = (y0-1)-floor(real(y0-1-1)/real(cgrid%n2))*cgrid%n2
  zi1 = (z0-1)-floor(real(z0-1-1)/real(cgrid%n3))*cgrid%n3
  
  ForcesJijsm(1,1) = &
    ForcesJijsm(1,1) &
    -2.*CoeffJijsm(2)*Fields(1,1,x0,y0,z1) &
    -2.*CoeffJijsm(2)*Fields(1,1,x0,y0,zi1) &
    -2.*CoeffJijsm(2)*Fields(1,1,x0,y1,z0) &
    -2.*CoeffJijsm(4)*Fields(1,1,x0,y1,z1) &
    -2.*CoeffJijsm(4)*Fields(1,1,x0,y1,zi1) &
    -2.*CoeffJijsm(2)*Fields(1,1,x0,yi1,z0) &
    -2.*CoeffJijsm(4)*Fields(1,1,x0,yi1,z1) &
    -2.*CoeffJijsm(4)*Fields(1,1,x0,yi1,zi1) &
    -2.*CoeffJijsm(1)*Fields(1,1,x1,y0,z0) &
    -2.*CoeffJijsm(3)*Fields(1,1,x1,y0,z1) &
    -2.*CoeffJijsm(3)*Fields(1,1,x1,y0,zi1) &
    -2.*CoeffJijsm(3)*Fields(1,1,x1,y1,z0) &
    -2.*CoeffJijsm(5)*Fields(1,1,x1,y1,z1) &
    -2.*CoeffJijsm(5)*Fields(1,1,x1,y1,zi1) &
    -2.*CoeffJijsm(3)*Fields(1,1,x1,yi1,z0) &
    -2.*CoeffJijsm(5)*Fields(1,1,x1,yi1,z1) &
    -2.*CoeffJijsm(5)*Fields(1,1,x1,yi1,zi1) &
    -2.*CoeffJijsm(1)*Fields(1,1,xi1,y0,z0) &
    -2.*CoeffJijsm(3)*Fields(1,1,xi1,y0,z1) &
    -2.*CoeffJijsm(3)*Fields(1,1,xi1,y0,zi1) &
    -2.*CoeffJijsm(3)*Fields(1,1,xi1,y1,z0) &
    -2.*CoeffJijsm(5)*Fields(1,1,xi1,y1,z1) &
    -2.*CoeffJijsm(5)*Fields(1,1,xi1,y1,zi1) &
    -2.*CoeffJijsm(3)*Fields(1,1,xi1,yi1,z0) &
    -2.*CoeffJijsm(5)*Fields(1,1,xi1,yi1,z1) &
    -2.*CoeffJijsm(5)*Fields(1,1,xi1,yi1,zi1) &
    -2.*CoeffJijsm(6)*Fields(2,1,x1,y1,z0) &
    -2.*CoeffJijsm(7)*Fields(2,1,x1,y1,z1) &
    -2.*CoeffJijsm(7)*Fields(2,1,x1,y1,zi1) &
    +2.*CoeffJijsm(6)*Fields(2,1,x1,yi1,z0) &
    +2.*CoeffJijsm(7)*Fields(2,1,x1,yi1,z1) &
    +2.*CoeffJijsm(7)*Fields(2,1,x1,yi1,zi1) &
    +2.*CoeffJijsm(6)*Fields(2,1,xi1,y1,z0) &
    +2.*CoeffJijsm(7)*Fields(2,1,xi1,y1,z1) &
    +2.*CoeffJijsm(7)*Fields(2,1,xi1,y1,zi1) &
    -2.*CoeffJijsm(6)*Fields(2,1,xi1,yi1,z0) &
    -2.*CoeffJijsm(7)*Fields(2,1,xi1,yi1,z1) &
    -2.*CoeffJijsm(7)*Fields(2,1,xi1,yi1,zi1) &
    -2.*CoeffJijsm(6)*Fields(3,1,x1,y0,z1) &
    +2.*CoeffJijsm(6)*Fields(3,1,x1,y0,zi1) &
    -2.*CoeffJijsm(7)*Fields(3,1,x1,y1,z1) &
    +2.*CoeffJijsm(7)*Fields(3,1,x1,y1,zi1) &
    -2.*CoeffJijsm(7)*Fields(3,1,x1,yi1,z1) &
    +2.*CoeffJijsm(7)*Fields(3,1,x1,yi1,zi1) &
    +2.*CoeffJijsm(6)*Fields(3,1,xi1,y0,z1) &
    -2.*CoeffJijsm(6)*Fields(3,1,xi1,y0,zi1) &
    +2.*CoeffJijsm(7)*Fields(3,1,xi1,y1,z1) &
    -2.*CoeffJijsm(7)*Fields(3,1,xi1,y1,zi1) &
    +2.*CoeffJijsm(7)*Fields(3,1,xi1,yi1,z1) &
    -2.*CoeffJijsm(7)*Fields(3,1,xi1,yi1,zi1)
  
  ForcesJijsm(2,1) = &
    ForcesJijsm(2,1) &
    -2.*CoeffJijsm(6)*Fields(1,1,x1,y1,z0) &
    -2.*CoeffJijsm(7)*Fields(1,1,x1,y1,z1) &
    -2.*CoeffJijsm(7)*Fields(1,1,x1,y1,zi1) &
    +2.*CoeffJijsm(6)*Fields(1,1,x1,yi1,z0) &
    +2.*CoeffJijsm(7)*Fields(1,1,x1,yi1,z1) &
    +2.*CoeffJijsm(7)*Fields(1,1,x1,yi1,zi1) &
    +2.*CoeffJijsm(6)*Fields(1,1,xi1,y1,z0) &
    +2.*CoeffJijsm(7)*Fields(1,1,xi1,y1,z1) &
    +2.*CoeffJijsm(7)*Fields(1,1,xi1,y1,zi1) &
    -2.*CoeffJijsm(6)*Fields(1,1,xi1,yi1,z0) &
    -2.*CoeffJijsm(7)*Fields(1,1,xi1,yi1,z1) &
    -2.*CoeffJijsm(7)*Fields(1,1,xi1,yi1,zi1) &
    -2.*CoeffJijsm(2)*Fields(2,1,x0,y0,z1) &
    -2.*CoeffJijsm(2)*Fields(2,1,x0,y0,zi1) &
    -2.*CoeffJijsm(1)*Fields(2,1,x0,y1,z0) &
    -2.*CoeffJijsm(3)*Fields(2,1,x0,y1,z1) &
    -2.*CoeffJijsm(3)*Fields(2,1,x0,y1,zi1) &
    -2.*CoeffJijsm(1)*Fields(2,1,x0,yi1,z0) &
    -2.*CoeffJijsm(3)*Fields(2,1,x0,yi1,z1) &
    -2.*CoeffJijsm(3)*Fields(2,1,x0,yi1,zi1) &
    -2.*CoeffJijsm(2)*Fields(2,1,x1,y0,z0) &
    -2.*CoeffJijsm(4)*Fields(2,1,x1,y0,z1) &
    -2.*CoeffJijsm(4)*Fields(2,1,x1,y0,zi1) &
    -2.*CoeffJijsm(3)*Fields(2,1,x1,y1,z0) &
    -2.*CoeffJijsm(5)*Fields(2,1,x1,y1,z1) &
    -2.*CoeffJijsm(5)*Fields(2,1,x1,y1,zi1) &
    -2.*CoeffJijsm(3)*Fields(2,1,x1,yi1,z0) &
    -2.*CoeffJijsm(5)*Fields(2,1,x1,yi1,z1) &
    -2.*CoeffJijsm(5)*Fields(2,1,x1,yi1,zi1) &
    -2.*CoeffJijsm(2)*Fields(2,1,xi1,y0,z0) &
    -2.*CoeffJijsm(4)*Fields(2,1,xi1,y0,z1) &
    -2.*CoeffJijsm(4)*Fields(2,1,xi1,y0,zi1) &
    -2.*CoeffJijsm(3)*Fields(2,1,xi1,y1,z0) &
    -2.*CoeffJijsm(5)*Fields(2,1,xi1,y1,z1) &
    -2.*CoeffJijsm(5)*Fields(2,1,xi1,y1,zi1) &
    -2.*CoeffJijsm(3)*Fields(2,1,xi1,yi1,z0) &
    -2.*CoeffJijsm(5)*Fields(2,1,xi1,yi1,z1) &
    -2.*CoeffJijsm(5)*Fields(2,1,xi1,yi1,zi1) &
    -2.*CoeffJijsm(6)*Fields(3,1,x0,y1,z1) &
    +2.*CoeffJijsm(6)*Fields(3,1,x0,y1,zi1) &
    +2.*CoeffJijsm(6)*Fields(3,1,x0,yi1,z1) &
    -2.*CoeffJijsm(6)*Fields(3,1,x0,yi1,zi1) &
    -2.*CoeffJijsm(7)*Fields(3,1,x1,y1,z1) &
    +2.*CoeffJijsm(7)*Fields(3,1,x1,y1,zi1) &
    +2.*CoeffJijsm(7)*Fields(3,1,x1,yi1,z1) &
    -2.*CoeffJijsm(7)*Fields(3,1,x1,yi1,zi1) &
    -2.*CoeffJijsm(7)*Fields(3,1,xi1,y1,z1) &
    +2.*CoeffJijsm(7)*Fields(3,1,xi1,y1,zi1) &
    +2.*CoeffJijsm(7)*Fields(3,1,xi1,yi1,z1) &
    -2.*CoeffJijsm(7)*Fields(3,1,xi1,yi1,zi1)
  
  ForcesJijsm(3,1) = &
    ForcesJijsm(3,1) &
    -2.*CoeffJijsm(6)*Fields(1,1,x1,y0,z1) &
    +2.*CoeffJijsm(6)*Fields(1,1,x1,y0,zi1) &
    -2.*CoeffJijsm(7)*Fields(1,1,x1,y1,z1) &
    +2.*CoeffJijsm(7)*Fields(1,1,x1,y1,zi1) &
    -2.*CoeffJijsm(7)*Fields(1,1,x1,yi1,z1) &
    +2.*CoeffJijsm(7)*Fields(1,1,x1,yi1,zi1) &
    +2.*CoeffJijsm(6)*Fields(1,1,xi1,y0,z1) &
    -2.*CoeffJijsm(6)*Fields(1,1,xi1,y0,zi1) &
    +2.*CoeffJijsm(7)*Fields(1,1,xi1,y1,z1) &
    -2.*CoeffJijsm(7)*Fields(1,1,xi1,y1,zi1) &
    +2.*CoeffJijsm(7)*Fields(1,1,xi1,yi1,z1) &
    -2.*CoeffJijsm(7)*Fields(1,1,xi1,yi1,zi1) &
    -2.*CoeffJijsm(6)*Fields(2,1,x0,y1,z1) &
    +2.*CoeffJijsm(6)*Fields(2,1,x0,y1,zi1) &
    +2.*CoeffJijsm(6)*Fields(2,1,x0,yi1,z1) &
    -2.*CoeffJijsm(6)*Fields(2,1,x0,yi1,zi1) &
    -2.*CoeffJijsm(7)*Fields(2,1,x1,y1,z1) &
    +2.*CoeffJijsm(7)*Fields(2,1,x1,y1,zi1) &
    +2.*CoeffJijsm(7)*Fields(2,1,x1,yi1,z1) &
    -2.*CoeffJijsm(7)*Fields(2,1,x1,yi1,zi1) &
    -2.*CoeffJijsm(7)*Fields(2,1,xi1,y1,z1) &
    +2.*CoeffJijsm(7)*Fields(2,1,xi1,y1,zi1) &
    +2.*CoeffJijsm(7)*Fields(2,1,xi1,yi1,z1) &
    -2.*CoeffJijsm(7)*Fields(2,1,xi1,yi1,zi1) &
    -2.*CoeffJijsm(1)*Fields(3,1,x0,y0,z1) &
    -2.*CoeffJijsm(1)*Fields(3,1,x0,y0,zi1) &
    -2.*CoeffJijsm(2)*Fields(3,1,x0,y1,z0) &
    -2.*CoeffJijsm(3)*Fields(3,1,x0,y1,z1) &
    -2.*CoeffJijsm(3)*Fields(3,1,x0,y1,zi1) &
    -2.*CoeffJijsm(2)*Fields(3,1,x0,yi1,z0) &
    -2.*CoeffJijsm(3)*Fields(3,1,x0,yi1,z1) &
    -2.*CoeffJijsm(3)*Fields(3,1,x0,yi1,zi1) &
    -2.*CoeffJijsm(2)*Fields(3,1,x1,y0,z0) &
    -2.*CoeffJijsm(3)*Fields(3,1,x1,y0,z1) &
    -2.*CoeffJijsm(3)*Fields(3,1,x1,y0,zi1) &
    -2.*CoeffJijsm(4)*Fields(3,1,x1,y1,z0) &
    -2.*CoeffJijsm(5)*Fields(3,1,x1,y1,z1) &
    -2.*CoeffJijsm(5)*Fields(3,1,x1,y1,zi1) &
    -2.*CoeffJijsm(4)*Fields(3,1,x1,yi1,z0) &
    -2.*CoeffJijsm(5)*Fields(3,1,x1,yi1,z1) &
    -2.*CoeffJijsm(5)*Fields(3,1,x1,yi1,zi1) &
    -2.*CoeffJijsm(2)*Fields(3,1,xi1,y0,z0) &
    -2.*CoeffJijsm(3)*Fields(3,1,xi1,y0,z1) &
    -2.*CoeffJijsm(3)*Fields(3,1,xi1,y0,zi1) &
    -2.*CoeffJijsm(4)*Fields(3,1,xi1,y1,z0) &
    -2.*CoeffJijsm(5)*Fields(3,1,xi1,y1,z1) &
    -2.*CoeffJijsm(5)*Fields(3,1,xi1,y1,zi1) &
    -2.*CoeffJijsm(4)*Fields(3,1,xi1,yi1,z0) &
    -2.*CoeffJijsm(5)*Fields(3,1,xi1,yi1,z1) &
    -2.*CoeffJijsm(5)*Fields(3,1,xi1,yi1,zi1)
  
  ForcesJijsm(1,2) = &
    ForcesJijsm(1,2) &
    +0.
  
  ForcesJijsm(2,2) = &
    ForcesJijsm(2,2) &
    +0.
  
  ForcesJijsm(3,2) = &
    ForcesJijsm(3,2) &
    +0.
  
  ForcesJijsm(1,3) = &
    ForcesJijsm(1,3) &
    +0.
  
  ForcesJijsm(2,3) = &
    ForcesJijsm(2,3) &
    +0.
  
  ForcesJijsm(3,3) = &
    ForcesJijsm(3,3) &
    +0.
  
  ForcesJijsm(1,4) = &
    ForcesJijsm(1,4) &
    +0.
  
  ForcesJijsm(2,4) = &
    ForcesJijsm(2,4) &
    +0.
  
  ForcesJijsm(3,4) = &
    ForcesJijsm(3,4) &
    +0.
  
  ForcesJijsm(4,4) = &
    ForcesJijsm(4,4) &
    +0.
  
  ForcesJijsm(5,4) = &
    ForcesJijsm(5,4) &
    +0.
  
  ForcesJijsm(6,4) = &
    ForcesJijsm(6,4) &
    +0.
  
  
End Function GetForcesJijsm
