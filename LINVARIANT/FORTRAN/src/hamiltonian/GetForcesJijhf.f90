Function GetForcesJijhf(x0, y0, z0, Fields, e0ij) Result(ForcesJijhf)
  
  Implicit none
  Real*8,  Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(in) :: e0ij(3,3)
  Integer, Intent(in) :: x0, y0, z0
  Real*8              :: eij(3,3), euij(3,3)
  Real*8              :: ForcesJijhf(Max(FieldDim, 6), NumField+1)
  
  Integer             :: x1,y1,z1,xi1,yi1,zi1
  
    
  ForcesJijhf = 0.0D0
    euij = StrainFromu(x0, y0, z0, Fields)
    eij = e0ij + euij
    
  x1 = (x0+1)-floor(real(x0+1-1)/real(cgrid%n1))*cgrid%n1
  y1 = (y0+1)-floor(real(y0+1-1)/real(cgrid%n2))*cgrid%n2
  z1 = (z0+1)-floor(real(z0+1-1)/real(cgrid%n3))*cgrid%n3
  xi1 = (x0-1)-floor(real(x0-1-1)/real(cgrid%n1))*cgrid%n1
  yi1 = (y0-1)-floor(real(y0-1-1)/real(cgrid%n2))*cgrid%n2
  zi1 = (z0-1)-floor(real(z0-1-1)/real(cgrid%n3))*cgrid%n3
  
  ForcesJijhf(1,1) = &
    ForcesJijhf(1,1) &
    +0.
  
  ForcesJijhf(2,1) = &
    ForcesJijhf(2,1) &
    +0.
  
  ForcesJijhf(3,1) = &
    ForcesJijhf(3,1) &
    +0.
  
  ForcesJijhf(1,2) = &
    ForcesJijhf(1,2) &
    -2.*CoeffJijhf(2)*Fields(1,2,x0,y0,z1) &
    -2.*CoeffJijhf(2)*Fields(1,2,x0,y0,zi1) &
    -2.*CoeffJijhf(2)*Fields(1,2,x0,y1,z0) &
    -2.*CoeffJijhf(4)*Fields(1,2,x0,y1,z1) &
    -2.*CoeffJijhf(4)*Fields(1,2,x0,y1,zi1) &
    -2.*CoeffJijhf(2)*Fields(1,2,x0,yi1,z0) &
    -2.*CoeffJijhf(4)*Fields(1,2,x0,yi1,z1) &
    -2.*CoeffJijhf(4)*Fields(1,2,x0,yi1,zi1) &
    -2.*CoeffJijhf(1)*Fields(1,2,x1,y0,z0) &
    -2.*CoeffJijhf(3)*Fields(1,2,x1,y0,z1) &
    -2.*CoeffJijhf(3)*Fields(1,2,x1,y0,zi1) &
    -2.*CoeffJijhf(3)*Fields(1,2,x1,y1,z0) &
    -2.*CoeffJijhf(5)*Fields(1,2,x1,y1,z1) &
    -2.*CoeffJijhf(5)*Fields(1,2,x1,y1,zi1) &
    -2.*CoeffJijhf(3)*Fields(1,2,x1,yi1,z0) &
    -2.*CoeffJijhf(5)*Fields(1,2,x1,yi1,z1) &
    -2.*CoeffJijhf(5)*Fields(1,2,x1,yi1,zi1) &
    -2.*CoeffJijhf(1)*Fields(1,2,xi1,y0,z0) &
    -2.*CoeffJijhf(3)*Fields(1,2,xi1,y0,z1) &
    -2.*CoeffJijhf(3)*Fields(1,2,xi1,y0,zi1) &
    -2.*CoeffJijhf(3)*Fields(1,2,xi1,y1,z0) &
    -2.*CoeffJijhf(5)*Fields(1,2,xi1,y1,z1) &
    -2.*CoeffJijhf(5)*Fields(1,2,xi1,y1,zi1) &
    -2.*CoeffJijhf(3)*Fields(1,2,xi1,yi1,z0) &
    -2.*CoeffJijhf(5)*Fields(1,2,xi1,yi1,z1) &
    -2.*CoeffJijhf(5)*Fields(1,2,xi1,yi1,zi1) &
    -2.*CoeffJijhf(6)*Fields(2,2,x1,y1,z0) &
    -2.*CoeffJijhf(7)*Fields(2,2,x1,y1,z1) &
    -2.*CoeffJijhf(7)*Fields(2,2,x1,y1,zi1) &
    +2.*CoeffJijhf(6)*Fields(2,2,x1,yi1,z0) &
    +2.*CoeffJijhf(7)*Fields(2,2,x1,yi1,z1) &
    +2.*CoeffJijhf(7)*Fields(2,2,x1,yi1,zi1) &
    +2.*CoeffJijhf(6)*Fields(2,2,xi1,y1,z0) &
    +2.*CoeffJijhf(7)*Fields(2,2,xi1,y1,z1) &
    +2.*CoeffJijhf(7)*Fields(2,2,xi1,y1,zi1) &
    -2.*CoeffJijhf(6)*Fields(2,2,xi1,yi1,z0) &
    -2.*CoeffJijhf(7)*Fields(2,2,xi1,yi1,z1) &
    -2.*CoeffJijhf(7)*Fields(2,2,xi1,yi1,zi1) &
    -2.*CoeffJijhf(6)*Fields(3,2,x1,y0,z1) &
    +2.*CoeffJijhf(6)*Fields(3,2,x1,y0,zi1) &
    -2.*CoeffJijhf(7)*Fields(3,2,x1,y1,z1) &
    +2.*CoeffJijhf(7)*Fields(3,2,x1,y1,zi1) &
    -2.*CoeffJijhf(7)*Fields(3,2,x1,yi1,z1) &
    +2.*CoeffJijhf(7)*Fields(3,2,x1,yi1,zi1) &
    +2.*CoeffJijhf(6)*Fields(3,2,xi1,y0,z1) &
    -2.*CoeffJijhf(6)*Fields(3,2,xi1,y0,zi1) &
    +2.*CoeffJijhf(7)*Fields(3,2,xi1,y1,z1) &
    -2.*CoeffJijhf(7)*Fields(3,2,xi1,y1,zi1) &
    +2.*CoeffJijhf(7)*Fields(3,2,xi1,yi1,z1) &
    -2.*CoeffJijhf(7)*Fields(3,2,xi1,yi1,zi1)
  
  ForcesJijhf(2,2) = &
    ForcesJijhf(2,2) &
    -2.*CoeffJijhf(6)*Fields(1,2,x1,y1,z0) &
    -2.*CoeffJijhf(7)*Fields(1,2,x1,y1,z1) &
    -2.*CoeffJijhf(7)*Fields(1,2,x1,y1,zi1) &
    +2.*CoeffJijhf(6)*Fields(1,2,x1,yi1,z0) &
    +2.*CoeffJijhf(7)*Fields(1,2,x1,yi1,z1) &
    +2.*CoeffJijhf(7)*Fields(1,2,x1,yi1,zi1) &
    +2.*CoeffJijhf(6)*Fields(1,2,xi1,y1,z0) &
    +2.*CoeffJijhf(7)*Fields(1,2,xi1,y1,z1) &
    +2.*CoeffJijhf(7)*Fields(1,2,xi1,y1,zi1) &
    -2.*CoeffJijhf(6)*Fields(1,2,xi1,yi1,z0) &
    -2.*CoeffJijhf(7)*Fields(1,2,xi1,yi1,z1) &
    -2.*CoeffJijhf(7)*Fields(1,2,xi1,yi1,zi1) &
    -2.*CoeffJijhf(2)*Fields(2,2,x0,y0,z1) &
    -2.*CoeffJijhf(2)*Fields(2,2,x0,y0,zi1) &
    -2.*CoeffJijhf(1)*Fields(2,2,x0,y1,z0) &
    -2.*CoeffJijhf(3)*Fields(2,2,x0,y1,z1) &
    -2.*CoeffJijhf(3)*Fields(2,2,x0,y1,zi1) &
    -2.*CoeffJijhf(1)*Fields(2,2,x0,yi1,z0) &
    -2.*CoeffJijhf(3)*Fields(2,2,x0,yi1,z1) &
    -2.*CoeffJijhf(3)*Fields(2,2,x0,yi1,zi1) &
    -2.*CoeffJijhf(2)*Fields(2,2,x1,y0,z0) &
    -2.*CoeffJijhf(4)*Fields(2,2,x1,y0,z1) &
    -2.*CoeffJijhf(4)*Fields(2,2,x1,y0,zi1) &
    -2.*CoeffJijhf(3)*Fields(2,2,x1,y1,z0) &
    -2.*CoeffJijhf(5)*Fields(2,2,x1,y1,z1) &
    -2.*CoeffJijhf(5)*Fields(2,2,x1,y1,zi1) &
    -2.*CoeffJijhf(3)*Fields(2,2,x1,yi1,z0) &
    -2.*CoeffJijhf(5)*Fields(2,2,x1,yi1,z1) &
    -2.*CoeffJijhf(5)*Fields(2,2,x1,yi1,zi1) &
    -2.*CoeffJijhf(2)*Fields(2,2,xi1,y0,z0) &
    -2.*CoeffJijhf(4)*Fields(2,2,xi1,y0,z1) &
    -2.*CoeffJijhf(4)*Fields(2,2,xi1,y0,zi1) &
    -2.*CoeffJijhf(3)*Fields(2,2,xi1,y1,z0) &
    -2.*CoeffJijhf(5)*Fields(2,2,xi1,y1,z1) &
    -2.*CoeffJijhf(5)*Fields(2,2,xi1,y1,zi1) &
    -2.*CoeffJijhf(3)*Fields(2,2,xi1,yi1,z0) &
    -2.*CoeffJijhf(5)*Fields(2,2,xi1,yi1,z1) &
    -2.*CoeffJijhf(5)*Fields(2,2,xi1,yi1,zi1) &
    -2.*CoeffJijhf(6)*Fields(3,2,x0,y1,z1) &
    +2.*CoeffJijhf(6)*Fields(3,2,x0,y1,zi1) &
    +2.*CoeffJijhf(6)*Fields(3,2,x0,yi1,z1) &
    -2.*CoeffJijhf(6)*Fields(3,2,x0,yi1,zi1) &
    -2.*CoeffJijhf(7)*Fields(3,2,x1,y1,z1) &
    +2.*CoeffJijhf(7)*Fields(3,2,x1,y1,zi1) &
    +2.*CoeffJijhf(7)*Fields(3,2,x1,yi1,z1) &
    -2.*CoeffJijhf(7)*Fields(3,2,x1,yi1,zi1) &
    -2.*CoeffJijhf(7)*Fields(3,2,xi1,y1,z1) &
    +2.*CoeffJijhf(7)*Fields(3,2,xi1,y1,zi1) &
    +2.*CoeffJijhf(7)*Fields(3,2,xi1,yi1,z1) &
    -2.*CoeffJijhf(7)*Fields(3,2,xi1,yi1,zi1)
  
  ForcesJijhf(3,2) = &
    ForcesJijhf(3,2) &
    -2.*CoeffJijhf(6)*Fields(1,2,x1,y0,z1) &
    +2.*CoeffJijhf(6)*Fields(1,2,x1,y0,zi1) &
    -2.*CoeffJijhf(7)*Fields(1,2,x1,y1,z1) &
    +2.*CoeffJijhf(7)*Fields(1,2,x1,y1,zi1) &
    -2.*CoeffJijhf(7)*Fields(1,2,x1,yi1,z1) &
    +2.*CoeffJijhf(7)*Fields(1,2,x1,yi1,zi1) &
    +2.*CoeffJijhf(6)*Fields(1,2,xi1,y0,z1) &
    -2.*CoeffJijhf(6)*Fields(1,2,xi1,y0,zi1) &
    +2.*CoeffJijhf(7)*Fields(1,2,xi1,y1,z1) &
    -2.*CoeffJijhf(7)*Fields(1,2,xi1,y1,zi1) &
    +2.*CoeffJijhf(7)*Fields(1,2,xi1,yi1,z1) &
    -2.*CoeffJijhf(7)*Fields(1,2,xi1,yi1,zi1) &
    -2.*CoeffJijhf(6)*Fields(2,2,x0,y1,z1) &
    +2.*CoeffJijhf(6)*Fields(2,2,x0,y1,zi1) &
    +2.*CoeffJijhf(6)*Fields(2,2,x0,yi1,z1) &
    -2.*CoeffJijhf(6)*Fields(2,2,x0,yi1,zi1) &
    -2.*CoeffJijhf(7)*Fields(2,2,x1,y1,z1) &
    +2.*CoeffJijhf(7)*Fields(2,2,x1,y1,zi1) &
    +2.*CoeffJijhf(7)*Fields(2,2,x1,yi1,z1) &
    -2.*CoeffJijhf(7)*Fields(2,2,x1,yi1,zi1) &
    -2.*CoeffJijhf(7)*Fields(2,2,xi1,y1,z1) &
    +2.*CoeffJijhf(7)*Fields(2,2,xi1,y1,zi1) &
    +2.*CoeffJijhf(7)*Fields(2,2,xi1,yi1,z1) &
    -2.*CoeffJijhf(7)*Fields(2,2,xi1,yi1,zi1) &
    -2.*CoeffJijhf(1)*Fields(3,2,x0,y0,z1) &
    -2.*CoeffJijhf(1)*Fields(3,2,x0,y0,zi1) &
    -2.*CoeffJijhf(2)*Fields(3,2,x0,y1,z0) &
    -2.*CoeffJijhf(3)*Fields(3,2,x0,y1,z1) &
    -2.*CoeffJijhf(3)*Fields(3,2,x0,y1,zi1) &
    -2.*CoeffJijhf(2)*Fields(3,2,x0,yi1,z0) &
    -2.*CoeffJijhf(3)*Fields(3,2,x0,yi1,z1) &
    -2.*CoeffJijhf(3)*Fields(3,2,x0,yi1,zi1) &
    -2.*CoeffJijhf(2)*Fields(3,2,x1,y0,z0) &
    -2.*CoeffJijhf(3)*Fields(3,2,x1,y0,z1) &
    -2.*CoeffJijhf(3)*Fields(3,2,x1,y0,zi1) &
    -2.*CoeffJijhf(4)*Fields(3,2,x1,y1,z0) &
    -2.*CoeffJijhf(5)*Fields(3,2,x1,y1,z1) &
    -2.*CoeffJijhf(5)*Fields(3,2,x1,y1,zi1) &
    -2.*CoeffJijhf(4)*Fields(3,2,x1,yi1,z0) &
    -2.*CoeffJijhf(5)*Fields(3,2,x1,yi1,z1) &
    -2.*CoeffJijhf(5)*Fields(3,2,x1,yi1,zi1) &
    -2.*CoeffJijhf(2)*Fields(3,2,xi1,y0,z0) &
    -2.*CoeffJijhf(3)*Fields(3,2,xi1,y0,z1) &
    -2.*CoeffJijhf(3)*Fields(3,2,xi1,y0,zi1) &
    -2.*CoeffJijhf(4)*Fields(3,2,xi1,y1,z0) &
    -2.*CoeffJijhf(5)*Fields(3,2,xi1,y1,z1) &
    -2.*CoeffJijhf(5)*Fields(3,2,xi1,y1,zi1) &
    -2.*CoeffJijhf(4)*Fields(3,2,xi1,yi1,z0) &
    -2.*CoeffJijhf(5)*Fields(3,2,xi1,yi1,z1) &
    -2.*CoeffJijhf(5)*Fields(3,2,xi1,yi1,zi1)
  
  ForcesJijhf(1,3) = &
    ForcesJijhf(1,3) &
    +0.
  
  ForcesJijhf(2,3) = &
    ForcesJijhf(2,3) &
    +0.
  
  ForcesJijhf(3,3) = &
    ForcesJijhf(3,3) &
    +0.
  
  ForcesJijhf(1,4) = &
    ForcesJijhf(1,4) &
    +0.
  
  ForcesJijhf(2,4) = &
    ForcesJijhf(2,4) &
    +0.
  
  ForcesJijhf(3,4) = &
    ForcesJijhf(3,4) &
    +0.
  
  ForcesJijhf(4,4) = &
    ForcesJijhf(4,4) &
    +0.
  
  ForcesJijhf(5,4) = &
    ForcesJijhf(5,4) &
    +0.
  
  ForcesJijhf(6,4) = &
    ForcesJijhf(6,4) &
    +0.
  
  
End Function GetForcesJijhf
