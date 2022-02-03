Function GetForcesDispu(x0, y0, z0, Fields, e0ij) Result(ForcesDispu)
  
  Implicit none
  Real*8,  Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(in) :: e0ij(3,3)
  Integer, Intent(in) :: x0, y0, z0
  Real*8              :: eij(3,3), euij(3,3)
  Real*8              :: ForcesDispu(Max(FieldDim, 6), NumField+1)
  
  Integer             :: x1,y1,z1,xi1,yi1,zi1
  
    
  ForcesDispu = 0.0D0
    euij = StrainFromu(x0, y0, z0, Fields)
    eij = e0ij + euij
    
  x1 = (x0+1)-floor(real(x0+1-1)/real(cgrid%n1))*cgrid%n1
  y1 = (y0+1)-floor(real(y0+1-1)/real(cgrid%n2))*cgrid%n2
  z1 = (z0+1)-floor(real(z0+1-1)/real(cgrid%n3))*cgrid%n3
  xi1 = (x0-1)-floor(real(x0-1-1)/real(cgrid%n1))*cgrid%n1
  yi1 = (y0-1)-floor(real(y0-1-1)/real(cgrid%n2))*cgrid%n2
  zi1 = (z0-1)-floor(real(z0-1-1)/real(cgrid%n3))*cgrid%n3
  
  ForcesDispu(1,1) = &
    ForcesDispu(1,1) &
    +0.
  
  ForcesDispu(2,1) = &
    ForcesDispu(2,1) &
    +0.
  
  ForcesDispu(3,1) = &
    ForcesDispu(3,1) &
    +0.
  
  ForcesDispu(1,2) = &
    ForcesDispu(1,2) &
    +0.
  
  ForcesDispu(2,2) = &
    ForcesDispu(2,2) &
    +0.
  
  ForcesDispu(3,2) = &
    ForcesDispu(3,2) &
    +0.
  
  ForcesDispu(1,3) = &
    ForcesDispu(1,3) &
    +0.25*CoeffDispu(5)*Fields(1,1,x0,y0,z0)**2 &
    +0.25*CoeffDispu(5)*Fields(1,1,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(5)*Fields(1,1,x0,yi1,z0)**2 &
    +0.25*CoeffDispu(5)*Fields(1,1,x0,yi1,zi1)**2 &
    -0.25*CoeffDispu(5)*Fields(1,1,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(5)*Fields(1,1,xi1,y0,zi1)**2 &
    -0.25*CoeffDispu(5)*Fields(1,1,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(5)*Fields(1,1,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(8)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
    +0.25*CoeffDispu(9)*Fields(1,2,x0,y0,z0)**2 &
    +0.25*CoeffDispu(8)*Fields(1,1,x0,y0,zi1)*Fields(1,2,x0,y0,zi1) &
    +0.25*CoeffDispu(9)*Fields(1,2,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(8)*Fields(1,1,x0,yi1,z0)*Fields(1,2,x0,yi1,z0) &
    +0.25*CoeffDispu(9)*Fields(1,2,x0,yi1,z0)**2 &
    +0.25*CoeffDispu(8)*Fields(1,1,x0,yi1,zi1)*Fields(1,2,x0,yi1,zi1) &
    +0.25*CoeffDispu(9)*Fields(1,2,x0,yi1,zi1)**2 &
    -0.25*CoeffDispu(8)*Fields(1,1,xi1,y0,z0)*Fields(1,2,xi1,y0,z0) &
    -0.25*CoeffDispu(9)*Fields(1,2,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(8)*Fields(1,1,xi1,y0,zi1)*Fields(1,2,xi1,y0,zi1) &
    -0.25*CoeffDispu(9)*Fields(1,2,xi1,y0,zi1)**2 &
    -0.25*CoeffDispu(8)*Fields(1,1,xi1,yi1,z0)*Fields(1,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(9)*Fields(1,2,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(8)*Fields(1,1,xi1,yi1,zi1)*Fields(1,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(9)*Fields(1,2,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,x0,y0,z0)*Fields(2,1,x0,y0,z0) &
    +0.25*CoeffDispu(2)*Fields(1,2,x0,y0,z0)*Fields(2,1,x0,y0,z0) &
    +0.25*CoeffDispu(4)*Fields(2,1,x0,y0,z0)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,x0,y0,zi1)*Fields(2,1,x0,y0,zi1) &
    +0.25*CoeffDispu(2)*Fields(1,2,x0,y0,zi1)*Fields(2,1,x0,y0,zi1) &
    +0.25*CoeffDispu(4)*Fields(2,1,x0,y0,zi1)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,x0,yi1,z0)*Fields(2,1,x0,yi1,z0) &
    -0.25*CoeffDispu(2)*Fields(1,2,x0,yi1,z0)*Fields(2,1,x0,yi1,z0) &
    +0.25*CoeffDispu(4)*Fields(2,1,x0,yi1,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,x0,yi1,zi1)*Fields(2,1,x0,yi1,zi1) &
    -0.25*CoeffDispu(2)*Fields(1,2,x0,yi1,zi1)*Fields(2,1,x0,yi1,zi1) &
    +0.25*CoeffDispu(4)*Fields(2,1,x0,yi1,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,xi1,y0,z0)*Fields(2,1,xi1,y0,z0) &
    +0.25*CoeffDispu(2)*Fields(1,2,xi1,y0,z0)*Fields(2,1,xi1,y0,z0) &
    -0.25*CoeffDispu(4)*Fields(2,1,xi1,y0,z0)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,xi1,y0,zi1)*Fields(2,1,xi1,y0,zi1) &
    +0.25*CoeffDispu(2)*Fields(1,2,xi1,y0,zi1)*Fields(2,1,xi1,y0,zi1) &
    -0.25*CoeffDispu(4)*Fields(2,1,xi1,y0,zi1)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,xi1,yi1,z0)*Fields(2,1,xi1,yi1,z0) &
    -0.25*CoeffDispu(2)*Fields(1,2,xi1,yi1,z0)*Fields(2,1,xi1,yi1,z0) &
    -0.25*CoeffDispu(4)*Fields(2,1,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,xi1,yi1,zi1)*Fields(2,1,xi1,yi1,zi1) &
    -0.25*CoeffDispu(2)*Fields(1,2,xi1,yi1,zi1)*Fields(2,1,xi1,yi1,zi1) &
    -0.25*CoeffDispu(4)*Fields(2,1,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +0.25*CoeffDispu(3)*Fields(1,2,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +0.25*CoeffDispu(6)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +0.25*CoeffDispu(7)*Fields(2,2,x0,y0,z0)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,x0,y0,zi1)*Fields(2,2,x0,y0,zi1) &
    +0.25*CoeffDispu(3)*Fields(1,2,x0,y0,zi1)*Fields(2,2,x0,y0,zi1) &
    +0.25*CoeffDispu(6)*Fields(2,1,x0,y0,zi1)*Fields(2,2,x0,y0,zi1) &
    +0.25*CoeffDispu(7)*Fields(2,2,x0,y0,zi1)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,x0,yi1,z0)*Fields(2,2,x0,yi1,z0) &
    -0.25*CoeffDispu(3)*Fields(1,2,x0,yi1,z0)*Fields(2,2,x0,yi1,z0) &
    +0.25*CoeffDispu(6)*Fields(2,1,x0,yi1,z0)*Fields(2,2,x0,yi1,z0) &
    +0.25*CoeffDispu(7)*Fields(2,2,x0,yi1,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,x0,yi1,zi1)*Fields(2,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(3)*Fields(1,2,x0,yi1,zi1)*Fields(2,2,x0,yi1,zi1) &
    +0.25*CoeffDispu(6)*Fields(2,1,x0,yi1,zi1)*Fields(2,2,x0,yi1,zi1) &
    +0.25*CoeffDispu(7)*Fields(2,2,x0,yi1,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,xi1,y0,z0)*Fields(2,2,xi1,y0,z0) &
    +0.25*CoeffDispu(3)*Fields(1,2,xi1,y0,z0)*Fields(2,2,xi1,y0,z0) &
    -0.25*CoeffDispu(6)*Fields(2,1,xi1,y0,z0)*Fields(2,2,xi1,y0,z0) &
    -0.25*CoeffDispu(7)*Fields(2,2,xi1,y0,z0)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,xi1,y0,zi1)*Fields(2,2,xi1,y0,zi1) &
    +0.25*CoeffDispu(3)*Fields(1,2,xi1,y0,zi1)*Fields(2,2,xi1,y0,zi1) &
    -0.25*CoeffDispu(6)*Fields(2,1,xi1,y0,zi1)*Fields(2,2,xi1,y0,zi1) &
    -0.25*CoeffDispu(7)*Fields(2,2,xi1,y0,zi1)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,xi1,yi1,z0)*Fields(2,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(3)*Fields(1,2,xi1,yi1,z0)*Fields(2,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(6)*Fields(2,1,xi1,yi1,z0)*Fields(2,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(7)*Fields(2,2,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,xi1,yi1,zi1)*Fields(2,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(3)*Fields(1,2,xi1,yi1,zi1)*Fields(2,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(6)*Fields(2,1,xi1,yi1,zi1)*Fields(2,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(7)*Fields(2,2,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
    +0.25*CoeffDispu(2)*Fields(1,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
    +0.25*CoeffDispu(4)*Fields(3,1,x0,y0,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,x0,y0,zi1)*Fields(3,1,x0,y0,zi1) &
    -0.25*CoeffDispu(2)*Fields(1,2,x0,y0,zi1)*Fields(3,1,x0,y0,zi1) &
    +0.25*CoeffDispu(4)*Fields(3,1,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,x0,yi1,z0)*Fields(3,1,x0,yi1,z0) &
    +0.25*CoeffDispu(2)*Fields(1,2,x0,yi1,z0)*Fields(3,1,x0,yi1,z0) &
    +0.25*CoeffDispu(4)*Fields(3,1,x0,yi1,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,x0,yi1,zi1)*Fields(3,1,x0,yi1,zi1) &
    -0.25*CoeffDispu(2)*Fields(1,2,x0,yi1,zi1)*Fields(3,1,x0,yi1,zi1) &
    +0.25*CoeffDispu(4)*Fields(3,1,x0,yi1,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,xi1,y0,z0)*Fields(3,1,xi1,y0,z0) &
    +0.25*CoeffDispu(2)*Fields(1,2,xi1,y0,z0)*Fields(3,1,xi1,y0,z0) &
    -0.25*CoeffDispu(4)*Fields(3,1,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,xi1,y0,zi1)*Fields(3,1,xi1,y0,zi1) &
    -0.25*CoeffDispu(2)*Fields(1,2,xi1,y0,zi1)*Fields(3,1,xi1,y0,zi1) &
    -0.25*CoeffDispu(4)*Fields(3,1,xi1,y0,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,xi1,yi1,z0)*Fields(3,1,xi1,yi1,z0) &
    +0.25*CoeffDispu(2)*Fields(1,2,xi1,yi1,z0)*Fields(3,1,xi1,yi1,z0) &
    -0.25*CoeffDispu(4)*Fields(3,1,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,xi1,yi1,zi1)*Fields(3,1,xi1,yi1,zi1) &
    -0.25*CoeffDispu(2)*Fields(1,2,xi1,yi1,zi1)*Fields(3,1,xi1,yi1,zi1) &
    -0.25*CoeffDispu(4)*Fields(3,1,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +0.25*CoeffDispu(3)*Fields(1,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +0.25*CoeffDispu(6)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +0.25*CoeffDispu(7)*Fields(3,2,x0,y0,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
    -0.25*CoeffDispu(3)*Fields(1,2,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
    +0.25*CoeffDispu(6)*Fields(3,1,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
    +0.25*CoeffDispu(7)*Fields(3,2,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
    +0.25*CoeffDispu(3)*Fields(1,2,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
    +0.25*CoeffDispu(6)*Fields(3,1,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
    +0.25*CoeffDispu(7)*Fields(3,2,x0,yi1,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(3)*Fields(1,2,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
    +0.25*CoeffDispu(6)*Fields(3,1,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
    +0.25*CoeffDispu(7)*Fields(3,2,x0,yi1,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
    +0.25*CoeffDispu(3)*Fields(1,2,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
    -0.25*CoeffDispu(6)*Fields(3,1,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
    -0.25*CoeffDispu(7)*Fields(3,2,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
    -0.25*CoeffDispu(3)*Fields(1,2,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
    -0.25*CoeffDispu(6)*Fields(3,1,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
    -0.25*CoeffDispu(7)*Fields(3,2,xi1,y0,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
    +0.25*CoeffDispu(3)*Fields(1,2,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(6)*Fields(3,1,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(7)*Fields(3,2,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(3)*Fields(1,2,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(6)*Fields(3,1,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(7)*Fields(3,2,xi1,yi1,zi1)**2
  
  ForcesDispu(2,3) = &
    ForcesDispu(2,3) &
    +0.25*CoeffDispu(4)*Fields(1,1,x0,y0,z0)**2 &
    +0.25*CoeffDispu(4)*Fields(1,1,x0,y0,zi1)**2 &
    -0.25*CoeffDispu(4)*Fields(1,1,x0,yi1,z0)**2 &
    -0.25*CoeffDispu(4)*Fields(1,1,x0,yi1,zi1)**2 &
    +0.25*CoeffDispu(4)*Fields(1,1,xi1,y0,z0)**2 &
    +0.25*CoeffDispu(4)*Fields(1,1,xi1,y0,zi1)**2 &
    -0.25*CoeffDispu(4)*Fields(1,1,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(4)*Fields(1,1,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(6)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
    +0.25*CoeffDispu(7)*Fields(1,2,x0,y0,z0)**2 &
    +0.25*CoeffDispu(6)*Fields(1,1,x0,y0,zi1)*Fields(1,2,x0,y0,zi1) &
    +0.25*CoeffDispu(7)*Fields(1,2,x0,y0,zi1)**2 &
    -0.25*CoeffDispu(6)*Fields(1,1,x0,yi1,z0)*Fields(1,2,x0,yi1,z0) &
    -0.25*CoeffDispu(7)*Fields(1,2,x0,yi1,z0)**2 &
    -0.25*CoeffDispu(6)*Fields(1,1,x0,yi1,zi1)*Fields(1,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(7)*Fields(1,2,x0,yi1,zi1)**2 &
    +0.25*CoeffDispu(6)*Fields(1,1,xi1,y0,z0)*Fields(1,2,xi1,y0,z0) &
    +0.25*CoeffDispu(7)*Fields(1,2,xi1,y0,z0)**2 &
    +0.25*CoeffDispu(6)*Fields(1,1,xi1,y0,zi1)*Fields(1,2,xi1,y0,zi1) &
    +0.25*CoeffDispu(7)*Fields(1,2,xi1,y0,zi1)**2 &
    -0.25*CoeffDispu(6)*Fields(1,1,xi1,yi1,z0)*Fields(1,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(7)*Fields(1,2,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(6)*Fields(1,1,xi1,yi1,zi1)*Fields(1,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(7)*Fields(1,2,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,x0,y0,z0)*Fields(2,1,x0,y0,z0) &
    +0.25*CoeffDispu(2)*Fields(1,2,x0,y0,z0)*Fields(2,1,x0,y0,z0) &
    +0.25*CoeffDispu(5)*Fields(2,1,x0,y0,z0)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,x0,y0,zi1)*Fields(2,1,x0,y0,zi1) &
    +0.25*CoeffDispu(2)*Fields(1,2,x0,y0,zi1)*Fields(2,1,x0,y0,zi1) &
    +0.25*CoeffDispu(5)*Fields(2,1,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,x0,yi1,z0)*Fields(2,1,x0,yi1,z0) &
    +0.25*CoeffDispu(2)*Fields(1,2,x0,yi1,z0)*Fields(2,1,x0,yi1,z0) &
    -0.25*CoeffDispu(5)*Fields(2,1,x0,yi1,z0)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,x0,yi1,zi1)*Fields(2,1,x0,yi1,zi1) &
    +0.25*CoeffDispu(2)*Fields(1,2,x0,yi1,zi1)*Fields(2,1,x0,yi1,zi1) &
    -0.25*CoeffDispu(5)*Fields(2,1,x0,yi1,zi1)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,xi1,y0,z0)*Fields(2,1,xi1,y0,z0) &
    -0.25*CoeffDispu(2)*Fields(1,2,xi1,y0,z0)*Fields(2,1,xi1,y0,z0) &
    +0.25*CoeffDispu(5)*Fields(2,1,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,xi1,y0,zi1)*Fields(2,1,xi1,y0,zi1) &
    -0.25*CoeffDispu(2)*Fields(1,2,xi1,y0,zi1)*Fields(2,1,xi1,y0,zi1) &
    +0.25*CoeffDispu(5)*Fields(2,1,xi1,y0,zi1)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,xi1,yi1,z0)*Fields(2,1,xi1,yi1,z0) &
    -0.25*CoeffDispu(2)*Fields(1,2,xi1,yi1,z0)*Fields(2,1,xi1,yi1,z0) &
    -0.25*CoeffDispu(5)*Fields(2,1,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,xi1,yi1,zi1)*Fields(2,1,xi1,yi1,zi1) &
    -0.25*CoeffDispu(2)*Fields(1,2,xi1,yi1,zi1)*Fields(2,1,xi1,yi1,zi1) &
    -0.25*CoeffDispu(5)*Fields(2,1,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +0.25*CoeffDispu(3)*Fields(1,2,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +0.25*CoeffDispu(8)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +0.25*CoeffDispu(9)*Fields(2,2,x0,y0,z0)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,x0,y0,zi1)*Fields(2,2,x0,y0,zi1) &
    +0.25*CoeffDispu(3)*Fields(1,2,x0,y0,zi1)*Fields(2,2,x0,y0,zi1) &
    +0.25*CoeffDispu(8)*Fields(2,1,x0,y0,zi1)*Fields(2,2,x0,y0,zi1) &
    +0.25*CoeffDispu(9)*Fields(2,2,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,x0,yi1,z0)*Fields(2,2,x0,yi1,z0) &
    +0.25*CoeffDispu(3)*Fields(1,2,x0,yi1,z0)*Fields(2,2,x0,yi1,z0) &
    -0.25*CoeffDispu(8)*Fields(2,1,x0,yi1,z0)*Fields(2,2,x0,yi1,z0) &
    -0.25*CoeffDispu(9)*Fields(2,2,x0,yi1,z0)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,x0,yi1,zi1)*Fields(2,2,x0,yi1,zi1) &
    +0.25*CoeffDispu(3)*Fields(1,2,x0,yi1,zi1)*Fields(2,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(8)*Fields(2,1,x0,yi1,zi1)*Fields(2,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(9)*Fields(2,2,x0,yi1,zi1)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,xi1,y0,z0)*Fields(2,2,xi1,y0,z0) &
    -0.25*CoeffDispu(3)*Fields(1,2,xi1,y0,z0)*Fields(2,2,xi1,y0,z0) &
    +0.25*CoeffDispu(8)*Fields(2,1,xi1,y0,z0)*Fields(2,2,xi1,y0,z0) &
    +0.25*CoeffDispu(9)*Fields(2,2,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,xi1,y0,zi1)*Fields(2,2,xi1,y0,zi1) &
    -0.25*CoeffDispu(3)*Fields(1,2,xi1,y0,zi1)*Fields(2,2,xi1,y0,zi1) &
    +0.25*CoeffDispu(8)*Fields(2,1,xi1,y0,zi1)*Fields(2,2,xi1,y0,zi1) &
    +0.25*CoeffDispu(9)*Fields(2,2,xi1,y0,zi1)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,xi1,yi1,z0)*Fields(2,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(3)*Fields(1,2,xi1,yi1,z0)*Fields(2,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(8)*Fields(2,1,xi1,yi1,z0)*Fields(2,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(9)*Fields(2,2,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,xi1,yi1,zi1)*Fields(2,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(3)*Fields(1,2,xi1,yi1,zi1)*Fields(2,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(8)*Fields(2,1,xi1,yi1,zi1)*Fields(2,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(9)*Fields(2,2,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
    +0.25*CoeffDispu(2)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
    +0.25*CoeffDispu(4)*Fields(3,1,x0,y0,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(2,1,x0,y0,zi1)*Fields(3,1,x0,y0,zi1) &
    -0.25*CoeffDispu(2)*Fields(2,2,x0,y0,zi1)*Fields(3,1,x0,y0,zi1) &
    +0.25*CoeffDispu(4)*Fields(3,1,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(2,1,x0,yi1,z0)*Fields(3,1,x0,yi1,z0) &
    +0.25*CoeffDispu(2)*Fields(2,2,x0,yi1,z0)*Fields(3,1,x0,yi1,z0) &
    -0.25*CoeffDispu(4)*Fields(3,1,x0,yi1,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(2,1,x0,yi1,zi1)*Fields(3,1,x0,yi1,zi1) &
    -0.25*CoeffDispu(2)*Fields(2,2,x0,yi1,zi1)*Fields(3,1,x0,yi1,zi1) &
    -0.25*CoeffDispu(4)*Fields(3,1,x0,yi1,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(2,1,xi1,y0,z0)*Fields(3,1,xi1,y0,z0) &
    +0.25*CoeffDispu(2)*Fields(2,2,xi1,y0,z0)*Fields(3,1,xi1,y0,z0) &
    +0.25*CoeffDispu(4)*Fields(3,1,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(2,1,xi1,y0,zi1)*Fields(3,1,xi1,y0,zi1) &
    -0.25*CoeffDispu(2)*Fields(2,2,xi1,y0,zi1)*Fields(3,1,xi1,y0,zi1) &
    +0.25*CoeffDispu(4)*Fields(3,1,xi1,y0,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(2,1,xi1,yi1,z0)*Fields(3,1,xi1,yi1,z0) &
    +0.25*CoeffDispu(2)*Fields(2,2,xi1,yi1,z0)*Fields(3,1,xi1,yi1,z0) &
    -0.25*CoeffDispu(4)*Fields(3,1,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(2,1,xi1,yi1,zi1)*Fields(3,1,xi1,yi1,zi1) &
    -0.25*CoeffDispu(2)*Fields(2,2,xi1,yi1,zi1)*Fields(3,1,xi1,yi1,zi1) &
    -0.25*CoeffDispu(4)*Fields(3,1,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(2,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +0.25*CoeffDispu(3)*Fields(2,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +0.25*CoeffDispu(6)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +0.25*CoeffDispu(7)*Fields(3,2,x0,y0,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(2,1,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
    -0.25*CoeffDispu(3)*Fields(2,2,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
    +0.25*CoeffDispu(6)*Fields(3,1,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
    +0.25*CoeffDispu(7)*Fields(3,2,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(2,1,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
    +0.25*CoeffDispu(3)*Fields(2,2,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
    -0.25*CoeffDispu(6)*Fields(3,1,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
    -0.25*CoeffDispu(7)*Fields(3,2,x0,yi1,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(2,1,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(3)*Fields(2,2,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(6)*Fields(3,1,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(7)*Fields(3,2,x0,yi1,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(2,1,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
    +0.25*CoeffDispu(3)*Fields(2,2,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
    +0.25*CoeffDispu(6)*Fields(3,1,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
    +0.25*CoeffDispu(7)*Fields(3,2,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(2,1,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
    -0.25*CoeffDispu(3)*Fields(2,2,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
    +0.25*CoeffDispu(6)*Fields(3,1,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
    +0.25*CoeffDispu(7)*Fields(3,2,xi1,y0,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(2,1,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
    +0.25*CoeffDispu(3)*Fields(2,2,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(6)*Fields(3,1,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(7)*Fields(3,2,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(2,1,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(3)*Fields(2,2,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(6)*Fields(3,1,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(7)*Fields(3,2,xi1,yi1,zi1)**2
  
  ForcesDispu(3,3) = &
    ForcesDispu(3,3) &
    +0.25*CoeffDispu(4)*Fields(1,1,x0,y0,z0)**2 &
    -0.25*CoeffDispu(4)*Fields(1,1,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(4)*Fields(1,1,x0,yi1,z0)**2 &
    -0.25*CoeffDispu(4)*Fields(1,1,x0,yi1,zi1)**2 &
    +0.25*CoeffDispu(4)*Fields(1,1,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(4)*Fields(1,1,xi1,y0,zi1)**2 &
    +0.25*CoeffDispu(4)*Fields(1,1,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(4)*Fields(1,1,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(6)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
    +0.25*CoeffDispu(7)*Fields(1,2,x0,y0,z0)**2 &
    -0.25*CoeffDispu(6)*Fields(1,1,x0,y0,zi1)*Fields(1,2,x0,y0,zi1) &
    -0.25*CoeffDispu(7)*Fields(1,2,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(6)*Fields(1,1,x0,yi1,z0)*Fields(1,2,x0,yi1,z0) &
    +0.25*CoeffDispu(7)*Fields(1,2,x0,yi1,z0)**2 &
    -0.25*CoeffDispu(6)*Fields(1,1,x0,yi1,zi1)*Fields(1,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(7)*Fields(1,2,x0,yi1,zi1)**2 &
    +0.25*CoeffDispu(6)*Fields(1,1,xi1,y0,z0)*Fields(1,2,xi1,y0,z0) &
    +0.25*CoeffDispu(7)*Fields(1,2,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(6)*Fields(1,1,xi1,y0,zi1)*Fields(1,2,xi1,y0,zi1) &
    -0.25*CoeffDispu(7)*Fields(1,2,xi1,y0,zi1)**2 &
    +0.25*CoeffDispu(6)*Fields(1,1,xi1,yi1,z0)*Fields(1,2,xi1,yi1,z0) &
    +0.25*CoeffDispu(7)*Fields(1,2,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(6)*Fields(1,1,xi1,yi1,zi1)*Fields(1,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(7)*Fields(1,2,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(4)*Fields(2,1,x0,y0,z0)**2 &
    -0.25*CoeffDispu(4)*Fields(2,1,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(4)*Fields(2,1,x0,yi1,z0)**2 &
    -0.25*CoeffDispu(4)*Fields(2,1,x0,yi1,zi1)**2 &
    +0.25*CoeffDispu(4)*Fields(2,1,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(4)*Fields(2,1,xi1,y0,zi1)**2 &
    +0.25*CoeffDispu(4)*Fields(2,1,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(4)*Fields(2,1,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(6)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +0.25*CoeffDispu(7)*Fields(2,2,x0,y0,z0)**2 &
    -0.25*CoeffDispu(6)*Fields(2,1,x0,y0,zi1)*Fields(2,2,x0,y0,zi1) &
    -0.25*CoeffDispu(7)*Fields(2,2,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(6)*Fields(2,1,x0,yi1,z0)*Fields(2,2,x0,yi1,z0) &
    +0.25*CoeffDispu(7)*Fields(2,2,x0,yi1,z0)**2 &
    -0.25*CoeffDispu(6)*Fields(2,1,x0,yi1,zi1)*Fields(2,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(7)*Fields(2,2,x0,yi1,zi1)**2 &
    +0.25*CoeffDispu(6)*Fields(2,1,xi1,y0,z0)*Fields(2,2,xi1,y0,z0) &
    +0.25*CoeffDispu(7)*Fields(2,2,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(6)*Fields(2,1,xi1,y0,zi1)*Fields(2,2,xi1,y0,zi1) &
    -0.25*CoeffDispu(7)*Fields(2,2,xi1,y0,zi1)**2 &
    +0.25*CoeffDispu(6)*Fields(2,1,xi1,yi1,z0)*Fields(2,2,xi1,yi1,z0) &
    +0.25*CoeffDispu(7)*Fields(2,2,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(6)*Fields(2,1,xi1,yi1,zi1)*Fields(2,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(7)*Fields(2,2,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
    +0.25*CoeffDispu(2)*Fields(1,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
    +0.25*CoeffDispu(1)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
    +0.25*CoeffDispu(2)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
    +0.25*CoeffDispu(5)*Fields(3,1,x0,y0,z0)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,x0,y0,zi1)*Fields(3,1,x0,y0,zi1) &
    +0.25*CoeffDispu(2)*Fields(1,2,x0,y0,zi1)*Fields(3,1,x0,y0,zi1) &
    +0.25*CoeffDispu(1)*Fields(2,1,x0,y0,zi1)*Fields(3,1,x0,y0,zi1) &
    +0.25*CoeffDispu(2)*Fields(2,2,x0,y0,zi1)*Fields(3,1,x0,y0,zi1) &
    -0.25*CoeffDispu(5)*Fields(3,1,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,x0,yi1,z0)*Fields(3,1,x0,yi1,z0) &
    +0.25*CoeffDispu(2)*Fields(1,2,x0,yi1,z0)*Fields(3,1,x0,yi1,z0) &
    -0.25*CoeffDispu(1)*Fields(2,1,x0,yi1,z0)*Fields(3,1,x0,yi1,z0) &
    -0.25*CoeffDispu(2)*Fields(2,2,x0,yi1,z0)*Fields(3,1,x0,yi1,z0) &
    +0.25*CoeffDispu(5)*Fields(3,1,x0,yi1,z0)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,x0,yi1,zi1)*Fields(3,1,x0,yi1,zi1) &
    +0.25*CoeffDispu(2)*Fields(1,2,x0,yi1,zi1)*Fields(3,1,x0,yi1,zi1) &
    -0.25*CoeffDispu(1)*Fields(2,1,x0,yi1,zi1)*Fields(3,1,x0,yi1,zi1) &
    -0.25*CoeffDispu(2)*Fields(2,2,x0,yi1,zi1)*Fields(3,1,x0,yi1,zi1) &
    -0.25*CoeffDispu(5)*Fields(3,1,x0,yi1,zi1)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,xi1,y0,z0)*Fields(3,1,xi1,y0,z0) &
    -0.25*CoeffDispu(2)*Fields(1,2,xi1,y0,z0)*Fields(3,1,xi1,y0,z0) &
    +0.25*CoeffDispu(1)*Fields(2,1,xi1,y0,z0)*Fields(3,1,xi1,y0,z0) &
    +0.25*CoeffDispu(2)*Fields(2,2,xi1,y0,z0)*Fields(3,1,xi1,y0,z0) &
    +0.25*CoeffDispu(5)*Fields(3,1,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,xi1,y0,zi1)*Fields(3,1,xi1,y0,zi1) &
    -0.25*CoeffDispu(2)*Fields(1,2,xi1,y0,zi1)*Fields(3,1,xi1,y0,zi1) &
    +0.25*CoeffDispu(1)*Fields(2,1,xi1,y0,zi1)*Fields(3,1,xi1,y0,zi1) &
    +0.25*CoeffDispu(2)*Fields(2,2,xi1,y0,zi1)*Fields(3,1,xi1,y0,zi1) &
    -0.25*CoeffDispu(5)*Fields(3,1,xi1,y0,zi1)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,xi1,yi1,z0)*Fields(3,1,xi1,yi1,z0) &
    -0.25*CoeffDispu(2)*Fields(1,2,xi1,yi1,z0)*Fields(3,1,xi1,yi1,z0) &
    -0.25*CoeffDispu(1)*Fields(2,1,xi1,yi1,z0)*Fields(3,1,xi1,yi1,z0) &
    -0.25*CoeffDispu(2)*Fields(2,2,xi1,yi1,z0)*Fields(3,1,xi1,yi1,z0) &
    +0.25*CoeffDispu(5)*Fields(3,1,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,xi1,yi1,zi1)*Fields(3,1,xi1,yi1,zi1) &
    -0.25*CoeffDispu(2)*Fields(1,2,xi1,yi1,zi1)*Fields(3,1,xi1,yi1,zi1) &
    -0.25*CoeffDispu(1)*Fields(2,1,xi1,yi1,zi1)*Fields(3,1,xi1,yi1,zi1) &
    -0.25*CoeffDispu(2)*Fields(2,2,xi1,yi1,zi1)*Fields(3,1,xi1,yi1,zi1) &
    -0.25*CoeffDispu(5)*Fields(3,1,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +0.25*CoeffDispu(3)*Fields(1,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +0.25*CoeffDispu(2)*Fields(2,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +0.25*CoeffDispu(3)*Fields(2,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +0.25*CoeffDispu(8)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +0.25*CoeffDispu(9)*Fields(3,2,x0,y0,z0)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
    +0.25*CoeffDispu(3)*Fields(1,2,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
    +0.25*CoeffDispu(2)*Fields(2,1,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
    +0.25*CoeffDispu(3)*Fields(2,2,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
    -0.25*CoeffDispu(8)*Fields(3,1,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
    -0.25*CoeffDispu(9)*Fields(3,2,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
    +0.25*CoeffDispu(3)*Fields(1,2,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
    -0.25*CoeffDispu(2)*Fields(2,1,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
    -0.25*CoeffDispu(3)*Fields(2,2,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
    +0.25*CoeffDispu(8)*Fields(3,1,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
    +0.25*CoeffDispu(9)*Fields(3,2,x0,yi1,z0)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
    +0.25*CoeffDispu(3)*Fields(1,2,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(2)*Fields(2,1,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(3)*Fields(2,2,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(8)*Fields(3,1,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(9)*Fields(3,2,x0,yi1,zi1)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
    -0.25*CoeffDispu(3)*Fields(1,2,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
    +0.25*CoeffDispu(2)*Fields(2,1,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
    +0.25*CoeffDispu(3)*Fields(2,2,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
    +0.25*CoeffDispu(8)*Fields(3,1,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
    +0.25*CoeffDispu(9)*Fields(3,2,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
    -0.25*CoeffDispu(3)*Fields(1,2,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
    +0.25*CoeffDispu(2)*Fields(2,1,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
    +0.25*CoeffDispu(3)*Fields(2,2,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
    -0.25*CoeffDispu(8)*Fields(3,1,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
    -0.25*CoeffDispu(9)*Fields(3,2,xi1,y0,zi1)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(3)*Fields(1,2,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(2)*Fields(2,1,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(3)*Fields(2,2,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
    +0.25*CoeffDispu(8)*Fields(3,1,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
    +0.25*CoeffDispu(9)*Fields(3,2,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(3)*Fields(1,2,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(2)*Fields(2,1,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(3)*Fields(2,2,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(8)*Fields(3,1,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(9)*Fields(3,2,xi1,yi1,zi1)**2
  
  ForcesDispu(1,4) = &
    ForcesDispu(1,4) &
    +0.
  
  ForcesDispu(2,4) = &
    ForcesDispu(2,4) &
    +0.
  
  ForcesDispu(3,4) = &
    ForcesDispu(3,4) &
    +0.
  
  ForcesDispu(4,4) = &
    ForcesDispu(4,4) &
    +0.
  
  ForcesDispu(5,4) = &
    ForcesDispu(5,4) &
    +0.
  
  ForcesDispu(6,4) = &
    ForcesDispu(6,4) &
    +0.
  
  
End Function GetForcesDispu
