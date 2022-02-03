Function GetForcesDisp(x0, y0, z0, Fields, e0ij) Result(ForcesDisp)
  
  Implicit none
  Real*8,  Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(in) :: e0ij(3,3)
  Integer, Intent(in) :: x0, y0, z0
  Real*8              :: eij(3,3), euij(3,3)
  Real*8              :: ForcesDisp(Max(FieldDim, 6), NumField+1)
  
  Integer             :: x1,y1,z1,xi1,yi1,zi1
  
    
  ForcesDisp = 0.0D0
    euij = StrainFromu(x0, y0, z0, Fields)
    eij = e0ij + euij
    
  x1 = (x0+1)-floor(real(x0+1-1)/real(cgrid%n1))*cgrid%n1
  y1 = (y0+1)-floor(real(y0+1-1)/real(cgrid%n2))*cgrid%n2
  z1 = (z0+1)-floor(real(z0+1-1)/real(cgrid%n3))*cgrid%n3
  xi1 = (x0-1)-floor(real(x0-1-1)/real(cgrid%n1))*cgrid%n1
  yi1 = (y0-1)-floor(real(y0-1-1)/real(cgrid%n2))*cgrid%n2
  zi1 = (z0-1)-floor(real(z0-1-1)/real(cgrid%n3))*cgrid%n3
  
  ForcesDisp(1,1) = &
    ForcesDisp(1,1) &
    -2.*CoeffDisp(1)*Fields(1,1,x0,y0,z0) &
    -4.*CoeffDisp(5)*Fields(1,1,x0,y0,z0)**3 &
    -(CoeffDisp(2)*Fields(1,2,x0,y0,z0)) &
    -3.*CoeffDisp(7)*Fields(1,1,x0,y0,z0)**2*Fields(1,2,x0,y0,z0) &
    -2.*CoeffDisp(10)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)**2 &
    -(CoeffDisp(13)*Fields(1,2,x0,y0,z0)**3) &
    -2.*CoeffDisp(4)*Fields(1,1,x0,y0,z0)*Fields(2,1,x0,y0,z0)**2 &
    -(CoeffDisp(6)*Fields(1,2,x0,y0,z0)*Fields(2,1,x0,y0,z0)**2) &
    -2.*CoeffDisp(6)*Fields(1,1,x0,y0,z0)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    -(CoeffDisp(8)*Fields(1,2,x0,y0,z0)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0)) &
    -2.*CoeffDisp(9)*Fields(1,1,x0,y0,z0)*Fields(2,2,x0,y0,z0)**2 &
    -(CoeffDisp(11)*Fields(1,2,x0,y0,z0)*Fields(2,2,x0,y0,z0)**2) &
    -2.*CoeffDisp(4)*Fields(1,1,x0,y0,z0)*Fields(3,1,x0,y0,z0)**2 &
    -(CoeffDisp(6)*Fields(1,2,x0,y0,z0)*Fields(3,1,x0,y0,z0)**2) &
    -2.*CoeffDisp(6)*Fields(1,1,x0,y0,z0)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    -(CoeffDisp(8)*Fields(1,2,x0,y0,z0)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0)) &
    -2.*CoeffDisp(9)*Fields(1,1,x0,y0,z0)*Fields(3,2,x0,y0,z0)**2 &
    -(CoeffDisp(11)*Fields(1,2,x0,y0,z0)*Fields(3,2,x0,y0,z0)**2)
  
  ForcesDisp(2,1) = &
    ForcesDisp(2,1) &
    -2.*CoeffDisp(1)*Fields(2,1,x0,y0,z0) &
    -2.*CoeffDisp(4)*Fields(1,1,x0,y0,z0)**2*Fields(2,1,x0,y0,z0) &
    -2.*CoeffDisp(6)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)*Fields(2,1,x0,y0,z0) &
    -2.*CoeffDisp(9)*Fields(1,2,x0,y0,z0)**2*Fields(2,1,x0,y0,z0) &
    -4.*CoeffDisp(5)*Fields(2,1,x0,y0,z0)**3 &
    -(CoeffDisp(2)*Fields(2,2,x0,y0,z0)) &
    -(CoeffDisp(6)*Fields(1,1,x0,y0,z0)**2*Fields(2,2,x0,y0,z0)) &
    -(CoeffDisp(8)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)*Fields(2,2,x0,y0,z0)) &
    -(CoeffDisp(11)*Fields(1,2,x0,y0,z0)**2*Fields(2,2,x0,y0,z0)) &
    -3.*CoeffDisp(7)*Fields(2,1,x0,y0,z0)**2*Fields(2,2,x0,y0,z0) &
    -2.*CoeffDisp(10)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0)**2 &
    -(CoeffDisp(13)*Fields(2,2,x0,y0,z0)**3) &
    -2.*CoeffDisp(4)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,y0,z0)**2 &
    -(CoeffDisp(6)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y0,z0)**2) &
    -2.*CoeffDisp(6)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    -(CoeffDisp(8)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0)) &
    -2.*CoeffDisp(9)*Fields(2,1,x0,y0,z0)*Fields(3,2,x0,y0,z0)**2 &
    -(CoeffDisp(11)*Fields(2,2,x0,y0,z0)*Fields(3,2,x0,y0,z0)**2)
  
  ForcesDisp(3,1) = &
    ForcesDisp(3,1) &
    -2.*CoeffDisp(1)*Fields(3,1,x0,y0,z0) &
    -2.*CoeffDisp(4)*Fields(1,1,x0,y0,z0)**2*Fields(3,1,x0,y0,z0) &
    -2.*CoeffDisp(6)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
    -2.*CoeffDisp(9)*Fields(1,2,x0,y0,z0)**2*Fields(3,1,x0,y0,z0) &
    -2.*CoeffDisp(4)*Fields(2,1,x0,y0,z0)**2*Fields(3,1,x0,y0,z0) &
    -2.*CoeffDisp(6)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
    -2.*CoeffDisp(9)*Fields(2,2,x0,y0,z0)**2*Fields(3,1,x0,y0,z0) &
    -4.*CoeffDisp(5)*Fields(3,1,x0,y0,z0)**3 &
    -(CoeffDisp(2)*Fields(3,2,x0,y0,z0)) &
    -(CoeffDisp(6)*Fields(1,1,x0,y0,z0)**2*Fields(3,2,x0,y0,z0)) &
    -(CoeffDisp(8)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)*Fields(3,2,x0,y0,z0)) &
    -(CoeffDisp(11)*Fields(1,2,x0,y0,z0)**2*Fields(3,2,x0,y0,z0)) &
    -(CoeffDisp(6)*Fields(2,1,x0,y0,z0)**2*Fields(3,2,x0,y0,z0)) &
    -(CoeffDisp(8)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0)*Fields(3,2,x0,y0,z0)) &
    -(CoeffDisp(11)*Fields(2,2,x0,y0,z0)**2*Fields(3,2,x0,y0,z0)) &
    -3.*CoeffDisp(7)*Fields(3,1,x0,y0,z0)**2*Fields(3,2,x0,y0,z0) &
    -2.*CoeffDisp(10)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0)**2 &
    -(CoeffDisp(13)*Fields(3,2,x0,y0,z0)**3)
  
  ForcesDisp(1,2) = &
    ForcesDisp(1,2) &
    -(CoeffDisp(2)*Fields(1,1,x0,y0,z0)) &
    -(CoeffDisp(7)*Fields(1,1,x0,y0,z0)**3) &
    -2.*CoeffDisp(3)*Fields(1,2,x0,y0,z0) &
    -2.*CoeffDisp(10)*Fields(1,1,x0,y0,z0)**2*Fields(1,2,x0,y0,z0) &
    -3.*CoeffDisp(13)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)**2 &
    -4.*CoeffDisp(14)*Fields(1,2,x0,y0,z0)**3 &
    -(CoeffDisp(6)*Fields(1,1,x0,y0,z0)*Fields(2,1,x0,y0,z0)**2) &
    -2.*CoeffDisp(9)*Fields(1,2,x0,y0,z0)*Fields(2,1,x0,y0,z0)**2 &
    -(CoeffDisp(8)*Fields(1,1,x0,y0,z0)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0)) &
    -2.*CoeffDisp(11)*Fields(1,2,x0,y0,z0)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    -(CoeffDisp(11)*Fields(1,1,x0,y0,z0)*Fields(2,2,x0,y0,z0)**2) &
    -2.*CoeffDisp(12)*Fields(1,2,x0,y0,z0)*Fields(2,2,x0,y0,z0)**2 &
    -(CoeffDisp(6)*Fields(1,1,x0,y0,z0)*Fields(3,1,x0,y0,z0)**2) &
    -2.*CoeffDisp(9)*Fields(1,2,x0,y0,z0)*Fields(3,1,x0,y0,z0)**2 &
    -(CoeffDisp(8)*Fields(1,1,x0,y0,z0)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0)) &
    -2.*CoeffDisp(11)*Fields(1,2,x0,y0,z0)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    -(CoeffDisp(11)*Fields(1,1,x0,y0,z0)*Fields(3,2,x0,y0,z0)**2) &
    -2.*CoeffDisp(12)*Fields(1,2,x0,y0,z0)*Fields(3,2,x0,y0,z0)**2
  
  ForcesDisp(2,2) = &
    ForcesDisp(2,2) &
    -(CoeffDisp(2)*Fields(2,1,x0,y0,z0)) &
    -(CoeffDisp(6)*Fields(1,1,x0,y0,z0)**2*Fields(2,1,x0,y0,z0)) &
    -(CoeffDisp(8)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)*Fields(2,1,x0,y0,z0)) &
    -(CoeffDisp(11)*Fields(1,2,x0,y0,z0)**2*Fields(2,1,x0,y0,z0)) &
    -(CoeffDisp(7)*Fields(2,1,x0,y0,z0)**3) &
    -2.*CoeffDisp(3)*Fields(2,2,x0,y0,z0) &
    -2.*CoeffDisp(9)*Fields(1,1,x0,y0,z0)**2*Fields(2,2,x0,y0,z0) &
    -2.*CoeffDisp(11)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    -2.*CoeffDisp(12)*Fields(1,2,x0,y0,z0)**2*Fields(2,2,x0,y0,z0) &
    -2.*CoeffDisp(10)*Fields(2,1,x0,y0,z0)**2*Fields(2,2,x0,y0,z0) &
    -3.*CoeffDisp(13)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0)**2 &
    -4.*CoeffDisp(14)*Fields(2,2,x0,y0,z0)**3 &
    -(CoeffDisp(6)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,y0,z0)**2) &
    -2.*CoeffDisp(9)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y0,z0)**2 &
    -(CoeffDisp(8)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0)) &
    -2.*CoeffDisp(11)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    -(CoeffDisp(11)*Fields(2,1,x0,y0,z0)*Fields(3,2,x0,y0,z0)**2) &
    -2.*CoeffDisp(12)*Fields(2,2,x0,y0,z0)*Fields(3,2,x0,y0,z0)**2
  
  ForcesDisp(3,2) = &
    ForcesDisp(3,2) &
    -(CoeffDisp(2)*Fields(3,1,x0,y0,z0)) &
    -(CoeffDisp(6)*Fields(1,1,x0,y0,z0)**2*Fields(3,1,x0,y0,z0)) &
    -(CoeffDisp(8)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)*Fields(3,1,x0,y0,z0)) &
    -(CoeffDisp(11)*Fields(1,2,x0,y0,z0)**2*Fields(3,1,x0,y0,z0)) &
    -(CoeffDisp(6)*Fields(2,1,x0,y0,z0)**2*Fields(3,1,x0,y0,z0)) &
    -(CoeffDisp(8)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y0,z0)) &
    -(CoeffDisp(11)*Fields(2,2,x0,y0,z0)**2*Fields(3,1,x0,y0,z0)) &
    -(CoeffDisp(7)*Fields(3,1,x0,y0,z0)**3) &
    -2.*CoeffDisp(3)*Fields(3,2,x0,y0,z0) &
    -2.*CoeffDisp(9)*Fields(1,1,x0,y0,z0)**2*Fields(3,2,x0,y0,z0) &
    -2.*CoeffDisp(11)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    -2.*CoeffDisp(12)*Fields(1,2,x0,y0,z0)**2*Fields(3,2,x0,y0,z0) &
    -2.*CoeffDisp(9)*Fields(2,1,x0,y0,z0)**2*Fields(3,2,x0,y0,z0) &
    -2.*CoeffDisp(11)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    -2.*CoeffDisp(12)*Fields(2,2,x0,y0,z0)**2*Fields(3,2,x0,y0,z0) &
    -2.*CoeffDisp(10)*Fields(3,1,x0,y0,z0)**2*Fields(3,2,x0,y0,z0) &
    -3.*CoeffDisp(13)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0)**2 &
    -4.*CoeffDisp(14)*Fields(3,2,x0,y0,z0)**3
  
  ForcesDisp(1,3) = &
    ForcesDisp(1,3) &
    +0.
  
  ForcesDisp(2,3) = &
    ForcesDisp(2,3) &
    +0.
  
  ForcesDisp(3,3) = &
    ForcesDisp(3,3) &
    +0.
  
  ForcesDisp(1,4) = &
    ForcesDisp(1,4) &
    +0.
  
  ForcesDisp(2,4) = &
    ForcesDisp(2,4) &
    +0.
  
  ForcesDisp(3,4) = &
    ForcesDisp(3,4) &
    +0.
  
  ForcesDisp(4,4) = &
    ForcesDisp(4,4) &
    +0.
  
  ForcesDisp(5,4) = &
    ForcesDisp(5,4) &
    +0.
  
  ForcesDisp(6,4) = &
    ForcesDisp(6,4) &
    +0.
  
  
End Function GetForcesDisp
