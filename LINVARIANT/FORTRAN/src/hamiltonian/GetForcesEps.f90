Function GetForcesEps(x0, y0, z0, Fields, e0ij) Result(ForcesEps)
  
  Implicit none
  Real*8,  Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(in) :: e0ij(3,3)
  Integer, Intent(in) :: x0, y0, z0
  Real*8              :: eij(3,3), euij(3,3)
  Real*8              :: ForcesEps(Max(FieldDim, 6), NumField+1)
  
  Integer             :: x1,y1,z1,xi1,yi1,zi1
  
    
  ForcesEps = 0.0D0
    euij = StrainFromu(x0, y0, z0, Fields)
    eij = e0ij + euij
    
  x1 = (x0+1)-floor(real(x0+1-1)/real(cgrid%n1))*cgrid%n1
  y1 = (y0+1)-floor(real(y0+1-1)/real(cgrid%n2))*cgrid%n2
  z1 = (z0+1)-floor(real(z0+1-1)/real(cgrid%n3))*cgrid%n3
  xi1 = (x0-1)-floor(real(x0-1-1)/real(cgrid%n1))*cgrid%n1
  yi1 = (y0-1)-floor(real(y0-1-1)/real(cgrid%n2))*cgrid%n2
  zi1 = (z0-1)-floor(real(z0-1-1)/real(cgrid%n3))*cgrid%n3
  
  ForcesEps(1,1) = &
    ForcesEps(1,1) &
    +0.
  
  ForcesEps(2,1) = &
    ForcesEps(2,1) &
    +0.
  
  ForcesEps(3,1) = &
    ForcesEps(3,1) &
    +0.
  
  ForcesEps(1,2) = &
    ForcesEps(1,2) &
    +0.
  
  ForcesEps(2,2) = &
    ForcesEps(2,2) &
    +0.
  
  ForcesEps(3,2) = &
    ForcesEps(3,2) &
    +0.
  
  ForcesEps(1,3) = &
    ForcesEps(1,3) &
    +0.
  
  ForcesEps(2,3) = &
    ForcesEps(2,3) &
    +0.
  
  ForcesEps(3,3) = &
    ForcesEps(3,3) &
    +0.
  
  ForcesEps(1,4) = &
    ForcesEps(1,4) &
    -CoeffEps(1) &
    -(CoeffEps(4)*e0ij(1,1)) &
    -(CoeffEps(3)*e0ij(2,2)) &
    -(CoeffEps(3)*e0ij(3,3)) 
  
  ForcesEps(2,4) = &
    ForcesEps(2,4) &
    -CoeffEps(1) &
    -(CoeffEps(3)*e0ij(1,1)) &
    -(CoeffEps(4)*e0ij(2,2)) &
    -(CoeffEps(3)*e0ij(3,3)) 
  
  ForcesEps(3,4) = &
    ForcesEps(3,4) &
    -CoeffEps(1) &
    -(CoeffEps(3)*e0ij(1,1)) &
    -(CoeffEps(3)*e0ij(2,2)) &
    -(CoeffEps(4)*e0ij(3,3)) 
  
  ForcesEps(4,4) = &
    ForcesEps(4,4) &
    -(CoeffEps(2)*e0ij(2,3)) 
  
  ForcesEps(5,4) = &
    ForcesEps(5,4) &
    -(CoeffEps(2)*e0ij(1,3)) 
  
  ForcesEps(6,4) = &
    ForcesEps(6,4) &
    -(CoeffEps(2)*e0ij(1,2)) 
  
  
End Function GetForcesEps
