Function GetForcesEpsu(x0, y0, z0, Fields, e0ij) Result(ForcesEpsu)
  
  Implicit none
  Real*8,  Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(in) :: e0ij(3,3)
  Integer, Intent(in) :: x0, y0, z0
  Real*8              :: eij(3,3), euij(3,3)
  Real*8              :: ForcesEpsu(Max(FieldDim, 6), NumField+1)
  
  Integer             :: x1,y1,z1,xi1,yi1,zi1
  
    
  ForcesEpsu = 0.0D0
    euij = StrainFromu(x0, y0, z0, Fields)
    eij = e0ij + euij
    
  x1 = (x0+1)-floor(real(x0+1-1)/real(cgrid%n1))*cgrid%n1
  y1 = (y0+1)-floor(real(y0+1-1)/real(cgrid%n2))*cgrid%n2
  z1 = (z0+1)-floor(real(z0+1-1)/real(cgrid%n3))*cgrid%n3
  xi1 = (x0-1)-floor(real(x0-1-1)/real(cgrid%n1))*cgrid%n1
  yi1 = (y0-1)-floor(real(y0-1-1)/real(cgrid%n2))*cgrid%n2
  zi1 = (z0-1)-floor(real(z0-1-1)/real(cgrid%n3))*cgrid%n3
  
  ForcesEpsu(1,1) = &
    ForcesEpsu(1,1) &
    +0.
  
  ForcesEpsu(2,1) = &
    ForcesEpsu(2,1) &
    +0.
  
  ForcesEpsu(3,1) = &
    ForcesEpsu(3,1) &
    +0.
  
  ForcesEpsu(1,2) = &
    ForcesEpsu(1,2) &
    +0.
  
  ForcesEpsu(2,2) = &
    ForcesEpsu(2,2) &
    +0.
  
  ForcesEpsu(3,2) = &
    ForcesEpsu(3,2) &
    +0.
  
  ForcesEpsu(1,3) = &
    ForcesEpsu(1,3) &
    +0.
  
  ForcesEpsu(2,3) = &
    ForcesEpsu(2,3) &
    +0.
  
  ForcesEpsu(3,3) = &
    ForcesEpsu(3,3) &
    +0.
  
  ForcesEpsu(1,4) = &
    ForcesEpsu(1,4) +0. 
  
  ForcesEpsu(2,4) = &
    ForcesEpsu(2,4) +0.
  
  ForcesEpsu(3,4) = &
    ForcesEpsu(3,4) +0.
  
  ForcesEpsu(4,4) = &
    ForcesEpsu(4,4) +0.
  
  ForcesEpsu(5,4) = &
    ForcesEpsu(5,4) +0.
  
  ForcesEpsu(6,4) = &
    ForcesEpsu(6,4) +0.
  
  
End Function GetForcesEpsu
