Function GetForcesu(x0, y0, z0, Fields, e0ij) Result(Forcesu)
  
  Implicit none
  Real*8,  Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(in) :: e0ij(3,3)
  Integer, Intent(in) :: x0, y0, z0
  Real*8              :: eij(3,3), euij(3,3)
  Real*8              :: Forcesu(Max(FieldDim, 6), NumField+1)
  
  Integer             :: x1,y1,z1,xi1,yi1,zi1
  
    
  Forcesu = 0.0D0
    euij = StrainFromu(x0, y0, z0, Fields)
    eij = e0ij + euij
    
  x1 = (x0+1)-floor(real(x0+1-1)/real(cgrid%n1))*cgrid%n1
  y1 = (y0+1)-floor(real(y0+1-1)/real(cgrid%n2))*cgrid%n2
  z1 = (z0+1)-floor(real(z0+1-1)/real(cgrid%n3))*cgrid%n3
  xi1 = (x0-1)-floor(real(x0-1-1)/real(cgrid%n1))*cgrid%n1
  yi1 = (y0-1)-floor(real(y0-1-1)/real(cgrid%n2))*cgrid%n2
  zi1 = (z0-1)-floor(real(z0-1-1)/real(cgrid%n3))*cgrid%n3
  
  Forcesu(1,1) = &
    Forcesu(1,1) &
    +0.
  
  Forcesu(2,1) = &
    Forcesu(2,1) &
    +0.
  
  Forcesu(3,1) = &
    Forcesu(3,1) &
    +0.
  
  Forcesu(1,2) = &
    Forcesu(1,2) &
    +0.
  
  Forcesu(2,2) = &
    Forcesu(2,2) &
    +0.
  
  Forcesu(3,2) = &
    Forcesu(3,2) &
    +0.
  
  Forcesu(1,3) = &
    Forcesu(1,3) &
    -4.*Coeffu(2)*Fields(1,3,x0,y0,z0) &
    -2.*Coeffu(4)*Fields(1,3,x0,y0,z0) &
    +Coeffu(2)*Fields(1,3,x0,y0,z1) &
    +Coeffu(2)*Fields(1,3,x0,y0,zi1) &
    +Coeffu(2)*Fields(1,3,x0,y1,z0) &
    +Coeffu(2)*Fields(1,3,x0,yi1,z0) &
    +Coeffu(4)*Fields(1,3,x1,y0,z0) &
    +Coeffu(4)*Fields(1,3,xi1,y0,z0) &
    +0.25*Coeffu(2)*Fields(2,3,x1,y1,z0) &
    +0.25*Coeffu(3)*Fields(2,3,x1,y1,z0) &
    -0.25*Coeffu(2)*Fields(2,3,x1,yi1,z0) &
    -0.25*Coeffu(3)*Fields(2,3,x1,yi1,z0) &
    -0.25*Coeffu(2)*Fields(2,3,xi1,y1,z0) &
    -0.25*Coeffu(3)*Fields(2,3,xi1,y1,z0) &
    +0.25*Coeffu(2)*Fields(2,3,xi1,yi1,z0) &
    +0.25*Coeffu(3)*Fields(2,3,xi1,yi1,z0) &
    +0.25*Coeffu(2)*Fields(3,3,x1,y0,z1) &
    +0.25*Coeffu(3)*Fields(3,3,x1,y0,z1) &
    -0.25*Coeffu(2)*Fields(3,3,x1,y0,zi1) &
    -0.25*Coeffu(3)*Fields(3,3,x1,y0,zi1) &
    -0.25*Coeffu(2)*Fields(3,3,xi1,y0,z1) &
    -0.25*Coeffu(3)*Fields(3,3,xi1,y0,z1) &
    +0.25*Coeffu(2)*Fields(3,3,xi1,y0,zi1) &
    +0.25*Coeffu(3)*Fields(3,3,xi1,y0,zi1)
  
  Forcesu(2,3) = &
    Forcesu(2,3) &
    +0.25*Coeffu(2)*Fields(1,3,x1,y1,z0) &
    +0.25*Coeffu(3)*Fields(1,3,x1,y1,z0) &
    -0.25*Coeffu(2)*Fields(1,3,x1,yi1,z0) &
    -0.25*Coeffu(3)*Fields(1,3,x1,yi1,z0) &
    -0.25*Coeffu(2)*Fields(1,3,xi1,y1,z0) &
    -0.25*Coeffu(3)*Fields(1,3,xi1,y1,z0) &
    +0.25*Coeffu(2)*Fields(1,3,xi1,yi1,z0) &
    +0.25*Coeffu(3)*Fields(1,3,xi1,yi1,z0) &
    -4.*Coeffu(2)*Fields(2,3,x0,y0,z0) &
    -2.*Coeffu(4)*Fields(2,3,x0,y0,z0) &
    +Coeffu(2)*Fields(2,3,x0,y0,z1) &
    +Coeffu(2)*Fields(2,3,x0,y0,zi1) &
    +Coeffu(4)*Fields(2,3,x0,y1,z0) &
    +Coeffu(4)*Fields(2,3,x0,yi1,z0) &
    +Coeffu(2)*Fields(2,3,x1,y0,z0) &
    +Coeffu(2)*Fields(2,3,xi1,y0,z0) &
    +0.25*Coeffu(2)*Fields(3,3,x0,y1,z1) &
    +0.25*Coeffu(3)*Fields(3,3,x0,y1,z1) &
    -0.25*Coeffu(2)*Fields(3,3,x0,y1,zi1) &
    -0.25*Coeffu(3)*Fields(3,3,x0,y1,zi1) &
    -0.25*Coeffu(2)*Fields(3,3,x0,yi1,z1) &
    -0.25*Coeffu(3)*Fields(3,3,x0,yi1,z1) &
    +0.25*Coeffu(2)*Fields(3,3,x0,yi1,zi1) &
    +0.25*Coeffu(3)*Fields(3,3,x0,yi1,zi1)
  
  Forcesu(3,3) = &
    Forcesu(3,3) &
    +0.25*Coeffu(2)*Fields(1,3,x1,y0,z1) &
    +0.25*Coeffu(3)*Fields(1,3,x1,y0,z1) &
    -0.25*Coeffu(2)*Fields(1,3,x1,y0,zi1) &
    -0.25*Coeffu(3)*Fields(1,3,x1,y0,zi1) &
    -0.25*Coeffu(2)*Fields(1,3,xi1,y0,z1) &
    -0.25*Coeffu(3)*Fields(1,3,xi1,y0,z1) &
    +0.25*Coeffu(2)*Fields(1,3,xi1,y0,zi1) &
    +0.25*Coeffu(3)*Fields(1,3,xi1,y0,zi1) &
    +0.25*Coeffu(2)*Fields(2,3,x0,y1,z1) &
    +0.25*Coeffu(3)*Fields(2,3,x0,y1,z1) &
    -0.25*Coeffu(2)*Fields(2,3,x0,y1,zi1) &
    -0.25*Coeffu(3)*Fields(2,3,x0,y1,zi1) &
    -0.25*Coeffu(2)*Fields(2,3,x0,yi1,z1) &
    -0.25*Coeffu(3)*Fields(2,3,x0,yi1,z1) &
    +0.25*Coeffu(2)*Fields(2,3,x0,yi1,zi1) &
    +0.25*Coeffu(3)*Fields(2,3,x0,yi1,zi1) &
    -4.*Coeffu(2)*Fields(3,3,x0,y0,z0) &
    -2.*Coeffu(4)*Fields(3,3,x0,y0,z0) &
    +Coeffu(4)*Fields(3,3,x0,y0,z1) &
    +Coeffu(4)*Fields(3,3,x0,y0,zi1) &
    +Coeffu(2)*Fields(3,3,x0,y1,z0) &
    +Coeffu(2)*Fields(3,3,x0,yi1,z0) &
    +Coeffu(2)*Fields(3,3,x1,y0,z0) &
    +Coeffu(2)*Fields(3,3,xi1,y0,z0)
  
  Forcesu(1,4) = &
    Forcesu(1,4) &
    +0.
  
  Forcesu(2,4) = &
    Forcesu(2,4) &
    +0.
  
  Forcesu(3,4) = &
    Forcesu(3,4) &
    +0.
  
  Forcesu(4,4) = &
    Forcesu(4,4) &
    +0.
  
  Forcesu(5,4) = &
    Forcesu(5,4) &
    +0.
  
  Forcesu(6,4) = &
    Forcesu(6,4) &
    +0.
  
  
End Function GetForcesu
