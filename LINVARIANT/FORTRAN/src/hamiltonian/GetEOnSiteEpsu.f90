Function GetEOnSiteEpsu(x0, y0, z0, Fields, e0ij) Result(EOnSiteEpsu)
  
  Implicit none
  Real*8,  Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(in) :: e0ij(3,3)
  Integer, Intent(in) :: x0, y0, z0
  Real*8              :: eij(3,3), euij(3,3)
  Real*8              :: EOnSiteEpsu
  
  Integer             :: x1,y1,z1,xi1,yi1,zi1
  
  EOnSiteEpsu = 0.0D0
  
  
  euij = StrainFromu(x0, y0, z0, Fields)
  eij = e0ij + euij
  
  x1 = (x0+1)-floor(real(x0+1-1)/real(cgrid%n1))*cgrid%n1
  y1 = (y0+1)-floor(real(y0+1-1)/real(cgrid%n2))*cgrid%n2
  z1 = (z0+1)-floor(real(z0+1-1)/real(cgrid%n3))*cgrid%n3
  xi1 = (x0-1)-floor(real(x0-1-1)/real(cgrid%n1))*cgrid%n1
  yi1 = (y0-1)-floor(real(y0-1-1)/real(cgrid%n2))*cgrid%n2
  zi1 = (z0-1)-floor(real(z0-1-1)/real(cgrid%n3))*cgrid%n3
  
  EOnSiteEpsu = &
    EOnSiteEpsu &
    -0.25*CoeffEpsu(4)*e0ij(1,1)*Fields(1,3,x0,y0,z0) &
    -0.25*CoeffEpsu(2)*e0ij(1,2)*Fields(1,3,x0,y0,z0) &
    -0.25*CoeffEpsu(2)*e0ij(1,3)*Fields(1,3,x0,y0,z0) &
    -0.25*CoeffEpsu(3)*e0ij(2,2)*Fields(1,3,x0,y0,z0) &
    -0.25*CoeffEpsu(3)*e0ij(3,3)*Fields(1,3,x0,y0,z0) &
    -0.25*CoeffEpsu(4)*e0ij(1,1)*Fields(1,3,x0,y0,z1) &
    -0.25*CoeffEpsu(2)*e0ij(1,2)*Fields(1,3,x0,y0,z1) &
    +0.25*CoeffEpsu(2)*e0ij(1,3)*Fields(1,3,x0,y0,z1) &
    -0.25*CoeffEpsu(3)*e0ij(2,2)*Fields(1,3,x0,y0,z1) &
    -0.25*CoeffEpsu(3)*e0ij(3,3)*Fields(1,3,x0,y0,z1) &
    -0.25*CoeffEpsu(4)*e0ij(1,1)*Fields(1,3,x0,y1,z0) &
    +0.25*CoeffEpsu(2)*e0ij(1,2)*Fields(1,3,x0,y1,z0) &
    -0.25*CoeffEpsu(2)*e0ij(1,3)*Fields(1,3,x0,y1,z0) &
    -0.25*CoeffEpsu(3)*e0ij(2,2)*Fields(1,3,x0,y1,z0) &
    -0.25*CoeffEpsu(3)*e0ij(3,3)*Fields(1,3,x0,y1,z0) &
    -0.25*CoeffEpsu(4)*e0ij(1,1)*Fields(1,3,x0,y1,z1) &
    +0.25*CoeffEpsu(2)*e0ij(1,2)*Fields(1,3,x0,y1,z1) &
    +0.25*CoeffEpsu(2)*e0ij(1,3)*Fields(1,3,x0,y1,z1) &
    -0.25*CoeffEpsu(3)*e0ij(2,2)*Fields(1,3,x0,y1,z1) &
    -0.25*CoeffEpsu(3)*e0ij(3,3)*Fields(1,3,x0,y1,z1) &
    +0.25*CoeffEpsu(4)*e0ij(1,1)*Fields(1,3,x1,y0,z0) &
    -0.25*CoeffEpsu(2)*e0ij(1,2)*Fields(1,3,x1,y0,z0) &
    -0.25*CoeffEpsu(2)*e0ij(1,3)*Fields(1,3,x1,y0,z0) &
    +0.25*CoeffEpsu(3)*e0ij(2,2)*Fields(1,3,x1,y0,z0) &
    +0.25*CoeffEpsu(3)*e0ij(3,3)*Fields(1,3,x1,y0,z0) &
    +0.25*CoeffEpsu(4)*e0ij(1,1)*Fields(1,3,x1,y0,z1) &
    -0.25*CoeffEpsu(2)*e0ij(1,2)*Fields(1,3,x1,y0,z1) &
    +0.25*CoeffEpsu(2)*e0ij(1,3)*Fields(1,3,x1,y0,z1) &
    +0.25*CoeffEpsu(3)*e0ij(2,2)*Fields(1,3,x1,y0,z1) &
    +0.25*CoeffEpsu(3)*e0ij(3,3)*Fields(1,3,x1,y0,z1) &
    +0.25*CoeffEpsu(4)*e0ij(1,1)*Fields(1,3,x1,y1,z0) &
    +0.25*CoeffEpsu(2)*e0ij(1,2)*Fields(1,3,x1,y1,z0) &
    -0.25*CoeffEpsu(2)*e0ij(1,3)*Fields(1,3,x1,y1,z0) &
    +0.25*CoeffEpsu(3)*e0ij(2,2)*Fields(1,3,x1,y1,z0) &
    +0.25*CoeffEpsu(3)*e0ij(3,3)*Fields(1,3,x1,y1,z0) &
    +0.25*CoeffEpsu(4)*e0ij(1,1)*Fields(1,3,x1,y1,z1) &
    +0.25*CoeffEpsu(2)*e0ij(1,2)*Fields(1,3,x1,y1,z1) &
    +0.25*CoeffEpsu(2)*e0ij(1,3)*Fields(1,3,x1,y1,z1) &
    +0.25*CoeffEpsu(3)*e0ij(2,2)*Fields(1,3,x1,y1,z1) &
    +0.25*CoeffEpsu(3)*e0ij(3,3)*Fields(1,3,x1,y1,z1) &
    -0.25*CoeffEpsu(3)*e0ij(1,1)*Fields(2,3,x0,y0,z0) &
    -0.25*CoeffEpsu(2)*e0ij(1,2)*Fields(2,3,x0,y0,z0) &
    -0.25*CoeffEpsu(4)*e0ij(2,2)*Fields(2,3,x0,y0,z0) &
    -0.25*CoeffEpsu(2)*e0ij(2,3)*Fields(2,3,x0,y0,z0) &
    -0.25*CoeffEpsu(3)*e0ij(3,3)*Fields(2,3,x0,y0,z0) &
    -0.25*CoeffEpsu(3)*e0ij(1,1)*Fields(2,3,x0,y0,z1) &
    -0.25*CoeffEpsu(2)*e0ij(1,2)*Fields(2,3,x0,y0,z1) &
    -0.25*CoeffEpsu(4)*e0ij(2,2)*Fields(2,3,x0,y0,z1) &
    +0.25*CoeffEpsu(2)*e0ij(2,3)*Fields(2,3,x0,y0,z1) &
    -0.25*CoeffEpsu(3)*e0ij(3,3)*Fields(2,3,x0,y0,z1) &
    +0.25*CoeffEpsu(3)*e0ij(1,1)*Fields(2,3,x0,y1,z0) &
    -0.25*CoeffEpsu(2)*e0ij(1,2)*Fields(2,3,x0,y1,z0) &
    +0.25*CoeffEpsu(4)*e0ij(2,2)*Fields(2,3,x0,y1,z0) &
    -0.25*CoeffEpsu(2)*e0ij(2,3)*Fields(2,3,x0,y1,z0) &
    +0.25*CoeffEpsu(3)*e0ij(3,3)*Fields(2,3,x0,y1,z0) &
    +0.25*CoeffEpsu(3)*e0ij(1,1)*Fields(2,3,x0,y1,z1) &
    -0.25*CoeffEpsu(2)*e0ij(1,2)*Fields(2,3,x0,y1,z1) &
    +0.25*CoeffEpsu(4)*e0ij(2,2)*Fields(2,3,x0,y1,z1) &
    +0.25*CoeffEpsu(2)*e0ij(2,3)*Fields(2,3,x0,y1,z1) &
    +0.25*CoeffEpsu(3)*e0ij(3,3)*Fields(2,3,x0,y1,z1) &
    -0.25*CoeffEpsu(3)*e0ij(1,1)*Fields(2,3,x1,y0,z0) &
    +0.25*CoeffEpsu(2)*e0ij(1,2)*Fields(2,3,x1,y0,z0) &
    -0.25*CoeffEpsu(4)*e0ij(2,2)*Fields(2,3,x1,y0,z0) &
    -0.25*CoeffEpsu(2)*e0ij(2,3)*Fields(2,3,x1,y0,z0) &
    -0.25*CoeffEpsu(3)*e0ij(3,3)*Fields(2,3,x1,y0,z0) &
    -0.25*CoeffEpsu(3)*e0ij(1,1)*Fields(2,3,x1,y0,z1) &
    +0.25*CoeffEpsu(2)*e0ij(1,2)*Fields(2,3,x1,y0,z1) &
    -0.25*CoeffEpsu(4)*e0ij(2,2)*Fields(2,3,x1,y0,z1) &
    +0.25*CoeffEpsu(2)*e0ij(2,3)*Fields(2,3,x1,y0,z1) &
    -0.25*CoeffEpsu(3)*e0ij(3,3)*Fields(2,3,x1,y0,z1) &
    +0.25*CoeffEpsu(3)*e0ij(1,1)*Fields(2,3,x1,y1,z0) &
    +0.25*CoeffEpsu(2)*e0ij(1,2)*Fields(2,3,x1,y1,z0) &
    +0.25*CoeffEpsu(4)*e0ij(2,2)*Fields(2,3,x1,y1,z0) &
    -0.25*CoeffEpsu(2)*e0ij(2,3)*Fields(2,3,x1,y1,z0) &
    +0.25*CoeffEpsu(3)*e0ij(3,3)*Fields(2,3,x1,y1,z0) &
    +0.25*CoeffEpsu(3)*e0ij(1,1)*Fields(2,3,x1,y1,z1) &
    +0.25*CoeffEpsu(2)*e0ij(1,2)*Fields(2,3,x1,y1,z1) &
    +0.25*CoeffEpsu(4)*e0ij(2,2)*Fields(2,3,x1,y1,z1) &
    +0.25*CoeffEpsu(2)*e0ij(2,3)*Fields(2,3,x1,y1,z1) &
    +0.25*CoeffEpsu(3)*e0ij(3,3)*Fields(2,3,x1,y1,z1) &
    -0.25*CoeffEpsu(3)*e0ij(1,1)*Fields(3,3,x0,y0,z0) &
    -0.25*CoeffEpsu(2)*e0ij(1,3)*Fields(3,3,x0,y0,z0) &
    -0.25*CoeffEpsu(3)*e0ij(2,2)*Fields(3,3,x0,y0,z0) &
    -0.25*CoeffEpsu(2)*e0ij(2,3)*Fields(3,3,x0,y0,z0) &
    -0.25*CoeffEpsu(4)*e0ij(3,3)*Fields(3,3,x0,y0,z0) &
    +0.25*CoeffEpsu(3)*e0ij(1,1)*Fields(3,3,x0,y0,z1) &
    -0.25*CoeffEpsu(2)*e0ij(1,3)*Fields(3,3,x0,y0,z1) &
    +0.25*CoeffEpsu(3)*e0ij(2,2)*Fields(3,3,x0,y0,z1) &
    -0.25*CoeffEpsu(2)*e0ij(2,3)*Fields(3,3,x0,y0,z1) &
    +0.25*CoeffEpsu(4)*e0ij(3,3)*Fields(3,3,x0,y0,z1) &
    -0.25*CoeffEpsu(3)*e0ij(1,1)*Fields(3,3,x0,y1,z0) &
    -0.25*CoeffEpsu(2)*e0ij(1,3)*Fields(3,3,x0,y1,z0) &
    -0.25*CoeffEpsu(3)*e0ij(2,2)*Fields(3,3,x0,y1,z0) &
    +0.25*CoeffEpsu(2)*e0ij(2,3)*Fields(3,3,x0,y1,z0) &
    -0.25*CoeffEpsu(4)*e0ij(3,3)*Fields(3,3,x0,y1,z0) &
    +0.25*CoeffEpsu(3)*e0ij(1,1)*Fields(3,3,x0,y1,z1) &
    -0.25*CoeffEpsu(2)*e0ij(1,3)*Fields(3,3,x0,y1,z1) &
    +0.25*CoeffEpsu(3)*e0ij(2,2)*Fields(3,3,x0,y1,z1) &
    +0.25*CoeffEpsu(2)*e0ij(2,3)*Fields(3,3,x0,y1,z1) &
    +0.25*CoeffEpsu(4)*e0ij(3,3)*Fields(3,3,x0,y1,z1) &
    -0.25*CoeffEpsu(3)*e0ij(1,1)*Fields(3,3,x1,y0,z0) &
    +0.25*CoeffEpsu(2)*e0ij(1,3)*Fields(3,3,x1,y0,z0) &
    -0.25*CoeffEpsu(3)*e0ij(2,2)*Fields(3,3,x1,y0,z0) &
    -0.25*CoeffEpsu(2)*e0ij(2,3)*Fields(3,3,x1,y0,z0) &
    -0.25*CoeffEpsu(4)*e0ij(3,3)*Fields(3,3,x1,y0,z0) &
    +0.25*CoeffEpsu(3)*e0ij(1,1)*Fields(3,3,x1,y0,z1) &
    +0.25*CoeffEpsu(2)*e0ij(1,3)*Fields(3,3,x1,y0,z1) &
    +0.25*CoeffEpsu(3)*e0ij(2,2)*Fields(3,3,x1,y0,z1) &
    -0.25*CoeffEpsu(2)*e0ij(2,3)*Fields(3,3,x1,y0,z1) &
    +0.25*CoeffEpsu(4)*e0ij(3,3)*Fields(3,3,x1,y0,z1) &
    -0.25*CoeffEpsu(3)*e0ij(1,1)*Fields(3,3,x1,y1,z0) &
    +0.25*CoeffEpsu(2)*e0ij(1,3)*Fields(3,3,x1,y1,z0) &
    -0.25*CoeffEpsu(3)*e0ij(2,2)*Fields(3,3,x1,y1,z0) &
    +0.25*CoeffEpsu(2)*e0ij(2,3)*Fields(3,3,x1,y1,z0) &
    -0.25*CoeffEpsu(4)*e0ij(3,3)*Fields(3,3,x1,y1,z0) &
    +0.25*CoeffEpsu(3)*e0ij(1,1)*Fields(3,3,x1,y1,z1) &
    +0.25*CoeffEpsu(2)*e0ij(1,3)*Fields(3,3,x1,y1,z1) &
    +0.25*CoeffEpsu(3)*e0ij(2,2)*Fields(3,3,x1,y1,z1) &
    +0.25*CoeffEpsu(2)*e0ij(2,3)*Fields(3,3,x1,y1,z1) &
    +0.25*CoeffEpsu(4)*e0ij(3,3)*Fields(3,3,x1,y1,z1)
  
  
End Function GetEOnSiteEpsu