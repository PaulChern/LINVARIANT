Function StrainFromu(x0, y0, z0, Fields) Result(euij)
  
  Use Parameters
  Implicit none
  Integer, Intent(in) :: x0, y0, z0
  Real*8,  Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8              :: euij(3,3)
  Integer             :: x1,y1,z1,xi1,yi1,zi1

  x1 = (x0+1)-floor(real(x0+1-1)/real(cgrid%n1))*cgrid%n1
  y1 = (y0+1)-floor(real(y0+1-1)/real(cgrid%n2))*cgrid%n2
  z1 = (z0+1)-floor(real(z0+1-1)/real(cgrid%n3))*cgrid%n3
  xi1 = (x0-1)-floor(real(x0-1-1)/real(cgrid%n1))*cgrid%n1
  yi1 = (y0-1)-floor(real(y0-1-1)/real(cgrid%n2))*cgrid%n2
  zi1 = (z0-1)-floor(real(z0-1-1)/real(cgrid%n3))*cgrid%n3
  
  euij = 0.0d0
  
  euij(1,1) = &
    euij(1,1) &
    -0.25*Fields(1,3,x0,y0,z0) &
    -0.25*Fields(1,3,x0,y0,z1) &
    -0.25*Fields(1,3,x0,y1,z0) &
    -0.25*Fields(1,3,x0,y1,z1) &
    +0.25*Fields(1,3,x1,y0,z0) &
    +0.25*Fields(1,3,x1,y0,z1) &
    +0.25*Fields(1,3,x1,y1,z0) &
    +0.25*Fields(1,3,x1,y1,z1)
  
  euij(1,2) = &
    euij(1,2) &
    -0.25*Fields(1,3,x0,y0,z0) &
    -0.25*Fields(1,3,x0,y0,z1) &
    +0.25*Fields(1,3,x0,y1,z0) &
    +0.25*Fields(1,3,x0,y1,z1) &
    -0.25*Fields(1,3,x1,y0,z0) &
    -0.25*Fields(1,3,x1,y0,z1) &
    +0.25*Fields(1,3,x1,y1,z0) &
    +0.25*Fields(1,3,x1,y1,z1) &
    -0.25*Fields(2,3,x0,y0,z0) &
    -0.25*Fields(2,3,x0,y0,z1) &
    -0.25*Fields(2,3,x0,y1,z0) &
    -0.25*Fields(2,3,x0,y1,z1) &
    +0.25*Fields(2,3,x1,y0,z0) &
    +0.25*Fields(2,3,x1,y0,z1) &
    +0.25*Fields(2,3,x1,y1,z0) &
    +0.25*Fields(2,3,x1,y1,z1)
  
  euij(1,3) = &
    euij(1,3) &
    -0.25*Fields(1,3,x0,y0,z0) &
    +0.25*Fields(1,3,x0,y0,z1) &
    -0.25*Fields(1,3,x0,y1,z0) &
    +0.25*Fields(1,3,x0,y1,z1) &
    -0.25*Fields(1,3,x1,y0,z0) &
    +0.25*Fields(1,3,x1,y0,z1) &
    -0.25*Fields(1,3,x1,y1,z0) &
    +0.25*Fields(1,3,x1,y1,z1) &
    -0.25*Fields(3,3,x0,y0,z0) &
    -0.25*Fields(3,3,x0,y0,z1) &
    -0.25*Fields(3,3,x0,y1,z0) &
    -0.25*Fields(3,3,x0,y1,z1) &
    +0.25*Fields(3,3,x1,y0,z0) &
    +0.25*Fields(3,3,x1,y0,z1) &
    +0.25*Fields(3,3,x1,y1,z0) &
    +0.25*Fields(3,3,x1,y1,z1)
  
  euij(2,2) = &
    euij(2,2) &
    -0.25*Fields(2,3,x0,y0,z0) &
    -0.25*Fields(2,3,x0,y0,z1) &
    +0.25*Fields(2,3,x0,y1,z0) &
    +0.25*Fields(2,3,x0,y1,z1) &
    -0.25*Fields(2,3,x1,y0,z0) &
    -0.25*Fields(2,3,x1,y0,z1) &
    +0.25*Fields(2,3,x1,y1,z0) &
    +0.25*Fields(2,3,x1,y1,z1)
  
  euij(2,3) = &
    euij(2,3) &
    -0.25*Fields(2,3,x0,y0,z0) &
    +0.25*Fields(2,3,x0,y0,z1) &
    -0.25*Fields(2,3,x0,y1,z0) &
    +0.25*Fields(2,3,x0,y1,z1) &
    -0.25*Fields(2,3,x1,y0,z0) &
    +0.25*Fields(2,3,x1,y0,z1) &
    -0.25*Fields(2,3,x1,y1,z0) &
    +0.25*Fields(2,3,x1,y1,z1) &
    -0.25*Fields(3,3,x0,y0,z0) &
    -0.25*Fields(3,3,x0,y0,z1) &
    +0.25*Fields(3,3,x0,y1,z0) &
    +0.25*Fields(3,3,x0,y1,z1) &
    -0.25*Fields(3,3,x1,y0,z0) &
    -0.25*Fields(3,3,x1,y0,z1) &
    +0.25*Fields(3,3,x1,y1,z0) &
    +0.25*Fields(3,3,x1,y1,z1)
  
  euij(3,3) = &
    euij(3,3) &
    -0.25*Fields(3,3,x0,y0,z0) &
    +0.25*Fields(3,3,x0,y0,z1) &
    -0.25*Fields(3,3,x0,y1,z0) &
    +0.25*Fields(3,3,x0,y1,z1) &
    -0.25*Fields(3,3,x1,y0,z0) &
    +0.25*Fields(3,3,x1,y0,z1) &
    -0.25*Fields(3,3,x1,y1,z0) &
    +0.25*Fields(3,3,x1,y1,z1)

  euij(2,1) = euij(1,2)
  euij(3,1) = euij(1,3)
  euij(3,2) = euij(2,3)
  
End Function StrainFromu
