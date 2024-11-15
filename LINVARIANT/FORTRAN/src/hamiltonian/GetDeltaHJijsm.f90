Function GetDeltaHJijsm(x0, y0, z0, Fields, e0ij, idelta, delta) Result(DeltaH)
  
  Implicit none
  Integer, Intent(in) :: x0, y0, z0, idelta
  Real*8,  Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(in) :: delta(FieldDimList(idelta))
  Real*8,  Intent(in) :: e0ij(3,3)
  Real*8              :: eij(3,3), euij(3,3)
  Real*8              :: DeltaH
  Real*8              :: DeltaHJijsm
  
  Integer             :: x1,y1,z1,xi1,yi1,zi1
  
  x1 = (x0+1)-floor(real(x0+1-1)/real(cgrid%n1))*cgrid%n1
  y1 = (y0+1)-floor(real(y0+1-1)/real(cgrid%n2))*cgrid%n2
  z1 = (z0+1)-floor(real(z0+1-1)/real(cgrid%n3))*cgrid%n3
  xi1 = (x0-1)-floor(real(x0-1-1)/real(cgrid%n1))*cgrid%n1
  yi1 = (y0-1)-floor(real(y0-1-1)/real(cgrid%n2))*cgrid%n2
  zi1 = (z0-1)-floor(real(z0-1-1)/real(cgrid%n3))*cgrid%n3
  
  DeltaH = 0.0D0
  DeltaHJijsm = 0.0D0
  euij = StrainFromu(x0, y0, z0, Fields)
  eij = e0ij + euij
  
  SELECT CASE (idelta)
  CASE (1)
    DeltaHJijsm = &
      DeltaHJijsm &
      +2.*CoeffJijsm(2)*delta(1)*Fields(1,1,x0,y0,z1) &
      +2.*CoeffJijsm(2)*delta(1)*Fields(1,1,x0,y0,zi1) &
      +2.*CoeffJijsm(2)*delta(1)*Fields(1,1,x0,y1,z0) &
      +2.*CoeffJijsm(4)*delta(1)*Fields(1,1,x0,y1,z1) &
      +2.*CoeffJijsm(4)*delta(1)*Fields(1,1,x0,y1,zi1) &
      +2.*CoeffJijsm(2)*delta(1)*Fields(1,1,x0,yi1,z0) &
      +2.*CoeffJijsm(4)*delta(1)*Fields(1,1,x0,yi1,z1) &
      +2.*CoeffJijsm(4)*delta(1)*Fields(1,1,x0,yi1,zi1) &
      +2.*CoeffJijsm(1)*delta(1)*Fields(1,1,x1,y0,z0) &
      +2.*CoeffJijsm(3)*delta(1)*Fields(1,1,x1,y0,z1) &
      +2.*CoeffJijsm(6)*delta(3)*Fields(1,1,x1,y0,z1) &
      +2.*CoeffJijsm(3)*delta(1)*Fields(1,1,x1,y0,zi1) &
      -2.*CoeffJijsm(6)*delta(3)*Fields(1,1,x1,y0,zi1) &
      +2.*CoeffJijsm(3)*delta(1)*Fields(1,1,x1,y1,z0) &
      +2.*CoeffJijsm(6)*delta(2)*Fields(1,1,x1,y1,z0) &
      +2.*CoeffJijsm(5)*delta(1)*Fields(1,1,x1,y1,z1) &
      +2.*CoeffJijsm(7)*delta(2)*Fields(1,1,x1,y1,z1) &
      +2.*CoeffJijsm(7)*delta(3)*Fields(1,1,x1,y1,z1) &
      +2.*CoeffJijsm(5)*delta(1)*Fields(1,1,x1,y1,zi1) &
      +2.*CoeffJijsm(7)*delta(2)*Fields(1,1,x1,y1,zi1) &
      -2.*CoeffJijsm(7)*delta(3)*Fields(1,1,x1,y1,zi1) &
      +2.*CoeffJijsm(3)*delta(1)*Fields(1,1,x1,yi1,z0) &
      -2.*CoeffJijsm(6)*delta(2)*Fields(1,1,x1,yi1,z0) &
      +2.*CoeffJijsm(5)*delta(1)*Fields(1,1,x1,yi1,z1) &
      -2.*CoeffJijsm(7)*delta(2)*Fields(1,1,x1,yi1,z1) &
      +2.*CoeffJijsm(7)*delta(3)*Fields(1,1,x1,yi1,z1) &
      +2.*CoeffJijsm(5)*delta(1)*Fields(1,1,x1,yi1,zi1) &
      -2.*CoeffJijsm(7)*delta(2)*Fields(1,1,x1,yi1,zi1) &
      -2.*CoeffJijsm(7)*delta(3)*Fields(1,1,x1,yi1,zi1) &
      +2.*CoeffJijsm(1)*delta(1)*Fields(1,1,xi1,y0,z0) &
      +2.*CoeffJijsm(3)*delta(1)*Fields(1,1,xi1,y0,z1) &
      -2.*CoeffJijsm(6)*delta(3)*Fields(1,1,xi1,y0,z1) &
      +2.*CoeffJijsm(3)*delta(1)*Fields(1,1,xi1,y0,zi1) &
      +2.*CoeffJijsm(6)*delta(3)*Fields(1,1,xi1,y0,zi1) &
      +2.*CoeffJijsm(3)*delta(1)*Fields(1,1,xi1,y1,z0) &
      -2.*CoeffJijsm(6)*delta(2)*Fields(1,1,xi1,y1,z0) &
      +2.*CoeffJijsm(5)*delta(1)*Fields(1,1,xi1,y1,z1) &
      -2.*CoeffJijsm(7)*delta(2)*Fields(1,1,xi1,y1,z1) &
      -2.*CoeffJijsm(7)*delta(3)*Fields(1,1,xi1,y1,z1) &
      +2.*CoeffJijsm(5)*delta(1)*Fields(1,1,xi1,y1,zi1) &
      -2.*CoeffJijsm(7)*delta(2)*Fields(1,1,xi1,y1,zi1) &
      +2.*CoeffJijsm(7)*delta(3)*Fields(1,1,xi1,y1,zi1) &
      +2.*CoeffJijsm(3)*delta(1)*Fields(1,1,xi1,yi1,z0) &
      +2.*CoeffJijsm(6)*delta(2)*Fields(1,1,xi1,yi1,z0) &
      +2.*CoeffJijsm(5)*delta(1)*Fields(1,1,xi1,yi1,z1) &
      +2.*CoeffJijsm(7)*delta(2)*Fields(1,1,xi1,yi1,z1) &
      -2.*CoeffJijsm(7)*delta(3)*Fields(1,1,xi1,yi1,z1) &
      +2.*CoeffJijsm(5)*delta(1)*Fields(1,1,xi1,yi1,zi1) &
      +2.*CoeffJijsm(7)*delta(2)*Fields(1,1,xi1,yi1,zi1) &
      +2.*CoeffJijsm(7)*delta(3)*Fields(1,1,xi1,yi1,zi1) &
      +2.*CoeffJijsm(2)*delta(2)*Fields(2,1,x0,y0,z1) &
      +2.*CoeffJijsm(2)*delta(2)*Fields(2,1,x0,y0,zi1) &
      +2.*CoeffJijsm(1)*delta(2)*Fields(2,1,x0,y1,z0) &
      +2.*CoeffJijsm(3)*delta(2)*Fields(2,1,x0,y1,z1) &
      +2.*CoeffJijsm(6)*delta(3)*Fields(2,1,x0,y1,z1) &
      +2.*CoeffJijsm(3)*delta(2)*Fields(2,1,x0,y1,zi1) &
      -2.*CoeffJijsm(6)*delta(3)*Fields(2,1,x0,y1,zi1) &
      +2.*CoeffJijsm(1)*delta(2)*Fields(2,1,x0,yi1,z0) &
      +2.*CoeffJijsm(3)*delta(2)*Fields(2,1,x0,yi1,z1) &
      -2.*CoeffJijsm(6)*delta(3)*Fields(2,1,x0,yi1,z1) &
      +2.*CoeffJijsm(3)*delta(2)*Fields(2,1,x0,yi1,zi1) &
      +2.*CoeffJijsm(6)*delta(3)*Fields(2,1,x0,yi1,zi1) &
      +2.*CoeffJijsm(2)*delta(2)*Fields(2,1,x1,y0,z0) &
      +2.*CoeffJijsm(4)*delta(2)*Fields(2,1,x1,y0,z1) &
      +2.*CoeffJijsm(4)*delta(2)*Fields(2,1,x1,y0,zi1) &
      +2.*CoeffJijsm(6)*delta(1)*Fields(2,1,x1,y1,z0) &
      +2.*CoeffJijsm(3)*delta(2)*Fields(2,1,x1,y1,z0) &
      +2.*CoeffJijsm(7)*delta(1)*Fields(2,1,x1,y1,z1) &
      +2.*CoeffJijsm(5)*delta(2)*Fields(2,1,x1,y1,z1) &
      +2.*CoeffJijsm(7)*delta(3)*Fields(2,1,x1,y1,z1) &
      +2.*CoeffJijsm(7)*delta(1)*Fields(2,1,x1,y1,zi1) &
      +2.*CoeffJijsm(5)*delta(2)*Fields(2,1,x1,y1,zi1) &
      -2.*CoeffJijsm(7)*delta(3)*Fields(2,1,x1,y1,zi1) &
      -2.*CoeffJijsm(6)*delta(1)*Fields(2,1,x1,yi1,z0) &
      +2.*CoeffJijsm(3)*delta(2)*Fields(2,1,x1,yi1,z0) &
      -2.*CoeffJijsm(7)*delta(1)*Fields(2,1,x1,yi1,z1) &
      +2.*CoeffJijsm(5)*delta(2)*Fields(2,1,x1,yi1,z1) &
      -2.*CoeffJijsm(7)*delta(3)*Fields(2,1,x1,yi1,z1) &
      -2.*CoeffJijsm(7)*delta(1)*Fields(2,1,x1,yi1,zi1) &
      +2.*CoeffJijsm(5)*delta(2)*Fields(2,1,x1,yi1,zi1) &
      +2.*CoeffJijsm(7)*delta(3)*Fields(2,1,x1,yi1,zi1) &
      +2.*CoeffJijsm(2)*delta(2)*Fields(2,1,xi1,y0,z0) &
      +2.*CoeffJijsm(4)*delta(2)*Fields(2,1,xi1,y0,z1) &
      +2.*CoeffJijsm(4)*delta(2)*Fields(2,1,xi1,y0,zi1) &
      -2.*CoeffJijsm(6)*delta(1)*Fields(2,1,xi1,y1,z0) &
      +2.*CoeffJijsm(3)*delta(2)*Fields(2,1,xi1,y1,z0) &
      -2.*CoeffJijsm(7)*delta(1)*Fields(2,1,xi1,y1,z1) &
      +2.*CoeffJijsm(5)*delta(2)*Fields(2,1,xi1,y1,z1) &
      +2.*CoeffJijsm(7)*delta(3)*Fields(2,1,xi1,y1,z1) &
      -2.*CoeffJijsm(7)*delta(1)*Fields(2,1,xi1,y1,zi1) &
      +2.*CoeffJijsm(5)*delta(2)*Fields(2,1,xi1,y1,zi1) &
      -2.*CoeffJijsm(7)*delta(3)*Fields(2,1,xi1,y1,zi1) &
      +2.*CoeffJijsm(6)*delta(1)*Fields(2,1,xi1,yi1,z0) &
      +2.*CoeffJijsm(3)*delta(2)*Fields(2,1,xi1,yi1,z0) &
      +2.*CoeffJijsm(7)*delta(1)*Fields(2,1,xi1,yi1,z1) &
      +2.*CoeffJijsm(5)*delta(2)*Fields(2,1,xi1,yi1,z1) &
      -2.*CoeffJijsm(7)*delta(3)*Fields(2,1,xi1,yi1,z1) &
      +2.*CoeffJijsm(7)*delta(1)*Fields(2,1,xi1,yi1,zi1) &
      +2.*CoeffJijsm(5)*delta(2)*Fields(2,1,xi1,yi1,zi1) &
      +2.*CoeffJijsm(7)*delta(3)*Fields(2,1,xi1,yi1,zi1) &
      +2.*CoeffJijsm(1)*delta(3)*Fields(3,1,x0,y0,z1) &
      +2.*CoeffJijsm(1)*delta(3)*Fields(3,1,x0,y0,zi1) &
      +2.*CoeffJijsm(2)*delta(3)*Fields(3,1,x0,y1,z0) &
      +2.*CoeffJijsm(6)*delta(2)*Fields(3,1,x0,y1,z1) &
      +2.*CoeffJijsm(3)*delta(3)*Fields(3,1,x0,y1,z1) &
      -2.*CoeffJijsm(6)*delta(2)*Fields(3,1,x0,y1,zi1) &
      +2.*CoeffJijsm(3)*delta(3)*Fields(3,1,x0,y1,zi1) &
      +2.*CoeffJijsm(2)*delta(3)*Fields(3,1,x0,yi1,z0) &
      -2.*CoeffJijsm(6)*delta(2)*Fields(3,1,x0,yi1,z1) &
      +2.*CoeffJijsm(3)*delta(3)*Fields(3,1,x0,yi1,z1) &
      +2.*CoeffJijsm(6)*delta(2)*Fields(3,1,x0,yi1,zi1) &
      +2.*CoeffJijsm(3)*delta(3)*Fields(3,1,x0,yi1,zi1) &
      +2.*CoeffJijsm(2)*delta(3)*Fields(3,1,x1,y0,z0) &
      +2.*CoeffJijsm(6)*delta(1)*Fields(3,1,x1,y0,z1) &
      +2.*CoeffJijsm(3)*delta(3)*Fields(3,1,x1,y0,z1) &
      -2.*CoeffJijsm(6)*delta(1)*Fields(3,1,x1,y0,zi1) &
      +2.*CoeffJijsm(3)*delta(3)*Fields(3,1,x1,y0,zi1) &
      +2.*CoeffJijsm(4)*delta(3)*Fields(3,1,x1,y1,z0) &
      +2.*CoeffJijsm(7)*delta(1)*Fields(3,1,x1,y1,z1) &
      +2.*CoeffJijsm(7)*delta(2)*Fields(3,1,x1,y1,z1) &
      +2.*CoeffJijsm(5)*delta(3)*Fields(3,1,x1,y1,z1) &
      -2.*CoeffJijsm(7)*delta(1)*Fields(3,1,x1,y1,zi1) &
      -2.*CoeffJijsm(7)*delta(2)*Fields(3,1,x1,y1,zi1) &
      +2.*CoeffJijsm(5)*delta(3)*Fields(3,1,x1,y1,zi1) &
      +2.*CoeffJijsm(4)*delta(3)*Fields(3,1,x1,yi1,z0) &
      +2.*CoeffJijsm(7)*delta(1)*Fields(3,1,x1,yi1,z1) &
      -2.*CoeffJijsm(7)*delta(2)*Fields(3,1,x1,yi1,z1) &
      +2.*CoeffJijsm(5)*delta(3)*Fields(3,1,x1,yi1,z1) &
      -2.*CoeffJijsm(7)*delta(1)*Fields(3,1,x1,yi1,zi1) &
      +2.*CoeffJijsm(7)*delta(2)*Fields(3,1,x1,yi1,zi1) &
      +2.*CoeffJijsm(5)*delta(3)*Fields(3,1,x1,yi1,zi1) &
      +2.*CoeffJijsm(2)*delta(3)*Fields(3,1,xi1,y0,z0) &
      -2.*CoeffJijsm(6)*delta(1)*Fields(3,1,xi1,y0,z1) &
      +2.*CoeffJijsm(3)*delta(3)*Fields(3,1,xi1,y0,z1) &
      +2.*CoeffJijsm(6)*delta(1)*Fields(3,1,xi1,y0,zi1) &
      +2.*CoeffJijsm(3)*delta(3)*Fields(3,1,xi1,y0,zi1) &
      +2.*CoeffJijsm(4)*delta(3)*Fields(3,1,xi1,y1,z0) &
      -2.*CoeffJijsm(7)*delta(1)*Fields(3,1,xi1,y1,z1) &
      +2.*CoeffJijsm(7)*delta(2)*Fields(3,1,xi1,y1,z1) &
      +2.*CoeffJijsm(5)*delta(3)*Fields(3,1,xi1,y1,z1) &
      +2.*CoeffJijsm(7)*delta(1)*Fields(3,1,xi1,y1,zi1) &
      -2.*CoeffJijsm(7)*delta(2)*Fields(3,1,xi1,y1,zi1) &
      +2.*CoeffJijsm(5)*delta(3)*Fields(3,1,xi1,y1,zi1) &
      +2.*CoeffJijsm(4)*delta(3)*Fields(3,1,xi1,yi1,z0) &
      -2.*CoeffJijsm(7)*delta(1)*Fields(3,1,xi1,yi1,z1) &
      -2.*CoeffJijsm(7)*delta(2)*Fields(3,1,xi1,yi1,z1) &
      +2.*CoeffJijsm(5)*delta(3)*Fields(3,1,xi1,yi1,z1) &
      +2.*CoeffJijsm(7)*delta(1)*Fields(3,1,xi1,yi1,zi1) &
      +2.*CoeffJijsm(7)*delta(2)*Fields(3,1,xi1,yi1,zi1) &
      +2.*CoeffJijsm(5)*delta(3)*Fields(3,1,xi1,yi1,zi1)
    
    DeltaH=DeltaHJijsm
  
  CASE (2)
    DeltaHJijsm = &
      DeltaHJijsm &
      +0.
    
    DeltaH=DeltaHJijsm
  
  CASE (3)
    DeltaHJijsm = &
      DeltaHJijsm &
      +0.
    
    DeltaH=DeltaHJijsm
  
  CASE (4)
    DeltaHJijsm = &
      DeltaHJijsm &
      +0.
    
    DeltaH=DeltaHJijsm
  
  CASE DEFAULT
      write(*,*) "mode out of range!"
      call abort
  END SELECT
  
End Function GetDeltaHJijsm
