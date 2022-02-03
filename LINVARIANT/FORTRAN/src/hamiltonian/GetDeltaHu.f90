Function GetDeltaHu(x0, y0, z0, Fields, e0ij, idelta, delta) Result(DeltaH)
  
  Implicit none
  Integer, Intent(in) :: x0, y0, z0, idelta
  Real*8,  Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(in) :: delta(FieldDimList(idelta))
  Real*8,  Intent(in) :: e0ij(3,3)
  Real*8              :: eij(3,3), euij(3,3)
  Real*8              :: DeltaH
  Real*8              :: DeltaHu
  
  Integer             :: x1,y1,z1,xi1,yi1,zi1
  
  x1 = (x0+1)-floor(real(x0+1-1)/real(cgrid%n1))*cgrid%n1
  y1 = (y0+1)-floor(real(y0+1-1)/real(cgrid%n2))*cgrid%n2
  z1 = (z0+1)-floor(real(z0+1-1)/real(cgrid%n3))*cgrid%n3
  xi1 = (x0-1)-floor(real(x0-1-1)/real(cgrid%n1))*cgrid%n1
  yi1 = (y0-1)-floor(real(y0-1-1)/real(cgrid%n2))*cgrid%n2
  zi1 = (z0-1)-floor(real(z0-1-1)/real(cgrid%n3))*cgrid%n3
  
  DeltaH = 0.0D0
  DeltaHu = 0.0D0
  euij = StrainFromu(x0, y0, z0, Fields)
  eij = e0ij + euij
  
  SELECT CASE (idelta)
  CASE (1)
    DeltaHu = &
      DeltaHu &
      +0.
    
    DeltaH=DeltaHu
  
  CASE (2)
    DeltaHu = &
      DeltaHu &
      +0.
    
    DeltaH=DeltaHu
  
  CASE (3)
    DeltaHu = &
      DeltaHu &
      +2.*Coeffu(2)*delta(1)**2 &
      +Coeffu(4)*delta(1)**2 &
      +2.*Coeffu(2)*delta(2)**2 &
      +Coeffu(4)*delta(2)**2 &
      +2.*Coeffu(2)*delta(3)**2 &
      +Coeffu(4)*delta(3)**2 &
      +4.*Coeffu(2)*delta(1)*Fields(1,3,x0,y0,z0) &
      +2.*Coeffu(4)*delta(1)*Fields(1,3,x0,y0,z0) &
      -(Coeffu(2)*delta(1)*Fields(1,3,x0,y0,z1)) &
      -(Coeffu(2)*delta(1)*Fields(1,3,x0,y0,zi1)) &
      -(Coeffu(2)*delta(1)*Fields(1,3,x0,y1,z0)) &
      -(Coeffu(2)*delta(1)*Fields(1,3,x0,yi1,z0)) &
      -(Coeffu(4)*delta(1)*Fields(1,3,x1,y0,z0)) &
      -0.25*Coeffu(2)*delta(3)*Fields(1,3,x1,y0,z1) &
      -0.25*Coeffu(3)*delta(3)*Fields(1,3,x1,y0,z1) &
      +0.25*Coeffu(2)*delta(3)*Fields(1,3,x1,y0,zi1) &
      +0.25*Coeffu(3)*delta(3)*Fields(1,3,x1,y0,zi1) &
      -0.25*Coeffu(2)*delta(2)*Fields(1,3,x1,y1,z0) &
      -0.25*Coeffu(3)*delta(2)*Fields(1,3,x1,y1,z0) &
      +0.25*Coeffu(2)*delta(2)*Fields(1,3,x1,yi1,z0) &
      +0.25*Coeffu(3)*delta(2)*Fields(1,3,x1,yi1,z0) &
      -(Coeffu(4)*delta(1)*Fields(1,3,xi1,y0,z0)) &
      +0.25*Coeffu(2)*delta(3)*Fields(1,3,xi1,y0,z1) &
      +0.25*Coeffu(3)*delta(3)*Fields(1,3,xi1,y0,z1) &
      -0.25*Coeffu(2)*delta(3)*Fields(1,3,xi1,y0,zi1) &
      -0.25*Coeffu(3)*delta(3)*Fields(1,3,xi1,y0,zi1) &
      +0.25*Coeffu(2)*delta(2)*Fields(1,3,xi1,y1,z0) &
      +0.25*Coeffu(3)*delta(2)*Fields(1,3,xi1,y1,z0) &
      -0.25*Coeffu(2)*delta(2)*Fields(1,3,xi1,yi1,z0) &
      -0.25*Coeffu(3)*delta(2)*Fields(1,3,xi1,yi1,z0) &
      +4.*Coeffu(2)*delta(2)*Fields(2,3,x0,y0,z0) &
      +2.*Coeffu(4)*delta(2)*Fields(2,3,x0,y0,z0) &
      -(Coeffu(2)*delta(2)*Fields(2,3,x0,y0,z1)) &
      -(Coeffu(2)*delta(2)*Fields(2,3,x0,y0,zi1)) &
      -(Coeffu(4)*delta(2)*Fields(2,3,x0,y1,z0)) &
      -0.25*Coeffu(2)*delta(3)*Fields(2,3,x0,y1,z1) &
      -0.25*Coeffu(3)*delta(3)*Fields(2,3,x0,y1,z1) &
      +0.25*Coeffu(2)*delta(3)*Fields(2,3,x0,y1,zi1) &
      +0.25*Coeffu(3)*delta(3)*Fields(2,3,x0,y1,zi1) &
      -(Coeffu(4)*delta(2)*Fields(2,3,x0,yi1,z0)) &
      +0.25*Coeffu(2)*delta(3)*Fields(2,3,x0,yi1,z1) &
      +0.25*Coeffu(3)*delta(3)*Fields(2,3,x0,yi1,z1) &
      -0.25*Coeffu(2)*delta(3)*Fields(2,3,x0,yi1,zi1) &
      -0.25*Coeffu(3)*delta(3)*Fields(2,3,x0,yi1,zi1) &
      -(Coeffu(2)*delta(2)*Fields(2,3,x1,y0,z0)) &
      -0.25*Coeffu(2)*delta(1)*Fields(2,3,x1,y1,z0) &
      -0.25*Coeffu(3)*delta(1)*Fields(2,3,x1,y1,z0) &
      +0.25*Coeffu(2)*delta(1)*Fields(2,3,x1,yi1,z0) &
      +0.25*Coeffu(3)*delta(1)*Fields(2,3,x1,yi1,z0) &
      -(Coeffu(2)*delta(2)*Fields(2,3,xi1,y0,z0)) &
      +0.25*Coeffu(2)*delta(1)*Fields(2,3,xi1,y1,z0) &
      +0.25*Coeffu(3)*delta(1)*Fields(2,3,xi1,y1,z0) &
      -0.25*Coeffu(2)*delta(1)*Fields(2,3,xi1,yi1,z0) &
      -0.25*Coeffu(3)*delta(1)*Fields(2,3,xi1,yi1,z0) &
      +4.*Coeffu(2)*delta(3)*Fields(3,3,x0,y0,z0) &
      +2.*Coeffu(4)*delta(3)*Fields(3,3,x0,y0,z0) &
      -(Coeffu(4)*delta(3)*Fields(3,3,x0,y0,z1)) &
      -(Coeffu(4)*delta(3)*Fields(3,3,x0,y0,zi1)) &
      -(Coeffu(2)*delta(3)*Fields(3,3,x0,y1,z0)) &
      -0.25*Coeffu(2)*delta(2)*Fields(3,3,x0,y1,z1) &
      -0.25*Coeffu(3)*delta(2)*Fields(3,3,x0,y1,z1) &
      +0.25*Coeffu(2)*delta(2)*Fields(3,3,x0,y1,zi1) &
      +0.25*Coeffu(3)*delta(2)*Fields(3,3,x0,y1,zi1) &
      -(Coeffu(2)*delta(3)*Fields(3,3,x0,yi1,z0)) &
      +0.25*Coeffu(2)*delta(2)*Fields(3,3,x0,yi1,z1) &
      +0.25*Coeffu(3)*delta(2)*Fields(3,3,x0,yi1,z1) &
      -0.25*Coeffu(2)*delta(2)*Fields(3,3,x0,yi1,zi1) &
      -0.25*Coeffu(3)*delta(2)*Fields(3,3,x0,yi1,zi1) &
      -(Coeffu(2)*delta(3)*Fields(3,3,x1,y0,z0)) &
      -0.25*Coeffu(2)*delta(1)*Fields(3,3,x1,y0,z1) &
      -0.25*Coeffu(3)*delta(1)*Fields(3,3,x1,y0,z1) &
      +0.25*Coeffu(2)*delta(1)*Fields(3,3,x1,y0,zi1) &
      +0.25*Coeffu(3)*delta(1)*Fields(3,3,x1,y0,zi1) &
      -(Coeffu(2)*delta(3)*Fields(3,3,xi1,y0,z0)) &
      +0.25*Coeffu(2)*delta(1)*Fields(3,3,xi1,y0,z1) &
      +0.25*Coeffu(3)*delta(1)*Fields(3,3,xi1,y0,z1) &
      -0.25*Coeffu(2)*delta(1)*Fields(3,3,xi1,y0,zi1) &
      -0.25*Coeffu(3)*delta(1)*Fields(3,3,xi1,y0,zi1)
    
    DeltaH=DeltaHu
  
  CASE (4)
    DeltaHu = &
      DeltaHu &
      +0.
    
    DeltaH=DeltaHu
  
  CASE DEFAULT
      write(*,*) "mode out of range!"
      call abort
  END SELECT
  
End Function GetDeltaHu
