Function GetDeltaHJijsmhf(x0, y0, z0, Fields, e0ij, idelta, delta) Result(DeltaH)
  
  Implicit none
  Integer, Intent(in) :: x0, y0, z0, idelta
  Real*8,  Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(in) :: delta(FieldDimList(idelta))
  Real*8,  Intent(in) :: e0ij(3,3)
  Real*8              :: eij(3,3), euij(3,3)
  Real*8              :: DeltaH
  Real*8              :: DeltaHJijsmhf
  
  Integer             :: x1,y1,z1,xi1,yi1,zi1

  x1 = (x0+1)-floor(real(x0+1-1)/real(cgrid%n1))*cgrid%n1
  y1 = (y0+1)-floor(real(y0+1-1)/real(cgrid%n2))*cgrid%n2
  z1 = (z0+1)-floor(real(z0+1-1)/real(cgrid%n3))*cgrid%n3
  xi1 = (x0-1)-floor(real(x0-1-1)/real(cgrid%n1))*cgrid%n1
  yi1 = (y0-1)-floor(real(y0-1-1)/real(cgrid%n2))*cgrid%n2
  zi1 = (z0-1)-floor(real(z0-1-1)/real(cgrid%n3))*cgrid%n3


  DeltaH = 0.0D0
  DeltaHJijsmhf = 0.0D0
  euij = StrainFromu(x0, y0, z0, Fields)
  eij = e0ij + euij
  
  SELECT CASE (idelta)
  CASE (1)
    DeltaHJijsmhf = &
      DeltaHJijsmhf &
      -(CoeffJijsmhf(1)*delta(3)*Fields(1,1,x0,y0,z1)*Fields(1,2,x0,y0,z0)) &
      +CoeffJijsmhf(1)*delta(3)*Fields(1,1,x0,y0,zi1)*Fields(1,2,x0,y0,z0) &
      -(CoeffJijsmhf(1)*delta(2)*Fields(1,1,x0,y1,z0)*Fields(1,2,x0,y0,z0)) &
      +CoeffJijsmhf(1)*delta(2)*Fields(1,1,x0,yi1,z0)*Fields(1,2,x0,y0,z0) &
      -(CoeffJijsmhf(2)*delta(2)*Fields(1,2,x0,y0,z0)*Fields(2,1,x1,y0,z0)) &
      +CoeffJijsmhf(2)*delta(2)*Fields(1,2,x0,y0,z0)*Fields(2,1,xi1,y0,z0) &
      -(CoeffJijsmhf(2)*delta(1)*Fields(1,1,x0,y1,z0)*Fields(2,2,x0,y0,z0)) &
      +CoeffJijsmhf(2)*delta(1)*Fields(1,1,x0,yi1,z0)*Fields(2,2,x0,y0,z0) &
      -(CoeffJijsmhf(1)*delta(3)*Fields(2,1,x0,y0,z1)*Fields(2,2,x0,y0,z0)) &
      +CoeffJijsmhf(1)*delta(3)*Fields(2,1,x0,y0,zi1)*Fields(2,2,x0,y0,z0) &
      -(CoeffJijsmhf(1)*delta(1)*Fields(2,1,x1,y0,z0)*Fields(2,2,x0,y0,z0)) &
      +CoeffJijsmhf(1)*delta(1)*Fields(2,1,xi1,y0,z0)*Fields(2,2,x0,y0,z0) &
      -(CoeffJijsmhf(2)*delta(3)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y1,z0)) &
      +CoeffJijsmhf(2)*delta(3)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,yi1,z0) &
      -(CoeffJijsmhf(2)*delta(3)*Fields(1,2,x0,y0,z0)*Fields(3,1,x1,y0,z0)) &
      +CoeffJijsmhf(2)*delta(3)*Fields(1,2,x0,y0,z0)*Fields(3,1,xi1,y0,z0) &
      -(CoeffJijsmhf(2)*delta(1)*Fields(1,1,x0,y0,z1)*Fields(3,2,x0,y0,z0)) &
      +CoeffJijsmhf(2)*delta(1)*Fields(1,1,x0,y0,zi1)*Fields(3,2,x0,y0,z0) &
      -(CoeffJijsmhf(2)*delta(2)*Fields(2,1,x0,y0,z1)*Fields(3,2,x0,y0,z0)) &
      +CoeffJijsmhf(2)*delta(2)*Fields(2,1,x0,y0,zi1)*Fields(3,2,x0,y0,z0) &
      -(CoeffJijsmhf(1)*delta(2)*Fields(3,1,x0,y1,z0)*Fields(3,2,x0,y0,z0)) &
      +CoeffJijsmhf(1)*delta(2)*Fields(3,1,x0,yi1,z0)*Fields(3,2,x0,y0,z0) &
      -(CoeffJijsmhf(1)*delta(1)*Fields(3,1,x1,y0,z0)*Fields(3,2,x0,y0,z0)) &
      +CoeffJijsmhf(1)*delta(1)*Fields(3,1,xi1,y0,z0)*Fields(3,2,x0,y0,z0)
    
    DeltaH=DeltaHJijsmhf
  
  CASE (2)
    DeltaHJijsmhf = &
      DeltaHJijsmhf &
      -(CoeffJijsmhf(2)*delta(3)*Fields(1,1,x0,y0,z0)*Fields(1,1,x0,y0,z1)) &
      +CoeffJijsmhf(2)*delta(3)*Fields(1,1,x0,y0,z0)*Fields(1,1,x0,y0,zi1) &
      -(CoeffJijsmhf(2)*delta(2)*Fields(1,1,x0,y0,z0)*Fields(1,1,x0,y1,z0)) &
      +CoeffJijsmhf(2)*delta(2)*Fields(1,1,x0,y0,z0)*Fields(1,1,x0,yi1,z0) &
      -(CoeffJijsmhf(1)*delta(1)*Fields(1,1,x0,y1,z0)*Fields(2,1,x0,y0,z0)) &
      +CoeffJijsmhf(1)*delta(1)*Fields(1,1,x0,yi1,z0)*Fields(2,1,x0,y0,z0) &
      -(CoeffJijsmhf(2)*delta(3)*Fields(2,1,x0,y0,z0)*Fields(2,1,x0,y0,z1)) &
      +CoeffJijsmhf(2)*delta(3)*Fields(2,1,x0,y0,z0)*Fields(2,1,x0,y0,zi1) &
      -(CoeffJijsmhf(1)*delta(2)*Fields(1,1,x0,y0,z0)*Fields(2,1,x1,y0,z0)) &
      -(CoeffJijsmhf(2)*delta(1)*Fields(2,1,x0,y0,z0)*Fields(2,1,x1,y0,z0)) &
      +CoeffJijsmhf(1)*delta(2)*Fields(1,1,x0,y0,z0)*Fields(2,1,xi1,y0,z0) &
      +CoeffJijsmhf(2)*delta(1)*Fields(2,1,x0,y0,z0)*Fields(2,1,xi1,y0,z0) &
      -(CoeffJijsmhf(1)*delta(1)*Fields(1,1,x0,y0,z1)*Fields(3,1,x0,y0,z0)) &
      +CoeffJijsmhf(1)*delta(1)*Fields(1,1,x0,y0,zi1)*Fields(3,1,x0,y0,z0) &
      -(CoeffJijsmhf(1)*delta(2)*Fields(2,1,x0,y0,z1)*Fields(3,1,x0,y0,z0)) &
      +CoeffJijsmhf(1)*delta(2)*Fields(2,1,x0,y0,zi1)*Fields(3,1,x0,y0,z0) &
      -(CoeffJijsmhf(1)*delta(3)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,y1,z0)) &
      -(CoeffJijsmhf(2)*delta(2)*Fields(3,1,x0,y0,z0)*Fields(3,1,x0,y1,z0)) &
      +CoeffJijsmhf(1)*delta(3)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,yi1,z0) &
      +CoeffJijsmhf(2)*delta(2)*Fields(3,1,x0,y0,z0)*Fields(3,1,x0,yi1,z0) &
      -(CoeffJijsmhf(1)*delta(3)*Fields(1,1,x0,y0,z0)*Fields(3,1,x1,y0,z0)) &
      -(CoeffJijsmhf(2)*delta(1)*Fields(3,1,x0,y0,z0)*Fields(3,1,x1,y0,z0)) &
      +CoeffJijsmhf(1)*delta(3)*Fields(1,1,x0,y0,z0)*Fields(3,1,xi1,y0,z0) &
      +CoeffJijsmhf(2)*delta(1)*Fields(3,1,x0,y0,z0)*Fields(3,1,xi1,y0,z0)
    
    DeltaH=DeltaHJijsmhf
  
  CASE (3)
    DeltaHJijsmhf = &
      DeltaHJijsmhf &
      +0.
    
    DeltaH=DeltaHJijsmhf
  
  CASE (4)
    DeltaHJijsmhf = &
      DeltaHJijsmhf &
      +0.
    
    DeltaH=DeltaHJijsmhf
  
  CASE DEFAULT
      write(*,*) "mode out of range!"
      call abort
  END SELECT
  
End Function GetDeltaHJijsmhf
