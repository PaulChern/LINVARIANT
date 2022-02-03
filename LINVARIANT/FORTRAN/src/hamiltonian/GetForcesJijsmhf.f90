Function GetForcesJijsmhf(Fields, e0ij) Result(ForcesJijsmhf)
  
  Implicit none
  Real*8,  Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(in) :: e0ij(3,3)
  Integer             :: x0, y0, z0
  Real*8              :: eij(3,3), euij(3,3)
  Real*8              :: ForcesJijsmhf(Max(FieldDim, 6), NumField+1, cgrid%n1, cgrid%n2, cgrid%n3)
  
  Integer             :: x1,y1,z1,xi1,yi1,zi1

  ForcesJijsmhf = 0.0D0
  !$OMP    PARALLEL DEFAULT(SHARED) PRIVATE(eij,euij,x0,y0,z0,x1,y1,z1,xi1,yi1,zi1)
  !$OMP    DO COLLAPSE(3)
  do z0 = 1, cgrid%n3
  do y0 = 1, cgrid%n2
  do x0 = 1, cgrid%n1
    
    euij = StrainFromu(x0, y0, z0, Fields)
    eij = e0ij + euij
    
    x1 = (x0+1)-floor(real(x0+1-1)/real(cgrid%n1))*cgrid%n1
    y1 = (y0+1)-floor(real(y0+1-1)/real(cgrid%n2))*cgrid%n2
    z1 = (z0+1)-floor(real(z0+1-1)/real(cgrid%n3))*cgrid%n3
    xi1 = (x0-1)-floor(real(x0-1-1)/real(cgrid%n1))*cgrid%n1
    yi1 = (y0-1)-floor(real(y0-1-1)/real(cgrid%n2))*cgrid%n2
    zi1 = (z0-1)-floor(real(z0-1-1)/real(cgrid%n3))*cgrid%n3
  
    ForcesJijsmhf(1,1,x0,y0,z0) = &
      ForcesJijsmhf(1,1,x0,y0,z0) &
      +CoeffJijsmhf(2)*Fields(1,1,x0,y1,z0)*Fields(2,2,x0,y0,z0) &
      -(CoeffJijsmhf(2)*Fields(1,1,x0,yi1,z0)*Fields(2,2,x0,y0,z0)) &
      +CoeffJijsmhf(1)*Fields(2,1,x1,y0,z0)*Fields(2,2,x0,y0,z0) &
      -(CoeffJijsmhf(1)*Fields(2,1,xi1,y0,z0)*Fields(2,2,x0,y0,z0)) &
      +CoeffJijsmhf(2)*Fields(1,1,x0,y0,z1)*Fields(3,2,x0,y0,z0) &
      -(CoeffJijsmhf(2)*Fields(1,1,x0,y0,zi1)*Fields(3,2,x0,y0,z0)) &
      +CoeffJijsmhf(1)*Fields(3,1,x1,y0,z0)*Fields(3,2,x0,y0,z0) &
      -(CoeffJijsmhf(1)*Fields(3,1,xi1,y0,z0)*Fields(3,2,x0,y0,z0))
    
    ForcesJijsmhf(2,1,x0,y0,z0) = &
      ForcesJijsmhf(2,1,x0,y0,z0) &
      +CoeffJijsmhf(1)*Fields(1,1,x0,y1,z0)*Fields(1,2,x0,y0,z0) &
      -(CoeffJijsmhf(1)*Fields(1,1,x0,yi1,z0)*Fields(1,2,x0,y0,z0)) &
      +CoeffJijsmhf(2)*Fields(1,2,x0,y0,z0)*Fields(2,1,x1,y0,z0) &
      -(CoeffJijsmhf(2)*Fields(1,2,x0,y0,z0)*Fields(2,1,xi1,y0,z0)) &
      +CoeffJijsmhf(2)*Fields(2,1,x0,y0,z1)*Fields(3,2,x0,y0,z0) &
      -(CoeffJijsmhf(2)*Fields(2,1,x0,y0,zi1)*Fields(3,2,x0,y0,z0)) &
      +CoeffJijsmhf(1)*Fields(3,1,x0,y1,z0)*Fields(3,2,x0,y0,z0) &
      -(CoeffJijsmhf(1)*Fields(3,1,x0,yi1,z0)*Fields(3,2,x0,y0,z0))
    
    ForcesJijsmhf(3,1,x0,y0,z0) = &
      ForcesJijsmhf(3,1,x0,y0,z0) &
      +CoeffJijsmhf(1)*Fields(1,1,x0,y0,z1)*Fields(1,2,x0,y0,z0) &
      -(CoeffJijsmhf(1)*Fields(1,1,x0,y0,zi1)*Fields(1,2,x0,y0,z0)) &
      +CoeffJijsmhf(1)*Fields(2,1,x0,y0,z1)*Fields(2,2,x0,y0,z0) &
      -(CoeffJijsmhf(1)*Fields(2,1,x0,y0,zi1)*Fields(2,2,x0,y0,z0)) &
      +CoeffJijsmhf(2)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y1,z0) &
      -(CoeffJijsmhf(2)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,yi1,z0)) &
      +CoeffJijsmhf(2)*Fields(1,2,x0,y0,z0)*Fields(3,1,x1,y0,z0) &
      -(CoeffJijsmhf(2)*Fields(1,2,x0,y0,z0)*Fields(3,1,xi1,y0,z0))
    
    ForcesJijsmhf(1,2,x0,y0,z0) = &
      ForcesJijsmhf(1,2,x0,y0,z0) &
      +CoeffJijsmhf(1)*Fields(1,1,x0,y1,z0)*Fields(2,1,x0,y0,z0) &
      -(CoeffJijsmhf(1)*Fields(1,1,x0,yi1,z0)*Fields(2,1,x0,y0,z0)) &
      +CoeffJijsmhf(2)*Fields(2,1,x0,y0,z0)*Fields(2,1,x1,y0,z0) &
      -(CoeffJijsmhf(2)*Fields(2,1,x0,y0,z0)*Fields(2,1,xi1,y0,z0)) &
      +CoeffJijsmhf(1)*Fields(1,1,x0,y0,z1)*Fields(3,1,x0,y0,z0) &
      -(CoeffJijsmhf(1)*Fields(1,1,x0,y0,zi1)*Fields(3,1,x0,y0,z0)) &
      +CoeffJijsmhf(2)*Fields(3,1,x0,y0,z0)*Fields(3,1,x1,y0,z0) &
      -(CoeffJijsmhf(2)*Fields(3,1,x0,y0,z0)*Fields(3,1,xi1,y0,z0))
    
    ForcesJijsmhf(2,2,x0,y0,z0) = &
      ForcesJijsmhf(2,2,x0,y0,z0) &
      +CoeffJijsmhf(2)*Fields(1,1,x0,y0,z0)*Fields(1,1,x0,y1,z0) &
      -(CoeffJijsmhf(2)*Fields(1,1,x0,y0,z0)*Fields(1,1,x0,yi1,z0)) &
      +CoeffJijsmhf(1)*Fields(1,1,x0,y0,z0)*Fields(2,1,x1,y0,z0) &
      -(CoeffJijsmhf(1)*Fields(1,1,x0,y0,z0)*Fields(2,1,xi1,y0,z0)) &
      +CoeffJijsmhf(1)*Fields(2,1,x0,y0,z1)*Fields(3,1,x0,y0,z0) &
      -(CoeffJijsmhf(1)*Fields(2,1,x0,y0,zi1)*Fields(3,1,x0,y0,z0)) &
      +CoeffJijsmhf(2)*Fields(3,1,x0,y0,z0)*Fields(3,1,x0,y1,z0) &
      -(CoeffJijsmhf(2)*Fields(3,1,x0,y0,z0)*Fields(3,1,x0,yi1,z0))
    
    ForcesJijsmhf(3,2,x0,y0,z0) = &
      ForcesJijsmhf(3,2,x0,y0,z0) &
      +CoeffJijsmhf(2)*Fields(1,1,x0,y0,z0)*Fields(1,1,x0,y0,z1) &
      -(CoeffJijsmhf(2)*Fields(1,1,x0,y0,z0)*Fields(1,1,x0,y0,zi1)) &
      +CoeffJijsmhf(2)*Fields(2,1,x0,y0,z0)*Fields(2,1,x0,y0,z1) &
      -(CoeffJijsmhf(2)*Fields(2,1,x0,y0,z0)*Fields(2,1,x0,y0,zi1)) &
      +CoeffJijsmhf(1)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,y1,z0) &
      -(CoeffJijsmhf(1)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,yi1,z0)) &
      +CoeffJijsmhf(1)*Fields(1,1,x0,y0,z0)*Fields(3,1,x1,y0,z0) &
      -(CoeffJijsmhf(1)*Fields(1,1,x0,y0,z0)*Fields(3,1,xi1,y0,z0))
    
    ForcesJijsmhf(1,3,x0,y0,z0) = &
      ForcesJijsmhf(1,3,x0,y0,z0) &
      +0.
    
    ForcesJijsmhf(2,3,x0,y0,z0) = &
      ForcesJijsmhf(2,3,x0,y0,z0) &
      +0.
    
    ForcesJijsmhf(3,3,x0,y0,z0) = &
      ForcesJijsmhf(3,3,x0,y0,z0) &
      +0.
    
    ForcesJijsmhf(1,4,x0,y0,z0) = &
      ForcesJijsmhf(1,4,x0,y0,z0) &
      +0.
    
    ForcesJijsmhf(2,4,x0,y0,z0) = &
      ForcesJijsmhf(2,4,x0,y0,z0) &
      +0.
    
    ForcesJijsmhf(3,4,x0,y0,z0) = &
      ForcesJijsmhf(3,4,x0,y0,z0) &
      +0.
    
    ForcesJijsmhf(4,4,x0,y0,z0) = &
      ForcesJijsmhf(4,4,x0,y0,z0) &
      +0.
    
    ForcesJijsmhf(5,4,x0,y0,z0) = &
      ForcesJijsmhf(5,4,x0,y0,z0) &
      +0.
    
    ForcesJijsmhf(6,4,x0,y0,z0) = &
      ForcesJijsmhf(6,4,x0,y0,z0) &
      +0.
    
  end do
  end do
  end do
  !$OMP    END DO
  !$OMP    END PARALLEL
  
End Function GetForcesJijsmhf
