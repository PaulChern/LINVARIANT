Function GetForces(Fields, e0ij) Result(Forces)
  
  Implicit none
  Real*8,  Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(in) :: e0ij(3,3)
  Real*8              :: eij(3,3), euij(3,3)
  Real*8              :: Forces(Max(FieldDim, 6), NumField+1, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8              :: ForcesDisp(Max(FieldDim, 6), NumField+1, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8              :: ForcesJijsm(Max(FieldDim, 6), NumField+1, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8              :: ForcesJijhf(Max(FieldDim, 6), NumField+1, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8              :: Forcesu(Max(FieldDim, 6), NumField+1, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8              :: ForcesDispu(Max(FieldDim, 6), NumField+1, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8              :: ForcesEpsDisp(Max(FieldDim, 6), NumField+1, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8              :: ForcesEpsu(Max(FieldDim, 6), NumField+1, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8              :: ForcesEps(Max(FieldDim, 6), NumField+1, cgrid%n1, cgrid%n2, cgrid%n3)
  
  Integer             :: x0, y0, z0
  Integer             :: x1,y1,z1,xi1,yi1,zi1
  
    
  Forces = 0.0D0
  ForcesDisp = 0.0D0
  ForcesJijsm = 0.0D0
  ForcesJijhf = 0.0D0
  Forcesu = 0.0D0
  ForcesDispu = 0.0D0
  ForcesEpsDisp = 0.0D0
  ForcesEpsu = 0.0D0
  ForcesEps = 0.0D0

!  !$OMP    PARALLEL  DEFAULT(SHARED) PRIVATE(x1,y1,z1,xi1,yi1,zi1,euij,eij)
!  !$OMP    DO COLLAPSE(1)
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
  
  ForcesDisp(1,1,x0,y0,z0) = &
    ForcesDisp(1,1,x0,y0,z0) &
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
  
  ForcesDisp(2,1,x0,y0,z0) = &
    ForcesDisp(2,1,x0,y0,z0) &
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
  
  ForcesDisp(3,1,x0,y0,z0) = &
    ForcesDisp(3,1,x0,y0,z0) &
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
  
  ForcesDisp(1,2,x0,y0,z0) = &
    ForcesDisp(1,2,x0,y0,z0) &
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
  
  ForcesDisp(2,2,x0,y0,z0) = &
    ForcesDisp(2,2,x0,y0,z0) &
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
  
  ForcesDisp(3,2,x0,y0,z0) = &
    ForcesDisp(3,2,x0,y0,z0) &
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
  
  ForcesDisp(1,3,x0,y0,z0) = &
    ForcesDisp(1,3,x0,y0,z0) &
    +0.
  
  ForcesDisp(2,3,x0,y0,z0) = &
    ForcesDisp(2,3,x0,y0,z0) &
    +0.
  
  ForcesDisp(3,3,x0,y0,z0) = &
    ForcesDisp(3,3,x0,y0,z0) &
    +0.
  
  ForcesDisp(1,4,x0,y0,z0) = &
    ForcesDisp(1,4,x0,y0,z0) &
    +0.
  
  ForcesDisp(2,4,x0,y0,z0) = &
    ForcesDisp(2,4,x0,y0,z0) &
    +0.
  
  ForcesDisp(3,4,x0,y0,z0) = &
    ForcesDisp(3,4,x0,y0,z0) &
    +0.
  
  ForcesDisp(4,4,x0,y0,z0) = &
    ForcesDisp(4,4,x0,y0,z0) &
    +0.
  
  ForcesDisp(5,4,x0,y0,z0) = &
    ForcesDisp(5,4,x0,y0,z0) &
    +0.
  
  ForcesDisp(6,4,x0,y0,z0) = &
    ForcesDisp(6,4,x0,y0,z0) &
    +0.
  
  ForcesJijsm(1,1,x0,y0,z0) = &
    ForcesJijsm(1,1,x0,y0,z0) &
    -2.*CoeffJijsm(2)*Fields(1,1,x0,y0,z1) &
    -2.*CoeffJijsm(2)*Fields(1,1,x0,y0,zi1) &
    -2.*CoeffJijsm(2)*Fields(1,1,x0,y1,z0) &
    -2.*CoeffJijsm(4)*Fields(1,1,x0,y1,z1) &
    -2.*CoeffJijsm(4)*Fields(1,1,x0,y1,zi1) &
    -2.*CoeffJijsm(2)*Fields(1,1,x0,yi1,z0) &
    -2.*CoeffJijsm(4)*Fields(1,1,x0,yi1,z1) &
    -2.*CoeffJijsm(4)*Fields(1,1,x0,yi1,zi1) &
    -2.*CoeffJijsm(1)*Fields(1,1,x1,y0,z0) &
    -2.*CoeffJijsm(3)*Fields(1,1,x1,y0,z1) &
    -2.*CoeffJijsm(3)*Fields(1,1,x1,y0,zi1) &
    -2.*CoeffJijsm(3)*Fields(1,1,x1,y1,z0) &
    -2.*CoeffJijsm(5)*Fields(1,1,x1,y1,z1) &
    -2.*CoeffJijsm(5)*Fields(1,1,x1,y1,zi1) &
    -2.*CoeffJijsm(3)*Fields(1,1,x1,yi1,z0) &
    -2.*CoeffJijsm(5)*Fields(1,1,x1,yi1,z1) &
    -2.*CoeffJijsm(5)*Fields(1,1,x1,yi1,zi1) &
    -2.*CoeffJijsm(1)*Fields(1,1,xi1,y0,z0) &
    -2.*CoeffJijsm(3)*Fields(1,1,xi1,y0,z1) &
    -2.*CoeffJijsm(3)*Fields(1,1,xi1,y0,zi1) &
    -2.*CoeffJijsm(3)*Fields(1,1,xi1,y1,z0) &
    -2.*CoeffJijsm(5)*Fields(1,1,xi1,y1,z1) &
    -2.*CoeffJijsm(5)*Fields(1,1,xi1,y1,zi1) &
    -2.*CoeffJijsm(3)*Fields(1,1,xi1,yi1,z0) &
    -2.*CoeffJijsm(5)*Fields(1,1,xi1,yi1,z1) &
    -2.*CoeffJijsm(5)*Fields(1,1,xi1,yi1,zi1) &
    -2.*CoeffJijsm(6)*Fields(2,1,x1,y1,z0) &
    -2.*CoeffJijsm(7)*Fields(2,1,x1,y1,z1) &
    -2.*CoeffJijsm(7)*Fields(2,1,x1,y1,zi1) &
    +2.*CoeffJijsm(6)*Fields(2,1,x1,yi1,z0) &
    +2.*CoeffJijsm(7)*Fields(2,1,x1,yi1,z1) &
    +2.*CoeffJijsm(7)*Fields(2,1,x1,yi1,zi1) &
    +2.*CoeffJijsm(6)*Fields(2,1,xi1,y1,z0) &
    +2.*CoeffJijsm(7)*Fields(2,1,xi1,y1,z1) &
    +2.*CoeffJijsm(7)*Fields(2,1,xi1,y1,zi1) &
    -2.*CoeffJijsm(6)*Fields(2,1,xi1,yi1,z0) &
    -2.*CoeffJijsm(7)*Fields(2,1,xi1,yi1,z1) &
    -2.*CoeffJijsm(7)*Fields(2,1,xi1,yi1,zi1) &
    -2.*CoeffJijsm(6)*Fields(3,1,x1,y0,z1) &
    +2.*CoeffJijsm(6)*Fields(3,1,x1,y0,zi1) &
    -2.*CoeffJijsm(7)*Fields(3,1,x1,y1,z1) &
    +2.*CoeffJijsm(7)*Fields(3,1,x1,y1,zi1) &
    -2.*CoeffJijsm(7)*Fields(3,1,x1,yi1,z1) &
    +2.*CoeffJijsm(7)*Fields(3,1,x1,yi1,zi1) &
    +2.*CoeffJijsm(6)*Fields(3,1,xi1,y0,z1) &
    -2.*CoeffJijsm(6)*Fields(3,1,xi1,y0,zi1) &
    +2.*CoeffJijsm(7)*Fields(3,1,xi1,y1,z1) &
    -2.*CoeffJijsm(7)*Fields(3,1,xi1,y1,zi1) &
    +2.*CoeffJijsm(7)*Fields(3,1,xi1,yi1,z1) &
    -2.*CoeffJijsm(7)*Fields(3,1,xi1,yi1,zi1)
  
  ForcesJijsm(2,1,x0,y0,z0) = &
    ForcesJijsm(2,1,x0,y0,z0) &
    -2.*CoeffJijsm(6)*Fields(1,1,x1,y1,z0) &
    -2.*CoeffJijsm(7)*Fields(1,1,x1,y1,z1) &
    -2.*CoeffJijsm(7)*Fields(1,1,x1,y1,zi1) &
    +2.*CoeffJijsm(6)*Fields(1,1,x1,yi1,z0) &
    +2.*CoeffJijsm(7)*Fields(1,1,x1,yi1,z1) &
    +2.*CoeffJijsm(7)*Fields(1,1,x1,yi1,zi1) &
    +2.*CoeffJijsm(6)*Fields(1,1,xi1,y1,z0) &
    +2.*CoeffJijsm(7)*Fields(1,1,xi1,y1,z1) &
    +2.*CoeffJijsm(7)*Fields(1,1,xi1,y1,zi1) &
    -2.*CoeffJijsm(6)*Fields(1,1,xi1,yi1,z0) &
    -2.*CoeffJijsm(7)*Fields(1,1,xi1,yi1,z1) &
    -2.*CoeffJijsm(7)*Fields(1,1,xi1,yi1,zi1) &
    -2.*CoeffJijsm(2)*Fields(2,1,x0,y0,z1) &
    -2.*CoeffJijsm(2)*Fields(2,1,x0,y0,zi1) &
    -2.*CoeffJijsm(1)*Fields(2,1,x0,y1,z0) &
    -2.*CoeffJijsm(3)*Fields(2,1,x0,y1,z1) &
    -2.*CoeffJijsm(3)*Fields(2,1,x0,y1,zi1) &
    -2.*CoeffJijsm(1)*Fields(2,1,x0,yi1,z0) &
    -2.*CoeffJijsm(3)*Fields(2,1,x0,yi1,z1) &
    -2.*CoeffJijsm(3)*Fields(2,1,x0,yi1,zi1) &
    -2.*CoeffJijsm(2)*Fields(2,1,x1,y0,z0) &
    -2.*CoeffJijsm(4)*Fields(2,1,x1,y0,z1) &
    -2.*CoeffJijsm(4)*Fields(2,1,x1,y0,zi1) &
    -2.*CoeffJijsm(3)*Fields(2,1,x1,y1,z0) &
    -2.*CoeffJijsm(5)*Fields(2,1,x1,y1,z1) &
    -2.*CoeffJijsm(5)*Fields(2,1,x1,y1,zi1) &
    -2.*CoeffJijsm(3)*Fields(2,1,x1,yi1,z0) &
    -2.*CoeffJijsm(5)*Fields(2,1,x1,yi1,z1) &
    -2.*CoeffJijsm(5)*Fields(2,1,x1,yi1,zi1) &
    -2.*CoeffJijsm(2)*Fields(2,1,xi1,y0,z0) &
    -2.*CoeffJijsm(4)*Fields(2,1,xi1,y0,z1) &
    -2.*CoeffJijsm(4)*Fields(2,1,xi1,y0,zi1) &
    -2.*CoeffJijsm(3)*Fields(2,1,xi1,y1,z0) &
    -2.*CoeffJijsm(5)*Fields(2,1,xi1,y1,z1) &
    -2.*CoeffJijsm(5)*Fields(2,1,xi1,y1,zi1) &
    -2.*CoeffJijsm(3)*Fields(2,1,xi1,yi1,z0) &
    -2.*CoeffJijsm(5)*Fields(2,1,xi1,yi1,z1) &
    -2.*CoeffJijsm(5)*Fields(2,1,xi1,yi1,zi1) &
    -2.*CoeffJijsm(6)*Fields(3,1,x0,y1,z1) &
    +2.*CoeffJijsm(6)*Fields(3,1,x0,y1,zi1) &
    +2.*CoeffJijsm(6)*Fields(3,1,x0,yi1,z1) &
    -2.*CoeffJijsm(6)*Fields(3,1,x0,yi1,zi1) &
    -2.*CoeffJijsm(7)*Fields(3,1,x1,y1,z1) &
    +2.*CoeffJijsm(7)*Fields(3,1,x1,y1,zi1) &
    +2.*CoeffJijsm(7)*Fields(3,1,x1,yi1,z1) &
    -2.*CoeffJijsm(7)*Fields(3,1,x1,yi1,zi1) &
    -2.*CoeffJijsm(7)*Fields(3,1,xi1,y1,z1) &
    +2.*CoeffJijsm(7)*Fields(3,1,xi1,y1,zi1) &
    +2.*CoeffJijsm(7)*Fields(3,1,xi1,yi1,z1) &
    -2.*CoeffJijsm(7)*Fields(3,1,xi1,yi1,zi1)
  
  ForcesJijsm(3,1,x0,y0,z0) = &
    ForcesJijsm(3,1,x0,y0,z0) &
    -2.*CoeffJijsm(6)*Fields(1,1,x1,y0,z1) &
    +2.*CoeffJijsm(6)*Fields(1,1,x1,y0,zi1) &
    -2.*CoeffJijsm(7)*Fields(1,1,x1,y1,z1) &
    +2.*CoeffJijsm(7)*Fields(1,1,x1,y1,zi1) &
    -2.*CoeffJijsm(7)*Fields(1,1,x1,yi1,z1) &
    +2.*CoeffJijsm(7)*Fields(1,1,x1,yi1,zi1) &
    +2.*CoeffJijsm(6)*Fields(1,1,xi1,y0,z1) &
    -2.*CoeffJijsm(6)*Fields(1,1,xi1,y0,zi1) &
    +2.*CoeffJijsm(7)*Fields(1,1,xi1,y1,z1) &
    -2.*CoeffJijsm(7)*Fields(1,1,xi1,y1,zi1) &
    +2.*CoeffJijsm(7)*Fields(1,1,xi1,yi1,z1) &
    -2.*CoeffJijsm(7)*Fields(1,1,xi1,yi1,zi1) &
    -2.*CoeffJijsm(6)*Fields(2,1,x0,y1,z1) &
    +2.*CoeffJijsm(6)*Fields(2,1,x0,y1,zi1) &
    +2.*CoeffJijsm(6)*Fields(2,1,x0,yi1,z1) &
    -2.*CoeffJijsm(6)*Fields(2,1,x0,yi1,zi1) &
    -2.*CoeffJijsm(7)*Fields(2,1,x1,y1,z1) &
    +2.*CoeffJijsm(7)*Fields(2,1,x1,y1,zi1) &
    +2.*CoeffJijsm(7)*Fields(2,1,x1,yi1,z1) &
    -2.*CoeffJijsm(7)*Fields(2,1,x1,yi1,zi1) &
    -2.*CoeffJijsm(7)*Fields(2,1,xi1,y1,z1) &
    +2.*CoeffJijsm(7)*Fields(2,1,xi1,y1,zi1) &
    +2.*CoeffJijsm(7)*Fields(2,1,xi1,yi1,z1) &
    -2.*CoeffJijsm(7)*Fields(2,1,xi1,yi1,zi1) &
    -2.*CoeffJijsm(1)*Fields(3,1,x0,y0,z1) &
    -2.*CoeffJijsm(1)*Fields(3,1,x0,y0,zi1) &
    -2.*CoeffJijsm(2)*Fields(3,1,x0,y1,z0) &
    -2.*CoeffJijsm(3)*Fields(3,1,x0,y1,z1) &
    -2.*CoeffJijsm(3)*Fields(3,1,x0,y1,zi1) &
    -2.*CoeffJijsm(2)*Fields(3,1,x0,yi1,z0) &
    -2.*CoeffJijsm(3)*Fields(3,1,x0,yi1,z1) &
    -2.*CoeffJijsm(3)*Fields(3,1,x0,yi1,zi1) &
    -2.*CoeffJijsm(2)*Fields(3,1,x1,y0,z0) &
    -2.*CoeffJijsm(3)*Fields(3,1,x1,y0,z1) &
    -2.*CoeffJijsm(3)*Fields(3,1,x1,y0,zi1) &
    -2.*CoeffJijsm(4)*Fields(3,1,x1,y1,z0) &
    -2.*CoeffJijsm(5)*Fields(3,1,x1,y1,z1) &
    -2.*CoeffJijsm(5)*Fields(3,1,x1,y1,zi1) &
    -2.*CoeffJijsm(4)*Fields(3,1,x1,yi1,z0) &
    -2.*CoeffJijsm(5)*Fields(3,1,x1,yi1,z1) &
    -2.*CoeffJijsm(5)*Fields(3,1,x1,yi1,zi1) &
    -2.*CoeffJijsm(2)*Fields(3,1,xi1,y0,z0) &
    -2.*CoeffJijsm(3)*Fields(3,1,xi1,y0,z1) &
    -2.*CoeffJijsm(3)*Fields(3,1,xi1,y0,zi1) &
    -2.*CoeffJijsm(4)*Fields(3,1,xi1,y1,z0) &
    -2.*CoeffJijsm(5)*Fields(3,1,xi1,y1,z1) &
    -2.*CoeffJijsm(5)*Fields(3,1,xi1,y1,zi1) &
    -2.*CoeffJijsm(4)*Fields(3,1,xi1,yi1,z0) &
    -2.*CoeffJijsm(5)*Fields(3,1,xi1,yi1,z1) &
    -2.*CoeffJijsm(5)*Fields(3,1,xi1,yi1,zi1)
  
  ForcesJijsm(1,2,x0,y0,z0) = &
    ForcesJijsm(1,2,x0,y0,z0) &
    +0.
  
  ForcesJijsm(2,2,x0,y0,z0) = &
    ForcesJijsm(2,2,x0,y0,z0) &
    +0.
  
  ForcesJijsm(3,2,x0,y0,z0) = &
    ForcesJijsm(3,2,x0,y0,z0) &
    +0.
  
  ForcesJijsm(1,3,x0,y0,z0) = &
    ForcesJijsm(1,3,x0,y0,z0) &
    +0.
  
  ForcesJijsm(2,3,x0,y0,z0) = &
    ForcesJijsm(2,3,x0,y0,z0) &
    +0.
  
  ForcesJijsm(3,3,x0,y0,z0) = &
    ForcesJijsm(3,3,x0,y0,z0) &
    +0.
  
  ForcesJijsm(1,4,x0,y0,z0) = &
    ForcesJijsm(1,4,x0,y0,z0) &
    +0.
  
  ForcesJijsm(2,4,x0,y0,z0) = &
    ForcesJijsm(2,4,x0,y0,z0) &
    +0.
  
  ForcesJijsm(3,4,x0,y0,z0) = &
    ForcesJijsm(3,4,x0,y0,z0) &
    +0.
  
  ForcesJijsm(4,4,x0,y0,z0) = &
    ForcesJijsm(4,4,x0,y0,z0) &
    +0.
  
  ForcesJijsm(5,4,x0,y0,z0) = &
    ForcesJijsm(5,4,x0,y0,z0) &
    +0.
  
  ForcesJijsm(6,4,x0,y0,z0) = &
    ForcesJijsm(6,4,x0,y0,z0) &
    +0.
  
  ForcesJijhf(1,1,x0,y0,z0) = &
    ForcesJijhf(1,1,x0,y0,z0) &
    +0.
  
  ForcesJijhf(2,1,x0,y0,z0) = &
    ForcesJijhf(2,1,x0,y0,z0) &
    +0.
  
  ForcesJijhf(3,1,x0,y0,z0) = &
    ForcesJijhf(3,1,x0,y0,z0) &
    +0.
  
  ForcesJijhf(1,2,x0,y0,z0) = &
    ForcesJijhf(1,2,x0,y0,z0) &
    -2.*CoeffJijhf(2)*Fields(1,2,x0,y0,z1) &
    -2.*CoeffJijhf(2)*Fields(1,2,x0,y0,zi1) &
    -2.*CoeffJijhf(2)*Fields(1,2,x0,y1,z0) &
    -2.*CoeffJijhf(4)*Fields(1,2,x0,y1,z1) &
    -2.*CoeffJijhf(4)*Fields(1,2,x0,y1,zi1) &
    -2.*CoeffJijhf(2)*Fields(1,2,x0,yi1,z0) &
    -2.*CoeffJijhf(4)*Fields(1,2,x0,yi1,z1) &
    -2.*CoeffJijhf(4)*Fields(1,2,x0,yi1,zi1) &
    -2.*CoeffJijhf(1)*Fields(1,2,x1,y0,z0) &
    -2.*CoeffJijhf(3)*Fields(1,2,x1,y0,z1) &
    -2.*CoeffJijhf(3)*Fields(1,2,x1,y0,zi1) &
    -2.*CoeffJijhf(3)*Fields(1,2,x1,y1,z0) &
    -2.*CoeffJijhf(5)*Fields(1,2,x1,y1,z1) &
    -2.*CoeffJijhf(5)*Fields(1,2,x1,y1,zi1) &
    -2.*CoeffJijhf(3)*Fields(1,2,x1,yi1,z0) &
    -2.*CoeffJijhf(5)*Fields(1,2,x1,yi1,z1) &
    -2.*CoeffJijhf(5)*Fields(1,2,x1,yi1,zi1) &
    -2.*CoeffJijhf(1)*Fields(1,2,xi1,y0,z0) &
    -2.*CoeffJijhf(3)*Fields(1,2,xi1,y0,z1) &
    -2.*CoeffJijhf(3)*Fields(1,2,xi1,y0,zi1) &
    -2.*CoeffJijhf(3)*Fields(1,2,xi1,y1,z0) &
    -2.*CoeffJijhf(5)*Fields(1,2,xi1,y1,z1) &
    -2.*CoeffJijhf(5)*Fields(1,2,xi1,y1,zi1) &
    -2.*CoeffJijhf(3)*Fields(1,2,xi1,yi1,z0) &
    -2.*CoeffJijhf(5)*Fields(1,2,xi1,yi1,z1) &
    -2.*CoeffJijhf(5)*Fields(1,2,xi1,yi1,zi1) &
    -2.*CoeffJijhf(6)*Fields(2,2,x1,y1,z0) &
    -2.*CoeffJijhf(7)*Fields(2,2,x1,y1,z1) &
    -2.*CoeffJijhf(7)*Fields(2,2,x1,y1,zi1) &
    +2.*CoeffJijhf(6)*Fields(2,2,x1,yi1,z0) &
    +2.*CoeffJijhf(7)*Fields(2,2,x1,yi1,z1) &
    +2.*CoeffJijhf(7)*Fields(2,2,x1,yi1,zi1) &
    +2.*CoeffJijhf(6)*Fields(2,2,xi1,y1,z0) &
    +2.*CoeffJijhf(7)*Fields(2,2,xi1,y1,z1) &
    +2.*CoeffJijhf(7)*Fields(2,2,xi1,y1,zi1) &
    -2.*CoeffJijhf(6)*Fields(2,2,xi1,yi1,z0) &
    -2.*CoeffJijhf(7)*Fields(2,2,xi1,yi1,z1) &
    -2.*CoeffJijhf(7)*Fields(2,2,xi1,yi1,zi1) &
    -2.*CoeffJijhf(6)*Fields(3,2,x1,y0,z1) &
    +2.*CoeffJijhf(6)*Fields(3,2,x1,y0,zi1) &
    -2.*CoeffJijhf(7)*Fields(3,2,x1,y1,z1) &
    +2.*CoeffJijhf(7)*Fields(3,2,x1,y1,zi1) &
    -2.*CoeffJijhf(7)*Fields(3,2,x1,yi1,z1) &
    +2.*CoeffJijhf(7)*Fields(3,2,x1,yi1,zi1) &
    +2.*CoeffJijhf(6)*Fields(3,2,xi1,y0,z1) &
    -2.*CoeffJijhf(6)*Fields(3,2,xi1,y0,zi1) &
    +2.*CoeffJijhf(7)*Fields(3,2,xi1,y1,z1) &
    -2.*CoeffJijhf(7)*Fields(3,2,xi1,y1,zi1) &
    +2.*CoeffJijhf(7)*Fields(3,2,xi1,yi1,z1) &
    -2.*CoeffJijhf(7)*Fields(3,2,xi1,yi1,zi1)
  
  ForcesJijhf(2,2,x0,y0,z0) = &
    ForcesJijhf(2,2,x0,y0,z0) &
    -2.*CoeffJijhf(6)*Fields(1,2,x1,y1,z0) &
    -2.*CoeffJijhf(7)*Fields(1,2,x1,y1,z1) &
    -2.*CoeffJijhf(7)*Fields(1,2,x1,y1,zi1) &
    +2.*CoeffJijhf(6)*Fields(1,2,x1,yi1,z0) &
    +2.*CoeffJijhf(7)*Fields(1,2,x1,yi1,z1) &
    +2.*CoeffJijhf(7)*Fields(1,2,x1,yi1,zi1) &
    +2.*CoeffJijhf(6)*Fields(1,2,xi1,y1,z0) &
    +2.*CoeffJijhf(7)*Fields(1,2,xi1,y1,z1) &
    +2.*CoeffJijhf(7)*Fields(1,2,xi1,y1,zi1) &
    -2.*CoeffJijhf(6)*Fields(1,2,xi1,yi1,z0) &
    -2.*CoeffJijhf(7)*Fields(1,2,xi1,yi1,z1) &
    -2.*CoeffJijhf(7)*Fields(1,2,xi1,yi1,zi1) &
    -2.*CoeffJijhf(2)*Fields(2,2,x0,y0,z1) &
    -2.*CoeffJijhf(2)*Fields(2,2,x0,y0,zi1) &
    -2.*CoeffJijhf(1)*Fields(2,2,x0,y1,z0) &
    -2.*CoeffJijhf(3)*Fields(2,2,x0,y1,z1) &
    -2.*CoeffJijhf(3)*Fields(2,2,x0,y1,zi1) &
    -2.*CoeffJijhf(1)*Fields(2,2,x0,yi1,z0) &
    -2.*CoeffJijhf(3)*Fields(2,2,x0,yi1,z1) &
    -2.*CoeffJijhf(3)*Fields(2,2,x0,yi1,zi1) &
    -2.*CoeffJijhf(2)*Fields(2,2,x1,y0,z0) &
    -2.*CoeffJijhf(4)*Fields(2,2,x1,y0,z1) &
    -2.*CoeffJijhf(4)*Fields(2,2,x1,y0,zi1) &
    -2.*CoeffJijhf(3)*Fields(2,2,x1,y1,z0) &
    -2.*CoeffJijhf(5)*Fields(2,2,x1,y1,z1) &
    -2.*CoeffJijhf(5)*Fields(2,2,x1,y1,zi1) &
    -2.*CoeffJijhf(3)*Fields(2,2,x1,yi1,z0) &
    -2.*CoeffJijhf(5)*Fields(2,2,x1,yi1,z1) &
    -2.*CoeffJijhf(5)*Fields(2,2,x1,yi1,zi1) &
    -2.*CoeffJijhf(2)*Fields(2,2,xi1,y0,z0) &
    -2.*CoeffJijhf(4)*Fields(2,2,xi1,y0,z1) &
    -2.*CoeffJijhf(4)*Fields(2,2,xi1,y0,zi1) &
    -2.*CoeffJijhf(3)*Fields(2,2,xi1,y1,z0) &
    -2.*CoeffJijhf(5)*Fields(2,2,xi1,y1,z1) &
    -2.*CoeffJijhf(5)*Fields(2,2,xi1,y1,zi1) &
    -2.*CoeffJijhf(3)*Fields(2,2,xi1,yi1,z0) &
    -2.*CoeffJijhf(5)*Fields(2,2,xi1,yi1,z1) &
    -2.*CoeffJijhf(5)*Fields(2,2,xi1,yi1,zi1) &
    -2.*CoeffJijhf(6)*Fields(3,2,x0,y1,z1) &
    +2.*CoeffJijhf(6)*Fields(3,2,x0,y1,zi1) &
    +2.*CoeffJijhf(6)*Fields(3,2,x0,yi1,z1) &
    -2.*CoeffJijhf(6)*Fields(3,2,x0,yi1,zi1) &
    -2.*CoeffJijhf(7)*Fields(3,2,x1,y1,z1) &
    +2.*CoeffJijhf(7)*Fields(3,2,x1,y1,zi1) &
    +2.*CoeffJijhf(7)*Fields(3,2,x1,yi1,z1) &
    -2.*CoeffJijhf(7)*Fields(3,2,x1,yi1,zi1) &
    -2.*CoeffJijhf(7)*Fields(3,2,xi1,y1,z1) &
    +2.*CoeffJijhf(7)*Fields(3,2,xi1,y1,zi1) &
    +2.*CoeffJijhf(7)*Fields(3,2,xi1,yi1,z1) &
    -2.*CoeffJijhf(7)*Fields(3,2,xi1,yi1,zi1)
  
  ForcesJijhf(3,2,x0,y0,z0) = &
    ForcesJijhf(3,2,x0,y0,z0) &
    -2.*CoeffJijhf(6)*Fields(1,2,x1,y0,z1) &
    +2.*CoeffJijhf(6)*Fields(1,2,x1,y0,zi1) &
    -2.*CoeffJijhf(7)*Fields(1,2,x1,y1,z1) &
    +2.*CoeffJijhf(7)*Fields(1,2,x1,y1,zi1) &
    -2.*CoeffJijhf(7)*Fields(1,2,x1,yi1,z1) &
    +2.*CoeffJijhf(7)*Fields(1,2,x1,yi1,zi1) &
    +2.*CoeffJijhf(6)*Fields(1,2,xi1,y0,z1) &
    -2.*CoeffJijhf(6)*Fields(1,2,xi1,y0,zi1) &
    +2.*CoeffJijhf(7)*Fields(1,2,xi1,y1,z1) &
    -2.*CoeffJijhf(7)*Fields(1,2,xi1,y1,zi1) &
    +2.*CoeffJijhf(7)*Fields(1,2,xi1,yi1,z1) &
    -2.*CoeffJijhf(7)*Fields(1,2,xi1,yi1,zi1) &
    -2.*CoeffJijhf(6)*Fields(2,2,x0,y1,z1) &
    +2.*CoeffJijhf(6)*Fields(2,2,x0,y1,zi1) &
    +2.*CoeffJijhf(6)*Fields(2,2,x0,yi1,z1) &
    -2.*CoeffJijhf(6)*Fields(2,2,x0,yi1,zi1) &
    -2.*CoeffJijhf(7)*Fields(2,2,x1,y1,z1) &
    +2.*CoeffJijhf(7)*Fields(2,2,x1,y1,zi1) &
    +2.*CoeffJijhf(7)*Fields(2,2,x1,yi1,z1) &
    -2.*CoeffJijhf(7)*Fields(2,2,x1,yi1,zi1) &
    -2.*CoeffJijhf(7)*Fields(2,2,xi1,y1,z1) &
    +2.*CoeffJijhf(7)*Fields(2,2,xi1,y1,zi1) &
    +2.*CoeffJijhf(7)*Fields(2,2,xi1,yi1,z1) &
    -2.*CoeffJijhf(7)*Fields(2,2,xi1,yi1,zi1) &
    -2.*CoeffJijhf(1)*Fields(3,2,x0,y0,z1) &
    -2.*CoeffJijhf(1)*Fields(3,2,x0,y0,zi1) &
    -2.*CoeffJijhf(2)*Fields(3,2,x0,y1,z0) &
    -2.*CoeffJijhf(3)*Fields(3,2,x0,y1,z1) &
    -2.*CoeffJijhf(3)*Fields(3,2,x0,y1,zi1) &
    -2.*CoeffJijhf(2)*Fields(3,2,x0,yi1,z0) &
    -2.*CoeffJijhf(3)*Fields(3,2,x0,yi1,z1) &
    -2.*CoeffJijhf(3)*Fields(3,2,x0,yi1,zi1) &
    -2.*CoeffJijhf(2)*Fields(3,2,x1,y0,z0) &
    -2.*CoeffJijhf(3)*Fields(3,2,x1,y0,z1) &
    -2.*CoeffJijhf(3)*Fields(3,2,x1,y0,zi1) &
    -2.*CoeffJijhf(4)*Fields(3,2,x1,y1,z0) &
    -2.*CoeffJijhf(5)*Fields(3,2,x1,y1,z1) &
    -2.*CoeffJijhf(5)*Fields(3,2,x1,y1,zi1) &
    -2.*CoeffJijhf(4)*Fields(3,2,x1,yi1,z0) &
    -2.*CoeffJijhf(5)*Fields(3,2,x1,yi1,z1) &
    -2.*CoeffJijhf(5)*Fields(3,2,x1,yi1,zi1) &
    -2.*CoeffJijhf(2)*Fields(3,2,xi1,y0,z0) &
    -2.*CoeffJijhf(3)*Fields(3,2,xi1,y0,z1) &
    -2.*CoeffJijhf(3)*Fields(3,2,xi1,y0,zi1) &
    -2.*CoeffJijhf(4)*Fields(3,2,xi1,y1,z0) &
    -2.*CoeffJijhf(5)*Fields(3,2,xi1,y1,z1) &
    -2.*CoeffJijhf(5)*Fields(3,2,xi1,y1,zi1) &
    -2.*CoeffJijhf(4)*Fields(3,2,xi1,yi1,z0) &
    -2.*CoeffJijhf(5)*Fields(3,2,xi1,yi1,z1) &
    -2.*CoeffJijhf(5)*Fields(3,2,xi1,yi1,zi1)
  
  ForcesJijhf(1,3,x0,y0,z0) = &
    ForcesJijhf(1,3,x0,y0,z0) &
    +0.
  
  ForcesJijhf(2,3,x0,y0,z0) = &
    ForcesJijhf(2,3,x0,y0,z0) &
    +0.
  
  ForcesJijhf(3,3,x0,y0,z0) = &
    ForcesJijhf(3,3,x0,y0,z0) &
    +0.
  
  ForcesJijhf(1,4,x0,y0,z0) = &
    ForcesJijhf(1,4,x0,y0,z0) &
    +0.
  
  ForcesJijhf(2,4,x0,y0,z0) = &
    ForcesJijhf(2,4,x0,y0,z0) &
    +0.
  
  ForcesJijhf(3,4,x0,y0,z0) = &
    ForcesJijhf(3,4,x0,y0,z0) &
    +0.
  
  ForcesJijhf(4,4,x0,y0,z0) = &
    ForcesJijhf(4,4,x0,y0,z0) &
    +0.
  
  ForcesJijhf(5,4,x0,y0,z0) = &
    ForcesJijhf(5,4,x0,y0,z0) &
    +0.
  
  ForcesJijhf(6,4,x0,y0,z0) = &
    ForcesJijhf(6,4,x0,y0,z0) &
    +0.
  
  Forcesu(1,1,x0,y0,z0) = &
    Forcesu(1,1,x0,y0,z0) &
    +0.
  
  Forcesu(2,1,x0,y0,z0) = &
    Forcesu(2,1,x0,y0,z0) &
    +0.
  
  Forcesu(3,1,x0,y0,z0) = &
    Forcesu(3,1,x0,y0,z0) &
    +0.
  
  Forcesu(1,2,x0,y0,z0) = &
    Forcesu(1,2,x0,y0,z0) &
    +0.
  
  Forcesu(2,2,x0,y0,z0) = &
    Forcesu(2,2,x0,y0,z0) &
    +0.
  
  Forcesu(3,2,x0,y0,z0) = &
    Forcesu(3,2,x0,y0,z0) &
    +0.
  
  Forcesu(1,3,x0,y0,z0) = &
    Forcesu(1,3,x0,y0,z0) &
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
  
  Forcesu(2,3,x0,y0,z0) = &
    Forcesu(2,3,x0,y0,z0) &
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
  
  Forcesu(3,3,x0,y0,z0) = &
    Forcesu(3,3,x0,y0,z0) &
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
  
  Forcesu(1,4,x0,y0,z0) = &
    Forcesu(1,4,x0,y0,z0) &
    +0.
  
  Forcesu(2,4,x0,y0,z0) = &
    Forcesu(2,4,x0,y0,z0) &
    +0.
  
  Forcesu(3,4,x0,y0,z0) = &
    Forcesu(3,4,x0,y0,z0) &
    +0.
  
  Forcesu(4,4,x0,y0,z0) = &
    Forcesu(4,4,x0,y0,z0) &
    +0.
  
  Forcesu(5,4,x0,y0,z0) = &
    Forcesu(5,4,x0,y0,z0) &
    +0.
  
  Forcesu(6,4,x0,y0,z0) = &
    Forcesu(6,4,x0,y0,z0) &
    +0.
  
  ForcesDispu(1,1,x0,y0,z0) = &
    ForcesDispu(1,1,x0,y0,z0) &
    +0.
  
  ForcesDispu(2,1,x0,y0,z0) = &
    ForcesDispu(2,1,x0,y0,z0) &
    +0.
  
  ForcesDispu(3,1,x0,y0,z0) = &
    ForcesDispu(3,1,x0,y0,z0) &
    +0.
  
  ForcesDispu(1,2,x0,y0,z0) = &
    ForcesDispu(1,2,x0,y0,z0) &
    +0.
  
  ForcesDispu(2,2,x0,y0,z0) = &
    ForcesDispu(2,2,x0,y0,z0) &
    +0.
  
  ForcesDispu(3,2,x0,y0,z0) = &
    ForcesDispu(3,2,x0,y0,z0) &
    +0.
  
  ForcesDispu(1,3,x0,y0,z0) = &
    ForcesDispu(1,3,x0,y0,z0) &
    +0.25*CoeffDispu(5)*Fields(1,1,x0,y0,z0)**2 &
    +0.25*CoeffDispu(5)*Fields(1,1,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(5)*Fields(1,1,x0,yi1,z0)**2 &
    +0.25*CoeffDispu(5)*Fields(1,1,x0,yi1,zi1)**2 &
    -0.25*CoeffDispu(5)*Fields(1,1,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(5)*Fields(1,1,xi1,y0,zi1)**2 &
    -0.25*CoeffDispu(5)*Fields(1,1,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(5)*Fields(1,1,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(8)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
    +0.25*CoeffDispu(9)*Fields(1,2,x0,y0,z0)**2 &
    +0.25*CoeffDispu(8)*Fields(1,1,x0,y0,zi1)*Fields(1,2,x0,y0,zi1) &
    +0.25*CoeffDispu(9)*Fields(1,2,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(8)*Fields(1,1,x0,yi1,z0)*Fields(1,2,x0,yi1,z0) &
    +0.25*CoeffDispu(9)*Fields(1,2,x0,yi1,z0)**2 &
    +0.25*CoeffDispu(8)*Fields(1,1,x0,yi1,zi1)*Fields(1,2,x0,yi1,zi1) &
    +0.25*CoeffDispu(9)*Fields(1,2,x0,yi1,zi1)**2 &
    -0.25*CoeffDispu(8)*Fields(1,1,xi1,y0,z0)*Fields(1,2,xi1,y0,z0) &
    -0.25*CoeffDispu(9)*Fields(1,2,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(8)*Fields(1,1,xi1,y0,zi1)*Fields(1,2,xi1,y0,zi1) &
    -0.25*CoeffDispu(9)*Fields(1,2,xi1,y0,zi1)**2 &
    -0.25*CoeffDispu(8)*Fields(1,1,xi1,yi1,z0)*Fields(1,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(9)*Fields(1,2,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(8)*Fields(1,1,xi1,yi1,zi1)*Fields(1,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(9)*Fields(1,2,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,x0,y0,z0)*Fields(2,1,x0,y0,z0) &
    +0.25*CoeffDispu(2)*Fields(1,2,x0,y0,z0)*Fields(2,1,x0,y0,z0) &
    +0.25*CoeffDispu(4)*Fields(2,1,x0,y0,z0)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,x0,y0,zi1)*Fields(2,1,x0,y0,zi1) &
    +0.25*CoeffDispu(2)*Fields(1,2,x0,y0,zi1)*Fields(2,1,x0,y0,zi1) &
    +0.25*CoeffDispu(4)*Fields(2,1,x0,y0,zi1)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,x0,yi1,z0)*Fields(2,1,x0,yi1,z0) &
    -0.25*CoeffDispu(2)*Fields(1,2,x0,yi1,z0)*Fields(2,1,x0,yi1,z0) &
    +0.25*CoeffDispu(4)*Fields(2,1,x0,yi1,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,x0,yi1,zi1)*Fields(2,1,x0,yi1,zi1) &
    -0.25*CoeffDispu(2)*Fields(1,2,x0,yi1,zi1)*Fields(2,1,x0,yi1,zi1) &
    +0.25*CoeffDispu(4)*Fields(2,1,x0,yi1,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,xi1,y0,z0)*Fields(2,1,xi1,y0,z0) &
    +0.25*CoeffDispu(2)*Fields(1,2,xi1,y0,z0)*Fields(2,1,xi1,y0,z0) &
    -0.25*CoeffDispu(4)*Fields(2,1,xi1,y0,z0)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,xi1,y0,zi1)*Fields(2,1,xi1,y0,zi1) &
    +0.25*CoeffDispu(2)*Fields(1,2,xi1,y0,zi1)*Fields(2,1,xi1,y0,zi1) &
    -0.25*CoeffDispu(4)*Fields(2,1,xi1,y0,zi1)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,xi1,yi1,z0)*Fields(2,1,xi1,yi1,z0) &
    -0.25*CoeffDispu(2)*Fields(1,2,xi1,yi1,z0)*Fields(2,1,xi1,yi1,z0) &
    -0.25*CoeffDispu(4)*Fields(2,1,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,xi1,yi1,zi1)*Fields(2,1,xi1,yi1,zi1) &
    -0.25*CoeffDispu(2)*Fields(1,2,xi1,yi1,zi1)*Fields(2,1,xi1,yi1,zi1) &
    -0.25*CoeffDispu(4)*Fields(2,1,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +0.25*CoeffDispu(3)*Fields(1,2,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +0.25*CoeffDispu(6)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +0.25*CoeffDispu(7)*Fields(2,2,x0,y0,z0)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,x0,y0,zi1)*Fields(2,2,x0,y0,zi1) &
    +0.25*CoeffDispu(3)*Fields(1,2,x0,y0,zi1)*Fields(2,2,x0,y0,zi1) &
    +0.25*CoeffDispu(6)*Fields(2,1,x0,y0,zi1)*Fields(2,2,x0,y0,zi1) &
    +0.25*CoeffDispu(7)*Fields(2,2,x0,y0,zi1)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,x0,yi1,z0)*Fields(2,2,x0,yi1,z0) &
    -0.25*CoeffDispu(3)*Fields(1,2,x0,yi1,z0)*Fields(2,2,x0,yi1,z0) &
    +0.25*CoeffDispu(6)*Fields(2,1,x0,yi1,z0)*Fields(2,2,x0,yi1,z0) &
    +0.25*CoeffDispu(7)*Fields(2,2,x0,yi1,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,x0,yi1,zi1)*Fields(2,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(3)*Fields(1,2,x0,yi1,zi1)*Fields(2,2,x0,yi1,zi1) &
    +0.25*CoeffDispu(6)*Fields(2,1,x0,yi1,zi1)*Fields(2,2,x0,yi1,zi1) &
    +0.25*CoeffDispu(7)*Fields(2,2,x0,yi1,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,xi1,y0,z0)*Fields(2,2,xi1,y0,z0) &
    +0.25*CoeffDispu(3)*Fields(1,2,xi1,y0,z0)*Fields(2,2,xi1,y0,z0) &
    -0.25*CoeffDispu(6)*Fields(2,1,xi1,y0,z0)*Fields(2,2,xi1,y0,z0) &
    -0.25*CoeffDispu(7)*Fields(2,2,xi1,y0,z0)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,xi1,y0,zi1)*Fields(2,2,xi1,y0,zi1) &
    +0.25*CoeffDispu(3)*Fields(1,2,xi1,y0,zi1)*Fields(2,2,xi1,y0,zi1) &
    -0.25*CoeffDispu(6)*Fields(2,1,xi1,y0,zi1)*Fields(2,2,xi1,y0,zi1) &
    -0.25*CoeffDispu(7)*Fields(2,2,xi1,y0,zi1)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,xi1,yi1,z0)*Fields(2,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(3)*Fields(1,2,xi1,yi1,z0)*Fields(2,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(6)*Fields(2,1,xi1,yi1,z0)*Fields(2,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(7)*Fields(2,2,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,xi1,yi1,zi1)*Fields(2,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(3)*Fields(1,2,xi1,yi1,zi1)*Fields(2,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(6)*Fields(2,1,xi1,yi1,zi1)*Fields(2,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(7)*Fields(2,2,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
    +0.25*CoeffDispu(2)*Fields(1,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
    +0.25*CoeffDispu(4)*Fields(3,1,x0,y0,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,x0,y0,zi1)*Fields(3,1,x0,y0,zi1) &
    -0.25*CoeffDispu(2)*Fields(1,2,x0,y0,zi1)*Fields(3,1,x0,y0,zi1) &
    +0.25*CoeffDispu(4)*Fields(3,1,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,x0,yi1,z0)*Fields(3,1,x0,yi1,z0) &
    +0.25*CoeffDispu(2)*Fields(1,2,x0,yi1,z0)*Fields(3,1,x0,yi1,z0) &
    +0.25*CoeffDispu(4)*Fields(3,1,x0,yi1,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,x0,yi1,zi1)*Fields(3,1,x0,yi1,zi1) &
    -0.25*CoeffDispu(2)*Fields(1,2,x0,yi1,zi1)*Fields(3,1,x0,yi1,zi1) &
    +0.25*CoeffDispu(4)*Fields(3,1,x0,yi1,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,xi1,y0,z0)*Fields(3,1,xi1,y0,z0) &
    +0.25*CoeffDispu(2)*Fields(1,2,xi1,y0,z0)*Fields(3,1,xi1,y0,z0) &
    -0.25*CoeffDispu(4)*Fields(3,1,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,xi1,y0,zi1)*Fields(3,1,xi1,y0,zi1) &
    -0.25*CoeffDispu(2)*Fields(1,2,xi1,y0,zi1)*Fields(3,1,xi1,y0,zi1) &
    -0.25*CoeffDispu(4)*Fields(3,1,xi1,y0,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,xi1,yi1,z0)*Fields(3,1,xi1,yi1,z0) &
    +0.25*CoeffDispu(2)*Fields(1,2,xi1,yi1,z0)*Fields(3,1,xi1,yi1,z0) &
    -0.25*CoeffDispu(4)*Fields(3,1,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,xi1,yi1,zi1)*Fields(3,1,xi1,yi1,zi1) &
    -0.25*CoeffDispu(2)*Fields(1,2,xi1,yi1,zi1)*Fields(3,1,xi1,yi1,zi1) &
    -0.25*CoeffDispu(4)*Fields(3,1,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +0.25*CoeffDispu(3)*Fields(1,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +0.25*CoeffDispu(6)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +0.25*CoeffDispu(7)*Fields(3,2,x0,y0,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
    -0.25*CoeffDispu(3)*Fields(1,2,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
    +0.25*CoeffDispu(6)*Fields(3,1,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
    +0.25*CoeffDispu(7)*Fields(3,2,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
    +0.25*CoeffDispu(3)*Fields(1,2,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
    +0.25*CoeffDispu(6)*Fields(3,1,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
    +0.25*CoeffDispu(7)*Fields(3,2,x0,yi1,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(3)*Fields(1,2,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
    +0.25*CoeffDispu(6)*Fields(3,1,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
    +0.25*CoeffDispu(7)*Fields(3,2,x0,yi1,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
    +0.25*CoeffDispu(3)*Fields(1,2,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
    -0.25*CoeffDispu(6)*Fields(3,1,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
    -0.25*CoeffDispu(7)*Fields(3,2,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
    -0.25*CoeffDispu(3)*Fields(1,2,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
    -0.25*CoeffDispu(6)*Fields(3,1,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
    -0.25*CoeffDispu(7)*Fields(3,2,xi1,y0,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
    +0.25*CoeffDispu(3)*Fields(1,2,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(6)*Fields(3,1,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(7)*Fields(3,2,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(3)*Fields(1,2,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(6)*Fields(3,1,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(7)*Fields(3,2,xi1,yi1,zi1)**2
  
  ForcesDispu(2,3,x0,y0,z0) = &
    ForcesDispu(2,3,x0,y0,z0) &
    +0.25*CoeffDispu(4)*Fields(1,1,x0,y0,z0)**2 &
    +0.25*CoeffDispu(4)*Fields(1,1,x0,y0,zi1)**2 &
    -0.25*CoeffDispu(4)*Fields(1,1,x0,yi1,z0)**2 &
    -0.25*CoeffDispu(4)*Fields(1,1,x0,yi1,zi1)**2 &
    +0.25*CoeffDispu(4)*Fields(1,1,xi1,y0,z0)**2 &
    +0.25*CoeffDispu(4)*Fields(1,1,xi1,y0,zi1)**2 &
    -0.25*CoeffDispu(4)*Fields(1,1,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(4)*Fields(1,1,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(6)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
    +0.25*CoeffDispu(7)*Fields(1,2,x0,y0,z0)**2 &
    +0.25*CoeffDispu(6)*Fields(1,1,x0,y0,zi1)*Fields(1,2,x0,y0,zi1) &
    +0.25*CoeffDispu(7)*Fields(1,2,x0,y0,zi1)**2 &
    -0.25*CoeffDispu(6)*Fields(1,1,x0,yi1,z0)*Fields(1,2,x0,yi1,z0) &
    -0.25*CoeffDispu(7)*Fields(1,2,x0,yi1,z0)**2 &
    -0.25*CoeffDispu(6)*Fields(1,1,x0,yi1,zi1)*Fields(1,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(7)*Fields(1,2,x0,yi1,zi1)**2 &
    +0.25*CoeffDispu(6)*Fields(1,1,xi1,y0,z0)*Fields(1,2,xi1,y0,z0) &
    +0.25*CoeffDispu(7)*Fields(1,2,xi1,y0,z0)**2 &
    +0.25*CoeffDispu(6)*Fields(1,1,xi1,y0,zi1)*Fields(1,2,xi1,y0,zi1) &
    +0.25*CoeffDispu(7)*Fields(1,2,xi1,y0,zi1)**2 &
    -0.25*CoeffDispu(6)*Fields(1,1,xi1,yi1,z0)*Fields(1,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(7)*Fields(1,2,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(6)*Fields(1,1,xi1,yi1,zi1)*Fields(1,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(7)*Fields(1,2,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,x0,y0,z0)*Fields(2,1,x0,y0,z0) &
    +0.25*CoeffDispu(2)*Fields(1,2,x0,y0,z0)*Fields(2,1,x0,y0,z0) &
    +0.25*CoeffDispu(5)*Fields(2,1,x0,y0,z0)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,x0,y0,zi1)*Fields(2,1,x0,y0,zi1) &
    +0.25*CoeffDispu(2)*Fields(1,2,x0,y0,zi1)*Fields(2,1,x0,y0,zi1) &
    +0.25*CoeffDispu(5)*Fields(2,1,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,x0,yi1,z0)*Fields(2,1,x0,yi1,z0) &
    +0.25*CoeffDispu(2)*Fields(1,2,x0,yi1,z0)*Fields(2,1,x0,yi1,z0) &
    -0.25*CoeffDispu(5)*Fields(2,1,x0,yi1,z0)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,x0,yi1,zi1)*Fields(2,1,x0,yi1,zi1) &
    +0.25*CoeffDispu(2)*Fields(1,2,x0,yi1,zi1)*Fields(2,1,x0,yi1,zi1) &
    -0.25*CoeffDispu(5)*Fields(2,1,x0,yi1,zi1)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,xi1,y0,z0)*Fields(2,1,xi1,y0,z0) &
    -0.25*CoeffDispu(2)*Fields(1,2,xi1,y0,z0)*Fields(2,1,xi1,y0,z0) &
    +0.25*CoeffDispu(5)*Fields(2,1,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,xi1,y0,zi1)*Fields(2,1,xi1,y0,zi1) &
    -0.25*CoeffDispu(2)*Fields(1,2,xi1,y0,zi1)*Fields(2,1,xi1,y0,zi1) &
    +0.25*CoeffDispu(5)*Fields(2,1,xi1,y0,zi1)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,xi1,yi1,z0)*Fields(2,1,xi1,yi1,z0) &
    -0.25*CoeffDispu(2)*Fields(1,2,xi1,yi1,z0)*Fields(2,1,xi1,yi1,z0) &
    -0.25*CoeffDispu(5)*Fields(2,1,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,xi1,yi1,zi1)*Fields(2,1,xi1,yi1,zi1) &
    -0.25*CoeffDispu(2)*Fields(1,2,xi1,yi1,zi1)*Fields(2,1,xi1,yi1,zi1) &
    -0.25*CoeffDispu(5)*Fields(2,1,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +0.25*CoeffDispu(3)*Fields(1,2,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +0.25*CoeffDispu(8)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +0.25*CoeffDispu(9)*Fields(2,2,x0,y0,z0)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,x0,y0,zi1)*Fields(2,2,x0,y0,zi1) &
    +0.25*CoeffDispu(3)*Fields(1,2,x0,y0,zi1)*Fields(2,2,x0,y0,zi1) &
    +0.25*CoeffDispu(8)*Fields(2,1,x0,y0,zi1)*Fields(2,2,x0,y0,zi1) &
    +0.25*CoeffDispu(9)*Fields(2,2,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,x0,yi1,z0)*Fields(2,2,x0,yi1,z0) &
    +0.25*CoeffDispu(3)*Fields(1,2,x0,yi1,z0)*Fields(2,2,x0,yi1,z0) &
    -0.25*CoeffDispu(8)*Fields(2,1,x0,yi1,z0)*Fields(2,2,x0,yi1,z0) &
    -0.25*CoeffDispu(9)*Fields(2,2,x0,yi1,z0)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,x0,yi1,zi1)*Fields(2,2,x0,yi1,zi1) &
    +0.25*CoeffDispu(3)*Fields(1,2,x0,yi1,zi1)*Fields(2,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(8)*Fields(2,1,x0,yi1,zi1)*Fields(2,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(9)*Fields(2,2,x0,yi1,zi1)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,xi1,y0,z0)*Fields(2,2,xi1,y0,z0) &
    -0.25*CoeffDispu(3)*Fields(1,2,xi1,y0,z0)*Fields(2,2,xi1,y0,z0) &
    +0.25*CoeffDispu(8)*Fields(2,1,xi1,y0,z0)*Fields(2,2,xi1,y0,z0) &
    +0.25*CoeffDispu(9)*Fields(2,2,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,xi1,y0,zi1)*Fields(2,2,xi1,y0,zi1) &
    -0.25*CoeffDispu(3)*Fields(1,2,xi1,y0,zi1)*Fields(2,2,xi1,y0,zi1) &
    +0.25*CoeffDispu(8)*Fields(2,1,xi1,y0,zi1)*Fields(2,2,xi1,y0,zi1) &
    +0.25*CoeffDispu(9)*Fields(2,2,xi1,y0,zi1)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,xi1,yi1,z0)*Fields(2,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(3)*Fields(1,2,xi1,yi1,z0)*Fields(2,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(8)*Fields(2,1,xi1,yi1,z0)*Fields(2,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(9)*Fields(2,2,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,xi1,yi1,zi1)*Fields(2,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(3)*Fields(1,2,xi1,yi1,zi1)*Fields(2,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(8)*Fields(2,1,xi1,yi1,zi1)*Fields(2,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(9)*Fields(2,2,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
    +0.25*CoeffDispu(2)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
    +0.25*CoeffDispu(4)*Fields(3,1,x0,y0,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(2,1,x0,y0,zi1)*Fields(3,1,x0,y0,zi1) &
    -0.25*CoeffDispu(2)*Fields(2,2,x0,y0,zi1)*Fields(3,1,x0,y0,zi1) &
    +0.25*CoeffDispu(4)*Fields(3,1,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(2,1,x0,yi1,z0)*Fields(3,1,x0,yi1,z0) &
    +0.25*CoeffDispu(2)*Fields(2,2,x0,yi1,z0)*Fields(3,1,x0,yi1,z0) &
    -0.25*CoeffDispu(4)*Fields(3,1,x0,yi1,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(2,1,x0,yi1,zi1)*Fields(3,1,x0,yi1,zi1) &
    -0.25*CoeffDispu(2)*Fields(2,2,x0,yi1,zi1)*Fields(3,1,x0,yi1,zi1) &
    -0.25*CoeffDispu(4)*Fields(3,1,x0,yi1,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(2,1,xi1,y0,z0)*Fields(3,1,xi1,y0,z0) &
    +0.25*CoeffDispu(2)*Fields(2,2,xi1,y0,z0)*Fields(3,1,xi1,y0,z0) &
    +0.25*CoeffDispu(4)*Fields(3,1,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(2,1,xi1,y0,zi1)*Fields(3,1,xi1,y0,zi1) &
    -0.25*CoeffDispu(2)*Fields(2,2,xi1,y0,zi1)*Fields(3,1,xi1,y0,zi1) &
    +0.25*CoeffDispu(4)*Fields(3,1,xi1,y0,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(2,1,xi1,yi1,z0)*Fields(3,1,xi1,yi1,z0) &
    +0.25*CoeffDispu(2)*Fields(2,2,xi1,yi1,z0)*Fields(3,1,xi1,yi1,z0) &
    -0.25*CoeffDispu(4)*Fields(3,1,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(2,1,xi1,yi1,zi1)*Fields(3,1,xi1,yi1,zi1) &
    -0.25*CoeffDispu(2)*Fields(2,2,xi1,yi1,zi1)*Fields(3,1,xi1,yi1,zi1) &
    -0.25*CoeffDispu(4)*Fields(3,1,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(2,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +0.25*CoeffDispu(3)*Fields(2,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +0.25*CoeffDispu(6)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +0.25*CoeffDispu(7)*Fields(3,2,x0,y0,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(2,1,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
    -0.25*CoeffDispu(3)*Fields(2,2,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
    +0.25*CoeffDispu(6)*Fields(3,1,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
    +0.25*CoeffDispu(7)*Fields(3,2,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(2,1,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
    +0.25*CoeffDispu(3)*Fields(2,2,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
    -0.25*CoeffDispu(6)*Fields(3,1,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
    -0.25*CoeffDispu(7)*Fields(3,2,x0,yi1,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(2,1,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(3)*Fields(2,2,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(6)*Fields(3,1,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(7)*Fields(3,2,x0,yi1,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(2,1,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
    +0.25*CoeffDispu(3)*Fields(2,2,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
    +0.25*CoeffDispu(6)*Fields(3,1,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
    +0.25*CoeffDispu(7)*Fields(3,2,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(2,1,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
    -0.25*CoeffDispu(3)*Fields(2,2,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
    +0.25*CoeffDispu(6)*Fields(3,1,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
    +0.25*CoeffDispu(7)*Fields(3,2,xi1,y0,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(2,1,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
    +0.25*CoeffDispu(3)*Fields(2,2,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(6)*Fields(3,1,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(7)*Fields(3,2,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(2,1,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(3)*Fields(2,2,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(6)*Fields(3,1,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(7)*Fields(3,2,xi1,yi1,zi1)**2
  
  ForcesDispu(3,3,x0,y0,z0) = &
    ForcesDispu(3,3,x0,y0,z0) &
    +0.25*CoeffDispu(4)*Fields(1,1,x0,y0,z0)**2 &
    -0.25*CoeffDispu(4)*Fields(1,1,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(4)*Fields(1,1,x0,yi1,z0)**2 &
    -0.25*CoeffDispu(4)*Fields(1,1,x0,yi1,zi1)**2 &
    +0.25*CoeffDispu(4)*Fields(1,1,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(4)*Fields(1,1,xi1,y0,zi1)**2 &
    +0.25*CoeffDispu(4)*Fields(1,1,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(4)*Fields(1,1,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(6)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
    +0.25*CoeffDispu(7)*Fields(1,2,x0,y0,z0)**2 &
    -0.25*CoeffDispu(6)*Fields(1,1,x0,y0,zi1)*Fields(1,2,x0,y0,zi1) &
    -0.25*CoeffDispu(7)*Fields(1,2,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(6)*Fields(1,1,x0,yi1,z0)*Fields(1,2,x0,yi1,z0) &
    +0.25*CoeffDispu(7)*Fields(1,2,x0,yi1,z0)**2 &
    -0.25*CoeffDispu(6)*Fields(1,1,x0,yi1,zi1)*Fields(1,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(7)*Fields(1,2,x0,yi1,zi1)**2 &
    +0.25*CoeffDispu(6)*Fields(1,1,xi1,y0,z0)*Fields(1,2,xi1,y0,z0) &
    +0.25*CoeffDispu(7)*Fields(1,2,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(6)*Fields(1,1,xi1,y0,zi1)*Fields(1,2,xi1,y0,zi1) &
    -0.25*CoeffDispu(7)*Fields(1,2,xi1,y0,zi1)**2 &
    +0.25*CoeffDispu(6)*Fields(1,1,xi1,yi1,z0)*Fields(1,2,xi1,yi1,z0) &
    +0.25*CoeffDispu(7)*Fields(1,2,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(6)*Fields(1,1,xi1,yi1,zi1)*Fields(1,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(7)*Fields(1,2,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(4)*Fields(2,1,x0,y0,z0)**2 &
    -0.25*CoeffDispu(4)*Fields(2,1,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(4)*Fields(2,1,x0,yi1,z0)**2 &
    -0.25*CoeffDispu(4)*Fields(2,1,x0,yi1,zi1)**2 &
    +0.25*CoeffDispu(4)*Fields(2,1,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(4)*Fields(2,1,xi1,y0,zi1)**2 &
    +0.25*CoeffDispu(4)*Fields(2,1,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(4)*Fields(2,1,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(6)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +0.25*CoeffDispu(7)*Fields(2,2,x0,y0,z0)**2 &
    -0.25*CoeffDispu(6)*Fields(2,1,x0,y0,zi1)*Fields(2,2,x0,y0,zi1) &
    -0.25*CoeffDispu(7)*Fields(2,2,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(6)*Fields(2,1,x0,yi1,z0)*Fields(2,2,x0,yi1,z0) &
    +0.25*CoeffDispu(7)*Fields(2,2,x0,yi1,z0)**2 &
    -0.25*CoeffDispu(6)*Fields(2,1,x0,yi1,zi1)*Fields(2,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(7)*Fields(2,2,x0,yi1,zi1)**2 &
    +0.25*CoeffDispu(6)*Fields(2,1,xi1,y0,z0)*Fields(2,2,xi1,y0,z0) &
    +0.25*CoeffDispu(7)*Fields(2,2,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(6)*Fields(2,1,xi1,y0,zi1)*Fields(2,2,xi1,y0,zi1) &
    -0.25*CoeffDispu(7)*Fields(2,2,xi1,y0,zi1)**2 &
    +0.25*CoeffDispu(6)*Fields(2,1,xi1,yi1,z0)*Fields(2,2,xi1,yi1,z0) &
    +0.25*CoeffDispu(7)*Fields(2,2,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(6)*Fields(2,1,xi1,yi1,zi1)*Fields(2,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(7)*Fields(2,2,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
    +0.25*CoeffDispu(2)*Fields(1,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
    +0.25*CoeffDispu(1)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
    +0.25*CoeffDispu(2)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
    +0.25*CoeffDispu(5)*Fields(3,1,x0,y0,z0)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,x0,y0,zi1)*Fields(3,1,x0,y0,zi1) &
    +0.25*CoeffDispu(2)*Fields(1,2,x0,y0,zi1)*Fields(3,1,x0,y0,zi1) &
    +0.25*CoeffDispu(1)*Fields(2,1,x0,y0,zi1)*Fields(3,1,x0,y0,zi1) &
    +0.25*CoeffDispu(2)*Fields(2,2,x0,y0,zi1)*Fields(3,1,x0,y0,zi1) &
    -0.25*CoeffDispu(5)*Fields(3,1,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,x0,yi1,z0)*Fields(3,1,x0,yi1,z0) &
    +0.25*CoeffDispu(2)*Fields(1,2,x0,yi1,z0)*Fields(3,1,x0,yi1,z0) &
    -0.25*CoeffDispu(1)*Fields(2,1,x0,yi1,z0)*Fields(3,1,x0,yi1,z0) &
    -0.25*CoeffDispu(2)*Fields(2,2,x0,yi1,z0)*Fields(3,1,x0,yi1,z0) &
    +0.25*CoeffDispu(5)*Fields(3,1,x0,yi1,z0)**2 &
    +0.25*CoeffDispu(1)*Fields(1,1,x0,yi1,zi1)*Fields(3,1,x0,yi1,zi1) &
    +0.25*CoeffDispu(2)*Fields(1,2,x0,yi1,zi1)*Fields(3,1,x0,yi1,zi1) &
    -0.25*CoeffDispu(1)*Fields(2,1,x0,yi1,zi1)*Fields(3,1,x0,yi1,zi1) &
    -0.25*CoeffDispu(2)*Fields(2,2,x0,yi1,zi1)*Fields(3,1,x0,yi1,zi1) &
    -0.25*CoeffDispu(5)*Fields(3,1,x0,yi1,zi1)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,xi1,y0,z0)*Fields(3,1,xi1,y0,z0) &
    -0.25*CoeffDispu(2)*Fields(1,2,xi1,y0,z0)*Fields(3,1,xi1,y0,z0) &
    +0.25*CoeffDispu(1)*Fields(2,1,xi1,y0,z0)*Fields(3,1,xi1,y0,z0) &
    +0.25*CoeffDispu(2)*Fields(2,2,xi1,y0,z0)*Fields(3,1,xi1,y0,z0) &
    +0.25*CoeffDispu(5)*Fields(3,1,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,xi1,y0,zi1)*Fields(3,1,xi1,y0,zi1) &
    -0.25*CoeffDispu(2)*Fields(1,2,xi1,y0,zi1)*Fields(3,1,xi1,y0,zi1) &
    +0.25*CoeffDispu(1)*Fields(2,1,xi1,y0,zi1)*Fields(3,1,xi1,y0,zi1) &
    +0.25*CoeffDispu(2)*Fields(2,2,xi1,y0,zi1)*Fields(3,1,xi1,y0,zi1) &
    -0.25*CoeffDispu(5)*Fields(3,1,xi1,y0,zi1)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,xi1,yi1,z0)*Fields(3,1,xi1,yi1,z0) &
    -0.25*CoeffDispu(2)*Fields(1,2,xi1,yi1,z0)*Fields(3,1,xi1,yi1,z0) &
    -0.25*CoeffDispu(1)*Fields(2,1,xi1,yi1,z0)*Fields(3,1,xi1,yi1,z0) &
    -0.25*CoeffDispu(2)*Fields(2,2,xi1,yi1,z0)*Fields(3,1,xi1,yi1,z0) &
    +0.25*CoeffDispu(5)*Fields(3,1,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(1)*Fields(1,1,xi1,yi1,zi1)*Fields(3,1,xi1,yi1,zi1) &
    -0.25*CoeffDispu(2)*Fields(1,2,xi1,yi1,zi1)*Fields(3,1,xi1,yi1,zi1) &
    -0.25*CoeffDispu(1)*Fields(2,1,xi1,yi1,zi1)*Fields(3,1,xi1,yi1,zi1) &
    -0.25*CoeffDispu(2)*Fields(2,2,xi1,yi1,zi1)*Fields(3,1,xi1,yi1,zi1) &
    -0.25*CoeffDispu(5)*Fields(3,1,xi1,yi1,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +0.25*CoeffDispu(3)*Fields(1,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +0.25*CoeffDispu(2)*Fields(2,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +0.25*CoeffDispu(3)*Fields(2,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +0.25*CoeffDispu(8)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +0.25*CoeffDispu(9)*Fields(3,2,x0,y0,z0)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
    +0.25*CoeffDispu(3)*Fields(1,2,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
    +0.25*CoeffDispu(2)*Fields(2,1,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
    +0.25*CoeffDispu(3)*Fields(2,2,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
    -0.25*CoeffDispu(8)*Fields(3,1,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
    -0.25*CoeffDispu(9)*Fields(3,2,x0,y0,zi1)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
    +0.25*CoeffDispu(3)*Fields(1,2,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
    -0.25*CoeffDispu(2)*Fields(2,1,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
    -0.25*CoeffDispu(3)*Fields(2,2,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
    +0.25*CoeffDispu(8)*Fields(3,1,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
    +0.25*CoeffDispu(9)*Fields(3,2,x0,yi1,z0)**2 &
    +0.25*CoeffDispu(2)*Fields(1,1,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
    +0.25*CoeffDispu(3)*Fields(1,2,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(2)*Fields(2,1,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(3)*Fields(2,2,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(8)*Fields(3,1,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
    -0.25*CoeffDispu(9)*Fields(3,2,x0,yi1,zi1)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
    -0.25*CoeffDispu(3)*Fields(1,2,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
    +0.25*CoeffDispu(2)*Fields(2,1,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
    +0.25*CoeffDispu(3)*Fields(2,2,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
    +0.25*CoeffDispu(8)*Fields(3,1,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
    +0.25*CoeffDispu(9)*Fields(3,2,xi1,y0,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
    -0.25*CoeffDispu(3)*Fields(1,2,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
    +0.25*CoeffDispu(2)*Fields(2,1,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
    +0.25*CoeffDispu(3)*Fields(2,2,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
    -0.25*CoeffDispu(8)*Fields(3,1,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
    -0.25*CoeffDispu(9)*Fields(3,2,xi1,y0,zi1)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(3)*Fields(1,2,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(2)*Fields(2,1,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
    -0.25*CoeffDispu(3)*Fields(2,2,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
    +0.25*CoeffDispu(8)*Fields(3,1,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
    +0.25*CoeffDispu(9)*Fields(3,2,xi1,yi1,z0)**2 &
    -0.25*CoeffDispu(2)*Fields(1,1,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(3)*Fields(1,2,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(2)*Fields(2,1,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(3)*Fields(2,2,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(8)*Fields(3,1,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
    -0.25*CoeffDispu(9)*Fields(3,2,xi1,yi1,zi1)**2
  
  ForcesDispu(1,4,x0,y0,z0) = &
    ForcesDispu(1,4,x0,y0,z0) &
    +0.
  
  ForcesDispu(2,4,x0,y0,z0) = &
    ForcesDispu(2,4,x0,y0,z0) &
    +0.
  
  ForcesDispu(3,4,x0,y0,z0) = &
    ForcesDispu(3,4,x0,y0,z0) &
    +0.
  
  ForcesDispu(4,4,x0,y0,z0) = &
    ForcesDispu(4,4,x0,y0,z0) &
    +0.
  
  ForcesDispu(5,4,x0,y0,z0) = &
    ForcesDispu(5,4,x0,y0,z0) &
    +0.
  
  ForcesDispu(6,4,x0,y0,z0) = &
    ForcesDispu(6,4,x0,y0,z0) &
    +0.
  
  ForcesEpsDisp(1,1,x0,y0,z0) = &
    ForcesEpsDisp(1,1,x0,y0,z0) &
    -2.*CoeffEpsDisp(5)*e0ij(1,1)*Fields(1,1,x0,y0,z0) &
    -2.*CoeffEpsDisp(4)*e0ij(2,2)*Fields(1,1,x0,y0,z0) &
    -2.*CoeffEpsDisp(4)*e0ij(3,3)*Fields(1,1,x0,y0,z0) &
    -2.*CoeffEpsDisp(5)*euij(1,1)*Fields(1,1,x0,y0,z0) &
    -2.*CoeffEpsDisp(4)*euij(2,2)*Fields(1,1,x0,y0,z0) &
    -2.*CoeffEpsDisp(4)*euij(3,3)*Fields(1,1,x0,y0,z0) &
    -(CoeffEpsDisp(8)*e0ij(1,1)*Fields(1,2,x0,y0,z0)) &
    -(CoeffEpsDisp(6)*e0ij(2,2)*Fields(1,2,x0,y0,z0)) &
    -(CoeffEpsDisp(6)*e0ij(3,3)*Fields(1,2,x0,y0,z0)) &
    -(CoeffEpsDisp(8)*euij(1,1)*Fields(1,2,x0,y0,z0)) &
    -(CoeffEpsDisp(6)*euij(2,2)*Fields(1,2,x0,y0,z0)) &
    -(CoeffEpsDisp(6)*euij(3,3)*Fields(1,2,x0,y0,z0)) &
    -(CoeffEpsDisp(1)*e0ij(1,2)*Fields(2,1,x0,y0,z0)) &
    -(CoeffEpsDisp(1)*euij(1,2)*Fields(2,1,x0,y0,z0)) &
    -(CoeffEpsDisp(2)*e0ij(1,2)*Fields(2,2,x0,y0,z0)) &
    -(CoeffEpsDisp(2)*euij(1,2)*Fields(2,2,x0,y0,z0)) &
    -(CoeffEpsDisp(1)*e0ij(1,3)*Fields(3,1,x0,y0,z0)) &
    -(CoeffEpsDisp(1)*euij(1,3)*Fields(3,1,x0,y0,z0)) &
    -(CoeffEpsDisp(2)*e0ij(1,3)*Fields(3,2,x0,y0,z0)) &
    -(CoeffEpsDisp(2)*euij(1,3)*Fields(3,2,x0,y0,z0))
  
  ForcesEpsDisp(2,1,x0,y0,z0) = &
    ForcesEpsDisp(2,1,x0,y0,z0) &
    -(CoeffEpsDisp(1)*e0ij(1,2)*Fields(1,1,x0,y0,z0)) &
    -(CoeffEpsDisp(1)*euij(1,2)*Fields(1,1,x0,y0,z0)) &
    -(CoeffEpsDisp(2)*e0ij(1,2)*Fields(1,2,x0,y0,z0)) &
    -(CoeffEpsDisp(2)*euij(1,2)*Fields(1,2,x0,y0,z0)) &
    -2.*CoeffEpsDisp(4)*e0ij(1,1)*Fields(2,1,x0,y0,z0) &
    -2.*CoeffEpsDisp(5)*e0ij(2,2)*Fields(2,1,x0,y0,z0) &
    -2.*CoeffEpsDisp(4)*e0ij(3,3)*Fields(2,1,x0,y0,z0) &
    -2.*CoeffEpsDisp(4)*euij(1,1)*Fields(2,1,x0,y0,z0) &
    -2.*CoeffEpsDisp(5)*euij(2,2)*Fields(2,1,x0,y0,z0) &
    -2.*CoeffEpsDisp(4)*euij(3,3)*Fields(2,1,x0,y0,z0) &
    -(CoeffEpsDisp(6)*e0ij(1,1)*Fields(2,2,x0,y0,z0)) &
    -(CoeffEpsDisp(8)*e0ij(2,2)*Fields(2,2,x0,y0,z0)) &
    -(CoeffEpsDisp(6)*e0ij(3,3)*Fields(2,2,x0,y0,z0)) &
    -(CoeffEpsDisp(6)*euij(1,1)*Fields(2,2,x0,y0,z0)) &
    -(CoeffEpsDisp(8)*euij(2,2)*Fields(2,2,x0,y0,z0)) &
    -(CoeffEpsDisp(6)*euij(3,3)*Fields(2,2,x0,y0,z0)) &
    -(CoeffEpsDisp(1)*e0ij(2,3)*Fields(3,1,x0,y0,z0)) &
    -(CoeffEpsDisp(1)*euij(2,3)*Fields(3,1,x0,y0,z0)) &
    -(CoeffEpsDisp(2)*e0ij(2,3)*Fields(3,2,x0,y0,z0)) &
    -(CoeffEpsDisp(2)*euij(2,3)*Fields(3,2,x0,y0,z0))
  
  ForcesEpsDisp(3,1,x0,y0,z0) = &
    ForcesEpsDisp(3,1,x0,y0,z0) &
    -(CoeffEpsDisp(1)*e0ij(1,3)*Fields(1,1,x0,y0,z0)) &
    -(CoeffEpsDisp(1)*euij(1,3)*Fields(1,1,x0,y0,z0)) &
    -(CoeffEpsDisp(2)*e0ij(1,3)*Fields(1,2,x0,y0,z0)) &
    -(CoeffEpsDisp(2)*euij(1,3)*Fields(1,2,x0,y0,z0)) &
    -(CoeffEpsDisp(1)*e0ij(2,3)*Fields(2,1,x0,y0,z0)) &
    -(CoeffEpsDisp(1)*euij(2,3)*Fields(2,1,x0,y0,z0)) &
    -(CoeffEpsDisp(2)*e0ij(2,3)*Fields(2,2,x0,y0,z0)) &
    -(CoeffEpsDisp(2)*euij(2,3)*Fields(2,2,x0,y0,z0)) &
    -2.*CoeffEpsDisp(4)*e0ij(1,1)*Fields(3,1,x0,y0,z0) &
    -2.*CoeffEpsDisp(4)*e0ij(2,2)*Fields(3,1,x0,y0,z0) &
    -2.*CoeffEpsDisp(5)*e0ij(3,3)*Fields(3,1,x0,y0,z0) &
    -2.*CoeffEpsDisp(4)*euij(1,1)*Fields(3,1,x0,y0,z0) &
    -2.*CoeffEpsDisp(4)*euij(2,2)*Fields(3,1,x0,y0,z0) &
    -2.*CoeffEpsDisp(5)*euij(3,3)*Fields(3,1,x0,y0,z0) &
    -(CoeffEpsDisp(6)*e0ij(1,1)*Fields(3,2,x0,y0,z0)) &
    -(CoeffEpsDisp(6)*e0ij(2,2)*Fields(3,2,x0,y0,z0)) &
    -(CoeffEpsDisp(8)*e0ij(3,3)*Fields(3,2,x0,y0,z0)) &
    -(CoeffEpsDisp(6)*euij(1,1)*Fields(3,2,x0,y0,z0)) &
    -(CoeffEpsDisp(6)*euij(2,2)*Fields(3,2,x0,y0,z0)) &
    -(CoeffEpsDisp(8)*euij(3,3)*Fields(3,2,x0,y0,z0))
  
  ForcesEpsDisp(1,2,x0,y0,z0) = &
    ForcesEpsDisp(1,2,x0,y0,z0) &
    -(CoeffEpsDisp(8)*e0ij(1,1)*Fields(1,1,x0,y0,z0)) &
    -(CoeffEpsDisp(6)*e0ij(2,2)*Fields(1,1,x0,y0,z0)) &
    -(CoeffEpsDisp(6)*e0ij(3,3)*Fields(1,1,x0,y0,z0)) &
    -(CoeffEpsDisp(8)*euij(1,1)*Fields(1,1,x0,y0,z0)) &
    -(CoeffEpsDisp(6)*euij(2,2)*Fields(1,1,x0,y0,z0)) &
    -(CoeffEpsDisp(6)*euij(3,3)*Fields(1,1,x0,y0,z0)) &
    -2.*CoeffEpsDisp(9)*e0ij(1,1)*Fields(1,2,x0,y0,z0) &
    -2.*CoeffEpsDisp(7)*e0ij(2,2)*Fields(1,2,x0,y0,z0) &
    -2.*CoeffEpsDisp(7)*e0ij(3,3)*Fields(1,2,x0,y0,z0) &
    -2.*CoeffEpsDisp(9)*euij(1,1)*Fields(1,2,x0,y0,z0) &
    -2.*CoeffEpsDisp(7)*euij(2,2)*Fields(1,2,x0,y0,z0) &
    -2.*CoeffEpsDisp(7)*euij(3,3)*Fields(1,2,x0,y0,z0) &
    -(CoeffEpsDisp(2)*e0ij(1,2)*Fields(2,1,x0,y0,z0)) &
    -(CoeffEpsDisp(2)*euij(1,2)*Fields(2,1,x0,y0,z0)) &
    -(CoeffEpsDisp(3)*e0ij(1,2)*Fields(2,2,x0,y0,z0)) &
    -(CoeffEpsDisp(3)*euij(1,2)*Fields(2,2,x0,y0,z0)) &
    -(CoeffEpsDisp(2)*e0ij(1,3)*Fields(3,1,x0,y0,z0)) &
    -(CoeffEpsDisp(2)*euij(1,3)*Fields(3,1,x0,y0,z0)) &
    -(CoeffEpsDisp(3)*e0ij(1,3)*Fields(3,2,x0,y0,z0)) &
    -(CoeffEpsDisp(3)*euij(1,3)*Fields(3,2,x0,y0,z0))
  
  ForcesEpsDisp(2,2,x0,y0,z0) = &
    ForcesEpsDisp(2,2,x0,y0,z0) &
    -(CoeffEpsDisp(2)*e0ij(1,2)*Fields(1,1,x0,y0,z0)) &
    -(CoeffEpsDisp(2)*euij(1,2)*Fields(1,1,x0,y0,z0)) &
    -(CoeffEpsDisp(3)*e0ij(1,2)*Fields(1,2,x0,y0,z0)) &
    -(CoeffEpsDisp(3)*euij(1,2)*Fields(1,2,x0,y0,z0)) &
    -(CoeffEpsDisp(6)*e0ij(1,1)*Fields(2,1,x0,y0,z0)) &
    -(CoeffEpsDisp(8)*e0ij(2,2)*Fields(2,1,x0,y0,z0)) &
    -(CoeffEpsDisp(6)*e0ij(3,3)*Fields(2,1,x0,y0,z0)) &
    -(CoeffEpsDisp(6)*euij(1,1)*Fields(2,1,x0,y0,z0)) &
    -(CoeffEpsDisp(8)*euij(2,2)*Fields(2,1,x0,y0,z0)) &
    -(CoeffEpsDisp(6)*euij(3,3)*Fields(2,1,x0,y0,z0)) &
    -2.*CoeffEpsDisp(7)*e0ij(1,1)*Fields(2,2,x0,y0,z0) &
    -2.*CoeffEpsDisp(9)*e0ij(2,2)*Fields(2,2,x0,y0,z0) &
    -2.*CoeffEpsDisp(7)*e0ij(3,3)*Fields(2,2,x0,y0,z0) &
    -2.*CoeffEpsDisp(7)*euij(1,1)*Fields(2,2,x0,y0,z0) &
    -2.*CoeffEpsDisp(9)*euij(2,2)*Fields(2,2,x0,y0,z0) &
    -2.*CoeffEpsDisp(7)*euij(3,3)*Fields(2,2,x0,y0,z0) &
    -(CoeffEpsDisp(2)*e0ij(2,3)*Fields(3,1,x0,y0,z0)) &
    -(CoeffEpsDisp(2)*euij(2,3)*Fields(3,1,x0,y0,z0)) &
    -(CoeffEpsDisp(3)*e0ij(2,3)*Fields(3,2,x0,y0,z0)) &
    -(CoeffEpsDisp(3)*euij(2,3)*Fields(3,2,x0,y0,z0))
  
  ForcesEpsDisp(3,2,x0,y0,z0) = &
    ForcesEpsDisp(3,2,x0,y0,z0) &
    -(CoeffEpsDisp(2)*e0ij(1,3)*Fields(1,1,x0,y0,z0)) &
    -(CoeffEpsDisp(2)*euij(1,3)*Fields(1,1,x0,y0,z0)) &
    -(CoeffEpsDisp(3)*e0ij(1,3)*Fields(1,2,x0,y0,z0)) &
    -(CoeffEpsDisp(3)*euij(1,3)*Fields(1,2,x0,y0,z0)) &
    -(CoeffEpsDisp(2)*e0ij(2,3)*Fields(2,1,x0,y0,z0)) &
    -(CoeffEpsDisp(2)*euij(2,3)*Fields(2,1,x0,y0,z0)) &
    -(CoeffEpsDisp(3)*e0ij(2,3)*Fields(2,2,x0,y0,z0)) &
    -(CoeffEpsDisp(3)*euij(2,3)*Fields(2,2,x0,y0,z0)) &
    -(CoeffEpsDisp(6)*e0ij(1,1)*Fields(3,1,x0,y0,z0)) &
    -(CoeffEpsDisp(6)*e0ij(2,2)*Fields(3,1,x0,y0,z0)) &
    -(CoeffEpsDisp(8)*e0ij(3,3)*Fields(3,1,x0,y0,z0)) &
    -(CoeffEpsDisp(6)*euij(1,1)*Fields(3,1,x0,y0,z0)) &
    -(CoeffEpsDisp(6)*euij(2,2)*Fields(3,1,x0,y0,z0)) &
    -(CoeffEpsDisp(8)*euij(3,3)*Fields(3,1,x0,y0,z0)) &
    -2.*CoeffEpsDisp(7)*e0ij(1,1)*Fields(3,2,x0,y0,z0) &
    -2.*CoeffEpsDisp(7)*e0ij(2,2)*Fields(3,2,x0,y0,z0) &
    -2.*CoeffEpsDisp(9)*e0ij(3,3)*Fields(3,2,x0,y0,z0) &
    -2.*CoeffEpsDisp(7)*euij(1,1)*Fields(3,2,x0,y0,z0) &
    -2.*CoeffEpsDisp(7)*euij(2,2)*Fields(3,2,x0,y0,z0) &
    -2.*CoeffEpsDisp(9)*euij(3,3)*Fields(3,2,x0,y0,z0)
  
  ForcesEpsDisp(1,3,x0,y0,z0) = &
    ForcesEpsDisp(1,3,x0,y0,z0) &
    +0.
  
  ForcesEpsDisp(2,3,x0,y0,z0) = &
    ForcesEpsDisp(2,3,x0,y0,z0) &
    +0.
  
  ForcesEpsDisp(3,3,x0,y0,z0) = &
    ForcesEpsDisp(3,3,x0,y0,z0) &
    +0.
  
  ForcesEpsDisp(1,4,x0,y0,z0) = &
    ForcesEpsDisp(1,4,x0,y0,z0) &
    -(CoeffEpsDisp(5)*Fields(1,1,x0,y0,z0)**2) &
    -(CoeffEpsDisp(8)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)) &
    -(CoeffEpsDisp(9)*Fields(1,2,x0,y0,z0)**2) &
    -(CoeffEpsDisp(4)*Fields(2,1,x0,y0,z0)**2) &
    -(CoeffEpsDisp(6)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0)) &
    -(CoeffEpsDisp(7)*Fields(2,2,x0,y0,z0)**2) &
    -(CoeffEpsDisp(4)*Fields(3,1,x0,y0,z0)**2) &
    -(CoeffEpsDisp(6)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0)) &
    -(CoeffEpsDisp(7)*Fields(3,2,x0,y0,z0)**2)
  
  ForcesEpsDisp(2,4,x0,y0,z0) = &
    ForcesEpsDisp(2,4,x0,y0,z0) &
    -(CoeffEpsDisp(4)*Fields(1,1,x0,y0,z0)**2) &
    -(CoeffEpsDisp(6)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)) &
    -(CoeffEpsDisp(7)*Fields(1,2,x0,y0,z0)**2) &
    -(CoeffEpsDisp(5)*Fields(2,1,x0,y0,z0)**2) &
    -(CoeffEpsDisp(8)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0)) &
    -(CoeffEpsDisp(9)*Fields(2,2,x0,y0,z0)**2) &
    -(CoeffEpsDisp(4)*Fields(3,1,x0,y0,z0)**2) &
    -(CoeffEpsDisp(6)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0)) &
    -(CoeffEpsDisp(7)*Fields(3,2,x0,y0,z0)**2)
  
  ForcesEpsDisp(3,4,x0,y0,z0) = &
    ForcesEpsDisp(3,4,x0,y0,z0) &
    -(CoeffEpsDisp(4)*Fields(1,1,x0,y0,z0)**2) &
    -(CoeffEpsDisp(6)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)) &
    -(CoeffEpsDisp(7)*Fields(1,2,x0,y0,z0)**2) &
    -(CoeffEpsDisp(4)*Fields(2,1,x0,y0,z0)**2) &
    -(CoeffEpsDisp(6)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0)) &
    -(CoeffEpsDisp(7)*Fields(2,2,x0,y0,z0)**2) &
    -(CoeffEpsDisp(5)*Fields(3,1,x0,y0,z0)**2) &
    -(CoeffEpsDisp(8)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0)) &
    -(CoeffEpsDisp(9)*Fields(3,2,x0,y0,z0)**2)
  
  ForcesEpsDisp(4,4,x0,y0,z0) = &
    ForcesEpsDisp(4,4,x0,y0,z0) &
    -(CoeffEpsDisp(1)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,y0,z0)) &
    -(CoeffEpsDisp(2)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y0,z0)) &
    -(CoeffEpsDisp(2)*Fields(2,1,x0,y0,z0)*Fields(3,2,x0,y0,z0)) &
    -(CoeffEpsDisp(3)*Fields(2,2,x0,y0,z0)*Fields(3,2,x0,y0,z0))
  
  ForcesEpsDisp(5,4,x0,y0,z0) = &
    ForcesEpsDisp(5,4,x0,y0,z0) &
    -(CoeffEpsDisp(1)*Fields(1,1,x0,y0,z0)*Fields(3,1,x0,y0,z0)) &
    -(CoeffEpsDisp(2)*Fields(1,2,x0,y0,z0)*Fields(3,1,x0,y0,z0)) &
    -(CoeffEpsDisp(2)*Fields(1,1,x0,y0,z0)*Fields(3,2,x0,y0,z0)) &
    -(CoeffEpsDisp(3)*Fields(1,2,x0,y0,z0)*Fields(3,2,x0,y0,z0))
  
  ForcesEpsDisp(6,4,x0,y0,z0) = &
    ForcesEpsDisp(6,4,x0,y0,z0) &
    -(CoeffEpsDisp(1)*Fields(1,1,x0,y0,z0)*Fields(2,1,x0,y0,z0)) &
    -(CoeffEpsDisp(2)*Fields(1,2,x0,y0,z0)*Fields(2,1,x0,y0,z0)) &
    -(CoeffEpsDisp(2)*Fields(1,1,x0,y0,z0)*Fields(2,2,x0,y0,z0)) &
    -(CoeffEpsDisp(3)*Fields(1,2,x0,y0,z0)*Fields(2,2,x0,y0,z0))
  
  ForcesEpsu(1,1,x0,y0,z0) = &
    ForcesEpsu(1,1,x0,y0,z0) &
    +0.
  
  ForcesEpsu(2,1,x0,y0,z0) = &
    ForcesEpsu(2,1,x0,y0,z0) &
    +0.
  
  ForcesEpsu(3,1,x0,y0,z0) = &
    ForcesEpsu(3,1,x0,y0,z0) &
    +0.
  
  ForcesEpsu(1,2,x0,y0,z0) = &
    ForcesEpsu(1,2,x0,y0,z0) &
    +0.
  
  ForcesEpsu(2,2,x0,y0,z0) = &
    ForcesEpsu(2,2,x0,y0,z0) &
    +0.
  
  ForcesEpsu(3,2,x0,y0,z0) = &
    ForcesEpsu(3,2,x0,y0,z0) &
    +0.
  
  ForcesEpsu(1,3,x0,y0,z0) = &
    ForcesEpsu(1,3,x0,y0,z0) &
    +0.
  
  ForcesEpsu(2,3,x0,y0,z0) = &
    ForcesEpsu(2,3,x0,y0,z0) &
    +0.

  ForcesEpsu(3,3,x0,y0,z0) = &
    ForcesEpsu(3,3,x0,y0,z0) &
    +0.

  ForcesEpsu(1,4,x0,y0,z0) = &
    ForcesEpsu(1,4,x0,y0,z0) &
    +0.
  
  ForcesEpsu(2,4,x0,y0,z0) = &
    ForcesEpsu(2,4,x0,y0,z0) &
    +0.
  
  ForcesEpsu(3,4,x0,y0,z0) = &
    ForcesEpsu(3,4,x0,y0,z0) &
    +0.
  
  ForcesEpsu(4,4,x0,y0,z0) = &
    ForcesEpsu(4,4,x0,y0,z0) &
    +0.
  
  ForcesEpsu(5,4,x0,y0,z0) = &
    ForcesEpsu(5,4,x0,y0,z0) &
    +0.
  
  ForcesEpsu(6,4,x0,y0,z0) = &
    ForcesEpsu(6,4,x0,y0,z0) &
    +0.
  
  ForcesEps(1,1,x0,y0,z0) = &
    ForcesEps(1,1,x0,y0,z0) &
    +0.
  
  ForcesEps(2,1,x0,y0,z0) = &
    ForcesEps(2,1,x0,y0,z0) &
    +0.
  
  ForcesEps(3,1,x0,y0,z0) = &
    ForcesEps(3,1,x0,y0,z0) &
    +0.
  
  ForcesEps(1,2,x0,y0,z0) = &
    ForcesEps(1,2,x0,y0,z0) &
    +0.
  
  ForcesEps(2,2,x0,y0,z0) = &
    ForcesEps(2,2,x0,y0,z0) &
    +0.
  
  ForcesEps(3,2,x0,y0,z0) = &
    ForcesEps(3,2,x0,y0,z0) &
    +0.
  
  ForcesEps(1,3,x0,y0,z0) = &
    ForcesEps(1,3,x0,y0,z0) &
    +0.
  
  ForcesEps(2,3,x0,y0,z0) = &
    ForcesEps(2,3,x0,y0,z0) &
    +0.
  
  ForcesEps(3,3,x0,y0,z0) = &
    ForcesEps(3,3,x0,y0,z0) &
    +0.
  
  ForcesEps(1,4,x0,y0,z0) = &
    ForcesEps(1,4,x0,y0,z0) &
    -CoeffEps(1) &
    -(CoeffEps(4)*eij(1,1)) &
    -(CoeffEps(3)*eij(2,2)) &
    -(CoeffEps(3)*eij(3,3)) 
  
  ForcesEps(2,4,x0,y0,z0) = &
    ForcesEps(2,4,x0,y0,z0) &
    -CoeffEps(1) &
    -(CoeffEps(3)*eij(1,1)) &
    -(CoeffEps(4)*eij(2,2)) &
    -(CoeffEps(3)*eij(3,3)) 

  
  ForcesEps(3,4,x0,y0,z0) = &
    ForcesEps(3,4,x0,y0,z0) &
    -CoeffEps(1) &
    -(CoeffEps(3)*eij(1,1)) &
    -(CoeffEps(3)*eij(2,2)) &
    -(CoeffEps(4)*eij(3,3)) 
  
  ForcesEps(4,4,x0,y0,z0) = &
    ForcesEps(4,4,x0,y0,z0) &
    -(CoeffEps(2)*eij(2,3)) 
  
  ForcesEps(5,4,x0,y0,z0) = &
    ForcesEps(5,4,x0,y0,z0) &
    -(CoeffEps(2)*eij(1,3)) 
  
  ForcesEps(6,4,x0,y0,z0) = &
    ForcesEps(6,4,x0,y0,z0) &
    -(CoeffEps(2)*eij(1,2)) 
  
  end do
  end do
  end do
!  !$OMP    END DO
!  !$OMP    END PARALLEL
  Forces=ForcesDisp+ForcesJijsm+ForcesJijhf+Forcesu+ForcesDispu+ForcesEpsDisp+ForcesEpsu+ForcesEps
  
End Function GetForces
