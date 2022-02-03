Function GetDeltaH(x0, y0, z0, Fields, e0ij, idelta, delta) Result(DeltaH)
  
  Implicit none
  Integer, Intent(in) :: x0, y0, z0, idelta
  Real*8,  Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(in) :: delta(FieldDimList(idelta))
  Real*8,  Intent(in) :: e0ij(3,3)
  Real*8              :: eij(3,3), euij(3,3)
  Real*8              :: DeltaH
  Real*8              :: DeltaHDisp
  Real*8              :: DeltaHJijsm
  Real*8              :: DeltaHJijhf
  Real*8              :: DeltaHu
  Real*8              :: DeltaHDispu
  Real*8              :: DeltaHEpsDisp
  Real*8              :: DeltaHEpsu
  Real*8              :: DeltaHEps
  
  Integer             :: x1,y1,z1,xi1,yi1,zi1
  
  x1 = (x0+1)-floor(real(x0+1-1)/real(cgrid%n1))*cgrid%n1
  y1 = (y0+1)-floor(real(y0+1-1)/real(cgrid%n2))*cgrid%n2
  z1 = (z0+1)-floor(real(z0+1-1)/real(cgrid%n3))*cgrid%n3
  xi1 = (x0-1)-floor(real(x0-1-1)/real(cgrid%n1))*cgrid%n1
  yi1 = (y0-1)-floor(real(y0-1-1)/real(cgrid%n2))*cgrid%n2
  zi1 = (z0-1)-floor(real(z0-1-1)/real(cgrid%n3))*cgrid%n3
  
  DeltaH = 0.0D0
  DeltaHDisp = 0.0D0
  DeltaHJijsm = 0.0D0
  DeltaHJijhf = 0.0D0
  DeltaHu = 0.0D0
  DeltaHDispu = 0.0D0
  DeltaHEpsDisp = 0.0D0
  DeltaHEpsu = 0.0D0
  DeltaHEps = 0.0D0
  euij = StrainFromu(x0, y0, z0, Fields)
  eij = e0ij + euij
  
  SELECT CASE (idelta)
  CASE (1)
    DeltaHDisp = &
      DeltaHDisp &
      +CoeffDisp(1)*delta(1)**2 &
      +CoeffDisp(5)*delta(1)**4 &
      +CoeffDisp(1)*delta(2)**2 &
      +CoeffDisp(4)*delta(1)**2*delta(2)**2 &
      +CoeffDisp(5)*delta(2)**4 &
      +CoeffDisp(1)*delta(3)**2 &
      +CoeffDisp(4)*delta(1)**2*delta(3)**2 &
      +CoeffDisp(4)*delta(2)**2*delta(3)**2 &
      +CoeffDisp(5)*delta(3)**4 &
      +2.*CoeffDisp(1)*delta(1)*Fields(1,1,x0,y0,z0) &
      +4.*CoeffDisp(5)*delta(1)**3*Fields(1,1,x0,y0,z0) &
      +2.*CoeffDisp(4)*delta(1)*delta(2)**2*Fields(1,1,x0,y0,z0) &
      +2.*CoeffDisp(4)*delta(1)*delta(3)**2*Fields(1,1,x0,y0,z0) &
      +6.*CoeffDisp(5)*delta(1)**2*Fields(1,1,x0,y0,z0)**2 &
      +CoeffDisp(4)*delta(2)**2*Fields(1,1,x0,y0,z0)**2 &
      +CoeffDisp(4)*delta(3)**2*Fields(1,1,x0,y0,z0)**2 &
      +4.*CoeffDisp(5)*delta(1)*Fields(1,1,x0,y0,z0)**3 &
      +CoeffDisp(2)*delta(1)*Fields(1,2,x0,y0,z0) &
      +CoeffDisp(7)*delta(1)**3*Fields(1,2,x0,y0,z0) &
      +CoeffDisp(6)*delta(1)*delta(2)**2*Fields(1,2,x0,y0,z0) &
      +CoeffDisp(6)*delta(1)*delta(3)**2*Fields(1,2,x0,y0,z0) &
      +3.*CoeffDisp(7)*delta(1)**2*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
      +CoeffDisp(6)*delta(2)**2*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
      +CoeffDisp(6)*delta(3)**2*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
      +3.*CoeffDisp(7)*delta(1)*Fields(1,1,x0,y0,z0)**2*Fields(1,2,x0,y0,z0) &
      +CoeffDisp(10)*delta(1)**2*Fields(1,2,x0,y0,z0)**2 &
      +CoeffDisp(9)*delta(2)**2*Fields(1,2,x0,y0,z0)**2 &
      +CoeffDisp(9)*delta(3)**2*Fields(1,2,x0,y0,z0)**2 &
      +2.*CoeffDisp(10)*delta(1)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)**2 &
      +CoeffDisp(13)*delta(1)*Fields(1,2,x0,y0,z0)**3 &
      +2.*CoeffDisp(1)*delta(2)*Fields(2,1,x0,y0,z0) &
      +2.*CoeffDisp(4)*delta(1)**2*delta(2)*Fields(2,1,x0,y0,z0) &
      +4.*CoeffDisp(5)*delta(2)**3*Fields(2,1,x0,y0,z0) &
      +2.*CoeffDisp(4)*delta(2)*delta(3)**2*Fields(2,1,x0,y0,z0) &
      +4.*CoeffDisp(4)*delta(1)*delta(2)*Fields(1,1,x0,y0,z0)*Fields(2,1,x0,y0,z0) &
      +2.*CoeffDisp(4)*delta(2)*Fields(1,1,x0,y0,z0)**2*Fields(2,1,x0,y0,z0) &
      +2.*CoeffDisp(6)*delta(1)*delta(2)*Fields(1,2,x0,y0,z0)*Fields(2,1,x0,y0,z0) &
      +2.*CoeffDisp(6)*delta(2)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)*Fields(2,1,x0,y0,z0) &
      +2.*CoeffDisp(9)*delta(2)*Fields(1,2,x0,y0,z0)**2*Fields(2,1,x0,y0,z0) &
      +CoeffDisp(4)*delta(1)**2*Fields(2,1,x0,y0,z0)**2 &
      +6.*CoeffDisp(5)*delta(2)**2*Fields(2,1,x0,y0,z0)**2 &
      +CoeffDisp(4)*delta(3)**2*Fields(2,1,x0,y0,z0)**2 &
      +2.*CoeffDisp(4)*delta(1)*Fields(1,1,x0,y0,z0)*Fields(2,1,x0,y0,z0)**2 &
      +CoeffDisp(6)*delta(1)*Fields(1,2,x0,y0,z0)*Fields(2,1,x0,y0,z0)**2 &
      +4.*CoeffDisp(5)*delta(2)*Fields(2,1,x0,y0,z0)**3 &
      +CoeffDisp(2)*delta(2)*Fields(2,2,x0,y0,z0) &
      +CoeffDisp(6)*delta(1)**2*delta(2)*Fields(2,2,x0,y0,z0) &
      +CoeffDisp(7)*delta(2)**3*Fields(2,2,x0,y0,z0) &
      +CoeffDisp(6)*delta(2)*delta(3)**2*Fields(2,2,x0,y0,z0) &
      +2.*CoeffDisp(6)*delta(1)*delta(2)*Fields(1,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      +CoeffDisp(6)*delta(2)*Fields(1,1,x0,y0,z0)**2*Fields(2,2,x0,y0,z0) &
      +CoeffDisp(8)*delta(1)*delta(2)*Fields(1,2,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      +CoeffDisp(8)*delta(2)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      +CoeffDisp(11)*delta(2)*Fields(1,2,x0,y0,z0)**2*Fields(2,2,x0,y0,z0) &
      +CoeffDisp(6)*delta(1)**2*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      +3.*CoeffDisp(7)*delta(2)**2*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      +CoeffDisp(6)*delta(3)**2*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      +2.*CoeffDisp(6)*delta(1)*Fields(1,1,x0,y0,z0)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      +CoeffDisp(8)*delta(1)*Fields(1,2,x0,y0,z0)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      +3.*CoeffDisp(7)*delta(2)*Fields(2,1,x0,y0,z0)**2*Fields(2,2,x0,y0,z0) &
      +CoeffDisp(9)*delta(1)**2*Fields(2,2,x0,y0,z0)**2 &
      +CoeffDisp(10)*delta(2)**2*Fields(2,2,x0,y0,z0)**2 &
      +CoeffDisp(9)*delta(3)**2*Fields(2,2,x0,y0,z0)**2 &
      +2.*CoeffDisp(9)*delta(1)*Fields(1,1,x0,y0,z0)*Fields(2,2,x0,y0,z0)**2 &
      +CoeffDisp(11)*delta(1)*Fields(1,2,x0,y0,z0)*Fields(2,2,x0,y0,z0)**2 &
      +2.*CoeffDisp(10)*delta(2)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0)**2 &
      +CoeffDisp(13)*delta(2)*Fields(2,2,x0,y0,z0)**3 &
      +2.*CoeffDisp(1)*delta(3)*Fields(3,1,x0,y0,z0) &
      +2.*CoeffDisp(4)*delta(1)**2*delta(3)*Fields(3,1,x0,y0,z0) &
      +2.*CoeffDisp(4)*delta(2)**2*delta(3)*Fields(3,1,x0,y0,z0) &
      +4.*CoeffDisp(5)*delta(3)**3*Fields(3,1,x0,y0,z0) &
      +4.*CoeffDisp(4)*delta(1)*delta(3)*Fields(1,1,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      +2.*CoeffDisp(4)*delta(3)*Fields(1,1,x0,y0,z0)**2*Fields(3,1,x0,y0,z0) &
      +2.*CoeffDisp(6)*delta(1)*delta(3)*Fields(1,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      +2.*CoeffDisp(6)*delta(3)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      +2.*CoeffDisp(9)*delta(3)*Fields(1,2,x0,y0,z0)**2*Fields(3,1,x0,y0,z0) &
      +4.*CoeffDisp(4)*delta(2)*delta(3)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      +2.*CoeffDisp(4)*delta(3)*Fields(2,1,x0,y0,z0)**2*Fields(3,1,x0,y0,z0) &
      +2.*CoeffDisp(6)*delta(2)*delta(3)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      +2.*CoeffDisp(6)*delta(3)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      +2.*CoeffDisp(9)*delta(3)*Fields(2,2,x0,y0,z0)**2*Fields(3,1,x0,y0,z0) &
      +CoeffDisp(4)*delta(1)**2*Fields(3,1,x0,y0,z0)**2 &
      +CoeffDisp(4)*delta(2)**2*Fields(3,1,x0,y0,z0)**2 &
      +6.*CoeffDisp(5)*delta(3)**2*Fields(3,1,x0,y0,z0)**2 &
      +2.*CoeffDisp(4)*delta(1)*Fields(1,1,x0,y0,z0)*Fields(3,1,x0,y0,z0)**2 &
      +CoeffDisp(6)*delta(1)*Fields(1,2,x0,y0,z0)*Fields(3,1,x0,y0,z0)**2 &
      +2.*CoeffDisp(4)*delta(2)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,y0,z0)**2 &
      +CoeffDisp(6)*delta(2)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y0,z0)**2 &
      +4.*CoeffDisp(5)*delta(3)*Fields(3,1,x0,y0,z0)**3 &
      +CoeffDisp(2)*delta(3)*Fields(3,2,x0,y0,z0) &
      +CoeffDisp(6)*delta(1)**2*delta(3)*Fields(3,2,x0,y0,z0) &
      +CoeffDisp(6)*delta(2)**2*delta(3)*Fields(3,2,x0,y0,z0) &
      +CoeffDisp(7)*delta(3)**3*Fields(3,2,x0,y0,z0) &
      +2.*CoeffDisp(6)*delta(1)*delta(3)*Fields(1,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +CoeffDisp(6)*delta(3)*Fields(1,1,x0,y0,z0)**2*Fields(3,2,x0,y0,z0) &
      +CoeffDisp(8)*delta(1)*delta(3)*Fields(1,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +CoeffDisp(8)*delta(3)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +CoeffDisp(11)*delta(3)*Fields(1,2,x0,y0,z0)**2*Fields(3,2,x0,y0,z0) &
      +2.*CoeffDisp(6)*delta(2)*delta(3)*Fields(2,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +CoeffDisp(6)*delta(3)*Fields(2,1,x0,y0,z0)**2*Fields(3,2,x0,y0,z0) &
      +CoeffDisp(8)*delta(2)*delta(3)*Fields(2,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +CoeffDisp(8)*delta(3)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +CoeffDisp(11)*delta(3)*Fields(2,2,x0,y0,z0)**2*Fields(3,2,x0,y0,z0) &
      +CoeffDisp(6)*delta(1)**2*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +CoeffDisp(6)*delta(2)**2*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +3.*CoeffDisp(7)*delta(3)**2*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +2.*CoeffDisp(6)*delta(1)*Fields(1,1,x0,y0,z0)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +CoeffDisp(8)*delta(1)*Fields(1,2,x0,y0,z0)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +2.*CoeffDisp(6)*delta(2)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +CoeffDisp(8)*delta(2)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +3.*CoeffDisp(7)*delta(3)*Fields(3,1,x0,y0,z0)**2*Fields(3,2,x0,y0,z0) &
      +CoeffDisp(9)*delta(1)**2*Fields(3,2,x0,y0,z0)**2 &
      +CoeffDisp(9)*delta(2)**2*Fields(3,2,x0,y0,z0)**2 &
      +CoeffDisp(10)*delta(3)**2*Fields(3,2,x0,y0,z0)**2 &
      +2.*CoeffDisp(9)*delta(1)*Fields(1,1,x0,y0,z0)*Fields(3,2,x0,y0,z0)**2 &
      +CoeffDisp(11)*delta(1)*Fields(1,2,x0,y0,z0)*Fields(3,2,x0,y0,z0)**2 &
      +2.*CoeffDisp(9)*delta(2)*Fields(2,1,x0,y0,z0)*Fields(3,2,x0,y0,z0)**2 &
      +CoeffDisp(11)*delta(2)*Fields(2,2,x0,y0,z0)*Fields(3,2,x0,y0,z0)**2 &
      +2.*CoeffDisp(10)*delta(3)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0)**2 &
      +CoeffDisp(13)*delta(3)*Fields(3,2,x0,y0,z0)**3
    
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
    
    DeltaHJijhf = &
      DeltaHJijhf &
      +0.
    
    DeltaHu = &
      DeltaHu &
      +0.
    
    DeltaHDispu = &
      DeltaHDispu &
      +0.
    
    DeltaHEpsDisp = &
      DeltaHEpsDisp &
      +CoeffEpsDisp(5)*delta(1)**2*e0ij(1,1) &
      +CoeffEpsDisp(4)*delta(2)**2*e0ij(1,1) &
      +CoeffEpsDisp(4)*delta(3)**2*e0ij(1,1) &
      +CoeffEpsDisp(1)*delta(1)*delta(2)*e0ij(1,2) &
      +CoeffEpsDisp(1)*delta(1)*delta(3)*e0ij(1,3) &
      +CoeffEpsDisp(4)*delta(1)**2*e0ij(2,2) &
      +CoeffEpsDisp(5)*delta(2)**2*e0ij(2,2) &
      +CoeffEpsDisp(4)*delta(3)**2*e0ij(2,2) &
      +CoeffEpsDisp(1)*delta(2)*delta(3)*e0ij(2,3) &
      +CoeffEpsDisp(4)*delta(1)**2*e0ij(3,3) &
      +CoeffEpsDisp(4)*delta(2)**2*e0ij(3,3) &
      +CoeffEpsDisp(5)*delta(3)**2*e0ij(3,3) &
      +CoeffEpsDisp(5)*delta(1)**2*euij(1,1) &
      +CoeffEpsDisp(4)*delta(2)**2*euij(1,1) &
      +CoeffEpsDisp(4)*delta(3)**2*euij(1,1) &
      +CoeffEpsDisp(1)*delta(1)*delta(2)*euij(1,2) &
      +CoeffEpsDisp(1)*delta(1)*delta(3)*euij(1,3) &
      +CoeffEpsDisp(4)*delta(1)**2*euij(2,2) &
      +CoeffEpsDisp(5)*delta(2)**2*euij(2,2) &
      +CoeffEpsDisp(4)*delta(3)**2*euij(2,2) &
      +CoeffEpsDisp(1)*delta(2)*delta(3)*euij(2,3) &
      +CoeffEpsDisp(4)*delta(1)**2*euij(3,3) &
      +CoeffEpsDisp(4)*delta(2)**2*euij(3,3) &
      +CoeffEpsDisp(5)*delta(3)**2*euij(3,3) &
      +2.*CoeffEpsDisp(5)*delta(1)*e0ij(1,1)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(1)*delta(2)*e0ij(1,2)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(1)*delta(3)*e0ij(1,3)*Fields(1,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(4)*delta(1)*e0ij(2,2)*Fields(1,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(4)*delta(1)*e0ij(3,3)*Fields(1,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(5)*delta(1)*euij(1,1)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(1)*delta(2)*euij(1,2)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(1)*delta(3)*euij(1,3)*Fields(1,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(4)*delta(1)*euij(2,2)*Fields(1,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(4)*delta(1)*euij(3,3)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(8)*delta(1)*e0ij(1,1)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(2)*e0ij(1,2)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(3)*e0ij(1,3)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(1)*e0ij(2,2)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(1)*e0ij(3,3)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(8)*delta(1)*euij(1,1)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(2)*euij(1,2)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(3)*euij(1,3)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(1)*euij(2,2)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(1)*euij(3,3)*Fields(1,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(4)*delta(2)*e0ij(1,1)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(1)*delta(1)*e0ij(1,2)*Fields(2,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(5)*delta(2)*e0ij(2,2)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(1)*delta(3)*e0ij(2,3)*Fields(2,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(4)*delta(2)*e0ij(3,3)*Fields(2,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(4)*delta(2)*euij(1,1)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(1)*delta(1)*euij(1,2)*Fields(2,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(5)*delta(2)*euij(2,2)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(1)*delta(3)*euij(2,3)*Fields(2,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(4)*delta(2)*euij(3,3)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(2)*e0ij(1,1)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(1)*e0ij(1,2)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(8)*delta(2)*e0ij(2,2)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(3)*e0ij(2,3)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(2)*e0ij(3,3)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(2)*euij(1,1)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(1)*euij(1,2)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(8)*delta(2)*euij(2,2)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(3)*euij(2,3)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(2)*euij(3,3)*Fields(2,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(4)*delta(3)*e0ij(1,1)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(1)*delta(1)*e0ij(1,3)*Fields(3,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(4)*delta(3)*e0ij(2,2)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(1)*delta(2)*e0ij(2,3)*Fields(3,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(5)*delta(3)*e0ij(3,3)*Fields(3,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(4)*delta(3)*euij(1,1)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(1)*delta(1)*euij(1,3)*Fields(3,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(4)*delta(3)*euij(2,2)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(1)*delta(2)*euij(2,3)*Fields(3,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(5)*delta(3)*euij(3,3)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(3)*e0ij(1,1)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(1)*e0ij(1,3)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(3)*e0ij(2,2)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(2)*e0ij(2,3)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(8)*delta(3)*e0ij(3,3)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(3)*euij(1,1)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(1)*euij(1,3)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(3)*euij(2,2)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(2)*euij(2,3)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(8)*delta(3)*euij(3,3)*Fields(3,2,x0,y0,z0)
    
    DeltaHEpsu = &
      DeltaHEpsu &
      +0.
    
    DeltaHEps = &
      DeltaHEps &
      +0.
    
    DeltaH=DeltaHDisp+DeltaHJijsm+DeltaHJijhf+DeltaHu+DeltaHDispu+DeltaHEpsDisp+DeltaHEpsu+DeltaHEps
  
  CASE (2)
    DeltaHDisp = &
      DeltaHDisp &
      +CoeffDisp(3)*delta(1)**2 &
      +CoeffDisp(14)*delta(1)**4 &
      +CoeffDisp(3)*delta(2)**2 &
      +CoeffDisp(12)*delta(1)**2*delta(2)**2 &
      +CoeffDisp(14)*delta(2)**4 &
      +CoeffDisp(3)*delta(3)**2 &
      +CoeffDisp(12)*delta(1)**2*delta(3)**2 &
      +CoeffDisp(12)*delta(2)**2*delta(3)**2 &
      +CoeffDisp(14)*delta(3)**4 &
      +CoeffDisp(2)*delta(1)*Fields(1,1,x0,y0,z0) &
      +CoeffDisp(13)*delta(1)**3*Fields(1,1,x0,y0,z0) &
      +CoeffDisp(11)*delta(1)*delta(2)**2*Fields(1,1,x0,y0,z0) &
      +CoeffDisp(11)*delta(1)*delta(3)**2*Fields(1,1,x0,y0,z0) &
      +CoeffDisp(10)*delta(1)**2*Fields(1,1,x0,y0,z0)**2 &
      +CoeffDisp(9)*delta(2)**2*Fields(1,1,x0,y0,z0)**2 &
      +CoeffDisp(9)*delta(3)**2*Fields(1,1,x0,y0,z0)**2 &
      +CoeffDisp(7)*delta(1)*Fields(1,1,x0,y0,z0)**3 &
      +2.*CoeffDisp(3)*delta(1)*Fields(1,2,x0,y0,z0) &
      +4.*CoeffDisp(14)*delta(1)**3*Fields(1,2,x0,y0,z0) &
      +2.*CoeffDisp(12)*delta(1)*delta(2)**2*Fields(1,2,x0,y0,z0) &
      +2.*CoeffDisp(12)*delta(1)*delta(3)**2*Fields(1,2,x0,y0,z0) &
      +3.*CoeffDisp(13)*delta(1)**2*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
      +CoeffDisp(11)*delta(2)**2*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
      +CoeffDisp(11)*delta(3)**2*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
      +2.*CoeffDisp(10)*delta(1)*Fields(1,1,x0,y0,z0)**2*Fields(1,2,x0,y0,z0) &
      +6.*CoeffDisp(14)*delta(1)**2*Fields(1,2,x0,y0,z0)**2 &
      +CoeffDisp(12)*delta(2)**2*Fields(1,2,x0,y0,z0)**2 &
      +CoeffDisp(12)*delta(3)**2*Fields(1,2,x0,y0,z0)**2 &
      +3.*CoeffDisp(13)*delta(1)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)**2 &
      +4.*CoeffDisp(14)*delta(1)*Fields(1,2,x0,y0,z0)**3 &
      +CoeffDisp(2)*delta(2)*Fields(2,1,x0,y0,z0) &
      +CoeffDisp(11)*delta(1)**2*delta(2)*Fields(2,1,x0,y0,z0) &
      +CoeffDisp(13)*delta(2)**3*Fields(2,1,x0,y0,z0) &
      +CoeffDisp(11)*delta(2)*delta(3)**2*Fields(2,1,x0,y0,z0) &
      +CoeffDisp(8)*delta(1)*delta(2)*Fields(1,1,x0,y0,z0)*Fields(2,1,x0,y0,z0) &
      +CoeffDisp(6)*delta(2)*Fields(1,1,x0,y0,z0)**2*Fields(2,1,x0,y0,z0) &
      +2.*CoeffDisp(11)*delta(1)*delta(2)*Fields(1,2,x0,y0,z0)*Fields(2,1,x0,y0,z0) &
      +CoeffDisp(8)*delta(2)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)*Fields(2,1,x0,y0,z0) &
      +CoeffDisp(11)*delta(2)*Fields(1,2,x0,y0,z0)**2*Fields(2,1,x0,y0,z0) &
      +CoeffDisp(9)*delta(1)**2*Fields(2,1,x0,y0,z0)**2 &
      +CoeffDisp(10)*delta(2)**2*Fields(2,1,x0,y0,z0)**2 &
      +CoeffDisp(9)*delta(3)**2*Fields(2,1,x0,y0,z0)**2 &
      +CoeffDisp(6)*delta(1)*Fields(1,1,x0,y0,z0)*Fields(2,1,x0,y0,z0)**2 &
      +2.*CoeffDisp(9)*delta(1)*Fields(1,2,x0,y0,z0)*Fields(2,1,x0,y0,z0)**2 &
      +CoeffDisp(7)*delta(2)*Fields(2,1,x0,y0,z0)**3 &
      +2.*CoeffDisp(3)*delta(2)*Fields(2,2,x0,y0,z0) &
      +2.*CoeffDisp(12)*delta(1)**2*delta(2)*Fields(2,2,x0,y0,z0) &
      +4.*CoeffDisp(14)*delta(2)**3*Fields(2,2,x0,y0,z0) &
      +2.*CoeffDisp(12)*delta(2)*delta(3)**2*Fields(2,2,x0,y0,z0) &
      +2.*CoeffDisp(11)*delta(1)*delta(2)*Fields(1,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      +2.*CoeffDisp(9)*delta(2)*Fields(1,1,x0,y0,z0)**2*Fields(2,2,x0,y0,z0) &
      +4.*CoeffDisp(12)*delta(1)*delta(2)*Fields(1,2,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      +2.*CoeffDisp(11)*delta(2)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      +2.*CoeffDisp(12)*delta(2)*Fields(1,2,x0,y0,z0)**2*Fields(2,2,x0,y0,z0) &
      +CoeffDisp(11)*delta(1)**2*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      +3.*CoeffDisp(13)*delta(2)**2*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      +CoeffDisp(11)*delta(3)**2*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      +CoeffDisp(8)*delta(1)*Fields(1,1,x0,y0,z0)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      +2.*CoeffDisp(11)*delta(1)*Fields(1,2,x0,y0,z0)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      +2.*CoeffDisp(10)*delta(2)*Fields(2,1,x0,y0,z0)**2*Fields(2,2,x0,y0,z0) &
      +CoeffDisp(12)*delta(1)**2*Fields(2,2,x0,y0,z0)**2 &
      +6.*CoeffDisp(14)*delta(2)**2*Fields(2,2,x0,y0,z0)**2 &
      +CoeffDisp(12)*delta(3)**2*Fields(2,2,x0,y0,z0)**2 &
      +CoeffDisp(11)*delta(1)*Fields(1,1,x0,y0,z0)*Fields(2,2,x0,y0,z0)**2 &
      +2.*CoeffDisp(12)*delta(1)*Fields(1,2,x0,y0,z0)*Fields(2,2,x0,y0,z0)**2 &
      +3.*CoeffDisp(13)*delta(2)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0)**2 &
      +4.*CoeffDisp(14)*delta(2)*Fields(2,2,x0,y0,z0)**3 &
      +CoeffDisp(2)*delta(3)*Fields(3,1,x0,y0,z0) &
      +CoeffDisp(11)*delta(1)**2*delta(3)*Fields(3,1,x0,y0,z0) &
      +CoeffDisp(11)*delta(2)**2*delta(3)*Fields(3,1,x0,y0,z0) &
      +CoeffDisp(13)*delta(3)**3*Fields(3,1,x0,y0,z0) &
      +CoeffDisp(8)*delta(1)*delta(3)*Fields(1,1,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      +CoeffDisp(6)*delta(3)*Fields(1,1,x0,y0,z0)**2*Fields(3,1,x0,y0,z0) &
      +2.*CoeffDisp(11)*delta(1)*delta(3)*Fields(1,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      +CoeffDisp(8)*delta(3)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      +CoeffDisp(11)*delta(3)*Fields(1,2,x0,y0,z0)**2*Fields(3,1,x0,y0,z0) &
      +CoeffDisp(8)*delta(2)*delta(3)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      +CoeffDisp(6)*delta(3)*Fields(2,1,x0,y0,z0)**2*Fields(3,1,x0,y0,z0) &
      +2.*CoeffDisp(11)*delta(2)*delta(3)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      +CoeffDisp(8)*delta(3)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      +CoeffDisp(11)*delta(3)*Fields(2,2,x0,y0,z0)**2*Fields(3,1,x0,y0,z0) &
      +CoeffDisp(9)*delta(1)**2*Fields(3,1,x0,y0,z0)**2 &
      +CoeffDisp(9)*delta(2)**2*Fields(3,1,x0,y0,z0)**2 &
      +CoeffDisp(10)*delta(3)**2*Fields(3,1,x0,y0,z0)**2 &
      +CoeffDisp(6)*delta(1)*Fields(1,1,x0,y0,z0)*Fields(3,1,x0,y0,z0)**2 &
      +2.*CoeffDisp(9)*delta(1)*Fields(1,2,x0,y0,z0)*Fields(3,1,x0,y0,z0)**2 &
      +CoeffDisp(6)*delta(2)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,y0,z0)**2 &
      +2.*CoeffDisp(9)*delta(2)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y0,z0)**2 &
      +CoeffDisp(7)*delta(3)*Fields(3,1,x0,y0,z0)**3 &
      +2.*CoeffDisp(3)*delta(3)*Fields(3,2,x0,y0,z0) &
      +2.*CoeffDisp(12)*delta(1)**2*delta(3)*Fields(3,2,x0,y0,z0) &
      +2.*CoeffDisp(12)*delta(2)**2*delta(3)*Fields(3,2,x0,y0,z0) &
      +4.*CoeffDisp(14)*delta(3)**3*Fields(3,2,x0,y0,z0) &
      +2.*CoeffDisp(11)*delta(1)*delta(3)*Fields(1,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +2.*CoeffDisp(9)*delta(3)*Fields(1,1,x0,y0,z0)**2*Fields(3,2,x0,y0,z0) &
      +4.*CoeffDisp(12)*delta(1)*delta(3)*Fields(1,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +2.*CoeffDisp(11)*delta(3)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +2.*CoeffDisp(12)*delta(3)*Fields(1,2,x0,y0,z0)**2*Fields(3,2,x0,y0,z0) &
      +2.*CoeffDisp(11)*delta(2)*delta(3)*Fields(2,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +2.*CoeffDisp(9)*delta(3)*Fields(2,1,x0,y0,z0)**2*Fields(3,2,x0,y0,z0) &
      +4.*CoeffDisp(12)*delta(2)*delta(3)*Fields(2,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +2.*CoeffDisp(11)*delta(3)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +2.*CoeffDisp(12)*delta(3)*Fields(2,2,x0,y0,z0)**2*Fields(3,2,x0,y0,z0) &
      +CoeffDisp(11)*delta(1)**2*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +CoeffDisp(11)*delta(2)**2*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +3.*CoeffDisp(13)*delta(3)**2*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +CoeffDisp(8)*delta(1)*Fields(1,1,x0,y0,z0)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +2.*CoeffDisp(11)*delta(1)*Fields(1,2,x0,y0,z0)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +CoeffDisp(8)*delta(2)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +2.*CoeffDisp(11)*delta(2)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +2.*CoeffDisp(10)*delta(3)*Fields(3,1,x0,y0,z0)**2*Fields(3,2,x0,y0,z0) &
      +CoeffDisp(12)*delta(1)**2*Fields(3,2,x0,y0,z0)**2 &
      +CoeffDisp(12)*delta(2)**2*Fields(3,2,x0,y0,z0)**2 &
      +6.*CoeffDisp(14)*delta(3)**2*Fields(3,2,x0,y0,z0)**2 &
      +CoeffDisp(11)*delta(1)*Fields(1,1,x0,y0,z0)*Fields(3,2,x0,y0,z0)**2 &
      +2.*CoeffDisp(12)*delta(1)*Fields(1,2,x0,y0,z0)*Fields(3,2,x0,y0,z0)**2 &
      +CoeffDisp(11)*delta(2)*Fields(2,1,x0,y0,z0)*Fields(3,2,x0,y0,z0)**2 &
      +2.*CoeffDisp(12)*delta(2)*Fields(2,2,x0,y0,z0)*Fields(3,2,x0,y0,z0)**2 &
      +3.*CoeffDisp(13)*delta(3)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0)**2 &
      +4.*CoeffDisp(14)*delta(3)*Fields(3,2,x0,y0,z0)**3
    
    DeltaHJijsm = &
      DeltaHJijsm &
      +0.
    
    DeltaHJijhf = &
      DeltaHJijhf &
      +2.*CoeffJijhf(2)*delta(1)*Fields(1,2,x0,y0,z1) &
      +2.*CoeffJijhf(2)*delta(1)*Fields(1,2,x0,y0,zi1) &
      +2.*CoeffJijhf(2)*delta(1)*Fields(1,2,x0,y1,z0) &
      +2.*CoeffJijhf(4)*delta(1)*Fields(1,2,x0,y1,z1) &
      +2.*CoeffJijhf(4)*delta(1)*Fields(1,2,x0,y1,zi1) &
      +2.*CoeffJijhf(2)*delta(1)*Fields(1,2,x0,yi1,z0) &
      +2.*CoeffJijhf(4)*delta(1)*Fields(1,2,x0,yi1,z1) &
      +2.*CoeffJijhf(4)*delta(1)*Fields(1,2,x0,yi1,zi1) &
      +2.*CoeffJijhf(1)*delta(1)*Fields(1,2,x1,y0,z0) &
      +2.*CoeffJijhf(3)*delta(1)*Fields(1,2,x1,y0,z1) &
      +2.*CoeffJijhf(6)*delta(3)*Fields(1,2,x1,y0,z1) &
      +2.*CoeffJijhf(3)*delta(1)*Fields(1,2,x1,y0,zi1) &
      -2.*CoeffJijhf(6)*delta(3)*Fields(1,2,x1,y0,zi1) &
      +2.*CoeffJijhf(3)*delta(1)*Fields(1,2,x1,y1,z0) &
      +2.*CoeffJijhf(6)*delta(2)*Fields(1,2,x1,y1,z0) &
      +2.*CoeffJijhf(5)*delta(1)*Fields(1,2,x1,y1,z1) &
      +2.*CoeffJijhf(7)*delta(2)*Fields(1,2,x1,y1,z1) &
      +2.*CoeffJijhf(7)*delta(3)*Fields(1,2,x1,y1,z1) &
      +2.*CoeffJijhf(5)*delta(1)*Fields(1,2,x1,y1,zi1) &
      +2.*CoeffJijhf(7)*delta(2)*Fields(1,2,x1,y1,zi1) &
      -2.*CoeffJijhf(7)*delta(3)*Fields(1,2,x1,y1,zi1) &
      +2.*CoeffJijhf(3)*delta(1)*Fields(1,2,x1,yi1,z0) &
      -2.*CoeffJijhf(6)*delta(2)*Fields(1,2,x1,yi1,z0) &
      +2.*CoeffJijhf(5)*delta(1)*Fields(1,2,x1,yi1,z1) &
      -2.*CoeffJijhf(7)*delta(2)*Fields(1,2,x1,yi1,z1) &
      +2.*CoeffJijhf(7)*delta(3)*Fields(1,2,x1,yi1,z1) &
      +2.*CoeffJijhf(5)*delta(1)*Fields(1,2,x1,yi1,zi1) &
      -2.*CoeffJijhf(7)*delta(2)*Fields(1,2,x1,yi1,zi1) &
      -2.*CoeffJijhf(7)*delta(3)*Fields(1,2,x1,yi1,zi1) &
      +2.*CoeffJijhf(1)*delta(1)*Fields(1,2,xi1,y0,z0) &
      +2.*CoeffJijhf(3)*delta(1)*Fields(1,2,xi1,y0,z1) &
      -2.*CoeffJijhf(6)*delta(3)*Fields(1,2,xi1,y0,z1) &
      +2.*CoeffJijhf(3)*delta(1)*Fields(1,2,xi1,y0,zi1) &
      +2.*CoeffJijhf(6)*delta(3)*Fields(1,2,xi1,y0,zi1) &
      +2.*CoeffJijhf(3)*delta(1)*Fields(1,2,xi1,y1,z0) &
      -2.*CoeffJijhf(6)*delta(2)*Fields(1,2,xi1,y1,z0) &
      +2.*CoeffJijhf(5)*delta(1)*Fields(1,2,xi1,y1,z1) &
      -2.*CoeffJijhf(7)*delta(2)*Fields(1,2,xi1,y1,z1) &
      -2.*CoeffJijhf(7)*delta(3)*Fields(1,2,xi1,y1,z1) &
      +2.*CoeffJijhf(5)*delta(1)*Fields(1,2,xi1,y1,zi1) &
      -2.*CoeffJijhf(7)*delta(2)*Fields(1,2,xi1,y1,zi1) &
      +2.*CoeffJijhf(7)*delta(3)*Fields(1,2,xi1,y1,zi1) &
      +2.*CoeffJijhf(3)*delta(1)*Fields(1,2,xi1,yi1,z0) &
      +2.*CoeffJijhf(6)*delta(2)*Fields(1,2,xi1,yi1,z0) &
      +2.*CoeffJijhf(5)*delta(1)*Fields(1,2,xi1,yi1,z1) &
      +2.*CoeffJijhf(7)*delta(2)*Fields(1,2,xi1,yi1,z1) &
      -2.*CoeffJijhf(7)*delta(3)*Fields(1,2,xi1,yi1,z1) &
      +2.*CoeffJijhf(5)*delta(1)*Fields(1,2,xi1,yi1,zi1) &
      +2.*CoeffJijhf(7)*delta(2)*Fields(1,2,xi1,yi1,zi1) &
      +2.*CoeffJijhf(7)*delta(3)*Fields(1,2,xi1,yi1,zi1) &
      +2.*CoeffJijhf(2)*delta(2)*Fields(2,2,x0,y0,z1) &
      +2.*CoeffJijhf(2)*delta(2)*Fields(2,2,x0,y0,zi1) &
      +2.*CoeffJijhf(1)*delta(2)*Fields(2,2,x0,y1,z0) &
      +2.*CoeffJijhf(3)*delta(2)*Fields(2,2,x0,y1,z1) &
      +2.*CoeffJijhf(6)*delta(3)*Fields(2,2,x0,y1,z1) &
      +2.*CoeffJijhf(3)*delta(2)*Fields(2,2,x0,y1,zi1) &
      -2.*CoeffJijhf(6)*delta(3)*Fields(2,2,x0,y1,zi1) &
      +2.*CoeffJijhf(1)*delta(2)*Fields(2,2,x0,yi1,z0) &
      +2.*CoeffJijhf(3)*delta(2)*Fields(2,2,x0,yi1,z1) &
      -2.*CoeffJijhf(6)*delta(3)*Fields(2,2,x0,yi1,z1) &
      +2.*CoeffJijhf(3)*delta(2)*Fields(2,2,x0,yi1,zi1) &
      +2.*CoeffJijhf(6)*delta(3)*Fields(2,2,x0,yi1,zi1) &
      +2.*CoeffJijhf(2)*delta(2)*Fields(2,2,x1,y0,z0) &
      +2.*CoeffJijhf(4)*delta(2)*Fields(2,2,x1,y0,z1) &
      +2.*CoeffJijhf(4)*delta(2)*Fields(2,2,x1,y0,zi1) &
      +2.*CoeffJijhf(6)*delta(1)*Fields(2,2,x1,y1,z0) &
      +2.*CoeffJijhf(3)*delta(2)*Fields(2,2,x1,y1,z0) &
      +2.*CoeffJijhf(7)*delta(1)*Fields(2,2,x1,y1,z1) &
      +2.*CoeffJijhf(5)*delta(2)*Fields(2,2,x1,y1,z1) &
      +2.*CoeffJijhf(7)*delta(3)*Fields(2,2,x1,y1,z1) &
      +2.*CoeffJijhf(7)*delta(1)*Fields(2,2,x1,y1,zi1) &
      +2.*CoeffJijhf(5)*delta(2)*Fields(2,2,x1,y1,zi1) &
      -2.*CoeffJijhf(7)*delta(3)*Fields(2,2,x1,y1,zi1) &
      -2.*CoeffJijhf(6)*delta(1)*Fields(2,2,x1,yi1,z0) &
      +2.*CoeffJijhf(3)*delta(2)*Fields(2,2,x1,yi1,z0) &
      -2.*CoeffJijhf(7)*delta(1)*Fields(2,2,x1,yi1,z1) &
      +2.*CoeffJijhf(5)*delta(2)*Fields(2,2,x1,yi1,z1) &
      -2.*CoeffJijhf(7)*delta(3)*Fields(2,2,x1,yi1,z1) &
      -2.*CoeffJijhf(7)*delta(1)*Fields(2,2,x1,yi1,zi1) &
      +2.*CoeffJijhf(5)*delta(2)*Fields(2,2,x1,yi1,zi1) &
      +2.*CoeffJijhf(7)*delta(3)*Fields(2,2,x1,yi1,zi1) &
      +2.*CoeffJijhf(2)*delta(2)*Fields(2,2,xi1,y0,z0) &
      +2.*CoeffJijhf(4)*delta(2)*Fields(2,2,xi1,y0,z1) &
      +2.*CoeffJijhf(4)*delta(2)*Fields(2,2,xi1,y0,zi1) &
      -2.*CoeffJijhf(6)*delta(1)*Fields(2,2,xi1,y1,z0) &
      +2.*CoeffJijhf(3)*delta(2)*Fields(2,2,xi1,y1,z0) &
      -2.*CoeffJijhf(7)*delta(1)*Fields(2,2,xi1,y1,z1) &
      +2.*CoeffJijhf(5)*delta(2)*Fields(2,2,xi1,y1,z1) &
      +2.*CoeffJijhf(7)*delta(3)*Fields(2,2,xi1,y1,z1) &
      -2.*CoeffJijhf(7)*delta(1)*Fields(2,2,xi1,y1,zi1) &
      +2.*CoeffJijhf(5)*delta(2)*Fields(2,2,xi1,y1,zi1) &
      -2.*CoeffJijhf(7)*delta(3)*Fields(2,2,xi1,y1,zi1) &
      +2.*CoeffJijhf(6)*delta(1)*Fields(2,2,xi1,yi1,z0) &
      +2.*CoeffJijhf(3)*delta(2)*Fields(2,2,xi1,yi1,z0) &
      +2.*CoeffJijhf(7)*delta(1)*Fields(2,2,xi1,yi1,z1) &
      +2.*CoeffJijhf(5)*delta(2)*Fields(2,2,xi1,yi1,z1) &
      -2.*CoeffJijhf(7)*delta(3)*Fields(2,2,xi1,yi1,z1) &
      +2.*CoeffJijhf(7)*delta(1)*Fields(2,2,xi1,yi1,zi1) &
      +2.*CoeffJijhf(5)*delta(2)*Fields(2,2,xi1,yi1,zi1) &
      +2.*CoeffJijhf(7)*delta(3)*Fields(2,2,xi1,yi1,zi1) &
      +2.*CoeffJijhf(1)*delta(3)*Fields(3,2,x0,y0,z1) &
      +2.*CoeffJijhf(1)*delta(3)*Fields(3,2,x0,y0,zi1) &
      +2.*CoeffJijhf(2)*delta(3)*Fields(3,2,x0,y1,z0) &
      +2.*CoeffJijhf(6)*delta(2)*Fields(3,2,x0,y1,z1) &
      +2.*CoeffJijhf(3)*delta(3)*Fields(3,2,x0,y1,z1) &
      -2.*CoeffJijhf(6)*delta(2)*Fields(3,2,x0,y1,zi1) &
      +2.*CoeffJijhf(3)*delta(3)*Fields(3,2,x0,y1,zi1) &
      +2.*CoeffJijhf(2)*delta(3)*Fields(3,2,x0,yi1,z0) &
      -2.*CoeffJijhf(6)*delta(2)*Fields(3,2,x0,yi1,z1) &
      +2.*CoeffJijhf(3)*delta(3)*Fields(3,2,x0,yi1,z1) &
      +2.*CoeffJijhf(6)*delta(2)*Fields(3,2,x0,yi1,zi1) &
      +2.*CoeffJijhf(3)*delta(3)*Fields(3,2,x0,yi1,zi1) &
      +2.*CoeffJijhf(2)*delta(3)*Fields(3,2,x1,y0,z0) &
      +2.*CoeffJijhf(6)*delta(1)*Fields(3,2,x1,y0,z1) &
      +2.*CoeffJijhf(3)*delta(3)*Fields(3,2,x1,y0,z1) &
      -2.*CoeffJijhf(6)*delta(1)*Fields(3,2,x1,y0,zi1) &
      +2.*CoeffJijhf(3)*delta(3)*Fields(3,2,x1,y0,zi1) &
      +2.*CoeffJijhf(4)*delta(3)*Fields(3,2,x1,y1,z0) &
      +2.*CoeffJijhf(7)*delta(1)*Fields(3,2,x1,y1,z1) &
      +2.*CoeffJijhf(7)*delta(2)*Fields(3,2,x1,y1,z1) &
      +2.*CoeffJijhf(5)*delta(3)*Fields(3,2,x1,y1,z1) &
      -2.*CoeffJijhf(7)*delta(1)*Fields(3,2,x1,y1,zi1) &
      -2.*CoeffJijhf(7)*delta(2)*Fields(3,2,x1,y1,zi1) &
      +2.*CoeffJijhf(5)*delta(3)*Fields(3,2,x1,y1,zi1) &
      +2.*CoeffJijhf(4)*delta(3)*Fields(3,2,x1,yi1,z0) &
      +2.*CoeffJijhf(7)*delta(1)*Fields(3,2,x1,yi1,z1) &
      -2.*CoeffJijhf(7)*delta(2)*Fields(3,2,x1,yi1,z1) &
      +2.*CoeffJijhf(5)*delta(3)*Fields(3,2,x1,yi1,z1) &
      -2.*CoeffJijhf(7)*delta(1)*Fields(3,2,x1,yi1,zi1) &
      +2.*CoeffJijhf(7)*delta(2)*Fields(3,2,x1,yi1,zi1) &
      +2.*CoeffJijhf(5)*delta(3)*Fields(3,2,x1,yi1,zi1) &
      +2.*CoeffJijhf(2)*delta(3)*Fields(3,2,xi1,y0,z0) &
      -2.*CoeffJijhf(6)*delta(1)*Fields(3,2,xi1,y0,z1) &
      +2.*CoeffJijhf(3)*delta(3)*Fields(3,2,xi1,y0,z1) &
      +2.*CoeffJijhf(6)*delta(1)*Fields(3,2,xi1,y0,zi1) &
      +2.*CoeffJijhf(3)*delta(3)*Fields(3,2,xi1,y0,zi1) &
      +2.*CoeffJijhf(4)*delta(3)*Fields(3,2,xi1,y1,z0) &
      -2.*CoeffJijhf(7)*delta(1)*Fields(3,2,xi1,y1,z1) &
      +2.*CoeffJijhf(7)*delta(2)*Fields(3,2,xi1,y1,z1) &
      +2.*CoeffJijhf(5)*delta(3)*Fields(3,2,xi1,y1,z1) &
      +2.*CoeffJijhf(7)*delta(1)*Fields(3,2,xi1,y1,zi1) &
      -2.*CoeffJijhf(7)*delta(2)*Fields(3,2,xi1,y1,zi1) &
      +2.*CoeffJijhf(5)*delta(3)*Fields(3,2,xi1,y1,zi1) &
      +2.*CoeffJijhf(4)*delta(3)*Fields(3,2,xi1,yi1,z0) &
      -2.*CoeffJijhf(7)*delta(1)*Fields(3,2,xi1,yi1,z1) &
      -2.*CoeffJijhf(7)*delta(2)*Fields(3,2,xi1,yi1,z1) &
      +2.*CoeffJijhf(5)*delta(3)*Fields(3,2,xi1,yi1,z1) &
      +2.*CoeffJijhf(7)*delta(1)*Fields(3,2,xi1,yi1,zi1) &
      +2.*CoeffJijhf(7)*delta(2)*Fields(3,2,xi1,yi1,zi1) &
      +2.*CoeffJijhf(5)*delta(3)*Fields(3,2,xi1,yi1,zi1)
    
    DeltaHu = &
      DeltaHu &
      +0.
    
    DeltaHDispu = &
      DeltaHDispu &
      +0.
    
    DeltaHEpsDisp = &
      DeltaHEpsDisp &
      +CoeffEpsDisp(9)*delta(1)**2*e0ij(1,1) &
      +CoeffEpsDisp(7)*delta(2)**2*e0ij(1,1) &
      +CoeffEpsDisp(7)*delta(3)**2*e0ij(1,1) &
      +CoeffEpsDisp(3)*delta(1)*delta(2)*e0ij(1,2) &
      +CoeffEpsDisp(3)*delta(1)*delta(3)*e0ij(1,3) &
      +CoeffEpsDisp(7)*delta(1)**2*e0ij(2,2) &
      +CoeffEpsDisp(9)*delta(2)**2*e0ij(2,2) &
      +CoeffEpsDisp(7)*delta(3)**2*e0ij(2,2) &
      +CoeffEpsDisp(3)*delta(2)*delta(3)*e0ij(2,3) &
      +CoeffEpsDisp(7)*delta(1)**2*e0ij(3,3) &
      +CoeffEpsDisp(7)*delta(2)**2*e0ij(3,3) &
      +CoeffEpsDisp(9)*delta(3)**2*e0ij(3,3) &
      +CoeffEpsDisp(9)*delta(1)**2*euij(1,1) &
      +CoeffEpsDisp(7)*delta(2)**2*euij(1,1) &
      +CoeffEpsDisp(7)*delta(3)**2*euij(1,1) &
      +CoeffEpsDisp(3)*delta(1)*delta(2)*euij(1,2) &
      +CoeffEpsDisp(3)*delta(1)*delta(3)*euij(1,3) &
      +CoeffEpsDisp(7)*delta(1)**2*euij(2,2) &
      +CoeffEpsDisp(9)*delta(2)**2*euij(2,2) &
      +CoeffEpsDisp(7)*delta(3)**2*euij(2,2) &
      +CoeffEpsDisp(3)*delta(2)*delta(3)*euij(2,3) &
      +CoeffEpsDisp(7)*delta(1)**2*euij(3,3) &
      +CoeffEpsDisp(7)*delta(2)**2*euij(3,3) &
      +CoeffEpsDisp(9)*delta(3)**2*euij(3,3) &
      +CoeffEpsDisp(8)*delta(1)*e0ij(1,1)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(2)*e0ij(1,2)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(3)*e0ij(1,3)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(1)*e0ij(2,2)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(1)*e0ij(3,3)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(8)*delta(1)*euij(1,1)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(2)*euij(1,2)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(3)*euij(1,3)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(1)*euij(2,2)*Fields(1,1,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(1)*euij(3,3)*Fields(1,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(9)*delta(1)*e0ij(1,1)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(2)*e0ij(1,2)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(3)*e0ij(1,3)*Fields(1,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(7)*delta(1)*e0ij(2,2)*Fields(1,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(7)*delta(1)*e0ij(3,3)*Fields(1,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(9)*delta(1)*euij(1,1)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(2)*euij(1,2)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(3)*euij(1,3)*Fields(1,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(7)*delta(1)*euij(2,2)*Fields(1,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(7)*delta(1)*euij(3,3)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(2)*e0ij(1,1)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(1)*e0ij(1,2)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(8)*delta(2)*e0ij(2,2)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(3)*e0ij(2,3)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(2)*e0ij(3,3)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(2)*euij(1,1)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(1)*euij(1,2)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(8)*delta(2)*euij(2,2)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(3)*euij(2,3)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(2)*euij(3,3)*Fields(2,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(7)*delta(2)*e0ij(1,1)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(1)*e0ij(1,2)*Fields(2,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(9)*delta(2)*e0ij(2,2)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(3)*e0ij(2,3)*Fields(2,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(7)*delta(2)*e0ij(3,3)*Fields(2,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(7)*delta(2)*euij(1,1)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(1)*euij(1,2)*Fields(2,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(9)*delta(2)*euij(2,2)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(3)*euij(2,3)*Fields(2,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(7)*delta(2)*euij(3,3)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(3)*e0ij(1,1)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(1)*e0ij(1,3)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(3)*e0ij(2,2)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(2)*e0ij(2,3)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(8)*delta(3)*e0ij(3,3)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(3)*euij(1,1)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(1)*euij(1,3)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(3)*euij(2,2)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(2)*euij(2,3)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(8)*delta(3)*euij(3,3)*Fields(3,1,x0,y0,z0) &
      +2.*CoeffEpsDisp(7)*delta(3)*e0ij(1,1)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(1)*e0ij(1,3)*Fields(3,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(7)*delta(3)*e0ij(2,2)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(2)*e0ij(2,3)*Fields(3,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(9)*delta(3)*e0ij(3,3)*Fields(3,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(7)*delta(3)*euij(1,1)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(1)*euij(1,3)*Fields(3,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(7)*delta(3)*euij(2,2)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(2)*euij(2,3)*Fields(3,2,x0,y0,z0) &
      +2.*CoeffEpsDisp(9)*delta(3)*euij(3,3)*Fields(3,2,x0,y0,z0)
    
    DeltaHEpsu = &
      DeltaHEpsu &
      +0.
    
    DeltaHEps = &
      DeltaHEps &
      +0.
    
    DeltaH=DeltaHDisp+DeltaHJijsm+DeltaHJijhf+DeltaHu+DeltaHDispu+DeltaHEpsDisp+DeltaHEpsu+DeltaHEps
  
  CASE (3)
    DeltaHDisp = &
      DeltaHDisp &
      +0.
    
    DeltaHJijsm = &
      DeltaHJijsm &
      +0.
    
    DeltaHJijhf = &
      DeltaHJijhf &
      +0.
    
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
    
    DeltaHDispu = &
      DeltaHDispu &
      +0. &
      +0.25*CoeffDispu(9)*delta(1)*Fields(1,2,xi1,yi1,zi1)**2 &
      +0.25*CoeffDispu(7)*delta(2)*Fields(1,2,xi1,yi1,zi1)**2 &
      +0.25*CoeffDispu(7)*delta(3)*Fields(1,2,xi1,yi1,zi1)**2 &
      +0.25*CoeffDispu(9)*delta(1)*Fields(1,2,xi1,yi1,z0)**2 &
      +0.25*CoeffDispu(7)*delta(2)*Fields(1,2,xi1,yi1,z0)**2 &
      -0.25*CoeffDispu(7)*delta(3)*Fields(1,2,xi1,yi1,z0)**2 &
      +0.25*CoeffDispu(9)*delta(1)*Fields(1,2,xi1,y0,zi1)**2 &
      -0.25*CoeffDispu(7)*delta(2)*Fields(1,2,xi1,y0,zi1)**2 &
      +0.25*CoeffDispu(7)*delta(3)*Fields(1,2,xi1,y0,zi1)**2 &
      +0.25*CoeffDispu(9)*delta(1)*Fields(1,2,xi1,y0,z0)**2 &
      -0.25*CoeffDispu(7)*delta(2)*Fields(1,2,xi1,y0,z0)**2 &
      -0.25*CoeffDispu(7)*delta(3)*Fields(1,2,xi1,y0,z0)**2 &
      -0.25*CoeffDispu(9)*delta(1)*Fields(1,2,x0,yi1,zi1)**2 &
      +0.25*CoeffDispu(7)*delta(2)*Fields(1,2,x0,yi1,zi1)**2 &
      +0.25*CoeffDispu(7)*delta(3)*Fields(1,2,x0,yi1,zi1)**2 &
      -0.25*CoeffDispu(9)*delta(1)*Fields(1,2,x0,yi1,z0)**2 &
      +0.25*CoeffDispu(7)*delta(2)*Fields(1,2,x0,yi1,z0)**2 &
      -0.25*CoeffDispu(7)*delta(3)*Fields(1,2,x0,yi1,z0)**2 &
      -0.25*CoeffDispu(9)*delta(1)*Fields(1,2,x0,y0,zi1)**2 &
      -0.25*CoeffDispu(7)*delta(2)*Fields(1,2,x0,y0,zi1)**2 &
      +0.25*CoeffDispu(7)*delta(3)*Fields(1,2,x0,y0,zi1)**2 &
      -0.25*CoeffDispu(9)*delta(1)*Fields(1,2,x0,y0,z0)**2 &
      -0.25*CoeffDispu(7)*delta(2)*Fields(1,2,x0,y0,z0)**2 &
      -0.25*CoeffDispu(7)*delta(3)*Fields(1,2,x0,y0,z0)**2 &
      +0.25*CoeffDispu(3)*delta(1)*Fields(1,2,xi1,yi1,zi1)*Fields(2,2,xi1,yi1,zi1) &
      +0.25*CoeffDispu(3)*delta(2)*Fields(1,2,xi1,yi1,zi1)*Fields(2,2,xi1,yi1,zi1) &
      +0.25*CoeffDispu(7)*delta(1)*Fields(2,2,xi1,yi1,zi1)**2 &
      +0.25*CoeffDispu(9)*delta(2)*Fields(2,2,xi1,yi1,zi1)**2 &
      +0.25*CoeffDispu(7)*delta(3)*Fields(2,2,xi1,yi1,zi1)**2 &
      +0.25*CoeffDispu(3)*delta(1)*Fields(1,2,xi1,yi1,z0)*Fields(2,2,xi1,yi1,z0) &
      +0.25*CoeffDispu(3)*delta(2)*Fields(1,2,xi1,yi1,z0)*Fields(2,2,xi1,yi1,z0) &
      +0.25*CoeffDispu(7)*delta(1)*Fields(2,2,xi1,yi1,z0)**2 &
      +0.25*CoeffDispu(9)*delta(2)*Fields(2,2,xi1,yi1,z0)**2 &
      -0.25*CoeffDispu(7)*delta(3)*Fields(2,2,xi1,yi1,z0)**2 &
      -0.25*CoeffDispu(3)*delta(1)*Fields(1,2,xi1,y0,zi1)*Fields(2,2,xi1,y0,zi1) &
      +0.25*CoeffDispu(3)*delta(2)*Fields(1,2,xi1,y0,zi1)*Fields(2,2,xi1,y0,zi1) &
      +0.25*CoeffDispu(7)*delta(1)*Fields(2,2,xi1,y0,zi1)**2 &
      -0.25*CoeffDispu(9)*delta(2)*Fields(2,2,xi1,y0,zi1)**2 &
      +0.25*CoeffDispu(7)*delta(3)*Fields(2,2,xi1,y0,zi1)**2 &
      -0.25*CoeffDispu(3)*delta(1)*Fields(1,2,xi1,y0,z0)*Fields(2,2,xi1,y0,z0) &
      +0.25*CoeffDispu(3)*delta(2)*Fields(1,2,xi1,y0,z0)*Fields(2,2,xi1,y0,z0) &
      +0.25*CoeffDispu(7)*delta(1)*Fields(2,2,xi1,y0,z0)**2 &
      -0.25*CoeffDispu(9)*delta(2)*Fields(2,2,xi1,y0,z0)**2 &
      -0.25*CoeffDispu(7)*delta(3)*Fields(2,2,xi1,y0,z0)**2 &
      +0.25*CoeffDispu(3)*delta(1)*Fields(1,2,x0,yi1,zi1)*Fields(2,2,x0,yi1,zi1) &
      -0.25*CoeffDispu(3)*delta(2)*Fields(1,2,x0,yi1,zi1)*Fields(2,2,x0,yi1,zi1) &
      -0.25*CoeffDispu(7)*delta(1)*Fields(2,2,x0,yi1,zi1)**2 &
      +0.25*CoeffDispu(9)*delta(2)*Fields(2,2,x0,yi1,zi1)**2 &
      +0.25*CoeffDispu(7)*delta(3)*Fields(2,2,x0,yi1,zi1)**2 &
      +0.25*CoeffDispu(3)*delta(1)*Fields(1,2,x0,yi1,z0)*Fields(2,2,x0,yi1,z0) &
      -0.25*CoeffDispu(3)*delta(2)*Fields(1,2,x0,yi1,z0)*Fields(2,2,x0,yi1,z0) &
      -0.25*CoeffDispu(7)*delta(1)*Fields(2,2,x0,yi1,z0)**2 &
      +0.25*CoeffDispu(9)*delta(2)*Fields(2,2,x0,yi1,z0)**2 &
      -0.25*CoeffDispu(7)*delta(3)*Fields(2,2,x0,yi1,z0)**2 &
      -0.25*CoeffDispu(3)*delta(1)*Fields(1,2,x0,y0,zi1)*Fields(2,2,x0,y0,zi1) &
      -0.25*CoeffDispu(3)*delta(2)*Fields(1,2,x0,y0,zi1)*Fields(2,2,x0,y0,zi1) &
      -0.25*CoeffDispu(7)*delta(1)*Fields(2,2,x0,y0,zi1)**2 &
      -0.25*CoeffDispu(9)*delta(2)*Fields(2,2,x0,y0,zi1)**2 &
      +0.25*CoeffDispu(7)*delta(3)*Fields(2,2,x0,y0,zi1)**2 &
      -0.25*CoeffDispu(3)*delta(1)*Fields(1,2,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      -0.25*CoeffDispu(3)*delta(2)*Fields(1,2,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      -0.25*CoeffDispu(7)*delta(1)*Fields(2,2,x0,y0,z0)**2 &
      -0.25*CoeffDispu(9)*delta(2)*Fields(2,2,x0,y0,z0)**2 &
      -0.25*CoeffDispu(7)*delta(3)*Fields(2,2,x0,y0,z0)**2 &
      +0.25*CoeffDispu(3)*delta(1)*Fields(1,2,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
      +0.25*CoeffDispu(3)*delta(3)*Fields(1,2,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
      +0.25*CoeffDispu(3)*delta(2)*Fields(2,2,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
      +0.25*CoeffDispu(3)*delta(3)*Fields(2,2,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
      +0.25*CoeffDispu(7)*delta(1)*Fields(3,2,xi1,yi1,zi1)**2 &
      +0.25*CoeffDispu(7)*delta(2)*Fields(3,2,xi1,yi1,zi1)**2 &
      +0.25*CoeffDispu(9)*delta(3)*Fields(3,2,xi1,yi1,zi1)**2 &
      -0.25*CoeffDispu(3)*delta(1)*Fields(1,2,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
      +0.25*CoeffDispu(3)*delta(3)*Fields(1,2,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
      -0.25*CoeffDispu(3)*delta(2)*Fields(2,2,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
      +0.25*CoeffDispu(3)*delta(3)*Fields(2,2,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
      +0.25*CoeffDispu(7)*delta(1)*Fields(3,2,xi1,yi1,z0)**2 &
      +0.25*CoeffDispu(7)*delta(2)*Fields(3,2,xi1,yi1,z0)**2 &
      -0.25*CoeffDispu(9)*delta(3)*Fields(3,2,xi1,yi1,z0)**2 &
      +0.25*CoeffDispu(3)*delta(1)*Fields(1,2,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
      +0.25*CoeffDispu(3)*delta(3)*Fields(1,2,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
      +0.25*CoeffDispu(3)*delta(2)*Fields(2,2,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
      -0.25*CoeffDispu(3)*delta(3)*Fields(2,2,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
      +0.25*CoeffDispu(7)*delta(1)*Fields(3,2,xi1,y0,zi1)**2 &
      -0.25*CoeffDispu(7)*delta(2)*Fields(3,2,xi1,y0,zi1)**2 &
      +0.25*CoeffDispu(9)*delta(3)*Fields(3,2,xi1,y0,zi1)**2 &
      -0.25*CoeffDispu(3)*delta(1)*Fields(1,2,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
      +0.25*CoeffDispu(3)*delta(3)*Fields(1,2,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
      -0.25*CoeffDispu(3)*delta(2)*Fields(2,2,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
      -0.25*CoeffDispu(3)*delta(3)*Fields(2,2,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
      +0.25*CoeffDispu(7)*delta(1)*Fields(3,2,xi1,y0,z0)**2 &
      -0.25*CoeffDispu(7)*delta(2)*Fields(3,2,xi1,y0,z0)**2 &
      -0.25*CoeffDispu(9)*delta(3)*Fields(3,2,xi1,y0,z0)**2 &
      +0.25*CoeffDispu(3)*delta(1)*Fields(1,2,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
      -0.25*CoeffDispu(3)*delta(3)*Fields(1,2,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
      +0.25*CoeffDispu(3)*delta(2)*Fields(2,2,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
      +0.25*CoeffDispu(3)*delta(3)*Fields(2,2,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
      -0.25*CoeffDispu(7)*delta(1)*Fields(3,2,x0,yi1,zi1)**2 &
      +0.25*CoeffDispu(7)*delta(2)*Fields(3,2,x0,yi1,zi1)**2 &
      +0.25*CoeffDispu(9)*delta(3)*Fields(3,2,x0,yi1,zi1)**2 &
      -0.25*CoeffDispu(3)*delta(1)*Fields(1,2,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
      -0.25*CoeffDispu(3)*delta(3)*Fields(1,2,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
      -0.25*CoeffDispu(3)*delta(2)*Fields(2,2,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
      +0.25*CoeffDispu(3)*delta(3)*Fields(2,2,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
      -0.25*CoeffDispu(7)*delta(1)*Fields(3,2,x0,yi1,z0)**2 &
      +0.25*CoeffDispu(7)*delta(2)*Fields(3,2,x0,yi1,z0)**2 &
      -0.25*CoeffDispu(9)*delta(3)*Fields(3,2,x0,yi1,z0)**2 &
      +0.25*CoeffDispu(3)*delta(1)*Fields(1,2,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
      -0.25*CoeffDispu(3)*delta(3)*Fields(1,2,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
      +0.25*CoeffDispu(3)*delta(2)*Fields(2,2,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
      -0.25*CoeffDispu(3)*delta(3)*Fields(2,2,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
      -0.25*CoeffDispu(7)*delta(1)*Fields(3,2,x0,y0,zi1)**2 &
      -0.25*CoeffDispu(7)*delta(2)*Fields(3,2,x0,y0,zi1)**2 &
      +0.25*CoeffDispu(9)*delta(3)*Fields(3,2,x0,y0,zi1)**2 &
      -0.25*CoeffDispu(3)*delta(1)*Fields(1,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      -0.25*CoeffDispu(3)*delta(3)*Fields(1,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      -0.25*CoeffDispu(3)*delta(2)*Fields(2,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      -0.25*CoeffDispu(3)*delta(3)*Fields(2,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      -0.25*CoeffDispu(7)*delta(1)*Fields(3,2,x0,y0,z0)**2 &
      -0.25*CoeffDispu(7)*delta(2)*Fields(3,2,x0,y0,z0)**2 &
      -0.25*CoeffDispu(9)*delta(3)*Fields(3,2,x0,y0,z0)**2 &
      +0.25*CoeffDispu(8)*delta(1)*Fields(1,1,xi1,yi1,zi1)*Fields(1,2,xi1,yi1,zi1) &
      +0.25*CoeffDispu(6)*delta(2)*Fields(1,1,xi1,yi1,zi1)*Fields(1,2,xi1,yi1,zi1) &
      +0.25*CoeffDispu(6)*delta(3)*Fields(1,1,xi1,yi1,zi1)*Fields(1,2,xi1,yi1,zi1) &
      +0.25*CoeffDispu(2)*delta(1)*Fields(1,1,xi1,yi1,zi1)*Fields(2,2,xi1,yi1,zi1) &
      +0.25*CoeffDispu(2)*delta(2)*Fields(1,1,xi1,yi1,zi1)*Fields(2,2,xi1,yi1,zi1) &
      +0.25*CoeffDispu(2)*delta(1)*Fields(1,1,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
      +0.25*CoeffDispu(2)*delta(3)*Fields(1,1,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
      +0.25*CoeffDispu(5)*delta(1)*Fields(1,1,xi1,yi1,zi1)**2 &
      +0.25*CoeffDispu(4)*delta(2)*Fields(1,1,xi1,yi1,zi1)**2 &
      +0.25*CoeffDispu(4)*delta(3)*Fields(1,1,xi1,yi1,zi1)**2 &
      +0.25*CoeffDispu(8)*delta(1)*Fields(1,1,xi1,yi1,z0)*Fields(1,2,xi1,yi1,z0) &
      +0.25*CoeffDispu(6)*delta(2)*Fields(1,1,xi1,yi1,z0)*Fields(1,2,xi1,yi1,z0) &
      -0.25*CoeffDispu(6)*delta(3)*Fields(1,1,xi1,yi1,z0)*Fields(1,2,xi1,yi1,z0) &
      +0.25*CoeffDispu(2)*delta(1)*Fields(1,1,xi1,yi1,z0)*Fields(2,2,xi1,yi1,z0) &
      +0.25*CoeffDispu(2)*delta(2)*Fields(1,1,xi1,yi1,z0)*Fields(2,2,xi1,yi1,z0) &
      -0.25*CoeffDispu(2)*delta(1)*Fields(1,1,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
      +0.25*CoeffDispu(2)*delta(3)*Fields(1,1,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
      +0.25*CoeffDispu(5)*delta(1)*Fields(1,1,xi1,yi1,z0)**2 &
      +0.25*CoeffDispu(4)*delta(2)*Fields(1,1,xi1,yi1,z0)**2 &
      -0.25*CoeffDispu(4)*delta(3)*Fields(1,1,xi1,yi1,z0)**2 &
      +0.25*CoeffDispu(8)*delta(1)*Fields(1,1,xi1,y0,zi1)*Fields(1,2,xi1,y0,zi1) &
      -0.25*CoeffDispu(6)*delta(2)*Fields(1,1,xi1,y0,zi1)*Fields(1,2,xi1,y0,zi1) &
      +0.25*CoeffDispu(6)*delta(3)*Fields(1,1,xi1,y0,zi1)*Fields(1,2,xi1,y0,zi1) &
      -0.25*CoeffDispu(2)*delta(1)*Fields(1,1,xi1,y0,zi1)*Fields(2,2,xi1,y0,zi1) &
      +0.25*CoeffDispu(2)*delta(2)*Fields(1,1,xi1,y0,zi1)*Fields(2,2,xi1,y0,zi1) &
      +0.25*CoeffDispu(2)*delta(1)*Fields(1,1,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
      +0.25*CoeffDispu(2)*delta(3)*Fields(1,1,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
      +0.25*CoeffDispu(5)*delta(1)*Fields(1,1,xi1,y0,zi1)**2 &
      -0.25*CoeffDispu(4)*delta(2)*Fields(1,1,xi1,y0,zi1)**2 &
      +0.25*CoeffDispu(4)*delta(3)*Fields(1,1,xi1,y0,zi1)**2 &
      +0.25*CoeffDispu(8)*delta(1)*Fields(1,1,xi1,y0,z0)*Fields(1,2,xi1,y0,z0) &
      -0.25*CoeffDispu(6)*delta(2)*Fields(1,1,xi1,y0,z0)*Fields(1,2,xi1,y0,z0) &
      -0.25*CoeffDispu(6)*delta(3)*Fields(1,1,xi1,y0,z0)*Fields(1,2,xi1,y0,z0) &
      -0.25*CoeffDispu(2)*delta(1)*Fields(1,1,xi1,y0,z0)*Fields(2,2,xi1,y0,z0) &
      +0.25*CoeffDispu(2)*delta(2)*Fields(1,1,xi1,y0,z0)*Fields(2,2,xi1,y0,z0) &
      -0.25*CoeffDispu(2)*delta(1)*Fields(1,1,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
      +0.25*CoeffDispu(2)*delta(3)*Fields(1,1,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
      +0.25*CoeffDispu(5)*delta(1)*Fields(1,1,xi1,y0,z0)**2 &
      -0.25*CoeffDispu(4)*delta(2)*Fields(1,1,xi1,y0,z0)**2 &
      -0.25*CoeffDispu(4)*delta(3)*Fields(1,1,xi1,y0,z0)**2 &
      -0.25*CoeffDispu(8)*delta(1)*Fields(1,1,x0,yi1,zi1)*Fields(1,2,x0,yi1,zi1) &
      +0.25*CoeffDispu(6)*delta(2)*Fields(1,1,x0,yi1,zi1)*Fields(1,2,x0,yi1,zi1) &
      +0.25*CoeffDispu(6)*delta(3)*Fields(1,1,x0,yi1,zi1)*Fields(1,2,x0,yi1,zi1) &
      +0.25*CoeffDispu(2)*delta(1)*Fields(1,1,x0,yi1,zi1)*Fields(2,2,x0,yi1,zi1) &
      -0.25*CoeffDispu(2)*delta(2)*Fields(1,1,x0,yi1,zi1)*Fields(2,2,x0,yi1,zi1) &
      +0.25*CoeffDispu(2)*delta(1)*Fields(1,1,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
      -0.25*CoeffDispu(2)*delta(3)*Fields(1,1,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
      -0.25*CoeffDispu(5)*delta(1)*Fields(1,1,x0,yi1,zi1)**2 &
      +0.25*CoeffDispu(4)*delta(2)*Fields(1,1,x0,yi1,zi1)**2 &
      +0.25*CoeffDispu(4)*delta(3)*Fields(1,1,x0,yi1,zi1)**2 &
      -0.25*CoeffDispu(8)*delta(1)*Fields(1,1,x0,yi1,z0)*Fields(1,2,x0,yi1,z0) &
      +0.25*CoeffDispu(6)*delta(2)*Fields(1,1,x0,yi1,z0)*Fields(1,2,x0,yi1,z0) &
      -0.25*CoeffDispu(6)*delta(3)*Fields(1,1,x0,yi1,z0)*Fields(1,2,x0,yi1,z0) &
      +0.25*CoeffDispu(2)*delta(1)*Fields(1,1,x0,yi1,z0)*Fields(2,2,x0,yi1,z0) &
      -0.25*CoeffDispu(2)*delta(2)*Fields(1,1,x0,yi1,z0)*Fields(2,2,x0,yi1,z0) &
      -0.25*CoeffDispu(2)*delta(1)*Fields(1,1,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
      -0.25*CoeffDispu(2)*delta(3)*Fields(1,1,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
      -0.25*CoeffDispu(5)*delta(1)*Fields(1,1,x0,yi1,z0)**2 &
      +0.25*CoeffDispu(4)*delta(2)*Fields(1,1,x0,yi1,z0)**2 &
      -0.25*CoeffDispu(4)*delta(3)*Fields(1,1,x0,yi1,z0)**2 &
      -0.25*CoeffDispu(8)*delta(1)*Fields(1,1,x0,y0,zi1)*Fields(1,2,x0,y0,zi1) &
      -0.25*CoeffDispu(6)*delta(2)*Fields(1,1,x0,y0,zi1)*Fields(1,2,x0,y0,zi1) &
      +0.25*CoeffDispu(6)*delta(3)*Fields(1,1,x0,y0,zi1)*Fields(1,2,x0,y0,zi1) &
      -0.25*CoeffDispu(2)*delta(1)*Fields(1,1,x0,y0,zi1)*Fields(2,2,x0,y0,zi1) &
      -0.25*CoeffDispu(2)*delta(2)*Fields(1,1,x0,y0,zi1)*Fields(2,2,x0,y0,zi1) &
      +0.25*CoeffDispu(2)*delta(1)*Fields(1,1,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
      -0.25*CoeffDispu(2)*delta(3)*Fields(1,1,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
      -0.25*CoeffDispu(5)*delta(1)*Fields(1,1,x0,y0,zi1)**2 &
      -0.25*CoeffDispu(4)*delta(2)*Fields(1,1,x0,y0,zi1)**2 &
      +0.25*CoeffDispu(4)*delta(3)*Fields(1,1,x0,y0,zi1)**2 &
      -0.25*CoeffDispu(8)*delta(1)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
      -0.25*CoeffDispu(6)*delta(2)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
      -0.25*CoeffDispu(6)*delta(3)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
      -0.25*CoeffDispu(2)*delta(1)*Fields(1,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      -0.25*CoeffDispu(2)*delta(2)*Fields(1,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      -0.25*CoeffDispu(2)*delta(1)*Fields(1,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      -0.25*CoeffDispu(2)*delta(3)*Fields(1,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      -0.25*CoeffDispu(5)*delta(1)*Fields(1,1,x0,y0,z0)**2 &
      -0.25*CoeffDispu(4)*delta(2)*Fields(1,1,x0,y0,z0)**2 &
      -0.25*CoeffDispu(4)*delta(3)*Fields(1,1,x0,y0,z0)**2 &
      +0.25*CoeffDispu(2)*delta(1)*Fields(1,2,xi1,yi1,zi1)*Fields(2,1,xi1,yi1,zi1) &
      +0.25*CoeffDispu(2)*delta(2)*Fields(1,2,xi1,yi1,zi1)*Fields(2,1,xi1,yi1,zi1) &
      +0.25*CoeffDispu(6)*delta(1)*Fields(2,1,xi1,yi1,zi1)*Fields(2,2,xi1,yi1,zi1) &
      +0.25*CoeffDispu(8)*delta(2)*Fields(2,1,xi1,yi1,zi1)*Fields(2,2,xi1,yi1,zi1) &
      +0.25*CoeffDispu(6)*delta(3)*Fields(2,1,xi1,yi1,zi1)*Fields(2,2,xi1,yi1,zi1) &
      +0.25*CoeffDispu(2)*delta(2)*Fields(2,1,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
      +0.25*CoeffDispu(2)*delta(3)*Fields(2,1,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
      +0.25*CoeffDispu(1)*delta(1)*Fields(1,1,xi1,yi1,zi1)*Fields(2,1,xi1,yi1,zi1) &
      +0.25*CoeffDispu(1)*delta(2)*Fields(1,1,xi1,yi1,zi1)*Fields(2,1,xi1,yi1,zi1) &
      +0.25*CoeffDispu(4)*delta(1)*Fields(2,1,xi1,yi1,zi1)**2 &
      +0.25*CoeffDispu(5)*delta(2)*Fields(2,1,xi1,yi1,zi1)**2 &
      +0.25*CoeffDispu(4)*delta(3)*Fields(2,1,xi1,yi1,zi1)**2 &
      +0.25*CoeffDispu(2)*delta(1)*Fields(1,2,xi1,yi1,z0)*Fields(2,1,xi1,yi1,z0) &
      +0.25*CoeffDispu(2)*delta(2)*Fields(1,2,xi1,yi1,z0)*Fields(2,1,xi1,yi1,z0) &
      +0.25*CoeffDispu(6)*delta(1)*Fields(2,1,xi1,yi1,z0)*Fields(2,2,xi1,yi1,z0) &
      +0.25*CoeffDispu(8)*delta(2)*Fields(2,1,xi1,yi1,z0)*Fields(2,2,xi1,yi1,z0) &
      -0.25*CoeffDispu(6)*delta(3)*Fields(2,1,xi1,yi1,z0)*Fields(2,2,xi1,yi1,z0) &
      -0.25*CoeffDispu(2)*delta(2)*Fields(2,1,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
      +0.25*CoeffDispu(2)*delta(3)*Fields(2,1,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
      +0.25*CoeffDispu(1)*delta(1)*Fields(1,1,xi1,yi1,z0)*Fields(2,1,xi1,yi1,z0) &
      +0.25*CoeffDispu(1)*delta(2)*Fields(1,1,xi1,yi1,z0)*Fields(2,1,xi1,yi1,z0) &
      +0.25*CoeffDispu(4)*delta(1)*Fields(2,1,xi1,yi1,z0)**2 &
      +0.25*CoeffDispu(5)*delta(2)*Fields(2,1,xi1,yi1,z0)**2 &
      -0.25*CoeffDispu(4)*delta(3)*Fields(2,1,xi1,yi1,z0)**2 &
      -0.25*CoeffDispu(2)*delta(1)*Fields(1,2,xi1,y0,zi1)*Fields(2,1,xi1,y0,zi1) &
      +0.25*CoeffDispu(2)*delta(2)*Fields(1,2,xi1,y0,zi1)*Fields(2,1,xi1,y0,zi1) &
      +0.25*CoeffDispu(6)*delta(1)*Fields(2,1,xi1,y0,zi1)*Fields(2,2,xi1,y0,zi1) &
      -0.25*CoeffDispu(8)*delta(2)*Fields(2,1,xi1,y0,zi1)*Fields(2,2,xi1,y0,zi1) &
      +0.25*CoeffDispu(6)*delta(3)*Fields(2,1,xi1,y0,zi1)*Fields(2,2,xi1,y0,zi1) &
      +0.25*CoeffDispu(2)*delta(2)*Fields(2,1,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
      -0.25*CoeffDispu(2)*delta(3)*Fields(2,1,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
      -0.25*CoeffDispu(1)*delta(1)*Fields(1,1,xi1,y0,zi1)*Fields(2,1,xi1,y0,zi1) &
      +0.25*CoeffDispu(1)*delta(2)*Fields(1,1,xi1,y0,zi1)*Fields(2,1,xi1,y0,zi1) &
      +0.25*CoeffDispu(4)*delta(1)*Fields(2,1,xi1,y0,zi1)**2 &
      -0.25*CoeffDispu(5)*delta(2)*Fields(2,1,xi1,y0,zi1)**2 &
      +0.25*CoeffDispu(4)*delta(3)*Fields(2,1,xi1,y0,zi1)**2 &
      -0.25*CoeffDispu(2)*delta(1)*Fields(1,2,xi1,y0,z0)*Fields(2,1,xi1,y0,z0) &
      +0.25*CoeffDispu(2)*delta(2)*Fields(1,2,xi1,y0,z0)*Fields(2,1,xi1,y0,z0) &
      +0.25*CoeffDispu(6)*delta(1)*Fields(2,1,xi1,y0,z0)*Fields(2,2,xi1,y0,z0) &
      -0.25*CoeffDispu(8)*delta(2)*Fields(2,1,xi1,y0,z0)*Fields(2,2,xi1,y0,z0) &
      -0.25*CoeffDispu(6)*delta(3)*Fields(2,1,xi1,y0,z0)*Fields(2,2,xi1,y0,z0) &
      -0.25*CoeffDispu(2)*delta(2)*Fields(2,1,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
      -0.25*CoeffDispu(2)*delta(3)*Fields(2,1,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
      -0.25*CoeffDispu(1)*delta(1)*Fields(1,1,xi1,y0,z0)*Fields(2,1,xi1,y0,z0) &
      +0.25*CoeffDispu(1)*delta(2)*Fields(1,1,xi1,y0,z0)*Fields(2,1,xi1,y0,z0) &
      +0.25*CoeffDispu(4)*delta(1)*Fields(2,1,xi1,y0,z0)**2 &
      -0.25*CoeffDispu(5)*delta(2)*Fields(2,1,xi1,y0,z0)**2 &
      -0.25*CoeffDispu(4)*delta(3)*Fields(2,1,xi1,y0,z0)**2 &
      +0.25*CoeffDispu(2)*delta(1)*Fields(1,2,x0,yi1,zi1)*Fields(2,1,x0,yi1,zi1) &
      -0.25*CoeffDispu(2)*delta(2)*Fields(1,2,x0,yi1,zi1)*Fields(2,1,x0,yi1,zi1) &
      -0.25*CoeffDispu(6)*delta(1)*Fields(2,1,x0,yi1,zi1)*Fields(2,2,x0,yi1,zi1) &
      +0.25*CoeffDispu(8)*delta(2)*Fields(2,1,x0,yi1,zi1)*Fields(2,2,x0,yi1,zi1) &
      +0.25*CoeffDispu(6)*delta(3)*Fields(2,1,x0,yi1,zi1)*Fields(2,2,x0,yi1,zi1) &
      +0.25*CoeffDispu(2)*delta(2)*Fields(2,1,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
      +0.25*CoeffDispu(2)*delta(3)*Fields(2,1,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
      +0.25*CoeffDispu(1)*delta(1)*Fields(1,1,x0,yi1,zi1)*Fields(2,1,x0,yi1,zi1) &
      -0.25*CoeffDispu(1)*delta(2)*Fields(1,1,x0,yi1,zi1)*Fields(2,1,x0,yi1,zi1) &
      -0.25*CoeffDispu(4)*delta(1)*Fields(2,1,x0,yi1,zi1)**2 &
      +0.25*CoeffDispu(5)*delta(2)*Fields(2,1,x0,yi1,zi1)**2 &
      +0.25*CoeffDispu(4)*delta(3)*Fields(2,1,x0,yi1,zi1)**2 &
      +0.25*CoeffDispu(2)*delta(1)*Fields(1,2,x0,yi1,z0)*Fields(2,1,x0,yi1,z0) &
      -0.25*CoeffDispu(2)*delta(2)*Fields(1,2,x0,yi1,z0)*Fields(2,1,x0,yi1,z0) &
      -0.25*CoeffDispu(6)*delta(1)*Fields(2,1,x0,yi1,z0)*Fields(2,2,x0,yi1,z0) &
      +0.25*CoeffDispu(8)*delta(2)*Fields(2,1,x0,yi1,z0)*Fields(2,2,x0,yi1,z0) &
      -0.25*CoeffDispu(6)*delta(3)*Fields(2,1,x0,yi1,z0)*Fields(2,2,x0,yi1,z0) &
      -0.25*CoeffDispu(2)*delta(2)*Fields(2,1,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
      +0.25*CoeffDispu(2)*delta(3)*Fields(2,1,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
      +0.25*CoeffDispu(1)*delta(1)*Fields(1,1,x0,yi1,z0)*Fields(2,1,x0,yi1,z0) &
      -0.25*CoeffDispu(1)*delta(2)*Fields(1,1,x0,yi1,z0)*Fields(2,1,x0,yi1,z0) &
      -0.25*CoeffDispu(4)*delta(1)*Fields(2,1,x0,yi1,z0)**2 &
      +0.25*CoeffDispu(5)*delta(2)*Fields(2,1,x0,yi1,z0)**2 &
      -0.25*CoeffDispu(4)*delta(3)*Fields(2,1,x0,yi1,z0)**2 &
      -0.25*CoeffDispu(2)*delta(1)*Fields(1,2,x0,y0,zi1)*Fields(2,1,x0,y0,zi1) &
      -0.25*CoeffDispu(2)*delta(2)*Fields(1,2,x0,y0,zi1)*Fields(2,1,x0,y0,zi1) &
      -0.25*CoeffDispu(6)*delta(1)*Fields(2,1,x0,y0,zi1)*Fields(2,2,x0,y0,zi1) &
      -0.25*CoeffDispu(8)*delta(2)*Fields(2,1,x0,y0,zi1)*Fields(2,2,x0,y0,zi1) &
      +0.25*CoeffDispu(6)*delta(3)*Fields(2,1,x0,y0,zi1)*Fields(2,2,x0,y0,zi1) &
      +0.25*CoeffDispu(2)*delta(2)*Fields(2,1,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
      -0.25*CoeffDispu(2)*delta(3)*Fields(2,1,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
      -0.25*CoeffDispu(1)*delta(1)*Fields(1,1,x0,y0,zi1)*Fields(2,1,x0,y0,zi1) &
      -0.25*CoeffDispu(1)*delta(2)*Fields(1,1,x0,y0,zi1)*Fields(2,1,x0,y0,zi1) &
      -0.25*CoeffDispu(4)*delta(1)*Fields(2,1,x0,y0,zi1)**2 &
      -0.25*CoeffDispu(5)*delta(2)*Fields(2,1,x0,y0,zi1)**2 &
      +0.25*CoeffDispu(4)*delta(3)*Fields(2,1,x0,y0,zi1)**2 &
      -0.25*CoeffDispu(2)*delta(1)*Fields(1,2,x0,y0,z0)*Fields(2,1,x0,y0,z0) &
      -0.25*CoeffDispu(2)*delta(2)*Fields(1,2,x0,y0,z0)*Fields(2,1,x0,y0,z0) &
      -0.25*CoeffDispu(6)*delta(1)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      -0.25*CoeffDispu(8)*delta(2)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      -0.25*CoeffDispu(6)*delta(3)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      -0.25*CoeffDispu(2)*delta(2)*Fields(2,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      -0.25*CoeffDispu(2)*delta(3)*Fields(2,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      -0.25*CoeffDispu(1)*delta(1)*Fields(1,1,x0,y0,z0)*Fields(2,1,x0,y0,z0) &
      -0.25*CoeffDispu(1)*delta(2)*Fields(1,1,x0,y0,z0)*Fields(2,1,x0,y0,z0) &
      -0.25*CoeffDispu(4)*delta(1)*Fields(2,1,x0,y0,z0)**2 &
      -0.25*CoeffDispu(5)*delta(2)*Fields(2,1,x0,y0,z0)**2 &
      -0.25*CoeffDispu(4)*delta(3)*Fields(2,1,x0,y0,z0)**2 &
      +0.25*CoeffDispu(2)*delta(1)*Fields(1,2,xi1,yi1,zi1)*Fields(3,1,xi1,yi1,zi1) &
      +0.25*CoeffDispu(2)*delta(3)*Fields(1,2,xi1,yi1,zi1)*Fields(3,1,xi1,yi1,zi1) &
      +0.25*CoeffDispu(2)*delta(2)*Fields(2,2,xi1,yi1,zi1)*Fields(3,1,xi1,yi1,zi1)
    
    DeltaHDispu = &
      DeltaHDispu &
      -0.25*CoeffDispu(1)*delta(1)*Fields(1,1,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      -0.25*CoeffDispu(1)*delta(3)*Fields(1,1,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      -0.25*CoeffDispu(2)*delta(1)*Fields(1,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      -0.25*CoeffDispu(2)*delta(3)*Fields(1,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      -0.25*CoeffDispu(1)*delta(2)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      -0.25*CoeffDispu(1)*delta(3)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      -0.25*CoeffDispu(2)*delta(2)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      -0.25*CoeffDispu(2)*delta(3)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      -0.25*CoeffDispu(4)*delta(1)*Fields(3,1,x0,y0,z0)**2 &
      -0.25*CoeffDispu(4)*delta(2)*Fields(3,1,x0,y0,z0)**2 &
      -0.25*CoeffDispu(5)*delta(3)*Fields(3,1,x0,y0,z0)**2 &
      +0.25*CoeffDispu(1)*delta(1)*Fields(1,1,x0,y0,zi1)*Fields(3,1,x0,y0,zi1) &
      -0.25*CoeffDispu(1)*delta(3)*Fields(1,1,x0,y0,zi1)*Fields(3,1,x0,y0,zi1) &
      +0.25*CoeffDispu(2)*delta(1)*Fields(1,2,x0,y0,zi1)*Fields(3,1,x0,y0,zi1) &
      -0.25*CoeffDispu(2)*delta(3)*Fields(1,2,x0,y0,zi1)*Fields(3,1,x0,y0,zi1) &
      +0.25*CoeffDispu(1)*delta(2)*Fields(2,1,x0,y0,zi1)*Fields(3,1,x0,y0,zi1) &
      -0.25*CoeffDispu(1)*delta(3)*Fields(2,1,x0,y0,zi1)*Fields(3,1,x0,y0,zi1) &
      +0.25*CoeffDispu(2)*delta(2)*Fields(2,2,x0,y0,zi1)*Fields(3,1,x0,y0,zi1) &
      -0.25*CoeffDispu(2)*delta(3)*Fields(2,2,x0,y0,zi1)*Fields(3,1,x0,y0,zi1) &
      -0.25*CoeffDispu(4)*delta(1)*Fields(3,1,x0,y0,zi1)**2 &
      -0.25*CoeffDispu(4)*delta(2)*Fields(3,1,x0,y0,zi1)**2 &
      +0.25*CoeffDispu(5)*delta(3)*Fields(3,1,x0,y0,zi1)**2 &
      -0.25*CoeffDispu(1)*delta(1)*Fields(1,1,x0,yi1,z0)*Fields(3,1,x0,yi1,z0) &
      -0.25*CoeffDispu(1)*delta(3)*Fields(1,1,x0,yi1,z0)*Fields(3,1,x0,yi1,z0) &
      -0.25*CoeffDispu(2)*delta(1)*Fields(1,2,x0,yi1,z0)*Fields(3,1,x0,yi1,z0) &
      -0.25*CoeffDispu(2)*delta(3)*Fields(1,2,x0,yi1,z0)*Fields(3,1,x0,yi1,z0) &
      -0.25*CoeffDispu(1)*delta(2)*Fields(2,1,x0,yi1,z0)*Fields(3,1,x0,yi1,z0) &
      +0.25*CoeffDispu(1)*delta(3)*Fields(2,1,x0,yi1,z0)*Fields(3,1,x0,yi1,z0) &
      -0.25*CoeffDispu(2)*delta(2)*Fields(2,2,x0,yi1,z0)*Fields(3,1,x0,yi1,z0) &
      +0.25*CoeffDispu(2)*delta(3)*Fields(2,2,x0,yi1,z0)*Fields(3,1,x0,yi1,z0) &
      -0.25*CoeffDispu(4)*delta(1)*Fields(3,1,x0,yi1,z0)**2 &
      +0.25*CoeffDispu(4)*delta(2)*Fields(3,1,x0,yi1,z0)**2 &
      -0.25*CoeffDispu(5)*delta(3)*Fields(3,1,x0,yi1,z0)**2 &
      +0.25*CoeffDispu(1)*delta(1)*Fields(1,1,x0,yi1,zi1)*Fields(3,1,x0,yi1,zi1) &
      -0.25*CoeffDispu(1)*delta(3)*Fields(1,1,x0,yi1,zi1)*Fields(3,1,x0,yi1,zi1) &
      +0.25*CoeffDispu(2)*delta(1)*Fields(1,2,x0,yi1,zi1)*Fields(3,1,x0,yi1,zi1) &
      -0.25*CoeffDispu(2)*delta(3)*Fields(1,2,x0,yi1,zi1)*Fields(3,1,x0,yi1,zi1) &
      +0.25*CoeffDispu(1)*delta(2)*Fields(2,1,x0,yi1,zi1)*Fields(3,1,x0,yi1,zi1) &
      +0.25*CoeffDispu(1)*delta(3)*Fields(2,1,x0,yi1,zi1)*Fields(3,1,x0,yi1,zi1) &
      +0.25*CoeffDispu(2)*delta(2)*Fields(2,2,x0,yi1,zi1)*Fields(3,1,x0,yi1,zi1) &
      +0.25*CoeffDispu(2)*delta(3)*Fields(2,2,x0,yi1,zi1)*Fields(3,1,x0,yi1,zi1) &
      -0.25*CoeffDispu(4)*delta(1)*Fields(3,1,x0,yi1,zi1)**2 &
      +0.25*CoeffDispu(4)*delta(2)*Fields(3,1,x0,yi1,zi1)**2 &
      +0.25*CoeffDispu(5)*delta(3)*Fields(3,1,x0,yi1,zi1)**2 &
      -0.25*CoeffDispu(1)*delta(1)*Fields(1,1,xi1,y0,z0)*Fields(3,1,xi1,y0,z0) &
      +0.25*CoeffDispu(1)*delta(3)*Fields(1,1,xi1,y0,z0)*Fields(3,1,xi1,y0,z0) &
      -0.25*CoeffDispu(2)*delta(1)*Fields(1,2,xi1,y0,z0)*Fields(3,1,xi1,y0,z0) &
      +0.25*CoeffDispu(2)*delta(3)*Fields(1,2,xi1,y0,z0)*Fields(3,1,xi1,y0,z0) &
      -0.25*CoeffDispu(1)*delta(2)*Fields(2,1,xi1,y0,z0)*Fields(3,1,xi1,y0,z0) &
      -0.25*CoeffDispu(1)*delta(3)*Fields(2,1,xi1,y0,z0)*Fields(3,1,xi1,y0,z0) &
      -0.25*CoeffDispu(2)*delta(2)*Fields(2,2,xi1,y0,z0)*Fields(3,1,xi1,y0,z0) &
      -0.25*CoeffDispu(2)*delta(3)*Fields(2,2,xi1,y0,z0)*Fields(3,1,xi1,y0,z0) &
      +0.25*CoeffDispu(4)*delta(1)*Fields(3,1,xi1,y0,z0)**2 &
      -0.25*CoeffDispu(4)*delta(2)*Fields(3,1,xi1,y0,z0)**2 &
      -0.25*CoeffDispu(5)*delta(3)*Fields(3,1,xi1,y0,z0)**2 &
      +0.25*CoeffDispu(1)*delta(1)*Fields(1,1,xi1,y0,zi1)*Fields(3,1,xi1,y0,zi1) &
      +0.25*CoeffDispu(1)*delta(3)*Fields(1,1,xi1,y0,zi1)*Fields(3,1,xi1,y0,zi1) &
      +0.25*CoeffDispu(2)*delta(1)*Fields(1,2,xi1,y0,zi1)*Fields(3,1,xi1,y0,zi1) &
      +0.25*CoeffDispu(2)*delta(3)*Fields(1,2,xi1,y0,zi1)*Fields(3,1,xi1,y0,zi1) &
      +0.25*CoeffDispu(1)*delta(2)*Fields(2,1,xi1,y0,zi1)*Fields(3,1,xi1,y0,zi1) &
      -0.25*CoeffDispu(1)*delta(3)*Fields(2,1,xi1,y0,zi1)*Fields(3,1,xi1,y0,zi1) &
      +0.25*CoeffDispu(2)*delta(2)*Fields(2,2,xi1,y0,zi1)*Fields(3,1,xi1,y0,zi1) &
      -0.25*CoeffDispu(2)*delta(3)*Fields(2,2,xi1,y0,zi1)*Fields(3,1,xi1,y0,zi1) &
      +0.25*CoeffDispu(4)*delta(1)*Fields(3,1,xi1,y0,zi1)**2 &
      -0.25*CoeffDispu(4)*delta(2)*Fields(3,1,xi1,y0,zi1)**2 &
      +0.25*CoeffDispu(5)*delta(3)*Fields(3,1,xi1,y0,zi1)**2 &
      -0.25*CoeffDispu(1)*delta(1)*Fields(1,1,xi1,yi1,z0)*Fields(3,1,xi1,yi1,z0) &
      +0.25*CoeffDispu(1)*delta(3)*Fields(1,1,xi1,yi1,z0)*Fields(3,1,xi1,yi1,z0) &
      -0.25*CoeffDispu(2)*delta(1)*Fields(1,2,xi1,yi1,z0)*Fields(3,1,xi1,yi1,z0) &
      +0.25*CoeffDispu(2)*delta(3)*Fields(1,2,xi1,yi1,z0)*Fields(3,1,xi1,yi1,z0) &
      -0.25*CoeffDispu(1)*delta(2)*Fields(2,1,xi1,yi1,z0)*Fields(3,1,xi1,yi1,z0) &
      +0.25*CoeffDispu(1)*delta(3)*Fields(2,1,xi1,yi1,z0)*Fields(3,1,xi1,yi1,z0) &
      -0.25*CoeffDispu(2)*delta(2)*Fields(2,2,xi1,yi1,z0)*Fields(3,1,xi1,yi1,z0) &
      +0.25*CoeffDispu(2)*delta(3)*Fields(2,2,xi1,yi1,z0)*Fields(3,1,xi1,yi1,z0) &
      +0.25*CoeffDispu(4)*delta(1)*Fields(3,1,xi1,yi1,z0)**2 &
      +0.25*CoeffDispu(4)*delta(2)*Fields(3,1,xi1,yi1,z0)**2 &
      -0.25*CoeffDispu(5)*delta(3)*Fields(3,1,xi1,yi1,z0)**2 &
      +0.25*CoeffDispu(1)*delta(1)*Fields(1,1,xi1,yi1,zi1)*Fields(3,1,xi1,yi1,zi1) &
      +0.25*CoeffDispu(1)*delta(3)*Fields(1,1,xi1,yi1,zi1)*Fields(3,1,xi1,yi1,zi1) &
      +0.25*CoeffDispu(1)*delta(2)*Fields(2,1,xi1,yi1,zi1)*Fields(3,1,xi1,yi1,zi1) &
      +0.25*CoeffDispu(1)*delta(3)*Fields(2,1,xi1,yi1,zi1)*Fields(3,1,xi1,yi1,zi1) &
      +0.25*CoeffDispu(2)*delta(3)*Fields(2,2,xi1,yi1,zi1)*Fields(3,1,xi1,yi1,zi1) &
      +0.25*CoeffDispu(4)*delta(1)*Fields(3,1,xi1,yi1,zi1)**2 &
      +0.25*CoeffDispu(4)*delta(2)*Fields(3,1,xi1,yi1,zi1)**2 &
      +0.25*CoeffDispu(5)*delta(3)*Fields(3,1,xi1,yi1,zi1)**2 &
      -0.25*CoeffDispu(6)*delta(1)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      -0.25*CoeffDispu(6)*delta(2)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      -0.25*CoeffDispu(8)*delta(3)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      -0.25*CoeffDispu(6)*delta(1)*Fields(3,1,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
      -0.25*CoeffDispu(6)*delta(2)*Fields(3,1,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
      +0.25*CoeffDispu(8)*delta(3)*Fields(3,1,x0,y0,zi1)*Fields(3,2,x0,y0,zi1) &
      -0.25*CoeffDispu(6)*delta(1)*Fields(3,1,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
      +0.25*CoeffDispu(6)*delta(2)*Fields(3,1,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
      -0.25*CoeffDispu(8)*delta(3)*Fields(3,1,x0,yi1,z0)*Fields(3,2,x0,yi1,z0) &
      -0.25*CoeffDispu(6)*delta(1)*Fields(3,1,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
      +0.25*CoeffDispu(6)*delta(2)*Fields(3,1,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
      +0.25*CoeffDispu(8)*delta(3)*Fields(3,1,x0,yi1,zi1)*Fields(3,2,x0,yi1,zi1) &
      +0.25*CoeffDispu(6)*delta(1)*Fields(3,1,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
      -0.25*CoeffDispu(6)*delta(2)*Fields(3,1,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
      -0.25*CoeffDispu(8)*delta(3)*Fields(3,1,xi1,y0,z0)*Fields(3,2,xi1,y0,z0) &
      +0.25*CoeffDispu(6)*delta(1)*Fields(3,1,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
      -0.25*CoeffDispu(6)*delta(2)*Fields(3,1,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
      +0.25*CoeffDispu(8)*delta(3)*Fields(3,1,xi1,y0,zi1)*Fields(3,2,xi1,y0,zi1) &
      +0.25*CoeffDispu(6)*delta(1)*Fields(3,1,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
      +0.25*CoeffDispu(6)*delta(2)*Fields(3,1,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
      -0.25*CoeffDispu(8)*delta(3)*Fields(3,1,xi1,yi1,z0)*Fields(3,2,xi1,yi1,z0) &
      +0.25*CoeffDispu(6)*delta(1)*Fields(3,1,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
      +0.25*CoeffDispu(6)*delta(2)*Fields(3,1,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1) &
      +0.25*CoeffDispu(8)*delta(3)*Fields(3,1,xi1,yi1,zi1)*Fields(3,2,xi1,yi1,zi1)
    
    DeltaHEpsDisp = &
      DeltaHEpsDisp &
      +0.
    
    DeltaHEpsu = &
      DeltaHEpsu &
      +0. 
!      -0.25*CoeffEpsu(4)*delta(1)*e0ij(1,1) &
!      -0.25*CoeffEpsu(3)*delta(2)*e0ij(1,1) &
!      -0.25*CoeffEpsu(3)*delta(3)*e0ij(1,1) &
!      -0.25*CoeffEpsu(2)*delta(1)*e0ij(1,2) &
!      -0.25*CoeffEpsu(2)*delta(2)*e0ij(1,2) &
!      -0.25*CoeffEpsu(2)*delta(1)*e0ij(1,3) &
!      -0.25*CoeffEpsu(2)*delta(3)*e0ij(1,3) &
!      -0.25*CoeffEpsu(3)*delta(1)*e0ij(2,2) &
!      -0.25*CoeffEpsu(4)*delta(2)*e0ij(2,2) &
!      -0.25*CoeffEpsu(3)*delta(3)*e0ij(2,2) &
!      -0.25*CoeffEpsu(2)*delta(2)*e0ij(2,3) &
!      -0.25*CoeffEpsu(2)*delta(3)*e0ij(2,3) &
!      -0.25*CoeffEpsu(3)*delta(1)*e0ij(3,3) &
!      -0.25*CoeffEpsu(3)*delta(2)*e0ij(3,3) &
!      -0.25*CoeffEpsu(4)*delta(3)*e0ij(3,3)
    
    DeltaHEps = &
      DeltaHEps &
      +0.
    
    DeltaH=DeltaHDisp+DeltaHJijsm+DeltaHJijhf+DeltaHu+DeltaHDispu+DeltaHEpsDisp+DeltaHEpsu+DeltaHEps
  
  CASE (4)
    DeltaHDisp = &
      DeltaHDisp &
      +0.
    
    DeltaHJijsm = &
      DeltaHJijsm &
      +0.
    
    DeltaHJijhf = &
      DeltaHJijhf &
      +0.
    
    DeltaHu = &
      DeltaHu &
      +0.
    
    DeltaHDispu = &
      DeltaHDispu &
      +0.
    
    DeltaHEpsDisp = &
      DeltaHEpsDisp &
      +CoeffEpsDisp(5)*delta(1)*Fields(1,1,x0,y0,z0)**2 &
      +CoeffEpsDisp(4)*delta(2)*Fields(1,1,x0,y0,z0)**2 &
      +CoeffEpsDisp(4)*delta(3)*Fields(1,1,x0,y0,z0)**2 &
      +CoeffEpsDisp(8)*delta(1)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(2)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(3)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
      +CoeffEpsDisp(9)*delta(1)*Fields(1,2,x0,y0,z0)**2 &
      +CoeffEpsDisp(7)*delta(2)*Fields(1,2,x0,y0,z0)**2 &
      +CoeffEpsDisp(7)*delta(3)*Fields(1,2,x0,y0,z0)**2 &
      +CoeffEpsDisp(1)*delta(6)*Fields(1,1,x0,y0,z0)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(6)*Fields(1,2,x0,y0,z0)*Fields(2,1,x0,y0,z0) &
      +CoeffEpsDisp(4)*delta(1)*Fields(2,1,x0,y0,z0)**2 &
      +CoeffEpsDisp(5)*delta(2)*Fields(2,1,x0,y0,z0)**2 &
      +CoeffEpsDisp(4)*delta(3)*Fields(2,1,x0,y0,z0)**2 &
      +CoeffEpsDisp(2)*delta(6)*Fields(1,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(6)*Fields(1,2,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(1)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(8)*delta(2)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(3)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
      +CoeffEpsDisp(7)*delta(1)*Fields(2,2,x0,y0,z0)**2 &
      +CoeffEpsDisp(9)*delta(2)*Fields(2,2,x0,y0,z0)**2 &
      +CoeffEpsDisp(7)*delta(3)*Fields(2,2,x0,y0,z0)**2 &
      +CoeffEpsDisp(1)*delta(5)*Fields(1,1,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(5)*Fields(1,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(1)*delta(4)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(4)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
      +CoeffEpsDisp(4)*delta(1)*Fields(3,1,x0,y0,z0)**2 &
      +CoeffEpsDisp(4)*delta(2)*Fields(3,1,x0,y0,z0)**2 &
      +CoeffEpsDisp(5)*delta(3)*Fields(3,1,x0,y0,z0)**2 &
      +CoeffEpsDisp(2)*delta(5)*Fields(1,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(5)*Fields(1,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(2)*delta(4)*Fields(2,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(3)*delta(4)*Fields(2,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(1)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(6)*delta(2)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(8)*delta(3)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
      +CoeffEpsDisp(7)*delta(1)*Fields(3,2,x0,y0,z0)**2 &
      +CoeffEpsDisp(7)*delta(2)*Fields(3,2,x0,y0,z0)**2 &
      +CoeffEpsDisp(9)*delta(3)*Fields(3,2,x0,y0,z0)**2
    
    DeltaHEpsu = &
      DeltaHEpsu &
      +0.
    
    DeltaHEps = &
      DeltaHEps &
      +CoeffEps(1)*delta(1) &
      +0.5*CoeffEps(4)*delta(1)**2 &
      +CoeffEps(1)*delta(2) &
      +CoeffEps(3)*delta(1)*delta(2) &
      +0.5*CoeffEps(4)*delta(2)**2 &
      +CoeffEps(1)*delta(3) &
      +CoeffEps(3)*delta(1)*delta(3) &
      +CoeffEps(3)*delta(2)*delta(3) &
      +0.5*CoeffEps(4)*delta(3)**2 &
      +0.5*CoeffEps(2)*delta(4)**2 &
      +0.5*CoeffEps(2)*delta(5)**2 &
      +0.5*CoeffEps(2)*delta(6)**2 &
      +CoeffEps(4)*delta(1)*e0ij(1,1) &
      +CoeffEps(3)*delta(2)*e0ij(1,1) &
      +CoeffEps(3)*delta(3)*e0ij(1,1) &
      +CoeffEps(2)*delta(6)*e0ij(1,2) &
      +CoeffEps(2)*delta(5)*e0ij(1,3) &
      +CoeffEps(3)*delta(1)*e0ij(2,2) &
      +CoeffEps(4)*delta(2)*e0ij(2,2) &
      +CoeffEps(3)*delta(3)*e0ij(2,2) &
      +CoeffEps(2)*delta(4)*e0ij(2,3) &
      +CoeffEps(3)*delta(1)*e0ij(3,3) &
      +CoeffEps(3)*delta(2)*e0ij(3,3) &
      +CoeffEps(4)*delta(3)*e0ij(3,3) 
    
    DeltaH=DeltaHDisp+DeltaHJijsm+DeltaHJijhf+DeltaHu+DeltaHDispu+DeltaHEpsDisp+DeltaHEpsu+DeltaHEps
  
  CASE DEFAULT
      write(*,*) "mode out of range!"
      call abort
  END SELECT
  
End Function GetDeltaH
