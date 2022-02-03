Function GetDeltaHDisp(x0, y0, z0, Fields, e0ij, idelta, delta) Result(DeltaH)
  
  Implicit none
  Integer, Intent(in) :: x0, y0, z0, idelta
  Real*8,  Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(in) :: delta(FieldDimList(idelta))
  Real*8,  Intent(in) :: e0ij(3,3)
  Real*8              :: eij(3,3), euij(3,3)
  Real*8              :: DeltaH
  Real*8              :: DeltaHDisp
  
  Integer             :: x1,y1,z1,xi1,yi1,zi1
  
  x1 = (x0+1)-floor(real(x0+1-1)/real(cgrid%n1))*cgrid%n1
  y1 = (y0+1)-floor(real(y0+1-1)/real(cgrid%n2))*cgrid%n2
  z1 = (z0+1)-floor(real(z0+1-1)/real(cgrid%n3))*cgrid%n3
  xi1 = (x0-1)-floor(real(x0-1-1)/real(cgrid%n1))*cgrid%n1
  yi1 = (y0-1)-floor(real(y0-1-1)/real(cgrid%n2))*cgrid%n2
  zi1 = (z0-1)-floor(real(z0-1-1)/real(cgrid%n3))*cgrid%n3
  
  DeltaH = 0.0D0
  DeltaHDisp = 0.0D0
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
    
    DeltaH=DeltaHDisp
  
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
    
    DeltaH=DeltaHDisp
  
  CASE (3)
    DeltaHDisp = &
      DeltaHDisp &
      +0.
    
    DeltaH=DeltaHDisp
  
  CASE (4)
    DeltaHDisp = &
      DeltaHDisp &
      +0.
    
    DeltaH=DeltaHDisp
  
  CASE DEFAULT
      write(*,*) "mode out of range!"
      call abort
  END SELECT
  
End Function GetDeltaHDisp
