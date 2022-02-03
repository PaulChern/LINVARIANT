Function GetEOnSite(x0, y0, z0, Fields, e0ij) Result(EOnSite)
  
  Implicit none
  Real*8,  Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(in) :: e0ij(3,3)
  Integer, Intent(in) :: x0, y0, z0
  Real*8              :: eij(3,3), euij(3,3)
  Real*8              :: EOnSite
  Real*8              :: EOnSiteDisp
  Real*8              :: EOnSiteJijsm
  Real*8              :: EOnSiteJijhf
  Real*8              :: EOnSiteEpsDisp
  Real*8              :: EOnSiteEps
  
  Integer             :: x1,y1,z1,xi1,yi1,zi1
  
  EOnSite = 0.0D0
  EOnSiteDisp = 0.0D0
  EOnSiteJijsm = 0.0D0
  EOnSiteJijhf = 0.0D0
  EOnSiteEpsDisp = 0.0D0
  EOnSiteEps = 0.0D0
  
  
  euij = StrainFromu(x0, y0, z0, Fields)
  eij = e0ij + euij
  
  x1 = (x0+1)-floor(real(x0+1-1)/real(cgrid%n1))*cgrid%n1
  y1 = (y0+1)-floor(real(y0+1-1)/real(cgrid%n2))*cgrid%n2
  z1 = (z0+1)-floor(real(z0+1-1)/real(cgrid%n3))*cgrid%n3
  xi1 = (x0-1)-floor(real(x0-1-1)/real(cgrid%n1))*cgrid%n1
  yi1 = (y0-1)-floor(real(y0-1-1)/real(cgrid%n2))*cgrid%n2
  zi1 = (z0-1)-floor(real(z0-1-1)/real(cgrid%n3))*cgrid%n3
  
  EOnSiteDisp = &
    EOnSiteDisp &
    +CoeffDisp(1)*Fields(1,1,x0,y0,z0)**2 &
    +CoeffDisp(5)*Fields(1,1,x0,y0,z0)**4 &
    +CoeffDisp(2)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
    +CoeffDisp(7)*Fields(1,1,x0,y0,z0)**3*Fields(1,2,x0,y0,z0) &
    +CoeffDisp(3)*Fields(1,2,x0,y0,z0)**2 &
    +CoeffDisp(10)*Fields(1,1,x0,y0,z0)**2*Fields(1,2,x0,y0,z0)**2 &
    +CoeffDisp(13)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)**3 &
    +CoeffDisp(14)*Fields(1,2,x0,y0,z0)**4 &
    +CoeffDisp(1)*Fields(2,1,x0,y0,z0)**2 &
    +CoeffDisp(4)*Fields(1,1,x0,y0,z0)**2*Fields(2,1,x0,y0,z0)**2 &
    +CoeffDisp(6)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)*Fields(2,1,x0,y0,z0)**2 &
    +CoeffDisp(9)*Fields(1,2,x0,y0,z0)**2*Fields(2,1,x0,y0,z0)**2 &
    +CoeffDisp(5)*Fields(2,1,x0,y0,z0)**4 &
    +CoeffDisp(2)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +CoeffDisp(6)*Fields(1,1,x0,y0,z0)**2*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +CoeffDisp(8)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +CoeffDisp(11)*Fields(1,2,x0,y0,z0)**2*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +CoeffDisp(7)*Fields(2,1,x0,y0,z0)**3*Fields(2,2,x0,y0,z0) &
    +CoeffDisp(3)*Fields(2,2,x0,y0,z0)**2 &
    +CoeffDisp(9)*Fields(1,1,x0,y0,z0)**2*Fields(2,2,x0,y0,z0)**2 &
    +CoeffDisp(11)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)*Fields(2,2,x0,y0,z0)**2 &
    +CoeffDisp(12)*Fields(1,2,x0,y0,z0)**2*Fields(2,2,x0,y0,z0)**2 &
    +CoeffDisp(10)*Fields(2,1,x0,y0,z0)**2*Fields(2,2,x0,y0,z0)**2 &
    +CoeffDisp(13)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0)**3 &
    +CoeffDisp(14)*Fields(2,2,x0,y0,z0)**4 &
    +CoeffDisp(1)*Fields(3,1,x0,y0,z0)**2 &
    +CoeffDisp(4)*Fields(1,1,x0,y0,z0)**2*Fields(3,1,x0,y0,z0)**2 &
    +CoeffDisp(6)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)*Fields(3,1,x0,y0,z0)**2 &
    +CoeffDisp(9)*Fields(1,2,x0,y0,z0)**2*Fields(3,1,x0,y0,z0)**2 &
    +CoeffDisp(4)*Fields(2,1,x0,y0,z0)**2*Fields(3,1,x0,y0,z0)**2 &
    +CoeffDisp(6)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y0,z0)**2 &
    +CoeffDisp(9)*Fields(2,2,x0,y0,z0)**2*Fields(3,1,x0,y0,z0)**2 &
    +CoeffDisp(5)*Fields(3,1,x0,y0,z0)**4 &
    +CoeffDisp(2)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +CoeffDisp(6)*Fields(1,1,x0,y0,z0)**2*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +CoeffDisp(8)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +CoeffDisp(11)*Fields(1,2,x0,y0,z0)**2*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +CoeffDisp(6)*Fields(2,1,x0,y0,z0)**2*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +CoeffDisp(8)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +CoeffDisp(11)*Fields(2,2,x0,y0,z0)**2*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +CoeffDisp(7)*Fields(3,1,x0,y0,z0)**3*Fields(3,2,x0,y0,z0) &
    +CoeffDisp(3)*Fields(3,2,x0,y0,z0)**2 &
    +CoeffDisp(9)*Fields(1,1,x0,y0,z0)**2*Fields(3,2,x0,y0,z0)**2 &
    +CoeffDisp(11)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0)*Fields(3,2,x0,y0,z0)**2 &
    +CoeffDisp(12)*Fields(1,2,x0,y0,z0)**2*Fields(3,2,x0,y0,z0)**2 &
    +CoeffDisp(9)*Fields(2,1,x0,y0,z0)**2*Fields(3,2,x0,y0,z0)**2 &
    +CoeffDisp(11)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0)*Fields(3,2,x0,y0,z0)**2 &
    +CoeffDisp(12)*Fields(2,2,x0,y0,z0)**2*Fields(3,2,x0,y0,z0)**2 &
    +CoeffDisp(10)*Fields(3,1,x0,y0,z0)**2*Fields(3,2,x0,y0,z0)**2 &
    +CoeffDisp(13)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0)**3 &
    +CoeffDisp(14)*Fields(3,2,x0,y0,z0)**4
  
  EOnSiteJijsm = &
    EOnSiteJijsm &
    +2.*CoeffJijsm(2)*Fields(1,1,x0,y0,z0)*Fields(1,1,x0,y0,z1) &
    +2.*CoeffJijsm(2)*Fields(1,1,x0,y0,z0)*Fields(1,1,x0,y0,zi1) &
    +2.*CoeffJijsm(2)*Fields(1,1,x0,y0,z0)*Fields(1,1,x0,y1,z0) &
    +2.*CoeffJijsm(4)*Fields(1,1,x0,y0,z0)*Fields(1,1,x0,y1,z1) &
    +2.*CoeffJijsm(4)*Fields(1,1,x0,y0,z0)*Fields(1,1,x0,y1,zi1) &
    +2.*CoeffJijsm(2)*Fields(1,1,x0,y0,z0)*Fields(1,1,x0,yi1,z0) &
    +2.*CoeffJijsm(4)*Fields(1,1,x0,y0,z0)*Fields(1,1,x0,yi1,z1) &
    +2.*CoeffJijsm(4)*Fields(1,1,x0,y0,z0)*Fields(1,1,x0,yi1,zi1) &
    +2.*CoeffJijsm(1)*Fields(1,1,x0,y0,z0)*Fields(1,1,x1,y0,z0) &
    +2.*CoeffJijsm(3)*Fields(1,1,x0,y0,z0)*Fields(1,1,x1,y0,z1) &
    +2.*CoeffJijsm(3)*Fields(1,1,x0,y0,z0)*Fields(1,1,x1,y0,zi1) &
    +2.*CoeffJijsm(3)*Fields(1,1,x0,y0,z0)*Fields(1,1,x1,y1,z0) &
    +2.*CoeffJijsm(5)*Fields(1,1,x0,y0,z0)*Fields(1,1,x1,y1,z1) &
    +2.*CoeffJijsm(5)*Fields(1,1,x0,y0,z0)*Fields(1,1,x1,y1,zi1) &
    +2.*CoeffJijsm(3)*Fields(1,1,x0,y0,z0)*Fields(1,1,x1,yi1,z0) &
    +2.*CoeffJijsm(5)*Fields(1,1,x0,y0,z0)*Fields(1,1,x1,yi1,z1) &
    +2.*CoeffJijsm(5)*Fields(1,1,x0,y0,z0)*Fields(1,1,x1,yi1,zi1) &
    +2.*CoeffJijsm(1)*Fields(1,1,x0,y0,z0)*Fields(1,1,xi1,y0,z0) &
    +2.*CoeffJijsm(3)*Fields(1,1,x0,y0,z0)*Fields(1,1,xi1,y0,z1) &
    +2.*CoeffJijsm(3)*Fields(1,1,x0,y0,z0)*Fields(1,1,xi1,y0,zi1) &
    +2.*CoeffJijsm(3)*Fields(1,1,x0,y0,z0)*Fields(1,1,xi1,y1,z0) &
    +2.*CoeffJijsm(5)*Fields(1,1,x0,y0,z0)*Fields(1,1,xi1,y1,z1) &
    +2.*CoeffJijsm(5)*Fields(1,1,x0,y0,z0)*Fields(1,1,xi1,y1,zi1) &
    +2.*CoeffJijsm(3)*Fields(1,1,x0,y0,z0)*Fields(1,1,xi1,yi1,z0) &
    +2.*CoeffJijsm(5)*Fields(1,1,x0,y0,z0)*Fields(1,1,xi1,yi1,z1) &
    +2.*CoeffJijsm(5)*Fields(1,1,x0,y0,z0)*Fields(1,1,xi1,yi1,zi1) &
    +2.*CoeffJijsm(6)*Fields(1,1,x1,y1,z0)*Fields(2,1,x0,y0,z0) &
    +2.*CoeffJijsm(7)*Fields(1,1,x1,y1,z1)*Fields(2,1,x0,y0,z0) &
    +2.*CoeffJijsm(7)*Fields(1,1,x1,y1,zi1)*Fields(2,1,x0,y0,z0) &
    -2.*CoeffJijsm(6)*Fields(1,1,x1,yi1,z0)*Fields(2,1,x0,y0,z0) &
    -2.*CoeffJijsm(7)*Fields(1,1,x1,yi1,z1)*Fields(2,1,x0,y0,z0) &
    -2.*CoeffJijsm(7)*Fields(1,1,x1,yi1,zi1)*Fields(2,1,x0,y0,z0) &
    -2.*CoeffJijsm(6)*Fields(1,1,xi1,y1,z0)*Fields(2,1,x0,y0,z0) &
    -2.*CoeffJijsm(7)*Fields(1,1,xi1,y1,z1)*Fields(2,1,x0,y0,z0) &
    -2.*CoeffJijsm(7)*Fields(1,1,xi1,y1,zi1)*Fields(2,1,x0,y0,z0) &
    +2.*CoeffJijsm(6)*Fields(1,1,xi1,yi1,z0)*Fields(2,1,x0,y0,z0) &
    +2.*CoeffJijsm(7)*Fields(1,1,xi1,yi1,z1)*Fields(2,1,x0,y0,z0) &
    +2.*CoeffJijsm(7)*Fields(1,1,xi1,yi1,zi1)*Fields(2,1,x0,y0,z0) &
    +2.*CoeffJijsm(2)*Fields(2,1,x0,y0,z0)*Fields(2,1,x0,y0,z1) &
    +2.*CoeffJijsm(2)*Fields(2,1,x0,y0,z0)*Fields(2,1,x0,y0,zi1) &
    +2.*CoeffJijsm(1)*Fields(2,1,x0,y0,z0)*Fields(2,1,x0,y1,z0) &
    +2.*CoeffJijsm(3)*Fields(2,1,x0,y0,z0)*Fields(2,1,x0,y1,z1) &
    +2.*CoeffJijsm(3)*Fields(2,1,x0,y0,z0)*Fields(2,1,x0,y1,zi1) &
    +2.*CoeffJijsm(1)*Fields(2,1,x0,y0,z0)*Fields(2,1,x0,yi1,z0) &
    +2.*CoeffJijsm(3)*Fields(2,1,x0,y0,z0)*Fields(2,1,x0,yi1,z1) &
    +2.*CoeffJijsm(3)*Fields(2,1,x0,y0,z0)*Fields(2,1,x0,yi1,zi1) &
    +2.*CoeffJijsm(2)*Fields(2,1,x0,y0,z0)*Fields(2,1,x1,y0,z0) &
    +2.*CoeffJijsm(4)*Fields(2,1,x0,y0,z0)*Fields(2,1,x1,y0,z1) &
    +2.*CoeffJijsm(4)*Fields(2,1,x0,y0,z0)*Fields(2,1,x1,y0,zi1) &
    +2.*CoeffJijsm(6)*Fields(1,1,x0,y0,z0)*Fields(2,1,x1,y1,z0) &
    +2.*CoeffJijsm(3)*Fields(2,1,x0,y0,z0)*Fields(2,1,x1,y1,z0) &
    +2.*CoeffJijsm(7)*Fields(1,1,x0,y0,z0)*Fields(2,1,x1,y1,z1) &
    +2.*CoeffJijsm(5)*Fields(2,1,x0,y0,z0)*Fields(2,1,x1,y1,z1) &
    +2.*CoeffJijsm(7)*Fields(1,1,x0,y0,z0)*Fields(2,1,x1,y1,zi1) &
    +2.*CoeffJijsm(5)*Fields(2,1,x0,y0,z0)*Fields(2,1,x1,y1,zi1) &
    -2.*CoeffJijsm(6)*Fields(1,1,x0,y0,z0)*Fields(2,1,x1,yi1,z0) &
    +2.*CoeffJijsm(3)*Fields(2,1,x0,y0,z0)*Fields(2,1,x1,yi1,z0) &
    -2.*CoeffJijsm(7)*Fields(1,1,x0,y0,z0)*Fields(2,1,x1,yi1,z1) &
    +2.*CoeffJijsm(5)*Fields(2,1,x0,y0,z0)*Fields(2,1,x1,yi1,z1) &
    -2.*CoeffJijsm(7)*Fields(1,1,x0,y0,z0)*Fields(2,1,x1,yi1,zi1) &
    +2.*CoeffJijsm(5)*Fields(2,1,x0,y0,z0)*Fields(2,1,x1,yi1,zi1) &
    +2.*CoeffJijsm(2)*Fields(2,1,x0,y0,z0)*Fields(2,1,xi1,y0,z0) &
    +2.*CoeffJijsm(4)*Fields(2,1,x0,y0,z0)*Fields(2,1,xi1,y0,z1) &
    +2.*CoeffJijsm(4)*Fields(2,1,x0,y0,z0)*Fields(2,1,xi1,y0,zi1) &
    -2.*CoeffJijsm(6)*Fields(1,1,x0,y0,z0)*Fields(2,1,xi1,y1,z0) &
    +2.*CoeffJijsm(3)*Fields(2,1,x0,y0,z0)*Fields(2,1,xi1,y1,z0) &
    -2.*CoeffJijsm(7)*Fields(1,1,x0,y0,z0)*Fields(2,1,xi1,y1,z1) &
    +2.*CoeffJijsm(5)*Fields(2,1,x0,y0,z0)*Fields(2,1,xi1,y1,z1) &
    -2.*CoeffJijsm(7)*Fields(1,1,x0,y0,z0)*Fields(2,1,xi1,y1,zi1) &
    +2.*CoeffJijsm(5)*Fields(2,1,x0,y0,z0)*Fields(2,1,xi1,y1,zi1) &
    +2.*CoeffJijsm(6)*Fields(1,1,x0,y0,z0)*Fields(2,1,xi1,yi1,z0) &
    +2.*CoeffJijsm(3)*Fields(2,1,x0,y0,z0)*Fields(2,1,xi1,yi1,z0) &
    +2.*CoeffJijsm(7)*Fields(1,1,x0,y0,z0)*Fields(2,1,xi1,yi1,z1) &
    +2.*CoeffJijsm(5)*Fields(2,1,x0,y0,z0)*Fields(2,1,xi1,yi1,z1) &
    +2.*CoeffJijsm(7)*Fields(1,1,x0,y0,z0)*Fields(2,1,xi1,yi1,zi1) &
    +2.*CoeffJijsm(5)*Fields(2,1,x0,y0,z0)*Fields(2,1,xi1,yi1,zi1) &
    +2.*CoeffJijsm(6)*Fields(1,1,x1,y0,z1)*Fields(3,1,x0,y0,z0) &
    -2.*CoeffJijsm(6)*Fields(1,1,x1,y0,zi1)*Fields(3,1,x0,y0,z0) &
    +2.*CoeffJijsm(7)*Fields(1,1,x1,y1,z1)*Fields(3,1,x0,y0,z0) &
    -2.*CoeffJijsm(7)*Fields(1,1,x1,y1,zi1)*Fields(3,1,x0,y0,z0) &
    +2.*CoeffJijsm(7)*Fields(1,1,x1,yi1,z1)*Fields(3,1,x0,y0,z0) &
    -2.*CoeffJijsm(7)*Fields(1,1,x1,yi1,zi1)*Fields(3,1,x0,y0,z0) &
    -2.*CoeffJijsm(6)*Fields(1,1,xi1,y0,z1)*Fields(3,1,x0,y0,z0) &
    +2.*CoeffJijsm(6)*Fields(1,1,xi1,y0,zi1)*Fields(3,1,x0,y0,z0) &
    -2.*CoeffJijsm(7)*Fields(1,1,xi1,y1,z1)*Fields(3,1,x0,y0,z0) &
    +2.*CoeffJijsm(7)*Fields(1,1,xi1,y1,zi1)*Fields(3,1,x0,y0,z0) &
    -2.*CoeffJijsm(7)*Fields(1,1,xi1,yi1,z1)*Fields(3,1,x0,y0,z0) &
    +2.*CoeffJijsm(7)*Fields(1,1,xi1,yi1,zi1)*Fields(3,1,x0,y0,z0) &
    +2.*CoeffJijsm(6)*Fields(2,1,x0,y1,z1)*Fields(3,1,x0,y0,z0) &
    -2.*CoeffJijsm(6)*Fields(2,1,x0,y1,zi1)*Fields(3,1,x0,y0,z0) &
    -2.*CoeffJijsm(6)*Fields(2,1,x0,yi1,z1)*Fields(3,1,x0,y0,z0) &
    +2.*CoeffJijsm(6)*Fields(2,1,x0,yi1,zi1)*Fields(3,1,x0,y0,z0) &
    +2.*CoeffJijsm(7)*Fields(2,1,x1,y1,z1)*Fields(3,1,x0,y0,z0) &
    -2.*CoeffJijsm(7)*Fields(2,1,x1,y1,zi1)*Fields(3,1,x0,y0,z0) &
    -2.*CoeffJijsm(7)*Fields(2,1,x1,yi1,z1)*Fields(3,1,x0,y0,z0) &
    +2.*CoeffJijsm(7)*Fields(2,1,x1,yi1,zi1)*Fields(3,1,x0,y0,z0) &
    +2.*CoeffJijsm(7)*Fields(2,1,xi1,y1,z1)*Fields(3,1,x0,y0,z0) &
    -2.*CoeffJijsm(7)*Fields(2,1,xi1,y1,zi1)*Fields(3,1,x0,y0,z0) &
    -2.*CoeffJijsm(7)*Fields(2,1,xi1,yi1,z1)*Fields(3,1,x0,y0,z0) &
    +2.*CoeffJijsm(7)*Fields(2,1,xi1,yi1,zi1)*Fields(3,1,x0,y0,z0) &
    +2.*CoeffJijsm(1)*Fields(3,1,x0,y0,z0)*Fields(3,1,x0,y0,z1) &
    +2.*CoeffJijsm(1)*Fields(3,1,x0,y0,z0)*Fields(3,1,x0,y0,zi1) &
    +2.*CoeffJijsm(2)*Fields(3,1,x0,y0,z0)*Fields(3,1,x0,y1,z0) &
    +2.*CoeffJijsm(6)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,y1,z1) &
    +2.*CoeffJijsm(3)*Fields(3,1,x0,y0,z0)*Fields(3,1,x0,y1,z1) &
    -2.*CoeffJijsm(6)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,y1,zi1) &
    +2.*CoeffJijsm(3)*Fields(3,1,x0,y0,z0)*Fields(3,1,x0,y1,zi1) &
    +2.*CoeffJijsm(2)*Fields(3,1,x0,y0,z0)*Fields(3,1,x0,yi1,z0) &
    -2.*CoeffJijsm(6)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,yi1,z1) &
    +2.*CoeffJijsm(3)*Fields(3,1,x0,y0,z0)*Fields(3,1,x0,yi1,z1) &
    +2.*CoeffJijsm(6)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,yi1,zi1) &
    +2.*CoeffJijsm(3)*Fields(3,1,x0,y0,z0)*Fields(3,1,x0,yi1,zi1) &
    +2.*CoeffJijsm(2)*Fields(3,1,x0,y0,z0)*Fields(3,1,x1,y0,z0) &
    +2.*CoeffJijsm(6)*Fields(1,1,x0,y0,z0)*Fields(3,1,x1,y0,z1) &
    +2.*CoeffJijsm(3)*Fields(3,1,x0,y0,z0)*Fields(3,1,x1,y0,z1) &
    -2.*CoeffJijsm(6)*Fields(1,1,x0,y0,z0)*Fields(3,1,x1,y0,zi1) &
    +2.*CoeffJijsm(3)*Fields(3,1,x0,y0,z0)*Fields(3,1,x1,y0,zi1) &
    +2.*CoeffJijsm(4)*Fields(3,1,x0,y0,z0)*Fields(3,1,x1,y1,z0) &
    +2.*CoeffJijsm(7)*Fields(1,1,x0,y0,z0)*Fields(3,1,x1,y1,z1) &
    +2.*CoeffJijsm(7)*Fields(2,1,x0,y0,z0)*Fields(3,1,x1,y1,z1) &
    +2.*CoeffJijsm(5)*Fields(3,1,x0,y0,z0)*Fields(3,1,x1,y1,z1) &
    -2.*CoeffJijsm(7)*Fields(1,1,x0,y0,z0)*Fields(3,1,x1,y1,zi1) &
    -2.*CoeffJijsm(7)*Fields(2,1,x0,y0,z0)*Fields(3,1,x1,y1,zi1) &
    +2.*CoeffJijsm(5)*Fields(3,1,x0,y0,z0)*Fields(3,1,x1,y1,zi1) &
    +2.*CoeffJijsm(4)*Fields(3,1,x0,y0,z0)*Fields(3,1,x1,yi1,z0) &
    +2.*CoeffJijsm(7)*Fields(1,1,x0,y0,z0)*Fields(3,1,x1,yi1,z1) &
    -2.*CoeffJijsm(7)*Fields(2,1,x0,y0,z0)*Fields(3,1,x1,yi1,z1) &
    +2.*CoeffJijsm(5)*Fields(3,1,x0,y0,z0)*Fields(3,1,x1,yi1,z1) &
    -2.*CoeffJijsm(7)*Fields(1,1,x0,y0,z0)*Fields(3,1,x1,yi1,zi1) &
    +2.*CoeffJijsm(7)*Fields(2,1,x0,y0,z0)*Fields(3,1,x1,yi1,zi1) &
    +2.*CoeffJijsm(5)*Fields(3,1,x0,y0,z0)*Fields(3,1,x1,yi1,zi1) &
    +2.*CoeffJijsm(2)*Fields(3,1,x0,y0,z0)*Fields(3,1,xi1,y0,z0) &
    -2.*CoeffJijsm(6)*Fields(1,1,x0,y0,z0)*Fields(3,1,xi1,y0,z1) &
    +2.*CoeffJijsm(3)*Fields(3,1,x0,y0,z0)*Fields(3,1,xi1,y0,z1) &
    +2.*CoeffJijsm(6)*Fields(1,1,x0,y0,z0)*Fields(3,1,xi1,y0,zi1) &
    +2.*CoeffJijsm(3)*Fields(3,1,x0,y0,z0)*Fields(3,1,xi1,y0,zi1) &
    +2.*CoeffJijsm(4)*Fields(3,1,x0,y0,z0)*Fields(3,1,xi1,y1,z0) &
    -2.*CoeffJijsm(7)*Fields(1,1,x0,y0,z0)*Fields(3,1,xi1,y1,z1) &
    +2.*CoeffJijsm(7)*Fields(2,1,x0,y0,z0)*Fields(3,1,xi1,y1,z1) &
    +2.*CoeffJijsm(5)*Fields(3,1,x0,y0,z0)*Fields(3,1,xi1,y1,z1) &
    +2.*CoeffJijsm(7)*Fields(1,1,x0,y0,z0)*Fields(3,1,xi1,y1,zi1) &
    -2.*CoeffJijsm(7)*Fields(2,1,x0,y0,z0)*Fields(3,1,xi1,y1,zi1) &
    +2.*CoeffJijsm(5)*Fields(3,1,x0,y0,z0)*Fields(3,1,xi1,y1,zi1) &
    +2.*CoeffJijsm(4)*Fields(3,1,x0,y0,z0)*Fields(3,1,xi1,yi1,z0) &
    -2.*CoeffJijsm(7)*Fields(1,1,x0,y0,z0)*Fields(3,1,xi1,yi1,z1) &
    -2.*CoeffJijsm(7)*Fields(2,1,x0,y0,z0)*Fields(3,1,xi1,yi1,z1) &
    +2.*CoeffJijsm(5)*Fields(3,1,x0,y0,z0)*Fields(3,1,xi1,yi1,z1) &
    +2.*CoeffJijsm(7)*Fields(1,1,x0,y0,z0)*Fields(3,1,xi1,yi1,zi1) &
    +2.*CoeffJijsm(7)*Fields(2,1,x0,y0,z0)*Fields(3,1,xi1,yi1,zi1) &
    +2.*CoeffJijsm(5)*Fields(3,1,x0,y0,z0)*Fields(3,1,xi1,yi1,zi1)
  
  EOnSiteJijhf = &
    EOnSiteJijhf &
    +2.*CoeffJijhf(2)*Fields(1,2,x0,y0,z0)*Fields(1,2,x0,y0,z1) &
    +2.*CoeffJijhf(2)*Fields(1,2,x0,y0,z0)*Fields(1,2,x0,y0,zi1) &
    +2.*CoeffJijhf(2)*Fields(1,2,x0,y0,z0)*Fields(1,2,x0,y1,z0) &
    +2.*CoeffJijhf(4)*Fields(1,2,x0,y0,z0)*Fields(1,2,x0,y1,z1) &
    +2.*CoeffJijhf(4)*Fields(1,2,x0,y0,z0)*Fields(1,2,x0,y1,zi1) &
    +2.*CoeffJijhf(2)*Fields(1,2,x0,y0,z0)*Fields(1,2,x0,yi1,z0) &
    +2.*CoeffJijhf(4)*Fields(1,2,x0,y0,z0)*Fields(1,2,x0,yi1,z1) &
    +2.*CoeffJijhf(4)*Fields(1,2,x0,y0,z0)*Fields(1,2,x0,yi1,zi1) &
    +2.*CoeffJijhf(1)*Fields(1,2,x0,y0,z0)*Fields(1,2,x1,y0,z0) &
    +2.*CoeffJijhf(3)*Fields(1,2,x0,y0,z0)*Fields(1,2,x1,y0,z1) &
    +2.*CoeffJijhf(3)*Fields(1,2,x0,y0,z0)*Fields(1,2,x1,y0,zi1) &
    +2.*CoeffJijhf(3)*Fields(1,2,x0,y0,z0)*Fields(1,2,x1,y1,z0) &
    +2.*CoeffJijhf(5)*Fields(1,2,x0,y0,z0)*Fields(1,2,x1,y1,z1) &
    +2.*CoeffJijhf(5)*Fields(1,2,x0,y0,z0)*Fields(1,2,x1,y1,zi1) &
    +2.*CoeffJijhf(3)*Fields(1,2,x0,y0,z0)*Fields(1,2,x1,yi1,z0) &
    +2.*CoeffJijhf(5)*Fields(1,2,x0,y0,z0)*Fields(1,2,x1,yi1,z1) &
    +2.*CoeffJijhf(5)*Fields(1,2,x0,y0,z0)*Fields(1,2,x1,yi1,zi1) &
    +2.*CoeffJijhf(1)*Fields(1,2,x0,y0,z0)*Fields(1,2,xi1,y0,z0) &
    +2.*CoeffJijhf(3)*Fields(1,2,x0,y0,z0)*Fields(1,2,xi1,y0,z1) &
    +2.*CoeffJijhf(3)*Fields(1,2,x0,y0,z0)*Fields(1,2,xi1,y0,zi1) &
    +2.*CoeffJijhf(3)*Fields(1,2,x0,y0,z0)*Fields(1,2,xi1,y1,z0) &
    +2.*CoeffJijhf(5)*Fields(1,2,x0,y0,z0)*Fields(1,2,xi1,y1,z1) &
    +2.*CoeffJijhf(5)*Fields(1,2,x0,y0,z0)*Fields(1,2,xi1,y1,zi1) &
    +2.*CoeffJijhf(3)*Fields(1,2,x0,y0,z0)*Fields(1,2,xi1,yi1,z0) &
    +2.*CoeffJijhf(5)*Fields(1,2,x0,y0,z0)*Fields(1,2,xi1,yi1,z1) &
    +2.*CoeffJijhf(5)*Fields(1,2,x0,y0,z0)*Fields(1,2,xi1,yi1,zi1) &
    +2.*CoeffJijhf(6)*Fields(1,2,x1,y1,z0)*Fields(2,2,x0,y0,z0) &
    +2.*CoeffJijhf(7)*Fields(1,2,x1,y1,z1)*Fields(2,2,x0,y0,z0) &
    +2.*CoeffJijhf(7)*Fields(1,2,x1,y1,zi1)*Fields(2,2,x0,y0,z0) &
    -2.*CoeffJijhf(6)*Fields(1,2,x1,yi1,z0)*Fields(2,2,x0,y0,z0) &
    -2.*CoeffJijhf(7)*Fields(1,2,x1,yi1,z1)*Fields(2,2,x0,y0,z0) &
    -2.*CoeffJijhf(7)*Fields(1,2,x1,yi1,zi1)*Fields(2,2,x0,y0,z0) &
    -2.*CoeffJijhf(6)*Fields(1,2,xi1,y1,z0)*Fields(2,2,x0,y0,z0) &
    -2.*CoeffJijhf(7)*Fields(1,2,xi1,y1,z1)*Fields(2,2,x0,y0,z0) &
    -2.*CoeffJijhf(7)*Fields(1,2,xi1,y1,zi1)*Fields(2,2,x0,y0,z0) &
    +2.*CoeffJijhf(6)*Fields(1,2,xi1,yi1,z0)*Fields(2,2,x0,y0,z0) &
    +2.*CoeffJijhf(7)*Fields(1,2,xi1,yi1,z1)*Fields(2,2,x0,y0,z0) &
    +2.*CoeffJijhf(7)*Fields(1,2,xi1,yi1,zi1)*Fields(2,2,x0,y0,z0) &
    +2.*CoeffJijhf(2)*Fields(2,2,x0,y0,z0)*Fields(2,2,x0,y0,z1) &
    +2.*CoeffJijhf(2)*Fields(2,2,x0,y0,z0)*Fields(2,2,x0,y0,zi1) &
    +2.*CoeffJijhf(1)*Fields(2,2,x0,y0,z0)*Fields(2,2,x0,y1,z0) &
    +2.*CoeffJijhf(3)*Fields(2,2,x0,y0,z0)*Fields(2,2,x0,y1,z1) &
    +2.*CoeffJijhf(3)*Fields(2,2,x0,y0,z0)*Fields(2,2,x0,y1,zi1) &
    +2.*CoeffJijhf(1)*Fields(2,2,x0,y0,z0)*Fields(2,2,x0,yi1,z0) &
    +2.*CoeffJijhf(3)*Fields(2,2,x0,y0,z0)*Fields(2,2,x0,yi1,z1) &
    +2.*CoeffJijhf(3)*Fields(2,2,x0,y0,z0)*Fields(2,2,x0,yi1,zi1) &
    +2.*CoeffJijhf(2)*Fields(2,2,x0,y0,z0)*Fields(2,2,x1,y0,z0) &
    +2.*CoeffJijhf(4)*Fields(2,2,x0,y0,z0)*Fields(2,2,x1,y0,z1) &
    +2.*CoeffJijhf(4)*Fields(2,2,x0,y0,z0)*Fields(2,2,x1,y0,zi1) &
    +2.*CoeffJijhf(6)*Fields(1,2,x0,y0,z0)*Fields(2,2,x1,y1,z0) &
    +2.*CoeffJijhf(3)*Fields(2,2,x0,y0,z0)*Fields(2,2,x1,y1,z0) &
    +2.*CoeffJijhf(7)*Fields(1,2,x0,y0,z0)*Fields(2,2,x1,y1,z1) &
    +2.*CoeffJijhf(5)*Fields(2,2,x0,y0,z0)*Fields(2,2,x1,y1,z1) &
    +2.*CoeffJijhf(7)*Fields(1,2,x0,y0,z0)*Fields(2,2,x1,y1,zi1) &
    +2.*CoeffJijhf(5)*Fields(2,2,x0,y0,z0)*Fields(2,2,x1,y1,zi1) &
    -2.*CoeffJijhf(6)*Fields(1,2,x0,y0,z0)*Fields(2,2,x1,yi1,z0) &
    +2.*CoeffJijhf(3)*Fields(2,2,x0,y0,z0)*Fields(2,2,x1,yi1,z0) &
    -2.*CoeffJijhf(7)*Fields(1,2,x0,y0,z0)*Fields(2,2,x1,yi1,z1) &
    +2.*CoeffJijhf(5)*Fields(2,2,x0,y0,z0)*Fields(2,2,x1,yi1,z1) &
    -2.*CoeffJijhf(7)*Fields(1,2,x0,y0,z0)*Fields(2,2,x1,yi1,zi1) &
    +2.*CoeffJijhf(5)*Fields(2,2,x0,y0,z0)*Fields(2,2,x1,yi1,zi1) &
    +2.*CoeffJijhf(2)*Fields(2,2,x0,y0,z0)*Fields(2,2,xi1,y0,z0) &
    +2.*CoeffJijhf(4)*Fields(2,2,x0,y0,z0)*Fields(2,2,xi1,y0,z1) &
    +2.*CoeffJijhf(4)*Fields(2,2,x0,y0,z0)*Fields(2,2,xi1,y0,zi1) &
    -2.*CoeffJijhf(6)*Fields(1,2,x0,y0,z0)*Fields(2,2,xi1,y1,z0) &
    +2.*CoeffJijhf(3)*Fields(2,2,x0,y0,z0)*Fields(2,2,xi1,y1,z0) &
    -2.*CoeffJijhf(7)*Fields(1,2,x0,y0,z0)*Fields(2,2,xi1,y1,z1) &
    +2.*CoeffJijhf(5)*Fields(2,2,x0,y0,z0)*Fields(2,2,xi1,y1,z1) &
    -2.*CoeffJijhf(7)*Fields(1,2,x0,y0,z0)*Fields(2,2,xi1,y1,zi1) &
    +2.*CoeffJijhf(5)*Fields(2,2,x0,y0,z0)*Fields(2,2,xi1,y1,zi1) &
    +2.*CoeffJijhf(6)*Fields(1,2,x0,y0,z0)*Fields(2,2,xi1,yi1,z0) &
    +2.*CoeffJijhf(3)*Fields(2,2,x0,y0,z0)*Fields(2,2,xi1,yi1,z0) &
    +2.*CoeffJijhf(7)*Fields(1,2,x0,y0,z0)*Fields(2,2,xi1,yi1,z1) &
    +2.*CoeffJijhf(5)*Fields(2,2,x0,y0,z0)*Fields(2,2,xi1,yi1,z1) &
    +2.*CoeffJijhf(7)*Fields(1,2,x0,y0,z0)*Fields(2,2,xi1,yi1,zi1) &
    +2.*CoeffJijhf(5)*Fields(2,2,x0,y0,z0)*Fields(2,2,xi1,yi1,zi1) &
    +2.*CoeffJijhf(6)*Fields(1,2,x1,y0,z1)*Fields(3,2,x0,y0,z0) &
    -2.*CoeffJijhf(6)*Fields(1,2,x1,y0,zi1)*Fields(3,2,x0,y0,z0) &
    +2.*CoeffJijhf(7)*Fields(1,2,x1,y1,z1)*Fields(3,2,x0,y0,z0) &
    -2.*CoeffJijhf(7)*Fields(1,2,x1,y1,zi1)*Fields(3,2,x0,y0,z0) &
    +2.*CoeffJijhf(7)*Fields(1,2,x1,yi1,z1)*Fields(3,2,x0,y0,z0) &
    -2.*CoeffJijhf(7)*Fields(1,2,x1,yi1,zi1)*Fields(3,2,x0,y0,z0) &
    -2.*CoeffJijhf(6)*Fields(1,2,xi1,y0,z1)*Fields(3,2,x0,y0,z0) &
    +2.*CoeffJijhf(6)*Fields(1,2,xi1,y0,zi1)*Fields(3,2,x0,y0,z0) &
    -2.*CoeffJijhf(7)*Fields(1,2,xi1,y1,z1)*Fields(3,2,x0,y0,z0) &
    +2.*CoeffJijhf(7)*Fields(1,2,xi1,y1,zi1)*Fields(3,2,x0,y0,z0) &
    -2.*CoeffJijhf(7)*Fields(1,2,xi1,yi1,z1)*Fields(3,2,x0,y0,z0) &
    +2.*CoeffJijhf(7)*Fields(1,2,xi1,yi1,zi1)*Fields(3,2,x0,y0,z0) &
    +2.*CoeffJijhf(6)*Fields(2,2,x0,y1,z1)*Fields(3,2,x0,y0,z0) &
    -2.*CoeffJijhf(6)*Fields(2,2,x0,y1,zi1)*Fields(3,2,x0,y0,z0) &
    -2.*CoeffJijhf(6)*Fields(2,2,x0,yi1,z1)*Fields(3,2,x0,y0,z0) &
    +2.*CoeffJijhf(6)*Fields(2,2,x0,yi1,zi1)*Fields(3,2,x0,y0,z0) &
    +2.*CoeffJijhf(7)*Fields(2,2,x1,y1,z1)*Fields(3,2,x0,y0,z0) &
    -2.*CoeffJijhf(7)*Fields(2,2,x1,y1,zi1)*Fields(3,2,x0,y0,z0) &
    -2.*CoeffJijhf(7)*Fields(2,2,x1,yi1,z1)*Fields(3,2,x0,y0,z0) &
    +2.*CoeffJijhf(7)*Fields(2,2,x1,yi1,zi1)*Fields(3,2,x0,y0,z0) &
    +2.*CoeffJijhf(7)*Fields(2,2,xi1,y1,z1)*Fields(3,2,x0,y0,z0) &
    -2.*CoeffJijhf(7)*Fields(2,2,xi1,y1,zi1)*Fields(3,2,x0,y0,z0) &
    -2.*CoeffJijhf(7)*Fields(2,2,xi1,yi1,z1)*Fields(3,2,x0,y0,z0) &
    +2.*CoeffJijhf(7)*Fields(2,2,xi1,yi1,zi1)*Fields(3,2,x0,y0,z0) &
    +2.*CoeffJijhf(1)*Fields(3,2,x0,y0,z0)*Fields(3,2,x0,y0,z1) &
    +2.*CoeffJijhf(1)*Fields(3,2,x0,y0,z0)*Fields(3,2,x0,y0,zi1) &
    +2.*CoeffJijhf(2)*Fields(3,2,x0,y0,z0)*Fields(3,2,x0,y1,z0) &
    +2.*CoeffJijhf(6)*Fields(2,2,x0,y0,z0)*Fields(3,2,x0,y1,z1) &
    +2.*CoeffJijhf(3)*Fields(3,2,x0,y0,z0)*Fields(3,2,x0,y1,z1) &
    -2.*CoeffJijhf(6)*Fields(2,2,x0,y0,z0)*Fields(3,2,x0,y1,zi1) &
    +2.*CoeffJijhf(3)*Fields(3,2,x0,y0,z0)*Fields(3,2,x0,y1,zi1) &
    +2.*CoeffJijhf(2)*Fields(3,2,x0,y0,z0)*Fields(3,2,x0,yi1,z0) &
    -2.*CoeffJijhf(6)*Fields(2,2,x0,y0,z0)*Fields(3,2,x0,yi1,z1) &
    +2.*CoeffJijhf(3)*Fields(3,2,x0,y0,z0)*Fields(3,2,x0,yi1,z1) &
    +2.*CoeffJijhf(6)*Fields(2,2,x0,y0,z0)*Fields(3,2,x0,yi1,zi1) &
    +2.*CoeffJijhf(3)*Fields(3,2,x0,y0,z0)*Fields(3,2,x0,yi1,zi1) &
    +2.*CoeffJijhf(2)*Fields(3,2,x0,y0,z0)*Fields(3,2,x1,y0,z0) &
    +2.*CoeffJijhf(6)*Fields(1,2,x0,y0,z0)*Fields(3,2,x1,y0,z1) &
    +2.*CoeffJijhf(3)*Fields(3,2,x0,y0,z0)*Fields(3,2,x1,y0,z1) &
    -2.*CoeffJijhf(6)*Fields(1,2,x0,y0,z0)*Fields(3,2,x1,y0,zi1) &
    +2.*CoeffJijhf(3)*Fields(3,2,x0,y0,z0)*Fields(3,2,x1,y0,zi1) &
    +2.*CoeffJijhf(4)*Fields(3,2,x0,y0,z0)*Fields(3,2,x1,y1,z0) &
    +2.*CoeffJijhf(7)*Fields(1,2,x0,y0,z0)*Fields(3,2,x1,y1,z1) &
    +2.*CoeffJijhf(7)*Fields(2,2,x0,y0,z0)*Fields(3,2,x1,y1,z1) &
    +2.*CoeffJijhf(5)*Fields(3,2,x0,y0,z0)*Fields(3,2,x1,y1,z1) &
    -2.*CoeffJijhf(7)*Fields(1,2,x0,y0,z0)*Fields(3,2,x1,y1,zi1) &
    -2.*CoeffJijhf(7)*Fields(2,2,x0,y0,z0)*Fields(3,2,x1,y1,zi1) &
    +2.*CoeffJijhf(5)*Fields(3,2,x0,y0,z0)*Fields(3,2,x1,y1,zi1) &
    +2.*CoeffJijhf(4)*Fields(3,2,x0,y0,z0)*Fields(3,2,x1,yi1,z0) &
    +2.*CoeffJijhf(7)*Fields(1,2,x0,y0,z0)*Fields(3,2,x1,yi1,z1) &
    -2.*CoeffJijhf(7)*Fields(2,2,x0,y0,z0)*Fields(3,2,x1,yi1,z1) &
    +2.*CoeffJijhf(5)*Fields(3,2,x0,y0,z0)*Fields(3,2,x1,yi1,z1) &
    -2.*CoeffJijhf(7)*Fields(1,2,x0,y0,z0)*Fields(3,2,x1,yi1,zi1) &
    +2.*CoeffJijhf(7)*Fields(2,2,x0,y0,z0)*Fields(3,2,x1,yi1,zi1) &
    +2.*CoeffJijhf(5)*Fields(3,2,x0,y0,z0)*Fields(3,2,x1,yi1,zi1) &
    +2.*CoeffJijhf(2)*Fields(3,2,x0,y0,z0)*Fields(3,2,xi1,y0,z0) &
    -2.*CoeffJijhf(6)*Fields(1,2,x0,y0,z0)*Fields(3,2,xi1,y0,z1) &
    +2.*CoeffJijhf(3)*Fields(3,2,x0,y0,z0)*Fields(3,2,xi1,y0,z1) &
    +2.*CoeffJijhf(6)*Fields(1,2,x0,y0,z0)*Fields(3,2,xi1,y0,zi1) &
    +2.*CoeffJijhf(3)*Fields(3,2,x0,y0,z0)*Fields(3,2,xi1,y0,zi1) &
    +2.*CoeffJijhf(4)*Fields(3,2,x0,y0,z0)*Fields(3,2,xi1,y1,z0) &
    -2.*CoeffJijhf(7)*Fields(1,2,x0,y0,z0)*Fields(3,2,xi1,y1,z1) &
    +2.*CoeffJijhf(7)*Fields(2,2,x0,y0,z0)*Fields(3,2,xi1,y1,z1) &
    +2.*CoeffJijhf(5)*Fields(3,2,x0,y0,z0)*Fields(3,2,xi1,y1,z1) &
    +2.*CoeffJijhf(7)*Fields(1,2,x0,y0,z0)*Fields(3,2,xi1,y1,zi1) &
    -2.*CoeffJijhf(7)*Fields(2,2,x0,y0,z0)*Fields(3,2,xi1,y1,zi1) &
    +2.*CoeffJijhf(5)*Fields(3,2,x0,y0,z0)*Fields(3,2,xi1,y1,zi1) &
    +2.*CoeffJijhf(4)*Fields(3,2,x0,y0,z0)*Fields(3,2,xi1,yi1,z0) &
    -2.*CoeffJijhf(7)*Fields(1,2,x0,y0,z0)*Fields(3,2,xi1,yi1,z1) &
    -2.*CoeffJijhf(7)*Fields(2,2,x0,y0,z0)*Fields(3,2,xi1,yi1,z1) &
    +2.*CoeffJijhf(5)*Fields(3,2,x0,y0,z0)*Fields(3,2,xi1,yi1,z1) &
    +2.*CoeffJijhf(7)*Fields(1,2,x0,y0,z0)*Fields(3,2,xi1,yi1,zi1) &
    +2.*CoeffJijhf(7)*Fields(2,2,x0,y0,z0)*Fields(3,2,xi1,yi1,zi1) &
    +2.*CoeffJijhf(5)*Fields(3,2,x0,y0,z0)*Fields(3,2,xi1,yi1,zi1)
  
  EOnSiteEpsDisp = &
    EOnSiteEpsDisp &
    +CoeffEpsDisp(5)*e0ij(1,1)*Fields(1,1,x0,y0,z0)**2 &
    +CoeffEpsDisp(4)*e0ij(2,2)*Fields(1,1,x0,y0,z0)**2 &
    +CoeffEpsDisp(4)*e0ij(3,3)*Fields(1,1,x0,y0,z0)**2 &
    +CoeffEpsDisp(5)*euij(1,1)*Fields(1,1,x0,y0,z0)**2 &
    +CoeffEpsDisp(4)*euij(2,2)*Fields(1,1,x0,y0,z0)**2 &
    +CoeffEpsDisp(4)*euij(3,3)*Fields(1,1,x0,y0,z0)**2 &
    +CoeffEpsDisp(8)*e0ij(1,1)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
    +CoeffEpsDisp(6)*e0ij(2,2)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
    +CoeffEpsDisp(6)*e0ij(3,3)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
    +CoeffEpsDisp(8)*euij(1,1)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
    +CoeffEpsDisp(6)*euij(2,2)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
    +CoeffEpsDisp(6)*euij(3,3)*Fields(1,1,x0,y0,z0)*Fields(1,2,x0,y0,z0) &
    +CoeffEpsDisp(9)*e0ij(1,1)*Fields(1,2,x0,y0,z0)**2 &
    +CoeffEpsDisp(7)*e0ij(2,2)*Fields(1,2,x0,y0,z0)**2 &
    +CoeffEpsDisp(7)*e0ij(3,3)*Fields(1,2,x0,y0,z0)**2 &
    +CoeffEpsDisp(9)*euij(1,1)*Fields(1,2,x0,y0,z0)**2 &
    +CoeffEpsDisp(7)*euij(2,2)*Fields(1,2,x0,y0,z0)**2 &
    +CoeffEpsDisp(7)*euij(3,3)*Fields(1,2,x0,y0,z0)**2 &
    +CoeffEpsDisp(1)*e0ij(1,2)*Fields(1,1,x0,y0,z0)*Fields(2,1,x0,y0,z0) &
    +CoeffEpsDisp(1)*euij(1,2)*Fields(1,1,x0,y0,z0)*Fields(2,1,x0,y0,z0) &
    +CoeffEpsDisp(2)*e0ij(1,2)*Fields(1,2,x0,y0,z0)*Fields(2,1,x0,y0,z0) &
    +CoeffEpsDisp(2)*euij(1,2)*Fields(1,2,x0,y0,z0)*Fields(2,1,x0,y0,z0) &
    +CoeffEpsDisp(4)*e0ij(1,1)*Fields(2,1,x0,y0,z0)**2 &
    +CoeffEpsDisp(5)*e0ij(2,2)*Fields(2,1,x0,y0,z0)**2 &
    +CoeffEpsDisp(4)*e0ij(3,3)*Fields(2,1,x0,y0,z0)**2 &
    +CoeffEpsDisp(4)*euij(1,1)*Fields(2,1,x0,y0,z0)**2 &
    +CoeffEpsDisp(5)*euij(2,2)*Fields(2,1,x0,y0,z0)**2 &
    +CoeffEpsDisp(4)*euij(3,3)*Fields(2,1,x0,y0,z0)**2 &
    +CoeffEpsDisp(2)*e0ij(1,2)*Fields(1,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +CoeffEpsDisp(2)*euij(1,2)*Fields(1,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +CoeffEpsDisp(3)*e0ij(1,2)*Fields(1,2,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +CoeffEpsDisp(3)*euij(1,2)*Fields(1,2,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +CoeffEpsDisp(6)*e0ij(1,1)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +CoeffEpsDisp(8)*e0ij(2,2)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +CoeffEpsDisp(6)*e0ij(3,3)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +CoeffEpsDisp(6)*euij(1,1)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +CoeffEpsDisp(8)*euij(2,2)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +CoeffEpsDisp(6)*euij(3,3)*Fields(2,1,x0,y0,z0)*Fields(2,2,x0,y0,z0) &
    +CoeffEpsDisp(7)*e0ij(1,1)*Fields(2,2,x0,y0,z0)**2 &
    +CoeffEpsDisp(9)*e0ij(2,2)*Fields(2,2,x0,y0,z0)**2 &
    +CoeffEpsDisp(7)*e0ij(3,3)*Fields(2,2,x0,y0,z0)**2 &
    +CoeffEpsDisp(7)*euij(1,1)*Fields(2,2,x0,y0,z0)**2 &
    +CoeffEpsDisp(9)*euij(2,2)*Fields(2,2,x0,y0,z0)**2 &
    +CoeffEpsDisp(7)*euij(3,3)*Fields(2,2,x0,y0,z0)**2 &
    +CoeffEpsDisp(1)*e0ij(1,3)*Fields(1,1,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
    +CoeffEpsDisp(1)*euij(1,3)*Fields(1,1,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
    +CoeffEpsDisp(2)*e0ij(1,3)*Fields(1,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
    +CoeffEpsDisp(2)*euij(1,3)*Fields(1,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
    +CoeffEpsDisp(1)*e0ij(2,3)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
    +CoeffEpsDisp(1)*euij(2,3)*Fields(2,1,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
    +CoeffEpsDisp(2)*e0ij(2,3)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
    +CoeffEpsDisp(2)*euij(2,3)*Fields(2,2,x0,y0,z0)*Fields(3,1,x0,y0,z0) &
    +CoeffEpsDisp(4)*e0ij(1,1)*Fields(3,1,x0,y0,z0)**2 &
    +CoeffEpsDisp(4)*e0ij(2,2)*Fields(3,1,x0,y0,z0)**2 &
    +CoeffEpsDisp(5)*e0ij(3,3)*Fields(3,1,x0,y0,z0)**2 &
    +CoeffEpsDisp(4)*euij(1,1)*Fields(3,1,x0,y0,z0)**2 &
    +CoeffEpsDisp(4)*euij(2,2)*Fields(3,1,x0,y0,z0)**2 &
    +CoeffEpsDisp(5)*euij(3,3)*Fields(3,1,x0,y0,z0)**2 &
    +CoeffEpsDisp(2)*e0ij(1,3)*Fields(1,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +CoeffEpsDisp(2)*euij(1,3)*Fields(1,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +CoeffEpsDisp(3)*e0ij(1,3)*Fields(1,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +CoeffEpsDisp(3)*euij(1,3)*Fields(1,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +CoeffEpsDisp(2)*e0ij(2,3)*Fields(2,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +CoeffEpsDisp(2)*euij(2,3)*Fields(2,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +CoeffEpsDisp(3)*e0ij(2,3)*Fields(2,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +CoeffEpsDisp(3)*euij(2,3)*Fields(2,2,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +CoeffEpsDisp(6)*e0ij(1,1)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +CoeffEpsDisp(6)*e0ij(2,2)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +CoeffEpsDisp(8)*e0ij(3,3)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +CoeffEpsDisp(6)*euij(1,1)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +CoeffEpsDisp(6)*euij(2,2)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +CoeffEpsDisp(8)*euij(3,3)*Fields(3,1,x0,y0,z0)*Fields(3,2,x0,y0,z0) &
    +CoeffEpsDisp(7)*e0ij(1,1)*Fields(3,2,x0,y0,z0)**2 &
    +CoeffEpsDisp(7)*e0ij(2,2)*Fields(3,2,x0,y0,z0)**2 &
    +CoeffEpsDisp(9)*e0ij(3,3)*Fields(3,2,x0,y0,z0)**2 &
    +CoeffEpsDisp(7)*euij(1,1)*Fields(3,2,x0,y0,z0)**2 &
    +CoeffEpsDisp(7)*euij(2,2)*Fields(3,2,x0,y0,z0)**2 &
    +CoeffEpsDisp(9)*euij(3,3)*Fields(3,2,x0,y0,z0)**2
  
  EOnSiteEps = &
    EOnSiteEps &
    +CoeffEps(1)*e0ij(1,1) &
    +0.5*CoeffEps(4)*e0ij(1,1)**2 &
    +0.5*CoeffEps(2)*e0ij(1,2)**2 &
    +0.5*CoeffEps(2)*e0ij(1,3)**2 &
    +CoeffEps(1)*e0ij(2,2) &
    +CoeffEps(3)*e0ij(1,1)*e0ij(2,2) &
    +0.5*CoeffEps(4)*e0ij(2,2)**2 &
    +0.5*CoeffEps(2)*e0ij(2,3)**2 &
    +CoeffEps(1)*e0ij(3,3) &
    +CoeffEps(3)*e0ij(1,1)*e0ij(3,3) &
    +CoeffEps(3)*e0ij(2,2)*e0ij(3,3) &
    +0.5*CoeffEps(4)*e0ij(3,3)**2 &
    +CoeffEps(1)*euij(1,1) &
    +0.5*CoeffEps(4)*euij(1,1)**2 &
    +0.5*CoeffEps(2)*euij(1,2)**2 &
    +0.5*CoeffEps(2)*euij(1,3)**2 &
    +CoeffEps(1)*euij(2,2) &
    +CoeffEps(3)*euij(1,1)*euij(2,2) &
    +0.5*CoeffEps(4)*euij(2,2)**2 &
    +0.5*CoeffEps(2)*euij(2,3)**2 &
    +CoeffEps(1)*euij(3,3) &
    +CoeffEps(3)*euij(1,1)*euij(3,3) &
    +CoeffEps(3)*euij(2,2)*euij(3,3) &
    +0.5*CoeffEps(4)*euij(3,3)**2
  
  EOnSite=EOnSiteDisp+EOnSiteJijsm+EOnSiteJijhf+EOnSiteEpsDisp+EOnSiteEps
  
End Function GetEOnSite
