Function GetEOnSiteEpsDisp(x0, y0, z0, Fields, e0ij) Result(EOnSiteEpsDisp)
  
  Implicit none
  Real*8,  Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(in) :: e0ij(3,3)
  Integer, Intent(in) :: x0, y0, z0
  Real*8              :: eij(3,3), euij(3,3)
  Real*8              :: EOnSiteEpsDisp
  
  Integer             :: x1,y1,z1,xi1,yi1,zi1
  
  EOnSiteEpsDisp = 0.0D0
  
  
  euij = StrainFromu(x0, y0, z0, Fields)
  eij = e0ij + euij
  
  x1 = (x0+1)-floor(real(x0+1-1)/real(cgrid%n1))*cgrid%n1
  y1 = (y0+1)-floor(real(y0+1-1)/real(cgrid%n2))*cgrid%n2
  z1 = (z0+1)-floor(real(z0+1-1)/real(cgrid%n3))*cgrid%n3
  xi1 = (x0-1)-floor(real(x0-1-1)/real(cgrid%n1))*cgrid%n1
  yi1 = (y0-1)-floor(real(y0-1-1)/real(cgrid%n2))*cgrid%n2
  zi1 = (z0-1)-floor(real(z0-1-1)/real(cgrid%n3))*cgrid%n3
  
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
  
  
End Function GetEOnSiteEpsDisp
