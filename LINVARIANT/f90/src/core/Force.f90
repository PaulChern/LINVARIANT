Module Force

  Use omp_lib
  Use Parameters
  Use Lattice
  Use Ewald
  Use Aux

  Implicit none
  Contains

  Include "ForcesExpr.f90"
 
  Include "GetForces.f90"
  Include "GetForcesEps.f90"
  Include "GetForcesG.f90"
  Include "GetForcesH.f90"
 
  Include "GetForcesGp.f90"
  Include "GetForcesGs.f90"
  Include "GetForcesGsp.f90"
 
  Include "GetForcesHp.f90"
  Include "GetForcesHu.f90"
  Include "GetForcesHup.f90"
  Include "GetForcesHsu.f90"
  Include "GetForcesHsp.f90"
  Include "GetForcesHsup.f90"

  Function GetGradient(Fields, e0ij) Result(g)
    Implicit none
    Real*8,  Intent(in)    :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8,  Intent(in)    :: e0ij(3,3)
    Real*8                 :: Forces(Max(FieldDim,6),NumField+1,cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8                 :: ForcesEta(6)
    Real*8                 :: g(FieldDim*NumField*cgrid%npts+6)
    Integer                :: i, ix, iy, iz

    Forces = 0.0D0
    If(DipoleQ) then
      Forces(1:3,:,:,:,:) = GetForcesEwald(Fields)
    End If
    Forces = -1.0D0*Forces - GetForces(Fields,e0ij)
    Do i = 1, 6
      ForcesEta(i) = Sum(Forces(i,4,:,:,:))
    End do
    
    g  = GetSysVector(Forces(1:FieldDim,1:NumField,:,:,:),eta2eij(ForcesEta))

  End Function GetGradient
End Module Force
