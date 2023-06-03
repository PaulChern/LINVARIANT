Subroutine GetPTSwapMap(TempList, EtotList, NumReplicas, Replicas)
  
  Use LINVARIANT
  Use Parameters
  Use Constants
  
  Implicit none
  Integer, Intent(in)    :: NumReplicas
  Real*8,  Intent(inout) :: TempList(NumReplicas)
  Integer, Intent(inout) :: Replicas(NumReplicas,5)
  Real*8,  Intent(inout) :: EtotList(NumReplicas)
  Integer                :: ireplica, TempReplica(3)
  Real*8                 :: dene, ProbTest, TempEtot
  Real*8                 :: rdum, AccProb
  
  do ireplica = 1, NumReplicas
    Replicas(ireplica,1) = ireplica
  end do
  
  do ireplica = 1, NumReplicas - 1
    dene = Hartree/k_bolt_ev*(1/TempList(ireplica+1)-1/TempList(ireplica))*(EtotList(ireplica+1)-EtotList(ireplica))
    AccProb = Min(1.0d0, Exp(dene))
    Call random_number(ProbTest)
    if(ProbTest .lt. AccProb) then
      TempReplica(:) = Replicas(ireplica,1:3) 
      Replicas(ireplica,1:3) = Replicas(ireplica+1,1:3)
      Replicas(ireplica+1,1:3) = TempReplica(:)
      TempEtot = EtotList(ireplica)
      EtotList(ireplica) = EtotList(ireplica+1)
      EtotList(ireplica+1) = TempEtot
    end if
  end do
  
  Replicas(1,3) = 1
  Replicas(NumReplicas,3) = -1
  do ireplica = 1, NumReplicas
    if(ireplica.ne.Replicas(ireplica,1)) then
      if(Replicas(ireplica,3).eq.1) then
        Replicas(ireplica,4) = Replicas(ireplica,4) + 1
      else if(Replicas(ireplica,3).eq.-1) then
        Replicas(ireplica,5) = Replicas(ireplica,5) + 1
      end if
    end if
  end do
  
End Subroutine GetPTSwapMap
  
  
Subroutine PTSwap(TempList, Fields, EwaldField, e0ij, dFieldsdt, de0ijdt, gm, NumReplicas, Replicas, ireplica, istep)
  
  Use LINVARIANT
  Use Constants
  Use Inputs
  use mpi
  
  Implicit none
  Real*8,  Intent(inout) :: TempList(NumReplicas)
  Integer, Intent(inout) :: Replicas(NumReplicas,5)
  Real*8,  Intent(inout) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(inout) :: dFieldsdt(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(inout) :: EwaldField(3, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8,  Intent(inout) :: e0ij(3,3)
  Real*8,  Intent(inout) :: de0ijdt(3,3)
  Real*8,  Intent(inout) :: gm(NumField+1)
  Integer, Intent(in)    :: NumReplicas, ireplica
  Integer                :: i, istep, DimFields1D 
  
  Real*8   :: MPI_EWALDFIELD_1D(cgrid%npts*NumField*FieldDim)
  Real*8   :: MPI_EWALDFIELD_1D_Recv(NumReplicas*cgrid%npts*NumField*FieldDim)
  Real*8   :: MPI_FIELD_1D(cgrid%npts*NumField*FieldDim)
  Real*8   :: MPI_FIELD_1D_Recv(NumReplicas*cgrid%npts*NumField*FieldDim)
  Real*8   :: MPI_ETA(6)
  Real*8   :: MPI_ETA_Recv(NumReplicas*6)
  Real*8   :: MPI_DFIELDDT_1D(cgrid%npts*NumField*FieldDim)
  Real*8   :: MPI_DFIELDDT_1D_Recv(NumReplicas*cgrid%npts*NumField*FieldDim)
  Real*8   :: MPI_DETADT(6)
  Real*8   :: MPI_DETADT_Recv(NumReplicas*6)
  Real*8   :: MPI_GM
  Real*8   :: MPI_GM_Recv(NumReplicas)
  Real*8   :: EtotRecv(NumReplicas)
  Real*8   :: Etot, Tk(2)
  
  DimFields1D = cgrid%npts*NumField*FieldDim
  
  if(CoolingSteps.eq.0)then
    Etot = GetEtot(Fields, e0ij)
    Tk = Thermometer(dFieldsdt, de0ijdt)
  else
    if(ireplica.eq.istep) then
      Etot = (ireplica-1)*1.0E6 
    elseif(ireplica.eq.istep-1)then
      Etot = (ireplica+1)*1.0E6
    else
      Etot = ireplica*1.0E6
    end if
  end if

  MPI_ETA = eij2eta(e0ij)
  MPI_FIELD_1D = AllFields1D(Fields)
  MPI_EWALDFIELD_1D = AllFields1D(EwaldField)
  !if(trim(Solver).eq."PTMD") then
  !  MPI_DETADT = eij2eta(de0ijdt)
  !  MPI_DFIELDDT_1D = AllFields1D(dFieldsdt)
  !end if

  if(mod(istep,TapeRate).eq.0) then
    do i = 1, NumReplicas
      if(i == ireplica) then
        write(*,'(I5,A1,I10,2F16.4,A3,F16.4,A12,F10.6)') i, "@", istep, Tk, " / ", TempList(i), "(K)  Epot: ", Etot/cgrid%npts
      end if
      call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    end do
  end if
  
  call MPI_GATHER(Etot,1,MPI_DOUBLE_PRECISION,&
  EtotRecv,1,MPI_DOUBLE_PRECISION,0,&
  MPI_COMM_WORLD,IERROR)
  call MPI_ALLGATHER(MPI_ETA,6,MPI_DOUBLE_PRECISION,&
  MPI_ETA_Recv,6,MPI_DOUBLE_PRECISION,&
  MPI_COMM_WORLD,IERROR)
  call MPI_ALLGATHER(MPI_FIELD_1D,DimFields1D,MPI_DOUBLE_PRECISION,&
  MPI_FIELD_1D_Recv,DimFields1D,MPI_DOUBLE_PRECISION,&
  MPI_COMM_WORLD,IERROR)
  !if(trim(Solver).eq."PTMD") then
  !  call MPI_ALLGATHER(MPI_DETADT,6,MPI_DOUBLE_PRECISION,&
  !  MPI_DETADT_Recv,6,MPI_DOUBLE_PRECISION,&
  !  MPI_COMM_WORLD,IERROR)
  !  call MPI_ALLGATHER(MPI_DFIELDDT_1D,DimFields1D,MPI_DOUBLE_PRECISION,&
  !  MPI_DFIELDDT_1D_Recv,DimFields1D,MPI_DOUBLE_PRECISION,&
  !  MPI_COMM_WORLD,IERROR)
  !  call MPI_ALLGATHER(MPI_GM,0,MPI_DOUBLE_PRECISION,&
  !  MPI_GM_Recv,0,MPI_DOUBLE_PRECISION,&
  !  MPI_COMM_WORLD,IERROR)
  !end if
  call MPI_ALLGATHER(MPI_EWALDFIELD_1D,DimFields1D,MPI_DOUBLE_PRECISION,&
  MPI_EWALDFIELD_1D_Recv,DimFields1D,MPI_DOUBLE_PRECISION,&
  MPI_COMM_WORLD,IERROR)
  
  call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
  
  if(ireplica.eq.1) then
    call GetPTSwapMap(TempList,EtotRecv,NumReplicas,Replicas)
  end if
  
  call MPI_BCAST(Replicas(:,1), NumReplicas, MPI_INTEGER, 0, MPI_COMM_WORLD, IERROR)
  call MPI_BCAST(TempList, NumReplicas, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERROR)
  call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
  e0ij=eta2eij(MPI_ETA_Recv(1+Replicas(ireplica,1)*6-6:&
  Replicas(ireplica,1)*6))
  Fields=AllFieldsFrom1D(MPI_FIELD_1D_Recv(1+Replicas(ireplica,1)*DimFields1D-DimFields1D:&
  Replicas(ireplica,1)*DimFields1D))
  !if(trim(Solver).eq."PTMD") then
  !  de0ijdt=eta2eij(MPI_DETADT_Recv(1+Replicas(ireplica,1)*6-6:&
  !  Replicas(ireplica,1)*6))
  !  dFieldsdt=AllFieldsFrom1D(MPI_DFIELDDT_1D_Recv(1+Replicas(ireplica,1)*DimFields1D-DimFields1D:&
  !  Replicas(ireplica,1)*DimFields1D))
  !  gm=MPI_GM_Recv(Replicas(ireplica,1))
  !end if
  EwaldField=AllFieldsFrom1D(MPI_EWALDFIELD_1D_Recv(1+Replicas(ireplica,1)*DimFields1D-DimFields1D:&
  1+Replicas(ireplica,1)*DimFields1D))

  if(ireplica.eq.NumReplicas) then
    e0ij = (1-MutationRatio)*e0ij &
         + MutationRatio*eta2eij(MPI_ETA_Recv(1+Replicas(1,1)*6-6:Replicas(1,1)*6))
    Fields = (1-MutationRatio)*Fields &
           + MutationRatio*AllFieldsFrom1D(MPI_FIELD_1D_Recv(1+Replicas(1,1)*DimFields1D-DimFields1D:Replicas(1,1)*DimFields1D))
    Call GetEwaldField(Fields, EwaldField)
  end if

  !if(ireplica.eq.NumReplicas) then
  !  Call GetInitConfig(Fields, e0ij, 'random', int2str5(ireplica-1))
  !  Call GetEwaldField(Fields, EwaldField)
  !end if

End Subroutine PTSwap
