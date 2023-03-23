#include "symbol.inc"
Program main
  Use LINVARIANT
  Use Invariant
  Use GreenFunction
  Use Wannier
  Use TB
  Use Parameters
  Use Constants
  Use timer
  Use Inputs
  Use Outputs
  Use fft
  Use Force
  Use Energy
  Use optimization
  Use Ewald
  Use Aux
  Use mpi
  
  Implicit none
  
  Integer                  :: i, j, istep, ireplica, ix, iy, iz, i_wl, ifield
  Real*8, allocatable      :: Fields(:,:,:,:,:)
  Real*8, allocatable      :: dFieldsdt(:,:,:,:,:)
  Real*8, allocatable      :: EwaldField(:,:,:,:,:)
  Real*8, allocatable      :: udamp(:), wl_h(:,:), wl_s(:,:)
  Real*8, allocatable      :: gm(:)
  Real*8, allocatable      :: TempList(:)
  Integer, allocatable     :: Replicas(:,:)
  character(len=10)        :: FileIndex
  Real*8                   :: etadamp, CGAlpha, wl_f
  Real*4                   :: tstart, tend
  Real*8                   :: e0ij(3,3), de0ijdt(3,3), detadt(6), Tk(2), dene
  Real*8                   :: Etot, Epot, Ekin, TimeNow
  Real(dp)                 :: rtest(10)
  complex(dp)              :: ctest(20), ltest(20)
  logical                  :: file_exists
  Integer                  :: wl_NumSteps

  IONODE = 0
  NODE_ME = 0
#ifdef MPI  
  call MPI_INIT(IERROR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, NCPU, IERROR)
  call MPI_COMM_RANK(MPI_COMM_WORLD, NODE_ME, IERROR)
#endif

  FileIndex = int2str5(NODE_ME)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Read control parameters and Coefficients!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  Call ReadParameters 
  Call ReadCoefficients
  do_io  Call fmkdir(trim(Solver)//".out")

  If(trim(Solver).eq."TRules") then
    Call MatrixRep
    Write(*,*) "Transformation rules generated."
    Call EXIT
  else if(trim(Solver).eq."JGreen") then
    Call ReadContour
    Call RectangleContour
    Call ReadWannier90('wannier90_hr', SpinDim)
    if(NumWann.ne.(SpinDim-1)*sum(SiteOrbInfo(2,:))) then
      write(*,*) "info in LINVARIANT.in and wannier_hr.dat don't match!"
      call abort
    end if 
    Call GreenJij
    Write(*,*) "Jij calculated."
    Call EXIT
  else if(trim(Solver).eq."FGreen") then
!    Call ReadContour
!    Call SemiCircleContour
    Call RectangleContour
    Call ReadWannier90('wannier90_hr', SpinDim)
    Call LoadWannFunc

    if(NumWann.ne.OrbMul*sum(SiteOrbInfo(2,:))) then
      write(*,*) "info in LINVARIANT.in and wannier_hr.dat don't match!"
      call abort
    end if
    Call GreenFij
    Write(*,*) "Fij calculated."
    Call EXIT
  else if(trim(Solver).eq."DFT") then
    Call DFT
    Call EXIT
  else if(trim(Solver).eq."test") then
    write(*,*) "This is a test"
    rtest = (/(i, i=1,10)/)
    call fft1d_r2c(1,-1,rtest,10,ctest,20)
    write(*,'(10F10.6)') rtest
    write(*,'(40F10.6)') ctest
    call fft1d_c2c(-1,1,ctest,20,ltest,20)
    write(*,'(40F10.6)') ltest
    write(*,'(40F10.6)') CSHIFT(ltest,1,DIM=1)
    Call EXIT
  End If

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Read or Write Ewald Matrix !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(EwaldMat(3,3,cgrid%n1,cgrid%n2,cgrid%n3))

#ifdef MPI
  do_io call EwaldMatrix(neighbourcut)
#endif

  call EwaldMatrix(neighbourcut)
  if(trim(Solver).eq."Ewald") Call EXIT

#ifdef MPI
  Call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
#endif

  allocate(EwaldHessian(3*cgrid%npts, 3*cgrid%npts, NumField))
  Call GetHessianEwald(EwaldMat, EwaldHessian)
  
  !!!!!!!!!!!!!!!!!!!!!
  ! Initialize Fields !
  !!!!!!!!!!!!!!!!!!!!!
  allocate(Fields(FieldDim,NumField,cgrid%n1,cgrid%n2,cgrid%n3))
  Call GetInitConfig(Fields, e0ij, RestartFields, FileIndex)

  allocate(dFieldsdt(FieldDim,NumField,cgrid%n1,cgrid%n2,cgrid%n3))
  Call GetInitConfig(dFieldsdt, de0ijdt, RestartVelocity, FileIndex)
 
  allocate(EwaldField(3,NumField,cgrid%n1,cgrid%n2,cgrid%n3))
  Call GetEwaldField(Fields, EwaldField)

  Call RemoveGridDrifts(Fields)
  Call RemoveGridDrifts(dFieldsdt)

  gm = 0.0D0
  TimeNow = 0.0D0
  detadt = 1.0D0
  do i = 1, 6
    If(CLAMPQ(i)) detadt(i) = 0.0D0
  end do
  de0ijdt = eta2eij(detadt)*de0ijdt
  Call AdaptTk(dFieldsdt, de0ijdt)
   
  !!!!!!!!!!!!!!!!!!!!!
  ! Start Simulations !
  !!!!!!!!!!!!!!!!!!!!!

#ifdef MPI
  Call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
#endif

  Call get_walltime(tend)
  allocate(udamp(NumField))
  allocate(gm(NumField+1))
  allocate(TempList(NCPU))
  allocate(Replicas(NCPU,5))

  If(trim(Solver).eq."Minimization") then
    !!!!!!!!!!!!!!!!!
    ! Initialize CG !
    !!!!!!!!!!!!!!!!!
    optdim = FieldDim*NumField*cgrid%npts+6
    !!!!!!!!!!!
    ! CG LOOP !
    !!!!!!!!!!!
    open(ifileno, &
    file=trim(Solver)//'.out/'//'trajectory_binary-'//trim(FileIndex)//'.dat', &
    form='unformatted',status='unknown')
    Call Minimization(Fields, e0ij)
    do_io Call WriteBinary(ifileno, Fields, dFieldsdt, e0ij, de0ijdt)
    Call Solver_finalize(Fields, e0ij, dFieldsdt, de0ijdt)
  else if(trim(Solver).eq."MD") then
    !!!!!!!!!!!!!!!!!
    ! Initialize MD !
    !!!!!!!!!!!!!!!!!
    open(ifileno, &
    file=trim(Solver)//'.out/'//'trajectory_binary-'//trim(FileIndex)//'.dat', &
    form='unformatted',status='unknown')
    Call RemoveGridDrifts(dFieldsdt)
    !!!!!!!!!!!
    ! MD LOOP !
    !!!!!!!!!!!
    do istep = 1, NumSteps
      do i = 1, NumField
        If(istep.gt.ThermoSteps.and.FrozenQ(i)) dFieldsdt(:,i,:,:,:) = 0.0D0
      end do
      Call TikTok(Fields, dFieldsdt, e0ij, de0ijdt, Temp, gm, TimeNow)
      if((mod(istep,TapeRate).eq.0).and.(istep.gt.ThermoSteps)) then
        tstart = tend
        Call get_walltime(tend)
        do_io Call MD_log(Fields, e0ij, dFieldsdt, de0ijdt, istep, tend-tstart, gm)
        do_io Call WriteBinary(ifileno, Fields, dFieldsdt, e0ij, de0ijdt)
      end if
    end do
    close(ifileno)
    Call Solver_finalize(Fields, e0ij, dFieldsdt, de0ijdt)
    Call Spectrum(NODE_ME)
  else if(trim(Solver).eq."WLMC") then
    !!!!!!!!!!!!!!!!!!!
    ! Initialize WLMC !
    !!!!!!!!!!!!!!!!!!!
    open(ifileno,file=trim(Solver)//'.out/wlmc_hs-'//trim(FileIndex)//'.dat',form='formatted',status='unknown')
    do i = 1, NumField
      udamp(i) = damp
    end do
    etadamp = damp
    Allocate(wl_s(2,ContourNPoints(2)))
    Allocate(wl_h(2,ContourNPoints(2)))
    wl_s = 0.0D0
    wl_h = 0
    wl_f = ee
    do i = 1, ContourNPoints(2)
      wl_s(1,i) = ContourMin + (i-1)*(ContourMax-ContourMin)/(ContourNPoints(2)-1)
      wl_h(1,i) = ContourMin + (i-1)*(ContourMax-ContourMin)/(ContourNPoints(2)-1)
    end do
    !!!!!!!!!!!!!
    ! WLMC LOOP !
    !!!!!!!!!!!!!
    Call random_seed()
    wl_NumSteps = 1
    ThermoSteps = Ceiling(wl_rate*(-8*DLog(10.0_dp)/DLog(wl_a)+2-wl_rate))
    do istep = 1, ThermoSteps + NumSteps
      if(mod(istep,wl_rate).eq.0) then
        tstart = tend
        Call get_walltime(tend)
        do_io Call WLMC_log(Fields, e0ij, istep, ThermoSteps, tend-tstart, udamp, etadamp, wl_f)
        wl_f = wl_f**(1/wl_a**(wl_rate-1))
      end if
      do i_wl = 1, wl_NumSteps
        Call WLMCStep(i_wl, Fields, e0ij, udamp, etadamp, wl_f)
        if(istep.gt.ThermoSteps) then
          write(ifileno, '(A8,I10)') "WL step:", istep-ThermoSteps
          do i = 1, ContourNPoints(2)
            write(ifileno, "(3E16.8)") wl_s(1,i), wl_h(2,i), wl_s(2,i)
            write(ifileno, "(3E25.15)") ((Sum(Fields(j,ifield,:,:,:))/cgrid%npts,j=1,FieldDim),ifield=1,NumField)
          end do 
        end if
      end do
      wl_f = wl_f**wl_a
      wl_NumSteps = Ceiling(wl_NumSteps*1/Sqrt(DLog(wl_f)))
    end do
    close(ifileno)
  else if(trim(Solver).eq."MCMC") then
    !!!!!!!!!!!!!!!!!!!
    ! Initialize MCMC !
    !!!!!!!!!!!!!!!!!!!
    open(ifileno, &
    file=trim(Solver)//'.out/'//'trajectory_binary-'//trim(FileIndex)//'.dat', &
    form='unformatted',status='unknown')
    dene = 0.0D0
    do i = 1, NumField
      udamp(i) = damp
    end do
    etadamp = damp
    !!!!!!!!!!!!!
    ! MCMC LOOP !
    !!!!!!!!!!!!!
    Call random_seed()
    do istep = 1, NumSteps
      Call MCMCStep(istep, Fields, EwaldField, e0ij, Temp, udamp, etadamp, dene)
      if((mod(istep,TapeRate).eq.0).and.(istep.gt.ThermoSteps)) then
        tstart = tend
        Call get_walltime(tend)
        do_io Call MCMC_log(Fields, e0ij, istep, tend-tstart, dene, udamp, etadamp)
        do_io Call WriteBinary(ifileno, Fields, dFieldsdt, e0ij, de0ijdt)
      end if
    end do
    close(ifileno)
    Call Solver_finalize(Fields, e0ij, dFieldsdt, de0ijdt)
  else if((trim(Solver).eq."PTMC").or.(trim(Solver).eq."PTMD")) then
    !!!!!!!!!!!!!!!!!
    ! Initialize PT !
    !!!!!!!!!!!!!!!!!
    io_begin
    INQUIRE(FILE='REPLICAS.dat', EXIST=file_exists)
    if (file_exists) then
      open(ifileno, file='REPLICAS.dat', status='old')
      do ireplica = 1, NCPU
        read(ifileno, *) TempList(ireplica)
      end do
      close(ifileno)
    else
      do ireplica = 1, NCPU
        TempList(ireplica) = ReplicaT0 + (ireplica-1)*(ReplicaTN-ReplicaT0)/(NCPU-1)
      end do
    end if
    io_end

#ifdef MPI
    call MPI_BCAST(TempList, NCPU, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, IERROR)
    call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
#endif

    do ireplica = 1, NCPU
      Replicas(ireplica,:) = (/ireplica, ireplica, 0, 0, 0/)
    end do

    dFieldsdt = Sqrt(TempList(NODE_ME+1)/Temp)*dFieldsdt
    de0ijdt = Sqrt(TempList(NODE_ME+1)/Temp)*de0ijdt

    do_io open(11,file=trim(Solver)//'.out/'//'REPLICAS.dat',status='unknown')
    do_io write(11, "("//trim(int2str5(NCPU))//"E15.6)") (TempList(i), i=1,NCPU)

    open(ifileno+100*NODE_ME,&
    file=trim(Solver)//'.out/'//'trajectory_binary-'//trim(FileIndex)//'.dat',&
    form='unformatted',status='unknown')
    Replicas(1,3) = 1
    Replicas(NCPU,3) = -1
    dene = 0.0D0
    do i = 1, NumField
      udamp(i) = damp
    end do
    etadamp = damp
    !!!!!!!!!!!
    ! PT LOOP !
    !!!!!!!!!!!
    Call random_seed()
    io_begin
      write(*,*) trim(Solver)//" solver runing in ", NCPU, " mpi processors."
      write(*,*) "Temperauture set: ", TempList
    io_end
    
    if(CoolingSteps.gt.0) then
      do ireplica = NCPU, 1, -1
        if(trim(Solver).eq."PTMC") then
          do istep = 1, CoolingSteps
            Call MCMCStep(istep, Fields, EwaldField, e0ij, &
            TempList(NODE_ME+1), udamp, etadamp, dene)
          end do
        else ! PTMD
          do istep = 1, CoolingSteps
            Call TikTok(Fields, dFieldsdt, e0ij, de0ijdt, TempList(NODE_ME+1), gm, TimeNow)
          end do
        end if

        if(ireplica.eq.NODE_ME+1) then
          Call WriteFinal('FinalConfig-'//trim(FileIndex)//'.dat', Fields, e0ij)
          Call WriteFinal('FinalVelocity-'//trim(FileIndex)//'.dat', dFieldsdt, de0ijdt)
        end if
        Call PTSwap(TempList, Fields, EwaldField, e0ij, dFieldsdt, de0ijdt, gm, NCPU, Replicas, NODE_ME+1, ireplica)
#ifdef MPI
        Call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
#endif
      end do
      RestartFields = 'FinalConfig'
      Call GetInitConfig(Fields, e0ij, RestartFields, FileIndex)
    end if

    CoolingSteps = 0
    gm = 0.0D0

    if(trim(Solver).eq."PTMC") then 
      do istep = 1, NumSteps
        Call MCMCStep(istep, Fields, EwaldField, e0ij, &
                      TempList(NODE_ME+1), udamp, etadamp, dene)
        if((mod(istep,TapeRate).eq.0).and.(istep.gt.ThermoSteps)) then
          io_begin
            tstart = tend
            Call get_walltime(tend)
            write(*,*) trim(Solver)//" step: ", istep, tend - tstart, "(second)"
          io_end
          Call WriteBinary(ifileno+100*NODE_ME, Fields, dFieldsdt, e0ij, de0ijdt)
        end if
        if((mod(istep,SwapRate).eq.0).and.(istep.lt.ThermoSteps)) then
          Call PTSwap(TempList, Fields, EwaldField, e0ij, dFieldsdt, de0ijdt, &
                      gm, NCPU, Replicas, NODE_ME+1, istep)
          io_begin 
            write(11, "("//trim(int2str5(NCPU))//"E15.6)") (TempList(i), i=1,NCPU)
            write(11, "("//trim(int2str5(NCPU))//"I8)") (Replicas(i,1), i=1,NCPU)
            write(11, "("//trim(int2str5(NCPU))//"I8)") (Replicas(i,2), i=1,NCPU)
            write(11, "("//trim(int2str5(NCPU))//"I8)") (Replicas(i,3), i=1,NCPU)
            write(11, "("//trim(int2str5(NCPU))//"I8)") (Replicas(i,4), i=1,NCPU)
            write(11, "("//trim(int2str5(NCPU))//"I8)") (Replicas(i,5), i=1,NCPU)
          io_end
        end if

#ifdef MPI
        Call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
#endif
      end do
    else ! PTMD
      do istep = 1, NumSteps
        Call TikTok(Fields, dFieldsdt, e0ij, de0ijdt, TempList(NODE_ME+1), gm, TimeNow)
        if((mod(istep,TapeRate).eq.0).and.(istep.gt.ThermoSteps)) then
          io_begin
            tstart = tend
            Call get_walltime(tend)
            write(*,*) trim(Solver)//" step: ", istep, tend - tstart, "(second)"
          io_end
          Epot = GetEtot(Fields, e0ij)
          Tk = Thermometer(dFieldsdt, de0ijdt)
          do i = 1, NCPU
            if(i == NODE_ME+1) then
              write(*,'(I5,A1,I10,2F10.4,A1,F10.4,A12,F10.6)') i, "@", istep, Tk, "/", &
                    TempList(i), "(K)  Epot: ", Epot
            end if

#ifdef MPI
            call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
#endif

          end do
          Call WriteBinary(ifileno+100*NODE_ME, Fields, dFieldsdt, e0ij, de0ijdt)
        end if
        if((mod(istep,SwapRate).eq.0).and.(istep.lt.ThermoSteps)) then
          Call PTSwap(TempList, Fields, EwaldField, e0ij, dFieldsdt, de0ijdt, gm, NCPU, Replicas, NODE_ME+1, istep)
          io_begin
            write(11, "("//trim(int2str5(NCPU))//"E15.6)") (TempList(i), i=1,NCPU)
            write(11, "("//trim(int2str5(NCPU))//"I8)") (Replicas(i,1), i=1,NCPU)
            write(11, "("//trim(int2str5(NCPU))//"I8)") (Replicas(i,2), i=1,NCPU)
            write(11, "("//trim(int2str5(NCPU))//"I8)") (Replicas(i,3), i=1,NCPU)
            write(11, "("//trim(int2str5(NCPU))//"I8)") (Replicas(i,4), i=1,NCPU)
            write(11, "("//trim(int2str5(NCPU))//"I8)") (Replicas(i,5), i=1,NCPU)
          io_end
        end if

#ifdef MPI
        Call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
#endif

      end do
    end if

    write(*,*) "PT process: ", NODE_ME, "of", NCPU, "done!"
    close(ifileno+100*NODE_ME)
    Call Solver_finalize(Fields, e0ij, dFieldsdt, de0ijdt)
    do_io  close(11)

#ifdef MPI
    call MPI_FINALIZE(IERROR)
#endif

  else
    do_io write(*,*) trim(Solver)//" is not implemented yet!"
    call abort
  end if
  
  !!!!!!!
  ! End !
  !!!!!!!
  deallocate(udamp)
  deallocate(TempList)
  deallocate(Replicas)
  deallocate(Fields)
  deallocate(EwaldMat)
  deallocate(EwaldHessian)
  deallocate(EwaldField)
  
  End Program main
