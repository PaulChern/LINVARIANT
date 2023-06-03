Module Inputs
  Use FileParser
  Use Parameters
  Use Constants
  Use LINVARIANT
  Use Aux

  implicit none
  
  logical :: sane_input = .true.
  
  Contains
  
  include "ReadCoefficients.f90"

  subroutine io_file_unit(funit)
    !===========================================
    !
    !! Returns an unused unit number
    !! so we can later open a file on that unit.
    !
    !===========================================
    implicit none
    integer,intent(inout) :: funit
    logical               :: file_open

    funit = funit - 1
    file_open = .true.
    do while (file_open)
      funit = funit + 1
      inquire(funit, OPENED=file_open)
    end do

  end subroutine io_file_unit
 
  Subroutine read_input
    implicit none
    integer             :: i_err
    logical             :: file_exists

    namelist/control/ &
      NameSim,        &
      potim,          &
      MutationRatio,  &
      RestartFields,  &
      RestartVelocity

    namelist/calculation/ &
      Solver

    namelist/system/ &
      Temp,          &
      Pressure,      &
      SpinDim

    namelist/grid/ &
      supercell,   &
      k_mesh,      &
      fft_grid,    &
      DeltaT

    namelist/Ewald/ &
      DipoleQ

    namelist/files/ &
      CoeffFile,    &
      MLFile

    namelist/constrain/ &
      ClampQ

    namelist/extfield/ &
      EfieldQ,         &
      TrainQ,          &
      TrainRate,       &
      EfieldType
  
    namelist/exchanges/ &
      WriteG,           &
      Jij_R

    namelist/mcmd/   &
      ReplicaT0,     &
      ReplicaTN,     &
      seed,          &
      CoolingSteps,  &
      NumSteps,      &
      ThermoSteps,   &
      ReservoirQ,    &
      ReservoirRate, &
      ReservoirRatio,&
      ReservoirTau,  &
      ReservoirMass, &
      TapeRate,      &
      SwapRate,      &
      damp,          &
      DampRatio,     &
      AcceptRatio

    INQUIRE(FILE="LINVARIANT.inp", EXIST=file_exists)
    if (file_exists) then
      open(ifileno, file='LINVARIANT.inp', status='old')

      read(ifileno, nml=control, iostat=i_err)
      rewind(ifileno)

      read(ifileno, nml=calculation, iostat=i_err)
      rewind(ifileno)

      read(ifileno, nml=system, iostat=i_err)
      rewind(ifileno)

      read(ifileno, nml=grid, iostat=i_err)
      rewind(ifileno)

      read(ifileno, nml=exchanges, iostat=i_err)
      rewind(ifileno)

      read(ifileno, nml=files, iostat=i_err)
      rewind(ifileno)

      read(ifileno, nml=ewald, iostat=i_err)
      rewind(ifileno)

      read(ifileno, nml=constrain, iostat=i_err)
      rewind(ifileno)

      read(ifileno, nml=extfield, iostat=i_err)
      rewind(ifileno)

      read(ifileno, nml=mcmd, iostat=i_err)
      rewind(ifileno)

      close(ifileno)

    end if

  End subroutine read_input
 
  Subroutine ReadParameters
    
    implicit none
    
    character(len=50)   :: keyword, cache 
    integer             :: rd_len,i_err,i,i_stat,i_errb, i_all, pos
    logical             :: comment
    logical             :: file_exists

    call set_input_defaults

    INQUIRE(FILE="LINVARIANT.inp", EXIST=file_exists)
    if (file_exists) then
      open(ifileno,file='LINVARIANT.inp')
      do
        10     continue
        ! Read file character for character until first whitespace

        call bytereader(keyword,rd_len,ifileno,i_errb)
        
        ! converting Capital letters
        call caps2small(keyword)
        ! check for comment markers (currently % and #)
        comment=(scan(trim(keyword),'%')==1).or.(scan(trim(keyword),'#')==1).or.&
        (scan(trim(keyword),'*')==1).or.(scan(trim(keyword),'=')==1.or.&
        (scan(trim(keyword),'!')==1))
        if (comment) then
          read(ifileno,*)
        else
          ! Parse keyword
          keyword=trim(keyword)
          select case(keyword)
            case('restartfields')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) RestartFields
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('restartvelocity')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) RestartVelocity
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('namesim')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) NameSim
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('solver')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) Solver
            TrajectoryFile = trim(Solver)
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('thermostate')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) ThermoState
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('numsteps')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) NumSteps
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('spindim')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) SpinDim
            OrbMul = 1
            if(SpinDim.eq.3) OrbMul = 2
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('thermosteps')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) ThermoSteps
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('coolingsteps')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) CoolingSteps
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('writeg')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) WriteG
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('mutationratio')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) MutationRatio
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('dipoleq')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) DipoleQ
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('efieldq')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) EfieldQ
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('trainq')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) TrainQ
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('efieldtype')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) EfieldType
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('reservoirq')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) ReservoirQ
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('reservoirrate')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) Reservoirrate
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('reservoirratio')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) Reservoirratio
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('reservoirtau')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) ReservoirTau
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('reservoirmass')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) ReservoirMass
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('clampq')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) clampq
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('frozenq')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) frozenq
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('deltat')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) DeltaT
            DeltaT = DeltaT/time_fs
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('replicat0')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) ReplicaT0
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('replicatn')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) ReplicaTN
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('temp')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) Temp
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('pressure')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) Pressure
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('taperate')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) TapeRate
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err             
          case('trainrate')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) TrainRate
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('swaprate')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) SwapRate
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('wl_rate')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) wl_rate
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('wl_a')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) wl_a
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('seed')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) seed
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('damp')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) damp
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('dampratio')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) DampRatio
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('acceptratio')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) AcceptRatio
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('coefffile')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) CoeffFile
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('mlfile')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) MLFile
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('supercell')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) cgrid%n1, cgrid%n2, cgrid%n3
            cgrid%npts = cgrid%n1*cgrid%n2*cgrid%n3
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('screening')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) screening
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('k_mesh')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) kgrid%n1, kgrid%n2, kgrid%n3
            kgrid%npts = kgrid%n1*kgrid%n2*kgrid%n3
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('fft_grid')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) rgrid%n1, rgrid%n2, rgrid%n3
            rgrid%npts = rgrid%n1*rgrid%n2*rgrid%n3
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('jij_r')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) Jij_R
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('potim')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) potim
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('optalgo')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) OptAlgo
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('opt_tol')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) opt_tol
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case('symprec')
            read(ifileno, '(A)', iostat=i_err) cache
            pos = scan(cache, '=')
            cache = trim(cache(pos+1:))
            read(cache, *, iostat=i_err) symprec
            if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          case default
            if(len(trim(keyword))>0) then
              read(ifileno,*)
            end if
            
          end select
        end if
        
        ! End of file
        if (i_errb==20) goto 20
        ! End of row
        if (i_errb==10) goto 10
      end do
      
      20  continue

      close(ifileno)

    end if
    
    return
  End Subroutine ReadParameters
  
  !---------------------------------------------------------------------------------
  ! SUBROUTINE: set_input_defaults
  !> Sets default values to input variables
  !---------------------------------------------------------------------------------
  Subroutine set_input_defaults
    Use FileParser
    
    implicit none
    
    real(dp) :: one=1.0_dp
    real(dp) :: zero=0.0_dp
    
    !Geometry and composition
    cgrid%n1          = 12
    cgrid%n2          = 12
    cgrid%n3          = 12
    cgrid%npts        = cgrid%n1*cgrid%n2*cgrid%n3

    rgrid%n1          = 50
    rgrid%n2          = 50
    rgrid%n3          = 50
    rgrid%npts        = rgrid%n1*rgrid%n2*rgrid%n3
    kgrid%n1          = 10
    kgrid%n2          = 10
    kgrid%n3          = 10
    kgrid%npts        = kgrid%n1*kgrid%n2*kgrid%n3
    DWq               = 0.0_dp
    symprec           = 1.0E-6_dp

    !TB
    NumWannSites      = 0
    Efermi            = 0.0_dp
    ContourMin        = -30.0_dp
    ContourMax        = 0.0_dp
    ContourHeight     = 0.01
    ContourNPoints    = [75,1500,75]
    Jij_R             = [1,0,0]
    potim             = 100
    WriteG            = .False.
    
    !Solvers
    Solver            = "MCMC"
    ThermoSteps       = 10000
    CoolingSteps      = 0
    NumSteps          = 10000
    TapeRate          = 1000
    TrainRate         = 10000000000000000
    MutationRatio     = 0.0
    DipoleQ           = .true.
    EfieldQ           = .false.
    TrainQ            = .false.
    EfieldType        = "dc"
    Screening         = 0.2
    CLAMPQ            = (/.false., .false., .false., .false., .false., .false./)
    FrozenQ           = (/.false., .false., .false./)
    call caps2small(Solver)
    
    !Molecular Dynamics
    DeltaT            = 1.0D-1
    ThermoState       = "Nose-Hoover"
    ReservoirQ        = (/.true., .true./)
    ReservoirRate     = 1
    ReservoirRatio    = 1.0D0
    ReservoirTau      = 1.0D0
    ReservoirMass     = 1.0D0
    Pressure          = 0.0D0
    Ephi              = (/200.0, 0.0, 0.0, 0.0/)
    call caps2small(ThermoState)
    
    !Markov Chain Monte Carlo
    seed              = 0
    damp              = 0.001
    dampRatio         = 1.5
    AcceptRatio       = 0.3
    AvrgInterval      = 100
    BuffMcAvrg        = 0
    BuffMcAvrg        = 0.001

    !Wang Landau Monte Carlo
    wl_a              = 0.5
    wl_rate           = 10
    
    !Parallel Tempering Monte Carlo
    ReplicaT0         = 0.0D0
    ReplicaTN         = 500.0D0
    SwapRate          = 1

    !Minimization
    OptAlgo           = "LBFGS"
    opt_tol           = 1.0D-5
    
    !Simulation
    RestartFields     = "random"
    RestartVelocity   = "random"
    NameSim           = "LINVARIANT"
    CoeffFile         = "Coefficients.dat"
    MLFile            = "ModelData.dat"
    TrajectoryFile    = trim(Solver)
    aunits            = "N"
    SpinDim           = 1
    OrbMul            = 1
  End Subroutine set_input_defaults
  
  Subroutine InitFromFile(filename, Fields, e0ij)
    Implicit None
    
    character(*), Intent(in) :: filename
    Integer                  :: FileHandle = 1111
    Integer                  :: ix, iy, iz, i, j, ifield
    Real(dp)                 :: qdotr(3), tmp(3)
    Real*8, Intent(inout)    :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8, Intent(inout)    :: e0ij(3,3)
    
    e0ij = 0.0D0
    
    open(FileHandle,file=filename,form='formatted',status='old')
    do ifield = 1, NumField
      Read(FileHandle, "(3E25.15)") ((DWq(j,i,ifield), j=1,3), i=1,FieldDim)
    end do
    Read(FileHandle, *) e0ij(1,1), e0ij(2,2), e0ij(3,3)
    Read(FileHandle, *) e0ij(2,3), e0ij(1,3), e0ij(1,2)
    do While (.True.)
!      Read(FileHandle, "(3I10)", END=999) ix, iy, iz
      Read(FileHandle, *, END=999) ix, iy, iz
      do ifield = 1, NumField
        do i = 1, FieldDim
          qdotr(i)=Real(Exp(2*pi*cmplx(0.0_dp,1.0_dp)*(DWq(1,i,ifield)*ix+DWq(2,i,ifield)*iy+DWq(3,i,ifield)*iz)))
        end do
!        Read(FileHandle, "(3E25.15)", END=999) tmp
        Read(FileHandle, *, END=999) tmp
        do i = 1, FieldDim
          Fields(i,ifield,ix,iy,iz) = FieldsBinary(i,ifield)*tmp(i)*tanh(100.0*qdotr(i))
        end do
      end do
    end do
    999  close(FileHandle)
    
  End Subroutine InitFromFile

  Subroutine open_input_file(unit,filename)

    integer, intent(inout)       :: unit
    character(len=*), intent(in) :: filename
    integer                      :: ios = 0

    call io_file_unit(unit)
    open(unit=unit, file=trim(filename), iostat=ios, status="old", action="read")
    if ( ios /= 0 ) then
      write(stdout,*) "Error opening file "//filename
      stop
    end if

  End subroutine open_input_file

  Subroutine ReadContour
    implicit none

    Real(dp)  ::  E_Re, E_Im
    Integer   ::  ioene=1111, NumEne
    Integer   ::  i, j

    Call open_input_file(ioene,"epath.dat")

    Read(ioene,*) NumEne
    Allocate(ContourPath(NumEne))
    Do i = 1, NumEne
      Read(ioene,*) E_Re, E_Im
      ContourPath(i) = cmplx(E_Re,E_Im)
    End do

  End subroutine ReadContour

  Function ReadXSF(filename) result(xcrysden)
    implicit none
    character(72),intent(in) :: filename
    real(dp)                 :: xcrysden(3*rgrid%n1,3*rgrid%n2,3*rgrid%n3)
    character(72)            :: cache
    integer                  :: xsf=2020, ios
    logical                  :: comment
    logical                  :: file_exists
    Integer                  :: nr1, nr2, nr3
    Integer                  :: iwann, is, i, ix, iy, iz

    call io_file_unit(xsf)
    INQUIRE(FILE=trim(filename), EXIST=file_exists)
    if (file_exists) then
      call open_input_file(xsf,trim(filename))
      do While (.True.)
        read(xsf,'(A)') cache
        comment = scan(adjustl(cache),'#')==1
        if(.not.comment) EXIT
      end do

      do i = 1, unitcell%nion+15
        Read(xsf,'(A)') cache
      end do

      Read(xsf,*) nr1, nr2, nr3
      If(nr1.neqv.3*rgrid%n1.and.nr2.neqv.3*rgrid%n2.and.nr3.neqv.3*rgrid%n3) then
        write(*,*) "rgrid doesn't match with xsf"
        Call abort
      End if

      do i = 1, 4
        Read(xsf,*) cache
      end do

      Read(xsf,*) (((xcrysden(ix,iy,iz),ix=1,nr1),iy=1,nr2),iz=1,nr3)
    end if

  End Function ReadXSF

  Include "GetInitConfig.f90"

End Module Inputs
