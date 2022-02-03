Module Parameters
  use Constants
  use common
  implicit none
  integer               :: NODE_ME, NCPU, IERROR, IONODE
  Integer, parameter    :: ifileno = 55                          !< File handle number for input fiiles
  Integer, parameter    :: ofileno = 66                          !< File handle number for output files
  Integer, parameter    :: stdin = 5
  Integer, parameter    :: stdout = 6
  Integer               :: NumField = 3, FieldDim = 3
  Real*8                :: epinf = 6.649D0
  Real*8, dimension(3)  :: FieldCharge = (/1.585541183242639, -0.5562034806932713, 0./)
  Integer, dimension(4) :: FieldDimList = (/3, 3, 3, 6/)
  Integer               :: OnSiteDim = 9
  Real*8                :: alat(3,3)
  Real*8, dimension (:,:,:,:,:), allocatable :: EwaldMat
  Real*8, dimension (:,:,:),     allocatable :: EwaldHessian
  
  Real*8                :: CoeffDisp(14)
  Real*8                :: CoeffJijsm(7)
  Real*8                :: CoeffJijhf(7)
  Real*8                :: CoeffJijsmhf(2) = 0.0
  Real*8                :: Coeffu(4)
  Real*8                :: CoeffDispu(9)
  Real*8                :: CoeffEpsDisp(9)
  Real*8                :: CoeffEpsu(4)
  Real*8                :: CoeffEps(4)
  Real*8                :: mass(4) 
  Real*8                :: NoseMass(4) 
  Real*8                :: EAmp(4) 
  Real*8                :: EPhi(4) 
  Real*8                :: GateField(4) 
  Real(dp)              :: DWq(3,3)

  !---------------------------------------------------------------------------------
  ! system
  !---------------------------------------------------------------------------------
  type(s_unitcell) :: unitcell
  Real(dp)         :: symprec

  !---------------------------------------------------------------------------------
  ! parallelization
  !---------------------------------------------------------------------------------
  integer        :: nproc_k, nproc_ob, nproc_rgrid

  !---------------------------------------------------------------------------------
  ! DFT
  !---------------------------------------------------------------------------------
  integer        :: iperiodic
  integer        :: layout_multipole


  !---------------------------------------------------------------------------------
  ! Geometry and composition
  !---------------------------------------------------------------------------------
  type(s_grid) :: cgrid
  type(s_grid) :: rgrid
  type(s_grid) :: kgrid
  Real(dp), allocatable  ::  kbz(:,:)

  !---------------------------------------------------------------------------------
  ! TB and Jij data
  !---------------------------------------------------------------------------------
  Integer                      :: SpinDim, NumWann, OrbMul
  Integer                      :: NumWannSites, ContourNPoints(3)
  Integer                      :: supercell(3), k_mesh(3), fft_grid(3), Jij_R(3)
  Integer,allocatable          :: NumRpts(:), WannSiteInd(:), WannBlockInd(:,:), WannDist(:,:,:)
  real(dp)                     :: ContourMin, ContourMax, ContourHeight
  Integer                      :: potim
  complex(dp), allocatable     :: WannFunc(:,:,:,:,:,:)
  complex(dp), allocatable     :: WannBasis(:,:,:,:,:)
  real(dp), allocatable        :: redcoord(:,:,:)       !< Coordinates for Heisenberg exchange couplings
  real(dp), allocatable        :: jc(:,:,:,:,:)         !< Exchange couplings
  Integer,  allocatable        :: SiteOrbInfo(:,:)
  Integer                      :: JijSites(2)
  complex(dp), allocatable     :: HRmn(:,:,:,:), ContourPath(:)
  integer, allocatable         :: Ti0(:,:,:)
  Real(dp),allocatable         :: Efermi(:)
  logical                      :: WriteG

  !---------------------------------------------------------------------------------
  ! Simulation paramters
  !---------------------------------------------------------------------------------
  character(len=40) :: RestartFields    !< restart file
  character(len=40) :: RestartVelocity  !< restart file
  character(len=10) :: NameSim          !< Name of simulation
  character(len=20) :: CoeffFile        !< Name of coefficient data file
  character(len=20) :: TrajectoryFile   !< Name of Config file
  character(len=1)  :: aunits           !< Atomic units to simulate model Hamiltonians (Y/N)

  !---------------------------------------------------------------------------------
  ! Solvers
  !---------------------------------------------------------------------------------
  character(len=20) :: Solver                 !< Model solver
  character(len=20) :: EfieldType             !< Electric field type
  integer           :: ThermoSteps            !< Thermo up
  integer           :: CoolingSteps           !< Thermo up
  integer           :: NumSteps               !< Number of Monte Carlo steps
  integer           :: TapeRate               !< Number of steps in measurement phase
  Logical           :: DipoleQ                !< dipole dipole interaction
  Logical           :: EfieldQ                !< dipole dipole interaction
  Logical           :: CLAMPQ(4)              !< strain frozen
  real*8            :: Temp                   !< Temperature
  real*8            :: Pressure               !< Pressure

  !---------------------------------------------------------------------------------
  ! Markov Chain Monte Carlo
  !---------------------------------------------------------------------------------
  Real*8  :: damp                   !< damp size, move scale
  Real*8  :: DampRatio              !< damp size, move scale
  Real*8  :: AcceptRatio            !< accept ratio
  integer :: seed                   !< random seed
  integer :: AvrgInterval           !< Sampling interval for averages
  integer :: BuffMcAvrg             !< Buffer size for averages

  !---------------------------------------------------------------------------------
  ! Wang Landau Monte Carlo
  !---------------------------------------------------------------------------------
  Real*8  :: wl_a                 !< reducing factor
  integer :: wl_rate              !< spiking rate

  !---------------------------------------------------------------------------------
  ! Parallel Tempering Monte Carlo
  !---------------------------------------------------------------------------------
  Real*8  :: ReplicaT0
  Real*8  :: ReplicaTN
  integer :: SwapRate               !< Number of sweeps for one swap

  !---------------------------------------------------------------------------------
  ! Molecular Dynamics
  !---------------------------------------------------------------------------------
  Real*8            :: deltaT
  character(len=20) :: ThermoState 
  Logical           :: ReservoirQ 
  Integer           :: ReservoirRate
  Real*8            :: ReservoirRatio
  Real*8            :: ReservoirTau
  Real*8            :: ReservoirMass

  !---------------------------------------------------------------------------------
  ! Minimization
  !---------------------------------------------------------------------------------
  Character(len=20) :: OptAlgo
  Real*8            :: opt_tol
  Integer           :: optdim
  Integer           :: opt_iter=0

End Module Parameters
