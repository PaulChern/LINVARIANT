module common
  use xc_f90_lib_m
  implicit none

  type s_grid
    integer :: n1, n2, n3
    Integer :: npts
    Integer :: ncells
  end type s_grid

  type s_rgrid
    integer              :: ndir,Nd                 ! ndir=3 --> dir=xx,yy,zz, ndir=6 --> dir=xx,yy,zz,yz,zx,xy
    integer,dimension(3) :: is,ie,num &             ! num=ie-is+1
                           ,is_overlap,ie_overlap & ! is_overlap=is-Nd, ie_overlap=ie+Nd
                           ,is_array,ie_array       ! allocate(array(is_array(1):ie_array(1), ...) )
    integer ,allocatable :: idx(:),idy(:),idz(:)    ! idx(is_overlap(1):ie_overlap(1))=is_array(1)~ie_array(1), ...
    integer ,allocatable :: is_all(:,:),ie_all(:,:) ! (1:3,0:nproc-1), is & ie for all MPI processes
    real(8) ,allocatable :: coordinate(:,:)         ! (minval(is_overlap):maxval(ie_overlap),1:3), coordinate of grids
  end type s_rgrid

  type s_reciprocal_grid
    logical,allocatable    :: if_Gzero(:,:,:)
    real(8),allocatable    :: vec_G(:,:,:,:)   ! G vector (reciprocal lattice vector)
    real(8),allocatable    :: coef(:,:,:)      ! 4*pi/|G|^2 (coefficient of the Poisson equation)
    real(8),allocatable    :: exp_ewald(:,:,:) ! exp(-|G|^2/(4*a_Ewald))
    complex(8),allocatable :: egx(:,:),egxc(:,:),egy(:,:),egyc(:,:),egz(:,:),egzc(:,:)
    complex(8),allocatable :: coef_nabla(:,:,:,:),coef_gxgy0(:,:,:),cos_cGdt(:,:,:),sin_cGdt(:,:,:) ! for single-scale Maxwell-TDDFT
  end type s_reciprocal_grid

  ! scalar field
  type s_scalar
    real(8),allocatable :: f(:,:,:) ! f(x,y,z)
  end type s_scalar

  ! vector field
  type s_vector
    real(8),allocatable :: v(:,:,:,:) ! v(1:3,x,y,z)
  end type s_vector

  type s_unitcell
    logical :: if_real_orbital
    integer :: ngrid,nspin,nband,nkpt,nion ! # of r-grid points, spin indices, orbitals, k points, and ions
    real(8) :: hvol,hgs(3),latt_a(3,3),vol,latt_b(3,3)
    real(8) :: rmatrix_a(3,3),rmatrix_b(3,3)
    real(8) :: mu
    real(8),allocatable :: vec_k(:,:)    ! (1:3,1:nk), k-vector
    real(8),allocatable :: wtk(:)        ! (1:nk), weight of k points
    real(8),allocatable :: rocc(:,:,:)   ! (1:no,1:nk,1:nspin), occupation rate
  ! atomic ...
    real(8),allocatable :: Mass(:)       ! (1:nelem), Atomic weight
    real(8),allocatable :: site(:,:)     ! (1:3,1:nion), atom position
    real(8),allocatable :: Velocity(:,:) ! (1:3,1:nion), atomic velocity
    real(8),allocatable :: Force(:,:)    ! (1:3,1:nion), force on atom
  ! external field
    real(8)        :: vec_Ac(3)       ! A/c (spatially averaged), A: vector potential, c: speed of light
    type(s_vector) :: Ac_micro        ! A/c (microscopic)      ! for single-scale Maxwell-TDDFT
    type(s_scalar) :: div_Ac          ! divergence of Ac_micro ! for single-scale Maxwell-TDDFT
    real(8)        :: vec_Ac_ext(3)   ! external vector potential for output
    real(8)        :: vec_E(3)        ! total electric field for output
    real(8)        :: vec_E_ext(3)    ! external electric potential for output
  end type s_unitcell

! Poisson equation
  type s_poisson
  ! for poisson_cg (conjugate-gradient method)
    integer :: iterVh                              ! iteration number for poisson_cg
    integer :: npole_partial                       ! number of multipoles calculated in each node
    integer :: npole_total                         ! total number of multipoles
    integer,allocatable :: ipole_tbl(:)            ! table for multipoles
    integer,allocatable :: ig_num(:)               ! number of grids for domains to which each multipole belongs
    integer,allocatable :: ig(:,:,:)               ! grid table for domains to which each multipole belongs
    integer,allocatable :: ig_bound(:,:,:)         ! grid table for boundaries
    real(8),allocatable :: wkbound(:), wkbound2(:) ! values on boundary represented in one-dimentional grid
    integer :: n_multipole_xyz(3)                  ! number of multipoles
  ! for Fourier transform
    complex(8),allocatable :: zrhoG_ele(:,:,:)     ! rho_ele(G): Fourier transform of the electron density
  ! for discrete Fourier transform (general)
    complex(8),allocatable :: ff1x(:,:,:),ff1y(:,:,:),ff1z(:,:,:),ff2x(:,:,:),ff2y(:,:,:),ff2z(:,:,:)
  ! for FFTE
    complex(8),allocatable :: a_ffte(:,:,:),b_ffte(:,:,:)
  end type s_poisson

! exchange-correlation functional
  type s_xc_functional
    integer :: xctype(3)
    integer :: ispin
    real(8) :: cval
    logical :: use_gradient
    logical :: use_laplacian
    logical :: use_kinetic_energy
    logical :: use_current
    TYPE(xc_f90_func_t)      :: func(3)
    TYPE(xc_f90_func_info_t) :: info(3)
  end type

! for persistent communication
  type s_pcomm_cache
    real(8), allocatable :: dbuf(:, :, :, :)
    complex(8), allocatable :: zbuf(:, :, :, :)
  end type s_pcomm_cache

  type s_parallel_info
  ! division of MPI processes (for orbital wavefunction)
    integer :: npk        ! k-points
    integer :: nporbital  ! orbital index
    integer :: nprgrid(3) ! r-space (x,y,z)
  ! parallelization of orbital wavefunction
    integer :: iaddress(5) ! address of MPI under orbital wavefunction(ix,iy,iz,io,ik)
    integer,allocatable :: imap(:,:,:,:,:) ! address map
    logical :: if_divide_rspace
    logical :: if_divide_orbit
    integer :: icomm_r,   id_r,   isize_r   ! communicator, process ID, & # of processes for r-space
    integer :: icomm_k,   id_k,   isize_k   ! communicator, process ID, & # of processes for k-space
    integer :: icomm_o,   id_o,   isize_o   ! communicator, process ID, & # of processes for orbital
    integer :: icomm_ro,  id_ro,  isize_ro  ! communicator, process ID, & # of processes for r-space & orbital
    integer :: icomm_ko,  id_ko,  isize_ko  ! communicator, process ID, & # of processes for k-space & orbital
    integer :: icomm_rko, id_rko, isize_rko ! communicator, process ID, & # of processes for r-space, k-space & orbital
    integer :: im_s,im_e,numm ! im=im_s,...,im_e, numm=im_e-im_s+1
    integer :: ik_s,ik_e,numk ! ik=ik_s,...,ik_e, numk=ik_e-ik_s+1
    integer :: io_s,io_e,numo ! io=io_s,...,io_e, numo=io_e-io_s+1
                              ! For calc_mode='RT' and temperature<0, these values are calculated from nelec.
                              ! In other cases, these are calculated from nstate.
  ! sub-communicators of icomm_r (r-space)
    integer :: icomm_x,id_x,isize_x ! x-axis
    integer :: icomm_y,id_y,isize_y ! y-axis
    integer :: icomm_z,id_z,isize_z ! z-axis
    integer :: icomm_xy,id_xy,isize_xy ! for singlescale FDTD
  ! for atom index #ia
    integer :: ia_s,ia_e ! ia=ia_s,...,ia_e
    integer :: nion_mg
    integer,allocatable :: ia_mg(:)
  ! for orbital index #io
    integer,allocatable :: irank_io(:) ! MPI rank of the orbital index #io
    integer,allocatable :: io_s_all(:) ! io_s for all orbital ranks
    integer,allocatable :: io_e_all(:) ! io_e for all orbital ranks
    integer,allocatable :: numo_all(:) ! numo for all orbital ranks
    integer :: numo_max ! max value of numo_all
#ifdef USE_SCALAPACK
    logical :: flag_blacs_gridinit
    integer :: icomm_sl ! for summation
    integer :: iam,nprocs
    integer,allocatable :: gridmap(:,:)
    integer :: nprow,npcol,myrow,mycol
    integer :: nrow_local,ncol_local,lda
    integer :: desca(9), descz(9)
    integer :: len_work  ! for PDSYEVD, PZHEEVD
    integer :: len_rwork ! for PZHEEVD
    integer,allocatable :: ndiv(:), i_tbl(:,:), j_tbl(:,:), iloc_tbl(:,:), jloc_tbl(:,:)
#endif
#ifdef USE_EIGENEXA
    logical :: flag_eigenexa_init
#endif
  end type s_parallel_info

! for update_overlap
  type s_sendrecv_grid
    ! Number of orbitals (4-th dimension of grid)
    integer :: nb
    ! Communicator
    integer :: icomm
    ! Neightboring MPI id (1:upside,2:downside, 1:x,2:y,3:z):
    integer :: neig(1:2, 1:3)
    ! Communication requests (1:send,2:recv, 1:upside,2:downside, 1:x,2:y,3:z):
    integer :: ireq_real8(1:2, 1:2, 1:3)
    integer :: ireq_complex8(1:2, 1:2, 1:3)
    ! PComm cache (1:src/2:dst, 1:upside,2:downside, 1:x,2:y,3:z)
    type(s_pcomm_cache) :: cache(1:2, 1:2, 1:3)
    ! Range (axis=1...3, 1:src/2:dst, dir=1:upside,2:downside, dim=1:x,2:y,3:z)
    integer :: is_block(1:3, 1:2, 1:2, 1:3)
    integer :: ie_block(1:3, 1:2, 1:2, 1:3)
    ! Initialization flags
    logical :: if_pcomm_real8_initialized
    logical :: if_pcomm_complex8_initialized
  end type s_sendrecv_grid

  type s_stencil
    logical :: if_orthogonal
    real(8) :: coef_lap0,coef_lap(4,3),coef_nab(4,3) ! (4,3) --> (Nd,3) (future work)
    real(8) :: coef_f(6) ! for non-orthogonal lattice
  end type s_stencil

! output files
  type s_ofile
    integer :: fh_eigen
    integer :: fh_rt
    integer :: fh_rt_energy
    integer :: fh_response
    integer :: fh_pulse
    integer :: fh_dft_md
    integer :: fh_ovlp,fh_nex
    character(100) :: file_eigen_data
    character(256) :: file_rt_data
    character(256) :: file_rt_energy_data
    character(256) :: file_response_data
    character(256) :: file_pulse_data
    character(256) :: file_dft_md
    character(256) :: file_ovlp,file_nex
    !
    character(256) :: dir_out_restart, dir_out_checkpoint
  end type s_ofile


end module common
