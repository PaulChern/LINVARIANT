!-------------------------------------------------------------------------------
! MODULE: Constants
!> Physical constants
!> The magnitude of the parameters has been chosen so that it mtaches
!> the values given in the NIST database with the CODATA recommended 2014 values
!> https://physics.nist.gov/cuu/Constants/
!-------------------------------------------------------------------------------
module Constants
  implicit none
  Integer, parameter    :: snglprec = selected_real_kind(6, 37)  !< define precision for single reals
  Integer, parameter    :: dp = selected_real_kind(15, 307) !< define precision for double reals
  !.. Scalar parameters
  complex(dp),parameter :: ii = (0.0d0, 1.0d0)            
  real(dp),parameter    :: pi = 3.141592653589793_dp
  real(dp),parameter    :: ee = 2.718281828459045_dp
  real(dp) :: gama         = 1.760859644d11     ! s^(-1)*T^(-1)
  real(dp) :: k_bolt       = 1.38064852d-23     ! J/K
  real(dp) :: k_bolt_ev    = 8.6173303d-5       ! eV/K
  real(dp) :: mub          = 5.7883818012d-5    ! eV/T
  real(dp) :: mu0          = 1.2566370614d-6    ! N/A^2
  real(dp) :: mry          = 2.179872325d-21    ! J
  real(dp) :: hbar_mev     = 6.582119514d-13    ! meV*s
  real(dp) :: hbar_ev      = 6.582119514d-16    ! eV*s
  real(dp) :: Joule_ev     = 6.241509126d18     ! eV
  real(dp) :: Hartree      = 27.211386245       ! eV
  real(dp) :: ry_ev        = 13.605693009_dp     ! eV
  real(dp) :: hbar         = 1.054571800e-34    ! J*s
  real(dp) :: ev           = 1.6021766208d-19   ! J
  real(dp) :: amu          = 1.660539040d-27    ! kg
  real(dp) :: ame          = 9.1093837015d-31   ! kg
  real(dp) :: mpme         = 1.83615267343d3    ! mp/me ratio
  real(dp) :: angstrom     = 1.0d-10            ! m
  real(dp) :: bohr         = 0.52917721067d-10  ! m
  real(dp) :: s2fs         = 1.0d12             ! fs
  real(dp) :: time_fs      = 2.4188843265857d-2 ! fs
  real(dp) :: g_e_abs      = 2.00231930436182
  real(dp) :: hbar2_2me    = 3.8099822830657697 ! eV*Angstrom^2
end module Constants
