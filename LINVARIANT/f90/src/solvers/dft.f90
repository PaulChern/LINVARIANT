Subroutine DFT
  Use Constants
  Use common
  Use lxc
  Use parallelization
!  Use KohnSham
  Implicit none

!  real(dp)              :: cval
!  type(s_xc_functional) :: xc_func
!  type(s_parallel_info) :: info
!  type(s_unitcell)      :: system
!  type(s_poisson)       :: poisson
!  type(s_stencil)       :: stencil
!  type(s_rgrid)         :: lg
!  type(s_rgrid)         :: mg
!  type(s_reciprocal_grid) :: fg
!  type(s_sendrecv_grid) :: srg, srg_scalar
!  type(s_ofile)         :: ofl
!  character(64)         :: spin
!  character(64)         :: xc
!  character(64)         :: xname
!  character(64)         :: cname
!
!  cval = 1.0d0
!  xc = 'pbe'
!  xname = 'none'
!  cname = 'none'
!  spin = 'unpolarized'

!  call init_xc(xc_func, spin, cval, xc, xname, cname)
!  call init_dft(nproc_group_global,info,lg,mg,system,stencil,fg,poisson,srg,srg_scalar,ofl)

End subroutine
