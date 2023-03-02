Subroutine Minimization(Fields, e0ij)
  Use Parameters
  Use Aux
  Use Energy
  Use Force
  Use Optimization

  Implicit none
  Real*8, Intent(inout)  :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
  Real*8, Intent(inout)  :: e0ij(3,3)
  Real(dp), allocatable  :: opt_lb(:), opt_ub(:), x(:)
  Real(dp)               :: opt_min
  integer*8              :: opt_obj
  integer                :: ires, istep, i
  include 'nlopt.f'

  Allocate(opt_lb(optdim))
  Allocate(opt_ub(optdim))
  Allocate(x(optdim))
  !!!!!!!!!!!
  ! CG LOOP !
  !!!!!!!!!!!

  if(trim(OptAlgo).eq."MMA") then
    call nlo_create(opt_obj, NLOPT_LD_MMA, optdim)
  else if(trim(OptAlgo).eq."CCSAQ") then
    call nlo_create(opt_obj, NLOPT_LD_CCSAQ, optdim)
  else if(trim(OptAlgo).eq."SLSQP") then
    call nlo_create(opt_obj, NLOPT_LD_SLSQP, optdim)
  else if(trim(OptAlgo).eq."PTNEWTON_RESTART") then
    call nlo_create(opt_obj, NLOPT_LD_TNEWTON_PRECOND_RESTART, optdim)
  else if(trim(OptAlgo).eq."PTNEWTON") then
    call nlo_create(opt_obj, NLOPT_LD_TNEWTON_PRECOND, optdim)
  else if(trim(OptAlgo).eq."TNEWTON_RESTART") then
    call nlo_create(opt_obj, NLOPT_LD_TNEWTON_RESTART, optdim)
  else if(trim(OptAlgo).eq."TNEWTON") then
    call nlo_create(opt_obj, NLOPT_LD_TNEWTON, optdim)
  else if(trim(OptAlgo).eq."SLBFGS2") then
    call nlo_create(opt_obj, NLOPT_LD_VAR2, optdim)
  else if(trim(OptAlgo).eq."SLBFGS1") then
    call nlo_create(opt_obj, NLOPT_LD_VAR1, optdim)
  else if(trim(OptAlgo).eq."LBFGS") then
    call nlo_create(opt_obj, NLOPT_LD_LBFGS, optdim)
!    call nlo_set_vector_storage(opt_obj, 10)
  else if(trim(OptAlgo).eq."COBYLA") then
    call nlo_create(opt_obj, NLOPT_LN_COBYLA, optdim)
  else if(trim(OptAlgo).eq."BOBYQA") then
    call nlo_create(opt_obj, NLOPT_LN_BOBYQA, optdim)
  else if(trim(OptAlgo).eq."NEWUOA") then
    call nlo_create(opt_obj, NLOPT_LN_NEWUOA, optdim)
  else if(trim(OptAlgo).eq."NEWUOA_BOUND") then
    call nlo_create(opt_obj, NLOPT_LN_NEWUOA_BOUND, optdim)
  else if(trim(OptAlgo).eq."PRAXIS") then
    call nlo_create(opt_obj, NLOPT_LN_PRAXIS, optdim)
  else if(trim(OptAlgo).eq."NELDERMEAD") then
    call nlo_create(opt_obj, NLOPT_LN_NELDERMEAD, optdim)
  else if(trim(OptAlgo).eq."Sbplx") then
    call nlo_create(opt_obj, NLOPT_LN_SBPLX, optdim)
  else if(trim(OptAlgo).eq."CRS") then
    call nlo_create(opt_obj, NLOPT_GN_CRS2_LM, optdim)
  else if(trim(OptAlgo).eq."StoGo") then
    call nlo_create(opt_obj, NLOPT_GD_STOGO, optdim)
  else if(trim(OptAlgo).eq."StoGoRand") then
    call nlo_create(opt_obj, NLOPT_GD_STOGO_RAND, optdim)
  else if(trim(OptAlgo).eq."AGS") then
    call nlo_create(opt_obj, NLOPT_GN_AGS, optdim)
  else if(trim(OptAlgo).eq."ISRES") then
    call nlo_create(opt_obj, NLOPT_GN_ISRES, optdim)
  else if(trim(OptAlgo).eq."ESCH") then
    call nlo_create(opt_obj, NLOPT_GN_ESCH, optdim)
  else
    call nlo_create(opt_obj, NLOPT_LD_LBFGS, optdim)
!    call nlo_set_vector_storage(opt_obj, 10)
  End if

  call nlo_get_lower_bounds(ires, opt_obj, opt_lb)
  call nlo_get_upper_bounds(ires, opt_obj, opt_ub)
  opt_lb = -10
  opt_ub = 10
  call nlo_set_lower_bounds(ires, opt_obj, opt_lb)
  call nlo_set_upper_bounds(ires, opt_obj, opt_ub)
  call nlo_set_min_objective(ires, opt_obj, OptEnergyFunc, 0)

  call nlo_remove_inequality_constraints(ires, opt_obj)
  call nlo_remove_equality_constraints(ires, opt_obj)

!  call nlo_add_inequality_constraint(ires, opt, myconstraint, d1, 1.D-8)
!  call nlo_add_inequality_constraint(ires, opt, myconstraint, d2, 1.D-8)

  call nlo_set_xtol_rel(ires, opt_obj, opt_tol)

  x = GetSysVector(Fields,e0ij)
  call nlo_optimize(ires, opt_obj, x, opt_min)
  if (ires.lt.0) then
     write(*,*) 'nlopt failed!'
  else
     write(*,*) 'found min at ', x(1), x(2)
     write(*,*) 'min val = ', opt_min
  endif

  call nlo_destroy(opt_obj)
  Call SysVector2Fields(x, Fields, e0ij)
  Deallocate(x)
  Deallocate(opt_lb)
  Deallocate(opt_ub)

End subroutine minimization
