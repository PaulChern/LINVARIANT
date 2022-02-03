Module KohnSham
  Implicit none
  Contains

  subroutine init_dft(comm,info,lg,mg,system,stencil,fg,poisson,srg,srg_scalar,ofile)
    use common
    use parameters
    use sendrecv_grid
    use init_communicator
    use init_poisson_sub
    use checkpoint_restart_sub, only: init_dir_out_restart
    use sym_rho_sub, only: init_sym_rho
    implicit none
    integer, intent(in)      :: comm
    type(s_parallel_info)    :: info
    type(s_rgrid)            :: lg,mg
    type(s_dft_system)       :: system
    type(s_stencil)          :: stencil
    type(s_reciprocal_grid)  :: fg
    type(s_poisson)          :: poisson
    type(s_sendrecv_grid)    :: srg,srg_scalar
    type(s_ofile)            :: ofile
    !
    integer,dimension(2,3)   :: neig
  
  ! electron system
    call init_dft_system(lg,system,stencil)
  
  ! process distribution
    info%npk       = nproc_k
    info%nporbital = nproc_ob
    info%nprgrid   = nproc_rgrid
    call init_process_distribution(system,comm,info)
    call init_communicator_dft(comm,info)
  
  ! parallelization
    call check_ffte_condition(info,lg)
    call init_grid_parallel(info,lg,mg) ! lg --> mg
    call init_parallel_dft(system,info)
    call create_sendrecv_neig(neig, info) ! neighboring node array
    ! sendrecv_grid object for wavefunction updates
    call init_sendrecv_grid(srg, mg, info%numo*info%numk*system%nspin, info%icomm_rko, neig)
    ! sendrecv_grid object for scalar field updates
    call init_sendrecv_grid(srg_scalar, mg, 1, info%icomm_rko, neig)
  
  ! symmetry
  
    call init_sym_rho( lg%num, mg%is, mg%ie, info%icomm_r )
  
  ! for Poisson equation
  
    poisson%iterVh = 0 ! Iteration counter
    select case(iperiodic)
    case(0)
      if(layout_multipole==2.or.layout_multipole==3) call make_corr_pole(lg,mg,system,poisson)
    case(3)
      call init_reciprocal_grid(lg,mg,fg,system,info,poisson)
    end select
    call set_ig_bound(lg,mg,poisson)
  
  ! output files
  
    call init_dir_out_restart(ofile)
  
  end subroutine init_dft
  
  !===================================================================================================================================
  
  subroutine init_dft_system(lg,system,stencil)
    use structures
    use lattice
    use salmon_global, only: al_vec1,al_vec2,al_vec3,al,spin,natom,nelem,nstate,iperiodic,num_kgrid,num_rgrid,dl, &
                             nproc_rgrid,Rion,Rion_red,nelec,calc_mode,temperature,projection_option,nelec_spin, &
                             iflag_atom_coor,ntype_atom_coor_reduced,epdir_re1,quiet
    use sym_sub, only: init_sym_sub
    use communication, only: comm_is_root
    use parallelization, only: nproc_id_global
    use occupation_so, only: SPIN_ORBIT_ON, init_occupation_so
    implicit none
    type(s_rgrid)      :: lg
    type(s_dft_system) :: system
    type(s_stencil)    :: stencil
    !
    integer :: ii,jj
    real(8) :: rsize(3),hgs(3),cnmat(0:12,12),bnmat(4,4)
  
    if(al_vec1(2)==0d0 .and. al_vec1(3)==0d0 .and. al_vec2(1)==0d0 .and. &
       al_vec2(3)==0d0 .and. al_vec3(1)==0d0 .and. al_vec3(2)==0d0) then
      stencil%if_orthogonal = .true.
      system%primitive_a = 0d0
      system%primitive_a(1,1) = al(1)
      system%primitive_a(2,2) = al(2)
      system%primitive_a(3,3) = al(3)
      rsize = al
    else
      stencil%if_orthogonal = .false.
      system%primitive_a(1:3,1) = al_vec1
      system%primitive_a(1:3,2) = al_vec2
      system%primitive_a(1:3,3) = al_vec3
      rsize(1) = sqrt(sum(al_vec1**2))
      rsize(2) = sqrt(sum(al_vec2**2))
      rsize(3) = sqrt(sum(al_vec3**2))
    end if
  
    if(sum(abs(dl)) == 0d0) then
      hgs(1:3) = rsize(1:3) / dble(num_rgrid(1:3))
    else
      hgs(1:3) = dl(1:3)
    end if
    call init_grid_whole(rsize,hgs,lg)
    system%hgs = hgs
    system%ngrid = lg%num(1) * lg%num(2) * lg%num(3)
  
    call init_lattice(system,stencil)
    call init_sym_sub( system%primitive_a, system%primitive_b, epdir_re1 )
    call init_kvector(num_kgrid,system)
  
    if(calc_mode=='RT') then
      system%if_real_orbital = .false.
    else if(calc_mode=='GS') then
      select case(iperiodic)
      case(0)
        system%if_real_orbital = .true.
      case(3)
        if(num_kgrid(1)*num_kgrid(2)*num_kgrid(3)==1 .and. stencil%if_orthogonal) then
          system%if_real_orbital = .true.
        else
          system%if_real_orbital = .false.
        end if
        if ( SPIN_ORBIT_ON ) system%if_real_orbital=.false.
      end select
    end if
    if ((.not. quiet) .and. comm_is_root(nproc_id_global)) then
       write(*,*) "  use of real value orbitals = ", system%if_real_orbital
    endif
  
    system%nion = natom
  
    if(spin=='unpolarized') then
      system%nspin=1
    else if(spin=='polarized') then
      system%nspin=2
    end if
  
    if(calc_mode=='RT'.and. temperature<-1.d-12 .and. projection_option=='no')then
      if( SPIN_ORBIT_ON )then
         system%no = nelec
      else if(system%nspin==2.and.sum(nelec_spin(:))>0)then
        system%no = maxval(nelec_spin(:))
      else
        if(mod(nelec,2)==0)then
          system%no = nelec/2
        else
          system%no = (nelec+1)/2
        end if
      end if
    else
      system%no = nstate
    end if
  
    allocate(system%mass(1:nelem))
    if ( allocated(system%Rion) ) deallocate(system%Rion)
    if ( allocated(system%rocc) ) deallocate(system%rocc)
    if ( allocated(system%Velocity) ) deallocate(system%Velocity)
    if ( allocated(system%Force) ) deallocate(system%Force)
    allocate(system%Rion(3,system%nion),system%rocc(system%no,system%nk,system%nspin))
    allocate(system%Velocity(3,system%nion),system%Force(3,system%nion))
    system%Velocity(:,:) =0d0
  
    if(iflag_atom_coor==ntype_atom_coor_reduced) then
      Rion = matmul(system%primitive_a,Rion_red) ! [ a1, a2, a3 ] * R_ion
    end if
    system%Rion = Rion
  
  ! initial value of occupation
    system%rocc = 0d0
    if ( SPIN_ORBIT_ON ) then
      call init_occupation_so( system%rocc, nelec )
    else
      select case(system%nspin)
      case(1)
        system%rocc(1:nelec/2,:,1) = 2d0
      case(2)
        if ( nelec > 0 ) then
          if ( mod(nelec,2) == 0 ) then
            system%rocc(1:nelec/2,:,1:2) = 1d0
          else
            system%rocc(1:(nelec-1)/2,:,1:2) = 1d0
            system%rocc((nelec-1)/2+1,:,1  ) = 1d0
          end if
        else if ( any(nelec_spin>0) ) then
          system%rocc(1:nelec_spin(1),:,1) = 1d0
          system%rocc(1:nelec_spin(2),:,2) = 1d0
        else
          write(*,*) "nelect or nelec_spin should be specified in input"
        end if
      end select
    end if
  
    call set_bn(bnmat)
    call set_cn(cnmat)
    if(stencil%if_orthogonal) then
      stencil%coef_lap0 = -0.5d0*cNmat(0,Nd)*(1.d0/Hgs(1)**2+1.d0/Hgs(2)**2+1.d0/Hgs(3)**2)
    else
      if(nproc_rgrid(1)*nproc_rgrid(2)*nproc_rgrid(3)/=1) &
        stop "error: nonorthogonal lattice and r-space parallelization"
      stencil%coef_lap0 = -0.5d0*cNmat(0,Nd)*  &
                        & ( stencil%coef_F(1)/Hgs(1)**2 + stencil%coef_F(2)/Hgs(2)**2 + stencil%coef_F(3)/Hgs(3)**2 )
    end if
    do jj=1,3
      do ii=1,4
        stencil%coef_lap(ii,jj) = cnmat(ii,4)/hgs(jj)**2
        stencil%coef_nab(ii,jj) = bnmat(ii,4)/hgs(jj)
      end do
    end do
  
    call set_gridcoordinate(lg,system)
  
    system%vec_Ac = 0d0     ! initial value
  
    system%vec_Ac_ext = 0d0 ! initial value for output
    system%vec_E = 0d0      ! initial value for output
    system%vec_E_ext = 0d0  ! initial value for output
  
    return
  end subroutine init_dft_system

End module 
