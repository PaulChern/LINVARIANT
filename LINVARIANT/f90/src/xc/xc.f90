Module lxc
  use Aux
  use common
  Contains

  subroutine lxc_test
    use xc_f90_lib_m
    implicit none
    TYPE(xc_f90_func_t) :: xc_func
    TYPE(xc_f90_func_info_t) :: xc_info
    integer(8), parameter :: npoints = 5
    real(8) :: rho(npoints) = (/0.1, 0.2, 0.3, 0.4, 0.5/)
    real(8) :: sigma(npoints) = (/0.2, 0.3, 0.4, 0.5, 0.6/)
    real(8) :: exc(npoints)
    integer :: i, vmajor, vminor, vmicro, func_id = 1
  
    ! Print out the version  
    call xc_f90_version(vmajor, vminor, vmicro)
    write(*,'("Libxc version: ",I1,".",I1,".",I1)') vmajor, vminor, vmicro
  
    ! Initialize the functional
    call xc_f90_func_init(xc_func, func_id, XC_UNPOLARIZED)
  
    ! Get the information
    xc_info = xc_f90_func_get_info(xc_func)
  
    ! Evaluate energy depending on the family
    select case (xc_f90_func_info_get_family(xc_info))
      case(XC_FAMILY_LDA)
        call xc_f90_lda_exc(xc_func, npoints, rho(1), exc(1))
      case(XC_FAMILY_GGA, XC_FAMILY_HYB_GGA)
        call xc_f90_gga_exc(xc_func, npoints, rho(1), sigma(1), exc(1))
    end select
  
    ! Print out density and energy density per particle
    do i = 1, npoints
      write(*,"(F8.6,1X,F9.6)") rho(i), exc(i)
    end do
  
    ! Deallocate memory
    call xc_f90_func_end(xc_func)
  
  end subroutine lxc_test

  subroutine init_xc(xc, spin, cval, xcname, xname, cname)
    implicit none
    type(s_xc_functional), intent(inout) :: xc
    character(*), intent(in)           :: spin
    real(8), intent(in)                :: cval
    character(*), intent(in), optional :: xcname
    character(*), intent(in), optional :: xname
    character(*), intent(in), optional :: cname

    ! Initialization of xc variable
    xc%xctype(1:3) = 0
    xc%cval = cval
    xc%use_gradient = .false.
    xc%use_laplacian = .false.
    xc%use_kinetic_energy = .false.
    xc%use_current = .false.

    if(spin=='unpolarized') then
      xc%ispin=0
    else if(spin=='polarized') then
      xc%ispin=1
    end if

    if( xcname=="none" .and. xname=="none" .and. cname=="none" ) then
       print '(A, A)', "Error! Exchange nad Correction functionals are not specified!"
       stop
    endif

    ! Exchange correlation
    if (present(xcname) .and. (len_trim(xcname) > 0)) then
      call setup_xcfunc(xcname)
    end if

    ! Exchange only
    if (present(xname) .and. (len_trim(xname) > 0)) then
      call setup_xfunc(xname)
    end if

    ! Correlation only
    if (present(cname) .and. (len_trim(cname) > 0)) then
      call setup_cfunc(cname)
    end if

    return
  contains

    subroutine setup_xcfunc(name)
      implicit none
      character(*), intent(in) :: name

      select case(lower(name))
      case('none')

        return
      
      case ('tbmbj')

        xc%xctype(1) = 4
        xc%use_gradient = .true.
        xc%use_laplacian = .true.
        xc%use_kinetic_energy = .true.
        xc%use_current = .true.
        return

      case ('bj_pw')

        xc%xctype(1) = 4 
        xc%cval = 1d0
        xc%use_gradient = .true.
        xc%use_laplacian = .true.
        xc%use_kinetic_energy = .true.
        xc%use_current = .true.
        return

      case('pz')
        xc%xctype(2) = 101
        xc%xctype(3) = 101
        call init_libxc('LDA_X', 2)
        call init_libxc('LDA_C_PZ', 3)
        return

      case('pzm')
        xc%xctype(2) = 101
        xc%xctype(3) = 101
        call init_libxc('LDA_X', 2)
        call init_libxc('LDA_C_PZ_MOD', 3)
        return

      case('pbe')
        xc%xctype(2) = 101
        xc%xctype(3) = 101
        call init_libxc('GGA_X_PBE', 2)
        call init_libxc('GGA_C_PBE', 3)
        return

      case default

        xc%xctype(1) = 101
        call init_libxc(name, 1)
        return
      end select
      return
    end subroutine



    subroutine setup_xfunc(name)
      implicit none
      character(*), intent(in) :: name


      select case(name)
      case('none')
        return

      case default

        xc%xctype(2) = 101
        call init_libxc(name, 2)
        return

      end select
      return
    end subroutine


    subroutine setup_cfunc(name)
      implicit none
      character(*), intent(in) :: name


      select case(name)
      case('none')
        return

      case default

        xc%xctype(3) = 101
        call init_libxc(name, 3)
        return

      end select
      return
    end subroutine



    subroutine init_libxc(libxc_name, ii)
      implicit none
      character(*), intent(in) :: libxc_name
      integer, intent(in) :: ii
      integer ::  ixc

      ixc = xc_f90_functional_get_number(trim(libxc_name))

      if (ixc <= 0) then
        print '(A, A)', "Undefined libxc functional:", trim(libxc_name)
        stop
      end if

      if (spin /= 'unpolarized') then
        print '(A)', "Spin polarized is not available"
        stop
      end if

      call xc_f90_func_init(xc%func(ii), ixc, XC_UNPOLARIZED)
      xc%info(ii) = xc_f90_func_get_info(xc%func(ii))

      select case (xc_f90_func_info_get_family(xc%info(ii)))
      case (XC_FAMILY_LDA)
      case (XC_FAMILY_GGA)
        xc%use_gradient = .true.
      case ( XC_FAMILY_MGGA)
        xc%use_gradient = .true.
        xc%use_laplacian = .true.
        xc%use_kinetic_energy = .true.
      case (XC_FAMILY_HYB_GGA, XC_FAMILY_HYB_MGGA)
        print '(A, A)', "Hybrid is not available:", trim(libxc_name)
        stop
      case default
        print '(A, A)', "Unknown Family:", trim(libxc_name)
        stop
      end select

      return
    end subroutine init_libxc

  end subroutine init_xc

End module lxc
