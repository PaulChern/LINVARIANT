Module timer
  Use constants
  Implicit none
  Contains

  subroutine get_date(cdate, ctime)
    !=======================================================
    !
    !! Returns two strings containing the date and the time
    !! in human-readable format. Uses a standard f90 call.
    !
    !=======================================================
    implicit none
    character(len=9), intent(out) :: cdate
    !! The date
    character(len=9), intent(out) :: ctime
    !! The time

    character(len=3), dimension(12) :: months
    data months/'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', &
      'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/
    integer date_time(8)
    !
    call date_and_time(values=date_time)
    !
    write (cdate, '(i2,a3,i4)') date_time(3), months(date_time(2)), date_time(1)
    write (ctime, '(i2.2,":",i2.2,":",i2.2)') date_time(5), date_time(6), date_time(7)

  end subroutine get_date

  Subroutine get_walltime(wctime)
    Implicit none
    Real*4, Intent(out) :: wctime
    Integer             :: r, c

    call system_clock(c, r)

    wctime = Real(c) / r

  End Subroutine get_walltime

  function get_wallclocktime()
    !==================================================================!
    !                                                                  !
    ! Returns elapsed wall clock time in seconds since its first call  !
    !                                                                  !
    !===================================================================
    implicit none

    real(dp) :: get_wallclocktime

    integer  :: c0, c1
    integer  :: rate
    logical  :: first = .true.
    save first, rate, c0

    if (first) then

      call system_clock(c0, rate)
      get_wallclocktime = 0.0_dp
      first = .false.
    else
      call system_clock(c1)
      get_wallclocktime = real(c1 - c0)/real(rate)
    endif
    return
  end function get_wallclocktime

  function get_time()
    !===========================================================
    !
    !! Returns elapsed CPU time in seconds since its first call.
    !! Uses standard f90 call
    !
    !===========================================================
    implicit none

    real(dp) :: get_time

    ! t0 contains the time of the first call
    ! t1 contains the present time
    real(kind=dp) :: t0, t1
    logical :: first = .true.
    save first, t0
    !
    call cpu_time(t1)
    !
    if (first) then
      t0 = t1
      get_time = 0.0_dp
      first = .false.
    else
      get_time = t1 - t0
    endif
    return
  end function get_time

End Module timer
