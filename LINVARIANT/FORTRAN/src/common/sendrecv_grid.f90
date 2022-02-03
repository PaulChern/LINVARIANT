!
!  Copyright 2019-2020 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!


module sendrecv_grid
  use common

  implicit none

  public :: init_sendrecv_grid
  public :: dealloc_cache
  public :: update_overlap
  public :: s_sendrecv_grid

  integer, parameter :: iside_up   = 1
  integer, parameter :: iside_down = 2
  integer, parameter :: itype_send = 1
  integer, parameter :: itype_recv = 2

  integer, public, parameter :: srg_initialization = 0
  integer, public, parameter :: srg_pack           = 1
  integer, public, parameter :: srg_unpack         = 2
  integer, public, parameter :: srg_communication  = 4
  integer, public, parameter :: srg_all            = 7

  interface update_overlap
  module procedure update_overlap_real8
  module procedure update_overlap_complex8
  end interface


  contains

  ! Flip 1 to 2, 2 to 1:
  integer function flip(i)
    implicit none
    integer, intent(in) :: i
    flip = 3 - i
    return
  end function flip

  ! Prepare unique MPI tag number (3~8) for send/recv based on
  ! the communication direction (idir, iside) from the sender:
  integer function get_tag(iside, idir)
    implicit none
    integer, intent(in) :: idir, iside
    get_tag = 2 * idir + iside
    return
  end function get_tag

  ! Initializing s_sendrecv_grid structure:
  ! NOTE:
  ! * This subroutine is commonly used for real type and complex type.
  ! * The cache region MUST be allocated after this initialization 
  !   by using `alloc_cache_real8/complex8`.
  subroutine init_sendrecv_grid(srg, rg, nb, icomm, neig)
    implicit none
    type(s_sendrecv_grid), intent(inout) :: srg
    type(s_rgrid), intent(in) :: rg
    integer, intent(in) ::  nb, icomm, neig(1:2, 1:3)

    integer :: idir, iaxis
    integer :: is_block(1:3, 1:2, 1:2, 1:3)
    integer :: ie_block(1:3, 1:2, 1:2, 1:3)
    
    ! Calculate shape (upper and lower bounds) of overlapped region:
    ! NOTE:
    !
    !  nd nd         nd nd    U/D: up/down side of `idir`-th direction
    ! +==+--+-------+--+==+   S/R: send/recv region, representing the
    ! |DR|DS|       |US|UR|        inner and outer side of local grid.
    ! +==+--+-------+--+==+   
    !    is            ie     
    !    <------------->      
    !    Local grid size
    !
    ! is/ie_block(idir, iside, itype, iaxis): 
    !   lower/upper bounds of each blocks
    !   * idir: direction (1:x, 2:y, 3:z)
    !   * iside: kind of region (1:upside, 2:downside)
    !   * itype: kind of region (1:send, 2:recv)
    !   * iaxis: axis (1:x, 2:y, 3:z)
    do idir = 1, 3 ! 1:x,2:y,3:z
      do iaxis = 1, 3 ! 1:x,2:y,3:z
        if (idir == iaxis) then
          ! upside-send (US) block:
          is_block(idir, itype_send, iside_up, idir) = rg%ie(idir) - rg%nd + 1
          ie_block(idir, itype_send, iside_up, idir) = rg%ie(idir)
          ! upside-recv (UR) block:
          is_block(idir, itype_recv, iside_up, idir) = rg%ie(idir) + 1
          ie_block(idir, itype_recv, iside_up, idir) = rg%ie(idir) + rg%nd
          ! downside-send (DS) block:
          is_block(idir, itype_send, iside_down, idir) = rg%is(idir)
          ie_block(idir, itype_send, iside_down, idir) = rg%is(idir) + rg%nd - 1
          ! downside-recv (DR) block:
          is_block(idir, itype_recv, iside_down, idir) = rg%is(idir) - rg%nd
          ie_block(idir, itype_recv, iside_down, idir) = rg%is(idir) - 1
        else
          is_block(iaxis, :, :, idir) = rg%is(iaxis)
          ie_block(iaxis, :, :, idir) = rg%ie(iaxis)
        end if
      end do
    end do

    ! Assign to s_sendrecv_grid structure:
    srg%nb = nb
    srg%neig(1:2, 1:3) = neig(1:2, 1:3)
    srg%is_block(:, :, :, :) = is_block
    srg%ie_block(:, :, :, :) = ie_block
    srg%icomm = icomm
    srg%ireq_real8 = -1
    srg%ireq_complex8 = -1
    ! Flag for persistent communication
    srg%if_pcomm_real8_initialized = .false.
    srg%if_pcomm_complex8_initialized = .false.

    return
  end subroutine init_sendrecv_grid

  subroutine dealloc_cache(srg)
    use communication, only: comm_get_groupinfo, comm_free_reqs, comm_proc_null
    implicit none
    type(s_sendrecv_grid), intent(inout) :: srg
    integer :: idir, iside, itype
    integer :: myrank, nprocs

    ! Obtain myrank in communication group:
    call comm_get_groupinfo(srg%icomm, myrank, nprocs)

    do idir = 1, 3
      do iside = 1, 2
        ! Release persistent communication requests
        if (srg%neig(iside, idir) /= comm_proc_null) then
          if (srg%neig(iside, idir) /= myrank) then
            if (srg%if_pcomm_real8_initialized) &
              call comm_free_reqs(srg%ireq_real8(:, iside, idir))
            if (srg%if_pcomm_complex8_initialized) &
              call comm_free_reqs(srg%ireq_complex8(:, iside, idir))
          end if
        end if

        do itype = 1, 2
          ! Release allocated cache regions
          if (allocated(srg%cache(itype, iside, idir)%dbuf)) &
            deallocate( srg%cache(itype, iside, idir)%dbuf)
          if (allocated(srg%cache(itype, iside, idir)%zbuf)) &
            deallocate( srg%cache(itype, iside, idir)%zbuf)
        end do
      end do
    end do

    srg%if_pcomm_real8_initialized = .false.
    srg%if_pcomm_complex8_initialized = .false.

    return
  end subroutine

  subroutine update_overlap_real8(srg, rg, data, istage)
    use communication, only: comm_get_groupinfo, comm_start_all, comm_wait_all, comm_proc_null
    implicit none
    type(s_sendrecv_grid), intent(inout) :: srg
    type(s_rgrid), intent(in) :: rg
    real(8), intent(inout) :: data( &
      rg%is_array(1):rg%ie_array(1), &
      rg%is_array(2):rg%ie_array(2), &
      rg%is_array(3):rg%ie_array(3), &
      1:srg%nb)
    integer, intent(in), optional :: istage
    integer :: idir, iside
    integer :: myrank, nprocs, iphase

    if (present(istage)) then
      iphase = istage
    else
      iphase = srg_all
    end if

    ! Obtain myrank in communication group:
    call comm_get_groupinfo(srg%icomm, myrank, nprocs)

    ! Exchange the overlap region with the neighboring node (or opposite side of itself).
    if (.not. srg%if_pcomm_real8_initialized) call alloc_cache()
    do idir = 1, 3 ! 1:x,2:y,3:z
      do iside = 1, 2 ! 1:up,2:down
        if (srg%neig(iside, idir) /= comm_proc_null) then
          if (srg%neig(iside, idir) /= myrank) then
            if (iand(iphase,srg_pack) > 0) then
              ! Store the overlap region into the cache
              call pack_cache(iside, idir)
            end if

            ! In the first call of this subroutine, setup the persistent communication:
            if (.not. srg%if_pcomm_real8_initialized) call init_pcomm(iside, idir)

            if (iand(iphase,srg_communication) > 0) then
              ! Start to communication
              call comm_start_all(srg%ireq_real8(:, iside, idir))
            end if
          else
            if (iand(iphase,srg_pack) > 0) then
              ! NOTE: If neightboring nodes are itself (periodic with single proc),
              !       a simple side-to-side copy is used instead of the MPI comm.
              call copy_self(iside, idir)
            end if
          end if
        end if
      end do
    end do

    do idir = 1, 3 ! 1:x,2:y,3:z
      do iside = 1, 2 ! 1:up,2:down
        if (srg%neig(iside, idir) /= comm_proc_null) then
          if (srg%neig(iside, idir) /= myrank) then
            if (iand(iphase,srg_communication) > 0) then
              ! Wait for recieving
              call comm_wait_all(srg%ireq_real8(:, iside, idir))
            end if

            if (iand(iphase,srg_unpack) > 0) then
              ! Write back the recieved cache
              call unpack_cache(iside, idir)
            end if
          end if
        end if
      end do
    end do

    if (.not. srg%if_pcomm_real8_initialized) &
      srg%if_pcomm_real8_initialized = .true. ! Update if_pcomm_initialized
    
    return
    contains

    subroutine alloc_cache()
      implicit none
      integer :: jdir, jside, jtype
      integer :: is_b(3), ie_b(3)
      do jdir = 1, 3
        do jside = 1, 2
          do jtype = 1, 2
            is_b(1:3) = srg%is_block(1:3, jtype, jside, jdir)
            ie_b(1:3) = srg%ie_block(1:3, jtype, jside, jdir)
            allocate(srg%cache(jtype, jside, jdir)%dbuf( &
              is_b(1):ie_b(1), is_b(2):ie_b(2), is_b(3):ie_b(3), 1:srg%nb))
          end do
        end do
      end do
    end subroutine

    subroutine init_pcomm(jside, jdir)
      use communication, only: comm_send_init, comm_recv_init
      implicit none
      integer, intent(in) :: jdir, jside
      ! Send (and initialize persistent communication)
      srg%ireq_real8(itype_send, jside, jdir) = comm_send_init( &
        srg%cache(itype_send, jside, jdir)%dbuf, &
        srg%neig(jside, jdir), &
        get_tag(jside, jdir), &
        srg%icomm)
      ! Recv (and initialize persistent communication)
      srg%ireq_real8(itype_recv, jside, jdir) = comm_recv_init( &
        srg%cache(itype_recv, jside, jdir)%dbuf, &
        srg%neig(jside, jdir), &
        get_tag(flip(jside), jdir), & ! `jside` in sender side
        srg%icomm)
    end subroutine init_pcomm

    subroutine pack_cache(jside, jdir)
      use pack_unpack, only: copy_data
      implicit none
      integer, intent(in) :: jdir, jside
      integer :: is_s(1:3), ie_s(1:3) ! src region
      is_s(1:3) = srg%is_block(1:3, itype_send, jside, jdir)
      ie_s(1:3) = srg%ie_block(1:3, itype_send, jside, jdir)
      call copy_data( &
        data(is_s(1):ie_s(1), is_s(2):ie_s(2), is_s(3):ie_s(3), 1:srg%nb), &
        srg%cache(itype_send, jside, jdir)%dbuf)
    end subroutine pack_cache

    subroutine unpack_cache(jside, jdir)
      use pack_unpack, only: copy_data
      implicit none
      integer, intent(in) :: jdir, jside
      integer :: is_d(1:3), ie_d(1:3) ! dst region
      is_d(1:3) = srg%is_block(1:3, itype_recv, jside, jdir)
      ie_d(1:3) = srg%ie_block(1:3, itype_recv, jside, jdir)
      call copy_data( &
        srg%cache(itype_recv, jside, jdir)%dbuf, &
        data(is_d(1):ie_d(1), is_d(2):ie_d(2), is_d(3):ie_d(3), 1:srg%nb))
    end subroutine unpack_cache

    subroutine copy_self(jside, jdir)
      use pack_unpack, only: copy_data
      integer, intent(in) :: jdir, jside
      integer :: is_s(1:3), ie_s(1:3) ! src region
      integer :: is_d(1:3), ie_d(1:3) ! dst region
      is_s(1:3) = srg%is_block(1:3, itype_send, flip(jside), jdir)
      ie_s(1:3) = srg%ie_block(1:3, itype_send, flip(jside), jdir)
      is_d(1:3) = srg%is_block(1:3, itype_recv, jside, jdir)
      ie_d(1:3) = srg%ie_block(1:3, itype_recv, jside, jdir)
      call copy_data( &
        data(is_s(1):ie_s(1), is_s(2):ie_s(2), is_s(3):ie_s(3), 1:srg%nb), &
        data(is_d(1):ie_d(1), is_d(2):ie_d(2), is_d(3):ie_d(3), 1:srg%nb))
    end subroutine copy_self

  end subroutine update_overlap_real8


  subroutine update_overlap_complex8(srg, rg, data, istage)
    use communication, only: comm_get_groupinfo, &
      & comm_start_all, comm_wait_all, comm_proc_null
    implicit none
    type(s_sendrecv_grid), intent(inout) :: srg
    type(s_rgrid), intent(in) :: rg
    complex(8), intent(inout) :: data( &
      rg%is_array(1):rg%ie_array(1), &
      rg%is_array(2):rg%ie_array(2), &
      rg%is_array(3):rg%ie_array(3), &
      1:srg%nb)
    integer, intent(in), optional :: istage
    integer :: idir, iside
    integer :: myrank, nprocs, iphase

    if (present(istage)) then
      iphase = istage
    else
      iphase = srg_all
    end if

    ! Obtain myrank in communication group:
    call comm_get_groupinfo(srg%icomm, myrank, nprocs)

    ! Exchange the overlap region with the neighboring node (or opposite side of itself).
    if (.not. srg%if_pcomm_complex8_initialized) call alloc_cache()
    do idir = 1, 3 ! 1:x,2:y,3:z
      do iside = 1, 2 ! 1:up,2:down
        if (srg%neig(iside, idir) /= comm_proc_null) then
          if (srg%neig(iside, idir) /= myrank) then
            if (iand(iphase,srg_pack) > 0) then
              ! Store the overlap region into the cache
              call pack_cache(iside, idir)
            end if

            ! In the first call of this subroutine, setup the persistent communication:
            if (.not. srg%if_pcomm_complex8_initialized) call init_pcomm(iside, idir)

            if (iand(iphase,srg_communication) > 0) then
              ! Start to communication
              call comm_start_all(srg%ireq_complex8(:, iside, idir))
            end if
          else
            if (iand(iphase,srg_pack) > 0) then
              ! NOTE: If neightboring nodes are itself (periodic with single proc),
              !       a simple side-to-side copy is used instead of the MPI comm.
              call copy_self(iside, idir)
            end if
          end if
        end if
      end do
    end do

    do idir = 1, 3 ! 1:x,2:y,3:z
      do iside = 1, 2 ! 1:up,2:down
        if (srg%neig(iside, idir) /= comm_proc_null) then
          if (srg%neig(iside, idir) /= myrank) then
            if (iand(iphase,srg_communication) > 0) then
              ! Wait for recieving
              call comm_wait_all(srg%ireq_complex8(:, iside, idir))
            end if
            if (iand(iphase,srg_unpack) > 0) then
              ! Write back the recieved cache
              call unpack_cache(iside, idir)
            end if
          end if
        end if
      end do
    end do

    if (.not. srg%if_pcomm_complex8_initialized) &
      srg%if_pcomm_complex8_initialized = .true. ! Update if_pcomm_initialized

    return
    contains

    subroutine alloc_cache()
      implicit none
      integer :: jdir, jside, jtype
      integer :: is_b(3), ie_b(3)
      do jdir = 1, 3
        do jside = 1, 2
          do jtype = 1, 2
            is_b(1:3) = srg%is_block(1:3, jtype, jside, jdir)
            ie_b(1:3) = srg%ie_block(1:3, jtype, jside, jdir)
            allocate(srg%cache(jtype, jside, jdir)%zbuf( &
              is_b(1):ie_b(1), is_b(2):ie_b(2), is_b(3):ie_b(3), 1:srg%nb))
          end do
        end do
      end do
    end subroutine

    subroutine init_pcomm(jside, jdir)
      use communication, only: comm_send_init, comm_recv_init
      implicit none
      integer, intent(in) :: jdir, jside
      ! Send (and initialize persistent communication)
      srg%ireq_complex8(itype_send, jside, jdir) = comm_send_init( &
        srg%cache(itype_send, jside, jdir)%zbuf, &
        srg%neig(jside, jdir), &
        get_tag(jside, jdir), &
        srg%icomm)
      ! Recv (and initialize persistent communication)
      srg%ireq_complex8(itype_recv, jside, jdir) = comm_recv_init( &
        srg%cache(itype_recv, jside, jdir)%zbuf, &
        srg%neig(jside, jdir), &
        get_tag(flip(jside), jdir), & ! `jside` in sender
        srg%icomm)
    end subroutine init_pcomm

    subroutine pack_cache(jside, jdir)
      use pack_unpack, only: copy_data
      implicit none
      integer, intent(in) :: jdir, jside
      integer :: is_s(1:3), ie_s(1:3) ! src region
      is_s(1:3) = srg%is_block(1:3, itype_send, jside, jdir)
      ie_s(1:3) = srg%ie_block(1:3, itype_send, jside, jdir)
      call copy_data( &
        data(is_s(1):ie_s(1), is_s(2):ie_s(2), is_s(3):ie_s(3), 1:srg%nb), &
        srg%cache(itype_send, jside, jdir)%zbuf)
    end subroutine pack_cache

    subroutine unpack_cache(jside, jdir)
      use pack_unpack, only: copy_data
      implicit none
      integer, intent(in) :: jdir, jside
      integer :: is_d(1:3), ie_d(1:3) ! dst region
      is_d(1:3) = srg%is_block(1:3, itype_recv, jside, jdir)
      ie_d(1:3) = srg%ie_block(1:3, itype_recv, jside, jdir)
      call copy_data( &
        srg%cache(itype_recv, jside, jdir)%zbuf, &
        data(is_d(1):ie_d(1), is_d(2):ie_d(2), is_d(3):ie_d(3), 1:srg%nb))
    end subroutine unpack_cache

    subroutine copy_self(jside, jdir)
      use pack_unpack, only: copy_data
      integer, intent(in) :: jdir, jside
      integer :: is_s(1:3), ie_s(1:3) ! src region
      integer :: is_d(1:3), ie_d(1:3) ! dst region
      is_s(1:3) = srg%is_block(1:3, itype_send, flip(jside), jdir)
      ie_s(1:3) = srg%ie_block(1:3, itype_send, flip(jside), jdir)
      is_d(1:3) = srg%is_block(1:3, itype_recv, jside, jdir)
      ie_d(1:3) = srg%ie_block(1:3, itype_recv, jside, jdir)
      call copy_data( &
        data(is_s(1):ie_s(1), is_s(2):ie_s(2), is_s(3):ie_s(3), 1:srg%nb), &
        data(is_d(1):ie_d(1), is_d(2):ie_d(2), is_d(3):ie_d(3), 1:srg%nb))
    end subroutine copy_self

  end subroutine update_overlap_complex8

  subroutine create_sendrecv_neig(neig, info)
    use parameters
    use common
    use communication, only: comm_proc_null
    implicit none
    integer, intent(out) :: neig(1:2, 1:3)
    type(s_parallel_info), intent(in) :: info
    !
    integer :: idir,iside,idisp,iret

    do idir=1,3
    do iside=1,2
      select case(iside)
        case(1); idisp = 1
        case(2); idisp = -1
      end select

      iret = get_neighbour_rank(info, idir, idisp)

      if (iret < 0 .and. iperiodic == 0) then
        neig(iside,idir) = comm_proc_null
      else
        neig(iside,idir) = iret
      end if
    end do
    end do
    contains
      ! convert 5-D address to rank
      function get_neighbour_rank(info, idir, idisp) result(irank)
        use parameters
        use common
        implicit none
        type(s_parallel_info), intent(in) :: info
        integer,               intent(in) :: idir, idisp
        integer                           :: iaddr(5),irank
        integer                           :: ishape(5)

        ishape(1:5) = [info%nprgrid(1:3), info%nporbital, info%npk]
        iaddr(1:5)  = info%iaddress(1:5)
        iaddr(idir) = iaddr(idir) + idisp

        if (iperiodic == 3) then
          ! periodic boundary
          if (iaddr(idir) <  0)            iaddr(idir) = iaddr(idir) + ishape(idir)
          if (iaddr(idir) >= ishape(idir)) iaddr(idir) = iaddr(idir) - ishape(idir)
        else
          ! dirichlet boundary
          if (iaddr(idir) <  0)            iaddr(idir) = -1
          if (iaddr(idir) >= ishape(idir)) iaddr(idir) = -1
        end if

        if (iaddr(idir) < 0) then
          irank = -1
        else
          irank = info%imap(iaddr(1),iaddr(2),iaddr(3),iaddr(4),iaddr(5))
        end if
      end function get_neighbour_rank

  end subroutine create_sendrecv_neig

end module sendrecv_grid


