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
module parallelization
  implicit none

  !!! multi process
  integer, public :: nproc_group_global
  integer, public :: nproc_id_global
  integer, public :: nproc_size_global

  !!! elapsed-time gap between compute node.
  real(8), private :: gap_time_system

  ! call once
  public :: setup_parallel
  public :: end_parallel

  ! util
  public :: adjust_elapse_time
  public :: get_thread_id
  public :: get_nthreads
  public :: is_distributed_parallel

  private :: find_elapse_time_gap

contains
  subroutine setup_parallel
    use communication
    implicit none
    call comm_init
    call comm_get_globalinfo(nproc_group_global, nproc_id_global, nproc_size_global)
    call find_elapse_time_gap
  end subroutine

  subroutine end_parallel
    use communication
    implicit none
    call comm_finalize
  end subroutine

  subroutine find_elapse_time_gap
    use communication
    use timer
    implicit none
    real(8) :: now, me
    call comm_sync_all
    now = get_time()
    me  = now
    call comm_sync_all
    call comm_get_min(now, nproc_group_global)
    gap_time_system = me - now
  end subroutine

  function adjust_elapse_time(time)
    implicit none
    real(8), intent(in) :: time
    real(8)             :: adjust_elapse_time
    adjust_elapse_time = time - gap_time_system
  end function

  function get_thread_id() result(nid)
#ifdef _OPENMP
    use omp_lib
#endif
    implicit none
    integer :: nid
#ifdef _OPENMP
    nid = omp_get_thread_num()
#else
    nid = 0
#endif
  end function

  function get_nthreads() result(nsize)
#ifdef _OPENMP
    use omp_lib
#endif
    implicit none
    integer :: nsize
#ifdef _OPENMP
    nsize = omp_get_max_threads()
#else
    nsize = 1
#endif
  end function

  function is_distributed_parallel() result(ret)
    implicit none
    logical :: ret
    ret = (nproc_size_global > 1)
  end function

end module
