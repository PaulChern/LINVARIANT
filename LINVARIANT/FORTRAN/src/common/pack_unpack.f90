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
module pack_unpack
  implicit none

  type array_shape
    integer :: nbeg
    integer :: nend
    integer :: nsize
  end type

  public :: array_shape
  public :: create_array_shape
  public :: pack_data, unpack_data
  public :: copy_data

  interface pack_data
    module procedure pack_data_3d_real8
    module procedure pack_data_3d_complex8
    module procedure pack_data_5d_real8
    module procedure pack_data_5d_complex8
  end interface

  interface unpack_data
    module procedure unpack_data_3d_real8
    module procedure unpack_data_3d_complex8
    module procedure unpack_data_5d_real8
    module procedure unpack_data_5d_complex8
  end interface

  interface copy_data
    module procedure copy_data_3d_real8
    module procedure copy_data_3d_complex8
    module procedure copy_data_4d_real8
    module procedure copy_data_4d_complex8
    module procedure copy_data_5d_real8
    module procedure copy_data_5d_complex8
    module procedure copy_data_6d_real8
    module procedure copy_data_6d_complex8
  end interface

contains
  function create_array_shape(nbeg,nend,nsize) result(s)
    implicit none
    integer, intent(in)           :: nbeg,nend
    integer, intent(in), optional :: nsize
    type(array_shape) :: s
    s%nbeg = nbeg
    s%nend = nend
    if (present(nsize)) then
      s%nsize = nsize
    else
      s%nsize = abs(nend - nbeg) + 1
    end if
  end function

  subroutine pack_data_3d_real8(nsrc_shape,ncopy_range,src,dst)
    implicit none
    type(array_shape), intent(in)  :: nsrc_shape(3)
    type(array_shape), intent(in)  :: ncopy_range(3)
    real(8),           intent(in)  :: src(nsrc_shape(1)%nbeg:nsrc_shape(1)%nend, &
                                          nsrc_shape(2)%nbeg:nsrc_shape(2)%nend, &
                                          nsrc_shape(3)%nbeg:nsrc_shape(3)%nend)
    real(8),           intent(out) :: dst(ncopy_range(1)%nsize, &
                                          ncopy_range(2)%nsize, &
                                          ncopy_range(3)%nsize)
    call copy_data(src(ncopy_range(1)%nbeg:ncopy_range(1)%nend,  &
                       ncopy_range(2)%nbeg:ncopy_range(2)%nend,  &
                       ncopy_range(3)%nbeg:ncopy_range(3)%nend), &
                   dst)
  end subroutine

  subroutine pack_data_3d_complex8(nsrc_shape,ncopy_range,src,dst)
    implicit none
    type(array_shape), intent(in)  :: nsrc_shape(3)
    type(array_shape), intent(in)  :: ncopy_range(3)
    complex(8),        intent(in)  :: src(nsrc_shape(1)%nbeg:nsrc_shape(1)%nend, &
                                          nsrc_shape(2)%nbeg:nsrc_shape(2)%nend, &
                                          nsrc_shape(3)%nbeg:nsrc_shape(3)%nend)
    complex(8),        intent(out) :: dst(ncopy_range(1)%nsize, &
                                          ncopy_range(2)%nsize, &
                                          ncopy_range(3)%nsize)
    call copy_data(src(ncopy_range(1)%nbeg:ncopy_range(1)%nend,  &
                       ncopy_range(2)%nbeg:ncopy_range(2)%nend,  &
                       ncopy_range(3)%nbeg:ncopy_range(3)%nend), &
                   dst)
  end subroutine

  subroutine pack_data_5d_real8(nsrc_shape,ncopy_range,src,dst)
    implicit none
    type(array_shape), intent(in)  :: nsrc_shape(5)
    type(array_shape), intent(in)  :: ncopy_range(5)
    real(8),           intent(in)  :: src(nsrc_shape(1)%nbeg:nsrc_shape(1)%nend, &
                                          nsrc_shape(2)%nbeg:nsrc_shape(2)%nend, &
                                          nsrc_shape(3)%nbeg:nsrc_shape(3)%nend, &
                                          nsrc_shape(4)%nbeg:nsrc_shape(4)%nend, &
                                          nsrc_shape(5)%nbeg:nsrc_shape(5)%nend)
    real(8),           intent(out) :: dst(ncopy_range(1)%nsize, &
                                          ncopy_range(2)%nsize, &
                                          ncopy_range(3)%nsize, &
                                          ncopy_range(4)%nsize, &
                                          ncopy_range(5)%nsize)
    call copy_data(src(ncopy_range(1)%nbeg:ncopy_range(1)%nend,  &
                       ncopy_range(2)%nbeg:ncopy_range(2)%nend,  &
                       ncopy_range(3)%nbeg:ncopy_range(3)%nend,  &
                       ncopy_range(4)%nbeg:ncopy_range(4)%nend,  &
                       ncopy_range(5)%nbeg:ncopy_range(5)%nend), &
                   dst)
  end subroutine

  subroutine pack_data_5d_complex8(nsrc_shape,ncopy_range,src,dst)
    implicit none
    type(array_shape), intent(in)  :: nsrc_shape(5)
    type(array_shape), intent(in)  :: ncopy_range(5)
    complex(8),        intent(in)  :: src(nsrc_shape(1)%nbeg:nsrc_shape(1)%nend, &
                                          nsrc_shape(2)%nbeg:nsrc_shape(2)%nend, &
                                          nsrc_shape(3)%nbeg:nsrc_shape(3)%nend, &
                                          nsrc_shape(4)%nbeg:nsrc_shape(4)%nend, &
                                          nsrc_shape(5)%nbeg:nsrc_shape(5)%nend)
    complex(8),        intent(out) :: dst(ncopy_range(1)%nsize, &
                                          ncopy_range(2)%nsize, &
                                          ncopy_range(3)%nsize, &
                                          ncopy_range(4)%nsize, &
                                          ncopy_range(5)%nsize)
    call copy_data(src(ncopy_range(1)%nbeg:ncopy_range(1)%nend,  &
                       ncopy_range(2)%nbeg:ncopy_range(2)%nend,  &
                       ncopy_range(3)%nbeg:ncopy_range(3)%nend,  &
                       ncopy_range(4)%nbeg:ncopy_range(4)%nend,  &
                       ncopy_range(5)%nbeg:ncopy_range(5)%nend), &
                   dst)
  end subroutine

  subroutine unpack_data_3d_real8(ndst_shape,ncopy_range,src,dst)
    implicit none
    type(array_shape), intent(in)  :: ndst_shape(3)
    type(array_shape), intent(in)  :: ncopy_range(3)
    real(8),           intent(in)  :: src(ncopy_range(1)%nsize, &
                                          ncopy_range(2)%nsize, &
                                          ncopy_range(3)%nsize)
    real(8),           intent(out) :: dst(ndst_shape(1)%nbeg:ndst_shape(1)%nend, &
                                          ndst_shape(2)%nbeg:ndst_shape(2)%nend, &
                                          ndst_shape(3)%nbeg:ndst_shape(3)%nend)
    call copy_data(src, &
                   dst(ncopy_range(1)%nbeg:ncopy_range(1)%nend, &
                       ncopy_range(2)%nbeg:ncopy_range(2)%nend, &
                       ncopy_range(3)%nbeg:ncopy_range(3)%nend))
  end subroutine

  subroutine unpack_data_3d_complex8(ndst_shape,ncopy_range,src,dst)
    implicit none
    type(array_shape), intent(in)  :: ndst_shape(3)
    type(array_shape), intent(in)  :: ncopy_range(3)
    complex(8),        intent(in)  :: src(ncopy_range(1)%nsize, &
                                          ncopy_range(2)%nsize, &
                                          ncopy_range(3)%nsize)
    complex(8),        intent(out) :: dst(ndst_shape(1)%nbeg:ndst_shape(1)%nend, &
                                          ndst_shape(2)%nbeg:ndst_shape(2)%nend, &
                                          ndst_shape(3)%nbeg:ndst_shape(3)%nend)
    call copy_data(src, &
                   dst(ncopy_range(1)%nbeg:ncopy_range(1)%nend, &
                       ncopy_range(2)%nbeg:ncopy_range(2)%nend, &
                       ncopy_range(3)%nbeg:ncopy_range(3)%nend))
  end subroutine

  subroutine unpack_data_5d_real8(ndst_shape,ncopy_range,src,dst)
    implicit none
    type(array_shape), intent(in)  :: ndst_shape(5)
    type(array_shape), intent(in)  :: ncopy_range(5)
    real(8),           intent(in)  :: src(ncopy_range(1)%nsize, &
                                          ncopy_range(2)%nsize, &
                                          ncopy_range(3)%nsize, &
                                          ncopy_range(4)%nsize, &
                                          ncopy_range(5)%nsize)
    real(8),           intent(out) :: dst(ndst_shape(1)%nbeg:ndst_shape(1)%nend, &
                                          ndst_shape(2)%nbeg:ndst_shape(2)%nend, &
                                          ndst_shape(3)%nbeg:ndst_shape(3)%nend, &
                                          ndst_shape(4)%nbeg:ndst_shape(4)%nend, &
                                          ndst_shape(5)%nbeg:ndst_shape(5)%nend)
    call copy_data(src, &
                   dst(ncopy_range(1)%nbeg:ncopy_range(1)%nend, &
                       ncopy_range(2)%nbeg:ncopy_range(2)%nend, &
                       ncopy_range(3)%nbeg:ncopy_range(3)%nend, &
                       ncopy_range(4)%nbeg:ncopy_range(4)%nend, &
                       ncopy_range(5)%nbeg:ncopy_range(5)%nend))
  end subroutine

  subroutine unpack_data_5d_complex8(ndst_shape,ncopy_range,src,dst)
    implicit none
    type(array_shape), intent(in)  :: ndst_shape(5)
    type(array_shape), intent(in)  :: ncopy_range(5)
    complex(8),        intent(in)  :: src(ncopy_range(1)%nsize, &
                                          ncopy_range(2)%nsize, &
                                          ncopy_range(3)%nsize, &
                                          ncopy_range(4)%nsize, &
                                          ncopy_range(5)%nsize)
    complex(8),        intent(out) :: dst(ndst_shape(1)%nbeg:ndst_shape(1)%nend, &
                                          ndst_shape(2)%nbeg:ndst_shape(2)%nend, &
                                          ndst_shape(3)%nbeg:ndst_shape(3)%nend, &
                                          ndst_shape(4)%nbeg:ndst_shape(4)%nend, &
                                          ndst_shape(5)%nbeg:ndst_shape(5)%nend)
    call copy_data(src, &
                   dst(ncopy_range(1)%nbeg:ncopy_range(1)%nend, &
                       ncopy_range(2)%nbeg:ncopy_range(2)%nend, &
                       ncopy_range(3)%nbeg:ncopy_range(3)%nend, &
                       ncopy_range(4)%nbeg:ncopy_range(4)%nend, &
                       ncopy_range(5)%nbeg:ncopy_range(5)%nend))
  end subroutine


  subroutine copy_data_3d_real8(src,dst)
    implicit none
    real(8), intent(in)  :: src(:,:,:)
    real(8), intent(out) :: dst(:,:,:)
    integer :: ix,iy,iz
    integer :: nx,ny,nz

    nz = size(src,3)
    ny = size(src,2)
    nx = size(src,1)

!$omp parallel do collapse(2) default(none) &
!$omp          private(ix,iy,iz) &
!$omp          firstprivate(nx,ny,nz) &
!$omp          shared(src,dst)
    do iz=1,nz
    do iy=1,ny
    do ix=1,nx
      dst(ix,iy,iz) = src(ix,iy,iz)
    end do
    end do
    end do
!$omp end parallel do
  end subroutine

  subroutine copy_data_3d_complex8(src,dst)
    implicit none
    complex(8), intent(in)  :: src(:,:,:)
    complex(8), intent(out) :: dst(:,:,:)
    integer :: ix,iy,iz
    integer :: nx,ny,nz

    nz = size(src,3)
    ny = size(src,2)
    nx = size(src,1)

!$omp parallel do collapse(2) default(none) &
!$omp          private(ix,iy,iz) &
!$omp          firstprivate(nx,ny,nz) &
!$omp          shared(src,dst)
    do iz=1,nz
    do iy=1,ny
    do ix=1,nx
      dst(ix,iy,iz) = src(ix,iy,iz)
    end do
    end do
    end do
!$omp end parallel do
  end subroutine

  subroutine copy_data_4d_real8(src,dst)
    implicit none
    real(8), intent(in)  :: src(:,:,:,:)
    real(8), intent(out) :: dst(:,:,:,:)
    integer :: nx,ny,nz,nw
    integer :: ix,iy,iz,iw

    nw = size(src,4)
    nz = size(src,3)
    ny = size(src,2)
    nx = size(src,1)

!$omp parallel do collapse(3) default(none) &
!$omp          private(ix,iy,iz,iw) &
!$omp          firstprivate(nx,ny,nz,nw) &
!$omp          shared(src,dst)
    do iw=1,nw
    do iz=1,nz
    do iy=1,ny
    do ix=1,nx
      dst(ix,iy,iz,iw) = src(ix,iy,iz,iw)
    end do
    end do
    end do
    end do
!$omp end parallel do
  end subroutine

  subroutine copy_data_4d_complex8(src,dst)
    implicit none
    complex(8), intent(in)  :: src(:,:,:,:)
    complex(8), intent(out) :: dst(:,:,:,:)
    integer :: nx,ny,nz,nw
    integer :: ix,iy,iz,iw

    nw = size(src,4)
    nz = size(src,3)
    ny = size(src,2)
    nx = size(src,1)

!$omp parallel do collapse(3) default(none) &
!$omp          private(ix,iy,iz,iw) &
!$omp          firstprivate(nx,ny,nz,nw) &
!$omp          shared(src,dst)
    do iw=1,nw
    do iz=1,nz
    do iy=1,ny
    do ix=1,nx
      dst(ix,iy,iz,iw) = src(ix,iy,iz,iw)
    end do
    end do
    end do
    end do
!$omp end parallel do
  end subroutine

  subroutine copy_data_5d_real8(src,dst)
    implicit none
    real(8), intent(in)  :: src(:,:,:,:,:)
    real(8), intent(out) :: dst(:,:,:,:,:)
    integer :: nx,ny,nz,nw,nl
    integer :: ix,iy,iz,iw,il

    nl = size(src,5)
    nw = size(src,4)
    nz = size(src,3)
    ny = size(src,2)
    nx = size(src,1)

!$omp parallel do collapse(4) default(none) &
!$omp          private(ix,iy,iz,iw,il) &
!$omp          firstprivate(nx,ny,nz,nw,nl) &
!$omp          shared(src,dst)
    do il=1,nl
    do iw=1,nw
    do iz=1,nz
    do iy=1,ny
    do ix=1,nx
      dst(ix,iy,iz,iw,il) = src(ix,iy,iz,iw,il)
    end do
    end do
    end do
    end do
    end do
!$omp end parallel do
  end subroutine

  subroutine copy_data_5d_complex8(src,dst)
    implicit none
    complex(8), intent(in)  :: src(:,:,:,:,:)
    complex(8), intent(out) :: dst(:,:,:,:,:)
    integer :: nx,ny,nz,nw,nl
    integer :: ix,iy,iz,iw,il

    nl = size(src,5)
    nw = size(src,4)
    nz = size(src,3)
    ny = size(src,2)
    nx = size(src,1)

!$omp parallel do collapse(4) default(none) &
!$omp          private(ix,iy,iz,iw,il) &
!$omp          firstprivate(nx,ny,nz,nw,nl) &
!$omp          shared(src,dst)
    do il=1,nl
    do iw=1,nw
    do iz=1,nz
    do iy=1,ny
    do ix=1,nx
      dst(ix,iy,iz,iw,il) = src(ix,iy,iz,iw,il)
    end do
    end do
    end do
    end do
    end do
!$omp end parallel do
  end subroutine

  subroutine copy_data_6d_real8(src,dst)
    implicit none
    real(8), intent(in)  :: src(:,:,:,:,:,:)
    real(8), intent(out) :: dst(:,:,:,:,:,:)
    integer :: nx,ny,nz,nw,nl,nm
    integer :: ix,iy,iz,iw,il,im

    nm = size(src,6)
    nl = size(src,5)
    nw = size(src,4)
    nz = size(src,3)
    ny = size(src,2)
    nx = size(src,1)

!$omp parallel do collapse(5) default(none) &
!$omp          private(ix,iy,iz,iw,il,im) &
!$omp          firstprivate(nx,ny,nz,nw,nl,nm) &
!$omp          shared(src,dst)
    do im=1,nm
    do il=1,nl
    do iw=1,nw
    do iz=1,nz
    do iy=1,ny
    do ix=1,nx
      dst(ix,iy,iz,iw,il,im) = src(ix,iy,iz,iw,il,im)
    end do
    end do
    end do
    end do
    end do
    end do
!$omp end parallel do
  end subroutine

  subroutine copy_data_6d_complex8(src,dst)
    implicit none
    complex(8), intent(in)  :: src(:,:,:,:,:,:)
    complex(8), intent(out) :: dst(:,:,:,:,:,:)
    integer :: nx,ny,nz,nw,nl,nm
    integer :: ix,iy,iz,iw,il,im

    nm = size(src,6)
    nl = size(src,5)
    nw = size(src,4)
    nz = size(src,3)
    ny = size(src,2)
    nx = size(src,1)

!$omp parallel do collapse(5) default(none) &
!$omp          private(ix,iy,iz,iw,il,im) &
!$omp          firstprivate(nx,ny,nz,nw,nl,nm) &
!$omp          shared(src,dst)
    do im=1,nm
    do il=1,nl
    do iw=1,nw
    do iz=1,nz
    do iy=1,ny
    do ix=1,nx
      dst(ix,iy,iz,iw,il,im) = src(ix,iy,iz,iw,il,im)
    end do
    end do
    end do
    end do
    end do
    end do
!$omp end parallel do
  end subroutine
end module
