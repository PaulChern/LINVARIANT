Module fft
  ! FourierParameters -> {1, 1}
  Use Inputs
  Use Aux
  Use Parameters
  Use Constants
  Contains

  Subroutine fft3d_c2c(a,b,field_r,n1,n2,n3,field_q,m1,m2,m3)
    implicit none
    include 'fftw3.f'
    Integer,Intent(in)      :: a, b
    Integer,Intent(in)      :: n1, n2, n3
    Integer,Intent(in)      :: m1, m2, m3
    complex(dp),Intent(in)  :: field_r(n1,n2,n3)
    complex(dp),Intent(out) :: field_q(m1,m2,m3)
    complex(dp)             :: work(n1,n2,n3)
    Integer                 :: lx1, lx2, ly1, ly2, lz1, lz2
    integer*8 plan

    call dfftw_plan_dft_3d(plan,n1,n2,n3,field_r,work,b,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, field_r, work)
    call dfftw_destroy_plan(plan)

    work = work/(n1*n2*n3)**((1.0-a)/2)

    If(n1.eq.m1.and.n2.eq.m2.and.n3.eq.m3) then
      field_q = work
    else
      lx1 = n1/2+1
      lx2 = m1-n1+n1/2+2
      ly1 = n2/2+1
      ly2 = m2-n2+n2/2+2
      lz1 = n3/2+1
      lz2 = m3-n3+n3/2+2
      field_q = cmplx(0.0_dp,0.0_dp)
      field_q(1:lx1,1:ly1,1:lz1) = work(1:lx1,1:ly1,1:lz1)             ! (0,0,0)
      field_q(lx2:m1,ly2:m2,lz2:m3) = work(lx1+1:n1,ly1+1:n2,lz1+1:n3) ! (1,1,1)
      field_q(lx2:m1,1:ly1,1:lz1) = work(lx1+1:n1,1:ly1,1:lz1)         ! (1,0,0)
      field_q(1:lx1,ly2:m2,1:lz1) = work(1:lx1,ly1+1:n2,1:lz1)         ! (0,1,0)
      field_q(1:lx1,1:ly1,lz2:m3) = work(1:lx1,1:ly1,lz1+1:n3)         ! (0,0,1)
      field_q(1:lx1,ly2:m2,lz2:m3) = work(1:lx1,ly1+1:n2,lz1+1:n3)     ! (0,1,1)
      field_q(lx2:m1,1:ly1,lz2:m3) = work(lx1+1:n1,1:ly1,lz1+1:n3)     ! (1,0,1)
      field_q(lx2:m1,ly2:m2,1:lz1) = work(lx1+1:n1,ly1+1:n2,1:lz1)     ! (1,1,0)
    End if
    field_q = ((1.0*m1*m2*m3)/(n1*n2*n3))*field_q

  End Subroutine fft3d_c2c

  Subroutine fft3d_r2c(a,b,field_r,n1,n2,n3,field_q,m1,m2,m3)
    implicit none
    include 'fftw3.f'
    Integer,Intent(in)      :: a, b
    Integer,Intent(in)      :: n1, n2, n3
    Integer,Intent(in)      :: m1, m2, m3
    real(dp),Intent(in)     :: field_r(n1,n2,n3)
    complex(dp),Intent(out) :: field_q(m1,m2,m3)
    complex(dp)             :: work(n1,n2,n3)
    Integer                 :: lx1, lx2, ly1, ly2, lz1, lz2
    integer*8 plan

    call dfftw_plan_dft_r2c_3d(plan,n1,n2,n3,field_r,work,b,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, field_r, work)
    call dfftw_destroy_plan(plan)

    work = work/(n1*n2*n3)**((1.0-a)/2)

    If(n1.eq.m1.and.n2.eq.m2.and.n3.eq.m3) then
      field_q = work
    else
      lx1 = n1/2+1
      lx2 = m1-n1+n1/2+2
      ly1 = n2/2+1
      ly2 = m2-n2+n2/2+2
      lz1 = n3/2+1
      lz2 = m3-n3+n3/2+2
      field_q = cmplx(0.0_dp,0.0_dp)
      field_q(1:lx1,1:ly1,1:lz1) = work(1:lx1,1:ly1,1:lz1)             ! (0,0,0)
      field_q(lx2:m1,ly2:m2,lz2:m3) = work(lx1+1:n1,ly1+1:n2,lz1+1:n3) ! (1,1,1)
      field_q(lx2:m1,1:ly1,1:lz1) = work(lx1+1:n1,1:ly1,1:lz1)         ! (1,0,0)
      field_q(1:lx1,ly2:m2,1:lz1) = work(1:lx1,ly1+1:n2,1:lz1)         ! (0,1,0)
      field_q(1:lx1,1:ly1,lz2:m3) = work(1:lx1,1:ly1,lz1+1:n3)         ! (0,0,1)
      field_q(1:lx1,ly2:m2,lz2:m3) = work(1:lx1,ly1+1:n2,lz1+1:n3)     ! (0,1,1)
      field_q(lx2:m1,1:ly1,lz2:m3) = work(lx1+1:n1,1:ly1,lz1+1:n3)     ! (1,0,1)
      field_q(lx2:m1,ly2:m2,1:lz1) = work(lx1+1:n1,ly1+1:n2,1:lz1)     ! (1,1,0)
    End if
    field_q = ((1.0*m1*m2*m3)/(n1*n2*n3))*field_q

  End Subroutine fft3d_r2c

  Subroutine fft1d_c2c(a,b,field_t,n1,field_w,m1)
    implicit none
    include 'fftw3.f'
    Integer,Intent(in)      :: a, b, n1, m1
    complex(dp),Intent(in)  :: field_t(n1)
    complex(dp),Intent(out) :: field_w(m1)
    complex(dp)             :: work(n1)
    integer*8 plan

    call dfftw_plan_dft_1d(plan,n1,field_t,work,b,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, field_t, work)
    call dfftw_destroy_plan(plan)

    work = work/n1**((1.0-a)/2)

    If(n1.eq.m1) then
      field_w = work
    else
      field_w = cmplx(0.0_dp,0.0_dp)
      field_w(1:n1/2+1) = work(1:n1/2+1)
      field_w(m1-n1+n1/2+2:m1) = work(n1/2+2:n1)
    End if
    field_w = ((1.0*m1)/n1)*field_w


  End Subroutine fft1d_c2c

  Subroutine fft1d_r2c(a,b,field_t,n1,field_w,m1)
    implicit none
    include 'fftw3.f'
    Integer,Intent(in)      :: a, b, m1, n1
    real(dp),Intent(in)     :: field_t(n1)
    complex(dp),Intent(out) :: field_w(m1)
    complex(dp)             :: work(n1)
    integer*8               :: plan
    Integer                 :: i

    call dfftw_plan_dft_r2c_1d(plan,n1,field_t,work,b,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, field_t, work)
    call dfftw_destroy_plan(plan)

    Do i = 1, n1/2
      if(i+1.ne.n1+1-i) work(n1+1-i) = Conjg(work(i+1))
    End do

    work = work/n1**((1.0-a)/2)

    If(n1.eq.m1) then
      field_w = work
    else
      field_w = cmplx(0.0_dp,0.0_dp)
      field_w(1:n1/2+1) = work(1:n1/2+1)
      field_w(m1-n1+n1/2+2:m1) = work(n1/2+2:n1)
    End if
    field_w = ((1.0*m1)/n1)*field_w

  End Subroutine fft1d_r2c

  Subroutine Spectrum(processor)
    Implicit None

    Integer         :: io_mode=11, io_out=12
    Integer         :: ix, iy, iz, imc, nmc, processor
    Integer         :: ifield, i, igrid
    Integer         :: qx, qy, qz
    Real*8          :: Etot, Epot, Ekin, norm
    Real*8          :: e11, e22, e33, e23, e13, e12
    Real(dp)        :: field_r(FieldDim,NumField,cgrid%n1,cgrid%n2,cgrid%n3,(NumSteps-ThermoSteps)/TapeRate)
    complex(dp)     :: field_q(cgrid%n1,cgrid%n2,cgrid%n3)
    complex(dp),allocatable :: spectra(:,:,:,:,:,:), wq(:)

    field_r = 0.0
    nmc = 0

    call io_file_unit(io_mode)
    Open(io_mode,file=trim(Solver)//'.out/'//'trajectory_binary-'//trim(int2str5(processor))//'.dat',&
    form='unformatted',status='old')

    do While (.True.)
      nmc = nmc + 1

      Read(io_mode, END=999) Etot, Epot, Ekin
      Read(io_mode, END=999) e11, e22, e33
      Read(io_mode, END=999) e23, e13, e12

      do igrid = 1, cgrid%npts
        Read(io_mode, END=999) ix, iy, iz
        do ifield = 1, NumField
          Read(io_mode, END=999) (field_r(i, ifield, ix, iy, iz, nmc), i=1,FieldDim)
        end do
      end do
    End do
    999 Close(io_mode)
    nmc = nmc - 1
  
    Allocate(spectra(nmc,cgrid%n1,cgrid%n2,cgrid%n3,3,NumField))
    Allocate(wq(nmc))

    Do imc = 1, nmc
      Do ifield = 1, NumField
        do i = 1, FieldDim
          field_q = cmplx(0.0_dp,0.0_dp)
          norm=sqrt(sum(field_r(i,ifield,:,:,:,imc)**2)/cgrid%npts)
          call fft3d_r2c(1,-1,field_r(i,ifield,:,:,:,imc),cgrid%n1,cgrid%n2,cgrid%n3,field_q,cgrid%n1,cgrid%n2,cgrid%n3)
          spectra(imc,:,:,:,i,ifield) = field_q/norm/cgrid%npts
        end do
      End do
    End do

    call io_file_unit(io_out)
    Open(io_out,file=trim(Solver)//'.out/FFTq-'//trim(int2str5(processor))//'.dat',status='unknown')

    do ifield = 1, NumField
    do i = 1, FieldDim
      do qz = 1, cgrid%n3
        do qy = 1, cgrid%n2
          do qx = 1, cgrid%n1
            norm=sqrt(sum(spectra(:,qx,qy,qz,i,ifield)**2)/nmc)
            wq = cmplx(0.0_dp,0.0_dp)
            call fft1d_c2c(0,-1,spectra(:,qx,qy,qz,i,ifield),nmc,wq,nmc)
            wq = wq/norm/nmc
            Write(io_out, "(5I10)") ifield, i, qx, qy, qz
            do imc = 1, nmc
              Write(io_out, "(2F25.15)") wq(imc)
            end do
          end do
        end do
      end do
    end do
    end do


    Close(io_out)
    Deallocate(spectra)
    Deallocate(wq)

  End subroutine Spectrum

End Module fft
