Module Outputs 
  Use LINVARIANT
  Use Parameters
  Use Lattice
  Use Energy
  Use Aux

  implicit none
  
  Contains
  
  subroutine fmkdir(newDirPath)
    implicit none
    
    character(len=*), intent(in) :: newDirPath
    character(len=256)           :: mkdirCmd
    logical                      :: dirExists
    
    ! Check if the directory exists first
    ! inquire(file=trim(newDirPath)//'/.', exist=dirExists)  ! Works with gfortran, but not ifort
    inquire(directory=newDirPath, exist=dirExists)         ! Works with ifort, but not gfortran
    
    if (dirExists) then
      write (*,*) "Directory already exists: '"//trim(newDirPath)//"'"
    else
      mkdirCmd = 'mkdir -p '//trim(newDirPath)
      write(*,'(a)') "Creating new directory: '"//trim(mkdirCmd)//"'"
      call system(mkdirCmd)
    endif
  end subroutine fmkdir
  
  Subroutine WriteBinary(FileHandle, Fields, dFieldsdt, e0ij, de0ijdt)
    Implicit None
    
    Integer        :: FileHandle
    Integer        :: ix, iy, iz, i, ifield
    Real*8         :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8         :: dFieldsdt(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8         :: Epot, Ekin(2), e0ij(3,3), de0ijdt(3,3)
    
    Epot = GetEtot(Fields, e0ij)
    Ekin = GetEkin(dFieldsdt, de0ijdt)
    
    write(FileHandle) Epot + Ekin(1), Epot, Ekin(1)
    write(FileHandle) e0ij(1,1), e0ij(2,2), e0ij(3,3)
    write(FileHandle) e0ij(2,3), e0ij(1,3), e0ij(1,2)
    do iz = 1, cgrid%n3
      do iy = 1, cgrid%n2
        do ix = 1, cgrid%n1
          write(FileHandle) ix, iy, iz
          !          write(FileHandle) EOnSite(ix,iy,iz,Fields)
          do ifield = 1, NumField
            write(FileHandle) (Fields(i,ifield,ix,iy,iz), i=1,FieldDimList(ifield))
          end do
        end do
      end do
    end do
  End Subroutine
  
  Subroutine WriteFinal(filename, Fields, e0ij)
    Implicit None
    
    character(*)   :: filename
    Integer        :: FileHandle = 1111
    Integer        :: ix, iy, iz, i, ifield
    Real*8         :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8         :: e0ij(3,3), euij(3,3)
    
    open(FileHandle,file=trim(Solver)//'.out/'//filename,form='formatted',status='unknown')
    write(FileHandle, "(3E25.15)") e0ij(1,1), e0ij(2,2), e0ij(3,3)
    write(FileHandle, "(3E25.15)") e0ij(2,3), e0ij(1,3), e0ij(1,2)
    do iz = 1, cgrid%n3
      do iy = 1, cgrid%n2
        do ix = 1, cgrid%n1
          write(FileHandle, "(3I10)") ix, iy, iz
          !          write(FileHandle, "(E25.15)")  EOnSite(ix,iy,iz,Fields)
          do ifield = 1, NumField
            write(FileHandle, "(3E25.15)") (Fields(i,ifield,ix,iy,iz), i=1,FieldDimList(ifield))
          end do
        end do
      end do
    end do
    close(FileHandle)

    open(FileHandle,file=trim(Solver)//'.out/Final_euij.dat',form='formatted',status='unknown')
    do iz = 1, cgrid%n3
      do iy = 1, cgrid%n2
        do ix = 1, cgrid%n1
          euij = StrainFromu(ix, iy, iz, Fields)
          write(FileHandle, "(3I10)") ix, iy, iz
          write(FileHandle, "(3E25.15)") (euij(i,:), i=1,3)
        end do
      end do
    end do
    close(FileHandle)    

  End Subroutine
  
  Subroutine BinaryToData(filename, processor)
    Implicit None
    
    character(*)    :: filename
    Integer         :: io_mode, io_out, io_avg, processor
    Integer         :: ix, iy, iz, nmc=0
    Integer         :: ifield, i, igrid
    Real*8          :: Etot, Epot, Ekin
    Real*8          :: e11, e22, e33, e23, e13, e12 
    Real*8          :: avgeij(6)
    Real*8          :: field(FieldDim, NumField)
    Real*8          :: avgconfig(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)

    
    io_mode = 11
    io_out = 12
    io_avg = 13
    
    Open(io_mode,file=trim(Solver)//'.out/'//'trajectory_binary-'//trim(int2str5(processor))//'.dat',&
    form='unformatted',status='old')
    Open(io_out,file=trim(Solver)//'.out/'//trim(filename)//'-'//trim(int2str5(processor))//'.dat',&
    status='unknown')
    Open(io_avg,file=trim(Solver)//'.out/AvgConfig-'//trim(int2str5(processor))//'.dat',&
    status='unknown')
   
    avgeij = 0.0d0
    avgconfig = 0.0d0
 
    do While (.True.)
      nmc = nmc + 1
      
      Read(io_mode, END=999) Etot, Epot, Ekin
      Read(io_mode, END=999) e11, e22, e33
      Read(io_mode, END=999) e23, e13, e12
      Write(io_out, "(3E25.15)") Etot, Epot, Ekin
      Write(io_out, "(3E25.15)") e11, e22, e33
      Write(io_out, "(3E25.15)") e23, e13, e12

      avgeij(1) = avgeij(1) + e11
      avgeij(2) = avgeij(2) + e22
      avgeij(3) = avgeij(3) + e33
      avgeij(4) = avgeij(4) + e23
      avgeij(5) = avgeij(5) + e13
      avgeij(6) = avgeij(6) + e12
      
      do igrid = 1, cgrid%n3*cgrid%n2*cgrid%n1
        Read(io_mode, END=999) ix, iy, iz
        do ifield = 1, NumField
          Read(io_mode, END=999) (field(i, ifield), i=1,FieldDimList(ifield))
          avgconfig(:,ifield,ix,iy,iz) = avgconfig(:,ifield,ix,iy,iz) + field(i,ifield)
        end do
        Write(io_out, "(3I10)") ix, iy, iz
        do ifield = 1, NumField
          Write(io_out, "(3E25.15)") (field(i, ifield), i=1,FieldDimList(ifield))
        end do
      end do
    End do
    
    999 Close(io_mode)
    Close(io_out)

    write(io_avg, "(3E25.15)") avgeij(1), avgeij(2), avgeij(3)
    write(io_avg, "(3E25.15)") avgeij(4), avgeij(5), avgeij(6)
    Do ix = 1, cgrid%n1
      Do iy = 1, cgrid%n2
        Do iz = 1, cgrid%n3
          Write(io_avg, "(3I10)") ix, iy, iz
          Do ifield = 1, NumField
            Write(io_avg, "(3E25.15)") (avgconfig(i, ifield, ix, iy, iz)/nmc, i=1,FieldDimList(ifield))
          End do
        end do
      end do
    end do

    Close(io_avg)
    
  END Subroutine BinaryToData
  
  Subroutine GetObservables(filename, order, processor)
    Implicit None
    
    character(*)    :: filename
    Integer         :: order, io_mode, io_out, io_phase, processor
    Integer         :: ix, iy, iz, imc=0, num
    Integer         :: ifield, i, j, igrid
    Real*8          :: Etot, Epot, Ekin
    Real*8          :: e11, e22, e33, e23, e13, e12
    Real*8          :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8          :: Observables(FieldDim, NumField)
    Real*8          :: phase(FieldDim, NumField), eta(6)
    
    io_mode = 11
    io_out = 12
    io_phase = 1112
    
    Open(io_mode,file=trim(Solver)//'.out/'//'trajectory_binary-'//trim(int2str5(processor))//'.dat',&
    form='unformatted',status='old')
    Open(io_out,file=trim(Solver)//'.out/'//trim(filename)//'-'//trim(int2str5(processor))//'.dat',&
    status='unknown')

    Observables = 0.0D0
    phase = 0.0D0
    eta = 0.0D0
 
    do While (.True.)
      imc = imc + 1
      
      Read(io_mode, END=999) Etot, Epot, Ekin
      Read(io_mode, END=999) e11, e22, e33
      Read(io_mode, END=999) e23, e13, e12
      Write(io_out, "(3E25.15)") Etot, Epot, Ekin
      Write(io_out, "(3E25.15)") e11, e22, e33
      Write(io_out, "(3E25.15)") e23, e13, e12
      eta = eta + (/e11, e22, e33, e23, e13, e12/)

      do igrid = 1, cgrid%n3*cgrid%n2*cgrid%n1
        Read(io_mode, END=999) ix, iy, iz
        do ifield = 1, NumField
          Read(io_mode, END=999) (Fields(i, ifield, ix, iy, iz), i=1,FieldDimList(ifield))
        end do
      end do
      
      do ifield = 1, NumField
        do j = 1, FieldDimList(ifield)
!          Observables(j, ifield) = Sum(Fields(j,ifield,:,:,:)**order)/cgrid%npts
          num = 0
          Do ix = 1, cgrid%n1
          Do iy = 1, cgrid%n2
          Do iz = 1, cgrid%n3
            if (ix > nint(GateField(1)) .and. ix <= cgrid%n1-nint(GateField(1)) .and. &
                iy > nint(GateField(1)) .and. iy <= cgrid%n2-nint(GateField(1)) .and. &
                iz > nint(GateField(1)) .and. iz <= cgrid%n3-nint(GateField(1))) then
              num = num + 1
              Observables(j, ifield) = Observables(j, ifield) + Fields(j,ifield,ix,iy,iz)
            end if
          End do
          End do
          End do
          Observables(j, ifield) = Observables(j, ifield)/num
        end do
        Write(io_out, "(3E25.15)") (Observables(i, ifield), i=1,FieldDimList(ifield))
      end do
      phase = phase + Observables
      
    End do
    
    999 Close(io_mode)
    Close(io_out)

    phase = phase/float(imc)
    eta   = eta/float(imc)
    Open(io_phase,file=trim(Solver)//'.out/phase-'//trim(int2str5(processor))//'.dat', status='unknown')
    Write(io_phase, "(3E25.15)") eta(1), eta(2), eta(3)
    Write(io_phase, "(3E25.15)") eta(4), eta(5), eta(6)
    Do ifield = 1, NumField
      Write(io_phase, "(3E25.15)") (phase(i, ifield), i=1,FieldDimList(ifield))
    End do
    Close(io_phase)
 
  END Subroutine GetObservables

  Subroutine MCMC_log(Fields, e0ij, istep, time, dene, udamp, etadamp)
    implicit none
    Real*8, Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8, Intent(in) :: e0ij(3,3), dene, udamp(NumField), etadamp
    Real*4, Intent(in) :: time
    Integer,Intent(in) :: istep
    Real*8             :: Epot, volume

    Epot = GetEtot(Fields, e0ij)
    volume = GetStrainVolume(e0ij)
    write(*,'(A9,I10,F10.2,A10)') trim(Solver)//" step: ", istep, time, "(second)"
    write(*,'(A9,4F10.6,A10,F8.4)') "damp: ", udamp, etadamp, "conv: ", dene
    write(*,'(A9,F10.6)') "Epot: ", Epot
    write(*,'(A9,F10.6,F10.6,F10.6)') "eii: ", e0ij(1,1), e0ij(2,2), e0ij(3,3)
    write(*,'(A9,F10.6,F10.6,F10.6)') "eij: ", e0ij(2,3), e0ij(1,3), e0ij(1,2)
    write(*,'(A9,F10.6)')"volume: ", volume

  End Subroutine MCMC_log

  Subroutine WLMC_log(Fields, e0ij, istep, ThermoSteps, time, udamp, etadamp, wl_f)
    implicit none
    Real*8, Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8, Intent(in) :: e0ij(3,3), udamp(NumField), etadamp, wl_f
    Real*4, Intent(in) :: time
    Integer,Intent(in) :: istep, ThermoSteps
    Real*8             :: Epot, volume

    Epot = GetEtot(Fields, e0ij)
    volume = GetStrainVolume(e0ij)
    If(istep.le.ThermoSteps) write(*,'(A15,I10)') "WL ThermoSteps:", ThermoSteps
    write(*,'(A9,I10,F10.2,A10)') trim(Solver)//" step: ", istep, time, "(second)"
    write(*,'(A9,4F10.6,A10,E8.4)') "damp: ", udamp, etadamp, "conv: ", Log(wl_f)
    write(*,'(A9,F10.6)') "Epot: ", Epot
    write(*,'(A9,F10.6,F10.6,F10.6)') "eii: ", e0ij(1,1), e0ij(2,2), e0ij(3,3)
    write(*,'(A9,F10.6,F10.6,F10.6)') "eij: ", e0ij(2,3), e0ij(1,3), e0ij(1,2)
    write(*,'(A9,F10.6)')"volume: ", volume

  End Subroutine WLMC_log

  Subroutine MD_log(Fields, e0ij, dFieldsdt, de0ijdt, istep, time, gm)
    Use Lattice
    implicit none
    Real*8, Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8, Intent(in) :: dFieldsdt(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8, Intent(in) :: e0ij(3,3), de0ijdt(3,3), gm(NumField+1)
    Real*4, Intent(in) :: time
    Integer,Intent(in) :: istep
    Real*8             :: Etot, Epot, Ekin(2), Tk(2), volume
    
    Epot = GetEtot(Fields, e0ij)
    Ekin = GetEkin(dFieldsdt, de0ijdt)
    Etot = Epot + Ekin(1)
    volume = GetStrainVolume(e0ij)
    Tk = Thermometer(dFieldsdt, de0ijdt)
    write(*,'(A9,I10,F6.2,A10,A4,2F10.4)') trim(Solver)//" step: ", istep, time, "(second)", "Tk:", Tk
    write(*,'(A9,4F10.6)')"gm: ", gm
    write(*,'(A9,3F10.6)') "Etot: ", Etot, Epot, Ekin(1)
    write(*,'(A9,4F10.6)') "Ekin: ", ResolveEkin(dFieldsdt, de0ijdt)
    write(*,'(A9,F10.6,F10.6,F10.6)') "eii: ", e0ij(1,1), e0ij(2,2), e0ij(3,3)
    write(*,'(A9,F10.6,F10.6,F10.6)') "eij: ", e0ij(2,3), e0ij(1,3), e0ij(1,2)
    write(*,'(A9,F10.6)')"volume: ", volume

  End Subroutine MD_log

  Subroutine Solver_finalize(Fields, e0ij, dFieldsdt, de0ijdt)
    implicit none
    Real*8, Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8, Intent(in) :: dFieldsdt(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8, Intent(in) :: e0ij(3,3), de0ijdt(3,3)
    character(len=10)  :: FileIndex
 
    FileIndex = int2str5(NODE_ME)
    
    Call WriteFinal('FinalConfig-'//trim(FileIndex)//'.dat', Fields, e0ij)
    Call WriteFinal('FinalVelocity-'//trim(FileIndex)//'.dat', dFieldsdt, de0ijdt)
    Call BinaryToData(trim(Solver), NODE_ME)
    Call GetObservables("Observables", 1, NODE_ME)
  End Subroutine Solver_finalize
  
End Module Outputs
