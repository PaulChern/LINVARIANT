Module Outputs 
  Use LINVARIANT
  Use Parameters
  Use Lattice
  Use Energy
  Use Aux

  implicit none
  
  Contains
    Include "ExportModelData.f90"
  
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
            write(FileHandle) (Fields(i,ifield,ix,iy,iz), i=1,FieldDim)
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
    Real*8         :: e0ij(3,3), euij(3,3,cgrid%n1, cgrid%n2, cgrid%n3)
    
    open(FileHandle,file=trim(Solver)//'.out/'//filename,form='formatted',status='unknown')
    Do ifield = 1, NumField
      Write(FileHandle, "("//trim(int2str(FieldDim))//"E25.15)") ((/0.0, 0.0, 0.0/), i=1,FieldDim)
    End do
    write(FileHandle, "(3E25.15)") e0ij(1,1), e0ij(2,2), e0ij(3,3)
    write(FileHandle, "(3E25.15)") e0ij(2,3), e0ij(1,3), e0ij(1,2)
    do iz = 1, cgrid%n3
      do iy = 1, cgrid%n2
        do ix = 1, cgrid%n1
          write(FileHandle, "(3I10)") ix, iy, iz, GridUnfold(2, ix, iy, iz)
          !          write(FileHandle, "(E25.15)")  EOnSite(ix,iy,iz,Fields)
          do ifield = 1, NumField
            write(FileHandle, "("//trim(int2str(FieldDim))//"E25.15)") (Fields(i,ifield,ix,iy,iz), i=1,FieldDim)
          end do
        end do
      end do
    end do

    close(FileHandle)

    open(FileHandle,file=trim(Solver)//'.out/Final_euij.dat',form='formatted',status='unknown')
    do iz = 1, cgrid%n3
      do iy = 1, cgrid%n2
        do ix = 1, cgrid%n1
          euij = GetHeterostructureStrain(Fields, ix, iy, iz, 0)
          write(FileHandle, "(3I10)") ix, iy, iz
          write(FileHandle, "(3E25.15)") (euij(i,:,ix,iy,iz), i=1,3)
        end do
      end do
    end do
    close(FileHandle)    

  End Subroutine
  
  Subroutine BinaryToData(filename, processor)
    Implicit None
    
    character(*)    :: filename
    Integer         :: io_mode, io_out, io_avg, processor
    Integer         :: ix, iy, iz, occ, nmc=0
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
      
      do igrid = 1, cgrid%npts
        Read(io_mode, END=999) ix, iy, iz, occ
        do ifield = 1, NumField
          Read(io_mode, END=999) (field(i, ifield), i=1,FieldDim)
          avgconfig(:,ifield,ix,iy,iz) = avgconfig(:,ifield,ix,iy,iz) + field(i,ifield)
        end do
        Write(io_out, "(3I10)") ix, iy, iz, occ
        do ifield = 1, NumField
          Write(io_out, "("//trim(int2str(FieldDim))//"E25.15)") (field(i, ifield), i=1,FieldDim)
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
          Write(io_avg, "(3I10)") ix, iy, iz, occ
          Do ifield = 1, NumField
            Write(io_avg, "("//trim(int2str(FieldDim))//"E25.15)") (avgconfig(i, ifield, ix, iy, iz)/nmc, i=1,FieldDim)
          End do
        end do
      end do
    end do

    Close(io_avg)
    
  END Subroutine BinaryToData

  Subroutine GetOpAvg(Fields, avg)
    Implicit None

    Real*8, Intent(in)    :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8, Intent(inout) :: avg(FieldDim, 2, NumField)
    Integer               :: ix, iy, iz, icell, num
    Integer               :: ifield, i, j

    avg = 0.0D0

    do ifield = 1, NumField
      do j = 1, FieldDim
        num = 0
        Do icell = 1, cgrid%ncells
          ix = GridFold(1, icell)
          iy = GridFold(2, icell)
          iz = GridFold(3, icell)
          if (ix > nint(GateField(1)) .and. ix <= cgrid%n1-nint(GateField(1)) .and. &
              iy > nint(GateField(1)) .and. iy <= cgrid%n2-nint(GateField(1)) .and. &
              iz > nint(GateField(1)) .and. iz <= cgrid%n3-nint(GateField(1))) then
            num = num + 1
              avg(j, 1, ifield) = avg(j, 1, ifield) + FieldsBinary(j,ifield)*Fields(j,ifield,ix,iy,iz)
              avg(j, 2, ifield) = avg(j, 2, ifield) + FieldsBinary(j,ifield)*Abs(Fields(j,ifield,ix,iy,iz))
          end if
        End do
        avg(j, :, ifield) = avg(j, :, ifield)/num
      end do
    end do

  END Subroutine GetOpAvg

  Subroutine GetAbsEuij(Fields, euij)
    Implicit None

    Real*8, Intent(in)    :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8, Intent(inout) :: euij(3, 3)
    Real*8                :: tmp(3,3,cgrid%n1, cgrid%n2, cgrid%n3)
    Integer               :: ix, iy, iz, icell

    euij = 0.0D0
    Do icell = 1, cgrid%ncells
      ix = GridFold(1, icell)
      iy = GridFold(2, icell)
      iz = GridFold(3, icell)
      tmp = GetHeterostructureStrain(Fields, ix, iy, iz, 0)
      euij = euij + Abs(tmp(:,:,ix,iy,iz))
    end do

    euij = euij/cgrid%npts

  END Subroutine GetAbsEuij

  Subroutine GetObservables(filename, order, processor)
    Implicit None
    
    character(*)    :: filename
    Integer         :: order, io_mode, io_out, io_phase, processor
    Integer         :: ix, iy, iz, icell, imc=0, num
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

    phase = 0.0D0
    eta = 0.0D0
 
    do While (.True.)
      Observables = 0.0D0
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
          Read(io_mode, END=999) (Fields(i, ifield, ix, iy, iz), i=1,FieldDim)
        end do
      end do
      
      do ifield = 1, NumField
        do j = 1, FieldDim
!          Observables(j, ifield) = Sum(Fields(j,ifield,:,:,:)**order)/cgrid%npts
          num = 0
        Do icell = 1, cgrid%ncells
          ix = GridFold(1, icell)
          iy = GridFold(2, icell)
          iz = GridFold(3, icell)
            if (ix > nint(GateField(1)) .and. ix <= cgrid%n1-nint(GateField(1)) .and. &
                iy > nint(GateField(1)) .and. iy <= cgrid%n2-nint(GateField(1)) .and. &
                iz > nint(GateField(1)) .and. iz <= cgrid%n3-nint(GateField(1))) then
              num = num + 1
              If (ifield .eq. NumField) then
                Observables(j, ifield) = Observables(j, ifield) + FieldsBinary(j,ifield)*Abs(Fields(j,ifield,ix,iy,iz))
              else
                Observables(j, ifield) = Observables(j, ifield) + FieldsBinary(j,ifield)*Fields(j,ifield,ix,iy,iz)
              End if
            end if
          End do
          Observables(j, ifield) = Observables(j, ifield)/num
        end do
        Write(io_out, "("//trim(int2str(FieldDim))//"E25.15)") (Observables(i, ifield), i=1,FieldDim)
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
      Write(io_phase, "("//trim(int2str(FieldDim))//"E25.15)") (phase(i, ifield), i=1,FieldDim)
    End do
    Do ifield = 1, NumField
      Write(io_phase, "("//trim(int2str(FieldDim))//"E25.15)") ((/0.0, 0.0, 0.0/), i=1,FieldDim)
    End do
    Close(io_phase)

    Call ExportMLModel(MLFile, Epot/cgrid%npts, phase, eta)
 
  END Subroutine GetObservables

  Subroutine MCMC_log(Fields, e0ij, istep, time, dene, udamp, etadamp)
    implicit none
    Real*8, Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8, Intent(in) :: e0ij(3,3), dene, udamp(NumField), etadamp
    Real*4, Intent(in) :: time
    Integer,Intent(in) :: istep
    Integer            :: ifield, i
    Real*8             :: avg(FieldDim, 2, NumField), abseuij(3,3)
    Real*8             :: Epot, volume

    Call GetOpAvg(Fields, avg)
    Call GetAbsEuij(Fields, abseuij)
    Epot = GetEtot(Fields, e0ij)/cgrid%npts
    volume = GetStrainVolume(e0ij)
    write(*,'(A15,I10,F10.2,A10)') trim(Solver)//" step: ", istep, time, "(second)"
    write(*, *) "------------------------------------------------------"
    write(*,'(A15,'//trim(int2str(size(udamp)))//'F10.6,F10.6)') "damp: ", udamp, etadamp
    write(*,'(A15,F10.6,A10,F8.4)') "Epot: ", Epot, "conv: ", dene
    do ifield = 1, NumField
      write(*,'(A12,I1,A2,'//trim(int2str(FieldDim))//'F10.6,A10,'//trim(int2str(FieldDim))//'F10.6)') "op", ifield, ": ", avg(:,1,ifield), "abs: ", avg(:,2,ifield)
    end do
    write(*, *) " "
    write(*,'(A15,3F10.6,A10,3F10.6)') "eii: ", e0ij(1,1), e0ij(2,2), e0ij(3,3), "abs: ", abseuij(1,1), abseuij(2,2), abseuij(3,3)
    write(*,'(A15,3F10.6,A10,3F10.6)') "eij: ", e0ij(2,3), e0ij(1,3), e0ij(1,2), "abs: ", abseuij(2,3), abseuij(1,3), abseuij(1,2)
    write(*,'(A15,F10.6)')  "volume: ", volume
    write(*, *) "------------------------------------------------------"

  End Subroutine MCMC_log

  Subroutine WLMC_log(Fields, e0ij, istep, ThermoSteps, time, udamp, etadamp, wl_f)
    implicit none
    Real*8, Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8, Intent(in) :: e0ij(3,3), udamp(NumField), etadamp, wl_f
    Real*4, Intent(in) :: time
    Integer,Intent(in) :: istep, ThermoSteps
    Integer            :: ifield, i
    Real*8             :: avg(FieldDim, 2, NumField)
    Real*8             :: Epot, volume

    Call GetOpAvg(Fields, avg)
    Epot = GetEtot(Fields, e0ij)
    volume = GetStrainVolume(e0ij)
    If(istep.le.ThermoSteps) write(*,'(A15,I10)') "WL ThermoSteps:", ThermoSteps
    write(*,'(A9,I10,F10.2,A10)') trim(Solver)//" step: ", istep, time, "(second)"
    write(*,'(A9,'//trim(int2str(size(udamp)))//'F10.6,F10.6,A10,F8.4)') "damp: ", udamp, etadamp, "conv: ", Log(wl_f)
    write(*,'(A9,F10.6)') "Epot: ", Epot
    do ifield = 1, NumField - 1
      write(*,'(A12,I1,A2,'//trim(int2str(FieldDim))//'F10.6,A10,'//trim(int2str(FieldDim))//'F10.6)') "op", ifield, ": ", avg(:,1,ifield), "abs: ", avg(:,2,ifield)
    end do
    write(*,'(A9,3F10.6)') "eii: ", e0ij(1,1), e0ij(2,2), e0ij(3,3)
    write(*,'(A9,3F10.6)') "eij: ", e0ij(2,3), e0ij(1,3), e0ij(1,2)
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
    Integer            :: ifield, i
    Real*8             :: avg(FieldDim, 2, NumField)
    Real*8             :: Etot, Epot, Ekin(2), Tk(2), volume
    
    Call GetOpAvg(Fields, avg)
    Epot = GetEtot(Fields, e0ij)
    Ekin = GetEkin(dFieldsdt, de0ijdt)
    Etot = Epot + Ekin(1)
    volume = GetStrainVolume(e0ij)
    Tk = Thermometer(dFieldsdt, de0ijdt)
    write(*,'(A9,I10,F6.2,A10,A4,2F10.4)') trim(Solver)//" step: ", istep, time, "(second)", "Tk:", Tk
    write(*,'(A9,'//trim(int2str(size(gm)))//'F10.6)')"gm: ", gm
    write(*,'(A9,3F10.6)') "Etot: ", Etot, Epot, Ekin(1)
    write(*,'(A9,4F10.6)') "Ekin: ", ResolveEkin(dFieldsdt, de0ijdt)
    do ifield = 1, NumField
      write(*,'(A12,I1,A2,'//trim(int2str(FieldDim))//'F10.6,A10,'//trim(int2str(FieldDim))//'F10.6)') "op", ifield, ": ", avg(:,1,ifield), "abs: ", avg(:,2,ifield)
    end do
    write(*,'(A9,3F10.6)') "eii: ", e0ij(1,1), e0ij(2,2), e0ij(3,3)
    write(*,'(A9,3F10.6)') "eij: ", e0ij(2,3), e0ij(1,3), e0ij(1,2)
    write(*,'(A9,F10.6)')"volume: ", volume

  End Subroutine MD_log

  Subroutine Solver_finalize(Fields, e0ij, dFieldsdt, de0ijdt)
    implicit none
    Real*8, Intent(in) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8, Intent(in) :: dFieldsdt(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8, Intent(in) :: e0ij(3,3), de0ijdt(3,3)
    character(len=10)  :: FileIndex
 
    FileIndex = int2str5(NODE_ME)
    
    Call WriteFinal('FinalVelocity-'//trim(FileIndex)//'.dat', dFieldsdt, de0ijdt)
    Call WriteFinal('FinalConfig-'//trim(FileIndex)//'.dat', Fields, e0ij)
    Call BinaryToData(trim(Solver), NODE_ME)
    Call GetObservables("Observables", 1, NODE_ME)
  End Subroutine Solver_finalize
  
End Module Outputs
