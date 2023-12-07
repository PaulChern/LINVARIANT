Module LINVARIANT
  
  Use Parameters
  Use Constants
  Use Energy
  Use Lattice
  Use Aux
  
  Implicit none
  Contains

  Subroutine LoadMesh
    Implicit none
    logical                  :: file_exists
    Integer                  :: NGx, NGy, NGz, id, ix, iy, iz, occ, scr

    INQUIRE(FILE='Supercell.inp', EXIST=file_exists)
    GridUnfold = 0
    GridFold   = 0

    If (file_exists) then
      open(ifileno, file='Supercell.inp', status='old')
      Read(ifileno, *) NGx, NGy, NGz
      If ((NGx.neqv.cgrid%n1).or.(NGy.neqv.cgrid%n2).or.(NGz.neqv.cgrid%n3)) then
        write(*,*) "Error: mismatch between Supercell.inp and supercell setup!"
        write(*,*) "cgrid: ", cgrid%n1, cgrid%n2, cgrid%n3
        write(*,*) "  BOX: ", NGx, NGy, NGz
        Call Abort
      End if
      do id = 1, cgrid%npts
        Read(ifileno, *) GridFold(:,id)
        ix = GridFold(1,id)
        iy = GridFold(2,id)
        iz = GridFold(3,id)
        GridUnfold(1, ix, iy, iz) = id
        GridUnfold(2:2+NumField,ix,iy,iz) = GridFold(4:4+NumField,id)
      end do
      close(ifileno)
    else
      do iz = 1, cgrid%n3
        do iy = 1, cgrid%n2
          do ix = 1, cgrid%n1
            id = ix + (iy - 1) * cgrid%n1 + (iz - 1) * cgrid%n2 * cgrid%n1
            GridFold(1, id) = ix
            GridFold(2, id) = iy
            GridFold(3, id) = iz
            GridFold(4:4+NumField, id) = 1
            GridUnfold(1, ix, iy, iz) = id
            GridUnfold(2:2+NumField, ix, iy, iz) = 1
          end do
        end do
      end do
    End if

    cgrid%ncells = Sum(GridFold(4,:))

  End Subroutine LoadMesh

  Subroutine ImposeBCS(Fields)
    Implicit none
    Real*8,  Intent(inout) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Integer                :: x0, y0, z0, ifield, id
    Integer                :: x1, y1, z1, xi1, yi1, zi1
    Integer                :: xn(3), yn(3), zn(3)
    Integer                :: i, j, k
    Real*8                 :: dist, occ

    !$OMP    PARALLEL DEFAULT(SHARED) PRIVATE(x0,y0,z0)
    !$OMP    DO
      do id = cgrid%ncells + 1, cgrid%npts
        x0 = GridFold(1, id)
        y0 = GridFold(2, id)
        z0 = GridFold(3, id)
        Fields(:, :, x0, y0, z0) = 0.0
      end do
    !$OMP    END DO
    !$OMP    END PARALLEL

    !$OMP    PARALLEL DEFAULT(SHARED) PRIVATE(x0,y0,z0,x1,y1,z1,xi1,yi1,zi1,xn,yn,zn)
    !$OMP    DO
      do id = cgrid%ncells + 1, cgrid%npts
        x0 = GridFold(1, id)
        y0 = GridFold(2, id)
        z0 = GridFold(3, id)

        x1 = (x0+1)-floor(real(x0+1-1)/real(cgrid%n1))*(cgrid%n1)
        y1 = (y0+1)-floor(real(y0+1-1)/real(cgrid%n2))*(cgrid%n2)
        z1 = (z0+1)-floor(real(z0+1-1)/real(cgrid%n3))*(cgrid%n3)
        xi1 = (x0-1)-floor(real(x0-1-1)/real(cgrid%n1))*(cgrid%n1)
        yi1 = (y0-1)-floor(real(y0-1-1)/real(cgrid%n2))*(cgrid%n2)
        zi1 = (z0-1)-floor(real(z0-1-1)/real(cgrid%n3))*(cgrid%n3)

        xn = (/xi1, x0, x1/)
        yn = (/yi1, y0, y1/)
        zn = (/zi1, z0, z1/)

        do ifield = 1, NumIRFields
          Fields(1, ifield, x0, y0, z0) = GridUnfold(2, x1, y0, z0)*Fields(1, ifield, x1, y0, z0) &
                                        + GridUnfold(2, xi1, y0, z0)*Fields(1, ifield, xi1, y0, z0)

          Fields(2, ifield, x0, y0, z0) = GridUnfold(2, x0, y1, z0)*Fields(2, ifield, x0, y1, z0) &
                                        + GridUnfold(2, x0, yi1, z0)*Fields(2, ifield, x0, yi1, z0)

          Fields(3, ifield, x0, y0, z0) = GridUnfold(2, x0, y0, z1)*Fields(3, ifield, x0, y0, z1) &
                                        + GridUnfold(2, x0, y0, zi1)*Fields(3, ifield, x0, y0, zi1)

          Fields(1:3, ifield, x0, y0, z0) = Fields(1:3, ifield, x0, y0, z0)/GridUnfold(2+ifield, x0, y0, z0)
        end do

        Fields(1, NumField, x0, y0, z0) = GridUnfold(2, x1, y0, z0)*Fields(1, NumField, x1, y0, z0) &
                                        + GridUnfold(2, xi1, y0, z0)*Fields(1, NumField, xi1, y0, z0)

        Fields(2, NumField, x0, y0, z0) = GridUnfold(2, x0, y1, z0)*Fields(2, NumField, x0, y1, z0) &
                                        + GridUnfold(2, x0, yi1, z0)*Fields(2, NumField, x0, yi1, z0)

        Fields(3, NumField, x0, y0, z0) = GridUnfold(2, x0, y0, z1)*Fields(3, NumField, x0, y0, z1) &
                                        + GridUnfold(2, x0, y0, zi1)*Fields(3, NumField, x0, y0, zi1)

        Fields(1:3, NumField, x0, y0, z0) = Fields(1:3, NumField, x0, y0, z0)/GridUnfold(2+NumField, x0, y0, z0)
      end do
    !$OMP    END DO
    !$OMP    END PARALLEL

  End Subroutine ImposeBCS

  Subroutine RemoveGridDrifts(Fields)
    Implicit none
    Real*8,  Intent(inout) :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8                 :: GridDriftx, GridDrifty, GridDriftz
    Integer                :: ix, iy, iz
  
    GridDriftx = Sum(Fields(1,NumField,1:cgrid%n1,1:cgrid%n2,1:cgrid%n3))/(cgrid%npts)
    GridDrifty = Sum(Fields(2,NumField,1:cgrid%n1,1:cgrid%n2,1:cgrid%n3))/(cgrid%npts)
    GridDriftz = Sum(Fields(3,NumField,1:cgrid%n1,1:cgrid%n2,1:cgrid%n3))/(cgrid%npts)
  
    do iz = 1, cgrid%n3
      do iy = 1, cgrid%n2
        do ix = 1, cgrid%n1
          Fields(1,NumField,ix,iy,iz) = Fields(1,NumField,ix,iy,iz) - GridDriftx
          Fields(2,NumField,ix,iy,iz) = Fields(2,NumField,ix,iy,iz) - GridDrifty
          Fields(3,NumField,ix,iy,iz) = Fields(3,NumField,ix,iy,iz) - GridDriftz
        end do
      end do
    end do
  
  End Subroutine RemoveGridDrifts
  
  Function AllFields1D(Fields) Result(Fields1D)
    Implicit none
    Real*8,  Intent(in)    :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8                 :: Fields1D(FieldDim*NumField*(cgrid%npts))
    Integer                :: ix, iy, iz, ifield, i, id

    do iz = 1, cgrid%n3
      do iy = 1, cgrid%n2
        do ix = 1, cgrid%n1
          do ifield = 1, NumField
            do i = 1, FieldDim
              id = (iz-1)*cgrid%n2*cgrid%n1*NumField*FieldDim &
              + (iy-1)*cgrid%n1*NumField*FieldDim &
              + (ix-1)*NumField*FieldDim &
              + (ifield-1)*FieldDim &
              + i
              Fields1D(id) = Fields(i,ifield,ix,iy,iz)
            end do
          end do
        end do
      end do
    end do
  
  End Function AllFields1D
  
  Function AllFieldsFrom1D(Fields1D) Result(Fields)
    Implicit none
    Real*8,  Intent(in)    :: Fields1D(FieldDim*NumField*(cgrid%npts))
    Real*8                 :: Fields(FieldDim,NumField,cgrid%n1, cgrid%n2, cgrid%n3)
  
    Fields = Reshape(Fields1D, (/FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3/))
  
  End Function AllFieldsFrom1D
  
  Function Thermometer(dFieldsdt, de0ijdt) Result(TK)
    Implicit none
    Real*8,  Intent(in) :: dFieldsdt(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8,  Intent(in) :: de0ijdt(3,3)
    Real*8              :: TK(2), Ekin(2) 
    Integer             :: i, ifield, ix, iy, iz, TotalDim, StrainDim
  
    TotalDim = 0
    do i = 1, NumField
      If(.Not.FrozenQ(i)) TotalDim = TotalDim+cgrid%ncells*Sum(FieldsBinary(:,i))
    end do
  
    StrainDim = 0
    do i = 1, 6
      If(.Not.CLAMPQ(i)) StrainDim = StrainDim+1
    end do
  
    Ekin = GetEkin(dFieldsdt, de0ijdt)
    TK(1) = 2*Ekin(1)*Hartree/k_bolt_ev/TotalDim
    TK(2) = 2*Ekin(2)*Hartree/k_bolt_ev/StrainDim
  
  End Function Thermometer
  
  Subroutine AdaptTk(dFieldsdt, de0ijdt)
    Implicit none
    Real*8,  Intent(inout) :: dFieldsdt(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8,  Intent(inout) :: de0ijdt(3,3)
    Real*8                 :: dFieldsdt_tmp(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8                 :: de0ijdt_tmp(3,3), detadt(6)
    Real*8                 :: TK, Ekin(NumField+1)
    Real*8                 :: rand
    Integer                :: i, ifield, ix, iy, iz, NumGrid, StrainDim
  
    NumGrid = cgrid%ncells
  
    Ekin = ResolveEkin(dFieldsdt, de0ijdt)
  
    Do i = 1, NumField
      TK = 2*Ekin(i)*Hartree/k_bolt_ev/(Sum(FieldsBinary(:,i))*NumGrid)
      dFieldsdt(:,i,:,:,:) = Sqrt(Temp/Tk)*dFieldsdt(:,i,:,:,:)
    End do
  
    StrainDim = 0
    do i = 1, 6
      If(.Not.CLAMPQ(i)) StrainDim = StrainDim+1
    end do
  
    TK = 2*Ekin(NumField+1)*Hartree/k_bolt_ev/StrainDim
    If (StrainDim.eq.0) then
      de0ijdt = 0.0
    Else If(TK .gt. 0.0D0) then
      de0ijdt = Sqrt(Temp/Tk)*de0ijdt
    Else
      do i = 1, 6
        Call random_number(rand)
        If(.Not.CLAMPQ(i)) then
          detadt(i) = rand-0.50d0
        else
          detadt(i) = 0.0D0
        end if
      end do
      de0ijdt = eta2eij(detadt)
      Ekin = ResolveEkin(dFieldsdt, de0ijdt)
      TK = 2*Ekin(NumField+1)*Hartree/k_bolt_ev/StrainDim
      de0ijdt = Sqrt(Temp/Tk)*de0ijdt
    End If

  
  End Subroutine AdaptTk
  
End Module LINVARIANT
