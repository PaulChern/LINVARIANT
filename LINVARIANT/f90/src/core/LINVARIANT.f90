Module LINVARIANT
  
  Use Parameters
  Use Constants
  Use Energy
  Use Lattice
  Use Aux
  
  Implicit none
  Contains
 
   Subroutine ThinFilmBCS(Fields)
     Implicit none
     Real*8,  Intent(inout) :: Fields(FieldDim, NumField, cgrid_a%n1+cgrid_b%n1, cgrid_a%n2+cgrid_b%n2, cgrid_a%n3+cgrid_b%n3)
     Integer                :: ix, iy, iz, ifield

     If (cgrid_b%n1.ne.0) then
       do iz = 1, cgrid_a%n3
         do iy = 1, cgrid_a%n2
           do ifield = 1, NumIRFields
             Fields(1, ifield, 1+cgrid_a%n1, iy, iz) = (1 - screening)*Fields(1, ifield, cgrid_a%n1, iy, iz)
             Fields(1, ifield, cgrid_a%n1+cgrid_b%n1, iy, iz) = (1 - screening)*Fields(1, ifield, 1, iy, iz)
           End do
         End do
       End do
     End if

     If (cgrid_b%n2.ne.0) then
       do iz = 1, cgrid_a%n3
         do ix = 1, cgrid_a%n1
           do ifield = 1, NumIRFields
             Fields(2, ifield, ix, 1+cgrid_a%n2, iz) = (1 - screening)*Fields(2, ifield, ix, cgrid_a%n2, iz)
             Fields(2, ifield, ix, cgrid_a%n2+cgrid_b%n2, iz) = (1 - screening)*Fields(2, ifield, ix, 1, iz)
           End do
         End do
       End do
     End If

     If (cgrid_b%n3.ne.0) then
       do iy = 1, cgrid_a%n2
         do ix = 1, cgrid_a%n1
           do ifield = 1, NumIRFields
             Fields(3, ifield, ix, iy, 1+cgrid_a%n3) = (1 - screening)*Fields(3, ifield, ix, iy, cgrid_a%n3)
             Fields(3, ifield, ix, iy, cgrid_a%n3+cgrid_b%n3) = (1 - screening)*Fields(3, ifield, ix, iy, 1)
           End do
         End do
       End do
     End if

  End Subroutine ThinFilmBCS

  Subroutine RemoveGridDrifts(Fields)
    Implicit none
    Real*8,  Intent(inout) :: Fields(FieldDim, NumField, cgrid_a%n1+cgrid_b%n1, cgrid_a%n2+cgrid_b%n2, cgrid_a%n3+cgrid_b%n3)
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
    Real*8,  Intent(in)    :: Fields(FieldDim, NumField, cgrid_a%n1+cgrid_b%n1, cgrid_a%n2+cgrid_b%n2, cgrid_a%n3+cgrid_b%n3)
    Real*8                 :: Fields1D(FieldDim*NumField*(cgrid_a%npts+cgrid_b%npts))
    Integer                :: ix, iy, iz, ifield, i, id
    Integer                :: NG1, NG2, NG3

    NG1 = cgrid_a%n1+cgrid_b%n1
    NG2 = cgrid_a%n2+cgrid_b%n2
    NG3 = cgrid_a%n3+cgrid_b%n3

    do iz = 1, NG3
      do iy = 1, NG2
        do ix = 1, NG1
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
    Real*8,  Intent(in)    :: Fields1D(FieldDim*NumField*(cgrid_a%npts+cgrid_b%npts))
    Real*8                 :: Fields(FieldDim,NumField,cgrid_a%n1+cgrid_b%n1, cgrid_a%n2+cgrid_b%n2, cgrid_a%n3+cgrid_b%n3)
  
    Fields = Reshape(Fields1D, (/FieldDim, NumField, cgrid_a%n1+cgrid_b%n1, cgrid_a%n2+cgrid_b%n2, cgrid_a%n3+cgrid_b%n3/))
  
  End Function AllFieldsFrom1D
  
  Function Thermometer(dFieldsdt, de0ijdt) Result(TK)
    Implicit none
    Real*8,  Intent(in) :: dFieldsdt(FieldDim, NumField, cgrid_a%n1+cgrid_b%n1, cgrid_a%n2+cgrid_b%n2, cgrid_a%n3+cgrid_b%n3)
    Real*8,  Intent(in) :: de0ijdt(3,3)
    Real*8              :: TK(2), Ekin(2) 
    Integer             :: i, ifield, ix, iy, iz, TotalDim, StrainDim
  
    TotalDim = 0
    do i = 1, NumField
      If(.Not.FrozenQ(i)) TotalDim = TotalDim+cgrid%npts*Sum(FieldsBinary(:,i))
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
    Real*8,  Intent(inout) :: dFieldsdt(FieldDim, NumField, cgrid_a%n1+cgrid_b%n1, cgrid_a%n2+cgrid_b%n2, cgrid_a%n3+cgrid_b%n3)
    Real*8,  Intent(inout) :: de0ijdt(3,3)
    Real*8                 :: dFieldsdt_tmp(FieldDim, NumField, cgrid_a%n1+cgrid_b%n1, cgrid_a%n2+cgrid_b%n2, cgrid_a%n3+cgrid_b%n3)
    Real*8                 :: de0ijdt_tmp(3,3), detadt(6)
    Real*8                 :: TK, Ekin(NumField+1)
    Real*8                 :: rand
    Integer                :: i, ifield, ix, iy, iz, NumGrid, StrainDim
  
    NumGrid = cgrid%npts
  
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
    If(TK .gt. 0.0D0) then
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
