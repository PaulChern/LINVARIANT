module symmetry
  Use Constants
  use spglib_f08
  Contains

  subroutine write_syminfo( max_num_sym, num_atom, &
       & lattice, symprec, atom_types, positions, &
       & mesh, shifted, is_time_reversal )
  
    implicit none
  
    ! Arguments ------------------------------------
    ! scalars
    integer, intent(in) :: num_atom, max_num_sym, is_time_reversal
    real(dp), intent(in) :: symprec
    ! arrays
    integer, intent(in), dimension(num_atom) :: atom_types
    integer, intent(in), dimension(3) :: mesh, shifted
    real(dp), intent(in), dimension(3, 3) :: lattice
    real(dp), intent(in), dimension(3, num_atom) :: positions
    ! Local variables-------------------------------
    ! scalars
    integer :: i, j, counter, weight, space_group, num_sym, indent
    integer :: num_ir_grid
    character(len=21) :: international
    character(len=7) :: schoenflies
    character(len=30) :: space
    character(len=128) :: fmt
    ! arrays
    integer, dimension(3, 3, max_num_sym) :: rotations
    integer, dimension(3, mesh(1)*mesh(2)*mesh(3)) :: grid_point
    integer, dimension(mesh(1)*mesh(2)*mesh(3)) :: map
    real(dp), dimension(3, max_num_sym) :: translations
  
    type(SpglibDataset) :: dset
    type(SpglibSpacegroupType) :: spgtype
  
    !**************************************************************************
  
    space = "                              "
  
  !   The allocatable components of dset get deallocated on scope exit
    dset = spg_get_dataset(lattice, positions, atom_types, num_atom, symprec)
  
    num_sym = dset % n_operations
  
    indent = 1
    if (dset % spacegroup_number /= 0) then
       spgtype = spg_get_spacegroup_type(dset % hall_number)
  
       print('(a, "space_group: ", i3)'), space(1:indent*2), dset % spacegroup_number
       print('(a, "international: ", a, a)' ), space(1:indent*2), trim(dset % international_symbol)
       print('(a, "schoenflies: ", a)'), space(1:indent*2), trim(spgtype % schoenflies)
    else
       print '("Space group could not be found,")'
    end if
  
    print ('(a, "atom-type:")'), space(1:indent*2)
    do i = 1, num_atom
       print('(a, "- { type: ", i3, "}")'), space(1:indent*2), atom_types(i)
    end do
    print('(a, "real-basis:")'), space(1:indent*2)
    do i = 1, 3
       print('(a, "- [", f19.14, ", ", f19.14, ", ", f19.14, "]")'), space(1:indent*2), lattice(i, :)
    end do
    print('(a, "position:")'), space(1:indent*2)
    do i = 1, num_atom
       print('(a, "- [", f17.14, ", ", f17.14, ", ", f17.14, "]")'), space(1:indent*2), positions(:, i)
    end do
    print('(a, "operation:")'), space(1:indent*2)
    do i = 1, num_sym
       print('(a, "- rotation: #", i4)'), space(1:indent*2), i
       do j = 1, 3
          print('(a, "  - [", i3,",", i3,",", i3,"]")'), space(1:indent*2), dset % rotations(:,j, i)
       end do
       print('(a, "  translation: [ ", f10.7,", ", f10.7,", ", f10.7,"]")'), space(1:indent*2), dset % translations(:,i)
    end do
  
  
    write(fmt, '(i0,a)') num_atom, "(x,i0,':',i0))"
    print "(a,'wyckoffs:         ', "//fmt,space(1:indent*2),(i, dset % wyckoffs(i), i = 1, dset % n_atoms)
    print "(a,'equivalent_atoms: ', "//fmt,space(1:indent*2),(i, dset % equivalent_atoms(i) + 1, i = 1, dset % n_atoms)
  
    print('(a, "reciprocal-mesh: [", i3, ", ", i3, ", ", i3, "]")'), space(1:indent*2), mesh(:)
    print('(a, "- shifted: [", i3, ", ", i3, ", ", i3, "]")'), space(1:indent*2), shifted(:)
    print('(a, "- is_time_reversal: ", i3)'), space(1:indent*2), is_time_reversal
  
    num_ir_grid = spg_get_ir_reciprocal_mesh( grid_point, map, &
        & mesh, shifted, is_time_reversal, lattice, positions, &
        & atom_types, num_atom, symprec )
  
    print('(a, "- num-ir-grid-point:", i4)'), space(1:indent*2), num_ir_grid
    print('(a, "- ir-grid-point:")'), space(1:indent*2)
    counter = 0
  
    do i = 1, product(mesh)
       if ( i-1 == map(i) ) then
          ! Ad-hoc and intuitive implementation of weight. Not optimal for very
          ! large size
          weight = count(map == i-1)
          counter = counter + 1
  
          print '(i5, 3(x,f9.6), 2x,i4,2x, f9.6, 3(x,i4))', counter, &
           & real(shifted + 2*grid_point(:, i))/real(2*mesh), map(i)+1, &
           & real(weight)/real(product(mesh)), grid_point(:, i)
       end if
    end do
  
  end subroutine write_syminfo
  
  subroutine LoadSymmetry(num_atom,lattice,positions,atom_types,kmesh,shifted,is_time_reversal)
    Use Parameters
  
    implicit none
  
    ! max_sym is the expected maximum number of symmetry operations.
    ! This can be very big in a supercell.
    integer :: max_num_sym=192
  
    integer,  intent(in)  :: num_atom
    real(dp), intent(in)  :: lattice(3,3)
    real(dp), intent(in)  :: positions(3,num_atom)
    integer,  intent(in)  :: atom_types(num_atom)
    integer,  intent(in)  :: kmesh(3)
    integer,  intent(in)  :: shifted(3)
    integer,  intent(in)  :: is_time_reversal

    type(SpglibDataset)   :: symdata
    type(SpglibSpacegroupType) :: spgtype
    integer               :: rotations(3,3,192)
    real(dp)              :: translations(3,192)
    character(len=128)    :: fmt
  
    integer               :: i, j, indent, counter, weight
    character(len=30)     :: space
    integer               :: grid_point(3,product(kmesh))
    integer               :: map(product(kmesh))
    integer               :: num_sym, num_ir_grid
  
    real(dp) :: a
  
    space = "                              "

    symdata = spg_get_dataset(lattice, positions, atom_types, num_atom, symprec)

    num_sym = symdata%n_operations

    indent = 1
    if (symdata%spacegroup_number.ne.0) then
       spgtype = spg_get_spacegroup_type(symdata%hall_number)
  
       print('(a, "space_group: ", i3)'), space(1:indent*2), symdata%spacegroup_number
       print('(a, "international: ", a, a)' ), space(1:indent*2), trim(symdata%international_symbol)
       print('(a, "schoenflies: ", a)'), space(1:indent*2), trim(spgtype%schoenflies)
    else
       print '("Space group could not be found,")'
    end if

    print('(a, "operation:")'), space(1:indent*2)
    do i = 1, num_sym
       print('(a, "- rotation: #", i4)'), space(1:indent*2), i
       do j = 1, 3
          print('(a, "  - [", i3,",", i3,",", i3,"]")'), space(1:indent*2), symdata%rotations(:,j, i)
       end do
       print('(a, "  translation: [ ", f10.7,", ", f10.7,", ", f10.7,"]")'), space(1:indent*2), symdata%translations(:,i)
    end do

    write(fmt, '(i0,a)') num_atom, "(x,i0,':',i0))"
    print "(a,'wyckoffs:         ', "//fmt,space(1:indent*2),(i, symdata%wyckoffs(i), i = 1, symdata%n_atoms)
    print "(a,'equivalent_atoms: ', "//fmt,space(1:indent*2),(i, symdata%equivalent_atoms(i) + 1, i = 1, symdata%n_atoms)

    print('(a, "reciprocal-kmesh: [", i3, ", ", i3, ", ", i3, "]")'), space(1:indent*2), kmesh(:)
    print('(a, "- shifted: [", i3, ", ", i3, ", ", i3, "]")'), space(1:indent*2), shifted(:)
    print('(a, "- is_time_reversal: ", i3)'), space(1:indent*2), is_time_reversal

    num_ir_grid = spg_get_ir_reciprocal_mesh( grid_point, map, &
        & kmesh, shifted, is_time_reversal, lattice, positions, &
        & atom_types, num_atom, symprec )

    print('(a, "- num-ir-grid-point:", i4)'), space(1:indent*2), num_ir_grid
    print('(a, "- ir-grid-point:")'), space(1:indent*2)
    counter = 0

    do i = 1, product(kmesh)
       if ( i-1 == map(i) ) then
          ! Ad-hoc and intuitive implementation of weight. Not optimal for very
          ! large size
          weight = count(map == i-1)
          counter = counter + 1

          print '(i5, 3(x,f9.6), 2x,i4,2x, f9.6, 3(x,i4))', counter, &
           & real(shifted + 2*grid_point(:, i))/real(2*kmesh), map(i)+1, &
           & real(weight)/real(product(kmesh)), grid_point(:, i)
       end if
    end do

  end subroutine LoadSymmetry
end module symmetry
