module optimization
  Use Energy
  Use Force
  Use Aux

  Implicit none
  Contains

  subroutine OptEnergyFunc(val, n, x, grad, need_gradient)
    Implicit none
    Real(dp)  :: val, x(n), grad(n)
    integer   :: n, need_gradient
    Real*8    :: Fields(FieldDim, NumField, cgrid%n1, cgrid%n2, cgrid%n3)
    Real*8    :: e0ij(3,3)
    Integer   :: ix, iy, iz

    opt_iter = opt_iter + 1

    Call SysVector2Fields(x, Fields, e0ij)
    if (need_gradient.ne.0) then
      grad  = GetGradient(Fields, e0ij)
    end if

    if (DipoleQ) then
      do iz = 1, cgrid%n3
        do iy = 1, cgrid%n2
          do ix = 1, cgrid%n1
            val = val + GetSiteEnergy(ix, iy, iz, Fields, e0ij) + GetSiteEnergyEwald(ix, iy, iz, Fields)
          end do
        end do
      end do
    else
      do iz = 1, cgrid%n3
        do iy = 1, cgrid%n2
          do ix = 1, cgrid%n1
            val = val + GetSiteEnergy(ix, iy, iz, Fields, e0ij)
          end do
        end do
      end do
    end if
    write(*,'(A20,I10,E25.15)') trim(OptAlgo)//":", opt_iter, val

  end subroutine OptEnergyFunc

end module optimization
