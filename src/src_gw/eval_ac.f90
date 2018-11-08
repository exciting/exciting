
subroutine eval_ac(ni, zi, fi, z, fz, dfz)
    use modinput
    use mod_aaa_approximant
    use mod_pade
    implicit none
    integer(4), intent(in)  :: ni
    complex(8), intent(in)  :: zi(ni)
    complex(8), intent(in)  :: fi(ni)
    complex(8), intent(in)  :: z
    complex(8), intent(out) :: fz
    complex(8), intent(out) :: dfz
    ! local
    type(aaa_approximant) :: aaa
    real(8)    :: tol, step
    complex(8) :: fz1, fz2

    select case( trim(input%gw%selfenergy%actype) )
        
        case('pade')
            if (dble(z) > 0.d0) then
                call pade_approximant(ni, zi, fi, z, fz, dfz)
            else
                call pade_approximant(ni, -zi, conjg(fi), z, fz, dfz)
            end if

        case('aaa')
            tol  = input%gw%selfenergy%tol
            if (dble(z) > 0.d0) then
                call set_aaa_approximant(aaa, -zi, fi, tol)
            else
                call set_aaa_approximant(aaa, zi, conjg(fi), tol)
            end if
            fz = get_aaa_approximant(aaa, z)
            ! Compute derivative dfz
            step = 0.01d0
            fz1 = get_aaa_approximant(aaa, z-step)
            fz2 = get_aaa_approximant(aaa, z+step)
            dfz = 0.5d0 * (fz2-fz1) / step
            call delete_aaa_approximant(aaa)

        case default
            stop 'Error(eval_ac): Unknown actype!'

    end select

end subroutine