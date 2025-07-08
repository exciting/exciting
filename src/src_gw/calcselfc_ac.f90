
!> Analytical continuation of the correlation self-energy from the complex to the real frequency axis
!> No need to tackle degeneracies, as the function in imaginary axis has proper degeneracy
subroutine calcselfc_ac()

    use modinput, only: input, getstructwgrid, emptynode
    use modgw, only: kset, ibgw, nbgw
    use mod_selfenergy, only: selfec, freq_selfc
    use mod_frequency, only: generate_freqgrid, delete_freqgrid
    use mod_aaa_approximant, only: aaa_approximant, set_aaa_approximant, &
                                   init_aaa_approximant, delete_aaa_approximant, &
                                   get_aaa_approximant
    use mod_pade, only: pade_approximant
    use constants, only: real_zero
    use precision, only: i32, dp

    implicit none

    ! local variables
    type(aaa_approximant) :: aaa_minus, aaa_plus
    integer(i32) :: iw, ik, ib, n_kpoints
    real(dp)    :: w
    complex(dp) :: sc, dsc
    complex(dp), allocatable :: zj(:), fj(:,:,:)

    ! imaginary frequency grid
    n_kpoints = size( selfec, 3 )
    allocate( fj(ibgw:nbgw, freq_selfc%nomeg, n_kpoints) )
    fj(:,:,:) = selfec(:,:,:)
    deallocate(selfec)
    allocate(zj(freq_selfc%nomeg))
    do iw = 1, freq_selfc%nomeg
        zj(iw) = cmplx(real_zero, freq_selfc%freqs(iw), dp)
    end do
    call delete_freqgrid(freq_selfc)

    ! real frequency grid
    if ( .not.associated(input%gw%selfenergy%wgrid) ) &
        input%gw%selfenergy%wgrid => getstructwgrid(emptynode)
    call generate_freqgrid(freq_selfc, &
                           input%gw%selfenergy%wgrid%type, &
                           'refreq', &
                           input%gw%selfenergy%wgrid%size, &
                           input%gw%selfenergy%wgrid%wmin, &
                           input%gw%selfenergy%wgrid%wmax)
    allocate( selfec(ibgw:nbgw, freq_selfc%nomeg, n_kpoints) )

    do ik = 1, n_kpoints
        do ib = ibgw, nbgw
            
            if (input%gw%selfenergy%actype == 'pade' ) then

                do iw = 1, freq_selfc%nomeg
                    w = freq_selfc%freqs(iw)
                    if (w < real_zero) then
                        call pade_approximant(size(zj), -zj, conjg(fj(ib, :, ik)), cmplx(w, real_zero, dp), sc, dsc)
                    else
                        call pade_approximant(size(zj), zj, fj(ib, :, ik), cmplx(w, real_zero, dp), sc, dsc)
                    end if
                    selfec(ib,iw,ik) = sc
                end do

            else if (input%gw%selfenergy%actype == 'aaa' ) then

                ! No idea why only this way works ...
                ! It's probably related to proper choice of the contour (causality)
                call set_aaa_approximant(aaa_plus, -zj, fj(ib,:,ik), input%gw%selfenergy%tol)
                call init_aaa_approximant(aaa_minus, aaa_plus%nj, -aaa_plus%zj, &
                                          conjg(aaa_plus%fj), conjg(aaa_plus%wj))
                do iw = 1, freq_selfc%nomeg
                    w = freq_selfc%freqs(iw)
                    if (w < real_zero) then
                        sc = get_aaa_approximant(aaa_minus, cmplx(w, real_zero, dp))
                    else
                        sc = get_aaa_approximant(aaa_plus, cmplx(w, real_zero, dp))
                    end if
                    selfec(ib,iw,ik) = sc
                end do

            end if
            
        end do ! ie
    end do ! ik

    if (input%gw%selfenergy%actype == 'aaa') then
        call delete_aaa_approximant(aaa_plus)
        call delete_aaa_approximant(aaa_minus)
    end if
      
end subroutine
