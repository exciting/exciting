
subroutine calcselfc_ac()

    use modinput
    use modmain
    use modgw
    use mod_frequency
    use mod_vxc
    use mod_aaa_approximant
    use mod_pade
    implicit none

    ! local variables
    type(aaa_approximant) :: aaa_minus, aaa_plus
    integer(4) :: iw, ik, ib
    real(8)    :: w
    complex(8) :: sc, dsc
    complex(8), allocatable :: zj(:), fj(:,:,:)

    ! imaginary frequency grid
    allocate(fj(ibgw:nbgw,freq_selfc%nomeg,kset%nkpt))
    fj(:,:,:) = selfec(:,:,:)
    deallocate(selfec)
    allocate(zj(freq_selfc%nomeg))
    do iw = 1, freq_selfc%nomeg
        zj(iw) = cmplx(0.d0, freq_selfc%freqs(iw), 8)
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
    allocate(selfec(ibgw:nbgw,freq_selfc%nomeg,kset%nkpt))

    do ik = 1, kset%nkpt
        do ib = ibgw, nbgw
            
            if (input%gw%selfenergy%actype == 'pade' ) then

                do iw = 1, freq_selfc%nomeg
                    w = freq_selfc%freqs(iw)
                    if (w < 0.d0) then
                        call pade_approximant(size(zj), -zj, conjg(fj(ib,:,ik)), cmplx(w,0.d0,8), sc, dsc)
                    else
                        call pade_approximant(size(zj), zj, fj(ib,:,ik), cmplx(w,0.d0,8), sc, dsc)
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
                    if (w < 0.d0) then
                        sc = get_aaa_approximant(aaa_minus, cmplx(w,0.d0,8))
                    else
                        sc = get_aaa_approximant(aaa_plus, cmplx(w,0.d0,8))
                    end if
                    selfec(ib,iw,ik) = sc
                end do

            end if
            
        end do ! ie
    end do ! ik

    deallocate(zj, fj)
    if (input%gw%selfenergy%actype == 'aaa') then
        call delete_aaa_approximant(aaa_plus)
        call delete_aaa_approximant(aaa_minus)
    end if
      
    return
end subroutine
