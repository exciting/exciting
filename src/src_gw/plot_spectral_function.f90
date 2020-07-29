
subroutine plot_spectral_function()
    use modinput
    use modmain,        only: pi, efermi
    use mod_vxc,        only: vxcnn
    use mod_selfenergy, only: selfex, selfec, freq_selfc, deltaE
    use modgw,          only: ibgw, nbgw, kset, evalfv
    implicit none
    integer(4) :: ik, ib, iw, n
    real(8) :: w, sRe, sIm, div
    complex(8) :: sc, dsc, sxc
    character(22) :: frmt
    real(8), allocatable :: sf(:)

    open(70,file='SpectralFunction.dat',form='FORMATTED',status='UNKNOWN',action='WRITE')
    open(71,file='Delta.dat',form='FORMATTED',status='UNKNOWN',action='WRITE')
    write(frmt, '("(",i8,"f14.6)")') 1 + nbgw-ibgw+1
    allocate(sf(ibgw:nbgw))
    do ik = 1, kset%nkpt
        write(70,'(a,i4,a,3f12.6,a,3f12.6,a,f12.6)') '# k-point: ik=', ik, '    vkl=', kset%vkl(:,ik), '    vkc=', kset%vkc(:,ik), '    wkpt=', kset%wkpt(ik) 
        write(71,'(a,i4,a,3f12.6,a,3f12.6,a,f12.6)') '# k-point: ik=', ik, '    vkl=', kset%vkl(:,ik), '    vkc=', kset%vkc(:,ik), '    wkpt=', kset%wkpt(ik) 
        do iw = 1, freq_selfc%nomeg
            w = freq_selfc%freqs(iw)
            ! compute spectral function
            do ib = ibgw, nbgw
                ! \Sigma_c(omega-\delta_E)
                call get_selfc(freq_selfc%nomeg, freq_selfc%freqs, selfec(ib,:,ik), &
                               w-deltaE, sc, dsc)
                sxc = selfex(ib,ik) + sc - vxcnn(ib,ik)
                sRe = dble(sxc)
                sIm = aimag(sxc) + input%gw%selfenergy%swidth
                div = (w-evalfv(ib,ik)-sRe)**2 + sIm**2
                sf(ib) = 1.d0/pi * abs(sIm) / div
            end do
            write(70,trim(frmt)) w, sf(ibgw:nbgw)
            write(71,trim(frmt)) w, dble(w-evalfv(ibgw:nbgw,ik)-selfex(ibgw:nbgw,ik)+vxcnn(ibgw:nbgw,ik))
        end do
        write(70,*); write(70,*)
        write(71,*); write(71,*)
    end do
    deallocate(sf)
    close(70)
    close(71)

end subroutine