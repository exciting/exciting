
subroutine plot_spectral_function()
    use modmain,        only: pi, evalsv, efermi
    use mod_vxc,        only: vxcnn
    use mod_selfenergy, only: selfex, selfec, freq_selfc
    use modgw,          only: ibgw, nbgw, kset
    implicit none
    integer(4) :: ik, ib, iw, n
    real(8) :: w, sRe, sIm, div
    complex(8) :: dsc
    character(22) :: frmt
    real(8), allocatable :: sf(:)

    open(70,file='SpectralFunction.dat',form='FORMATTED',status='UNKNOWN',action='WRITE')
    open(71,file='Delta.dat',form='FORMATTED',status='UNKNOWN',action='WRITE')
    write(frmt, '("(",i8,"f14.6)")') 1 + nbgw-ibgw+1
    allocate(sf(ibgw:nbgw))
    do ik = 1, kset%nkpt
        write(70,*) '# ik = ', ik
        write(71,*) '# ik = ', ik
        do iw = 1, freq_selfc%nomeg
            w = freq_selfc%freqs(iw)
            ! compute spectral function
            do ib = ibgw, nbgw
                dsc = selfex(ib,ik) + selfec(ib,iw,ik) - vxcnn(ib,ik)
                sRe = dble(dsc)
                sIm = aimag(dsc)
                div = (w - (evalsv(ib,ik)-efermi) - sRe)**2 + sIm**2
                sf(ib) = 1.d0/pi * abs(sIm) / div
            end do
            write(70,trim(frmt)) w+efermi, sf(:)
            write(71,trim(frmt)) w, dble(w-evalsv(ibgw:nbgw,ik)-selfex(ibgw:nbgw,ik)+vxcnn(ibgw:nbgw,ik))
        end do
        write(70,*); write(70,*)
        write(71,*); write(71,*)
    end do
    deallocate(sf)
    close(70)
    close(71)

end subroutine