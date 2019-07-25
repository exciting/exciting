subroutine calchead(ik, iomstart, iomend,ndim)
    !
    ! This subroutine calculate the head of the dielectric matrix
    !
    use modinput
    use modmain, only : zzero, zone, pi, evalcr, idxas
    use modgw
    implicit none

    ! input/output
    integer(4), intent(in) :: ik
    integer(4), intent(in) :: iomstart, iomend
    integer(4), intent(in) :: ndim

    ! local
    integer(4) :: ic, icg
    integer(4) :: ia, is, ias
    integer(4) :: ie1, ie2
    integer(4) :: ikp, iop, jop
    integer(4) :: iom
    real(8) :: edif    ! energy difference
    real(8) :: tstart, tend
    complex(8) :: coefh
    complex(8) :: pnm, zsum

    if (input%gw%debug) write(fdebug,*) ' ---- calchead started ----'
    call timesec(tstart)

    ! position in the non-reducied grid
    ikp = kset%ik2ikp(ik)

    ! constant prefactor
    coefh = cmplx(4.d0*pi*vi, 0.d0, 8)

    ! loop over tensor components
    do jop = 1, 3
    do iop = 1, 3

        ! loop over frequencies
        do iom = iomstart, iomend

            ! Inter-band contribution
            zsum = zzero
            do ie2 = numin, nstdf
            do ie1 = 1, ndim
                if (ie1 <= nomax) then
                    ! valence-valence
                    edif = evalfv(ie2,ikp)-evalfv(ie1,ikp)
                    if (dabs(edif) > 1.d-6) then
                        pnm = pmatvv(ie1,ie2,iop)*conjg(pmatvv(ie1,ie2,jop))
                        zsum = zsum + fnm(ie1,ie2,iom,ik)*pnm/(edif*edif)
                    end if
                else
                    ! core-valence
                    icg = ie1-nomax
                    is = corind(icg,1)
                    ia = corind(icg,2)
                    ias = idxas(ia,is)
                    ic = corind(icg,3)
                    edif = evalfv(ie2,ikp)-evalcr(ic,ias)
                    if (dabs(edif) > 1.0d-6) then
                        pnm = pmatcv(icg,ie2,iop)*conjg(pmatcv(icg,ie2,jop))
                        zsum = zsum + fnm(ie1,ie2,iom,ik)*pnm/(edif*edif)
                    end if
                end if
            end do ! ie2
            end do ! ie1
            epsh(iop,jop,iom) = epsh(iop,jop,iom) + coefh*zsum

            !-------------------------
            ! Intra-band contribution
            !-------------------------
            if (metallic) then
                zsum = zzero
                do ie1 = numin, nomax
                    pnm = pmatvv(ie1,ie1,iop)*conjg(pmatvv(ie1,ie1,jop))
                    zsum = zsum + kwfer(ie1,ik)*pnm
                enddo
                ! for imaginary frequency, a negative sign is needed
                if (freq%fconv == 'imfreq') zsum = -zsum
                epsh(iop,jop,iom) = epsh(iop,jop,iom) + &
                                    coefh*zsum/(freq%freqs(iom)**2)
            end if ! metallic

        end do ! iom

    end do ! iop
    end do ! jop

    ! timing
    call timesec(tend)
    time_dfhead = time_dfhead+tend-tstart
    if (input%gw%debug) then
        write(fdebug,*) ' ---- calchead finished ----'
        write(fdebug,*) ' ik = ', ik
        do iom = iomstart, iomend
            write(fdebug,*) iom, epsh(1,1,iom)
        end do
    end if

    return
end subroutine
