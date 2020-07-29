
subroutine qdepwsum(iq,iomstart,iomend,ndim)

    use modinput
    use modmain
    use modgw

    implicit none
    integer(4), intent(in) :: iq
    integer(4), intent(in) :: iomstart, iomend
    integer(4), intent(in) :: ndim

    integer(4) :: iom, n, m
    integer(4) :: ik, jk, ikp, jkp
    integer(4) :: ia, is, ias, ic, icg
    real(8)    :: de, wkp, ene, occ, eta, ff
    complex(8) :: z1, z2, sfact
    real(8)    :: tstart, tend
    complex(8), allocatable :: om(:)

    call timesec(tstart)

    ! spin degeneracy
    sfact = 2.d0

    if (allocated(fnm)) deallocate(fnm)
    allocate(fnm(1:ndim,numin:nstdf,iomstart:iomend,1:kqset%nkpt))
    fnm(:,:,:,:) = zzero

    if (allocated(om)) deallocate(om)
    allocate(om(iomstart:iomend))

    select case (freq%fconv)
      case('refreq')
        om(iomstart:iomend) = freq%freqs(iomstart:iomend)
        eta = input%gw%freqgrid%eta
      case('imfreq')
        om(iomstart:iomend) = zi*freq%freqs(iomstart:iomend)
        eta = 0.d0
      case default
        stop 'Not supported option!'
    end select

    wkp = 1.0d0 / dble(kqset%nkpt)

    do ik = 1, kqset%nkpt
      jk  = kqset%kqid(ik,iq)
      ikp = kset%ik2ikp(ik)
      jkp = kset%ik2ikp(jk)

      do n = 1, ndim

        if (n <= nomax) then
          ene = evalfv(n,ikp)
          occ = occfv(n,ikp)/2.d0
        else
          icg = n - nomax
          is  = corind(icg,1)
          ia  = corind(icg,2)
          ic  = corind(icg,3)
          ias = idxas(ia,is)
          ene = evalcr(ic,ias)
          occ = 1.d0
        end if

        do m = numin, nstdf

          do iom = iomstart, iomend
            ff = occ * ( 1.d0 - occfv(m,jkp)/2.d0 )
            de = evalfv(m,jkp) - ene
            z1 = om(iom) - de + zi*eta
            z2 = om(iom) + de - zi*eta
            fnm(n,m,iom,ik) = sfact * ff * (1.d0/z1 - 1.d0/z2) * wkp
          end do ! iom

        end do ! m

      end do ! n

    end do ! ik

    deallocate(om)

    if (.false.) then
      iom = 1
      do ik = 1, kqset%nkpt
        write(*,*) 'iq, ik = ', iq, ik
        do n = 1, nomax
        do m = numin, nstdf
          write(*,'(2i4,2f12.6)') n, m, fnm(n,m,iom,ik)
        end do
        end do
      end do
    end if

    call timesec(tend)
    time_bzinit = time_bzinit+tend-tstart

end subroutine
