
subroutine qdepwsum(iq,iomstart,iomend,ndim)
    use constants, only: zzero, zi
    use mod_atoms, only: idxas
    use mod_bands, only: numin, nomax, nstdf, evalfv, occfv
    use mod_core_states, only: corind
    use mod_corestate, only: evalcr
    use modinput, only: input
    use modgw, only: fnm, kqset, freq, kqset, kset, time_bzinit
    use precision, only: dp, i32

    implicit none

    integer(i32), intent(in) :: iq
    integer(i32), intent(in) :: iomstart, iomend
    integer(i32), intent(in) :: ndim

    integer(i32) :: iom, n, m
    integer(i32) :: ik, jk, ikp, jkp
    integer(i32) :: ia, is, ias, ic, icg
    real(dp)     :: de, wkp, ene, occ, eta, ff
    complex(dp)  :: z1, z2, sfact
    real(dp)     :: tstart, tend
    complex(dp), allocatable :: om(:)

    call timesec(tstart)

    ! spin degeneracy
    sfact = 2.0_dp

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

    wkp = 1.0_dp / dble(kqset%nkpt)

    do ik = 1, kqset%nkpt
      jk  = kqset%kqid(ik,iq)
      ikp = kset%ik2ikp(ik)
      jkp = kset%ik2ikp(jk)

      do n = 1, ndim

        if (n <= nomax) then
          ene = evalfv(n,ikp)
          occ = occfv(n,ikp)/2.0_dp
        else
          icg = n - nomax
          is  = corind(icg,1)
          ia  = corind(icg,2)
          ic  = corind(icg,3)
          ias = idxas(ia,is)
          ene = evalcr(ic,ias)
          occ = 1.0_dp
        end if

        do m = numin, nstdf

          do iom = iomstart, iomend
            ff = occ * ( 1.0_dp - occfv(m,jkp)/2.0_dp )
            de = evalfv(m,jkp) - ene
            z1 = om(iom) - de + zi*eta
            z2 = om(iom) + de - zi*eta
            fnm(n,m,iom,ik) = sfact * ff * (1.0_dp/z1 - 1.0_dp/z2) * wkp
          end do ! iom

        end do ! m

      end do ! n

    end do ! ik

    deallocate(om)

    call timesec(tend)
    time_bzinit = time_bzinit+tend-tstart

end subroutine
