
subroutine calcselfx(iq)
!
! Calculate the q-dependent self-energy contribution
!
    use modinput
    use modmain, only: nstfv, apwordmax, lmmaxapw, natmtot, &
                       pi, idxas, zzero, nmatmax, zone
    use modgw
    use mod_coulomb_potential
    use mod_mpi_gw, only : myrank
    use mod_hdf5
    implicit none

    ! input/output
    integer(4), intent(in) :: iq

    ! local
    integer(4) :: ik, ikp, jk
    integer(4) :: mdim, nmdim
    real(8)    :: tstart, tend, t0, t1
    integer(4) :: ie1, ie2, im
    integer(4) :: ia, is, ias, ic, icg
    real(8)    :: wkq, sxs2, fnk
    complex(8) :: sx, vc
    complex(8) :: mvm     ! Sum_ij{M^i*V^c_{ij}*conjg(M^j)}
    complex(8), allocatable :: evecfv(:,:)

    integer :: k, l, ispn, ist, jst
    complex(8) :: zsum

    ! external routines
    complex(8), external :: zdotc

    call timesec(tstart)

    ! singular term prefactor (q->0)
    sxs2 = 4.d0*pi*vi

    !----------------------------------------
    ! Set v-diagonal mixed product basis set
    !----------------------------------------
    if (vccut) then
        mbsiz = matsiz
        if (allocated(barc)) deallocate(barc)
        allocate(barc(matsiz,mbsiz))
        barc(:,:) = zzero
        do im = 1, matsiz
            if (barcev(im) > 0.d0) then
                vc = cmplx(barcev(im),0.d0,8)
                barc(:,im) = vmat(:,im)*sqrt(vc)
            end if
        end do
    else
        call setbarcev(0.d0)
    end if

    !--------------------------------------------------
    ! total number of states (n->m + n->c transisions)
    !--------------------------------------------------
    if ((input%gw%coreflag=='all').or. &
    &   (input%gw%coreflag=='xal')) then
      mdim = nomax+ncg
    else
      mdim = nomax
    end if
    nmdim = (nbgw-ibgw+1)*mdim

    allocate(eveckalm(nstfv,apwordmax,lmmaxapw,natmtot))
    allocate(eveckpalm(nstfv,apwordmax,lmmaxapw,natmtot))
    allocate(eveck(nmatmax,nstfv))
    allocate(eveckp(nmatmax,nstfv))

    allocate(minmmat(mbsiz,ibgw:nbgw,1:mdim))
    minmmat(:,:,:) = zzero
    ! msize = sizeof(minmmat)*b2mb
    ! write(*,'(" calcselfx: rank, size(minmmat) (Mb):",i4,f12.2)') myrank, msize

    !================================
    ! loop over irreducible k-points
    !================================
    ! write(*,*)
    do ikp = 1, kset%nkpt
      ! write(*,*) 'calcselfx: rank, (iq, ikp):', myrank, iq, ikp

      ! k vector
      ik = kset%ikp2ik(ikp)
      ! k-q vector
      jk = kqset%kqid(ik,iq)

      ! get KS eigenvectors
      allocate(evecfv(nmatmax,nstfv))
      call get_evec_gw(kqset%vkl(:,jk), Gkqset%vgkl(:,:,:,jk), evecfv)
      eveckp = conjg(evecfv)
      call get_evec_gw(kqset%vkl(:,ik), Gkqset%vgkl(:,:,:,ik), evecfv)
      eveck = evecfv
      deallocate(evecfv)

      call expand_evec(ik, 't')
      call expand_evec(jk, 'c')

      ! Calculate M^i_{nm}+M^i_{cm}
      call expand_products(ik, iq, ibgw, nbgw, -1, 1, mdim, nomax, minmmat)

      !========================================================
      ! Calculate the contribution to the exchange self-energy
      !========================================================
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ie1,ie2,mvm,icg,is,ia,ias,ic,fnk,sx)
!$OMP DO
#endif
      do ie1 = ibgw, nbgw

        ! sum over occupied states
        sx = zzero
        do ie2 = 1, mdim
          !=======================
          ! Valence contribution
          !=======================
          if (ie2 <= nomax) then
            mvm = zdotc(mbsiz, minmmat(:,ie1,ie2), 1, minmmat(:,ie1,ie2), 1)
            sx = sx - kiw(ie2,jk)*mvm
          else
            !=============================
            ! Core electron contribution
            !=============================
            icg = ie2-nomax
            is = corind(icg,1)
            ia = corind(icg,2)
            ias = idxas(ia,is)
            ic = corind(icg,3)
            mvm = zdotc(mbsiz, minmmat(:,ie1,ie2), 1, minmmat(:,ie1,ie2), 1)
            sx = sx - ciw(ic,ias)*mvm
          end if ! occupied states
        end do ! ie2

        ! add singular term (q->0)
        if ((Gamma) .and. (ie1 <= nomax)) then
          ! occupation number
          fnk = kiw(ie1,ik) * kqset%nkpt
          sx  = sx - sxs2*fnk*singc2
        end if

        selfex(ie1,ikp) = selfex(ie1,ikp) + sx

      end do ! ie1
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif

      ! debugging info
      if (input%gw%debug) then
        write(fdebug,*) 'EXCHANGE SELF-ENERGY: iq=', iq, ' ikp=', ikp
        write(fdebug,*) 'state   Sigma_x'
        do ie1 = ibgw, nbgw
          write(fdebug,*) ie1, selfex(ie1,ikp)
        end do
        write(fdebug,*)
      end if

    end do ! ikp

    deallocate(minmmat)
    deallocate(eveck)
    deallocate(eveckp)
    deallocate(eveckalm)
    deallocate(eveckpalm)

    ! timing
    call timesec(tend)
    time_selfx = time_selfx+tend-tstart

    return
end subroutine
!EOC
