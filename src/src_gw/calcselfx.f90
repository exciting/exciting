!BOP
!!ROUTINE: calcselfx
!
!!INTERFACE: 
!
subroutine calcselfx(iq)
!      
!!DESCRIPTION:
!
! This subroutine calculates the q-dependent self energy contribution
!
!!USES:
    use modinput
    use modmain,    only : nstfv, apwordmax, lmmaxapw, natmtot, nspnfv, &
    &                      pi, idxas, zzero, nmatmax, zone, occsv
    use modgw
    use mod_mpi_gw, only : myrank
    use mod_hdf5
      
!!INPUT PARAMETERS:
    implicit none
    integer(4), intent(in) :: iq
    
!!LOCAL VARIABLES:            
    integer(4) :: ik, ikp, jk
    integer(4) :: mdim, nmdim
    real(8)    :: tstart, tend, t0, t1
    integer(4) :: ie1, ie2, im
    integer(4) :: ia, is, ias, ic, icg
    real(8)    :: wkq, sxs2
    complex(8) :: sx, vc
    complex(8) :: mvm     ! Sum_ij{M^i*V^c_{ij}*conjg(M^j)}
    complex(8), allocatable :: evecsv(:,:,:)
    
    integer :: k, l, ispn, ist, jst
    complex(8) :: zsum
    
    ! external routines 
    complex(8), external :: zdotc    
    
!!REVISION HISTORY:
!
! Created Nov 2013 by DIN
!
!EOP
!BOC
    ! if (myrank==0) then
    !   write(*,*)
    !   write(*,*) ' ---- calcselfx started ----'
    !   write(*,*)
    ! end if
    call timesec(tstart)
   
    if (vccut) then
      sxs2 = 0.d0
      !----------------------------------------
      ! Set v-diagonal mixed product basis set
      !----------------------------------------
      mbsiz = matsiz
      if (allocated(barc)) deallocate(barc)
      allocate(barc(matsiz,mbsiz))
      barc(:,:) = zzero
      do im = 1, matsiz
        vc = cmplx(barcev(im),0.d0,8)
        barc(:,im) = vmat(:,im)*sqrt(vc)
      end do
    else
      ! singular term prefactor (q->0)
      sxs2 = -4.d0*pi*vi
      !----------------------------------------
      ! Set v-diagonal mixed product basis set
      !----------------------------------------
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
    
    allocate(eveckalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveckpalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveck(nmatmax,nstsv))
    allocate(eveckp(nmatmax,nstsv))
    
    allocate(minmmat(mbsiz,ibgw:nbgw,1:mdim))
    minmmat(:,:,:) = zzero
    msize = sizeof(minmmat)*b2mb
    ! write(*,'(" calcselfx: rank, size(minmmat) (Mb):",i4,f12.2)') myrank, msize
    
    !================================
    ! loop over irreducible k-points
    !================================
    do ispn = 1, nspinor
    do ikp = 1, kset%nkpt
    
      ! write(*,*)
      ! write(*,*) '(calcselfx): k-point loop ikp=', ikp
    
      ! k vector
      ik = kset%ikp2ik(ikp)
      ! k-q vector 
      jk = kqset%kqid(ik,iq)
      
      ! get KS eigenvectors
      call timesec(t0)
      allocate(evecsv(nmatmax,nstsv,nspinor))
      call getevecsvgw_new('GW_EVECSV.OUT',jk,kqset%vkl(:,jk),nmatmax,nstsv,nspinor,evecsv)
      !call getevecfv(kqset%vkl(:,jk),Gkset%vgkl(:,:,:,jk),evecsv)
      eveckp = conjg(evecsv(:,:,ispn))
      
      call getevecsvgw_new('GW_EVECSV.OUT',ik,kqset%vkl(:,ik),nmatmax,nstsv,nspinor,evecsv)
      !call getevecfv(kqset%vkl(:,ik),Gkset%vgkl(:,:,:,ik),evecsv)
      eveck = evecsv(:,:,ispn)
      deallocate(evecsv)
      call timesec(t1)
      time_io = time_io+t1-t0
      
      ! Calculate M^i_{nm}+M^i_{cm}
      call expand_evec(ik,'t')
      call expand_evec(jk,'c')
      call expand_products(ik,iq,ibgw,nbgw,-1,1,mdim,nomax,minmmat)
        
      !========================================================
      ! Calculate the contribution to the exchange self-energy
      !========================================================
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ie1,sx,ie2,mvm,icg,ia,is,ic,ias)
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
            mvm = zdotc(mbsiz,minmmat(:,ie1,ie2),1,minmmat(:,ie1,ie2),1)
            sx = sx-kiw(ie2,jk)*mvm
          else
            !============================= 
            ! Core electron contribution
            !============================= 
            icg = ie2-nomax
            is = corind(icg,1)
            ia = corind(icg,2)
            ias = idxas(ia,is)
            ic = corind(icg,3)
            mvm = zdotc(mbsiz,minmmat(:,ie1,ie2),1,minmmat(:,ie1,ie2),1)
            sx = sx-ciw(ic,ias)*mvm
          end if ! occupied states
        end do ! ie2
        
        ! add singular term (q->0)
        if (Gamma.and.(ie1<=nomax)) sx = sx+sxs2*singc2*kiw(ie1,ik)*kqset%nkpt
        !if (Gamma.and.(dabs(kiw(ie1,ik))>1.d-6)) sx = sx+sxs2*singc2*kiw(ie1,ik)*kqset%nkpt
        
        selfex(ie1,ikp) = selfex(ie1,ikp)+sx
        
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
    end do ! ispn

    deallocate(minmmat)
    deallocate(eveck)
    deallocate(eveckp)
    deallocate(eveckalm)
    deallocate(eveckpalm)

    ! timing
    call timesec(tend)
    time_selfx = time_selfx+tend-tstart
    
    ! if (myrank==0) then
    !   write(*,*)
    !   write(*,*) ' ---- calcselfx ended ----'
    !   write(*,*)
    ! end if
    
    return
end subroutine
!EOC
