
subroutine calcselfcSR(iq,iomstart,iomend)

    use modmain
    use modgw
    use mod_mpi_gw, only : myrank
    implicit none
    integer, intent(in) :: iq
    integer, intent(in) :: iomstart, iomend
    
    integer(4) :: i, ipw, ipw0, npw
    real(8) :: lambda, vc
    real(8) :: gpq(3), gpq2
    complex(8), allocatable :: tmat(:,:)
    
    integer :: lwork, info
    complex(8), allocatable :: work(:)
    real(8),    allocatable :: rwork(:)

    real(8) :: beta
    
    integer(4) :: ik, ikp, jk, ispn
    integer(4) :: ie, iom
    integer(4) :: mdim
    integer(4) :: fid
    character(120) :: fname_mwm
    real(8) :: tstart, tend, t0, t1
    complex(8), allocatable :: evecsv(:,:,:)
    character(len=10), external :: int2str
    
    if (myrank==0) then
      write(*,*)
      write(*,*) ' ---- calcselfcSR started ----'
      write(*,*)
    end if
    call timesec(tstart)
    
    ! range separation parameter
    lambda = 0.11d0 
    
    !----------------------
    ! PW-MB overlap matrix
    !----------------------
    call calcmpwmix(iq)

    !================================================
    ! Calculate the short-ranged Coulomb potential
    !================================================
    mbsiz = matsiz
    if (allocated(barc)) deallocate(barc)
    allocate(barc(matsiz,matsiz))
    barc(:,:) = 0.d0

    npw = Gqbarc%ngk(1,iq)
    allocate(tmat(matsiz,npw))
    tmat(:,:) = zzero
    
    ipw0 = 1
    if (Gamma) then
      ipw0 = 2
      beta = (6.d0*pi/(omega*kqset%nkpt))**(1.d0/3.d0)
      i_sz = 16.d0*pi*pi*kqset%nkpt/omega*(beta-lambda*dsqrt(pi)*derf(beta/(2.d0*lambda)))
      tmat(:,1) = i_sz*mpwmix(:,1)
    end if
    
    !------------------
    ! Loop over G+q
    !------------------
    do ipw = ipw0, npw
      gpq(1:3) = Gset%vgc(1:3,Gqbarc%igkig(ipw,1,iq))+kqset%vqc(1:3,iq)
      gpq2 = gpq(1)*gpq(1)+gpq(2)*gpq(2)+gpq(3)*gpq(3)
      vc = 4.d0*pi/gpq2*(1.d0-dexp(-gpq2/(4.d0*lambda*lambda)))
      tmat(:,ipw) = vc*mpwmix(:,ipw)
    end do ! ipw
    
    ! Coulomb potential in MB representation
    call zgemm('n','c',matsiz,matsiz,npw, &
    &          zone, &
    &          tmat,matsiz, &
    &          mpwmix,matsiz, &
    &          zzero, &
    &          barc,matsiz)
    deallocate(tmat)
    
    lwork = 2*matsiz
    allocate(work(lwork),rwork(3*matsiz))
    call zheev( 'v','u',matsiz,barc,matsiz, &
    &           barcev,work,lwork,rwork,info)
    call errmsg(info.ne.0,'CALCBARCMB',"Fail to diag. barc by zheev !!!")
    deallocate(work,rwork)
    
    !========================================
    ! Calculate the linear response function
    !========================================
    call init_dielectric_function(mbsiz,iomstart,iomend,Gamma)
    call calcepsilon(iq,iomstart,iomend)
    
    !======================================
    ! Calculate the correlation selfenergy
    !======================================
    
    ! total number of active states
    if (input%gw%coreflag=='all') then
      mdim = nstse+ncg
    else
      mdim = nstse
    end if
    
    !-------------------------------------------
    ! products M*W^c*M
    !-------------------------------------------
    allocate(mwm(ibgw:nbgw,1:mdim,1:freq%nomeg))
    msize = sizeof(mwm)*b2mb
    write(*,'(" calcselfcSR: size(mwm) (Mb):",f12.2)') msize
    
    allocate(eveckalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveckpalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveck(nmatmax,nstsv))
    allocate(eveckp(nmatmax,nstsv))
    allocate(minmmat(1:mbsiz,ibgw:nbgw,1:mdim))
    msize = sizeof(minmmat)*b2mb
    write(*,'(" calcselfcSR: rank, size(minmmat) (Mb):",i4,f12.2)') myrank, msize
    
    !================================
    ! loop over irreducible k-points
    !================================
    do ispn = 1, nspinor
    do ikp = 1, kset%nkpt
    
      write(*,*)
      write(*,*) 'calcselfcSR: k-point loop ikp=', ikp
    
      ! k vector
      ik = kset%ikp2ik(ikp)
      ! k-q vector 
      jk = kqset%kqid(ik,iq)
      
      ! get KS eigenvectors
      call timesec(t0)
      allocate(evecsv(nmatmax,nstsv,nspinor))
      call getevecsvgw_new('GW_EVECSV.OUT',jk,kqset%vkl(:,jk),nmatmax,nstsv,nspinor,evecsv)
      eveckp = conjg(evecsv(:,:,ispn))
      call getevecsvgw_new('GW_EVECSV.OUT',ik,kqset%vkl(:,ik),nmatmax,nstsv,nspinor,evecsv)
      eveck = evecsv(:,:,ispn)
      deallocate(evecsv)
      call timesec(t1)
      time_io = time_io+t1-t0
        
      ! Calculate M^i_{nm}+M^i_{cm}
      call expand_evec(ik,'t')
      call expand_evec(jk,'c')
      call expand_products(ik,iq,ibgw,nbgw,-1,1,mdim,nstse,minmmat)

      call calcselfcSR_k(ikp,iq,mdim)
      
    end do ! ikp
    end do ! ispn
    
    deallocate(minmmat)
    deallocate(eveck)
    deallocate(eveckp)
    deallocate(eveckalm)
    deallocate(eveckpalm)
    deallocate(mwm)
    
    call delete_dielectric_function(Gamma)
    if (allocated(kcw)) deallocate(kcw)
    if (allocated(unw)) deallocate(unw)
    
    ! timing
    call timesec(tend)
    time_selfc = time_selfc+tend-tstart
    
    if (myrank==0) then
      write(*,*)
      write(*,*) ' ---- calcselfcSR ended ----'
      write(*,*)
    end if
    
contains

    subroutine calcselfcSR_k(ikp,iq,mdim)
        implicit none
        integer, intent(in) :: ikp
        integer, intent(in) :: iq
        integer, intent(in) :: mdim
    
        integer :: iom
        integer :: ie1, ie2, nmdim
        integer :: ik, jk, jkp
        integer :: ia, is, ias, ic, icg
        real(8) :: wkq, enk 
        complex(8) :: sc
        complex(8) :: xnm(1:freq%nomeg)
        complex(8) :: wm(mbsiz,ibgw:nbgw,1:mdim)
    
        complex(8), external :: freqconv
        complex(8), external :: zdotu, zdotc
        external zhemm
    
        wkq = 1.d0/dble(kqset%nkpt)
        
        do iom = 1, freq%nomeg
          ! calculate \sum_{j} W^c_{ij}M^j_{nm}
          nmdim = (nbgw-ibgw+1)*mdim
          call zhemm('l','u',mbsiz,nmdim, &
          &          zone,vPv(:,:,iom),mbsiz,minmmat,mbsiz, &
          &          zzero,wm,mbsiz)
          
          do ie2 = 1, mdim
            do ie1 = ibgw, nbgw
              mwm(ie1,ie2,iom) = wkq*zdotc(mbsiz,minmmat(:,ie1,ie2),1,wm(:,ie1,ie2),1)
            end do
          end do
        end do ! iom
    
        !===========================
        ! Frequency convolution
        !===========================
        ! k-point
        ik = kset%ikp2ik(ikp)
        ! k-q point
        jk = kqset%kqid(ik,iq)
        jkp = kset%ik2ikp(jk)

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(iom,ie1,sc,ie2,enk,xnm,icg,is,ia,ic,ias)
!$OMP DO
#endif    
        do iom = 1, freq%nomeg
          ! loop over states
          do ie1 = ibgw, nbgw
      
            ! sum over states
            sc = zzero
            do ie2 = 1, mdim
        
              if (ie2<=nstsv) then
                !============================= 
                ! Valence electron contribution
                !============================= 
                enk = evalsv(ie2,jkp)-efermi
                xnm(1:freq%nomeg) = mwm(ie1,ie2,1:freq%nomeg)
              else
                !============================= 
                ! Core electron contribution
                !=============================
                icg = ie2-nstsv
                is = corind(icg,1)
                ia = corind(icg,2)
                ic = corind(icg,3)
                ias = idxas(ia,is)
                enk = evalcr(ic,ias)-efermi
                xnm(1:freq%nomeg) = mwm(ie1,ie2,1:freq%nomeg)
              end if ! val/cor
          
              !------------------------
              ! Frequency convolution
              !------------------------
              sc = sc+freqconv(iom,freq%nomeg,freq%freqs(iom), &
              &                enk,xnm,freq%freqs,freq%womeg)
                  
            end do ! ie2
    
            selfecSR(ie1,iom,ikp) = selfecSR(ie1,iom,ikp)+sc

          end do ! ie1
        end do ! iom
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif        
    
        return
    end subroutine

end subroutine
