      
subroutine test_vcoul

    use modinput
    use modmain,               only: zzero, evalsv, efermi
    use modgw
    use mod_coulomb_potential, only: barc, delete_coulomb_potential
    use mod_mpi_gw
    use m_getunit
    use mod_hdf5
            
    implicit none
    integer(4) :: iq, ik, jk, ispn, iom
    integer(4) :: i, j, mdim, nmdim, im
    complex(8) :: vc
    complex(8), allocatable :: evecsv(:,:,:)
    complex(8), allocatable :: tmat(:,:), tvec(:)

    complex(8), external :: zdotc

    !===========================================================================
    ! Initialization
    !===========================================================================
    
    ! prepare GW global data
    call init_gw
    
    ! clean not used anymore global exciting variables
    call clean_gndstate

    if (.not.input%gw%rpmat) then
      !========================================================
      ! calculate momentum matrix elements and store to a file
      !========================================================
      call calcpmatgw
    end if
    
    ! occupancy dependent BZ integration weights
    call kintw

    !--------------------------------------------------
    ! total number of states (n->m + n->c transisions)
    !--------------------------------------------------
    if ((input%gw%coreflag=='all').or. &
    &   (input%gw%coreflag=='xal')) then
      mdim = nstsv+ncg
    else
      mdim = nstsv
    end if
    nmdim = (nbgw-ibgw+1)*mdim
    
    ispn = 1
    ik = 1
    iq = 1
    jk = kqset%kqid(ik,iq)

    Gamma = gammapoint(kqset%vqc(:,iq))

    !========================================
    ! Calculate interstitial basis functions
    !========================================
    matsiz = locmatsiz+Gqset%ngk(1,iq)
    call diagsgi(iq)
    call calcmpwipw(iq)

    !======================================
    ! Calculate the bare Coulomb potential
    !======================================
    call calcbarcmb(iq)
    call setbarcev(input%gw%barecoul%barcevtol)
    call delete_coulomb_potential

    !===================================
    ! Calculate the dielectric function
    !===================================
    call init_dielectric_function(mbsiz,1,freq%nomeg,Gamma)
    call calcepsilon(iq,1,freq%nomeg)
    ! write(*,*) "epsilon = ", sum(epsilon(:,:,2))

    !========================================
    ! Calculate M^I_{ij} coefficients
    !======================================== 
    ! arrays to store products of KS eigenvectors with the matching coefficients
    allocate(eveckalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveckpalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveck(nmatmax,nstsv))
    allocate(eveckp(nmatmax,nstsv))
    allocate(evecsv(nmatmax,nstsv,nspinor))
    call getevecsvgw('GW_EVECSV.OUT',jk,kqset%vkl(:,jk),nmatmax,nstsv,nspinor,evecsv)
    eveckp = conjg(evecsv(:,:,ispn))
    call getevecsvgw('GW_EVECSV.OUT',ik,kqset%vkl(:,ik),nmatmax,nstsv,nspinor,evecsv)
    eveck = evecsv(:,:,ispn)
    deallocate(evecsv)
    call expand_evec(ik,'t')
    call expand_evec(jk,'c')
    allocate(minmmat(mbsiz,ibgw:nbgw,1:mdim))
    call expand_products(ik,iq,ibgw,nbgw,-1,1,mdim,nstsv,minmmat)
  
    !=================================================
    ! Calculate \sum_{IJ} (M^I_{ii})* v_{IJ} M^J_{jj}
    !=================================================
    allocate(tmat(ibgw:nbgw,ibgw:nbgw))
    tmat = 0.d0
    do i = ibgw, nbgw
      do j = ibgw, nbgw
        tmat(i,j) = zdotc(mbsiz,minmmat(:,i,i),1,minmmat(:,j,j),1)
      enddo
      write(8,'(15f16.8)') abs(tmat(i,:))
    enddo

    !------------------------------------------------
    ! \sum_{IJ} (M^I_{ii})* \epsilon_{IJ} M^J_{jj}
    !------------------------------------------------
    iom = 1
    
    do im = 1, mbsiz
      epsilon(im,im,iom) = epsilon(im,im,iom)-zone
    end do ! im

    tmat = 0.d0
    allocate(tvec(mbsiz))
    do i = ibgw, nbgw
      do j = ibgw, nbgw
        tvec = 0.d0
        call zgemv('n',mbsiz,mbsiz,zone,epsilon(:,:,iom),mbsiz, &
        &           minmmat(:,j,j),1,zzero,tvec,1)
        tmat(i,j) = zdotc(mbsiz,minmmat(:,i,i),1,tvec,1)
      enddo
      write(7,'(15f16.8)') abs(tmat(i,:))
    enddo
    deallocate(tvec,tmat)

    call delete_dielectric_function(Gamma)
    if (allocated(kcw)) deallocate(kcw)
    if (allocated(unw)) deallocate(unw)
    !------------------------------------------------
          
    
    ! clean unused data
    if (allocated(mpwipw)) deallocate(mpwipw)
    if (allocated(barc)) deallocate(barc)
    if (allocated(evalsv)) deallocate(evalsv)
    call delete_freqgrid(freq)
    call delete_k_vectors(kset)
    call delete_G_vectors(Gset)
    call delete_Gk_vectors(Gkset)
    call delete_kq_vectors(kqset)
    call delete_Gk_vectors(Gqset)
    call delete_Gk_vectors(Gqbarc)
    
    return
end subroutine
