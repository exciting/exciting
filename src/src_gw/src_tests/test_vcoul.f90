      
subroutine test_vcoul

    use modinput
    use modmain,               only : zzero, evalsv, efermi
    use modgw
    use mod_mpi_gw
    use m_getunit
    use mod_hdf5
            
    implicit none
    integer(4) :: iq, ik, jk, ispn
    integer(4) :: i, j, mdim, nmdim, im
    complex(8) :: vc
    complex(8), allocatable :: evecsv(:,:,:)

    complex(8), external :: zdotc

    !===========================================================================
    ! Initialization
    !===========================================================================
    
    ! prepare GW global data
    call init_gw
    
    ! clean not used anymore global exciting variables
    call clean_gndstate

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
    ! call setbarcev(input%gw%barecoul%barcevtol)
    ! call delete_coulomb_potential
    mbsiz = matsiz
    if (allocated(barc)) deallocate(barc)
    allocate(barc(matsiz,mbsiz))
    barc(:,:) = zzero
    do im = 1, matsiz
      vc = cmplx(barcev(im),0.d0,8)
      barc(:,im) = vmat(:,im)*sqrt(vc)
    end do

    !========================================
    ! Calculate M^I_{ij} coefficients
    !======================================== 
    ! arrays to store products of KS eigenvectors with the matching coefficients
    allocate(eveckalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveckpalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveck(nmatmax,nstsv))
    allocate(eveckp(nmatmax,nstsv))
    allocate(evecsv(nmatmax,nstsv,nspinor))
    call getevecsvgw_new('GW_EVECSV.OUT',jk,kqset%vkl(:,jk),nmatmax,nstsv,nspinor,evecsv)
    eveckp = conjg(evecsv(:,:,ispn))
    call getevecsvgw_new('GW_EVECSV.OUT',ik,kqset%vkl(:,ik),nmatmax,nstsv,nspinor,evecsv)
    eveck = evecsv(:,:,ispn)
    deallocate(evecsv)
    call expand_evec(ik,'t')
    call expand_evec(jk,'c')
    allocate(minmmat(mbsiz,ibgw:nbgw,1:mdim))
    call expand_products(ik,iq,ibgw,nbgw,-1,1,mdim,nomax,minmmat)
  
    !=================================================
    ! Calculate \sum_{IJ} (M^I_{ii})* v_{IJ} M^J_{jj}
    !=================================================

    do i = 1, 4
    do j = 1, 4
      vc = zdotc(mbsiz,minmmat(:,i,i),1,minmmat(:,j,j),1)
      write(*,'(2i,f16.8)') i, j, dble(vc)
    enddo
    print*, ""
    enddo


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
