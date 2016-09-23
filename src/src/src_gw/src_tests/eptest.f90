subroutine eptest(ik,jk,iq)

    use modmain
    use modgw

    implicit none
    integer(4), intent(in) :: ik     ! index of the k-point.
    integer(4), intent(in) :: jk     ! index of the k'-points.
    integer(4), intent(in) :: iq     ! index of the q-point
    integer :: ispn
    complex(8), allocatable :: evecsv(:,:,:)
    
    call diagsgi(iq)
    call calcmpwipw(iq)
    
    matsiz = locmatsiz+Gqset%ngk(1,iq)
    write(*,101) iq, locmatsiz, Gqset%ngk(1,iq), matsiz
    
    allocate(eveckalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveckpalm(nstsv,apwordmax,lmmaxapw,natmtot))
    allocate(eveck(nmatmax,nstsv))
    allocate(eveckp(nmatmax,nstsv))

    ! get KS eigenvectors
    allocate(evecsv(nmatmax,nstsv,nspinor))
    call getevecsvgw_new('GW_EVECSV.OUT',jk,kqset%vkl(:,jk),nmatmax,nstsv,nspinor,evecsv)
    eveckp = conjg(evecsv(:,:,1))
    call getevecsvgw_new('GW_EVECSV.OUT',ik,kqset%vkl(:,ik),nmatmax,nstsv,nspinor,evecsv)
    eveck = evecsv(:,:,1)
    deallocate(evecsv)

    call expand_evec(ik,'t')
    call expand_evec(jk,'c')
    
    allocate(minmmat(matsiz,nstfv,nstfv))
    call calcminm(ik,iq,1,nstfv,1,nstfv,minmmat)

    deallocate(eveck)
    deallocate(eveckp)
    deallocate(eveckalm)
    deallocate(eveckpalm)

101 format(/,'Data for q-point nr.:',i4,/,4x, &
    &        'Mixed basis:',/,4x, &
    &        'Number of atomic basis functions:       ',i4,/,4x, &
    &        'Number of interstitial basis functions: ',i4,/,4x, &
    &        'Total number of basis functions:        ',i4,/)

    return      
end subroutine

