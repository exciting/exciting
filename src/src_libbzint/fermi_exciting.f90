
subroutine fermi_exciting(lspin,nel,nb,nik,eband,ntet,tetc,wtet,vt,efer,eg,df)
    use order
    implicit none
    ! input parameters
    logical, intent(in)    :: lspin         ! flag for spin-polarized calculations
    real(8), intent(in)    :: nel           ! number of electrons
    integer(4), intent(in) :: nb            ! Maximum number of bands
    integer(4), intent(in) :: nik           ! Number of irreducible k-points
    real(8),    intent(in) :: eband(nb,nik) ! Band energies
    integer(4), intent(in) :: ntet          ! Number of tetrahedra
    integer(4), intent(in) :: tetc(4,*)     ! id. numbers of the corners of the tetrahedr
    integer(4), intent(in) :: wtet(*)       ! weight of each tetrahedron
    real(8), intent(in)    :: vt            ! the volume of the tetrahedra
    ! output variables
    real(8), intent(out)   :: efer          ! Efermi
    real(8), intent(out)   :: eg            ! band gap energy
    real(8), intent(out)   :: df            ! dos at Efermi
    ! local variables
    integer :: ik, ib, it, vbm(2), cbm(2)
    integer :: nvm
    real(8) :: evbm, ecbm, eint
    real(8) :: ocmin, ocmax, ocint, sfact
    integer :: nitmax = 200
    real(8) :: eps = 1.d-8
    real(8) :: eb(nb,nik)

    logical :: debug=.false.
    
    real(8), external :: idos_exciting
    real(8), external :: dostet_exciting
    
    !------------------
    ! Initialization
    !------------------
    if (lspin) then
      sfact = 1.d0
    else
      sfact = 2.d0
    end if
    
    ! Important to use sorted energies
    eb(:,:) = eband(:,:)
    do ik = 1, nik
      call sort(nb,eb(:,ik))
    end do
    
    ! nvm is the number of bands for an insulating system 
    ! since for a system with gap, the procedure to determine the
    ! band gap can be unstable, just try first whether it is an
    ! insulating system, but such a simplistic way to determine the Fermi energy
    ! is valid only for spin un-polarized cases 
    if (.not.lspin) then 
      nvm  = nint(nel/2.d0) 
      evbm = maxval(eb(nvm,:))
      ecbm = minval(eb(nvm+1,:))
      eint = (evbm+ecbm)/2.d0
      ocint = sfact*idos_exciting(nb,nik,eb,ntet,tetc,wtet,vt,eint)
      if ((ecbm>=evbm).and.(abs(ocint-nel)<eps)) then 
        efer = eint
        eg = ecbm-evbm
        return 
      end if 
    end if
    
    ! find the minimal and maximal band energy  
    evbm = minval(eb)
    ecbm = maxval(eb,mask=eb<1.d3)
    if (debug) write(6,*) 'evbm, ecbm = ', evbm, ecbm
    
    ocmin = sfact*idos_exciting(nb,nik,eb,ntet,tetc,wtet,vt,evbm)
    ocmax = sfact*idos_exciting(nb,nik,eb,ntet,tetc,wtet,vt,ecbm)
    if (debug) write(6,*) 'ocmin, ocmax = ', ocmin, ocmax
        
    if (ocmax<=nel) then
      write(6,'(a)')  'ERROR(libbzint::fermi_exciting): Not enough bands!'
      write(6,'(a,f10.4,2f10.2)') '  emax, ocmax, nel = ', ecbm, ocmax, nel 
      stop
    end if 

    if (debug) write(6,1) '#it',"evbm","ecbm","eint","ocmin","ocmax","ocint"

    !--------------------------------------------------------------------------
    ! Use bisection method to determine solver the equation N( Efermi ) = Nel
    !--------------------------------------------------------------------------
    do it = 1, nitmax
      eint = evbm+0.5d0*(ecbm-evbm)
      ocint = sfact*idos_exciting(nb,nik,eb,ntet,tetc,wtet,vt,eint)
      if (debug) write(6,2) it, evbm, ecbm, eint, ocmin, ocmax, ocint
      if (abs(ocint-nel)<eps) exit  
      if (ocint>nel) then
        ecbm = eint
        ocmax = ocint
      else
        evbm = eint
        ocmin = ocint
      end if
    end do
    if (it>=nitmax) then
      write(6,'(a)')  'ERROR(libbzint::fermi_exciting): Failed to converge!'
      stop
    end if

    df = sfact*dostet_exciting(nb,nik,eb,ntet,tetc,wtet,vt,eint)
    
    if (df<1.0d-4) then
      !-------------------------
      ! System with a band gap
      !-------------------------
      evbm = maxval(eb,mask=eb<eint)
      ecbm = minval(eb,mask=eb>eint)
      efer = 0.5d0*(evbm+ecbm)
      eg = ecbm-evbm
      if (debug) then
        write(6,*)
        write(6,*) 'Insulator:'
        write(6,*) '  E_VBM = ',  evbm
        write(6,*) '  E_CBM = ',  ecbm
        write(6,*) '  Eg = ', eg
        write(6,*) '  E_Fermi = ', efer
      end if

    else
      !------------------------
      ! Metal
      !------------------------
      efer = eint
      if (debug) then
        write(6,*) 'Metal:'
        write(6,*) '  DOS at E_Fermi = ', df
        write(6,*) '  E_Fermi = ', efer
      end if
      eg = 0.d0
    end if
    
    1 format(A5,6A12)
    2 format(I5,6f12.6)
end subroutine
      
