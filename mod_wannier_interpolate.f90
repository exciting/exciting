module m_wannier_interpolate_density
  use modmain
  use mod_wannier
  use m_plotmat
  implicit none
  contains

subroutine wannier_interpolate_density( int_kset)
  use m_wannier_interpolate_eigsys
  use m_wannier_interpolate_fft

  type( k_set), intent( in) :: int_kset

  integer :: ik, ik2, iq, ist, jst, is, ia, ias, tp, ir, irc, igk, ifg, ig, iv(3)
  real(8) :: efermiint, wgt, rc(3), s(3,3), si(3,3), v1(3)
  complex(8) :: z1
  real(8), allocatable :: evalfv(:,:), eval1(:,:), evalint(:,:), occint(:,:)
  complex(8), allocatable :: evecint(:,:,:), w(:,:)
  complex(8), allocatable :: evecfv(:,:,:),  evecsv(:,:), apwalm(:,:,:,:,:), wfmtlm(:,:), wfmttp(:,:,:)
  complex(8), allocatable :: zfft(:,:), zfft2(:,:)
  real(8), allocatable :: phase(:,:), rmttp(:,:), rhomt_int(:,:,:), rhoir_int(:)
  real(8), allocatable :: rlfft(:,:)
  complex(8), allocatable :: wfmttpq(:,:,:,:), wfirq(:,:)
  type( k_set) :: tmp_kset
  real(8), allocatable :: rfft(:,:)
  integer, allocatable :: gfft(:,:)

  character( 256) :: fname
  
  allocate( rfft( 3, wf_Gset%ngrtot))
  allocate( gfft( 3, wf_Gset%ngrtot))
  !allocate( zfft( wf_Gset%ngrtot, wf_fst:wf_lst))
  !allocate( zfft2( wf_Gset%ngrtot, wf_fst:wf_lst))
  !allocate( evecfv( nmatmax, nstfv, nspnfv))
  do ig = 1, wf_Gset%ngrtot
    ifg = igfft( ig)
    rfft( :, ifg) = dble( wf_Gset%ivg( :, ig))/ngrid
    gfft( :, ifg) = wf_Gset%ivg( :, ig)
  end do

  !do ik = 1, wf_kset%nkpt
  !  call findkpt( wf_kset%vkl( :, ik), is, ia)
  !  if( ia .eq. 2) write(*,*) ik
  !end do
  !return
  !ik = 14
  !call findkpt( wf_kset%vkl( :, ik), is, ia)
  !if( ia .eq. 2) then
  !  s = dble( symlat( :, :, lsplsymc( is)))
  !  v1 = vtlsymc( :, is)
  !  write(*,'(3F13.6)') s(1,:)
  !  write(*,'(3F13.6)') s(2,:)
  !  write(*,'(3F13.6)') s(3,:)
  !  write(*,*)
  !  write(*,'(3F13.6)') v1
  !  write(*,*)
  !  ! read eigenvector      
  !  if( input%properties%wannier%input .eq. "groundstate") then
  !    call getevecfv( wf_kset%vkl( :, ik), wf_Gkset%vgkl( :, :, :, ik), evecfv)
  !  else if( input%properties%wannier%input .eq. "hybrid") then
  !    call getevecfv( wf_kset%vkl( :, ik), wf_Gkset%vgkl( :, :, :, ik), evecfv)
  !  else if( input%properties%wannier%input .eq. "gw") then
  !    call getevecsvgw_new( "GW_EVECSV.OUT", ik, wf_kset%vkl( :, ik), nmatmax, nstfv, nspnfv, evecfv)
  !  else
  !    call terminate
  !  end if
  !  zfft = zzero
  !  do ist = wf_fst, wf_lst
  !    do igk = 1, wf_Gkset%ngk( 1, ik)
  !      ifg = igfft( wf_Gkset%igkig( igk, 1, ik))
  !      zfft( ifg, ist) = evecfv( igk, ist, 1)
  !    end do
  !    call zfftifc( 3, ngrid, 1, zfft( :, ist))
  !  end do
  !  do ir = 1, wf_Gset%ngrtot
  !    write( *, '(3F13.6,3I,SP,F13.6)') rfft( :, ir), gfft( :, ir), dble( zfft( ir, 3))**2 + aimag( zfft( ir, 3))**2
  !  end do
  !end if
  !return
      

  allocate( eval1( wf_fst:wf_lst, wf_kset%nkpt))
  allocate( evalfv( nstfv, nspnfv))
  ! read FV eigenvalues
  do ik = 1, wf_kset%nkpt
    call getevalfv( wf_kset%vkl( :, ik), evalfv)
    eval1( :, ik) = evalfv( :, 1)
  end do

  allocate( evalint( wf_fst:wf_lst, int_kset%nkpt))
  allocate( evecint( wf_fst:wf_lst, wf_fst:wf_lst, int_kset%nkpt))
  allocate( phase( wf_kset%nkpt, int_kset%nkpt))
  ! interpolate eigenenergies
  call wannier_interpolate_eigsys( eval1( wf_fst:wf_lst, :), int_kset, evalint, evecint, phase)

  allocate( occint( wf_fst:wf_lst, int_kset%nkpt))
  ! interpolate occupancies and Fermi energy
  call wannier_interpolate_occupancy( int_kset, evalint, occint, efermiint)
  write(*,'(F23.16)') efermiint

  !return

  call readstate
  call readfermi
  call linengy
  call genapwfr
  call genlofr
  call olprad

  allocate( evecfv( nmatmax, nstfv, nspnfv))
  allocate( apwalm( ngkmax, apwordmax, lmmaxapw, natmtot, nspnfv))
  allocate( wfmttp( lmmaxvr, nrcmtmax, wf_fst:wf_lst))
  allocate( wfmtlm( lmmaxvr, nrcmtmax))
  allocate( zfft( wf_Gset%ngrtot, wf_fst:wf_lst))

  !********************************
  ! store wavefunctions in file
  !********************************

  do ik = 1, wf_kset%nkpt
    ! get matching coefficients
    call match( wf_Gkset%ngk( 1, ik), wf_Gkset%gkc( :, 1, ik), wf_Gkset%tpgkc( :, :, 1, ik), wf_Gkset%sfacgk( :, :, 1, ik), apwalm( :, :, :, :, 1))
      
    ! read eigenvector      
    if( input%properties%wannier%input .eq. "groundstate") then
      call getevecfv( wf_kset%vkl( :, ik), wf_Gkset%vgkl( :, :, :, ik), evecfv)
    else if( input%properties%wannier%input .eq. "hybrid") then
      call getevecfv( wf_kset%vkl( :, ik), wf_Gkset%vgkl( :, :, :, ik), evecfv)
    else if( input%properties%wannier%input .eq. "gw") then
      call getevecsvgw_new( "GW_EVECSV.OUT", ik, wf_kset%vkl( :, ik), nmatmax, nstfv, nspnfv, evecfv)
    else
      call terminate
    end if

    ! build muffin-tin wavefunctions
    do is = 1, nspecies
      do ia = 1, natoms( is)
        ias = idxas( ia, is)
        ! build muffin-tin wavefunctions
        ! and transform to spherical coordinates
        wfmttp = zzero
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ist, wfmtlm)
!$OMP DO SCHEDULE( guided) 
#endif
        do ist = wf_fst, wf_lst
          call wavefmt( input%groundstate%lradstep, input%groundstate%lmaxvr, is, ia, wf_Gkset%ngk( 1, ik), &
               apwalm, evecfv( :, ist, 1), lmmaxvr, wfmtlm( :, :))
          call zgemm( 'N', 'N', lmmaxvr, nrcmt( is), lmmaxvr, zone, &
               zbshtvr, lmmaxvr, &
               wfmtlm( :, 1:nrcmt( is)), lmmaxvr, zzero, &
               wfmttp( :, 1:nrcmt( is), ist), lmmaxvr)
        end do
#ifdef USEOMP                
!$OMP END DO  
!$OMP END PARALLEL
#endif
      end do
    end do

    ! write muffin-tin wavefunction to file
    call wannier_putwfmt( ik, wfmttp)

    ! build interstitial wavefunction
    zfft = zzero
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ist, igk, ifg)
!$OMP DO  
#endif
    do ist = wf_fst, wf_lst
      do igk = 1, wf_Gkset%ngk( 1, ik)
        ifg = igfft( wf_Gkset%igkig( igk, 1, ik))
        zfft( ifg, ist) = evecfv( igk, ist, 1)
      end do
      call zfftifc( 3, ngrid, 1, zfft( :, ist))
    end do
#ifdef USEOMP                
!$OMP END DO  
!$OMP END PARALLEL
#endif
    
    ! write interstitial wavefunction to file
    call wannier_putwfir( ik, zfft)

  end do
  
  deallocate( apwalm, evecfv, wfmtlm)
  write(*,*) "wavefunctions written"
    

  allocate( w( wf_fst:wf_lst, wf_fst:wf_lst))
  allocate( rmttp( lmmaxvr, nrcmtmax))
  ! interpolated wavefunction
  allocate( wfmttpq( lmmaxvr, nrcmtmax, wf_fst:wf_lst, natmtot))
  allocate( wfirq( wf_Gset%ngrtot, wf_fst:wf_lst))
  ! interpolated density
  allocate( rhomt_int( lmmaxvr, nrmtmax, natmtot))
  allocate( rhoir_int( wf_Gset%ngrtot))

  rhomt_int = 0.d0
  rhoir_int = 0.d0

  write(*,*) "start interpolation"
  write(*,*) ngrid
  write(*,*) wf_Gset%intgv(1,:)
  write(*,*) wf_Gset%intgv(2,:)
  write(*,*) wf_Gset%intgv(3,:)

  do iq = 1, int_kset%nkpt
    write(*,*) iq
    wfmttpq = zzero
    wfirq = zzero
    do ik = 1, wf_kset%nkpt

      ! build mixing matrix
      call zgemm( 'N', 'N', wf_nst, wf_nst, wf_nst, cmplx( phase( ik, iq), 0, 8), &
           wf_transform( :, :, ik), wf_nst, &
           evecint( :, :, iq), wf_nst, zzero, &
           w, wf_nst)

      call wannier_getwfmt( ik, wfmtlm)

      do is = 1, nspecies
        do ia = 1, natoms( is)
          ias = idxas( ia, is)

#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ir)
!$OMP DO SCHEDULE( guided) 
#endif
          do ir = 1, nrcmt( is)
            call zgemm( 'N', 'N', lmmaxvr, wf_nst, wf_nst, zone, &
                 wfmttp( :, ir, :), lmmaxvr, &
                 w, wf_nst, zone, &
                 wfmttpq( :, ir, :, ias), lmmaxvr)
          end do
#ifdef USEOMP                
!$OMP END DO  
!$OMP END PARALLEL
#endif

        end do
      end do

      call wannier_getwfir( ik, zfft)

      call zgemm( 'N', 'N', wf_Gset%ngrtot, wf_nst, wf_nst, zone, &
           zfft, wf_Gset%ngrtot, &
           w, wf_nst, zone, &
           wfirq, wf_Gset%ngrtot)

    end do

    !********************************
    ! sum up density contributions
    !********************************
    do is = 1, nspecies
      do ia = 1, natoms( is)
        ias = idxas( ia, is)

        do ist = wf_fst, wf_lst
          rmttp = 0.d0
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( irc, tp, z1)
!$OMP DO  
#endif
          do irc = 1, nrcmt( is)
            do tp = 1, lmmaxvr
              z1 = wfmttpq( tp, irc, ist, ias)
              rmttp( tp, irc) = rmttp( tp, irc) + int_kset%wkpt( iq)*occint( ist, iq)*(dble( z1)**2 + aimag( z1)**2)
            end do
          end do
#ifdef USEOMP                
!$OMP END DO  
!$OMP END PARALLEL
#endif
          
          irc = 0
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ir, irc, tp, z1)
!$OMP DO  
#endif
          do ir = 1, nrmt( is), input%groundstate%lradstep
            irc = 1 + (ir-1)/input%groundstate%lradstep
            call dgemv( 'N', lmmaxvr, lmmaxvr, 1.d0, &
                 rfshtvr, lmmaxvr, &
                 rmttp( :, irc), 1, 1.d0, &
                 rhomt_int( :, ir, ias), 1)
          end do
#ifdef USEOMP                
!$OMP END DO  
!$OMP END PARALLEL
#endif
        end do

      end do
    end do


#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( ist, ir, z1) reduction(+: rhoir_int)
!$OMP DO  
#endif
    do ist = wf_fst, wf_lst
      do ir = 1, wf_Gset%ngrtot
        z1 = wfirq( ir, ist)
        rhoir_int( ir) = rhoir_int( ir) + int_kset%wkpt( iq)*occint( ist, iq)*(dble( z1)**2 + aimag( z1)**2)/omega
      end do
    end do
#ifdef USEOMP                
!$OMP END DO  
!$OMP END PARALLEL
#endif

  end do

  deallocate( w, wfmttp, rmttp, zfft, wfirq, wfmttpq)

  call wannier_destroywf
  !stop

  rhomt = rhomt_int
  rhoir = rhoir_int

  call charge
  write(*,'(100F13.6)') chgmt
  write(*,'(F13.6)') chgir
  write(*,'(F13.6)') chgcalc
  write(*,'(F13.6)') chgtot
  
  write(*,*) "start calculation"
  rhomt = 0.d0
  rhoir = 0.d0

  call generate_k_vectors( tmp_kset, bvec, input%groundstate%ngridk, input%groundstate%vkloff, input%groundstate%reducek)

  allocate( evecsv( nstsv, nstsv))
  deallocate( evalint, evecint, phase, occint)
  allocate( evalint( wf_fst:wf_lst, tmp_kset%nkpt))
  allocate( evecint( wf_fst:wf_lst, wf_fst:wf_lst, tmp_kset%nkpt))
  allocate( phase( wf_kset%nkpt, tmp_kset%nkpt))
  allocate( occint( wf_fst:wf_lst, int_kset%nkpt))

  ! interpolate eigenenergies
  call wannier_interpolate_eigsys( eval1( wf_fst:wf_lst, :), tmp_kset, evalint, evecint, phase)
  ! interpolate occupancies and Fermi energy
  call wannier_interpolate_occupancy( tmp_kset, evalint, occint, efermiint)
  occsv = 0.d0
  occsv( 1:26, :) = 2.d0
  !occsv( wf_fst:wf_lst, :) = occint

  do ik = 1, nkpt
    call getevecfv( vkl( :, ik), vgkl( :, :, :, ik), evecfv)
    call getevecsv( vkl( :, ik), evecsv)
    call rhovalk( ik, evecfv, evecsv)
    call genrhoir( ik, evecfv, evecsv)
  end do

  call charge
  write(*,'(100F13.6)') chgmt
  write(*,'(F13.6)') chgir
  write(*,'(F13.6)') chgcalc
  write(*,'(F13.6)') chgtot

  !do is = 1, nspecies
  !  do ia = 1, natoms( is)
  !    ias = idxas( ia, is)
  !    do tp = 1, lmmaxvr
  !      do ir = 1, nrmt( is)
  !        write( *, '(3(I3.3,3x),SP,3F23.6)') ias, tp, ir, rhomt( tp, ir, ias), rhomt_int( tp, ir, ias), rhomt_int( tp, ir, ias)/rhomt( tp, ir, ias) - 1.d0
  !      end do
  !    end do
  !    write(*,*)
  !  end do
  !end do

  do ir = 1, wf_Gset%ngrtot
    write( *, '(3F13.6,3x,SP,3F23.6)') rfft( :, ir), rhoir( ir), rhoir_int( ir), rhoir_int( ir)/rhoir( ir) - 1.d0
  end do

  return
end subroutine wannier_interpolate_density

subroutine wannier_interpolate_occupancy( int_kset, evalint, occint, efermiint)
  use mod_eigenvalue_occupancy, only: occmax
  use mod_charge_and_moment, only: chgval

  type( k_set), intent( in) :: int_kset
  real(8), intent( in) :: evalint( wf_fst:wf_lst, int_kset%nkpt)
  real(8), intent( out) :: occint( wf_fst:wf_lst, int_kset%nkpt)
  real(8), intent( out) :: efermiint

  integer, parameter :: maxit = 1000

  integer :: iq, ist, nvm, it
  real(8) :: e0, e1, fermidos, chg, x, t1
  
  real(8) :: sdelta, stheta

  if( input%groundstate%stypenumber .ge. 0 ) then
    t1 = 1.d0/input%groundstate%swidth
    nvm = nint( chgval/occmax)
    if( (wf_fst .ne. 1) .or. (nvm .ge. wf_lst)) then
      write( *, '("Error (wanner_interpolate_occupancy): Please wannierize all occupied and at least one unoccupied band for the calculation of occupation-numbers.")')
      call terminate
    end if
    ! check for insulator or semiconductor
    e0 = maxval( evalint( nvm, :))
    e1 = minval( evalint( nvm+1, :))
    efermiint = 0.5*(e0 + e1)

    fermidos = 0.d0
    chg = 0.d0
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iq, ist, x) reduction(+: chg, fermidos)
!$OMP DO  
#endif
    do iq = 1, int_kset%nkpt
      do ist = wf_fst, wf_lst
        x = (evalint( ist, iq) - efermiint)*t1
        fermidos = fermidos + int_kset%wkpt( iq)*sdelta( input%groundstate%stypenumber, x)*t1
        occint( ist, iq) = occmax*stheta( input%groundstate%stypenumber, -x)
        chg = chg + int_kset%wkpt( iq)*occint( ist, iq)
      end do
    end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
    fermidos = fermidos+occmax
    if( (e1 .ge. e0) .and. (abs( chg - chgval) .lt. input%groundstate%epsocc)) then
    !  write( *, '("Info (wannier_interpolate_occupancy): System has gap. Simplistic method used in determining efermi and occupation")')
    else
    ! metal found
      e0 = evalint( 1, 1)
      e1 = e0
      do ist = wf_fst, wf_lst
        e0 = min( e0, minval( evalint( ist, :)))
        e1 = max( e1, maxval( evalint( ist, :)))
      end do

      it = 0
      do while( it .lt. maxit)
        efermiint = 0.5*(e0 + e1)
        chg = 0.d0
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iq, ist, x) reduction(+: chg)
!$OMP DO  
#endif
        do iq = 1, int_kset%nkpt
          do ist = wf_fst, wf_lst
            x = (efermiint - evalint( ist, iq))*t1
            occint( ist, iq) = occmax*stheta( input%groundstate%stypenumber, x)
            chg = chg + int_kset%wkpt( iq)*occint( ist, iq)
          end do
        end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
        if( chg .lt. chgval) then
          e0 = efermiint
        else
          e1 = efermiint
        end if
        if( (e1-e0) .lt. input%groundstate%epsocc) then
          it = maxit+1
        else
          it = it + 1
        end if
      end do
    end if
    if( it .eq. maxit) then
      write( *, '("Error (wannier_interpolate_occupancy): Fermi energy could not be found.")')
      call terminate
    end if
  else
    write( *, '("Error (wannier_interpolate_occupancy): Not implemented for this stype.")')
    call terminate
  end if

  return
end subroutine wannier_interpolate_occupancy

end module m_wannier_interpolate_density
