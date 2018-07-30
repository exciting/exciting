module mod_wannier_helper
  use mod_wannier_variables

  use mod_spin,                 only: nspinor, nspnfv
  use mod_eigensystem,          only: nmatmax_ptr, nmat_ptr, nmat, nmatmax, npmat 
  use mod_APW_LO,               only: nlotot
  use mod_lattice,              only: bvec
  use mod_Gvector,              only: intgv
  use mod_Gkvector!,             only: ngk_ptr, vgkl_ptr, vgkl, vgkc
  use mod_eigenvalue_occupancy, only: nstfv, occmax
  use mod_charge_and_moment,    only: chgval
  use mod_misc,                 only: filext
  use mod_atoms,                only: natmtot
  use modinput
  use mod_kpoint
  use m_getunit

  implicit none

! methods
  contains

    subroutine wannier_setkpts
      integer :: ik, ispn
      !write(*,*) "k set"
      select case (input%properties%wannier%input)
        case( "gs")
          call generate_k_vectors( wf_kset, bvec, input%groundstate%ngridk, input%groundstate%vkloff, .false.)
          if( input%groundstate%reducek) then
            wf_kset%vkl = vklnr
            wf_kset%vkc = vkcnr
          else
            wf_kset%vkl = vkl
            wf_kset%vkc = vkc
          end if
          call generate_k_vectors( wf_kset_red, bvec, input%groundstate%ngridk, input%groundstate%vkloff, .true., .false.)
          nstfv = int( chgval/2.d0) + input%groundstate%nempty + 1
        case( "gw")
          call generate_k_vectors( wf_kset, bvec, input%gw%ngridq, input%gw%vqloff, input%gw%reduceq)
          vkl = wf_kset%vkl
          vkc = wf_kset%vkc
          nkpt = wf_kset%nkpt
          call generate_k_vectors( wf_kset, bvec, input%gw%ngridq, input%gw%vqloff, .false.)
          vklnr = wf_kset%vkl
          vkcnr = wf_kset%vkc
          nkptnr = wf_kset%nkpt
          call generate_k_vectors( wf_kset_red, bvec, input%gw%ngridq, input%gw%vqloff, .true.)
          nstfv = int( chgval/2.d0) + input%gw%nempty + 1
        case( "qsgw")
          call generate_k_vectors( wf_kset, bvec, input%gw%ngridq, input%gw%vqloff, input%gw%reduceq)
          if( allocated( vkl)) deallocate( vkl)
          allocate( vkl( 3, wf_kset%nkpt))
          vkl_ptr => vkl
          if( allocated( vkc)) deallocate( vkc)
          allocate( vkc( 3, wf_kset%nkpt))
          vkl = wf_kset%vkl
          vkc = wf_kset%vkc
          nkpt = wf_kset%nkpt
          call generate_k_vectors( wf_kset, bvec, input%gw%ngridq, input%gw%vqloff, .false.)
          vklnr = wf_kset%vkl
          vkcnr = wf_kset%vkc
          call generate_k_vectors( wf_kset_red, bvec, input%gw%ngridq, input%gw%vqloff, .true.)
          nstfv = int( chgval/2.d0) + input%gw%nempty + 1
        case( "hybrid")
          call readkpts
          !call generate_k_vectors( wf_kset, bvec, input%groundstate%ngridk, input%groundstate%vkloff, input%groundstate%reducek)
          !vkl = wf_kset%vkl
          !vkc = wf_kset%vkc
          !nkpt = wf_kset%nkpt
          !! recalculate G+k-vectors
          !do ik = 1, nkpt
          !  do is = 1, nspnfv
          !    vl (:) = vkl (:, ik)
          !    vc (:) = vkc (:, ik)
          !    call gengpvec( vl, vc, ngk( is, ik), igkig( :, is, ik), vgkl( :, :, is, ik), vgkc( :, :, is, ik), gkc( :, is, ik), tpgkc( :, :, is, ik))
          !    call gensfacgp( ngk(is, ik), vgkc( :, :, is, ik), ngkmax, sfacgk( :, :, is, ik))
          !  end do
          !end do
          call generate_k_vectors( wf_kset, bvec, input%groundstate%ngridk, input%groundstate%vkloff, .false.)
          call generate_k_vectors( wf_kset_red, bvec, input%groundstate%ngridk, input%groundstate%vkloff, .true., .true.)
          nstfv = int( chgval/2.d0) + input%groundstate%nempty + 1
        case default
          write(*,*)
          write(*, '("Error (wannier_setkpt): ",a," is not a valid input.")') input%properties%wannier%input
          stop
      end select
      call generate_G_vectors( wf_Gset, bvec, intgv, input%groundstate%gmaxvr)
      call generate_Gk_vectors( wf_Gkset, wf_kset, wf_Gset, gkmax)

      ! we may have changed the global vkl, so we have to update related variables
      if( allocated(ngk)) deallocate( ngk)
      allocate( ngk(nspnfv, nkpt))
      ngk_ptr => ngk
      if( allocated(igkig)) deallocate( igkig)
      allocate( igkig(ngkmax, nspnfv, nkpt))
      if( allocated(vgkl)) deallocate( vgkl)
      allocate( vgkl(3, ngkmax, nspnfv, nkpt))
      vgkl_ptr => vgkl
      if( allocated(vgkc)) deallocate( vgkc)
      allocate( vgkc(3, ngkmax, nspnfv, nkpt))
      if( allocated(gkc)) deallocate( gkc)
      allocate( gkc(ngkmax, nspnfv, nkpt))
      if( allocated(tpgkc)) deallocate( tpgkc)
      allocate( tpgkc(2, ngkmax, nspnfv, nkpt))
      if( allocated(sfacgk)) deallocate( sfacgk)
      allocate( sfacgk(ngkmax, natmtot, nspnfv, nkpt))
      if( allocated(nmat)) deallocate( nmat)
      allocate( nmat(nspnfv, nkpt))
      nmat_ptr => nmat
      if( allocated(npmat)) deallocate( npmat)
      allocate( npmat(nspnfv, nkpt))
      nmatmax = 0
      nmatmax_ptr => nmatmax
      do ik = 1, nkpt
        do ispn = 1, nspnfv
          call gengpvec( vkl( :, ik), vkc( :, ik), &
                  ngk( ispn, ik), &
                  igkig( :, ispn, ik), &
                  vgkl( :, :, ispn, ik), &
                  vgkc( :, :, ispn, ik), &
                  gkc( :, ispn, ik), &
                  tpgkc( :, :, ispn, ik))
          call getngkmax
          call gensfacgp( ngk( ispn, ik), vgkc( :, :, ispn, ik), ngkmax, &
                  sfacgk( :, :, ispn, ik))
          nmat( ispn, ik) = ngk( ispn, ik) + nlotot
          nmatmax = max( nmatmax, nmat( ispn, ik))
          npmat( ispn, ik) = (nmat( ispn, ik)*(nmat( ispn, ik) + 1))/2
        end do
      end do

      return
    end subroutine wannier_setkpts

    subroutine wannier_genradfun
      character(256) :: fname

      fname = filext

      if( input%properties%wannier%input .eq. 'hybrid') then
        filext = '_PBE.OUT'
      else
        filext = '.OUT'
      end if
      
      call readstate
      call linengy
      call genapwfr     ! APW radial functions
      call genlofr      ! LO radial functions
      call olprad

      filext = fname

      return
    end subroutine wannier_genradfun

    subroutine wannier_getevec( ik, evec)
      integer, intent( in) :: ik
      complex(8), intent( out) :: evec( nmatmax_ptr, nstfv, nspinor)

      integer :: ist
      real(8) :: phase

      if( (input%properties%wannier%input .eq. "gs") .or. (input%properties%wannier%input .eq. "qsgw")) then
        call getevecfv( wf_kset%vkl(:, ik), wf_Gkset%vgkl( :, :, :, ik), evec)
      else if( input%properties%wannier%input .eq. "hybrid") then
        call getevecfv( wf_kset%vkl(:, ik), wf_GKset%vgkl( :, :, :, ik), evec)
      else if( input%properties%wannier%input .eq. "gw") then
        call getevecsvgw_new( "GW_EVECSV.OUT", ik, wf_kset%vkl( :, ik), nmatmax_ptr, nstfv, nspinor, evec)
      else
        stop
      end if

      !call plotmat( evec( :, :, 1))
      !write(*,*)

      ! phase correction
      !if( allocated( wf_evecphase)) then
      !  do ist = wf_fst, wf_lst
      !    phase = atan2( dble( aimag( evec( wf_evecphase( ist, ik), ist, 1))), dble( evec( wf_evecphase( ist, ik), ist, 1)))
      !    evec( :, ist, 1) = evec( :, ist, 1)*cmplx( cos( phase), sin( -phase), 8)
      !  end do
      !end if

      return
    end subroutine wannier_getevec

    subroutine wannier_geteval( eval, fst, lst)
      real(8), allocatable, intent( out) :: eval(:,:)
      integer, intent( out) :: fst, lst

      integer :: ik, ist, un, recl, nkpqp, fstqp, lstqp, nk, isym, iq, nkequi, isymequi( wf_kset%nkpt), ikequi( wf_kset%nkpt)
      real(8) :: vl(3), efermiqp, efermiks
      character(256) :: fname
      logical :: exist

      real(8), allocatable :: evalqp(:), evalks(:), evalfv(:,:)

      if( input%properties%wannier%input .eq. "gw") then
        call getunit( un)
        write( fname, '("EVALQP.OUT")')
        inquire( file=trim( fname), exist=exist)
        if( .not. exist) then
          write(*,*)
          write( *, '("Error (wfint_init): File EVALQP.OUT does not exist!")')
          stop
        end if
        inquire( iolength=recl) nkpqp, fstqp, lstqp
        open( un, file=trim( fname), action='read', form='unformatted', access='direct', recl=recl)
        read( un, rec=1), nkpqp, fstqp, lstqp
        close( un)
        allocate( evalqp( fstqp:lstqp))
        allocate( evalks( fstqp:lstqp))
        allocate( eval( fstqp:lstqp, wf_kset%nkpt))
        fst = fstqp
        lst = lstqp
        inquire( iolength=recl) nkpqp, fstqp, lstqp, vl, evalqp, evalks, efermiqp, efermiks
        open( un, file=trim( fname), action='read', form='unformatted', access='direct', recl=recl)
        do ik = 1, nkpqp
          read( un, rec=ik) nkpqp, fstqp, lstqp, vl, evalqp, evalks, efermiqp, efermiks
          call findequivkpt( vl, wf_kset, nkequi, isymequi, ikequi)
          do iq = 1, nkequi
            eval( :, ikequi( iq)) = evalqp(:)
          end do
        end do
        close( un)
        deallocate( evalqp, evalks)
      else
        allocate( evalfv( nstfv, nspnfv))
        allocate( eval( nstfv, wf_kset%nkpt))
        fst = 1
        lst = nstfv
        do ik = 1, wf_kset%nkpt
          call getevalfv( wf_kset%vkl( :, ik), evalfv)
          eval( :, ik) = evalfv( :, 1)
        end do
        deallocate( evalfv)
      end if

      return
    end subroutine wannier_geteval

    subroutine wannier_occupy( kset, eval, fst, lst, efermi, occ, usetetra_)
      use mod_opt_tetra

      type( k_set), intent( in) :: kset 
      integer, intent( in) :: fst, lst
      real(8), intent( in) :: eval( fst:lst, kset%nkpt)
      logical, optional, intent( in) :: usetetra_
      real(8), intent( out) :: efermi
      real(8), optional, intent( out) :: occ( fst:lst, kset%nkpt)

      integer, parameter :: maxit = 1000
      
      integer :: iq, ist, nvm, it
      logical :: usetetra
      real(8) :: e0, e1, chg, x, t1, df, occ_tmp( lst, kset%nkpt)
      
      real(8) :: sdelta, stheta
      
      usetetra = .true.
      occ_tmp = 0.d0
      df = 1.d-2

      nvm = nint( chgval/occmax)
      if( (fst .ne. 1) .and. (fst .le. nvm)) then
        write(*,*)
        write( *, '("Warning (wannier_occupy): The lowest band given is ",I3,". All bands below are considered to be fully occupied.")') fst
        occ_tmp( 1:(fst-1), :) = occmax
      end if
      if( fst .gt. nvm) then
        write(*,*)
        write( *, '("Warning (wannier_occupy): No valence bands given. All are considered to be unoccupied. Fermi energy set to lowest energy given.")')
        if( present( occ)) occ = 0.d0
        efermi = minval( eval( fst, :))
        return
      end if
      if( (lst .le. nvm)) then
        write(*,*)
        write( *, '("Warning (wannier_occupy): At least one conduction band has to be given in order to determine occupancies. All bands given are considered to be fully occupied. Fermi energy set to highest energy given.")')
        if( present( occ)) occ = occmax
        efermi = maxval( eval( lst, :))
        return
      end if
      ! check for insulator or semiconductor
      e0 = maxval( eval( nvm, :))
      e1 = minval( eval( nvm+1, :))
      efermi = 0.5*(e0 + e1)
    
      chg = 0.d0
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iq, ist) reduction(+: chg)
!$OMP DO  
#endif
      do iq = 1, kset%nkpt
        do ist = 1, fst-1
          chg = chg + kset%wkpt( iq)*occ_tmp( ist, iq)
        end do
        do ist = fst, lst
          if( eval( ist, iq) .le. efermi) occ_tmp( ist, iq) = occmax
          chg = chg + kset%wkpt( iq)*occ_tmp( ist, iq)
        end do
      end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
      if( (e1 .ge. e0) .and. (abs( chg - chgval) .lt. input%groundstate%epsocc)) then
        write(*,*)
        write( *, '("Info (wannier_occupy): System has gap. Fermi level set to the middle of the gap.")')
        usetetra = .false.
      else
        ! metal found
        if( input%groundstate%stypenumber .ge. 0 ) then
          t1 = 1.d0/input%groundstate%swidth
          it = 0
          e0 = eval( fst, 1)
          e1 = e0
          do ist = fst, lst
            e0 = min( e0, minval( eval( ist, :)))
            e1 = max( e1, maxval( eval( ist, :)))
          end do
    
          do while( it .lt. maxit)
            efermi = 0.5*(e0 + e1)
            chg = 0.d0
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iq, ist, x) reduction(+: chg)
!$OMP DO  
#endif
            do iq = 1, kset%nkpt
              do ist = 1, fst-1
                chg = chg + kset%wkpt( iq)*occ_tmp( ist, iq)
              end do
              do ist = fst, lst
                x = (efermi - eval( ist, iq))*t1
                occ_tmp( ist, iq) = occmax*stheta( input%groundstate%stypenumber, x)
                chg = chg + kset%wkpt( iq)*occ_tmp( ist, iq)
              end do
            end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            if( chg .lt. chgval) then
              e0 = efermi
            else
              e1 = efermi
            end if
            if( (e1-e0) .lt. input%groundstate%epsocc) then
              it = maxit+1
            else
              it = it + 1
            end if
          end do
          usetetra = .true.
          if( present( usetetra_)) usetetra = usetetra_
          if( it .eq. maxit) then
            write(*,*)
            write( *, '("Error (wannier_occupy): Fermi energy could not be found.")')
            stop
          end if
        else
          usetetra = .true.
          df = 1.d0
          if( present( usetetra_)) usetetra = usetetra_
          if( .not. usetetra) then
            write(*,*)
            write( *, '("Error (wannier_occupy): Not implemented for this stype.")')
            stop
          end if
        end if
      end if

      if( usetetra) then
        write(*,*)
        write( *, '("Info (wannier_occupy): Use tetrahedron method in determining efermi and occupation")')
        call opt_tetra_init( 2, kset%bvec, kset%ngridk, kset%nkpt, kset%ikmap)
        call opt_tetra_efermi( chgval/dble( occmax)-fst+1, kset%nkpt, lst-fst+1, eval( fst:lst, :), efermi, occ_tmp( fst:lst, :), ef0=efermi, df0=df)
      end if

      if( present( occ)) occ(:,:) = occ_tmp( fst:lst, :)
      return
    end subroutine wannier_occupy

    subroutine wannier_getefermi( efermi)
      real(8), intent( out) :: efermi

      integer :: fst, lst
      real(8), allocatable :: evalfv(:,:)

      call wannier_geteval( evalfv, fst, lst)
      call wannier_occupy( wf_kset, evalfv, fst, lst, efermi)

      return
    end subroutine wannier_getefermi

end module mod_wannier_helper
