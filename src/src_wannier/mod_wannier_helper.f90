! Collection of helper functions and wrappers commonly used 
! at various occasions
module mod_wannier_helper
  use mod_wannier_variables

  use mod_large_io,             only: inquire_large, open_direct_unformatted_large
  use mod_spin,                  only: nspinor, nspnfv
  use svlo,                      only: get_num_of_basis_functions_sv
  use mod_eigensystem,           only: nmatmax_ptr, nmatmax, nmat_ptr, nmat, npmat
  use mod_APW_LO,                only: nlotot
  use mod_lattice,               only: bvec
  use mod_Gvector,               only: intgv
  use mod_Gkvector!,              only: ngk_ptr, vgkl_ptr, vgkl, vgkc
  use mod_eigenvalue_occupancy,  only: nstfv, nstsv, occmax
  use mod_charge_and_moment,     only: chgval
  use mod_potential_and_density, only: xctype
  use mod_misc,                  only: filext
  use mod_atoms,                 only: natmtot
  use sorting,                   only: sort_index_1d  
  use mod_gen_lo,                only: genlofr
  use mod_kpoint
  use precision,                 only: long_int

  implicit none

! methods
  contains

    !=====================================================================================
    ! set the k-point set compatible with the global k-points
    ! in the respective context
    subroutine wfhelp_setkpts
      integer :: ik, ispn, nempty, ngridk(3)
      logical :: reducek
      real(8) :: vkloff(3)
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
          nstfv = min( minval( nmat_ptr), int( chgval/2.d0) + input%groundstate%nempty + 1)
        case( "gw")
          input%groundstate%stypenumber = -1 ! turn on LIBBZINT
          if (xctype(1) >= 400) then
            ! GW@hybrids
            nstfv = min( minval( nmat_ptr), int( chgval/2.d0) + input%groundstate%nempty + 1)
            call init1
            call generate_k_vectors( wf_kset, bvec, input%groundstate%ngridk, input%groundstate%vkloff, .false.)
            call generate_k_vectors( wf_kset_red, bvec, input%groundstate%ngridk, input%groundstate%vkloff, .true.)
          else
            ! GW@DFT
            reducek = input%groundstate%reducek
            input%groundstate%reducek = .false.
            ngridk = input%groundstate%ngridk
            input%groundstate%ngridk = input%gw%ngridq
            vkloff = input%groundstate%vkloff
            input%groundstate%vkloff = input%gw%vqloff
            nempty = input%groundstate%nempty
            input%groundstate%nempty = input%gw%nempty
            call init1
            input%groundstate%reducek = reducek
            input%groundstate%ngridk = ngridk
            input%groundstate%vkloff = vkloff
            input%groundstate%nempty = nempty
            call generate_k_vectors( wf_kset, bvec, input%gw%ngridq, input%gw%vqloff, .false.)
            call generate_k_vectors( wf_kset_red, bvec, input%gw%ngridq, input%gw%vqloff, .true.)
          end if
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
          nstfv = min( minval( nmat_ptr), int( chgval/2.d0) + input%groundstate%nempty + 1)
        case( "hybrid")
          nstfv = min( minval( nmat_ptr), int( chgval/2.d0) + input%groundstate%nempty + 1)
          input%groundstate%stypenumber = -1 ! turn on LIBBZINT
          call generate_k_vectors( wf_kset, bvec, input%groundstate%ngridk, input%groundstate%vkloff, .false.)
          call generate_k_vectors( wf_kset_red, bvec, input%groundstate%ngridk, input%groundstate%vkloff, .true.)
          call init1
        case default
          if( mpiglobal%rank .eq. 0) then
            write(*,*)
            write(*, '("Error (wfhelp_setkpt): ",a," is not a valid input.")') input%properties%wannier%input
          end if
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
      call getngkmax
      do ik = 1, nkpt
        do ispn = 1, nspnfv
          call gengpvec( vkl( :, ik), vkc( :, ik), &
                  ngk( ispn, ik), &
                  igkig( :, ispn, ik), &
                  vgkl( :, :, ispn, ik), &
                  vgkc( :, :, ispn, ik), &
                  gkc( :, ispn, ik), &
                  tpgkc( :, :, ispn, ik))
          call gensfacgp( ngk( ispn, ik), vgkc( :, :, ispn, ik), ngkmax, &
                  sfacgk( :, :, ispn, ik))
          nmat( ispn, ik) = ngk( ispn, ik) + nlotot
          nmatmax = max( nmatmax, nmat( ispn, ik))
          npmat( ispn, ik) = (nmat( ispn, ik)*(nmat( ispn, ik) + 1))/2
        end do
      end do

      return
    end subroutine wfhelp_setkpts

    !=====================================================================================
    ! generate radial functions and related quantities
    subroutine wfhelp_genradfun
      character(256) :: fxt

      fxt = filext

      if( input%properties%wannier%input .eq. 'hybrid') then
        filext = '_PBE.OUT'
      else if ( input%properties%wannier%input .eq. 'gw') then
        if (xctype(1) >= 400) then
          ! case of GW@hybrids
          filext = '_PBE.OUT'
        end if
      end if

      call readstate
      call linengy
      call genapwfr     ! APW radial functions
      call genlofr      ! LO radial functions
      call olprad

      filext = trim( fxt)

      return
    end subroutine wfhelp_genradfun

    !=====================================================================================
    !> wrapper to fetch eigenenvectors in the respective context
    subroutine wfhelp_getevec( ik, evec, sort)
      use constants, only: zone, zzero
      use svlo, only: get_num_of_basis_functions_sv

      integer, intent( in)           :: ik
      logical, intent(in), optional  :: sort
      complex(8), intent( out)       :: evec( nmatmax_ptr, nstsv, nspinor)

      integer :: ispn, ist, num_of_basis_funs_sv, nmatk
      complex(8), allocatable :: auxmat(:,:)
      complex(8), allocatable :: evecfv(:,:,:) ! EVECFV.OUT eigenvectors
      complex(8), allocatable :: evecsv(:,:)   ! EVECSV.OUT eigenvectors
      logical :: sort_
      character(22) :: filext0

      sort_= .true.
      if( present( sort)) sort_ = sort
      allocate( evecfv( nmatmax_ptr, nstfv, nspnfv ) )
      allocate( evecsv( nstsv, nstsv ) )

      ! load EVECFV, EVECSV from file
      if( (input%properties%wannier%input .eq. "gs") .or. (input%properties%wannier%input .eq. "qsgw")) then
        call getevecfv( wf_kset%vkl( :, ik ), wf_Gkset%vgkl( :, :, :, ik ), evecfv )
        call wfhelp_getevecsv( ik, evecsv, sort=.false. )
      else if( input%properties%wannier%input .eq. "hybrid") then
        call getevecfv( wf_kset%vkl( :, ik ), wf_Gkset%vgkl( :, :, :, ik ), evecfv )
        call wfhelp_getevecsv( ik, evecsv, sort=.false. )
      else if( input%properties%wannier%input .eq. "gw") then
        if (xctype(1) >= 400) then
          call getevecfv( wf_kset%vkl( :, ik ), wf_Gkset%vgkl( :, :, :, ik ), evecfv )
          call wfhelp_getevecsv( ik, evecsv, sort=.false. )
        else
          filext0 = filext
          filext  = "_GW.OUT"
          call getevecfv( wf_kset%vkl( :, ik ), wf_Gkset%vgkl( :, :, :, ik ), evecfv )
          call wfhelp_getevecsv( ik, evecsv, sort=.false. )
          filext = filext0
        end if
      else
        stop
      end if

      ! set evec (FV), or convert EVECFV/EVECSV from file to evec LAPW+lo-expansion (SV)
      if ( .not. associated(input%groundstate%spin) ) then
        ! no spin
        evec = evecfv
      else
        ! spin
        if ( .not. issvlo() ) then
          ! SV
          ! sigma 1 component
          call zgemm( 'n', 'n', nmatmax_ptr, nstsv, nstfv, zone, &
                evecfv( 1, 1, 1 ), nmatmax_ptr, &
                evecsv( 1, 1 ), nstsv, zzero, &
                evec( 1, 1, 1 ), nmatmax_ptr)
          ! sigma 2 component
          call zgemm( 'n', 'n', nmatmax_ptr, nstsv, nstfv, zone, & 
                evecfv( 1, 1, 1 ), nmatmax_ptr, &
                evecsv( nstfv+1, 1 ), nstsv, zzero, &
                evec( 1, 1, 2), nmatmax_ptr)
        else 
          ! SVLO
          evec = zzero
          num_of_basis_funs_sv = get_num_of_basis_functions_sv()
          ! FVAPW times SVFV (APW) components
          ! sigma 1 component
          call zgemm( 'n', 'n', wf_Gkset%ngk( 1, ik ), nstsv, nstfv, zone, &
                evecfv( 1, 1, 1 ), nmatmax_ptr, &
                evecsv( 1, 1 ), nstsv, zzero, &
                evec( 1, 1, 1 ), nmatmax_ptr )
          ! sigma 2 component
          call zgemm( 'n', 'n', wf_Gkset%ngk( 1, ik ), nstsv, nstfv, zone, &
                evecfv( 1, 1, 1 ), nmatmax_ptr, &
                evecsv( num_of_basis_funs_sv + 1, 1 ), nstsv, zzero, &
                evec( 1, 1, 2 ), nmatmax_ptr )
          ! SVLO (LO) components
          ! sigma 1 component
          nmatk = wf_Gkset%ngk( 1, ik) + nlotot
          evec( wf_Gkset%ngk( 1, ik )+1 : nmatk, 1:nstsv, 1 ) = evecsv( nstfv+1 : nstfv+nlotot, 1:nstsv)
          ! sigma 2 component
          evec( wf_Gkset%ngk( 1, ik )+1 : nmatk, 1:nstsv, 2 ) = evecsv( num_of_basis_funs_sv+nstfv+1 : num_of_basis_funs_sv+nstfv+nlotot, 1:nstsv)
        end if
      end if
      deallocate( evecfv )
      deallocate( evecsv )
      
      ! sort evec according to eval-sorting
      if ( sort_ ) then
        if ( allocated( wf_index_map ) ) then
          allocate( auxmat(size(evec, dim=1), size(evec, dim=2)) )
          do ispn = 1, nspinor
            auxmat = evec(:, :, ispn)
            evec(:, :, ispn) = auxmat(:, wf_index_map(:, ik))
          end do
          deallocate( auxmat )
        else
          if( mpiglobal%rank .eq. 0) then
            write(*,*)
            write( *,'("Error (wfhelp_getevec): No index-map sorting eigenvalues/eigenvectors allocated!")')
          end if
          stop
        end if
      end if


      ! phase correction !NOTE feature unused currently. Adapt for spinor evec
      ! if( wf_fixphases) then
      !   !write(*,*) 'phase correction'
      !   allocate( auxmat( nmatmax_ptr, nstfv))
      !   auxmat = evec( :, :, 1)
      !   call zgemm( 'n', 'n', nmatmax_ptr, nstfv, nstfv, zone, &
      !          auxmat, nmatmax_ptr, &
      !          wf_evecphases( :, :, ik), nstfv, zzero, &
      !          evec( :, :, 1), nmatmax_ptr)
      !   deallocate( auxmat)
      ! end if

      return
    end subroutine wfhelp_getevec

    !=====================================================================================
    ! wrapper to fetch second variation eigenenvectors in the respective context
    subroutine wfhelp_getevecsv( ik, evecsv, sort)
      integer, intent(in)           :: ik
      logical, intent(in), optional  :: sort
      complex(8), intent(out)       :: evecsv(nstsv, nstsv) !NOTE for future svlo nstfv has to be replaced by get_num_of_basis_functions_sv

      complex(8), allocatable :: auxmat(:,:)
      logical :: sort_
      character(22) :: filext0

      sort_= .true.
      if( present( sort ) ) sort_ = sort

      ! SVLO failsafe. Getevecsv only used in mod_wannier_spin for S_z calculation in wannier gauge
      if( issvlo() ) then
        if( mpiglobal%rank .eq. 0) then
          write(*,*)
          write( *,'("Error (wfhelp_getevecsv): no svlo capability yet!")')
          stop
        end if
      end if

      ! load EVECSV from file
      if( (input%properties%wannier%input .eq. "gs") .or. (input%properties%wannier%input .eq. "qsgw")) then
        call getevecsv( wf_kset%vkl( :, ik ), evecsv )
      else if( input%properties%wannier%input .eq. "hybrid") then
        call getevecsv( wf_kset%vkl( :, ik ), evecsv )
      else if( input%properties%wannier%input .eq. "gw") then
        if (xctype(1) >= 400) then
          call getevecsv( wf_kset%vkl( :, ik ), evecsv )
        else
          filext0 = filext
          filext  = "_GW.OUT"
          call getevecsv( wf_kset%vkl( :, ik ), evecsv )
          filext = filext0
        end if
      else
        stop
      end if

      ! sort evec according to eval-sorting
      if ( sort_ ) then
        if ( allocated( wf_index_map ) ) then
          allocate( auxmat, source=evecsv )
          evecsv(:, :) = auxmat(:, wf_index_map(:, ik))
          deallocate( auxmat )
        else
          if( mpiglobal%rank .eq. 0) then
            write(*,*)
            write( *,'("Error (wfhelp_getevecsv): No index-map sorting  eigenvalues/eigenvectors allocated!")')
          end if
          stop
        end if
      end if
    end subroutine wfhelp_getevecsv

    !=====================================================================================
    ! wrapper to fetch eigenenergies in the respective context
    subroutine wfhelp_geteval( eval, fst, lst, mode, reduce, sort)
      real(8), allocatable, intent( out)          :: eval(:,:)
      integer, intent( out)                       :: fst, lst
      character(*), optional, intent( in)         :: mode
      logical, optional, intent( in)              :: reduce
      logical, optional, intent( in)              :: sort
      
     
      integer :: ik, ikk, ist, un, nkpqp, fstqp, lstqp, nk, iq, nkequi, isymequi( wf_kset%nkpt), ikequi( wf_kset%nkpt)
      integer(long_int) :: recl
      real(8) :: vl(3), efermiqp, efermiks
      character(256) :: mode_, fname, fxt
      logical :: reduce_, exist
      logical :: sort_
      type( k_set) :: kset

      real(8), allocatable :: evalqp(:), evalks(:), eval_file(:,:), evalsv(:)

      mode_ = input%properties%wannier%input
      if( present( mode)) mode_ = trim( mode)
      sort_= .true.
      if( present( sort)) sort_ = sort
      reduce_ = .false.
      if( present( reduce)) reduce_ = reduce
      if( allocated( eval)) deallocate( eval)
        ! KS energies on GS grid ('gs')
        ! generalized KS energies on GS grid ('hybrid')
        ! KS energies on BSE grid ('bse')
        if( (mode_ .eq. 'gs') .or. (mode_ .eq. 'hybrid') .or. (mode_ .eq. 'bse')) then
          fxt = filext
          nk = nstsv
          if( mode_ .eq. 'bse') then
            call generate_k_vectors( kset, wf_kset%bvec, input%xs%ngridk, input%xs%vkloff, reduce_)
            nstsv = int( chgval/2.d0) + input%xs%nempty + 1
          else
            filext = '.OUT'
            call generate_k_vectors( kset, wf_kset%bvec, input%groundstate%ngridk, input%groundstate%vkloff, reduce_)
          end if
          fst = 1
          lst = nstsv
          allocate( eval( fst:lst, kset%nkpt))
          allocate( eval_file( nstsv, nspnfv))
          do ik = 1, kset%nkpt
            call getevalsv( kset%vkl( :, ik), eval_file(:, 1))
            ikk = ik
            if( mode_ .ne. 'bse') call findkptinset( kset%vkl( :, ik), wf_kset, ist, ikk)
            eval(:, ikk) = eval_file(:, 1)
          end do
          deallocate( eval_file)
          filext = trim( fxt)
          nstsv = nk
        ! QP energies on GW grid ('gw')
        ! KS energies on GW grid ('gwks')
        else if( (mode_ .eq. 'gw') .or. (mode_ .eq. 'gwks')) then
          call generate_k_vectors( kset, wf_kset%bvec, input%gw%ngridq, input%gw%vqloff, reduce_)
          if ( nspinor == 1 ) then
            write( fname, '("EVALQP.OUT")')
          else
            ! spinor case
            write( fname, '("EVALQPSV.OUT")')
          end if
          inquire( file=trim( fname), exist=exist)
          if( .not. exist) then
            if( mpiglobal%rank .eq. 0) then
              write(*,*)
              write( *, '("Error (wfhelp_geteval): File EVALQP.OUT does not exist!")')
            end if
            stop
          end if
          call inquire_large( recl, [nkpqp, fstqp, lstqp] )
          call open_direct_unformatted_large( un, trim( fname), "read", recl, "old" )
          read( un, rec=1) nkpqp, fstqp, lstqp
          close( un)
          allocate( evalqp( fstqp:lstqp))
          allocate( evalks( fstqp:lstqp))
          fst = min( fstqp, nstsv )
          lst = min( lstqp, nstsv )
          allocate( eval( fst:lst, kset%nkpt))
          call inquire_large( recl, [nkpqp, fstqp, lstqp], vl, evalqp, evalks, [efermiqp, efermiks] )
          call open_direct_unformatted_large( un, trim( fname), "read", recl, "old" )
          do ik = 1, nkpqp
            read( un, rec=ik) nkpqp, fstqp, lstqp, vl, evalqp, evalks, efermiqp, efermiks
            call findequivkpt( vl, kset, nkequi, isymequi, ikequi)
            do iq = 1, nkequi
              if( mode_ .eq. 'gwks') then
                eval( :, ikequi( iq)) = evalks(fst:lst)
              else
                eval( :, ikequi( iq)) = evalqp(fst:lst)
              end if
            end do
          end do
          close( un)
          deallocate( evalqp, evalks)
        else
          if( mpiglobal%rank .eq. 0) then
            write(*,*)
            write( *, '("Error (wfhelp_geteval): Given mode not supported.")')
          end if
          stop
        end if
        ! sort evals in ascending order and sort evals if sort given
        if ( .not. allocated( wf_index_map ) ) then
          ! generate wf_index_map
          allocate( wf_index_map(nstsv, wf_kset%nkpt) )
          do ik = 1, wf_kset%nkpt
            wf_index_map(:, ik) = [(ist, ist=1, nstsv)]
            wf_index_map(fst:lst, ik) = fst + sort_index_1d( lst-fst+1, eval(:, ik) ) - 1
          end do
        end if
        if ( sort_ ) then
          do ik = 1, wf_kset%nkpt
            eval(:, ik) = eval(wf_index_map(fst:lst, ik), ik)
          end do
        end if
      
      return
    end subroutine wfhelp_geteval

    !=====================================================================================
    ! generic routine to calculate occupations and Fermi energy
    ! for a given set of eigenenergies
    subroutine wfhelp_occupy( kset, eval, fst, lst, efermi, occ, tetra)
      use mod_opt_tetra
      use mod_occupy, only: find_fermi

      type( k_set), intent( in)           :: kset
      integer, intent( in)                :: fst, lst
      real(8), intent( in)                :: eval(fst:lst, kset%nkpt)
      real(8), intent( out)               :: efermi
      real(8), intent( out)               :: occ(fst:lst, kset%nkpt)
      type( t_set), optional, intent( in) :: tetra

      integer, parameter :: maxit = 1000

      integer :: iq, ist, nvm, it, nst
      logical :: usetetra
      real(8) :: e0, e1, chg, chg0, x, t1, df
      
      real(8), external :: stheta

      nst = lst - fst + 1
      usetetra = .false.
      if( present( tetra)) usetetra = .true.
      df = 1.d-1
      chg0 = 0.d0; occ = 0.d0

      nvm = nint( chgval/occmax)
      if( (fst .ne. 1) .and. (fst .le. nvm)) then
        if( mpiglobal%rank .eq. 0) then
          write(*,*)
          write( *, '("Warning (wfhelp_occupy): The lowest band given is ",I3,". All bands below are considered to be fully occupied.")') fst
        end if
        chg0 = (fst-1)*occmax
      end if
      if( fst .gt. nvm) then
        if( mpiglobal%rank .eq. 0) then
          write(*,*)
          write( *, '("Warning (wfhelp_occupy): No valence bands given. All bands are considered to be unoccupied. Fermi energy set to lowest energy given.")')
        end if
        occ = 0.d0
        efermi = minval( eval(fst,:))
        return
      end if
      if( (lst .le. nvm)) then
        if( mpiglobal%rank .eq. 0) then
          write(*,*)
          write( *, '("Warning (wfhelp_occupy): At least one conduction band has to be given in order to determine occupancies. All bands given are considered to be fully occupied. Fermi energy set to highest energy given.")')
        end if
        occ = occmax
        efermi = maxval( eval(lst,:))
        return
      end if
      ! check for insulator or semiconductor
      e0 = maxval( eval( nvm, :))
      e1 = minval( eval( nvm+1, :))
      efermi = 0.5*(e0 + e1)

      chg = chg0
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iq, ist) reduction(+: chg)
!$OMP DO
#endif
      do iq = 1, kset%nkpt
        do ist = fst, lst
          if( eval( ist, iq) .le. efermi) occ( ist, iq) = occmax
          chg = chg + kset%wkpt( iq)*occ( ist, iq)
        end do
      end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif

      ! system has gap
      if( (e1 .ge. e0) .and. (abs( chg - chgval) .lt. input%groundstate%epsocc)) &
        return

      ! system is metallic
      ! use improved tetrahedron integration
      if( usetetra .or. input%groundstate%stypenumber == -2) then
        if( mpiglobal%rank .eq. 0) then
          write(*,*)
          write( *, '("Info (wfhelp_occupy): Use tetrahedron method in determining efermi and occupation")')
        end if
        call opt_tetra_efermi( tetra, chgval/dble( occmax)-fst+1, kset%nkpt, nst, eval, efermi, occ, ef0=efermi, df0=df)
        do iq = 1, kset%nkpt
          occ(:,iq) = occmax*occ(:,iq)/kset%wkpt(iq)
        end do
      ! use standard smearing method
      else if( input%groundstate%stypenumber >= 0 ) then
        call find_fermi( kset%nkpt, kset%wkpt, nst, eval, chgval-chg0, occmax, &
          input%groundstate%stypenumber, input%groundstate%swidth, input%groundstate%epsocc, &
          efermi, occ )
      ! use libbzint
      else if( input%groundstate%stypenumber == -1) then
        ! Calculate the Fermi energy
        call fermi_exciting( .false., chgval-chg0, nst, kset%nkpt, eval, &
                             kset%ntet, kset%tnodes, kset%wtet, kset%tvol, &
                             efermi, e0, e1)                               
        ! Calculate state occupation numbers
        call tetiw( kset%nkpt, kset%ntet, nst, eval, kset%tnodes, kset%wtet, kset%tvol, efermi, occ)
        do iq = 1, kset%nkpt
          do ist = fst, lst
            occ(ist,iq) = occmax/kset%wkpt(iq)*occ(ist,iq)
          end do
        end do
      end if

      return
    end subroutine wfhelp_occupy

    !=====================================================================================
    ! calculate the Fermi energy
    ! (dependent on context)
    subroutine wfhelp_getefermi( efermi, tetra)
      use mod_opt_tetra
      real(8), intent(out) :: efermi
      type(t_set), optional, intent(in) :: tetra

      integer :: fst, lst
      real(8), allocatable :: evalfv(:,:), occ(:,:)

      call wfhelp_geteval( evalfv, fst, lst)
      allocate( occ( fst:lst, wf_kset%nkpt))
      if( present( tetra)) then
        call wfhelp_occupy( wf_kset, evalfv, fst, lst, efermi, occ, tetra=tetra)
      else
        call wfhelp_occupy( wf_kset, evalfv, fst, lst, efermi, occ)
      end if
      deallocate( occ)

      return
    end subroutine wfhelp_getefermi

    !=====================================================================================
    ! initialize set of tetrahedra for tetrahedron integration
    subroutine wfhelp_init_tetra( tetra, kset, ttype, reduce)
      use mod_opt_tetra

      type( t_set), intent( inout) :: tetra
      type( k_set), intent( in)        :: kset
      integer, optional, intent( in)   :: ttype
      logical, optional, intent( in)   :: reduce

      integer :: ttype_ = 1
      logical :: reduce_ = .false.

      if( present( ttype)) ttype_ = ttype
      if( present( reduce)) reduce_ = reduce

      if( .true. .or. .not. tetra%initialized) call opt_tetra_init( tetra, kset, ttype_, reduce_)
      return
    end subroutine wfhelp_init_tetra

    !=====================================================================================
    ! this routine aims to fix the phases of the eigenvectors
    ! such that the transformation matrices U also work for
    ! eigenvectors from another (but sufficiently close) DFT run
    subroutine wfhelp_fixphases
      use constants, only: zone, zzero

      integer :: ik, ist, jst, nmatp
      real(8), allocatable :: eval(:,:)
      complex(8), allocatable :: evec(:,:,:), phase(:,:)

      if( wf_fixphases) then
        wf_fixphases = .false.
        if( allocated( wf_evecphases)) deallocate( wf_evecphases)
        call wfhelp_geteval( eval, ist, jst)
        allocate( wf_evecphases( nstfv, nstfv, wf_kset%nkpt))
        allocate( phase( ist:jst, ist:jst))
        wf_evecphases = zzero
        do ik = 1, nstfv
          wf_evecphases( ik, ik, :) = zone
        end do
        allocate( evec( nmatmax_ptr, nstfv, nspinor))
        do ik = 1, wf_kset%nkpt
          nmatp = wf_Gkset%ngk( 1, ik) + nlotot
          call wfhelp_getevec( ik, evec)
          call wfhelp_getphases( evec( 1:nmatp, ist:jst, 1), eval( ist:jst, ik), nmatp, jst-ist+1, phase)
          wf_evecphases( ist:jst, ist:jst, ik) = phase
        end do
        deallocate( evec, eval, phase)
        wf_fixphases = .true.
      end if
      return
    end subroutine wfhelp_fixphases

    subroutine wfhelp_getphases( evec, eval, ngp, nst, phases)
      use constants, only: zone, zzero
      use m_linalg,      only: zhediag

      complex(8), intent( in)    :: evec( ngp, nst)
      real(8), intent( in)       :: eval( nst)
      integer, intent( in)       :: ngp, nst
      complex(8), intent( out)   :: phases( nst, nst)

      integer :: ist, jst, kst, ndeg, igp
      real(8) :: epse, epsp
      complex(8) :: phase
      complex(8), allocatable :: pert(:,:), auxmat(:,:), p(:,:)
      real(8), allocatable :: pval(:)
      
      epse = 1.d-4
      epsp = 1.d-1

      phases = zzero
      do ist = 1, nst
        phases( ist, ist) = zone
      end do

      allocate( pert( ngp, ngp), auxmat( ngp, ngp), p( nst, nst), pval( nst))
      ! fix degeneracies
      ndeg = 1
      do jst = 1, nst-1
        if( (eval( jst+1) - eval( jst) .ge. epse) .and. (ndeg .gt. 1)) then
          ist = jst - ndeg + 1
          ! construct ficticious hermitian perturbation
          auxmat = zzero
          do kst = ist, jst
            auxmat( :, kst) = auxmat( :, kst) + evec( :, kst)
          end do
          call zgemm( 'c', 'n', ngp, ngp, ngp, zone, &
                 auxmat, ngp, &
                 auxmat, ngp, zzero, &
                 pert, ngp)
          call zgemm( 'c', 'n', ndeg, ngp, ngp, zone, &
                 evec( :, ist:jst), ngp, &
                 pert, ngp, zzero, &
                 auxmat( 1:ndeg, :), ndeg)
          call zgemm( 'n', 'n', ndeg, ndeg, ngp, zone, &
                 auxmat( 1:ndeg, :), ndeg, &
                 evec( :, ist:jst), ngp, zzero, &
                 p, nst)
          ! find unitary mixing of degenerate states
          pval = 0.d0
          call zhediag( p( 1:ndeg, 1:ndeg), pval( 1:ndeg), evec=phases( ist:jst, ist:jst))
        end if
        if( eval( jst+1) - eval( jst) .lt. epse) then
          ndeg = ndeg + 1
        else
          ndeg = 1
        end if
      end do
      ! fix global phase
      call zgemm( 'n', 'n', ngp, nst, nst, zone, &
             evec, ngp, &
             phases, nst, zzero, &
             auxmat( :, 1:nst), ngp)
      do ist = 1, nst
        phase = zzero
        do igp = 1, ngp
          phase = phase + auxmat( igp, ist)
          if( abs( phase) .gt. epsp) exit
        end do
        if( abs( phase) .gt. epsp) then
          phases( :, ist) = phases( :, ist)*conjg( phase)/abs( phase)
        else
          if( mpiglobal%rank .eq. 0) then
            write(*,*) 'global phase not found for state ', ist
          end if
        end if
      end do

      deallocate( pert, auxmat, p, pval)
      return
    end subroutine wfhelp_getphases

end module mod_wannier_helper
