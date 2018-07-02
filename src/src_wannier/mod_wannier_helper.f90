module mod_wannier_helper
  use mod_wannier_variables

  use mod_spin,                 only: nspinor, nspnfv
  use mod_eigensystem,          only: nmatmax_ptr, nmat, npmat 
  use mod_APW_LO,               only: nlotot
  use mod_lattice,              only: bvec
  use mod_Gvector,              only: intgv
  use mod_Gkvector
  use mod_eigenvalue_occupancy, only: nstfv
  use mod_charge_and_moment,    only: chgval
  use mod_misc,                 only: filext
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
          write(*, '(" ERROR (wannier_setkpt): ",a," is not a valid input.")') input%properties%wannier%input
          stop
      end select
      call generate_G_vectors( wf_Gset, bvec, intgv, input%groundstate%gmaxvr)
      call generate_Gk_vectors( wf_Gkset, wf_kset, wf_Gset, gkmax)

      ! we may have changed the global vkl, so we have to update related variables
      do ik = 1, nkpt
        do ispn = 1, nspnfv
          call gengpvec( vkl( :, ik), vkc( :, ik), &
                  ngk( ispn, ik), &
                  igkig( :, ispn, ik), &
                  vgkl( :, :, ispn, ik), &
                  vgkc( :, :, ispn, ik), &
                  gkc( :, ispn, ik), &
                  tpgkc( :, :, ispn, ik))
          !call getngkmax
          call gensfacgp( ngk( ispn, ik), vgkc( :, :, ispn, ik), ngkmax, &
                  sfacgk( :, :, ispn, ik))
          nmat( ispn, ik) = ngk( ispn, ik) + nlotot
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

end module mod_wannier_helper
