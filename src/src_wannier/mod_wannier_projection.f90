module mod_wannier_projection
  use mod_wannier_variables
  use mod_wannier_helper

  use mod_atoms
  use mod_eigensystem
  use mod_APW_LO
  use mod_muffin_tin
  use mod_eigenvalue_occupancy
  use constants,                 only: zzero, zone, y00, twopi
  use mod_Gkvector,              only: ngkmax_ptr
  use mod_spin,                  only: nspinor
  use mod_potential_and_density, only: veffmt

  implicit none

  integer :: wfpro_noccmax = 10    ! maximum principle quantum number of occupied states
  integer :: wfpro_nunocc = 10     ! number of additional unoccupied states
  integer :: wfpro_dordmax = 1
  integer :: wfpro_nproj = 0
  real(8) :: wfpro_epsld = 1.d-3

  integer :: wfpro_nradfun
  integer :: wfpro_nprojtot, wfpro_nprojtot_n
  integer :: wfpro_lmax, wfpro_nmax
  integer :: wfpro_un
  
  integer, allocatable :: wfpro_nst(:), wfpro_states(:,:,:), wfpro_projst(:,:), wfpro_projused(:), wfpro_nn(:,:), wfpro_vn(:,:,:,:)
  real(8), allocatable :: wfpro_radfun(:,:)
  complex(8), allocatable :: wfpro_proj(:,:,:)

! methods
  contains

    ! finds atoms in neighboring unit cells such that
    ! all bonds are included in the set of atoms
    subroutine wfpro_neighcells
      integer :: i, j, k, is, ia, ias, jas, igroup, vi(3)
      real(8) :: a(3,3), b(3,3), d, vac1(3), val1(3), vac2(3), val2(3), vd(3), vr1(3), vr2(3), vr(3), aposc( 3, natmtot)
      integer(4) :: natom, nbond, nlbond, atoms( 27*natmtot), acell( 3, 27*natmtot), bonds( 2, 27*natmtot), lbonds( 27*natmtot)
      real(8) :: btol, d0, dist( natmtot, natmtot), bondv( 3, 27*natmtot), lbondv( 3, 27*natmtot)
      logical :: added

      btol = 0.2d0

      if( allocated( wfpro_nn)) deallocate( wfpro_nn)
      allocate( wfpro_nn( natmtot, wf_ngroups))
      if( allocated( wfpro_vn)) deallocate( wfpro_vn)
      allocate( wfpro_vn( 3, 27*natmtot, natmtot, wf_ngroups))

      wfpro_vn = 0
      wfpro_nn = 0

      added = .false.
      do igroup = 1, wf_ngroups
        if( wf_groups(igroup)%neighcells) then
          added = .true.
        else
          wfpro_nn(:,igroup) = 1
        end if
      end do
      if( .not. added) return

      a = input%structure%crystal%basevect
      call r3minv( a, b)
      do is = 1, nspecies
        do ia = 1, natoms( is)
          ias = idxas( ia, is)
          aposc( :, ias) = atposc( :, ia, is)
        end do
      end do

      ! shortest lattice vector
      vd = a(:,1)
      d0 = norm2( vd)
      do i = -3, 3
        vr1 = dble( i)*a(:,1)
        do j = -3, 3
          vr2 = vr1 + dble( j)*a(:,2)
          do k = -3, 3
            vr = vr2 + dble( k)*a(:,3)
            d = norm2( vd + vr)
            if( d .lt. input%structure%epslat) cycle
            d0 = min( d0, d)
          end do
        end do
      end do

      ! distance matrix
      do ias = 1, natmtot
        vac1 = aposc( :, ias)
        call r3mv( b, vac1, val1)
        do jas = 1, natmtot
          vac2 = aposc( :, jas)
          call r3mv( b, vac2, val2)
          vr = val2 - val1
          call r3ws( input%structure%epslat, a, vr, vi)
          call r3mv( a, vr, vd)
          d = norm2( vd)
          if( d .lt. input%structure%epslat) d = d0
          dist( ias, jas) = d
        end do
      end do

      ! add atoms
      natom = 0
      atoms = 0
      acell = 0
      nbond = 0
      bonds = 0
      bondv = 0.d0

      ! loop over all atoms in unit cell
      do ias = 1, natmtot
        vac1 = aposc( :, ias)
        d0 = (1.d0+btol)*minval( dist(:,ias))
        added = .false.
        do i = 1, natom
          if( atoms(i) .eq. ias .and. sum( abs( acell(:,i))) .eq. 0) then
            added = .true.
            exit
          end if
        end do
        nlbond = 0
        ! find local bonds
        do jas = 1, natmtot
          if( dist( jas, ias) .gt. d0) cycle
          vac2 = aposc( :, jas)
          vd = vac2 - vac1
          do i = -3, 3
            vr1 = dble( i)*a(:,1)
            do j = -3, 3
              vr2 = vr1 + dble( j)*a(:,2)
              do k = -3, 3
                vr = vr2 + dble( k)*a(:,3)
                if( abs( norm2( vd+vr) - dist( jas, ias)) .le. input%structure%epslat) then
                  nlbond = nlbond + 1
                  lbonds( nlbond) = jas
                  lbondv( :, nlbond) = vd + vr
                end if
              end do
            end do
          end do
        end do
        ! add bonds that are not yet included
        do i = 1, nlbond
          call r3mv( b, lbondv(:,i) + vac1 - aposc( :, lbonds(i)), vr)
          do j = 1, nbond
            if( (bonds(1,j) .eq. ias .and. bonds(2,j) .eq. lbonds(i) .and. norm2( bondv(:,j)-lbondv(:,i)) .lt. input%structure%epslat) .or. &
                (bonds(1,j) .eq. lbonds(i) .and. bonds(2,j) .eq. ias .and. norm2( bondv(:,j)+lbondv(:,i)) .lt. input%structure%epslat)) exit
          end do
          if( j .gt. nbond) then
            if( .not. added) then
              natom = natom + 1
              atoms( natom) = ias
              added = .true.
            end if
            nbond = nbond + 1
            bonds( :, nbond) = (/ias, lbonds(i)/)
            bondv( :, nbond) = lbondv(:,i)
            do j = 1, natom
              if( atoms(j) .eq. lbonds(i) .and. sum( abs( acell(:,j)-nint( vr))) .eq. 0) exit
            end do
            if( j .gt. natom) then
              natom = natom + 1
              atoms( natom) = lbonds(i)
              acell( :, natom) = nint( vr)
            end if
          end if
        end do
      end do

      do igroup = 1, wf_ngroups
        do ias = 1, natmtot
          do i = 1, natom
            if( (atoms(i) .eq. ias) .and. wf_groups( igroup)%neighcells) then
              wfpro_nn( ias, igroup) = wfpro_nn( ias, igroup) + 1
              wfpro_vn( :, wfpro_nn( ias, igroup), ias, igroup) = acell(:,i)
            end if
          end do
          wfpro_nn( ias, igroup) = max( 1, wfpro_nn( ias, igroup))
        end do
      end do

      return
    end subroutine wfpro_neighcells

    ! finds linearization energies for local orbital radial functions
    subroutine wfpro_getline( nmax, lmax, line)
      integer, intent( in) :: nmax, lmax
      real(8), intent( out) :: line( 0:nmax, 0:lmax, nspecies)

      integer :: is, nr, ias, l, n, nn
      real(8) :: vr( nrmtmax), p0s( nrmtmax), q0s( nrmtmax), q1s( nrmtmax), hp0( nrmtmax)
      real(8) :: ens( 0:nmax), elo, ehi, flo, fhi, emi, fmi

      ! copied from genlofr
      line = 0.d0
      do is = 1, nspecies
        nr = nrmt( is)
        ias = idxas( 1, is)
        vr(1:nr) = veffmt(1, 1:nr, ias) * y00
        do l = 0, lmax
          do n = 0, nmax
            ens( n) = 0.d0
            call rdirac( 0, n+l+1, l, l+1, nr, spr( :, is), vr, ens( n), p0s, q0s, .false., .false.)
          end do
          do n = 0, nmax
            ehi = ens( n)
            call rschroddme( 0, l, 0, ehi, nr, spr( :, is), vr, nn, p0s, hp0, q0s, q1s)
            fhi = hp0( nr)
            if( p0s( nr) .eq. p0s( nr-1)) then
              elo = ehi
            else
              if( n .eq. 0) then 
                elo = 2*ens(0) - ens(1) ! assuming lowest eigenenergy is negative
              else
                elo = ens( n-1)
              endif
              call rschroddme( 0, l, 0, elo, nr, spr( :, is), vr, nn, p0s, hp0, q0s, q1s) 
              flo = hp0( nr)
              if( ehi .lt. elo) then
                write(*,*)
                write(*, '("Error (wfpro_getline): Oops! This was not supposed to happen.")')
                stop
              endif
              do while( ehi - elo .gt. 1d-6)
                emi = 0.5d0*(ehi + elo)
                call rschroddme( 0, l, 0, emi, nr, spr( :, is), vr, nn, p0s, hp0, q0s, q1s)
                fmi = hp0( nr)
                if( fmi*fhi .lt. 0) then
                  flo = fmi
                  elo = emi
                else
                  fhi = fmi
                  ehi = emi
                endif
              end do
            end if
            line( n, l, is) = 0.5d0*(ens( n) + 0.5d0*(ehi + elo))
          end do
        end do
      end do

      return
    end subroutine wfpro_getline

    ! costructs a set of atomic states described by (n,l,m,o)
    ! according to the aufbau principle
    subroutine wfpro_getstates
      integer :: is, ia, ist, l, lmax, m, n, nl, nemax, ntotmax, nunocc
      integer, allocatable :: nnl(:,:,:), nmin(:)

      if( allocated( wfpro_nst)) deallocate( wfpro_nst)
      allocate( wfpro_nst( nspecies))
      if( allocated( wfpro_states)) deallocate( wfpro_states)
      allocate( wfpro_states( 2, 2*wfpro_noccmax+wfpro_nunocc-1, nspecies))

      ntotmax = ceiling( sqrt( 0.25d0 + wfpro_noccmax*(wfpro_noccmax+1) + 2.d0*wfpro_nunocc) - 0.5d0)
      allocate( nnl( 0:ntotmax, ntotmax, nspecies), nmin( nspecies))

      nnl = 0
      nmin = 1
      wfpro_nst = 0
      wfpro_states = 0
      wfpro_lmax = 0
      wfpro_nmax = 0

      ! get occupancy per state
      do is = 1, nspecies
        do ist = 1, spnst( is)
          n = spn( ist, is)
          l = spl( ist, is)
          if( spcore( ist, is)) cycle
          nnl( l, n, is) = nnl( l, n, is) + nint( spocc( ist, is))
        end do
        do n = 1, ntotmax
          if( maxval( nnl( :, n, is)) .gt. 0) exit
        end do
        nmin( is) = n
      end do

      ! select states
      do is = 1, nspecies
        nunocc = 0
        do nl = 1, 2*ntotmax-1
          lmax = ceiling( 0.5d0*nl) - 1
          do l = lmax, 0, -1
            n = nl - l
            if( (n .gt. ntotmax) .or. (n .lt. nmin( is))) cycle
            nemax = 2*(2*l+1)
            if( nunocc .lt. wfpro_nunocc) then
            !if( nnl( l, n, is) .eq. nemax) then
              wfpro_nst( is) = wfpro_nst( is) + 1
              wfpro_states( :, wfpro_nst( is), is) = (/n, l/)
              wfpro_lmax = max( wfpro_lmax, l)
              wfpro_nmax = max( wfpro_nmax, n)
              if( nnl( l, n, is) .eq. 0) nunocc = nunocc + 1
            end if
          end do
        end do
      end do

      wfpro_nprojtot = 0
      wfpro_nradfun = 0
      do is = 1, nspecies
        do ia = 1, natoms( is)
          do ist = 1, wfpro_nst( is)
            wfpro_nprojtot = wfpro_nprojtot + 2*wfpro_states( 2, ist, is) + 1
            wfpro_nradfun = wfpro_nradfun + 1
          end do
        end do
      end do
      wfpro_nprojtot = wfpro_dordmax*wfpro_nprojtot
      wfpro_nradfun = wfpro_dordmax*wfpro_nradfun
      !write(*,*) 'wfpro_nprojtot', wfpro_nprojtot
      !write(*,*) 'wfpro_nradfun', wfpro_nradfun

      if( allocated( wfpro_projst)) deallocate( wfpro_projst)
      allocate( wfpro_projst( 7, wfpro_nprojtot))
      if( allocated( wfpro_projused)) deallocate( wfpro_projused)
      allocate( wfpro_projused( wfpro_nprojtot))
        
      wfpro_projst = 0
      
      nl = 0
      lmax = 0
      do is = 1, nspecies
        do ia = 1, natoms( is)
          do n = 1, wfpro_dordmax
            do ist = 1, wfpro_nst( is)
              lmax = lmax + 1
              l = wfpro_states( 2, ist, is)
              do m = -l, l
                nl = nl + 1
                wfpro_projst( 1, nl) = is
                wfpro_projst( 2, nl) = ia
                wfpro_projst( 3, nl) = lmax                       ! idx of radial function
                wfpro_projst( 4, nl) = wfpro_states( 1, ist, is)  ! principal quantum number
                wfpro_projst( 5, nl) = l                          ! angular quantum number
                wfpro_projst( 6, nl) = m                          ! magnetic quantum number
                wfpro_projst( 7, nl) = n                          ! energy derivative order
              end do
            end do
          end do
        end do
      end do

      return
    end subroutine wfpro_getstates

    ! generates local orbital radial functions
    subroutine wfpro_genradfun
      integer :: is, ia, ias, ist, l, n, nr, io1, io2, ord(2), dord, np, pord
      integer :: iradfun, nn, ir, j, info
      real(8) :: line( 0:wfpro_nmax, 0:wfpro_lmax, nspecies), t1

      real(8) :: vr( nrmtmax), fr( nrmtmax), gr( nrmtmax), cf( 3, nrmtmax)
      real(8) :: p0( nrmtmax, 2), p1( nrmtmax, 2)
      real(8) :: q0( nrmtmax, 2), q1( nrmtmax, 2)
      real(8) :: p0s( nrmtmax), p1s( nrmtmax) ,q0s( nrmtmax), q1s( nrmtmax)
      real(8) :: hp0( nrmtmax)

      integer, allocatable :: ipiv(:)
      real(8), allocatable :: xa(:), ya(:)
      real(8), allocatable :: a(:,:), b(:), c(:)
! external functions
      real(8) :: polynom
      external polynom
      
      call wfpro_getline( wfpro_nmax, wfpro_lmax, line)

      ! generate radial functions (copied from genlofr)
      if( allocated( wfpro_radfun)) deallocate( wfpro_radfun)
      allocate( wfpro_radfun( nrmtmax, wfpro_nradfun))
      np = 4
      pord = 2
      allocate( ipiv( np))
      allocate( xa( np), ya( np), c( np))
      allocate( a( np, np), b( np))
      iradfun = 0
      do is = 1, nspecies
        nr = nrmt( is)
        do ia = 1, natoms( is)
          ias = idxas( ia, is)
          
          vr( 1:nr) = veffmt( 1, 1:nr, ias)*y00
          do dord = 1, wfpro_dordmax
            ord(1) = 0
            ord(2) = dord
            do ist = 1, wfpro_nst( is)
              a = 0.d0
              iradfun = iradfun + 1
              n = wfpro_states( 1, ist, is)
              l = wfpro_states( 2, ist, is)
              do io2 = 1, pord
                call rschroddme( ord( io2), l, 0, line( n-l-1, l, is), nr, spr( :, is), vr, nn, p0( :, io2), p1( :, io2), q0( :, io2), q1( :, io2))
                do ir = 1, nr
                  fr( ir) = p0( ir, io2)**2
                end do
                call fderiv( -1, nr, spr( :, is), fr, gr, cf)
                t1 = 1.d0/dsqrt( abs( gr( nr)))
                p0( 1:nr, io2) = t1*p0( 1:nr, io2)
                p1( 1:nr, io2) = t1*p1( 1:nr, io2)
                q0( 1:nr, io2) = t1*q0( 1:nr, io2)
                q1( 1:nr, io2) = t1*q1( 1:nr, io2)
                do j = 1, np
                  ir = nr - np + j
                  xa( j) = spr( ir, is)
                  ya( j) = p0( ir, io2)/spr( ir, is)
                end do
                do io1 = 1, pord
                  a( io1, io2) = polynom( io1-1, np, xa, ya, c, rmt( is))
                end do
              end do
              b = 0.d0
              b( pord) = 1.d0
              call dgesv( pord, 1, a, np, ipiv, b, np, info)
              if( info .ne. 0) then
                write(*,*)
                write(*,'("Error (wfpro_genradfun): DGESV returned non zero error code:",i4)') info
                stop
              end if
              p0s = 0.d0
              p1s = 0.d0
              q0s = 0.d0
              q1s = 0.d0
              do io1 = 1, pord
                t1 = b( io1)
                p0s( 1:nr) = p0s( 1:nr) + t1*p0( 1:nr, io1)
                p1s( 1:nr) = p1s( 1:nr) + t1*p1( 1:nr, io1)
                q0s( 1:nr) = q0s( 1:nr) + t1*q0( 1:nr, io1)
                q1s( 1:nr) = q1s( 1:nr) + t1*q1( 1:nr, io1)
              end do
              do ir = 1, nr
                fr( ir) = p0s( ir)**2
              end do
              call fderiv( -1, nr, spr( :, is), fr, gr, cf)
              t1 = 1.d0/dsqrt( abs( gr( nr)))
              p0s( 1:nr) = t1*p0s( 1:nr)
              p1s( 1:nr) = t1*p1s( 1:nr)
              q0s( 1:nr) = t1*q0s( 1:nr)
              q1s( 1:nr) = t1*q1s( 1:nr)
              call rschrodapp( l, nr, spr( :, is), vr, p0s, q0s, q1s, hp0)
              do ir = 1, nr
                t1 = 1.d0/spr( ir, is)
                wfpro_radfun( ir, iradfun) = t1*p0s( ir)
              end do
            end do
          end do

        end do
      end do
            
      return
    end subroutine wfpro_genradfun

    subroutine wfpro_projection
      use m_getunit
      integer :: ik, iproj, is, ias, io, ilo, l, m, lm, ig, ir, nr
      real(8) :: t0, t1
      real(8) :: fr( nrmtmax), gr( nrmtmax), cf( 3, nrmtmax)

      real(8), allocatable :: rolpi(:,:)
      complex(8), allocatable :: evecfv(:,:,:), apwalm(:,:,:,:,:), auxmat(:,:)

      if( input%properties%wannier%printproj) then
        call getunit( wfpro_un)
        open( wfpro_un, file=trim( wf_filename)//"_PROJECTION"//trim(filext), action='write', form='formatted')
      end if
        
      call timesec( t0)

      if( associated( input%properties%wannier%projection)) then
        wfpro_nunocc = input%properties%wannier%projection%nunocc
        wfpro_dordmax = max( 1, input%properties%wannier%projection%dordmax)
        wfpro_nproj = input%properties%wannier%projection%nprojtot
        wfpro_epsld = input%properties%wannier%projection%epsld
      end if

      call wfpro_neighcells
      call wfpro_getstates
      call wfpro_genradfun
      call wfpro_removeldprojf( wfpro_epsld)
            
      if( mpiglobal%rank .eq. 0) then
        write( wf_info, '(" calculate projection overlap matrices...")')
      end if
      !write(*,*) "rolpi"
      allocate( rolpi( wf_nprojtot, apwordmax+nlomax))
      rolpi(:,:) = 0.d0
      do iproj = 1, wf_nprojtot
        is = wf_projst( 1, iproj)
        ias = idxas( wf_projst( 2, iproj), is)
        nr = nrmt( is)
        l = wf_projst( 5, iproj)
        lm = idxlm( l, wf_projst( 6, iproj))
        do io = 1, apword( l, is)
          fr = 0.d0
          do ir = 1, nr
            fr( ir) = wfpro_radfun( ir, wf_projst( 3, iproj))*apwfr( ir, 1, io, l, ias)*spr( ir, is)**2
          end do
          call fderiv( -1, nr, spr( :, is), fr, gr, cf)
          rolpi( iproj, io) = gr( nr)
        end do
        do ilo = 1, nlorb( is)
          if( l .eq. lorbl( ilo, is)) then
            fr = 0.d0
            do ir = 1, nr
              fr( ir) = wfpro_radfun( ir, wf_projst( 3, iproj))*lofr( ir, 1, ilo, ias)*spr( ir, is)**2
            end do
            call fderiv( -1, nr, spr( :, is), fr, gr, cf)
            rolpi( iproj, apwordmax+ilo) = gr( nr)
          end if
        end do
      end do

      if( allocated( wfpro_proj)) deallocate( wfpro_proj)
      allocate( wfpro_proj( wf_fst:wf_lst, wf_nprojtot, wf_kset%nkpt))

      allocate( evecfv( nmatmax_ptr, nstfv, nspinor))
      allocate( apwalm( ngkmax_ptr, apwordmax, lmmaxapw, natmtot, nspinor))
      allocate( auxmat( nmatmax_ptr, wf_nprojtot))

      auxmat = zzero
      do ik = firstofset( mpiglobal%rank, wf_kset%nkpt), lastofset( mpiglobal%rank, wf_kset%nkpt)   
        call wfhelp_getevec( ik, evecfv)
        call match( wf_Gkset%ngk( 1, ik), wf_Gkset%gkc( :, 1, ik), wf_Gkset%tpgkc( :, :, 1, ik), wf_Gkset%sfacgk( :, :, 1, ik), apwalm(:, :, :, :, 1))

        ! projection matrices 
        auxmat(:,:) = zzero
        do iproj = 1, wf_nprojtot
          is = wf_projst( 1, iproj)
          ias = idxas( wf_projst( 2, iproj), is)
          lm = idxlm( wf_projst( 5, iproj), wf_projst( 6, iproj))
          do io = 1, apword( wf_projst( 5, iproj), is)
            auxmat( 1:wf_Gkset%ngk( 1, ik), iproj) = auxmat( 1:wf_Gkset%ngk( 1, ik), iproj) + &
                conjg( apwalm( 1:wf_Gkset%ngk( 1, ik), io, lm, ias, 1))*cmplx( rolpi( iproj, io), 0.d0, 8)
          end do
          do ilo = 1, nlorb( is)
            l = lorbl( ilo, is)
            if( l .eq. wf_projst( 5, iproj)) then
              do m = -l, l
                if( m .eq. wf_projst( 6, iproj)) then
                  ig = wf_Gkset%ngk( 1, ik) + idxlo( lm, ilo, ias)
                  auxmat( ig, iproj) = cmplx( rolpi( iproj, apwordmax+ilo), 0.d0, 8)
                end if
              end do
            end if
          end do
        end do

        call zgemm( 'c', 'n', wf_nst, wf_nprojtot, wf_Gkset%ngk( 1, ik)+nlotot, zone, &
               evecfv( 1, wf_fst, 1), nmatmax_ptr, &
               auxmat, nmatmax_ptr, zzero, &
               wfpro_proj( wf_fst, 1, ik), wf_nst)

      end do 
      call mpi_allgatherv_ifc( set=wf_kset%nkpt, rlen=wf_nst*wf_nprojtot, zbuf=wfpro_proj)
      call barrier

      deallocate( evecfv, apwalm, auxmat)
      deallocate( rolpi)

      call wfpro_selectprojf

      call timesec( t1)
      if( mpiglobal%rank .eq. 0) then
        write( wf_info, '(5x,"duration (seconds): ",T40,3x,F10.1)') t1-t0
        write( wf_info, '(5x,"#k-points: ",T40,7x,I6)') wf_kset%nkpt
        write( wf_info, '(5x,"#states: ",T40,7x,I6)') wf_nst
        write( wf_info, '(5x,"#projection functions: ",T40,7x,I6)') wf_nprojtot
        write( wf_info, *)
        call flushifc( wf_info)
      end if

      if( allocated( wfpro_nst     )) deallocate( wfpro_nst     )
      if( allocated( wfpro_states  )) deallocate( wfpro_states  )
      if( allocated( wfpro_projst  )) deallocate( wfpro_projst  )
      if( allocated( wfpro_projused)) deallocate( wfpro_projused)
      if( allocated( wfpro_nn      )) deallocate( wfpro_nn      )
      if( allocated( wfpro_vn      )) deallocate( wfpro_vn      )
      if( allocated( wfpro_radfun  )) deallocate( wfpro_radfun  )
      if( allocated( wfpro_proj    )) deallocate( wfpro_proj    )

      !write(*,*) "projection done"
      if( input%properties%wannier%printproj) then
        call wfpro_writepro_lo( wfpro_un)
        close( wfpro_un)
      end if

      return
    end subroutine wfpro_projection

    subroutine wfpro_removeldprojf( eps)
      use m_linalg
      real(8), intent( in), optional :: eps

      integer :: igroup, i, is, ias, lm, j, ias2, lm2, iproj, iproj2, ir, nr
      real(8) :: tmp, eps_
      real(8) :: fr( nrmtmax), gr( nrmtmax), cf( 3, nrmtmax)

      real(8), allocatable :: projolp(:,:), epslist(:), auxmat(:,:)

      eps_ = 1.d-3
      if( present( eps)) eps_ = eps

      ! overlap between all projection functions
      allocate( projolp( wfpro_nprojtot, wfpro_nprojtot))
      projolp(:,:) = 0.d0
      do iproj = 1, wfpro_nprojtot
        projolp( iproj, iproj) = 1.d0
        is = wfpro_projst( 1, iproj)
        ias = idxas( wfpro_projst( 2, iproj), is)
        nr = nrmt( is)
        lm = idxlm( wfpro_projst( 5, iproj), wfpro_projst( 6, iproj))
        do iproj2 = iproj+1, wfpro_nprojtot
          ias2 = idxas( wfpro_projst( 2, iproj2), wfpro_projst( 1, iproj2))
          lm2 = idxlm( wfpro_projst( 5, iproj2), wfpro_projst( 6, iproj2))
          if( (lm .eq. lm2) .and. (ias .eq. ias2)) then
            fr = 0.d0
            do ir = 1, nr
              fr( ir) = wfpro_radfun( ir, wfpro_projst( 3, iproj))*wfpro_radfun( ir, wfpro_projst( 3, iproj2))*spr( ir, is)**2
            end do
            call fderiv( -1, nr, spr( :, is), fr, gr, cf)
            projolp( iproj, iproj2) = gr( nr)
            projolp( iproj2, iproj) = gr( nr)
          end if
        end do
      end do

      ! find linearly independent ones (using QR decomposition)
      allocate( auxmat( wfpro_nprojtot, wfpro_nprojtot))
      allocate( epslist( wfpro_nprojtot))
      call rqr( projolp, r=auxmat)
      tmp = 0.d0
      do iproj = 1, wfpro_nprojtot
        tmp = max( tmp, abs( auxmat( iproj, iproj)))
      end do
      
      epslist = -1.d0
      do iproj = 1, wfpro_nprojtot
        epslist( iproj) = abs( auxmat( iproj, iproj))/tmp
      end do
      deallocate( auxmat)

      wf_nprojtot = 0
      wfpro_nprojtot_n = 0
      wfpro_projused = 0
      if( wfpro_nproj .eq. 0) wfpro_nproj = wfpro_nprojtot
      do iproj = 1, wfpro_nprojtot
        nr = maxloc( epslist, 1)
        if( (epslist( nr) .gt. eps_) .and. (wf_nprojtot .lt. wfpro_nproj)) then
          wf_nprojtot = wf_nprojtot + 1
          wfpro_nprojtot_n = wfpro_nprojtot_n + &
            maxval( wfpro_nn( idxas( wfpro_projst(2,iproj), wfpro_projst(1,iproj)), :))
          wfpro_projused( nr) = 1
        end if
        epslist( nr) = -1.d0
      end do

      if( allocated( wf_projolp)) deallocate( wf_projolp)
      allocate( wf_projolp( wf_nprojtot, wf_nprojtot))
      if( allocated( wf_projst)) deallocate( wf_projst)
      allocate( wf_projst( 7, wf_nprojtot))

      wf_projolp = 0.d0
      i = 1
      do iproj = 1, wfpro_nprojtot
        if( wfpro_projused( iproj) .eq. 0) cycle
        wf_projst( 1:7, i) = wfpro_projst( :, iproj)
        j = 1
        do iproj2 = 1, wfpro_nprojtot
          if( wfpro_projused( iproj2) .eq. 0) cycle
          wf_projolp( i, j) = projolp( iproj, iproj2)
          j = j + 1
        end do
        i = i + 1
      end do

      deallocate( projolp, epslist)

      do igroup = 1, wf_ngroups
        if( allocated( wf_groups( igroup)%projused)) deallocate( wf_groups( igroup)%projused)
        allocate( wf_groups( igroup)%projused( wf_nprojtot))
      end do

      return
    end subroutine wfpro_removeldprojf

    ! construct a set of orthonormal projection functions from local orbitals
    subroutine wfpro_selectprojf
      use m_linalg, only: zhegdiag
      integer :: igroup, nproj1, nproj2, iproj1, iproj2, ik, i, j, k, l, n, fst, lst, nwf, is, ia, ias, nadd, in, iat
      real(8) :: s

      integer, allocatable :: npat(:), apdat(:,:,:)
      real(8), allocatable :: eval(:), apchar(:,:,:)
      complex(8), allocatable :: zmat(:,:), pmat(:,:,:), smat(:,:), omat(:,:), cmat(:,:) 
      complex(8), allocatable :: auxmat(:,:,:,:), aux1(:,:), aux2(:,:)

      if( input%properties%wannier%printproj) then
        allocate( apchar( 0:4, wfpro_nprojtot_n, wf_ngroups))
        allocate( apdat( 5, wfpro_nprojtot_n, wf_ngroups))
        apchar = 0.d0
      end if
      do igroup = 1, wf_ngroups
        
        nproj1 = size( input%properties%wannier%grouparray( igroup)%group%projectorarray, 1)
        nproj2 = input%properties%wannier%grouparray( igroup)%group%nproj
        if( nproj2 .eq. 0) then
          wf_groups( igroup)%nproj = wf_nprojtot
        else
          wf_groups( igroup)%nproj = nproj2
        end if
        if( nproj1 .gt. 0) then
          wf_groups( igroup)%nproj = nproj1
        end if
        wf_groups( igroup)%projused = 0

        ! if special projectors are defined use only them
        if( (wf_groups( igroup)%method .eq. "pro") .or. (wf_groups( igroup)%method .eq. "promax")) then
          if( wf_groups( igroup)%nproj .ne. wf_groups( igroup)%nwf) then
            if( mpiglobal%rank .eq. 0) then
              write(*,*)
              write( *, '("Error (wfpro_selectprojf): The number of projectors must be equal to the number of Wannier functions for group ",I2,".")') igroup
            end if
            stop
          end if
          do iproj1 = 1, nproj1
            if( (input%properties%wannier%grouparray( igroup)%group%projectorarray( iproj1)%projector%nr .lt. 1) .or. &
                (input%properties%wannier%grouparray( igroup)%group%projectorarray( iproj1)%projector%nr .gt. wf_nprojtot)) then
              if( mpiglobal%rank .eq. 0) then
                write(*,*)
                write( *, '("Error (wfpro_selectprojf): ",I4," is not a valid index for projection local-orbitals in group ",I2,".")') &
                    input%properties%wannier%grouparray( igroup)%group%projectorarray( iproj1)%projector%nr, igroup
                write(*, '(" Here is a list of local-orbitals that can be used for projection:")')
                call wfpro_writepro_lo( 6)
              end if
              stop
            end if
            wf_groups( igroup)%projused( input%properties%wannier%grouparray( igroup)%group%projectorarray( iproj1)%projector%nr) = 1
          end do
        else if( wf_groups( igroup)%method .eq. 'opf' .or. &
                 wf_groups( igroup)%method .eq. 'opfmax' .or. &
                 wf_groups( igroup)%method .eq. 'disSMV' .or. &
                 wf_groups( igroup)%method .eq. 'disFull') then
          if( wf_groups( igroup)%nproj .lt. wf_groups( igroup)%nwf) then
            if( mpiglobal%rank .eq. 0) then
              write(*,*)
              write( *, '("Error (wfpro_selectprojf): The number of projectors must be greater or equal to the number of Wannier functions for group ",I2,".")') igroup
            end if
            stop
          end if
          do iproj1 = 1, nproj1
            if( (input%properties%wannier%grouparray( igroup)%group%projectorarray( iproj1)%projector%nr .lt. 1) .or. &
                (input%properties%wannier%grouparray( igroup)%group%projectorarray( iproj1)%projector%nr .gt. wf_nprojtot)) then
              if( mpiglobal%rank .eq. 0) then
                write(*,*)
                write( *, '("Error (wfpro_selectprojf): ",I4," is not a valid index for projection local-orbitals in group ",I2,".")') &
                    input%properties%wannier%grouparray( igroup)%group%projectorarray( iproj1)%projector%nr, igroup
                write(*, '(" Here is a list of local-orbitals that can be used for projection:")')
                call wfpro_writepro_lo( 6)
              end if
              stop
            end if
            wf_groups( igroup)%projused( input%properties%wannier%grouparray( igroup)%group%projectorarray( iproj1)%projector%nr) = 1
          end do
        end if

        ! construct set of orthonormal atomic projectors
        if( nproj1 .le. 0) then
          wf_groups( igroup)%nproj = 0
          wf_groups( igroup)%projused = 1
          fst = wf_groups( igroup)%fst
          lst = wf_groups( igroup)%lst
          nwf = wf_groups( igroup)%nwf

          nadd = 0
          do is = 1, nspecies
            do ia = 1, natoms( is)
              n = 0
              do i = 1, wf_nprojtot
                if( wf_groups( igroup)%projused( i) .eq. 0 .or. &
                    wf_projst( 1, i) .ne. is .or. &
                    wf_projst( 2, i) .ne. ia) cycle
                n = n+1
              end do
              nadd = max( nadd, n)
            end do
          end do
          allocate( npat( sum( wfpro_nn(:,igroup))))
          allocate( auxmat( fst:lst, nadd, wf_kset%nkpt, sum( wfpro_nn(:,igroup))))
          auxmat = zzero

          iat = 0
          do is = 1, nspecies
            do ia = 1, natoms( is)
              ias = idxas( ia, is)

              n = 0
              do i = 1, wf_nprojtot
                if( wf_groups( igroup)%projused( i) .eq. 0 .or. &
                    wf_projst( 1, i) .ne. is .or. &
                    wf_projst( 2, i) .ne. ia) cycle
                n = n+1
              end do

              allocate( zmat(n,n), pmat( fst:lst, n, wf_kset%nkpt), smat(nwf,n), omat(n,n), cmat(n,n), eval(n))
              allocate( aux1(n,n), aux2(n,n))
              do in = 1, wfpro_nn( ias, igroup)
                iat = iat + 1
                ! projection on (rotated) KS states
                zmat = zzero
                pmat = zzero
                do ik = 1, wf_kset%nkpt
                  s = -twopi*dot_product( wf_kset%vkl(:,ik), dble( wfpro_vn(:,in,ias,igroup)))
                  j = 0
                  do i = 1, wf_nprojtot
                    if( wf_groups( igroup)%projused( i) .eq. 0 .or. &
                        wf_projst( 1, i) .ne. is .or. &
                        wf_projst( 2, i) .ne. ia) cycle
                    j = j + 1
                    pmat( :, j, ik) = cmplx( cos( s), sin( s), 8)*wfpro_proj( fst:lst, i, ik)
                  end do
                  call zgemm( 'c', 'n', nwf, n, lst-fst+1, zone, &
                         wf_transform( fst, wf_groups( igroup)%fwf, ik), wf_nst, &
                         pmat( fst, 1, ik), lst-fst+1, zzero, &
                         smat, nwf)
                  call zgemm( 'c', 'n', n, n, nwf, zone, &
                         smat, nwf, &
                         smat, nwf, zone, &
                         zmat, n)
                end do
                zmat = zmat/dble( wf_kset%nkpt)

                ! projection function overlaps
                i = 0
                do k = 1, wf_nprojtot
                  if( wf_groups( igroup)%projused( k) .eq. 0 .or. &
                      wf_projst( 1, k) .ne. is .or. &
                      wf_projst( 2, k) .ne. ia) cycle
                  i = i + 1
                  j = 0
                  do l = 1, wf_nprojtot
                    if( wf_groups( igroup)%projused( l) .eq. 0 .or. &
                        wf_projst( 1, l) .ne. is .or. &
                        wf_projst( 2, l) .ne. ia) cycle
                    j = j + 1
                    omat(j,i) = cmplx( wf_projolp(l,k), 0.d0, 8)
                  end do
                end do

                ! construct projection matrix
                call zhegdiag( zmat, omat, eval, evec=cmat)

                npat( iat) = 0
                if( nproj2 .eq. 0) then
                  do while( .true.)
                    if( npat( iat) .ge. n) exit
                    if( eval( n-npat( iat)) .le. input%properties%wannier%grouparray( igroup)%group%epsproj) exit
                    npat( iat) = npat( iat) + 1
                  end do
                else
                  npat( iat) = nproj2
                end if
                npat( iat) = max( 1, npat( iat))

                do ik = 1, wf_kset%nkpt
                  call zgemm( 'n', 'n', lst-fst+1, npat( iat), n, zone, &
                         pmat( fst, 1, ik), lst-fst+1, &
                         cmat( 1, n-npat( iat)+1), n, zzero, &
                         auxmat( fst, 1, ik, iat), lst-fst+1)
                end do

                ! compute atomic projector angular character
                if( input%properties%wannier%printproj) then
                  do i = 1, npat( iat)
                    apdat( 1, wf_groups( igroup)%nproj+i, igroup) = is
                    apdat( 2, wf_groups( igroup)%nproj+i, igroup) = ia
                    apdat( 3:5, wf_groups( igroup)%nproj+i, igroup) = wfpro_vn(:,in,ias,igroup)
                    j = 0
                    do iproj1 = 1, wf_nprojtot
                      if( wf_groups( igroup)%projused( iproj1) .eq. 0 .or. &
                          wf_projst( 1, iproj1) .ne. is .or. &
                          wf_projst( 2, iproj1) .ne. ia) cycle
                      j = j+1
                      l = wf_projst( 5, iproj1)
                      if( l .gt. 4) cycle
                      k = 0
                      do iproj2 = 1, wf_nprojtot
                        if( wf_groups( igroup)%projused( iproj2) .eq. 0 .or. &
                            wf_projst( 1, iproj2) .ne. is .or. &
                            wf_projst( 2, iproj2) .ne. ia) cycle
                        k = k+1
                        if( wf_projst( 5, iproj2) .ne. l) cycle
                        apchar( l, wf_groups( igroup)%nproj+i, igroup) = apchar( l, wf_groups( igroup)%nproj+i, igroup) + &
                          dble( conjg( cmat( j, n-npat( iat)+i))*cmat( k, n-npat( iat)+i)*omat(j,k))
                      end do
                    end do
                  end do
                end if

                wf_groups( igroup)%nproj = wf_groups( igroup)%nproj + npat( iat)
              end do
              deallocate( zmat, pmat, smat, omat, cmat, eval)
              deallocate( aux1, aux2)
            end do
          end do

          if( wf_groups( igroup)%nproj .lt. wf_groups( igroup)%nwf) then
            if( mpiglobal%rank .eq. 0) then
              write(*,*)
              write(*,'("Error (wfpro_selectprojf): The number of projection functions is smaller than the number of Wannier functions for group ",i2,".")') igroup
              write(*,'(4x,"projection functions: ",i6)') wf_groups( igroup)%nproj
              write(*,'(4x,"Wannier functions:    ",i6)') wf_groups( igroup)%nwf
              write(*,'("Consider the following options:")')
              write(*,'(4x,"- smaller eigenvalue cutoff (attribute epsproj)")')
              write(*,'(4x,"- larger number of unoccpied states per atom (attribute nunocc)")')
              write(*,'(4x,"- include atoms in neighboring cells (attribute neighcells)")')
            end if
            stop
          end if
          if( allocated( wf_groups( igroup)%projection)) deallocate( wf_groups( igroup)%projection)
          allocate( wf_groups( igroup)%projection( fst:lst, wf_groups( igroup)%nproj, wf_kset%nkpt))
          nadd = 0
          do iat = 1, sum( wfpro_nn(:,igroup))
            do ik = 1, wf_kset%nkpt
              wf_groups( igroup)%projection( :, (nadd+1):(nadd+npat( iat)), ik) = auxmat( :, 1:npat( iat), ik, iat)
            end do
            nadd = nadd + npat( iat)
          end do
          deallocate( npat, auxmat)
        
        ! if special projectors are given, use them
        else
          if( allocated( wf_groups( igroup)%projection)) deallocate( wf_groups( igroup)%projection)
          allocate( wf_groups( igroup)%projection( fst:lst, wf_groups( igroup)%nproj, wf_kset%nkpt))
          j = 0
          do i = 1, wf_nprojtot
            if( wf_groups( igroup)%projused( i) .eq. 0) cycle
            j = j + 1
            do ik = 1, wf_kset%nkpt
              wf_groups( igroup)%projection( :, j, ik) = wfpro_proj( fst:lst, i, ik)
            end do
          end do
        end if
      end do

      if( input%properties%wannier%printproj) then
        call wfpro_writepro_ap( wfpro_un, apdat, apchar)
        deallocate( apdat, apchar)
      end if

      return
    end subroutine wfpro_selectprojf

    ! print out projection local-orbitals
    subroutine wfpro_writepro_lo( un)
      integer, intent( in) :: un
      
      ! local variables
      integer :: iproj, igroup

      call printbox( un, '*', "Local-orbitals for projection")
      write( un, *)
      write( un, '(80("-"))')
      write( un, '(7x,"#",6x,"atom",4x,"  n /  l /  m",4x,"dord",4x,"used")')
      write( un, '(80("-"))')
      do iproj = 1, wf_nprojtot
        write( un, '(4x,i4,4x,a2,x,i3,4x,2(i3," /"),i3,4x,i4,4x,a1)', advance='no') &
            iproj, &
            spsymb( wf_projst( 1, iproj)), &
            wf_projst( 2, iproj), &
            wf_projst( 4:7, iproj)
        do igroup = 1, wf_ngroups
          if( wf_groups( igroup)%projused( iproj) .eq. 1) then
            write( un, '(i2,x)', advance='no') igroup
          else
            write( un, '(3x)', advance='no')
          end if
        end do
        write( un, *)
      end do
      write( un, *)
      return
    end subroutine wfpro_writepro_lo

    ! print out atomic projectors
    subroutine wfpro_writepro_ap( un, apdat, apchar)
      integer, intent( in) :: un, apdat( 5, wfpro_nprojtot_n, wf_ngroups)
      real(8), intent( in) :: apchar( 0:4, wfpro_nprojtot_n, wf_ngroups)

      integer :: igroup, iproj

      call printbox( un, '*', "Atomic projection functions")
      write( un, *)
      write( un, '(80("-"))')
      write( un, '(7x,"#",6x,"atom",10x,"cell",4x,2x,"s",4x,"p",4x,"d",4x,"f",4x,"g")')
      write( un, '(80("-"))')
      do igroup = 1, wf_ngroups
        do iproj = 1, wf_groups( igroup)%nproj
          write( un, '(4x,i4,4x,a2,x,i3,4x,"(",2(i2,","),i2,")",4x,4(i3,2x),i3)') &
            iproj, &
            spsymb( apdat( 1, iproj, igroup)), &
            apdat( 2, iproj, igroup), &
            apdat( 3:5, iproj, igroup), &
            nint( 1.d2*apchar( :, iproj, igroup))
        end do
        write( un, '(80("-"))')
      end do
      write( un, *)
      return
    end subroutine wfpro_writepro_ap

end module mod_wannier_projection
