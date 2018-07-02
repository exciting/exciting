module mod_wannier_projection
  use mod_wannier_variables
  use mod_wannier_helper

  use mod_atoms
  use mod_eigensystem
  use mod_APW_LO
  use mod_muffin_tin
  use mod_eigenvalue_occupancy
  use mod_constants
  use mod_Gkvector,              only: ngkmax_ptr
  use mod_spin,                  only: nspinor
  use mod_potential_and_density, only: veffmt

  implicit none

  integer :: wfpro_noccmax = 10
  integer :: wfpro_nunocc = 20
  integer :: wfpro_dordmax = 2
  integer :: wfpro_nproj = 0
  real(8) :: wfpro_epsld = 1.d-3

  integer :: wfpro_nradfun
  integer :: wfpro_nprojtot
  integer :: wfpro_lmax
  
  integer, allocatable :: wfpro_nst(:), wfpro_states(:,:,:), wfpro_projst(:,:), wfpro_projused(:)
  real(8), allocatable :: wfpro_radfun(:,:)
  real(8), allocatable :: wfpro_radint(:,:)

! methods
  contains

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
        vr( 1:nr) = veffmt( 1, 1:nr, ias)*y00
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
                elo = 2*ens(0) - ens(1) ! assuming lowest eigenenergy is  negative
              else
                elo = ens( n-1)
              endif
              call rschroddme( 0, l, 0, elo, nr, spr( :, is), vr, nn, p0s, hp0, q0s, q1s) 
              flo = hp0( nr)
              if( ehi .lt. elo) then
                write(*,*) 'warning'
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

    subroutine wfpro_getstates
      integer :: is, ia, ist, l, lmax, m, n, nl, nemax, nunocc, nnl( 2*wfpro_noccmax+wfpro_nunocc-1, nspecies)

      allocate( wfpro_nst( nspecies))
      allocate( wfpro_states( 2, 2*wfpro_noccmax+wfpro_nunocc-1, nspecies))

      nnl = 0
      wfpro_nst = 0
      wfpro_states = 0
      wfpro_lmax = 0

      ! get occupancy per state
      do is = 1, nspecies
        do ist = 1, spnst( is)
          n = spn( ist, is)
          l = spl( ist, is)
          nnl( n+l, is) = nnl( n+l, is) + nint( spocc( ist, is))
        end do
      end do

      ! select fully occupied states
      do is = 1, nspecies
        do n = 1, wfpro_noccmax
          do l = 0, n-1
            nl = n + l
            nemax = 2*ceiling( 0.5d0*nl)**2
            do ist = 1, spnst( is)
              if( (n .eq. spn( ist, is)) .and. (l .eq. spl( ist, is)) .and..not. spcore( ist, is) .and. (nnl( nl, is)/nemax .eq. 1)) then
                wfpro_nst( is) = wfpro_nst( is) + 1
                wfpro_states( :, wfpro_nst( is), is) = (/n, l/)
                wfpro_lmax = max( wfpro_lmax, l)
                exit
              end if
            end do
          end do
        end do
      end do

      ! select partially occupied and unoccupied states
      do is = 1, nspecies
        nunocc = 0
        do nl = 1, 2*wfpro_noccmax+wfpro_nunocc-1
          lmax = ceiling( 0.5d0*nl) - 1
          nemax = 2*ceiling( 0.5d0*nl)**2
          if( nnl( nl, is)/nemax .lt. 1) then
            m = 0
            do l = lmax, 0, -1
              if( nunocc .lt. wfpro_nunocc) then
                wfpro_nst( is) = wfpro_nst( is) + 1
                if( nnl( nl, is) .le. m) nunocc = nunocc + 1
                wfpro_states( :, wfpro_nst( is), is) = (/nl-l, l/)
                wfpro_lmax = max( wfpro_lmax, l)
              end if
              m = m + 2*(2*l + 1)
            end do
          end if
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

      allocate( wfpro_projst( 5, wfpro_nprojtot))
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
                wfpro_projst( 3, nl) = lmax ! idx of radial function
                wfpro_projst( 4, nl) = l
                wfpro_projst( 5, nl) = m
              end do
            end do
          end do
        end do
      end do

      return
    end subroutine wfpro_getstates

    subroutine wfpro_genradfun
      integer :: is, ia, ias, ist, l, n, nr, io1, io2, ord(2), dord, np
      integer :: iradfun, nn, ir, j, info
      real(8) :: line( 0:(2*wfpro_noccmax+wfpro_nunocc-1), 0:wfpro_lmax, nspecies), t1

      real(8) :: vr( nrmtmax), fr( nrmtmax), gr( nrmtmax), cf( 3, nrmtmax)
      real(8) :: p0( nrmtmax, maxlorbord), p1( nrmtmax, maxlorbord)
      real(8) :: q0( nrmtmax, maxlorbord), q1( nrmtmax, maxlorbord)
      real(8) :: p0s( nrmtmax), p1s( nrmtmax) ,q0s( nrmtmax), q1s( nrmtmax)
      real(8) :: hp0( nrmtmax)

      integer, allocatable :: ipiv(:)
      real(8), allocatable :: xa(:), ya(:)
      real(8), allocatable :: a(:,:), b(:), c(:)
! external functions
      real(8) :: polynom
      external polynom
      
      call wfpro_getline( (2*wfpro_noccmax+wfpro_nunocc-1), wfpro_lmax, line)

      ! generate radial functions (copied from genlofr)
      allocate( wfpro_radfun( nrmtmax, wfpro_nradfun))
      np = 5
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
            ord = (/0, dord/)
            do ist = 1, wfpro_nst( is)
              iradfun = iradfun + 1
              n = wfpro_states( 1, ist, is)
              l = wfpro_states( 2, ist, is)
              do io2 = 1, 2
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
                do io1 = 1, 2
                  a( io1, io2) = polynom( io1-1, np, xa, ya, c, rmt( is))
                end do
              end do
              b = 0.d0
              b(2) = 1.d0
              call dgesv( 2, 1, a, np, ipiv, b, np, info)
              if( info .ne. 0) then
                write(*,*) 'Ooops!'
                stop
              end if
              p0s = 0.d0
              p1s = 0.d0
              q0s = 0.d0
              q1s = 0.d0
              do io1 = 1, 2
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
              if (input%groundstate%SymmetricKineticEnergy) then
                do ir = 1, nr
                  t1 = 1.d0/spr( ir, is)
                  wfpro_radfun( ir, iradfun) = t1*p0s( ir)
                end do
              else
                call rschrodapp( l, nr, spr( :, is), vr, p0s, q0s, q1s, hp0)
                do ir = 1, nr
                  t1 = 1.d0/spr( ir, is)
                  wfpro_radfun( ir, iradfun) = t1*p0s( ir)
                end do
              end if
            end do
          end do

        end do
      end do
            
      return
    end subroutine wfpro_genradfun
    
    ! print out projection local-orbitals
    subroutine wfpro_showproj
      ! local variables
      integer :: iproj

      write(*,*)
      write(*,*) 'Local-orbitals for generating Wannier functions via projection.'
      write(*,*) '     nr   species   atom   l    m'
      do iproj = 1, wf_nprojtot
        write( *, '(5x,i3,3x)', advance='no') iproj
        write( *, '(2x,a5,3x)', advance='no') spsymb( wf_projst( 1, iproj))
        write( *, '(i4,3x)', advance='no') wf_projst( 2, iproj)
        write( *, '(i2,3x)', advance='no') wf_projst( 4, iproj)
        write( *, '(i2,3x)', advance='no') wf_projst( 5, iproj)
        write(*,*)
      end do
      write(*,*)
      return
    end subroutine wfpro_showproj

    subroutine wfpro_projection

      integer :: ik, iproj, ia, is, ias, io, ilo, l, m, lm, ig, i, ir, nr
      real(8) :: t0, t1
      real(8) :: fr( nrmtmax), gr( nrmtmax), cf( 3, nrmtmax)

      real(8), allocatable :: rolpi(:,:)
      complex(8), allocatable :: evecfv(:,:,:), apwalm(:,:,:,:,:), auxmat(:,:)

      call timesec( t0)

      if( associated( input%properties%wannier%projection)) then
        wfpro_nunocc = input%properties%wannier%projection%nunocc
        wfpro_dordmax = input%properties%wannier%projection%dordmax
        wfpro_nproj = input%properties%wannier%projection%nprojtot
        wfpro_epsld = input%properties%wannier%projection%epsld
      end if

      call wfpro_getstates
      call wfpro_genradfun
      call wfpro_removeldprojf( wfpro_epsld)
            
      write( wf_info, '(" calculate projection overlap matrices...")')
      allocate( rolpi( wf_nprojtot, apwordmax+nlomax))
      rolpi(:,:) = 0.d0
      do iproj = 1, wf_nprojtot
        is = wf_projst( 1, iproj)
        ias = idxas( wf_projst( 2, iproj), is)
        nr = nrmt( is)
        l = wf_projst( 4, iproj)
        lm = idxlm( l, wf_projst( 5, iproj))
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

      if( allocated( wf_projection)) deallocate( wf_projection)
      allocate( wf_projection( wf_fst:wf_lst, wf_nprojtot, wf_kset%nkpt))

      allocate( evecfv( nmatmax_ptr, nstfv, nspinor))
      allocate( apwalm( ngkmax_ptr, apwordmax, lmmaxapw, natmtot, nspinor))
      allocate( auxmat( nmatmax_ptr, wf_nprojtot))

      do ik = 1, wf_kset%nkpt   
        call wannier_getevec( ik, evecfv)
        call match( wf_Gkset%ngk( 1, ik), wf_Gkset%gkc( :, 1, ik), wf_Gkset%tpgkc( :, :, 1, ik), wf_Gkset%sfacgk( :, :, 1, ik), apwalm(:, :, :, :, 1))

        ! projection matrices 
        auxmat(:,:) = zzero
        do iproj = 1, wf_nprojtot
          is = wf_projst( 1, iproj)
          ias = idxas( wf_projst( 2, iproj), is)
          lm = idxlm( wf_projst( 4, iproj), wf_projst( 5, iproj))
          do io = 1, apword( wf_projst( 4, iproj), is)
            auxmat( 1:wf_Gkset%ngk( 1, ik), iproj) = auxmat( 1:wf_Gkset%ngk( 1, ik), iproj) + &
                conjg( apwalm( 1:wf_Gkset%ngk( 1, ik), io, lm, ias, 1))*cmplx( rolpi( iproj, io), 0, 8)
          end do
          do ilo = 1, nlorb( is)
            l = lorbl( ilo, is)
            if( l .eq. wf_projst( 4, iproj)) then
              do m = -l, l
                if( m .eq. wf_projst( 5, iproj)) then
                  ig = wf_Gkset%ngk( 1, ik) + idxlo( lm, ilo, ias)
                  auxmat( ig, iproj) = cmplx( rolpi( iproj, apwordmax+ilo), 0, 8)
                end if
              end do
            end if
          end do
        end do

        call zgemm( 'c', 'n', wf_nst, wf_nprojtot, wf_Gkset%ngk( 1, ik)+nlotot, zone, &
               evecfv( :, wf_fst:wf_lst, 1), nmatmax_ptr, &
               auxmat, nmatmax_ptr, zzero, &
               wf_projection( :, :, ik), wf_nst)

        !write(*,'(I,3F13.6)') ik, wf_kset%vkl( :, ik)
        !call plotmat( wf_projection( :, :, ik))
        !write(*,*)
      end do 

      deallocate( evecfv, apwalm, auxmat)
      deallocate( rolpi)

      call wfpro_selectprojf

      !call writematlab( wf_projection(:,:,12), 'pmfull')
      call timesec( t1)
      write( wf_info, '(5x,"duration (seconds): ",T40,3x,F10.1)') t1-t0
      write( wf_info, '(5x,"#k-points: ",T40,7x,I6)') wf_kset%nkpt
      write( wf_info, '(5x,"#states: ",T40,7x,I6)') wf_nst
      write( wf_info, '(5x,"#projection functions: ",T40,7x,I6)') wf_nprojtot
      write( wf_info, *)
      call flushifc( wf_info)

      return
    end subroutine wfpro_projection

    subroutine wfpro_removeldprojf( eps)
      real(8), intent( in), optional :: eps

      integer :: igroup, i, is, ias, lm, j, ias2, lm2, iproj, iproj2, ir, nr, lapack_info, lapack_lwork
      real(8) :: tmp, eps_
      real(8) :: fr( nrmtmax), gr( nrmtmax), cf( 3, nrmtmax)

      real(8), allocatable :: projolp(:,:), epslist(:)
      real(8), allocatable :: lapack_tau(:), lapack_work(:)

      eps_ = 1.d-3
      if( present( eps)) eps_ = eps

      allocate( projolp( wfpro_nprojtot, wfpro_nprojtot))
      allocate( epslist( wfpro_nprojtot))
      projolp(:,:) = 0.d0
      i = 0
      do iproj = 1, wfpro_nprojtot
        i = i + 1
        is = wfpro_projst( 1, iproj)
        ias = idxas( wfpro_projst( 2, iproj), is)
        nr = nrmt( is)
        lm = idxlm( wfpro_projst( 4, iproj), wfpro_projst( 5, iproj))
        j = 0
        do iproj2 = 1, wfpro_nprojtot
          j = j + 1
          ias2 = idxas( wfpro_projst( 2, iproj2), wfpro_projst( 1, iproj2))
          lm2 = idxlm( wfpro_projst( 4, iproj2), wfpro_projst( 5, iproj2))
          if( (lm .eq. lm2) .and. (ias .eq. ias2)) then
            fr = 0.d0
            do ir = 1, nr
              fr( ir) = wfpro_radfun( ir, wfpro_projst( 3, iproj))*wfpro_radfun( ir, wfpro_projst( 3, iproj2))*spr( ir, is)**2
            end do
            call fderiv( -1, nr, spr( :, is), fr, gr, cf)
            projolp( i, j) = gr( nr)
          end if
        end do
      end do
      allocate( lapack_tau( wfpro_nprojtot))
      allocate( lapack_work( 1))
      call dgeqrf( wfpro_nprojtot, wfpro_nprojtot, projolp, wfpro_nprojtot, lapack_tau, lapack_work, -1, lapack_info)
      lapack_lwork = nint( lapack_work( 1))
      deallocate( lapack_work)
      allocate( lapack_work( lapack_lwork))
      !call writematlab( cmplx( projolp, 0, 8), 'olp')
      call dgeqrf( wfpro_nprojtot, wfpro_nprojtot, projolp, wfpro_nprojtot, lapack_tau, lapack_work, lapack_lwork, lapack_info)
      !call writematlab( cmplx( projolp, 0, 8), 'qr')
      tmp = 0.d0
      do iproj = 1, wfpro_nprojtot
        tmp = max( tmp, abs( projolp( iproj, iproj)))
      end do
      
      epslist = -1.d0
      do iproj = 1, wfpro_nprojtot
        epslist( iproj) = abs( projolp( iproj, iproj))/tmp
      end do

      wf_nprojtot = 0
      wfpro_projused = 0
      if( wfpro_nproj .eq. 0) wfpro_nproj = wfpro_nprojtot
      do iproj = 1, wfpro_nprojtot
        nr = maxloc( epslist, 1)
        if( (epslist( nr) .gt. eps_) .and. (wf_nprojtot .lt. wfpro_nproj)) then
          wf_nprojtot = wf_nprojtot + 1
          wfpro_projused( nr) = 1
        end if
        epslist( nr) = -1.d0
      end do
      deallocate( projolp, lapack_tau, lapack_work, epslist)

      if( allocated( wf_projst)) deallocate( wf_projst)
      allocate( wf_projst( 5, wf_nprojtot))

      i = 0
      do iproj = 1, wfpro_nprojtot
        if( wfpro_projused( iproj) .eq. 1) then
          i = i + 1
          wf_projst( :, i) = wfpro_projst( :, iproj)
        end if
      end do

      do igroup = 1, wf_ngroups
        if( allocated( wf_groups( igroup)%projused)) deallocate( wf_groups( igroup)%projused)
        allocate( wf_groups( igroup)%projused( wf_nprojtot))
        wf_groups( igroup)%projused = 1
        wf_groups( igroup)%nprojused = wf_nprojtot
      end do

      return
    end subroutine wfpro_removeldprojf

    subroutine wfpro_selectprojf
      integer :: igroup, nproj1, nproj2, iproj, ik, i, n
      real(8) :: loscore( wf_nprojtot), escore( wf_fst:wf_lst, wf_kset%nkpt), win(2), fac

      real(8), allocatable :: eval(:,:)

      call wannier_geteval( eval, i, n)
      fac = 4.d0
      
      do igroup = 1, wf_ngroups
        
        nproj1 = size( input%properties%wannier%grouparray( igroup)%group%projectorarray, 1)
        nproj2 = input%properties%wannier%grouparray( igroup)%group%nproj
        if( nproj2 .eq. 0) nproj2 = wf_nprojtot
        if( nproj1 .eq. 0) then
          wf_groups( igroup)%nprojused = min( wf_nprojtot, nproj2)
        else
          wf_groups( igroup)%nprojused = nproj1
        end if
        wf_groups( igroup)%projused = 0

        ! if special projectors are defined use them
        if( (wf_groups( igroup)%method .eq. "pro") .or. (wf_groups( igroup)%method .eq. "promax")) then
          if( wf_groups( igroup)%nprojused .ne. wf_groups( igroup)%nwf) then
            write( *, '(" Error (wfpro_selectprojf): The number of projectors must be equal to the number of Wannier functions for group ",I2,".")') igroup
            stop
          end if
          do iproj = 1, nproj1
            if( (input%properties%wannier%grouparray( igroup)%group%projectorarray( iproj)%projector%nr .lt. 1) .or. &
                (input%properties%wannier%grouparray( igroup)%group%projectorarray( iproj)%projector%nr .gt. wf_nprojtot)) then
              write( *, '(" Error (wfpro_selectprojf): ",I4," is not a valid index for projection local-orbitals in group ",I2,".")') &
                  input%properties%wannier%grouparray( igroup)%group%projectorarray( iproj)%projector%nr, igroup
              write(*, '(" Here is a list of local-orbitals that can be used for projection:")')
              call wfpro_showproj
              stop
            end if
            wf_groups( igroup)%projused( input%properties%wannier%grouparray( igroup)%group%projectorarray( iproj)%projector%nr) = 1
          end do
        else if( (wf_groups( igroup)%method .eq. "opf") &
            .or. (wf_groups( igroup)%method .eq. "opfmax") &
            .or. (wf_groups( igroup)%method .eq. "scdm") &
            .or. (wf_groups( igroup)%method .eq. "scdmmax") &
            .or. (wf_groups( igroup)%method .eq. "disentangle")) then
          if( wf_groups( igroup)%nprojused .lt. wf_groups( igroup)%nwf) then
            write( *, '(" Error (wfpro_selectprojf): The number of projectors must be greater or equal to the number of Wannier functions for group ",I2,".")') igroup
            stop
          end if
          do iproj = 1, nproj1
            if( (input%properties%wannier%grouparray( igroup)%group%projectorarray( iproj)%projector%nr .lt. 1) .or. &
                (input%properties%wannier%grouparray( igroup)%group%projectorarray( iproj)%projector%nr .gt. wf_nprojtot)) then
              write( *, '(" Error (wfpro_selectprojf): ",I4," is not a valid index for projection local-orbitals in group ",I2,".")') &
                  input%properties%wannier%grouparray( igroup)%group%projectorarray( iproj)%projector%nr, igroup
              write(*, '(" Here is a list of local-orbitals that can be used for projection:")')
              call wfpro_showproj
              stop
            end if
            wf_groups( igroup)%projused( input%properties%wannier%grouparray( igroup)%group%projectorarray( iproj)%projector%nr) = 1
          end do
        end if

        ! otherwise select best fitting ones
        if( wf_groups( igroup)%method .eq. "disentangle") then
          escore = 0.d0
          win = wf_groups( igroup)%win_i
          if( sum( abs( win)) .lt. 1.d-23) then
            win = wf_groups( igroup)%win_o
          end if
          do ik = 1, wf_kset%nkpt
            do i = 1, wf_groups( igroup)%win_ni( ik)
              n = wf_groups( igroup)%win_ii( i, ik)
              escore( n, ik) = escore( n, ik) + &
                1.d0/( 1.d0 + exp( -fac*( eval( n, ik) - win(1))/( win(2) - win(1)))) * &
                1.d0/( 1.d0 + exp(  fac*( eval( n, ik) - win(2))/( win(2) - win(1))))
            end do
            do i = 1, wf_groups( igroup)%win_no( ik)
              n = wf_groups( igroup)%win_io( i, ik)
              escore( n, ik) = escore( n, ik) + &
                1.d0/( 1.d0 + exp( -fac*( eval( n, ik) - win(1))/( win(2) - win(1)))) * &
                1.d0/( 1.d0 + exp(  fac*( eval( n, ik) - win(2))/( win(2) - win(1))))
            end do
          end do
        else
          escore = 1.d0
        end if
        if( nproj1 .eq. 0) then
          loscore = 0.d0
          do iproj = 1, wf_nprojtot
            do ik = 1, wf_kset%nkpt
              do i = wf_groups( igroup)%fst, wf_groups( igroup)%lst
                loscore( iproj) = loscore( iproj) + abs( wf_projection( i, iproj, ik))*escore( i, ik)
              end do
            end do
          end do
          n = 0
          do iproj = 1, wf_nprojtot
            i = maxloc( loscore, 1)
            if( n .lt. wf_groups( igroup)%nprojused) then
              n = n + 1
              wf_groups( igroup)%projused( i) = 1
            end if
            loscore( i) = -1.d0
          end do
        end if

      end do

      return
    end subroutine wfpro_selectprojf

end module mod_wannier_projection
