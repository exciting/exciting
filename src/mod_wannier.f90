module mod_wannier
  use mod_APW_LO
  use mod_atoms 
  use mod_kpoint
  use mod_constants
  use mod_muffin_tin
  use mod_Gkvector
  use mod_eigensystem
  use mod_spin
  use mod_eigenvalue_occupancy
  use mod_lattice
  use mod_symmetry
  use m_ematqk
  use m_plotmat

  use mod_Gvector

  implicit none

! variable
  integer :: wf_nprojtot, wf_nprojused, wf_bandstart, wf_nband
  
  integer, allocatable :: wf_projst(:,:), wf_projused(:)
  complex(8), allocatable :: wf_transform(:,:,:)

! methods
  contains
    !BOP
    ! !ROUTINE: wfinit
    ! !INTERFACE:
    !
    subroutine wfinit
      ! !USES:
      ! !DESCRIPTION:
      !   Reads local-orbitals from species files which are indicated to be used
      !   as projection functions for the generation of Wannier functions to the
      !   module {\tt mod\_wannier}.
      !
      ! !REVISION HISTORY:
      !   Created September 2016 (Sebastian)
      !EOP
      !BOC

      ! local variables
      integer :: is, ia, ilo, l, m
      ! read local-orbitals for projection
      allocate( wf_projst( nlotot, 5))
      wf_projst(:,:) = 0
      wf_nprojtot = 0
      do is = 1, nspecies
        do ia = 1, natoms( is)
          do ilo = 1, nlorb( is)
            if( lorbwfproj( ilo, is)) then
              l = lorbl( ilo, is)
              do m = -l, l
                wf_nprojtot = wf_nprojtot + 1
                wf_projst( wf_nprojtot, 1) = is
                wf_projst( wf_nprojtot, 2) = ia
                wf_projst( wf_nprojtot, 3) = ilo
                wf_projst( wf_nprojtot, 4) = l
                wf_projst( wf_nprojtot, 5) = m 
              end do
            end if
          end do
        end do
      end do
      if( wf_nprojtot .eq. 0) then
        write(*,*) 'ERROR (wf_init): No local-orbitals found for projection.'
        stop
      end if
      return
    end subroutine wfinit
    !EOC

    !BOP
    ! !ROUTINE: genwf
    ! !INTERFACE:
    !
    subroutine genwf( bandstart, nband, nproj, loproj)
      ! !USES:
      ! !INPUT/OUTPUT PARAMETERS:
      !   bandstart : n1, lowest band index of the band range used for generation of
      !               Wannier functions (in,integer)
      !   nband     : N, number of bands used for generation (in,integer)
      !   nproj     : N, number of projection orbitals used for generation (in,integer)
      !   loproj    : indices of local-orbitals used as projection functions
      !               (in, integer(nproj))
      ! !DESCRIPTION:
      !   Generates unitary transformation matrices $U_{\vec{k}}$ via the
      !   projection method from which a set of smooth Bloch-like states 
      !   $$ \left|\Psi_{m,\vec{k}}^W\right\rangle = \sum\limits_{n=n_1}^{n_2}
      !   \left(U_{\vec{k}}\right)_{nm} \left|\Psi_{n,\vec{k}}^H\right\rangle $$
      !   out of the $N=n_2-n_1+1$ Hamiltonian eigenstates $\Psi_{n,\vec{k}}^H$
      !   in the band range $[n_1,n_2]$ can be constructed and used for the
      !   generation of localized Wannier functions
      !   $$ \left| \vec{R}n\right\rangle = \frac{1}{N_k} \sum\limits_{\vec{k}}
      !   e^{-i\vec{k}\cdot\vec{R}} \left| \Psi_{n,\vec{k}}^W\right\rangle $$
      !   The matrices are stored in the module variable {\tt wf\_transform
      !   (complex(nband,nproj,nkptnr))}.
      !
      ! !REVISION HISTORY:
      !   Created September 2016 (Sebastian)
      !EOP
      !BOC

      integer, intent( in) :: bandstart, nband, nproj, loproj( nproj) 

      ! local variables
      integer :: iproj, is, js, ja, ias, nr, l, lm, io, ilo, ir, ig, ik, i, m
      integer :: isym, iknr, ngknr 
      integer, allocatable :: igkignr(:)
      real(8), allocatable :: vgklnr(:,:,:), vgkcnr(:,:,:), gkcnr(:), tpgkcnr(:,:)
      complex(8), allocatable :: sfacgknr(:,:)
      complex(8), allocatable :: evecfv(:,:,:), apwalm(:,:,:,:,:), olpm(:,:), projm(:,:), projmhyb(:,:), auxmat(:,:), auxmat2(:,:), evec(:,:)
      real(8), allocatable :: rolpi(:,:), uf(:), gf(:), cf(:,:), eval(:)

      if( .not. allocated( wf_projst)) call wfinit

      if( allocated( wf_projused)) deallocate( wf_projused)
      allocate( wf_projused( nproj))
      wf_bandstart = bandstart
      wf_nband = nband
      wf_projused = 0
      wf_nprojused = 0
      do iproj = 1, nproj
        if( (loproj( iproj) .lt. 1) .or. (loproj( iproj) .gt. wf_nprojtot)) then
          write(*,'(a,i2,a)') ' ERROR (genwf): ', loproj( iproj), ' is not a valid index for projection local-orbitals.'
          write(*,*) 'Here is a list of local-orbitals that can be used for projection:'
          call wfshowproj
          return
        else
          wf_projused( iproj) = loproj( iproj)
          wf_nprojused = wf_nprojused + 1
        end if
      end do
         
      ! radial overlap integrals
      call readstate
      call linengy
      call genapwfr     ! APW radial functions
      call genlofr      ! LO radial functions

      allocate( rolpi( wf_nprojused, apwordmax+maxlorb))
      allocate( uf( nrmtmax), gf( nrmtmax), cf( 3, nrmtmax))
      rolpi(:,:) = zzero
      do iproj = 1, wf_nprojused
        uf(:) = 0d0
        is = wf_projst( wf_projused( iproj), 1)
        nr = nrmt( is)
        ias = idxas( wf_projst( wf_projused( iproj), 2), is)
        lm = idxlm( wf_projst( wf_projused( iproj), 4), wf_projst( wf_projused( iproj), 5))
        do io = 1, apword( wf_projst( wf_projused( iproj), 4), is)
          do ir = 1, nr
            uf( ir) = apwfr( ir, 1, io, wf_projst( wf_projused( iproj), 4), ias)*lofr( ir, 1, wf_projst( wf_projused( iproj), 3), ias)*spr( ir, is)**2
          end do
          call fderiv( -1, nr, spr(:, is), uf, gf, cf)
          rolpi( iproj, io) = gf( nr)
        end do
        do ilo = 1, nlorb( is)
          l = lorbl( ilo, is)
          if( l .eq. wf_projst( wf_projused( iproj), 4)) then
            do ir = 1, nr
              uf( ir) = lofr( ir, 1, ilo, ias)*lofr( ir, 1, wf_projst( wf_projused( iproj), 3), ias)*spr( ir, is)**2
            end do
            call fderiv( -1, nr, spr( :, is), uf, gf, cf)
            rolpi( iproj, apwordmax+ilo) = gf( nr)
          end if
        end do
      end do

      allocate( evecfv( nmatmax, nstfv, nspnfv))
      allocate( apwalm( ngkmax, apwordmax, lmmaxapw, natmtot, nspnfv))
      allocate( projm( wf_nband, wf_nprojused), &
                olpm( wf_nprojused, wf_nprojused), &
                projmhyb( wf_nband, wf_nprojused))
      allocate( eval( wf_nprojused), evec( wf_nprojused, wf_nprojused), auxmat( ngkmax+nlotot, wf_nprojused), auxmat2( wf_nprojused, wf_nprojused))
      if( allocated( wf_transform)) deallocate( wf_transform)
      allocate( wf_transform( wf_nband, wf_nprojused, nkptnr))
      allocate( igkignr( ngkmax))
      allocate( vgklnr( 3, ngkmax, nspnfv), vgkcnr( 3, ngkmax, nspnfv), gkcnr( ngkmax), tpgkcnr( 2, ngkmax))
      allocate( sfacgknr( ngkmax, natmtot))
#ifdef USEOMP
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, evecfv, vgklnr, isym, ik, ngknr, igkignr, vgkcnr, gkcnr, tpgkcnr, sfacgknr, apwalm, iproj, is, ias, lm, ig, auxmat, io, js, ja, ilo, l, m, projm, olpm, auxmat2, eval, evec)
!!$OMP DO
#endif
      write(*,*) wf_nband, nkptnr, wf_nprojused
      do iknr = 1, nkptnr
        evecfv = zzero
        vgklnr = zzero
        call findkpt( vklnr( :, iknr), isym, ik)
        ! find G+k-vectors for non-reduced k-point
        Call gengpvec( vklnr( :, iknr), vkcnr( :, iknr), ngknr, igkignr, vgklnr(:,:,1), vgkcnr(:,:,1), gkcnr, tpgkcnr)
        ! find structure factors
        Call gensfacgp( ngknr, vgkcnr, ngkmax, sfacgknr)
        ! get basis function coefficients and matching coefficients
        call getevecfv( vklnr(:, iknr), vgklnr, evecfv)
        call match( ngknr, gkcnr, tpgkcnr, sfacgknr, apwalm(:, :, :, :, 1))
        
        ! projection matrix elements and overlap matrix elements
        do iproj = 1, wf_nprojused
          is = wf_projst( wf_projused( iproj), 1)
          ias = idxas( wf_projst( wf_projused( iproj), 2), is)
          lm = idxlm( wf_projst( wf_projused( iproj), 4), wf_projst( wf_projused( iproj), 5))
          do ig = 1, ngknr
            auxmat( ig, iproj) = zzero
            do io = 1, apword( wf_projst( wf_projused( iproj), 4), is)
              auxmat( ig, iproj) = auxmat( ig, iproj) + conjg( apwalm( ig, io, lm, ias, 1))*rolpi( iproj, io)
            end do
          end do
          do js = 1, nspecies
            do ja = 1, natoms( js)
              do ilo = 1, nlorb( js)
                l = lorbl( ilo, js)
                do m = -l, l
                  ig = ngknr + idxlo( idxlm( l, m), ilo, idxas( ja, js))
                  auxmat( ig, iproj) = zzero
                  if( (idxas( ja, js) .eq. ias) .and. (idxlm( l, m) .eq. lm)) then
                    auxmat( ig, iproj) = cmplx( rolpi( iproj, apwordmax+ilo), 0, 8)
                  end if
                end do
              end do
            end do
          end do
        end do
        call ZGEMM( 'C', 'N', wf_nband, wf_nprojused, ngknr+nlotot, zone, &
             evecfv( 1:(ngknr+nlotot), wf_bandstart:(wf_bandstart+wf_nband-1), 1), ngknr+nlotot, &
             auxmat( 1:(ngknr+nlotot), :), ngknr+nlotot, zzero, &
             projm, wf_nband)

        ! XYZ.amn
        !do ig = 1, wf_nprojused
        !  do i = 1, wf_nband
        !    write( *, '(I,I,I,F23.16,F23.16)') i, ig, iknr, projm( i, ig)
        !  end do
        !end do
        !call plotmat( projm)

        ! hybridization
        !projmhyb( :, 1) = 1*projm( :, 1) + 1*projm( :, 2) + 1*projm( :, 3) + 1*projm( :, 4)
        !projmhyb( :, 2) = 1*projm( :, 1) + 1*projm( :, 2) - 1*projm( :, 3) - 1*projm( :, 4)
        !projmhyb( :, 3) = 1*projm( :, 1) - 1*projm( :, 2) - 1*projm( :, 3) + 1*projm( :, 4)
        !projmhyb( :, 4) = 1*projm( :, 1) - 1*projm( :, 2) + 1*projm( :, 3) - 1*projm( :, 4)

        !projm = projmhyb

        olpm(:,:) = zzero
        call ZGEMM( 'C', 'N', wf_nprojused, wf_nprojused, wf_nband, zone, &
             projm, wf_nband, &
             projm, wf_nband, zzero, &
             olpm, wf_nprojused)

        ! transformation matrices
        call diaghermat( wf_nprojused, olpm, eval, evec)
        auxmat2(:,:) = zzero
        i = 1
        do io = 1, wf_nprojused
          !write(*,*) io, eval( io)
          if( abs( eval( io)) .le. 1.d-10) then
            i = i + 1
            !write(*,*) "ERROR (genwf): 0 eigenvalue ", iknr, io
          else
            auxmat2( io, io) = zone/sqrt( cmplx( eval( io), 0, 8)) 
          end if
        end do

        call ZGEMM( 'N', 'N', wf_nprojused, wf_nprojused-i+1, wf_nprojused-i+1, zone, &
             evec( :, i:wf_nprojused), wf_nprojused, &
             auxmat2( i:wf_nprojused, i:wf_nprojused), wf_nprojused-i+1, zzero, &
             olpm( :, i:wf_nprojused), wf_nprojused)
        call ZGEMM( 'N', 'C', wf_nprojused, wf_nprojused, wf_nprojused-i+1, zone, &
             olpm( :, i:wf_nprojused), wf_nprojused, &
             evec( :, i:wf_nprojused), wf_nprojused, zzero, &
             auxmat2, wf_nprojused)
        call ZGEMM( 'N', 'N', wf_nband, wf_nprojused, wf_nprojused, zone, &
             projm, wf_nband, &
             auxmat2, wf_nprojused, zzero, &
             wf_transform( :, :, iknr), wf_nband)
      end do 
#ifdef USEOMP
!!$OMP END DO
!!$OMP END PARALLEL
#endif
      deallocate( rolpi, projm, olpm, auxmat, auxmat2, eval, evec, apwalm, evecfv, uf, gf, cf)
      return
    end subroutine genwf
    !EOC
    
    !BOP
    ! !ROUTINE: genmlwf
    ! !INTERFACE:
    !
    subroutine genmlwf( bandstart, nband, nproj, loproj)
      ! !USES:
      ! !INPUT/OUTPUT PARAMETERS:
      !   bandstart : n1, lowest band index of the band range used for generation of
      !               Wannier functions (in,integer)
      !   nband     : N, number of bands used for generation (in,integer)
      !   loproj    : indices of local-orbitals used as projection functions
      !               (in, integer(nband))
      ! !DESCRIPTION:
      !   Does the same thing as {\tt genwf} does but the transformation
      !   matrices are used for the generation of maximally localized Wannier
      !   functions. The matrices are computed in a self consistent loop. 
      !
      ! !REVISION HISTORY:
      !   Created September 2016 (Sebastian)
      !EOP
      !BOC

      integer, intent( in) :: bandstart, nband, nproj, loproj( nproj) 
      
      ! local variables
      integer :: ngridk(3), nvec(3), ix, iy, iz, i, ib, ik, iknr, ngknr
      character :: read_file
      logical :: file_exists, reading_succes
      real(8) :: bwgtsm, mixing, maxdiff, vec1(3), vec2(3)

      ! allocatable arrays
      integer, allocatable :: idxn(:,:), nn(:), diffk(:)
      integer, allocatable :: igkignr(:)
      real(8), allocatable :: vgklnr(:,:,:), vgkcnr(:,:,:), gkcnr(:), tpgkcnr(:,:)
      real(8), allocatable :: ndist(:), nvl(:,:,:), nvc(:,:,:), bwgt(:), ravg(:,:), diffval(:,:)
      complex(8), allocatable :: mlwf_emat(:,:), auxmat(:,:), mlwf_m0(:,:,:,:), mlwf_m(:,:,:,:), eval(:), evec(:,:), mlwf_r(:,:), mlwf_t(:,:), mlwf_dw(:,:), mlwf_transform(:,:,:)
      complex(8), allocatable :: evecfv1(:,:,:), evecfv2(:,:,:)

      ! write next neighbours of each k-point to array
      allocate( ndist( nkptnr), nn( nkptnr), &
                nvl( 3, nkptnr, nkptnr), &
                nvc( 3, nkptnr, nkptnr)) 
      call wfneighbors( ndist, nn, nvl, nvc) 
      allocate( idxn( nn(1), nkptnr))
      ngridk = input%groundstate%ngridk
      do iz = 0, ngridk( 3)-1
        do iy = 0, ngridk( 2)-1
          do ix = 0, ngridk( 1)-1
            iknr = modulo( iz, ngridk( 3))*ngridk( 2)*ngridk( 1) + &
                   modulo( iy, ngridk( 2))*ngridk( 1) + modulo( ix, ngridk(1))+1
            do i = 1, nn(1)
              idxn( i, iknr) = modulo( iz + nint( nvl( 3, 1, i)*ngridk(3)), ngridk( 3))*ngridk( 2)*ngridk( 1) + &
                               modulo( iy + nint( nvl( 2, 1, i)*ngridk(2)), ngridk( 2))*ngridk( 1) + &
                               modulo( ix + nint( nvl( 1, 1, i)*ngridk(1)), ngridk(1)) + 1
            end do
          end do
        end do
      end do
      
      read_file = 'n'
      ! check wether WANNIER_MLWF.OUT exists and ask for use
      inquire( file='WANNIER_MLWF.OUT', exist=file_exists)
      if( file_exists) then
          write(*,*) 'There is already a file called WANNIER_MLWF.OUT.'
          write(*,*) 'Should transformation matrices for MLWF be read from this? (y/n)'
          read( *, '(1a)') read_file
      end if

      if( read_file .eq. 'y') then
        ! read transformation matrices from file
        call wfinit
        call wfreadfile( 'WANNIER_MLWF.OUT', reading_succes)
        if( .not. reading_succes) then
          write(*,*) 'WANNIER_MLWF.OUT contains invalid content. Reading aborted. Start recalculation.'
          read_file = 'n'
        end if
      end if
      if( read_file .eq. 'n') then
        ! generate transformation matrices from projection as initial input
        call genwf( bandstart, nband, nproj, loproj)
        ! calculating inner products <u(m,k)|u(n,k+b)>
        allocate( mlwf_emat( wf_nband, wf_nband))
        !allocate( mlwf_emat( nstfv, nstfv))
        allocate( auxmat( wf_nband, wf_nprojused))
        allocate( mlwf_m0( wf_nprojused, wf_nprojused, nn(1), nkptnr))
        write(*,*) 'Computing inner products...'
        ! call before parallel loop in order to initialize save variables in emat_wannier

        allocate( evecfv1( nmatmax, nstfv, nspnfv), evecfv2( nmatmax, nstfv, nspnfv))
        allocate( igkignr( ngkmax))
        allocate( vgklnr( 3, ngkmax, nspnfv), vgkcnr( 3, ngkmax, nspnfv), gkcnr( ngkmax), tpgkcnr( 2, ngkmax))

        do ib = 1, nn(1) !!ATTENTION!!nn
          call emat_init( nvl( :, 1, ib), (/0d0, 0d0, 0d0/), input%groundstate%lmaxapw, 8)
#ifdef USEOMP
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, vkcnr, ngknr, igkignr, vgklnr, vgkcnr, gkcnr, tpgkcnr, evecfv1, evecfv2, auxmat, mlwf_emat)
!!$OMP DO
#endif
          do iknr = 1, nkptnr
            
            ! find G+k-vectors and eigenvectors for non-reduced k-point k
            Call gengpvec( vklnr( :, iknr), vkcnr( :, iknr), ngknr, igkignr, vgklnr(:,:,1), vgkcnr(:,:,1), gkcnr, tpgkcnr)
            call getevecfv( vklnr( :, iknr), vgklnr, evecfv1(:,:,1))
            ! find G+k-vectors and eigenvectors for non-reduced k-point k+b
            Call gengpvec( vklnr( :, idxn( ib, iknr)), vkcnr( :, idxn( ib, iknr)), ngknr, igkignr, vgklnr(:,:,1), vgkcnr(:,:,1), gkcnr, tpgkcnr)
            call getevecfv( vklnr( :, idxn( ib, iknr)), vgklnr, evecfv2(:,:,1))

            call emat_genemat( iknr, idxn( ib, iknr), wf_bandstart, wf_nband, wf_bandstart, wf_nband, evecfv1(:,:,1), evecfv2(:,:,1), mlwf_emat)

            call ZGEMM( 'N', 'N', wf_nband, wf_nprojused, wf_nband, zone, &
                 mlwf_emat, wf_nband, &
                 wf_transform( :, :, idxn( ib, iknr)), wf_nband, zzero, &
                 auxmat, wf_nband)
            call ZGEMM( 'C', 'N', wf_nprojused, wf_nprojused, wf_nband, zone, &
                 wf_transform( :, :, iknr), wf_nband, &
                 auxmat, wf_nband, zzero, &
                 mlwf_m0( :, :, ib, iknr), wf_nprojused)
            ix = ix + 1
            ! print out progress
            write( 6, '(a,10x,i3,a)', advance='no') achar( 13), nint( 100.0d0*ix/nkptnr/nn(1)), '%'
            flush( 6)
            !write( *, '(5I6)') iknr, idxn( ib, iknr), nint( vklnr( :, iknr) + nvl( :, 1, ib) - vklnr( :, idxn( ib, iknr)))
            !do iy = 1, wf_nband
            !  do iz = 1, wf_nband
            !    write( *, '(F23.16,F23.16)') mlwf_emat( iz, iy)
            !  end do
            !end do
          end do
#ifdef USEOMP
!!$OMP END DO
!!$OMP END PARALLEL
#endif
        end do
        write(6,*)

        !stop

        ! calculating weights of b-vectors
        allocate( bwgt( nn(1)))
        bwgtsm = 0d0
        do ib = 1, nn(1)
          bwgt( ib) = 3d0/(nn(1)*norm2( nvc( :, 1, ib))**2)
          bwgtsm = bwgtsm + bwgt( ib)
        end do
        
        ! mixing parameter for self consistent minimization
        mixing = dble( 0.5)
    
        ! initialize transformation matrices
        allocate( mlwf_transform( wf_nprojused, wf_nprojused, nkptnr))
        mlwf_transform = zzero
        do ix = 1, wf_nprojused
          mlwf_transform( ix, ix, :) = zone
        end do

        ! start minimization loop
        allocate( evec( wf_nprojused, wf_nprojused), eval( wf_nprojused))
        allocate( mlwf_r( wf_nprojused, wf_nprojused), mlwf_t( wf_nprojused, wf_nprojused), mlwf_dw( wf_nprojused, wf_nprojused))
        allocate( mlwf_m( wf_nprojused, wf_nprojused, nn(1), nkptnr))
        allocate( ravg( 3, wf_nprojused))
        allocate( diffval( 0:5000, nkptnr), diffk( 0:5000))
        mlwf_m = mlwf_m0
        ! maximum norm of difference between old and new transformation matrices
        diffval = 1.0d0
        iz = 0
        write(*,*) 'Waiting for MLWF to converge...'
        do while( (iz .lt. 5000) .and. (maxval( diffval( iz, :)) .ge. dble( 1.0E-8)))
          iz = iz + 1
          diffval( iz, :) = 0.d0
          ravg = 0.d0
          do ix = 1, wf_nprojused
            do iknr = 1, nkptnr
              do ib = 1, nn(1)
                ravg( :, ix) = ravg( :, ix) - bwgt( ib)/nkptnr*real( aimag( log( mlwf_m( ix, ix, ib, iknr))))*nvc( :, 1, ib)
              end do
            end do
          end do
      
#ifdef USEOMP
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, ib, mlwf_r, mlwf_t, mlwf_dw, eval, evec, ix, auxmat, maxdiff)
!!$OMP DO
#endif
          do iknr = 1, nkptnr
            mlwf_dw(:,:) = zzero
            do ib = 1, nn(1)
              ! calculating R and T
              do ix = 1, wf_nprojused
                mlwf_r( :, ix) = mlwf_m( :, ix, ib, iknr)*conjg( mlwf_m( ix, ix, ib, iknr))
                mlwf_t( :, ix) = mlwf_m( :, ix, ib, iknr)/mlwf_m( ix, ix, ib, iknr)*(real( aimag( log( mlwf_m( ix, ix, ib, iknr)))) + dot_product( nvc( :, 1, ib), ravg( :, ix)))
              end do
        
              ! calculating dW
              mlwf_r = mlwf_r - conjg( transpose( mlwf_r))
              mlwf_t = mlwf_t + conjg( transpose( mlwf_t))
              mlwf_dw(:,:) = mlwf_dw(:,:) + bwgt( ib)*( 0.5*mlwf_r(:,:) + 0.5*zi*mlwf_t(:,:))
            end do
            mlwf_dw(:,:) = mlwf_dw(:,:)*mixing/bwgtsm

            ! updating transformation matrices
            call diaggenmat( wf_nprojused, mlwf_dw, eval, evec)
            auxmat(:,:) = zzero
            do ix = 1, wf_nprojused
              auxmat( ix, ix) = exp( eval( ix)) 
            end do
            call ZGEMM( 'N', 'N', wf_nprojused, wf_nprojused, wf_nprojused, zone, &
                 evec, wf_nprojused, &
                 auxmat, wf_nprojused, zzero, &
                 mlwf_dw, wf_nprojused)
            call ZGEMM( 'N', 'C', wf_nprojused, wf_nprojused, wf_nprojused, zone, &
                 mlwf_dw, wf_nprojused, &
                 evec, wf_nprojused, zzero, &
                 auxmat, wf_nprojused)
            call ZGEMM( 'N', 'N', wf_nprojused, wf_nprojused, wf_nprojused, zone, &
                 mlwf_transform( :, :, iknr), wf_nprojused, &
                 auxmat, wf_nprojused, zzero, &
                 mlwf_dw, wf_nprojused)
            maxdiff = norm2( abs( reshape( mlwf_transform( :, :, iknr) - mlwf_dw, (/wf_nprojused**2/))))!/norm2( abs( reshape( mlwf_transform( :, :, ik), (/wf_nprojused**2/))))
            if( maxdiff .ge. maxval( diffval( iz, :))) diffk( iz) = iknr
            diffval( iz, iknr) = maxdiff

            mlwf_transform( :, :, iknr) = mlwf_dw(:,:)
          end do
#ifdef USEOMP
!!$OMP END DO
!!$OMP END PARALLEL
#endif
    
          ! updating M
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, ib, auxmat)
!$OMP DO
#endif
          do iknr = 1, nkptnr
            do ib = 1, nn(1)
              call ZGEMM( 'N', 'N', wf_nprojused, wf_nprojused, wf_nprojused, zone, &
                   mlwf_m0( :, :, ib, iknr), wf_nprojused, &
                   mlwf_transform( :, :, idxn( ib, iknr)), wf_nprojused, zzero, &
                   auxmat, wf_nprojused)
              call ZGEMM( 'C', 'N', wf_nprojused, wf_nprojused, wf_nprojused, zone, &
                   mlwf_transform( :, :, iknr), wf_nprojused, &
                   auxmat, wf_nprojused, zzero, &
                   mlwf_m( :, :, ib, iknr), wf_nprojused)
            end do
          end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif

          ! print out progress
          write( 6, '(a,10x,a,i4,a,E9.2)', advance='no') achar( 13), '[', iz, ']dU =', maxdiff
          flush( 6)
        end do
        write(*,*)
        if( iz .ge. 10000) then
          write(*,*) 'ERROR: Not converged after 10000 cycles.'
        else
          write(*,'(a35,i4,a8)') ' SUCCES: Convergence reached after ', iz, ' cycles.' 
        end if
        
        !do iz = 1, 10000
        !  write(*,'(I6)', advance='no') iz
        !  do ik = 1, nkpt
        !    write(*,'(3x,F16.12)', advance='no') diffval( iz, ik)
        !  end do
        !  write(*,*)
        !end do

        ! generating final transformation matrices for Hamiltonian eigenstates
        do iknr = 1, nkptnr
          call ZGEMM( 'N', 'N', wf_nprojused, wf_nprojused, wf_nprojused, zone, &
               wf_transform( :, :, iknr), wf_nprojused, &
               mlwf_transform( :, :, iknr), wf_nprojused, zzero, &
               auxmat, wf_nprojused)
          wf_transform( :, :, iknr) = auxmat(:,:)
        end do

        ! write resulting transformation matrices to file
        call wfwritefile( 'WANNIER_MLWF.OUT')
  
        deallocate( auxmat, eval, evec, mlwf_m0, mlwf_m, mlwf_r, mlwf_t, mlwf_dw, mlwf_transform)
      end if


    end subroutine genmlwf
    !EOC
    
    ! print out projection local-orbitals
    subroutine wfshowproj
      ! local variables
      integer :: iproj

      write(*,*)
      write(*,*) 'Local-orbitals for generating Wannier functions via projection.'
      write(*,*) '     nr   species   atom   l    m    used'
      do iproj = 1, wf_nprojtot
        write( *, '(5x,i2,3x)', advance='no') iproj
        write( *, '(2x,a5,3x)', advance='no') spsymb( wf_projst( iproj, 1))
        write( *, '(i4,3x)', advance='no') wf_projst( iproj, 2)
        write( *, '(i2,3x)', advance='no') wf_projst( iproj, 4)
        write( *, '(i2,3x)', advance='no') wf_projst( iproj, 5)
        if( any( wf_projused .eq. iproj)) write( *, '(1x,a)', advance='no') '*'
        write(*,*)
      end do
      write(*,*)
      return
    end subroutine wfshowproj
    
    ! reads transformation matrices from file
    subroutine wfreadfile( filename, succes)
      character(*), intent( in) :: filename
      logical, intent( out) :: succes

      ! local variables
      integer :: ik, ix, iy

      succes = .true.
      open( 50, file=trim( filename), action='READ', form='FORMATTED')
      read( 50, '(i2)') wf_bandstart
      read( 50, '(i2)') wf_nprojused
      if( wf_nprojused .gt. wf_nprojtot) then
        write(*,*) 'ERROR: Content in file does not fit to content in species files. Check local-orbitals.'
        succes = .false.
        return
      end if
      if( allocated( wf_projused)) deallocate( wf_projused)
      allocate( wf_projused( wf_nprojused))
      if( allocated( wf_transform)) deallocate( wf_transform)
      allocate( wf_transform( wf_nprojused, wf_nprojused, nkpt))
      do ix = 1, wf_nprojused
        read( 50, '(i2,2x)', advance='no') wf_projused( ix)
      end do
      read( 50, *)
      outerloop: do ik = 1, nkpt
        read( 50, '(i)') ix
        if( .not.( ix .eq. ik)) then
          succes = .false.
          exit outerloop
        end if
        do iy = 1, wf_nprojused
          do ix = 1, wf_nprojused
            read( 50, '(E23.16,3x,E23.16)') wf_transform( ix, iy, ik)
          end do
        end do
      end do outerloop
      close( 50)
      if( succes) write(*,*) 'Transformation matrices succesfully read.'
      return
    end subroutine wfreadfile
    
    ! writes transformation matrices to file
    subroutine wfwritefile( filename)
      character(*), intent( in) :: filename

      ! local variables
      integer :: ik, ix, iy

      open( 50, file=trim( filename), action='WRITE', form='FORMATTED')
      write( 50, '(i2,3x,a)') wf_bandstart, 'starting band'
      write( 50, '(i2,3x,a)') wf_nprojused, 'number of bands'
      do ix = 1, wf_nprojused
        write( 50, '(i2,2x)', advance='no') wf_projused( ix)
      end do
      write( 50, *) 'projectors used'
      do ik = 1, nkpt
        write( 50, '(i)') ik
        do iy = 1, wf_nprojused
          do ix = 1, wf_nprojused
            write( 50, '(E23.16,3x,E23.16)') wf_transform( ix, iy, ik)
          end do
        end do
      end do
      close( 50)
      write( *, '(a,a)') ' Transformation matrices written to file ', trim( filename)
      return
    end subroutine wfwritefile  

    subroutine angchar( iknr, lmax, nmin, nmax)
      integer, intent( in) :: iknr, lmax, nmin, nmax

      integer :: n, l, m, lm, io, ilo, ia, is, ias, ngknr
      integer, allocatable :: igkignr(:)
      real(8), allocatable :: vgklnr(:,:,:), vgkcnr(:,:,:), gkcnr(:), tpgkcnr(:,:)
      complex(8), allocatable :: sfacgknr(:,:), evecfv1(:,:,:)
      complex(8), allocatable :: evecfv(:,:,:), apwalm(:,:,:,:,:)
      complex(8), allocatable :: output(:,:,:)

      allocate( evecfv( nmatmax, nstfv, nspnfv))
      allocate( apwalm( ngkmax, apwordmax, lmmaxapw, natmtot, nspnfv))
      allocate( igkignr( ngkmax))
      allocate( vgklnr( 3, ngkmax, nspnfv), vgkcnr( 3, ngkmax, nspnfv), gkcnr( ngkmax), tpgkcnr( 2, ngkmax))
      allocate( sfacgknr( ngkmax, natmtot))
      allocate( output( natmtot, nstfv, lmmaxapw))

      ! find G+k-vectors for non-reduced k-point
      Call gengpvec( vklnr( :, iknr), vkcnr( :, iknr), ngknr, igkignr, vgklnr(:,:,1), vgkcnr(:,:,1), gkcnr, tpgkcnr)
      ! find structure factors
      Call gensfacgp( ngknr, vgkcnr, ngkmax, sfacgknr)
      ! get basis function coefficients and matching coefficients
      call getevecfv( vklnr(:, iknr), vgklnr, evecfv)
      call match( ngknr, gkcnr, tpgkcnr, sfacgknr, apwalm(:, :, :, :, 1))
      output(:,:,:) = zzero
      do is = 1, nspecies
        do ia = 1, natoms( is)
          ias = idxas( ia, is)
          do io = 1, apword( l, is)
            call ZGEMM( 'T', 'N', nstfv, lmmaxapw, ngknr, zone, &
                 evecfv( 1:ngknr, :, 1), ngknr, &
                 apwalm( 1:ngknr, io, 1:lmmaxapw, ias, 1), ngknr, zone, &
                 output( ias, :, :), nstfv)
          end do
          !do ilo = 1, nlorb( is)
          !  l = lorbl( ilo, is)
          !  do m = -l, l
          !    lm = idxlm( l, m)
          !    output( ias, :, lm) = output( ias, :, lm) + evecfv( ngknr+idxlo( lm, ilo, ias), 1:nmax, 1) 
          !  end do
          !end do
        end do
      end do

      write( *, '("k-point: ",I2,3F13.6)') iknr, vklnr( :, iknr)
      write( *, '("band    l     m       atoms")')
      do n = nmin, nmax
        do l = 0, lmax
          do m = -l, l
            lm = idxlm( l, m)
            write( *, '(3(I3,3x))', advance='no') n, l, m
            do is = 1, nspecies
              do ia = 1, natoms( is)
                ias = idxas( ia, is)
                write( *, '(F8.1,SP,F8.1,3x)', advance='no') 100*abs( output( ias, n, lm))/dot_product( sqrt( abs( output( ias, n, :))), sqrt( abs( output( ias, n, :)))), atan2( real( aimag( output( ias, n, lm))), real( real( output( ias, n, lm))))/twopi*360
              end do
            end do
            write(*,*)
          end do
        end do
        write(*,*) 
      end do
      return
    end subroutine angchar

    subroutine wfneighbors( d, m, nvl, nvc)
      integer, intent( out) :: m( nkptnr)
      real(8), intent( out) :: d( nkptnr), nvl( 3, nkptnr, nkptnr), nvc( 3, nkptnr, nkptnr)
      integer :: ix, iy, iz, ngridk(3), i, j, n 
      integer :: mt( nkptnr)
      real(8) :: vl(3), vc(3), dist
      real(8) :: nvlt( 3, nkptnr, nkptnr), nvct( 3, nkptnr, nkptnr), dt( nkptnr)
      
      dt = 0.d0
      mt = 0
      nvlt = 0.d0
      nvct = 0.d0
      ngridk = input%groundstate%ngridk
      i = 0
      do iz = 0, ngridk(3)-1
        do iy = 0, ngridk(2)-1
          do ix = 0, ngridk(1)-1
            vl = (/dble( ix)/ngridk(1), dble( iy)/ngridk(2), dble( iz)/ngridk(3)/)
            call r3mv( bvec, vl, vc)
            dist = norm2( vc)
            if( minval( abs( dt(:) - dist)) .gt. input%structure%epslat) then
              i = i + 1
              dt( i) = dist
            end if
          end do
        end do
      end do
      n = i
      do iz = ngridk(3)-1, -ngridk(3)+1, -1
        do iy = ngridk(2)-1, -ngridk(2)+1, -1
          do ix = ngridk(1)-1, -ngridk(1)+1, -1
            vl = (/dble( ix)/ngridk(1), dble( iy)/ngridk(2), dble( iz)/ngridk(3)/)
            call r3mv( bvec, vl, vc)
            dist = norm2( vc)
            j = minloc( abs( dt(:) - dist), 1)
            if( abs( dt( j) - dist) .lt. input%structure%epslat) then
              mt( j) = mt( j) + 1
              nvlt( :, j, mt( j)) = vl
              nvct( :, j, mt( j)) = vc
            end if
          end do
        end do
      end do
      m = 0
      d = 0.d0
      nvl = 0.d0
      nvc = 0.d0
      dist = minval( dt)
      do i = 1, n
        j = minloc( abs( dt( 1:n) - dist), 1)
        d( i) = dt( j)
        m( i) = mt( j)
        nvl( :, i, :) = nvlt( :, j, :)
        nvc( :, i, :) = nvct( :, j, :)
        dt( j) = 1.d6 
        dist = minval( dt)
      end do
      return
    end subroutine wfneighbors
end module mod_wannier
