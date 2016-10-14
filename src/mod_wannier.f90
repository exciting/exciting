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
  use m_ematqk

  use mod_Gvector

  implicit none

! variable
  integer :: wf_nprojtot, wf_nprojused, wf_bandstart
  
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
      !   Created September 2016 (ST)
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
    subroutine genwf( bandstart, nband, loproj)
      ! !USES:
      ! !INPUT/OUTPUT PARAMETERS:
      !   bandstart : n1, lowest band index of the band range used for generation of
      !               Wannier functions (in,integer)
      !   nband     : N, number of bands used for generation (in,integer)
      !   loproj    : indices of local-orbitals used as projection functions
      !               (in, integer(nband))
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
      !   (complex(nband,nband,nkpt))}.
      !
      ! !REVISION HISTORY:
      !   Created September 2016 (ST)
      !EOP
      !BOC

      integer, intent( in) :: bandstart, nband, loproj( nband) 

      ! local variables
      integer :: iband, iproj, is, ias, nr, l, lm, io, ir, ig, ik
      complex(8) :: auxc
      complex(8), allocatable :: evecfv(:,:,:), apwalm(:,:,:,:,:), olpm(:,:), projm(:,:), auxmat(:,:), evec(:,:)
      real(8), allocatable :: rolpi(:,:), uf(:), gf(:), cf(:,:), eval(:)

      if( .not. allocated( wf_projst)) call wfinit

      if( allocated( wf_projused)) deallocate( wf_projused)
      allocate( wf_projused( nband))
      wf_bandstart = bandstart
      wf_projused = 0
      wf_nprojused = 0
      do iband = 1, nband
        if( (loproj( iband) .lt. 1) .or. (loproj( iband) .gt. wf_nprojtot)) then
          write(*,'(a,i2,a)') ' ERROR (genwf): ', loproj( iband), ' is not a valid index for projection local-orbitals.'
          write(*,*) 'Here is a list of local-orbitals that can be used for projection:'
          call wfshowproj
          return
        else
          wf_projused( iband) = loproj( iband)
          wf_nprojused = wf_nprojused + 1
        end if
      end do
         
      ! radial overlap integrals
      call readstate
      call linengy
      call genapwfr
      call genlofr
      call olprad

      allocate( rolpi( wf_nprojused, apwordmax))
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
          call fderiv(-1, nr, spr(:, is), uf, gf, cf)
          rolpi( iproj, io) = gf( nr)
        end do
      end do

      allocate( evecfv( nmatmax, nstfv, nspnfv))
      allocate( apwalm( ngkmax, apwordmax, lmmaxapw, natmtot, nspnfv))
      allocate( projm( wf_nprojused, wf_nprojused), olpm( wf_nprojused, wf_nprojused))
      allocate( eval( wf_nprojused), evec( wf_nprojused, wf_nprojused), auxmat( wf_nprojused, wf_nprojused))
      if( allocated( wf_transform)) deallocate( wf_transform)
      allocate( wf_transform( wf_nprojused, wf_nprojused, nkpt))
      do ik = 1, nkpt
        ! get basis function coefficients and matching coefficients
        call getevecfv( vkl(:, ik), vgkl(:, :, :, ik), evecfv)
        call match( ngk(1, ik), gkc(:, 1, ik), tpgkc(:, :, &
             & 1, ik), sfacgk(:, :, 1, ik), apwalm(:, :, :, :, &
             & 1))

        ! projection matrix elements and overlap matrix elements
        do iproj = 1, wf_nprojused
          is = wf_projst( wf_projused( iproj), 1)
          nr = nrmt( is)
          ias = idxas( wf_projst( wf_projused( iproj), 2), is)
          lm = idxlm( wf_projst( wf_projused( iproj), 4), wf_projst( wf_projused( iproj), 5))
          do iband = 1, wf_nprojused
            projm( iband, iproj) = zzero
            do ig = 1, ngk( 1, ik)
              auxc = zzero
              do io = 1, apword( wf_projst( wf_projused( iproj), 4), is)
                auxc = auxc + conjg( apwalm( ig, io, lm, ias, 1))*rolpi( iproj, io)
              end do
              projm( iband, iproj) = projm( iband, iproj) + conjg( evecfv( ig, bandstart+iband-1, 1))*auxc
            end do
            do io = 1, nlorb( is)
              l = lorbl( io, is)
              if( l .eq. wf_projst( wf_projused( iproj), 4)) then
                projm( iband, iproj) = projm( iband, iproj) + conjg( evecfv( ngk( 1, ik)+idxlo( lm, io, ias), bandstart+iband-1, 1))*ololo( io, wf_projst( wf_projused( iproj), 3), ias)
              end if
            end do
          end do
        end do 

        olpm(:,:) = zzero
        call ZGEMM( 'C', 'N', wf_nprojused, wf_nprojused, wf_nprojused, zone, projm, wf_nprojused, projm, wf_nprojused, zzero, olpm, wf_nprojused)

        ! transformation matrices
        call diaghermat( wf_nprojused, olpm, eval, evec)
        auxmat(:,:) = zzero
        do io = 1, wf_nprojused
          if( eval( io) .eq. 0d0) then
            write(*,*) "ERROR (genwf): 0 eigenvalue ", ik, io
          else
            auxmat( io, io) = zone/sqrt( cmplx( eval( io), 0, 8)) 
          end if
        end do

        call ZGEMM( 'N', 'N', wf_nprojused, wf_nprojused, wf_nprojused, zone, evec, wf_nprojused, auxmat, wf_nprojused, zzero, olpm, wf_nprojused)
        call ZGEMM( 'N', 'C', wf_nprojused, wf_nprojused, wf_nprojused, zone, olpm, wf_nprojused, evec, wf_nprojused, zzero, auxmat, wf_nprojused)
        call ZGEMM( 'N', 'N', wf_nprojused, wf_nprojused, wf_nprojused, zone, projm, wf_nprojused, auxmat, wf_nprojused, zzero, wf_transform( :, :, ik), wf_nprojused)
      end do 
      deallocate( rolpi, projm, olpm, auxmat, eval, evec, apwalm, evecfv, uf, gf, cf)
      return
    end subroutine genwf
    !EOC
    
    !BOP
    ! !ROUTINE: genmlwf
    ! !INTERFACE:
    !
    subroutine genmlwf( bandstart, nband, loproj)
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
      !   Created September 2016 (ST)
      !EOP
      !BOC

      integer, intent( in) :: bandstart, nband, loproj( nband) 
      
      ! local variables
      integer :: ngridk(3), nvec(3), nn, ix, iy, iz, i, ib, ik
      character :: read_file
      logical :: file_exists, reading_succes
      real(8) :: bwgtsm, mixing, maxdiff, vec1(3), vec2(3)

      ! allocatable arrays
      integer, allocatable :: idxn(:,:), diffk(:)
      real(8), allocatable :: nvl(:,:), nvc(:,:), bwgt(:), ravg(:,:), diffval(:,:)
      complex(8), allocatable :: mlwf_emat(:,:), auxmat(:,:), mlwf_m0(:,:,:,:), mlwf_m(:,:,:,:), eval(:), evec(:,:), mlwf_r(:,:), mlwf_t(:,:), mlwf_dw(:,:), mlwf_transform(:,:,:)
      complex(8), allocatable :: evecfv1(:,:,:), evecfv2(:,:,:)

      ! write next neighbours of each k-point to array
      ngridk = input%groundstate%ngridk
      nn = 0
      do ix = 1, 3
        nn = nn + 2 - 2/ngridk( ix)
      end do
      allocate( nvl( 3, nn), nvc( 3, nn), idxn( nn, nkpt))
      nvl(:,:) = 0d0
      i = 1
      do ix = 1, 3
        if( ngridk( ix) .eq. 2) then
          nvl( ix, i) = 0.5d0
          i = i + 1
        else if( ngridk( ix) .gt. 2) then
          nvl( ix, i) = 1d0/ngridk( ix)
          nvl( ix, i+1) = -1d0/ngridk( ix)
          i = i + 2
        end if
      end do
      do iz = 0, ngridk( 3)-1
        do iy = 0, ngridk( 2)-1
          do ix = 0, ngridk( 1)-1
            nvec(:) = (/ix, iy, iz/)
            ik = modulo( iz, ngridk( 3))*ngridk( 2)*ngridk( 1) + &
              modulo( iy, ngridk( 2))*ngridk( 1) + modulo( ix, ngridk(1))+1
            ib = 1
            do i = 1, 3
              if( ngridk( i) .eq. 2) then
                nvec( i) = nvec( i) + 1
                idxn( ib, ik) = modulo( nvec(3), ngridk( 3))*ngridk( 2)*ngridk( 1) + &
                  modulo( nvec(2), ngridk( 2))*ngridk( 1) + modulo( nvec(1), ngridk(1))+1
                nvec( i) = nvec( i) - 1
                ib = ib + 1
              else if( ngridk( i) .gt. 2) then
                nvec( i) = nvec( i) + 1
                idxn( ib, ik) = modulo( nvec(3), ngridk( 3))*ngridk( 2)*ngridk( 1) + &
                  modulo( nvec(2), ngridk( 2))*ngridk( 1) + modulo( nvec(1), ngridk(1))+1
                ib = ib + 1
                nvec( i) = nvec( i) - 2
                idxn( ib, ik) = modulo( nvec(3), ngridk( 3))*ngridk( 2)*ngridk( 1) + &
                  modulo( nvec(2), ngridk( 2))*ngridk( 1) + modulo( nvec(1), ngridk(1))+1
                nvec( i) = nvec( i) + 1
                ib = ib + 1
              end if
            end do
          end do
        end do
      end do
      do ib = 1, nn
        call r3mv( bvec, nvl( :, ib), nvc( :, ib))
      end do
      
      read_file = 'n'
      ! check wether WANNIER_MLWF.OUT exists and ask for use
      inquire( file='WANNIER_MLWF.OUT', exist=file_exists)
      if( .not. file_exists) then
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
        call genwf( bandstart, nband, loproj)
        ! calculating inner products <u(m,k)|u(n,k+b)>
        allocate( mlwf_emat( wf_nprojused, wf_nprojused))
        !allocate( mlwf_emat( nstfv, nstfv))
        allocate( auxmat( wf_nprojused, wf_nprojused))
        allocate( mlwf_m0( wf_nprojused, wf_nprojused, nn, nkpt))
        write(*,*) 'Computing inner products...'
        ! call before parallel loop in order to initialize save variables in emat_wannier

        allocate( evecfv1( nmatmax, nstfv, nspnfv), evecfv2( nmatmax, nstfv, nspnfv))

        ix = 0
        do ib = 1, nn !!ATTENTION!!nn
          call emat_init( nvl( :, ib), (/0d0, 0d0, 0d0/), input%groundstate%lmaxapw, 8)
          do ik = 1, nkpt
            
            call getevecfv( vkl( :, ik), vgkl( :, :, :, ik), evecfv1(:,:,1))
            call getevecfv( vkl( :, idxn( ib, ik)), vgkl( :, :, :, idxn( ib, ik)), evecfv2(:,:,1))
            !call emat_genemat( ik, idxn( ib, ik), 1, 4, 5, 6, evecfv1(:,:,1), evecfv2(:,:,1), mlwf_emat)
            call emat_genemat( ik, idxn( ib, ik), wf_bandstart, wf_nprojused, wf_bandstart, wf_nprojused, evecfv1(:,:,1), evecfv2(:,:,1), mlwf_emat)

            !call emat_wannier( ik, nvl( :, ib), wf_bandstart, wf_nprojused, wf_bandstart, wf_nprojused, mlwf_emat)
            !mlwf_emat = conjg( mlwf_emat)

            !write(*,*) '------------------------'
            !write(*,*) ib, ik
            !do iy = 1, wf_nprojused
            !  do iz = 1, wf_nprojused
            !    write(*,'(2(I3,3x),SP,E23.16,E23.16,"i")') iy, iz, mlwf_emat( iy, iz)
            !  end do
            !end do
            !call emat_wannier( ik, nvl( :, ib), 1, nstfv, 1, nstfv, mlwf_emat)
            
            call ZGEMM( 'N', 'N', wf_nprojused, wf_nprojused, wf_nprojused, zone, mlwf_emat, wf_nprojused, wf_transform( :, :, idxn( ib, ik)), wf_nprojused, zzero, auxmat, wf_nprojused)
            call ZGEMM( 'C', 'N', wf_nprojused, wf_nprojused, wf_nprojused, zone, wf_transform( :, :, ik), wf_nprojused, auxmat, wf_nprojused, zzero, mlwf_m0( :, :, ib, ik), wf_nprojused)
            ix = ix + 1
            ! print out progress
            write( 6, '(a,10x,i3,a)', advance='no') achar( 13), nint( 100.0d0*ix/nkpt/nn), '%'
            flush( 6)
          end do
          !write(*,*) ib
        end do
        write(6,*)

        !stop

        ! calculating weights of b-vectors
        allocate( bwgt( nn))
        bwgtsm = 0d0
        do ib = 1, nn
          bwgt( ib) = 3d0/(nn*norm2( nvc( :, ib))**2)
          bwgtsm = bwgtsm + bwgt( ib)
        end do
        
        ! mixing parameter for self consistent minimization
        mixing = dble( 0.5)
    
        ! initialize transformation matrices
        allocate( mlwf_transform( wf_nprojused, wf_nprojused, nkpt))
        mlwf_transform = zzero
        do ix = 1, wf_nprojused
          mlwf_transform( ix, ix, :) = zone
        end do

        ! start minimization loop
        allocate( evec( wf_nprojused, wf_nprojused), eval( wf_nprojused))
        allocate( mlwf_r( wf_nprojused, wf_nprojused), mlwf_t( wf_nprojused, wf_nprojused), mlwf_dw( wf_nprojused, wf_nprojused))
        allocate( mlwf_m( wf_nprojused, wf_nprojused, nn, nkpt))
        allocate( ravg( 3, wf_nprojused))
        allocate( diffval( 0:5000, nkpt), diffk( 0:5000))
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
            do ik = 1, nkpt
              do ib = 1, nn
                ravg( :, ix) = ravg( :, ix) - bwgt( ib)/nkpt*real( aimag( log( mlwf_m( ix, ix, ib, ik))))*nvc( :, ib)
              end do
            end do
          end do
      
          do ik = 1, nkpt
            mlwf_dw(:,:) = zzero
            do ib = 1, nn
              ! calculating R and T
              do ix = 1, wf_nprojused
                mlwf_r( :, ix) = mlwf_m( :, ix, ib, ik)*conjg( mlwf_m( ix, ix, ib, ik))
                mlwf_t( :, ix) = mlwf_m( :, ix, ib, ik)/mlwf_m( ix, ix, ib, ik)*(real( aimag( log( mlwf_m( ix, ix, ib, ik)))) + dot_product( nvc( :, ib), ravg( :, ix)))
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
            call ZGEMM( 'N', 'N', wf_nprojused, wf_nprojused, wf_nprojused, zone, evec, wf_nprojused, auxmat, wf_nprojused, zzero, mlwf_dw, wf_nprojused)
            call ZGEMM( 'N', 'C', wf_nprojused, wf_nprojused, wf_nprojused, zone, mlwf_dw, wf_nprojused, evec, wf_nprojused, zzero, auxmat, wf_nprojused)
            call ZGEMM( 'N', 'N', wf_nprojused, wf_nprojused, wf_nprojused, zone, mlwf_transform( :, :, ik), wf_nprojused, auxmat, wf_nprojused, zzero, mlwf_dw, wf_nprojused)
            maxdiff = norm2( abs( reshape( mlwf_transform( :, :, ik) - mlwf_dw, (/wf_nprojused**2/))))!/norm2( abs( reshape( mlwf_transform( :, :, ik), (/wf_nprojused**2/))))
            if( maxdiff .ge. maxval( diffval( iz, :))) diffk( iz) = ik
            diffval( iz, ik) = maxdiff

            mlwf_transform( :, :, ik) = mlwf_dw(:,:)
          end do
    
          ! updating M
          do ik = 1, nkpt
            do ib = 1, nn
              call ZGEMM( 'N', 'N', wf_nprojused, wf_nprojused, wf_nprojused, zone, mlwf_m0( :, :, ib, ik), wf_nprojused, mlwf_transform( :, :, idxn( ib, ik)), wf_nprojused, zzero, auxmat, wf_nprojused)
              call ZGEMM( 'C', 'N', wf_nprojused, wf_nprojused, wf_nprojused, zone, mlwf_transform( :, :, ik), wf_nprojused, auxmat, wf_nprojused, zzero, mlwf_m( :, :, ib, ik), wf_nprojused)
            end do
          end do

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
        do ik = 1, nkpt
          call ZGEMM( 'N', 'N', wf_nprojused, wf_nprojused, wf_nprojused, zone, wf_transform( :, :, ik), wf_nprojused, mlwf_transform( :, :, ik), wf_nprojused, zzero, auxmat, wf_nprojused)
          wf_transform( :, :, ik) = auxmat(:,:)
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
end module mod_wannier
