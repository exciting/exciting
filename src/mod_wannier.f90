module mod_wannier
  use mod_APW_LO
  use mod_atoms 
  use mod_kpoint
  use mod_constants
  use mod_muffin_tin
  use mod_Gvector
  use mod_Gkvector
  use mod_eigensystem, only: nmatmax, idxlo
  use mod_spin, only: nspnfv
  use mod_lattice, only: bvec
  use mod_eigenvalue_occupancy, only: nstfv
  !use mod_symmetry
  use m_ematqk
  use m_plotmat
  use modmpi
  use mod_misc

  implicit none

! module variables
  integer :: wf_nprojtot, wf_bandstart, wf_nband
  integer :: wf_n_usedshells                            ! number of shells used for gradient calculation
  logical :: wf_initialized = .false.
  
  integer, allocatable :: wf_projst(:,:), wf_projused(:)
  complex(8), allocatable :: wf_transform(:,:,:)        ! unitary transformation matrices
  complex(8), allocatable :: wf_projection(:,:,:)       ! full overlap of wafefunctions and local-orbitals
  complex(8), allocatable :: wf_opf(:,:)                ! expansion coefficients for optimal projection functions
  complex(8), allocatable :: wf_emat(:,:,:,:,:)         ! plane-wave matrix elements for neighboring k-points
  real(8), allocatable :: wf_centers(:,:)               ! centers of Wannier functions
  real(8), allocatable :: wf_omega(:)                   ! localization functional for each Wannier function
  real(8), allocatable :: wf_n_dist(:)                  ! distance of neighbors in given shell
  integer, allocatable :: wf_n_n(:)                     ! number of neighbors in given shell
  integer, allocatable :: wf_n_ik(:,:,:)                ! k-point index of given neighbor in given shell for given k-point
  real(8), allocatable :: wf_n_vl(:,:,:), wf_n_vc(:,:,:)! vector of given neighbor in given shell in lattice and cartesian coordinates
  real(8), allocatable :: wf_n_wgt(:)                   ! geometric weight of given shell

! methods
  contains
    !BOP
    ! !ROUTINE: wannier_init
    ! !INTERFACE:
    !
    subroutine wannier_init
      ! !USES:
      ! !DESCRIPTION:
      !   Reads local-orbitals from species files which are indicated to be used
      !   as projection functions for the generation of Wannier functions to the
      !   module {\tt mod\_wannier} and constructs necessary geometry (nearest 
      !   neighbors and geometric weights) as well as the overlap of the
      !   wavefunctions and the local-orbitals.
      !
      ! !REVISION HISTORY:
      !   Created September 2016 (Tillack)
      !EOP
      !BOC

      ! local variables
      integer :: iproj, is, ia, js, ja, ias, nr, l, m, lm, io, ilo, ir, ig
      integer :: iknr, ngknr
      integer, allocatable :: igkignr(:)
      real(8), allocatable :: vgklnr(:,:,:), vgkcnr(:,:,:), gkcnr(:), tpgkcnr(:,:), sval(:)
      complex(8), allocatable :: sfacgknr(:,:)
      complex(8), allocatable :: evecfv(:,:,:), apwalm(:,:,:,:,:), auxmat(:,:), lsvec(:,:), rsvec(:,:) 
      real(8), allocatable :: rolpi(:,:), uf(:), gf(:), cf(:,:)
      
      call init0
      call init1

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
        write(*,*) 'ERROR (wannier_init): No local-orbitals found for projection.'
        stop
      end if

      ! find geometry
      call wannier_geometry
      
      ! radial overlap integrals
      call readstate
      call linengy
      call genapwfr     ! APW radial functions
      call genlofr      ! LO radial functions

      allocate( rolpi( wf_nprojtot, apwordmax+maxlorb))
      allocate( uf( nrmtmax), gf( nrmtmax), cf( 3, nrmtmax))
      rolpi(:,:) = zzero
      do iproj = 1, wf_nprojtot
        uf(:) = 0d0
        is = wf_projst( iproj, 1)
        nr = nrmt( is)
        ias = idxas( wf_projst( iproj, 2), is)
        lm = idxlm( wf_projst( iproj, 4), wf_projst( iproj, 5))
        do io = 1, apword( wf_projst( iproj, 4), is)
          do ir = 1, nr
            uf( ir) = apwfr( ir, 1, io, wf_projst( iproj, 4), ias)*lofr( ir, 1, wf_projst( iproj, 3), ias)*spr( ir, is)**2
          end do
          call fderiv( -1, nr, spr(:, is), uf, gf, cf)
          rolpi( iproj, io) = gf( nr)
        end do
        do ilo = 1, nlorb( is)
          l = lorbl( ilo, is)
          if( l .eq. wf_projst( iproj, 4)) then
            do ir = 1, nr
              uf( ir) = lofr( ir, 1, ilo, ias)*lofr( ir, 1, wf_projst( iproj, 3), ias)*spr( ir, is)**2
            end do
            call fderiv( -1, nr, spr( :, is), uf, gf, cf)
            rolpi( iproj, apwordmax+ilo) = gf( nr)
          end if
        end do
      end do

      allocate( evecfv( nmatmax, nstfv, nspnfv))
      allocate( apwalm( ngkmax, apwordmax, lmmaxapw, natmtot, nspnfv))
      allocate( wf_projection( nstfv, wf_nprojtot, nkptnr))
      allocate( auxmat( ngkmax+nlotot, wf_nprojtot))
      allocate( sval( min( nstfv, wf_nprojtot)), &
                lsvec( nstfv, nstfv), &
                rsvec( wf_nprojtot, wf_nprojtot))
      allocate( igkignr( ngkmax))
      allocate( vgklnr( 3, ngkmax, nspnfv), vgkcnr( 3, ngkmax, nspnfv), gkcnr( ngkmax), tpgkcnr( 2, ngkmax))
      allocate( sfacgknr( ngkmax, natmtot))

      do iknr = 1, nkptnr
        vgklnr = zzero
        ! find G+k-vectors for non-reduced k-point
        call gengpvec( vklnr( :, iknr), vkcnr( :, iknr), ngknr, igkignr, vgklnr(:,:,1), vgkcnr(:,:,1), gkcnr, tpgkcnr)
        ! find structure factors
        call gensfacgp( ngknr, vgkcnr, ngkmax, sfacgknr)
        ! get basis function coefficients and matching coefficients
        call getevecfv( vklnr(:, iknr), vgklnr, evecfv)
        call match( ngknr, gkcnr, tpgkcnr, sfacgknr, apwalm(:, :, :, :, 1))

        ! projection matrices 
        do iproj = 1, wf_nprojtot
          is = wf_projst( iproj, 1)
          ias = idxas( wf_projst( iproj, 2), is)
          lm = idxlm( wf_projst( iproj, 4), wf_projst( iproj, 5))
          do ig = 1, ngknr
            auxmat( ig, iproj) = zzero
            do io = 1, apword( wf_projst( iproj, 4), is)
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

        call ZGEMM( 'C', 'N', nstfv, wf_nprojtot, ngknr+nlotot, zone, &
             evecfv( 1:(ngknr+nlotot), :, 1), ngknr+nlotot, &
             auxmat( 1:(ngknr+nlotot), :), ngknr+nlotot, zzero, &
             wf_projection( :, :, iknr), nstfv)
      end do 

      wf_initialized = .true.
      call wannier_showproj
      deallocate( auxmat, rolpi, apwalm, evecfv, uf, gf, cf, sfacgknr, vgklnr, vgkcnr, gkcnr, tpgkcnr, igkignr, sval, lsvec, rsvec)
      return
    end subroutine wannier_init
    !EOC

    !BOP
    ! !ROUTINE: wannier_emat
    ! !INTERFACE:
    !
    subroutine wannier_emat
      ! !USES:
      ! !DESCRIPTION:
      !   Generates plane-wave matrix elements for neighboring k-points
      !   needed in order to calculate well/maximally localized Wannier
      !   functions.
      ! !REVISION HISTORY:
      !   Created September 2016 (Tillack)
      !EOP
      !BOC

      ! local variables
      integer :: is, n, k1, k2, iknr, ngknr

      ! allocatable arrays
      integer, allocatable :: igkignr(:)
      real(8), allocatable :: vgklnr(:,:,:), vgkcnr(:,:,:), gkcnr(:), tpgkcnr(:,:)
      complex(8), allocatable :: evecfv1(:,:,:), evecfv2(:,:,:)
      
      write(*,*) "wannier_emat..."

      if( allocated( wf_emat)) deallocate( wf_emat)
      allocate( wf_emat( nstfv, nstfv, maxval( wf_n_n( 1:wf_n_usedshells)), wf_n_usedshells, nkptnr))
      wf_emat = zzero

      allocate( igkignr( ngkmax))
      allocate( vgklnr( 3, ngkmax, nspnfv), vgkcnr( 3, ngkmax, nspnfv), gkcnr( ngkmax), tpgkcnr( 2, ngkmax))
      allocate( evecfv1( nmatmax, nstfv, nspnfv), evecfv2( nmatmax, nstfv, nspnfv))

      do is = 1, wf_n_usedshells
        do n = 1, wf_n_n( is)     
          call emat_init( wf_n_vl( :, n, is), (/0, 0, 0/), input%groundstate%lmaxapw, 8)
          k1 = 1
          k2 = nkptnr
#ifdef MPI
          k1 = firstofset( rank, nkptnr)
          k2 = lastofset( rank, nkptnr)
#endif
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, ngknr, igkignr, vgklnr, vgkcnr, gkcnr, tpgkcnr, evecfv1, evecfv2)
!$OMP DO  
#endif
          do iknr = k1, k2   
            ! find G+k-vecto rs and eigenvectors for non-reduced k-point k
            call gengpvec( vklnr( :, iknr), vkcnr( :, iknr), ngknr, igkignr, vgklnr(:,:,1), vgkcnr(:,:,1), gkcnr, tpgkcnr)
            call getevecfv( vklnr( :, iknr), vgklnr, evecfv1)
            ! find G+k+b-vectors and eigenvectors for non-reduced k-point k+b
            call gengpvec( vklnr( :, wf_n_ik( n, is, iknr)), vkcnr( :, wf_n_ik( n, is, iknr)), ngknr, igkignr, vgklnr(:,:,1), vgkcnr(:,:,1), gkcnr, tpgkcnr)
            call getevecfv( vklnr( :, wf_n_ik( n, is, iknr)), vgklnr, evecfv2)
            ! generate plane-wave matrix elements
            call emat_genemat( iknr, 1, nstfv, 1, nstfv, evecfv1(:,:,1), evecfv2(:,:,1), wf_emat( :, :, n, is, iknr))
          end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
#ifdef MPI
        call mpi_allgatherv_ifc( nkptnr, wf_nband**2, zbuf=mlwf_m0(:,:,:,ib))
        call mpi_allreduce( omegai, omegaod, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
        end do
      end do
      call emat_destroy
      return
      !EOC
    end subroutine wannier_emat
    !EOP

    !BOP
    ! !ROUTINE: wannier_gen_pro
    ! !INTERFACE:
    !
    subroutine wannier_gen_pro( bandstart, nband, loproj)
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
      !   Created September 2016 (Tillack)
      !EOP
      !BOC

      integer, intent( in) :: bandstart, nband, loproj( nband)

      ! local variables
      integer :: iproj, iknr, i
      real(8), allocatable :: sval(:), ravg(:,:), omega(:)
      complex(8), allocatable :: projm(:,:), lsvec(:,:), rsvec(:,:)

      if( .not. wf_initialized) call wannier_init
      wf_nband = nband
      wf_bandstart = bandstart
      if( allocated( wf_projused)) deallocate( wf_projused)
      allocate( wf_projused( nband))

      ! check projection orbitals
      do iproj = 1, nband
        if( (loproj( iproj) .lt. 1) .or. (loproj( iproj) .gt. wf_nprojtot)) then
          write(*,'(a,i2,a)') ' ERROR (wannier_gen_pro): ', loproj( iproj), ' is not a valid index for projection local-orbitals.'
          write(*,*) 'Here is a list of local-orbitals that can be used for projection:'
          call wannier_showproj
          return
        end if
        wf_projused( iproj) = loproj( iproj)
      end do

      allocate( projm( wf_nband, wf_nband), sval( wf_nband), lsvec( wf_nband, wf_nband), rsvec( wf_nband, wf_nband))
      if( allocated( wf_transform)) deallocate( wf_transform)
      allocate( wf_transform( wf_nband, wf_nband, nkptnr))
         
#ifdef USEOMP                
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, i, projm, sval, lsvec, rsvec)
!!$OMP DO  
#endif
      do iknr = 1, nkptnr
        ! build selected projection matrix
        do i = 1, wf_nband
          projm( :, i) = wf_projection( wf_bandstart:(wf_bandstart+wf_nband-1), loproj( i), iknr)
        end do
        ! build transformation matrices
        call zgesdd_wrapper( projm, wf_nband, wf_nband, sval, lsvec, rsvec)
        if( minval( sval) .lt. input%structure%epslat) then
          write(*,*) ' WARNING (wannier_gen_pro): The local-orbitals selected for projection may not describe the angular character of the selected bands sufficiently. Add local-orbitals with different angular momentum.'
        end if
        call ZGEMM( 'N', 'N', wf_nband, wf_nband, wf_nband, zone, &
             lsvec, wf_nband, &
             rsvec, wf_nband, zzero, &
             wf_transform( :, :, iknr), wf_nband)
      end do 
#ifdef USEOMP
!!$OMP END DO
!!$OMP END PARALLEL
#endif
      allocate( omega( wf_nband), ravg( 3, wf_nband))
      call wannier_loc
      write(*,'("Omega: ",F13.6)') sum( wf_omega)
      write(*,'(100F13.6)') wf_omega
      write(*,'(100F13.6)') wf_centers( 1, :)
      write(*,'(100F13.6)') wf_centers( 2, :)
      write(*,'(100F13.6)') wf_centers( 3, :)
      call wannier_writefile( 'WANNIER')
      deallocate( projm, sval, lsvec, rsvec)
      return
    end subroutine wannier_gen_pro
    !EOC
    
    !BOP
    ! !ROUTINE: wannier_gen_opf
    ! !INTERFACE:
    !
    subroutine wannier_gen_opf( bandstart, nband)
      ! !USES:
      ! !INPUT/OUTPUT PARAMETERS:
      !   bandstart : n1, lowest band index of the band range used for generation of
      !               Wannier functions (in,integer)
      !   nband     : N, number of bands used for generation (in,integer)
      ! !DESCRIPTION:
      !   Does the same thing as {\tt genwf} does but the transformation
      !   matrices are used for the generation of maximally localized Wannier
      !   functions. The matrices are computed in a self consistent loop. 
      !
      ! !REVISION HISTORY:
      !   Created September 2016 (Tillack)
      !EOP
      !BOC
      integer, intent( in) :: bandstart, nband

      ! local variables
      integer :: iknr, is, n, nx, iproj, i, j, l, ix, loproj( nband), minit
      real(8) :: lambda, phi, phimean, uncertainty, theta, evalmin, omegamin, limit
      real(8) :: opf_min_q(3,3), opf_min_p(3,1), opf_min_c, opf_min_sol(3), eval(3), revec(3,3), opf_min_big(6,6), opf_min_aux(3,3), tp(2)
      complex(8) :: opf_min_z(3,1), evec(3,3), evec_big(6,6), eval_big(6), r1, r2

      ! allocatable arrays
      complex(8), allocatable :: opf_transform(:,:,:), lsvec(:,:), rsvec(:,:), opf_mixing(:,:), opf_x(:,:,:), opf_x0(:,:,:)
      complex(8), allocatable :: auxmat(:,:), auxmat2(:,:), projm(:,:), lsvec2(:,:), rsvec2(:,:)
      real(8), allocatable :: sval(:), opf_t(:), ravg(:,:), omega(:), sval2(:), phihist(:), phimeanhist(:)

      if( .not. wf_initialized) call wannier_init
      wf_nband = nband
      wf_bandstart = bandstart
      if( allocated( wf_projused)) deallocate( wf_projused)
      allocate( wf_projused( wf_nprojtot))

      do iproj = 1, wf_nprojtot
        wf_projused( iproj) = iproj
      end do

      ! build rectangular transformation matrices
      allocate( opf_transform( wf_nband, wf_nprojtot, nkptnr))
      allocate( lsvec( wf_nband, wf_nband), rsvec( wf_nprojtot, wf_nprojtot), sval( wf_nband))
      do iknr = 1, nkptnr
        call zgesdd_wrapper( wf_projection( wf_bandstart:(wf_bandstart+wf_nband-1), :, iknr), wf_nband, wf_nprojtot, sval, lsvec, rsvec)
        call zgemm( 'N', 'N', wf_nband, wf_nprojtot, wf_nband, zone, &
             lsvec, wf_nband, &
             rsvec( 1:wf_nband, :), wf_nband, zzero, &
             opf_transform( :, :, iknr), wf_nband)
      end do

      ! generate plane-wave matrix-elements
      !write(*,*) 'Computing inner products...'
      call wannier_emat
      write(*,*) "done!"

      lambda = 0.5d0*dot_product( dble( wf_n_n( 1:wf_n_usedshells)), wf_n_wgt( 1:wf_n_usedshells))
      minit = 100
      limit = min( 1.d-3, 1.d3*input%properties%wannier%uncertainty)

      !write(*,*) 'Preparation ...'
      allocate( opf_x( wf_nprojtot, wf_nprojtot, nkptnr*(1+sum( wf_n_n( 1:wf_n_usedshells)))))
      allocate( opf_x0( wf_nprojtot, wf_nprojtot, nkptnr*(1+sum( wf_n_n( 1:wf_n_usedshells)))))
      allocate( opf_t( nkptnr*(1+sum( wf_n_n( 1:wf_n_usedshells)))))
      ! build enlarged overlaps
      allocate( auxmat( wf_nband, wf_nprojtot), auxmat2( wf_nprojtot, wf_nprojtot))
      i = 0
      do iknr = 1, nkptnr
        do is = 1, wf_n_usedshells
          do n = 1, wf_n_n( is)
            i = i + 1
            call zgemm( 'N', 'N', wf_nband, wf_nprojtot, wf_nband, zone, &
                 wf_emat( wf_bandstart:(wf_bandstart+wf_nband-1), wf_bandstart:(wf_bandstart+wf_nband-1), n, is, iknr), wf_nband, &
                 opf_transform( :, :, wf_n_ik( n, is, iknr)), wf_nband, zzero, &
                 auxmat, wf_nband)
            call zgemm( 'C', 'N', wf_nprojtot, wf_nprojtot, wf_nband, zone, &
                 opf_transform( :, :, iknr), wf_nband, &
                 auxmat, wf_nband, zzero, &
                 opf_x0( :, :, i), wf_nprojtot)
            opf_t( i) = -wf_n_wgt( is)
          end do
        end do
      end do
      ! build constraint matrices
      do iknr = 1, nkptnr
        i = i + 1
        call zgemm( 'C', 'N', wf_nprojtot, wf_nprojtot, wf_nband, zone, &
             wf_projection( wf_bandstart:(wf_bandstart+wf_nband-1), :, iknr), wf_nband, &
             wf_projection( wf_bandstart:(wf_bandstart+wf_nband-1), :, iknr), wf_nband, zzero, &
             opf_x0( :, :, i), wf_nprojtot)
        do is = 1, wf_nprojtot
          opf_x0( is, is, i) = opf_x0( is, is, i) - zone
        end do
        opf_t( i) = lambda
      end do
      nx = i

      allocate( wf_opf( wf_nprojtot, wf_nband))
      allocate( opf_mixing( wf_nprojtot, wf_nprojtot))
      allocate( omega( wf_nband), ravg( 3, wf_nband))
      allocate( projm( wf_nband, wf_nband), sval2( wf_nband), lsvec2( wf_nband, wf_nband), rsvec2( wf_nband, wf_nband))
      allocate( wf_transform( wf_nband, wf_nband, iknr))
      allocate( phihist( minit), phimeanhist( minit))
      deallocate( lsvec, rsvec)
      allocate( lsvec(3,3), rsvec(3,3))
      omegamin = 1.d13
      opf_mixing = zzero
      opf_x = opf_x0
      do i = 1, wf_nprojtot
        opf_mixing( i, i) = zone
      end do
      l = 1
      i = 1
      n = 0
      phihist = 0.d0
      phimeanhist = 0.d0
      uncertainty = 1.d0
      ! start minimization
      do while( (n .lt. minit) .or. ((n .lt. 500000) .and. (uncertainty .gt. limit)))
        i = 1
        do while( (i .le. wf_nband) .and. ((n .lt. minit) .or. (uncertainty .gt. limit)))
          j = i + 1
          do while( (j .le. wf_nprojtot) .and. ((n .lt. minit) .or. (uncertainty .gt. limit)))
            n = n + 1
            ! calculate Lagrangian
            phi = 0.d0
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( is, iknr) reduction(+:phi)
!$OMP DO  
#endif
            do is = 1, nx
              do iknr = 1, wf_nband
                phi = phi + opf_t( is)*opf_x( iknr, iknr, is)*conjg( opf_x( iknr, iknr, is))
              end do
            end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            phi = phi/nkptnr
            phi = phi + wf_nband*dot_product( wf_n_n( 1:wf_n_usedshells), wf_n_wgt( 1:wf_n_usedshells))

            ! convergence analysis
            if( n .eq. 1) then
              phihist(:) = phi
              phimeanhist(:) = phi
            end if
            phihist = cshift( phihist, -1)
            phimeanhist = cshift( phimeanhist, -1)
            phihist(1) = phi
            phimean = sum( phihist(:))/minit
            phimeanhist(1) = phimean
            uncertainty = sqrt( sum( (phihist(:)-phimeanhist(:))**2)/(minit-1))/phi

            omega = 0.d0
            ! calculate omega
!#ifdef USEOMP                
!!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, projm, sval2, lsvec2, rsvec2)
!!$OMP DO  
!#endif
!            do iknr = 1, nkptnr
!              call zgemm( 'N', 'N', wf_nband, wf_nband, wf_nprojtot, zone, &
!                   wf_projection( wf_bandstart:(wf_bandstart+wf_nband-1), :, iknr), wf_nband, &
!                   opf_mixing( :, 1:wf_nband), wf_nprojtot, zzero, &
!                   projm, wf_nband)
!              call zgesdd_wrapper( projm, wf_nband, wf_nband, sval2, lsvec2, rsvec2)
!              call ZGEMM( 'N', 'N', wf_nband, wf_nband, wf_nband, zone, &
!                   lsvec2, wf_nband, &
!                   rsvec2, wf_nband, zzero, &
!                   wf_transform( :, :, iknr), wf_nband)
!            end do
!#ifdef USEOMP
!!$OMP END DO
!!$OMP END PARALLEL
!#endif
!            call wannier_loc
!            write(*,'(4(I4,3x),3(F13.6,3x))') n, l, i, j, phi, sum( wf_omega), uncertainty
             write(*,'(4(I4,3x),2(F13.6,3x))') n, l, i, j, phi, uncertainty
!            if( sum( wf_omega) .lt. omegamin) then
!              omegamin = sum( wf_omega)
!              wf_opf = opf_mixing( :, 1:wf_nband)
!            end if


            opf_min_q = 0.d0
            opf_min_p = 0.d0
            opf_min_c = 0.d0
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( is, opf_min_z) reduction(+: opf_min_q, opf_min_p)
!$OMP DO  
#endif
            do is = 1, nx
              opf_min_z(1,1) = 0.5d0*(opf_x( i, i, is) - opf_x( j, j, is))
              opf_min_z(2,1) = -0.5d0*(opf_x( i, j, is) + opf_x( j, i, is))
              opf_min_z(3,1) = 0.5d0*zi*(opf_x( i, j, is) - opf_x( j, i, is))
              opf_min_q = opf_min_q + opf_t( is)*dble( matmul( opf_min_z, transpose( conjg( opf_min_z))))
              opf_min_p = opf_min_p + opf_t( is)*dble( conjg( opf_x( i, i, is) + opf_x( j, j, is))*opf_min_z)
              !opf_min_c = 0.25d0*opf_t( is)*abs( opf_x( i, i, is) + opf_x( j, j, is))**2
            end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            call diagsymmat( 3, opf_min_q, eval, revec)
            if( j .le. wf_nband) then
              opf_min_sol(:) = revec( :, 1)
              if( opf_min_sol(1) .lt. 0.d0) opf_min_sol = -opf_min_sol
              !write(*,*) "Case 1"
            else
              opf_min_big = 0.d0
              opf_min_big( 1:3, 1:3) = 2*opf_min_q
              opf_min_big( 1:3, 4:6) = 0.25d0*matmul( opf_min_p, transpose( opf_min_p)) - matmul( opf_min_q, opf_min_q)
              opf_min_big( 4:6, :) = 0.d0
              opf_min_big( 4, 1) = 1.d0
              opf_min_big( 5, 2) = 1.d0
              opf_min_big( 6, 3) = 1.d0
              call diaggenmat( 6, cmplx( opf_min_big, 0, 8), eval_big, evec_big)
              evalmin = maxval( abs( eval_big))
              do is = 1, 6
                if( abs( aimag( eval_big( is))) .lt. input%structure%epslat) evalmin = min( evalmin, real( eval_big( is)))
              end do
              do is = 1, 3
                opf_min_q( is, is) = opf_min_q( is, is) - evalmin
              end do
              if( minval( eval(:) - evalmin) .lt. 1.d-6) then
                call zgesdd_wrapper( cmplx( opf_min_q, 0.d0, 8), 3, 3, eval, lsvec, rsvec)
                do is = 1, 3
                  if( eval( is) .gt. 1.d-6) then
                    rsvec( is, :) = rsvec( is, :)/eval( is)
                  else
                    rsvec( is, :) = 0.d0
                  end if
                end do
                opf_min_aux = real( matmul( transpose( rsvec), transpose( lsvec)))
                call r3mv( opf_min_aux, -0.5d0*opf_min_p(:,1), opf_min_sol)
                !write(*,*) "Case 2"
                call r3mv( opf_min_q, opf_min_sol, opf_min_aux(:,1))
                if( (norm2( opf_min_aux(:,1) + 0.5d0*opf_min_p(:,1)) .ge. 1.d-6) .or. (norm2( opf_min_sol) .gt. 1.d0)) then
                  opf_min_sol = 10.d0
                else
                  opf_min_sol = opf_min_sol + revec(:,1)/norm2( revec(:,1))*sqrt( 1.d0 - norm2( opf_min_sol))
                end if
              else
                !write(*,*) "Case 3"
                call r3minv( opf_min_q, opf_min_aux)
                call r3mv( opf_min_aux, -0.5d0*opf_min_p(:,1), opf_min_sol)
              end if
            end if
            if( abs( 1.d0 - norm2( opf_min_sol)) .le. input%structure%epslat) then
              call sphcrd( (/opf_min_sol(2), opf_min_sol(3), opf_min_sol(1)/), r1, tp)
              phi = tp(2)
              theta = tp(1)/2
              r1 = cmplx( cos( theta), 0, 8)
              r2 = exp( zi*phi)*sin( theta)
              ! update mixing matrix
              auxmat(1,:) = opf_mixing(:,i)
              auxmat(2,:) = opf_mixing(:,j)
              opf_mixing(:,i) = auxmat(1,:)*r1 - auxmat(2,:)*conjg( r2)
              opf_mixing(:,j) = auxmat(2,:)*r1 + auxmat(1,:)*r2
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( is, auxmat)
!$OMP DO  
#endif
              ! update X-matrices
              do is = 1, nx
                auxmat(1,:) = opf_x( :, i, is)
                auxmat(2,:) = opf_x( :, j, is)
                opf_x( :, i, is) = auxmat(1,:)*r1 - auxmat(2,:)*conjg( r2)
                opf_x( :, j, is) = auxmat(2,:)*r1 + auxmat(1,:)*r2
                auxmat(1,:) = opf_x( i, :, is)
                auxmat(2,:) = opf_x( j, :, is)
                opf_x( i, :, is) = auxmat(1,:)*r1 - auxmat(2,:)*r2
                opf_x( j, :, is) = auxmat(2,:)*r1 + auxmat(1,:)*conjg( r2)
              end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            else
              !write(*,*) "Doh!"
            end if
            j = j + 1
          end do
          i = i + 1
        end do
        l = l + 1
      end do
      wf_opf = opf_mixing( :, 1:wf_nband)

      do iknr = 1, nkptnr
        call zgemm( 'N', 'N', wf_nband, wf_nband, wf_nprojtot, zone, &
             wf_projection( wf_bandstart:(wf_bandstart+wf_nband-1), :, iknr), wf_nband, &
             wf_opf, wf_nprojtot, zzero, &
             projm, wf_nband)
        call zgesdd_wrapper( projm, wf_nband, wf_nband, sval2, lsvec2, rsvec2)
        call ZGEMM( 'N', 'N', wf_nband, wf_nband, wf_nband, zone, &
             lsvec2, wf_nband, &
             rsvec2, wf_nband, zzero, &
             wf_transform( :, :, iknr), wf_nband)
      end do
      call wannier_loc
      write(*,*) sum( wf_omega)
      write(*,'(100F13.6)') wf_omega
      write(*,'(100F13.6)') wf_centers( 1, :)
      write(*,'(100F13.6)') wf_centers( 2, :)
      write(*,'(100F13.6)') wf_centers( 3, :)
      call wannier_writefile( 'WANNIER')
      return
      !EOC
    end subroutine wannier_gen_opf
    !EOP

    !BOP
    ! !ROUTINE: wannier_gen_max
    ! !INTERFACE:
    !
    subroutine wannier_gen_max( bandstart, nband, loproj, fromfile)
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
      !   Created September 2016 (Tillack)
      !EOP
      !BOC

      integer, intent( in) :: bandstart, nband
      integer, optional, intent( in) :: loproj( nband)
      logical, optional, intent( in) :: fromfile
      
      ! local variables
      integer :: ngridk(3), ix, iy, iz, i, n, is, iknr, ngknr, k1, k2, z0, z1
      real(8) :: mixing, alpha, grad
      real(8) :: omegastart, omega, omegamean, uncertainty
      real(8), allocatable :: omegai(:), omegad(:), omegaod(:), omega2(:)
      integer :: minit
      logical :: succes

      ! allocatable arrays
      integer, allocatable :: idxn(:,:), nn(:), integerlist(:)
      integer, allocatable :: igkignr(:)
      real(8), allocatable :: vgklnr(:,:,:), vgkcnr(:,:,:), gkcnr(:), tpgkcnr(:,:)
      real(8), allocatable :: ndist(:), nvl(:,:,:), nvc(:,:,:), bwgt(:), ravg(:,:)
      real(8), allocatable :: omegahist(:), omegameanhist(:), eval(:)
      complex(8), allocatable :: auxmat(:,:), mlwf_m0(:,:,:,:,:), mlwf_m(:,:,:,:,:), evec(:,:), mlwf_r(:,:), mlwf_t(:,:), mlwf_dw(:,:), mlwf_transform(:,:,:)
      complex(8), allocatable :: auxmatcpy(:,:), mlwf_dwcpy(:,:)
      complex(8), allocatable :: evecfv1(:,:,:), evecfv2(:,:,:)
      complex(8), allocatable :: auxmatread(:,:)

      ! generate transformation matrices from projection as initial input
      !write(*,*) 'Performing initial projection step...'
      succes = .true.
      if( present( fromfile)) then
        if( fromfile) then
          call wannier_init
          call wannier_readfile( 'WANNIER', succes)
          write(*,*) wf_nband
          write(*,*) wf_bandstart
          write(*,*) wf_nprojtot
        end if
      else if( present( loproj)) then
        call wannier_gen_pro( bandstart, nband, loproj)
      else
        call wannier_gen_opf( bandstart, nband)
      end if
      if( .not. succes) then
        call wannier_gen_opf( bandstart, nband)
      end if

      allocate( auxmat( wf_nband, wf_nband))
      allocate( mlwf_m0( wf_nband, wf_nband, maxval( wf_n_n( 1:wf_n_usedshells)), wf_n_usedshells, nkptnr))
      allocate( omegai( wf_nband))
      if( .not. allocated( wf_emat)) then
        write(*,*) 'Computing inner products...'
        call wannier_emat
      endif

      omegai = 0.d0

      do is = 1, wf_n_usedshells
        do n = 1, wf_n_n( is)     
          k1 = 1
          k2 = nkptnr
#ifdef MPI
          k1 = firstofset( rank, nkptnr)
          k2 = lastofset( rank, nkptnr)
#endif
#ifdef USEOMP                
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, ix, ngknr, auxmat) reduction(+:omegai)
!$OMP DO  
#endif
          do iknr = k1, k2   
                             
            call ZGEMM( 'N', 'N', wf_nband, wf_nband, wf_nband, zone, &
                 wf_emat( wf_bandstart:(wf_bandstart+wf_nband-1), wf_bandstart:(wf_bandstart+wf_nband-1), n, is, iknr), wf_nband, &
                 wf_transform( :, :, wf_n_ik( n, is, iknr)), wf_nband, zzero, &
                 auxmat, wf_nband)
            call ZGEMM( 'C', 'N', wf_nband, wf_nband, wf_nband, zone, &
                 wf_transform( :, :, iknr), wf_nband, &
                 auxmat, wf_nband, zzero, &
                 mlwf_m0( :, :, n, is, iknr), wf_nband)
  
            ! independent part of localization functional
            do ix = 1, wf_nband
              omegai( ix) = omegai( ix) + wf_n_wgt( is)*(1d0 - sum( abs( mlwf_m0( :, ix, n, is, iknr))**2))
            end do

          end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
#ifdef MPI
        call mpi_allgatherv_ifc( nkptnr, wf_nband**2, zbuf=mlwf_m0(:,:,:,ib))
        call mpi_allreduce( omegai, omegaod, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif
        end do
      end do
      
      omegai = omegai/nkptnr
#ifdef MPI
      omegai = omegaod/nkptnr
#endif
      !write(*,'("Omega I: ",F23.16)') sum( omegai)
      
      ! mixing parameter for self consistent minimization
      mixing = dble( 0.5)
      alpha = mixing
      ! minimum number of iterations
      minit = 100
#ifdef MPI
      if( rank .eq. 0) then
#endif
        ! initialize transformation matrices
        allocate( mlwf_transform( wf_nband, wf_nband, nkptnr))
        mlwf_transform = zzero
        do ix = 1, wf_nband
          mlwf_transform( ix, ix, :) = zone
        end do
        !mlwf_transform = wf_transform

        ! start minimization loop
        allocate( evec( wf_nband, wf_nband), eval( wf_nband))
        allocate( mlwf_r( wf_nband, wf_nband), &
                  mlwf_t( wf_nband, wf_nband), &
                  mlwf_dw( wf_nband, wf_nband), &
                  mlwf_dwcpy( wf_nband, wf_nband), &
                  auxmatcpy( wf_nband, wf_nband))
        allocate( mlwf_m( wf_nband, wf_nband, maxval( wf_n_n( 1:wf_n_usedshells)), wf_n_usedshells, nkptnr))
        allocate( ravg( 3, wf_nband))
        allocate( omegahist( minit), omegameanhist( minit), omegaod( wf_nband), omegad( wf_nband), omega2( wf_nband))
        allocate( integerlist( minit))
        do iz = 1, minit
          integerlist( iz) = iz
        end do
        mlwf_m = mlwf_m0
        iz = 0
        omegahist = 0.d0
        omegameanhist = 0.d0
        grad = 1.d0

        !write(*,*) 'Minimize localization functional...'
        do while( (iz .lt. minit) .or. ((iz .lt. 5000) .and. ((uncertainty .gt. input%properties%wannier%uncertainty) .or. grad .lt. 0.d0) .and. (omegastart .gt. omega)))
          iz = iz + 1
          !alpha = mixing/(1.d0+uncertainty+1.d-3*iz)
          !write(*,*) iz
          ! centers of Wannier functions
          ravg = 0.d0
          do ix = 1, wf_nband
            do is = 1, wf_n_usedshells
              do n = 1, wf_n_n( is)
                ravg( :, ix) = ravg( :, ix) - wf_n_wgt( is)/nkptnr*sum( real( aimag( log( mlwf_m( ix, ix, n, is, :)))))*wf_n_vc( :, n, is)
              end do
            end do
          end do
    
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, is, n, mlwf_r, mlwf_t, mlwf_dw, eval, evec, ix, auxmat, mlwf_dwcpy, auxmatcpy)
!$OMP DO
#endif
          do iknr = 1, nkptnr
            mlwf_dw(:,:) = zzero
            do is = 1, wf_n_usedshells
              do n = 1, wf_n_n( is)
                ! calculating R and T
                do ix = 1, wf_nband
                  mlwf_r( :, ix) = mlwf_m( :, ix, n, is, iknr)*conjg( mlwf_m( ix, ix, n, is, iknr))
                  mlwf_t( :, ix) = mlwf_m( :, ix, n, is, iknr)/mlwf_m( ix, ix, n, is, iknr)*(real( aimag( log( mlwf_m( ix, ix, n, is, iknr)))) + dot_product( wf_n_vc( :, n, is), ravg( :, ix)))
                end do
        
                ! calculating dW
                mlwf_r = mlwf_r - conjg( transpose( mlwf_r))
                mlwf_t = mlwf_t + conjg( transpose( mlwf_t))
                mlwf_dw(:,:) = mlwf_dw(:,:) + wf_n_wgt( is)*( 0.5*mlwf_r(:,:) + 0.5*zi*mlwf_t(:,:))
              end do
            end do
            mlwf_dw(:,:) = mlwf_dw(:,:)*alpha/dot_product( dble( wf_n_n( 1:wf_n_usedshells)), wf_n_wgt( 1:wf_n_usedshells))

            ! updating transformation matrices
            call diaghermat( wf_nband, zi*mlwf_dw, eval, evec)
            grad = norm2( eval)
            !if( uncertainty .gt. 5.d-3) then
            !  alpha = 1.d0
            !else
            !  alpha = 1.d0
            !endif
            do ix = 1, wf_nband
              mlwf_dw( :, ix) = exp( -zi*eval( ix))*evec( :, ix) 
            end do
            call ZGEMM( 'N', 'C', wf_nband, wf_nband, wf_nband, zone, &
                 mlwf_dw, wf_nband, &
                 evec, wf_nband, zzero, &
                 auxmat, wf_nband)
            call ZGEMM( 'N', 'N', wf_nband, wf_nband, wf_nband, zone, &
                 mlwf_transform( :, :, iknr), wf_nband, &
                 auxmat, wf_nband, zzero, &
                 mlwf_dw, wf_nband)
            mlwf_transform( :, :, iknr) = mlwf_dw(:,:)
          end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
    
          ! gauge-dependent part of localization functional
          omegaod = 0.d0
          omegad = 0.d0
          omega2 = 0.d0

          ! updating M
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, n, is, ix, auxmat) reduction(+:omegaod, omegad, omega2)
!$OMP DO  
#endif
          do iknr = 1, nkptnr
            do is = 1, wf_n_usedshells 
              do n = 1, wf_n_n( is)
                call ZGEMM( 'N', 'N', wf_nband, wf_nband, wf_nband, zone, &
                     mlwf_m0( :, :, n, is, iknr), wf_nband, &
                     mlwf_transform( :, :, wf_n_ik( n, is, iknr)), wf_nband, zzero, &
                     auxmat, wf_nband)
                call ZGEMM( 'C', 'N', wf_nband, wf_nband, wf_nband, zone, &
                     mlwf_transform( :, :, iknr), wf_nband, &
                     auxmat, wf_nband, zzero, &
                     mlwf_m( :, :, n, is, iknr), wf_nband)
                !omegaod = omegaod + bwgt( ib)*sum( abs( mlwf_m( :, :, iknr, ib))**2)
                do ix = 1, wf_nband
                  omegad( ix) = omegad( ix) + wf_n_wgt( is)*(real( aimag( log( mlwf_m( ix, ix, n, is, iknr)))) &
                                  + dot_product( wf_n_vc( :, n, is), ravg( :, ix)))**2
                  do iy = 1, wf_nband
                    if( iy .ne. ix) then
                      omegaod( ix) = omegaod( ix) + wf_n_wgt( is)*abs( mlwf_m( iy, ix, n, is, iknr))**2
                    end if
                  end do
                  omega2( ix) = omega2( ix) + wf_n_wgt( is)*(1.d0 - abs( mlwf_m( ix, ix, n, is, iknr))**2 + &
                      & ( real( aimag( log( mlwf_m( ix, ix, n, is, iknr))))&
                      &   + dot_product( wf_n_vc( :, n, is), ravg( :, ix)))**2 )
                end do
              end do
            end do
          end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
          ! convergence analysis
          omegaod = omegaod/nkptnr
          omegad = omegad/nkptnr
          omega2 = omega2/nkptnr

          omega = sum( omegai) + sum( omegad) + sum( omegaod)

          if( iz .eq. 1) then
            omegahist(:) = 0.d0
            omegameanhist(:) = 0.d0
            omegastart = omega
          end if
          omegahist = cshift( omegahist, -1)
          omegameanhist = cshift( omegameanhist, -1)
          omegahist(1) = omega
          omegamean = sum( omegahist(:))/min( minit, iz)
          omegameanhist(1) = omegamean
          if( iz .eq. 1) then
            uncertainty = 1.d0
            grad = 1.d0
          else
            uncertainty = sqrt( sum( (omegahist(:)-omegameanhist(:))**2)/(min( minit, iz)-1))/omega
            grad = dot_product( dble( integerlist( 1:min( minit, iz))), omegahist( 1:min( minit, iz))-omegamean) - (min( minit, iz)+1)*0.5d0*sum( omegahist( 1:min( minit, iz))-omegamean)
            grad = grad/sum( (dble( integerlist( 1:min( minit, iz)))-(min( minit, iz)+1)*0.5d0)**2)/omega
          end if
          !if( (iz .gt. minit) .and. (grad .lt. 1.d-5)) then
          !  z0 = iz
          !  if( iz .eq. z1+1) mixing = alpha
          !  alpha = mixing*exp( -0.8*2.d-3*(iz-z1))
          !else if( (iz .gt. minit) .and. (grad .lt. 1.d-3)) then
          !  z1 = iz
          !  if( iz .eq. z0+1) mixing = alpha
          !  alpha = mixing + 2*(0.8d0 - mixing)/pi*atan( 2.d-3*(iz-z0))
          !else
          !  z0 = iz
          !  z1 = iz
          !end if

          ! print out progress
          !write( 6, '(a,10x,"[",I4,"] omega = ",F13.6,5x,"uncertainty: ",F13.6)', advance='no') achar( 13), iz, omega, uncertainty 
          !flush( 6)
          !write(*,'(F13.6)') omega
          write(*,'(I4,4F23.16)') iz, omega, uncertainty, grad, alpha!, omega2 !omegai+omegad+omegaod
        end do

        write(*,*)
        if( omega .gt. omegastart) then
          write(*, '("ERROR (genmlwf): Localization functional diverged. Procedure aborted after ",I4," loops.")') iz
        else if( iz .ge. 5000) then
          write(*,*) 'ERROR (genmlwf): Not converged after 5000 cycles.'
        else
          write(*,'(" SUCCES: Convergence reached after ",I4," cycles.")') iz 
          write(*,'(" Localization gain: ",I3,"%")') nint( 100d0*(omegastart-omega)/omega)
        end if
        
        ! generating final transformation matrices for Hamiltonian eigenstates
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, auxmat)
!$OMP DO  
#endif
        do iknr = 1, nkptnr
          call ZGEMM( 'N', 'N', wf_nband, wf_nband, wf_nband, zone, &
               wf_transform( :, :, iknr), wf_nband, &
               mlwf_transform( :, :, iknr), wf_nband, zzero, &
               auxmat, wf_nband)
          wf_transform( :, :, iknr) = auxmat(:,:)
        end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
        call wannier_loc
        write(*,*) sum( wf_omega)
        write(*,'(100F13.6)') wf_omega
        write(*,'(100F13.6)') wf_centers( 1, :)
        write(*,'(100F13.6)') wf_centers( 2, :)
        write(*,'(100F13.6)') wf_centers( 3, :)
        ! write resulting transformation matrices to file
        call wannier_writefile( 'WANNIER')
#ifdef MPI
        call barrier
      else
        call barrier
      end if
#endif
  
      deallocate( auxmat, eval, evec, mlwf_m0, mlwf_m, mlwf_r, mlwf_t, mlwf_dw, mlwf_transform, omegahist, omegameanhist)
      return
      !EOC
    end subroutine wannier_gen_max
    !EOP
    
    subroutine wannier_gen_fromfile()
      logical :: succes
      real(8), allocatable :: ravg(:,:), omega(:)

      call wannier_init
      call wannier_readfile( 'WANNIER', succes)
      if( .not. succes) then
        call terminate
      end if

      !allocate( ravg( 3, wf_nband), omega( wf_nband))
      !call wannier_loc( omega, centers=ravg)
      !write(*,*) sum( omega)
      !write(*,'(100F13.6)') omega
      !write(*,'(100F13.6)') ravg( 1, :)
      !write(*,'(100F13.6)') ravg( 2, :)
      !write(*,'(100F13.6)') ravg( 3, :)
      return
    end subroutine wannier_gen_fromfile

    subroutine wannier_loc()
      integer :: iknr, is, n, i
      complex(8), allocatable :: loc_m(:,:,:,:,:), auxmat(:,:)
      
      allocate( auxmat( wf_nband, wf_nband))
      allocate( loc_m( wf_nband, wf_nband, maxval( wf_n_n( 1:wf_n_usedshells)), wf_n_usedshells, nkptnr))

      if( .not. wf_initialized) call wannier_init
      if( .not. allocated( wf_emat)) call wannier_emat
      if( .not. allocated( wf_centers)) allocate( wf_centers( 3, wf_nband))
      if( .not. allocated( wf_omega)) allocate( wf_omega( wf_nband))

#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, is, n, auxmat)
!$OMP DO
#endif
      do iknr = 1, nkptnr
        do is = 1, wf_n_usedshells 
          do n = 1, wf_n_n( is)
            call ZGEMM( 'N', 'N', wf_nband, wf_nband, wf_nband, zone, &
                 wf_emat( wf_bandstart:(wf_bandstart+wf_nband-1), wf_bandstart:(wf_bandstart+wf_nband-1), n, is, iknr), wf_nband, &
                 wf_transform( :, :, wf_n_ik( n, is, iknr)), wf_nband, zzero, &
                 auxmat, wf_nband)
            call ZGEMM( 'C', 'N', wf_nband, wf_nband, wf_nband, zone, &
                 wf_transform( :, :, iknr), wf_nband, &
                 auxmat, wf_nband, zzero, &
                 loc_m( :, :, n, is, iknr), wf_nband)
          end do
        end do
      end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
      wf_centers = 0.d0
      do i = 1, wf_nband
        do is = 1, wf_n_usedshells
          do n = 1, wf_n_n( is)
            wf_centers( :, i) = wf_centers( :, i) - wf_n_wgt( is)/nkptnr*sum( real( aimag( log( loc_m( i, i, n, is, :)))))*wf_n_vc( :, n, is)
          end do
        end do
      end do
      wf_omega = 0.d0
#ifdef USEOMP
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( iknr, is, n, i, auxmat) reduction(+:wf_omega)
!$OMP DO
#endif
      do iknr = 1, nkptnr
        do is = 1, wf_n_usedshells 
          do n = 1, wf_n_n( is)
            do i = 1, wf_nband
              wf_omega( i) = wf_omega( i) + wf_n_wgt( is)*(1.d0 - abs( loc_m( i, i, n, is, iknr))**2 + &
                  & ( real( aimag( log( loc_m( i, i, n, is, iknr))))&
                  &   + dot_product( wf_n_vc( :, n, is), wf_centers( :, i)))**2 )
            end do
          end do
        end do
      end do
#ifdef USEOMP
!$OMP END DO
!$OMP END PARALLEL
#endif
      wf_omega = wf_omega/nkptnr

      deallocate( loc_m, auxmat)
      return
    end subroutine wannier_loc

    ! print out projection local-orbitals
    subroutine wannier_showproj
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
    end subroutine wannier_showproj
    
    ! reads transformation matrices from file
    subroutine wannier_readfile( filename, succes)
      character(*), intent( in) :: filename
      logical, intent( out) :: succes

      ! local variables
      integer :: ik, ix, iy

      succes = .true.
      inquire( file=trim( filename)//trim( filext), exist=succes)
      if( .not. succes) then
        write(*,*) 'ERROR (wannier_readfile): File '//trim( filename)//trim( filext)//' does not exist.'
        return
      end if
      open( 50, file=trim( filename)//trim( filext), action='READ', form='UNFORMATTED', status='OLD')
      read( 50) wf_bandstart
      read( 50) wf_nband
      if( wf_nband .gt. wf_nprojtot) then
        write(*,*) 'WARNING (wannier_readfile): Content in file does not fit to content in species files. Check local-orbitals.'
      end if
      if( allocated( wf_projused)) deallocate( wf_projused)
      allocate( wf_projused( wf_nband))
      if( allocated( wf_projection)) deallocate( wf_projection)
      allocate( wf_centers( 3, wf_nband))
      if( allocated( wf_omega)) deallocate( wf_omega)
      allocate( wf_omega( wf_nband))
      if( allocated( wf_transform)) deallocate( wf_transform)
      allocate( wf_transform( wf_nband, wf_nband, nkptnr))
      do ix = 1, wf_nband
        read( 50) wf_projused( ix)
      end do
      do ix = 1, wf_nband
        read( 50) wf_omega( ix)
      end do
      do iy = 1, wf_nband
        do ix = 1, 3
          read( 50) wf_centers( ix, iy)
        end do
      end do
      outerloop: do ik = 1, nkptnr
        read( 50) ix
        if( .not.( ix .eq. ik)) then
          write(*,*) 'ERROR (wannier_readfile): File contains invalid or non fitting content.' 
          succes = .false.
          exit outerloop
        end if
        do iy = 1, wf_nband
          do ix = 1, wf_nband
            read( 50) wf_transform( ix, iy, ik)
          end do
        end do
      end do outerloop
      close( 50)
      if( succes) write(*,*) 'Transformation matrices succesfully read.'
      return
    end subroutine wannier_readfile
    
    ! writes transformation matrices to file
    subroutine wannier_writefile( filename)
      character(*), intent( in) :: filename
      ! local variables
      integer :: ik, ix, iy

      open( 50, file=trim( filename)//trim( filext), action='WRITE', form='UNFORMATTED')
      write( 50) wf_bandstart
      write( 50) wf_nband
      do ix = 1, wf_nband
        write( 50) wf_projused( ix)
      end do
      do ix = 1, wf_nband
        write( 50) wf_omega( ix)
      end do
      do iy = 1, wf_nband
        do ix = 1, 3
          write( 50) wf_centers( ix, iy)
        end do
      end do
      do ik = 1, nkptnr
        write( 50) ik
        do iy = 1, wf_nband
          do ix = 1, wf_nband
            write( 50) wf_transform( ix, iy, ik)
          end do
        end do
      end do
      close( 50)
      write( *, '(a,a)') ' Transformation matrices written to file ', trim( filename)//trim( filext)
      return
    end subroutine wannier_writefile  

    subroutine angchar( iknr, lmax, nmin, nmax)
      integer, intent( in) :: iknr, lmax, nmin, nmax

      integer :: n, l, m, lm, io, ia, is, ias, ngknr
      integer, allocatable :: igkignr(:)
      real(8), allocatable :: vgklnr(:,:,:), vgkcnr(:,:,:), gkcnr(:), tpgkcnr(:,:)
      complex(8), allocatable :: sfacgknr(:,:)
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

    !BOP
    ! !ROUTINE: wannier_geometry
    ! !INTERFACE:
    !
    subroutine wannier_geometry
      ! !USES:
      ! !INPUT/OUTPUT PARAMETERS:
      ! !DESCRIPTION:
      !   Sets geometry environment for Wannier functions. Finds near neighbors
      !   and corresponding geometric weights as well as ne number of shells
      !   needed for gradient calculation.
      !
      ! !REVISION HISTORY:
      !   Created September 2016 (Tillack)
      !EOP
      !BOC

      integer :: ix, iy, iz, ngridk(3), i, j, nshell, nkmax, iknr
      integer :: mt( nkptnr)
      real(8) :: vl(3), vc(3), dist
      real(8) :: nvlt( 3, nkptnr, nkptnr), nvct( 3, nkptnr, nkptnr), dt( nkptnr)
      real(8) :: coeff(6,6), coeffcpy(6,6), right(6), sval(6), lsvec(6,6), rsvec(6,6), coeffinv(6,6)
      logical :: stopshell

      integer :: lwork, info
      real(8), allocatable :: work(:)
      integer, allocatable :: iwork(:)
      
      ! find possible distances (shells)
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
      
      ! find all possible neighbors
      nshell = i
      nkmax = 0
      do iz = ngridk(3)-1, -ngridk(3)+1, -1
        do iy = ngridk(2)-1, -ngridk(2)+1, -1
          do ix = ngridk(1)-1, -ngridk(1)+1, -1
            vl = (/dble( ix)/ngridk(1), dble( iy)/ngridk(2), dble( iz)/ngridk(3)/)
            call r3mv( bvec, vl, vc)
            dist = norm2( vc)
            j = minloc( abs( dt(:) - dist), 1)
            if( abs( dt( j) - dist) .lt. input%structure%epslat) then
              mt( j) = mt( j) + 1
              nkmax = max( nkmax, mt(j))
              nvlt( :, mt( j), j) = vl
              nvct( :, mt( j), j) = vc
            end if
          end do
        end do
      end do
      
      ! allocating geometry arrays
      if( .not. allocated( wf_n_dist)) allocate( wf_n_dist( nshell))
      if( .not. allocated( wf_n_n)) allocate( wf_n_n( nshell))
      if( .not. allocated( wf_n_ik)) allocate( wf_n_ik( nkmax, nshell, nkptnr))
      if( .not. allocated( wf_n_vl)) allocate( wf_n_vl( 3, nkmax, nshell))
      if( .not. allocated( wf_n_vc)) allocate( wf_n_vc( 3, nkmax, nshell))
      if( .not. allocated( wf_n_wgt)) allocate( wf_n_wgt( nshell))

      ! sort everything by distance
      wf_n_n = 0
      wf_n_dist = 0.d0
      wf_n_vl = 0.d0
      wf_n_vc = 0.d0
      dist = minval( dt)
      do i = 1, nshell
        j = minloc( abs( dt( 1:nshell) - dist), 1)
        wf_n_dist( i) = dt( j)
        wf_n_n( i) = mt( j)
        wf_n_vl( :, :, i) = nvlt( :, :, j)
        wf_n_vc( :, :, i) = nvct( :, :, j)
        dt( j) = 1.d13 
        dist = minval( dt)
      end do
      
      ! find k-point indices of neighbors
      do iz = 0, ngridk( 3)-1
        do iy = 0, ngridk( 2)-1
          do ix = 0, ngridk( 1)-1
            iknr = modulo( iz, ngridk( 3))*ngridk( 2)*ngridk( 1) + &
                   modulo( iy, ngridk( 2))*ngridk( 1) + modulo( ix, ngridk(1))+1
            do j = 1, nshell
              do i = 1, wf_n_n( j)
                wf_n_ik( i, j, iknr) = modulo( iz + nint( wf_n_vl( 3, i, j)*ngridk(3)), ngridk( 3))*ngridk( 2)*ngridk( 1) + &
                                       modulo( iy + nint( wf_n_vl( 2, i, j)*ngridk(2)), ngridk( 2))*ngridk( 1) + &
                                       modulo( ix + nint( wf_n_vl( 1, i, j)*ngridk(1)), ngridk(1)) + 1
              end do
            end do
          end do
        end do
      end do
      
      ! find number of shells needed for gradient calculation and geometric weights
      j = 1
      stopshell = .false.
      coeff = 0.d0
      right(:) = (/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
      wf_n_wgt = 0.d0
      allocate( work( 1), iwork( 8*j))
      do while( (j .le. 6) .and. (.not. stopshell))
        do i = 1, wf_n_n( j)
          coeff( 1, j) = coeff( 1, j) + wf_n_vc( 1, i, j)**2
          coeff( 2, j) = coeff( 2, j) + wf_n_vc( 2, i, j)**2
          coeff( 3, j) = coeff( 3, j) + wf_n_vc( 3, i, j)**2
          coeff( 4, j) = coeff( 4, j) + wf_n_vc( 1, i, j)*wf_n_vc( 2, i, j)
          coeff( 5, j) = coeff( 5, j) + wf_n_vc( 1, i, j)*wf_n_vc( 3, i, j)
          coeff( 6, j) = coeff( 6, j) + wf_n_vc( 2, i, j)*wf_n_vc( 3, i, j)
        end do
        coeffcpy = coeff
        ! find pseudo-inverse of coeff
        call dgesdd( 'A', 6, j, coeffcpy( :, 1:j), 6, sval( 1:j), lsvec, 6, rsvec( 1:j, 1:j), j, work, -1, iwork, info)
        lwork = work(1)
        if( allocated( work)) deallocate( work)
        allocate( work( lwork))
        sval = 0.d0
        rsvec = 0.d0
        call dgesdd( 'A', 6, j, coeffcpy( :, 1:j), 6, sval( 1:j), lsvec, 6, rsvec( 1:j, 1:j), j, work, lwork, iwork, info)
        do i = 1, j
          if( sval( i) .gt. input%structure%epslat) then
            rsvec( i, :) = rsvec( i, 1:j)/sval( i)
          else
            rsvec( i, :) = 0.d0
          end if
        end do
        do i = j+1, 6
          rsvec( i, :) = 0.d0
        end do
        coeffinv( 1:j, :) = matmul( transpose( rsvec( :, 1:j)), transpose( lsvec))
        sval = matmul( matmul( coeff( :, 1:j), coeffinv( 1:j, :)), right)
        if( (sum( abs( sval - right)) .lt. input%structure%epslat) .or. (j .eq. 6)) then
          stopshell = .true.
          wf_n_wgt( 1:j) = matmul( coeffinv( 1:j, :), right)
          wf_n_usedshells = j
        end if
        j = j+1
      end do
      return
      !EOC
    end subroutine wannier_geometry
    !EOP

    !subroutine wfplotwannier( band)
    !  integer, intent( in) :: band
    !
    !  character(64) :: prefix
    !  character(64) :: filename

    !  integer :: iknr, iband, is, ia

    !  if( band .lt. 10) then
    !    write( prefix, '("WANNIER_N0",I1)') band
    !  else
    !    write( prefix, '("WANNIER_N",I2)') band
    !  end if

    !  input%properties%wfplot%plot3d%usesym=.false.
    !  labels=>create_plotlablels("Probability Density","WANNIER3D",3)
    !  call set_plotlabel_axis(labels,1,"a","lattice coordinate","graceunit")
    !  call set_plotlabel_axis(labels,2,"b","lattice coordinate","graceunit")
    !  call set_plotlabel_axis(labels,3,"c","lattice coordinate","graceunit")
    !  call set_plotlabel_axis(labels,4,"Wannier Function Norm Squared","","graceunit")
    !  Call plot3d(labels, 1, input%groundstate%lmaxapw, lmmaxapw, &
    !      & rhomt, rhoir, input%properties%wfplot%plot3d)
    !  call destroy_plotlablels(labels)
    !  if (rank==0) then
    !    Write(*,*)
    !    Write(*, '("Info(wfplotwannier):")')
    !    Write(*, '(" 3D Wannier function modulus squared written to WANNIER3D.xml")')
    !  end if
    !end subroutine wfplotwannier
end module mod_wannier
