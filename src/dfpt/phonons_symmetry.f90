!> This module contains variable types describing and procedures to find
!> special symmetry adapted displacement patterns for phonon calculations,
!> so called *irreducible representations*.
!> 
! See 'phonons_symmetry.md' for further documentation and explanations.
!>{!../src/dfpt/phonons_symmetry.md!}
module phonons_symmetry
  use precision, only: dp
  use modmpi, only: terminate_if_false

  implicit none
  private

  ! DERIVED DATA TYPES
  !> variable type describing an irreducible representation \(I\)
  type, public :: irrep
    !> dimension of the irrep \(d_I\) (number of members)
    integer :: dim = 0
    !> displacement patterns \(p^{I \mu}_{\alpha i}({\bf q})\)
    complex(dp), allocatable :: pat(:,:,:)
    !> matrix representation \(\texttt{S}^{I}({\bf q})\) of each symmetry operation
    complex(dp), allocatable :: symmat(:,:,:)
    contains
      !> free memory
      procedure :: free => free_irrep
  end type

  !> variable type describing the symmetry adapted basis set for
  !> the expansion of the phonons with given \({\bf q}\)
  type, public :: irrep_basis
    !> number of symmetries in the small group of \({\bf q}\)
    integer :: nsym = 0
    !> global indices of symmetries in small group of \({\bf q}\)
    integer, allocatable :: isym(:)
    !> lattice vectors that map \({\bf \rm S}{\bf q}\) back to 1st BZ
    integer, allocatable :: ivsym(:,:)
    !> number of irreps
    integer :: nirrep = 0
    !> set of irreps
    type(irrep), allocatable :: irreps(:)
    contains
      !> free memory
      procedure :: free => free_irrep_basis
  end type

  public :: ph_sym_find_irreps

  contains

    !> This subroutine finds the irreps for a given wavevector \({\bf q}\).
    !>
    !> This is done by first finding the small group of \({\bf q}\), \(\mathcal{G}_{\bf q}\).
    !> Then, a random dynamical matrix is set up and symmetrized with the symmetries
    !> of \(\mathcal{G}_{\bf q}\). The eigenvectors of this random symmetrized dynamical
    !> matrix describe the displacement patterns \(p^{I \mu}_{\alpha i}({\bf q})\) and
    !> all eigenvectors corresponding to degenerate eigenvalues form the \(d_I\) members
    !> of an irrep \(I\). In the end, for each symmetry operation from \(\mathcal{G}_{\bf q}\),
    !> its matrix representation \({\bf \texttt{S}}^I({\bf q})\) in the basis of the irreps
    !> is found.
    !>
    !> When canonical displacement patterns are assumed, \(3N_{\rm at}\) one-dimensional
    !> irreps are returned and the respective displacement patterns are set to
    !> \(p^{I \mu=1}_{\alpha i}({\bf q}) = 1\) for \(I = 3(\alpha - 1) + i\) and 0 otherwise.
    !> The only symmetry considered in this case is the identity.
    !>
    !> When the unperturbed system is assumed, \(N_{\rm at}\) three-dimensional 
    !> irreps are returned and the respective displacement patterns are set to
    !> \(p^{I \mu}_{\alpha i}({\bf q}) = \delta_{i \mu}\, \delta_{\alpha I}\).
    !> The number of symmetries is set to the total number of crystal symmetries in
    !> \(\mathcal{G}\) and the matrix representations are set to the \(3 \times 3\) 
    !> rotation matrices \(\texttt{S}^I_{\mu \nu}({\bf 0}) = {\rm S}_{\mu \nu}\).
    subroutine ph_sym_find_irreps( vql, basis, &
        canonical, unperturbed)
      use phonons_util, only: ph_util_symmetrize_dyn, ph_util_diag_dynmat, ph_util_symmetry_T
      use math_utils, only: get_degeneracies, all_zero
      use constants, only: zzero, zone
      use mod_lattice, only: bvec
      use mod_symmetry, only: nsymcrys, symlat, lsplsymc, ieqatom, symlatc
      use mod_atoms, only: natmtot, nspecies, natoms, idxas
      use modinput
      !> wavevector \({\bf q}\) in lattice coordinates
      real(dp), intent(in) :: vql(3)
      !> basis of displacement patterns
      type(irrep_basis), intent(out) :: basis
      !> use canonical displacement patterns (default: `.false.`)
      logical, optional, intent(in) :: canonical
      !> use symmetries of unperturbed system (default: `.false.`), 
      !> has to be used in combination with `vql=[0.0,0.0,0.0]`
      logical, optional, intent(in) :: unperturbed

      real(dp), parameter :: eps = 1e-8_dp     ! epsilon for degenerate eigenvalues

      logical :: canon, unpert
      integer :: i, j, k, is, ia, ias, ja, jas
      real(dp) :: v(3), arand(3,3), brand(3,3)
      character(256) :: errmsg

      integer, allocatable :: deg(:,:)
      real(dp), allocatable :: rauxdyn(:,:,:)
      complex(dp), allocatable :: auxdyn(:,:), dynsym(:,:), symmat(:,:)

      ! use canonical displacement patterns
      canon = .false.
      if( present( canonical ) ) canon = canonical
      ! use symmetries of unperturbed system
      unpert = .false.
      if( present( unperturbed) ) unpert = (unperturbed .and. (sum( abs( vql ) ) < 1e-12_dp))

      ! allocate arrays
      call basis%free
      allocate( basis%isym(nsymcrys) )
      allocate( basis%ivsym(3, nsymcrys) )

      ! find small group of q
      call findgroupq( .false., vql, input%structure%epslat, bvec, symlat, nsymcrys, lsplsymc, &
             basis%nsym, basis%isym, basis%ivsym )

      ! find lattice vectors G with S^T q = q + G
      do i = 1, basis%nsym
        call r3mtv( dble( symlat(:, :, lsplsymc(basis%isym(i))) ), vql, v )
        basis%ivsym(:, i) = nint( v - vql )
        call terminate_if_false( sum( abs( v - vql - basis%ivsym(:, i) ) ) < input%structure%epslat, '(ph_sym_find_irreps) &
          Error in finding lattice vectors for small group of q.' )
      end do
      
      ! generate random dynamical matrix
      ! We want the irreps to be identical in each run, on each machine and for all compilers.
      ! Therefore, we need the same random numbers everytime (see lcgrand) and we need
      ! to fix the gauge (phase) of the irreps (eigenvectors).
      allocate( rauxdyn(3*natmtot, 3*natmtot, 2) )
      allocate( auxdyn(3*natmtot, 3*natmtot), source=zzero )
      allocate( dynsym(3*natmtot, 3*natmtot), source=zzero )
      allocate( symmat(3*natmtot, 3*natmtot), source=zzero )
      call lcgrand( rauxdyn, size( rauxdyn ), 12 )
      do is = 1, nspecies
        do ia = 1, natoms(is)
          ias = idxas(ia, is)
          i = (ias - 1) * 3
          inner: do ja = 1, natoms(is)
            jas = idxas(ja, is)
            j = (jas - 1) * 3
            arand = rauxdyn(i+1:i+3, j+1:j+3, 1)
            brand = rauxdyn(i+1:i+3, j+1:j+3, 2)
            if( norm2( vql ) < input%structure%epslat ) brand = 0.0_dp
            do k = 1, basis%nsym
              if( ja == ieqatom(ia, is, basis%isym(k)) ) then
                dynsym(i+1:i+3, j+1:j+3) = cmplx( arand,  brand, dp )
                cycle inner
              end if
            end do
          end do inner
        end do
      end do
      dynsym = dynsym + conjg( transpose( dynsym ) )

      call ph_util_symmetrize_dyn( vql, dynsym, basis%nsym, basis%isym )

      ! diagonalize dynamical matrix
      if( canon ) then
        ! unity matrix for eigenvectors (= irreps) for canonical displacements
        basis%nsym = 1
        basis%isym(1) = 1
        basis%ivsym = 0
        do i = 1, 3*natmtot
          rauxdyn(i, 1, 1) = dble( i )
          auxdyn(i, i) = zone
        end do
      else
        call ph_util_diag_dynmat( dynsym, rauxdyn(:, 1, 1), auxdyn )
      end if

      ! find irreps
      if( unpert ) then
        basis%nirrep = natmtot
        allocate( basis%irreps(basis%nirrep) )
        do i = 1, basis%nirrep
          basis%irreps(i)%dim = 3
          allocate( basis%irreps(i)%pat(3, natmtot, 3), source=zzero )
          do j = 1, 3
            basis%irreps(i)%pat(j, :, j) = zone
          end do
          allocate( basis%irreps(i)%symmat(basis%irreps(i)%dim, basis%irreps(i)%dim, basis%nsym) )
        end do
      else
        deg = get_degeneracies( rauxdyn(:, 1, 1), eps )
        basis%nirrep = size( deg, dim=2 )
        allocate( basis%irreps(basis%nirrep) )
        do k = 1, basis%nirrep
          basis%irreps(k)%dim = deg(3, k)
          i = deg(1, k); j = deg(2, k)
          allocate( basis%irreps(k)%pat(3, natmtot, basis%irreps(k)%dim), &
                    source=reshape( auxdyn(:, i:j), [3, natmtot, basis%irreps(k)%dim] ) )
          allocate( basis%irreps(k)%symmat(basis%irreps(k)%dim, basis%irreps(k)%dim, basis%nsym) )
        end do
      end if

      ! find matrix representation of each symmetry in irrep basis
      do i = 1, basis%nsym
        symmat = ph_util_symmetry_T( basis%isym(i), vql )
        call zgemm( 'c', 'n', 3*natmtot, 3*natmtot, 3*natmtot, zone, &
               auxdyn, 3*natmtot, &
               symmat, 3*natmtot, zzero, &
               dynsym, 3*natmtot )
        call zgemm( 'n', 'n', 3*natmtot, 3*natmtot, 3*natmtot, zone, &
               dynsym, 3*natmtot, &
               auxdyn, 3*natmtot, zzero, &
               symmat, 3*natmtot )
        k = 0
        do j = 1, basis%nirrep
          if( unpert) then
            k = lsplsymc(basis%isym(i))
            basis%irreps(j)%symmat(:, :, i) = cmplx( symlatc(:, :, k), 0, dp )
          else
            basis%irreps(j)%symmat(:, :, i) = symmat(k+1:k+basis%irreps(j)%dim, k+1:k+basis%irreps(j)%dim)
            k = k + basis%irreps(j)%dim
            ! check if symmetry matrix is unitary
            call zgemm( 'c', 'n', basis%irreps(j)%dim, basis%irreps(j)%dim, basis%irreps(j)%dim, zone, &
                   basis%irreps(j)%symmat(1, 1, i), basis%irreps(j)%dim, &
                   basis%irreps(j)%symmat(1, 1, i), basis%irreps(j)%dim, zzero, &
                   dynsym, 3*natmtot )
            do ja = 1, basis%irreps(j)%dim
              dynsym(ja, ja) = dynsym(ja, ja) - zone
            end do
            write( errmsg, '("(ph_sym_find_irreps) Non-unitary symmetry matrix for q=",3f13.6,", irrep ",i3," and symmetry ",i2,".")' ) vql, j, i
            call terminate_if_false( all_zero( dynsym(1:basis%irreps(j)%dim, 1:basis%irreps(j)%dim), tol=eps ), trim( errmsg ) )
          end if
        end do
      end do

      deallocate( rauxdyn, auxdyn, dynsym, symmat )
    end subroutine ph_sym_find_irreps

    !> destroy irrep
    subroutine free_irrep( this )
      class(irrep), intent(inout) :: this
      if( allocated( this%pat ) ) deallocate( this%pat )
      if( allocated( this%symmat ) ) deallocate( this%symmat )
      this%dim = 0
    end subroutine free_irrep

    !> destroy irrep basis
    subroutine free_irrep_basis( this )
      class(irrep_basis), intent(inout) :: this
      integer :: i
      if( allocated( this%isym ) ) deallocate( this%isym )
      if( allocated( this%ivsym ) ) deallocate( this%ivsym )
      if( allocated( this%irreps ) ) then
        do i = 1, this%nirrep
          call this%irreps(i)%free
        end do
        deallocate( this%irreps )
      end if
      this%nsym = 0
      this%nirrep = 0
    end subroutine free_irrep_basis

end module phonons_symmetry
