!> Evaluate the square root of the Matrix with the bare Coulomb potential, whose
!> eigenvalues and eigenvectors are stored as `barcev` and `vmat`
subroutine setbarcev(eig_tol, remove_g_equal_zero)
#include "asserts.fpp"
    use constants, only: zzero
    use mod_coulomb_potential, only: barc, barcev, vmat
    use mod_product_basis, only: mbsiz, matsiz
    use modinput, only: input
    use modgw, only: fdebug
    use modmpi, only: rank
    use precision, only: dp, i32
    use xlapack, only: matrix_multiply
#include "offload.fpp"

    implicit none

    !> Eigenvalues smaller than `eig_tol` are discarded
    real(dp), intent(in) :: eig_tol
    !> If `.true.`, remove the eigenvector corresponding to a plane-wave with \(G=0\)
    logical, intent(in) :: remove_g_equal_zero

    integer(i32) :: i, j, imax, n_eigs, dim
    logical, allocatable :: keep(:) ! indicate which barc eigenvectors are kept
    complex(dp), allocatable :: wi0(:)
    complex(dp), allocatable :: wi0new(:)

!!REVISION HISTORY:
! Created July 31,2009 by Hong Jiang
! Readjusted Jan, 2012 by DIN
! Reformulated Jun 2024 by Ronaldo

    keep = ( barcev >= eig_tol )
    dim = size( vmat, 1 )
    n_eigs = size( vmat, 2 )
    CALL_ASSERT( dim == matsiz, 'dim must be equal to the global variable matsiz' )
    CALL_ASSERT( n_eigs == size( barcev ), 'barcev must have n_eigs elements' )
    
    if (remove_g_equal_zero) then
      allocate( wi0(dim), source=zzero )
      call calcwmix0(wi0)

      allocate( wi0new(n_eigs), source=zzero )
      call matrix_multiply(vmat, wi0, wi0new, trans_A='C')

      ! find the index of the diagonalized barc eigenvector that has maximal
      ! overlap with G=0 (constant) plane wave
      imax = maxloc( abs(wi0new), dim=1 )
      if (input%gw%debug .and. rank==0) then
        write(fdebug,*)'- Maximum singular eigenvector ###'
        write(fdebug,'("immax, max(wi0new), barcev(immax): ",i4,4x,f12.6,4x,f12.6)') imax, abs(wi0new(imax))**2, barcev(imax)
      end if
      ! exclude this matrix element
      keep(imax) = .false.
    end if

    if (input%gw%debug) then
      if (mbsiz < matsiz) then
        if (rank==0) then
          write(fdebug,*) "Info(setbarcev): Product basis size has been changed"
          write(fdebug,*) " - Old basis set size =", matsiz
          write(fdebug,*) " - New basis set size =", mbsiz
        end if
      end if
    end if

    ! Build the trasformation matrix
    mbsiz = count( keep )
    if (allocated(barc)) then
        OMP_OFFLOAD target exit data map(delete: barc)
        deallocate(barc)
    end if
    allocate( barc(matsiz, mbsiz), source=zzero )

    i = 0
    do j = 1, n_eigs
      if ( keep(j) ) then
        i = i + 1
        barc(:, i) = vmat(:, j)*sqrt( cmplx( barcev(j), kind=dp ) )
      end if
    end do

    OMP_OFFLOAD target enter data map(always, to: barc)

end subroutine
!EOC
