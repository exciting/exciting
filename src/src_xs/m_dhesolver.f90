module m_dhesolver
  use modmpi
  use herm_eigensolver
  use modscl
  use mod_hdf5
  implicit none

  contains

    !BOP
    ! !ROUTINE: dhesolver
    ! !INTERFACE:
    subroutine dhesolver(ham, eval, binfo, evec, i1, i2, v1, v2, found, eecs)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    !   type(blacsinfo) :: binfo       ! Info type describing the BLACS grid
    !   integer(4), optional :: i1, i2 ! Index range of eigen-solutions
    !   real(8), optional :: v1, v2    ! Range of eigenvalues to search for
    !   integer(4), optional :: eecs   ! Estimate for eigenvalue clustering.
    !                                  ! Apriori not known but needed for propper
    !                                  ! orthogonalization of eigenvectors (default = 3)
    ! IN/OUT:
    !   type(dzmat) :: ham         ! 2D block cyclic distributed hermitian matrix
    !   real(8) :: eval(ham%nrows) ! Real valued eigenvalues in ascending order
    !   type(dzmat), optional :: evec   ! 2D block cyclic distributed eigenvector matrix
    ! OUT:
    !   integer(4), optional :: found ! How many solutions were found
    !
    ! !DESCRIPTION:
    !   Takes the upper triangular part of an distributed complex matrix matrix,
    !   assumed to be hermitian and finds eigenvalues and eigenvectors using
    !   the appropriate solver (ELPA, ScaLAPACK, or LAPACK).
    !
    ! !REVISION HISTORY:
    !   Created 2016 (Aurich)
    !EOP
    !BOC

      implicit none

      ! Arguments
      type(dzmat), intent(inout) :: ham
      real(8), intent(inout) :: eval(:)
      type(blacsinfo), intent(in) :: binfo
      ! Optional arguments
      type(dzmat), intent(inout), optional :: evec
      integer(4), intent(in), optional :: eecs
      integer(4), intent(in), optional :: i1, i2
      real(8), intent(in), optional :: v1, v2
      integer(4), intent(out), optional :: found

      call he_eigensolver_wrapper(ham, eval, binfo, evec, i1, i2, v1, v2, found, eecs)

    end subroutine dhesolver
    !EOC
end module m_dhesolver