module m_dhesolver
  use modmpi
  use modinput, only: input
  use herm_eigensolver
  use modscl
  use mod_hdf5
  implicit none

  contains

    !BOP
    ! !ROUTINE: dhesolver
    ! !INTERFACE:
    subroutine dhesolver(ham, eval, binfo, evec, i1, i2, v1, v2, found)
    ! !INPUT/OUTPUT PARAMETERS:
    ! IN:
    !   type(blacsinfo) :: binfo       ! Info type describing the BLACS grid
    !   integer(4), optional :: i1, i2 ! Index range of eigen-solutions
    !   real(8), optional :: v1, v2    ! Range of eigenvalues to search for
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
      !> ham:Hermitian matrix to diagonalize
      type(dzmat), intent(inout) :: ham
      !> eval:Eigenvalues (output)
      real(8), intent(inout) :: eval(:)
      !> binfo:BLACS context information
      type(blacsinfo), intent(in) :: binfo
      ! Optional arguments
      !> evec:Eigenvectors (optional)
      type(dzmat), intent(inout), optional :: evec
      !> i1,i2:Index bounds for eigenvalue subset (optional)
      integer(4), intent(in), optional :: i1, i2
      !> v1,v2:Value bounds for eigenvalue subset (optional)
      real(8), intent(in), optional :: v1, v2
      !> found:Number of eigenvalues found (optional)
      integer(4), intent(out), optional :: found

      logical :: distributed, sane

      if (present(evec)) then
        distributed = ham%isdistributed .and. evec%isdistributed
        sane = (ham%isdistributed .eqv. evec%isdistributed)
      else
        distributed = ham%isdistributed
        sane = .true.
      end if

      call terminate_if_false(sane, 'Error(dhesolver): Inconsistent matrix distribution')

      select case( input%xs%bse%bsesolver )
        case( 'elpa1StageSolver' )
          call terminate_if_false(check_lower_index_and_energy_selection(i1, v1, v2), &
                  &'Error(dhesolver): Requested subset selection not available for ELPA 1stage solver.&
                  & Use a different solver or compute all eigenvalues and select yourself.')
          call elpa_eigensolver(ham, eval, binfo, '1', evec, i2, found)
        case( 'elpa2StageSolver' )
          call terminate_if_false(check_lower_index_and_energy_selection(i1, v1, v2), &
                  &'Error(dhesolver): Requested subset selection not available for ELPA 2stage solver.&
                  & Use a different solver or compute all eigenvalues and select yourself.')
          call elpa_eigensolver(ham, eval, binfo, '2', evec, i2, found)
        case( 'scalapackPzheevx' )
          call scalapack_eigensolver_pzheevx(ham, eval, binfo, evec, i1, i2, v1, v2, found, input%xs%bse%eecs)
        case( 'scalapackPzheevd' )
          call terminate_if_false(check_lower_index_and_energy_selection(i1, v1, v2), &
                  &'Error(dhesolver): Requested subset selection not available for ScaLAPACK pzheevd solver.&
                  & Use a different solver or compute all eigenvalues and select yourself.')
          call terminate_if_false(.not. check_upper_index_subset_selection(ham, i2), &
                  &'Error(dhesolver): Requested upper index selection not available for ScaLAPACK pzheevd solver.&
                  & Use a different solver or compute all eigenvalues and select yourself.')
          call terminate_if_false(present(evec), 'Error(dhesolver): The ScaLAPACK pzheevd solver needs eigenvectors.')
          call scalapack_eigensolver_pzheevd(ham, eval, binfo, evec, found)
        case( 'lapack' )
          call lapack_eigensolver(ham, eval, binfo, evec, i1, i2, v1, v2, found)
        case default
          call terminate('Error(dhesolver): Unknown bsesolver. This is most likely an implementation error.')
      end select

    end subroutine dhesolver

    !> Check if lower index selection is 1 and no energy selection is given
    function check_lower_index_and_energy_selection(i1, v1, v2) result(is_valid_subset_selection)
      !> i1:Lower index bound for eigenvalue subset (optional)
      integer(4), intent(in), optional :: i1
      !> v1,v2:Value bounds for eigenvalue subset (optional)
      real(8), intent(in), optional :: v1, v2
      !> is_valid_subset_selection:check if lower index selection is 1 and no energy selection is given
      logical :: is_valid_subset_selection

      is_valid_subset_selection = .true.
      if (present(i1)) then
        if (i1 /= 1) then
          is_valid_subset_selection = .false.
        end if
      end if
      if (present(v1) .or. present(v2)) then
        is_valid_subset_selection = .false.
      end if

    end function check_lower_index_and_energy_selection

    !> Check if upper index selection is requesting a real subset
    function check_upper_index_subset_selection(ham, i2) result(is_upper_index_subset_selection)
      !> ham:Hermitian matrix to diagonalize
      type(dzmat), intent(in) :: ham
      !> i2:Upper index bound for eigenvalue subset (optional)
      integer(4), intent(in), optional :: i2
      !> is_upper_index_subset_selection:check if upper index selection is requesting a real subset
      logical :: is_upper_index_subset_selection

      is_upper_index_subset_selection = .false.
      if (present(i2)) then
        if (i2 /= ham%nrows) then
          is_upper_index_subset_selection = .true.
        end if
      end if

    end function check_upper_index_subset_selection
    !EOC
end module m_dhesolver