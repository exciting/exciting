!> This module contains routines that detect and filter ghost states.
!> Ghost states can occur, for example, in calculations with semicore or
!> Dirac-type local orbitals.
module ghost_band_filter

   use precision, only: dp
   use modmpi, only: terminate

   implicit none

   private
   public :: detect_ghost_bands, filter_ghost_bands, reshape_sv_arrays, restore_shape_sv_arrays, &
            &set_default_ghost_band_parameters
 
 contains
   !> Return a large artificial eigenvalue that is guaranteed to lie above the
   !> current spectrum, including the case where all eigenvalues are negative.
   function get_shifted_ghost_state_eval(evals) result(ghost_state_eval)

      !> eigenvalues
      real(dp), intent(in) :: evals(:, :)

      real(dp) :: ghost_state_eval
      real(dp) :: largest_eval
      real(dp) :: spectrum_scale

      largest_eval = maxval(evals)
      spectrum_scale = max(1.0_dp, maxval(abs(evals)))

      ghost_state_eval = largest_eval + 10.0_dp * spectrum_scale

   end function get_shifted_ghost_state_eval

   !> This subroutine identifies ghost states.
   !> A state is classified as a ghost state if its eigenenergy is lower than the lowest linearization energy 
   !> by more than a tolerance specified in the input file. 
   !> A warning is printed if ghost states are detected. 
   subroutine detect_ghost_bands(lorbe0, apwe0, n_states, evals, tolerance_smallest_allowed_eval, n_ghost_states, emit_warnings)
        
      !> local orbital energies
      real(dp), intent(in) :: lorbe0(:, :, :)
      !> APW energies
      real(dp), intent(in) :: apwe0(:, :, :)
      !> number of states
      integer, intent(in) :: n_states
      !> eigenvalues 
      real(dp), intent(in) :: evals(:, :)
      !> allowed tolerance for eigenvalues to be smaller than the smallest linearization energy
      real(dp), intent(in) :: tolerance_smallest_allowed_eval
      !> number of ghost states
      integer, intent(out) :: n_ghost_states
      !> if false, warning message is suppresed
      logical, intent(in), optional :: emit_warnings

      ! local variables
      real(dp) :: smallest_allowed_eval
      ! state index
      integer :: ist
      ! error message
      character(len=200) :: msg
      character(len=20) :: tol_str
      logical :: do_warn

      do_warn = .true.
      if (present(emit_warnings)) do_warn = emit_warnings
      
      ! detect ghost bands
      smallest_allowed_eval = merge(minval(lorbe0), minval(apwe0), minval(lorbe0) < minval(apwe0))
      
      n_ghost_states = 0
      ist = n_ghost_states

      do while (any(evals(ist+1, :) < smallest_allowed_eval - tolerance_smallest_allowed_eval))
         ist = ist + 1
         if (ist > n_states-1) then
            call terminate('Error(detect_ghost_bands): Something has gone terribly wrong with the KS spectrum.&
            & All eigenvalues are smaller than the smallest allowed eigenvalue.')
         endif
      enddo

      ! Print warning if ghost states were detected
      if (ist > 0) then
         n_ghost_states = ist

         if (do_warn) then
            write(tol_str, '(F6.3)') tolerance_smallest_allowed_eval
            msg = 'Bands with eigenvalues more than ' // trim(tol_str) // ' &
                 Ha below the lowest linearization energy were detected.'
            call warning('Warning(detect_ghost_bands):')
            call warning(msg)
            call warning('This can indicate the existence of ghost bands.')
         end if 
      end if
       
   end subroutine detect_ghost_bands

   !> This subroutine processes ghost states.
   !> Ghost states are assigned a large artificial eigenvalue above the current
   !> spectrum, and the eigenvalue and eigenvector arrays are then reordered
   !> accordingly.
   !> Note that this can lead to errors, e.g. when summing over all eigenvalues.
   subroutine filter_ghost_bands(n_states, n_ghost_states, evals, evecs, emit_warnings)

      !> number of states
      integer, intent(in) :: n_states
      !> number of ghost states
      integer, intent(in) :: n_ghost_states
      !> eigenvalues 
      real(dp), intent(inout) :: evals(:, :)
      !> eigenvectors 
      complex(dp), intent(inout) :: evecs(:, :, :)
      !> if false, warning message is suppresed     
      logical, intent(in), optional :: emit_warnings

      ! local variables
      ! state index
      integer :: ist
      ! temporary evec array for reordering
      complex(dp), allocatable :: evecs_temp(:, :, :)
      ! shape of evecs
      integer :: shape_evecs(3)
      ! large artificial eigenvalue for ghost states
      real(dp) :: ghost_state_eval
      ! warning message
      character(len=200) :: msg
      ! first band index occupied by shifted ghost states
      integer :: first_shifted_ghost_band
      logical :: do_warn

      do_warn = .true.
      if (present(emit_warnings)) do_warn = emit_warnings


      ghost_state_eval = get_shifted_ghost_state_eval(evals)
      first_shifted_ghost_band = n_states - n_ghost_states + 1

      if (do_warn) then
         call warning('Warning(ghost_band_filter): Ghost band filter is activated.')
         call warning('To keep the original matrix sizes, the ghost states are not completely removed,')
         call warning('but are assigned a large artificial eigenvalue above the current spectrum.')
         call warning('This can lead to problems, for example when summing over all eigenvalues.')
         write(msg, '("Shifted ghost bands start at band index ", I0, ".")') first_shifted_ghost_band
         call warning(trim(msg))
      end if

      ! save eigenvectors of ghost states
      shape_evecs = shape(evecs)
      allocate( evecs_temp(shape_evecs(1), n_ghost_states, shape_evecs(3)) )
      evecs_temp(:, 1:n_ghost_states, :) = evecs(:, 1:n_ghost_states, :)

      ! shift non-ghost states to the beginning of the arrays
      do ist = n_ghost_states+1, n_states
         evals(ist-n_ghost_states, :) = evals(ist, :)
         evecs(:, ist-n_ghost_states, :) = evecs(:, ist, :)
      enddo

      ! append ghost states to the end of the arrays
      evals(n_states-n_ghost_states+1:n_states, :) = ghost_state_eval
      evecs(:, n_states-n_ghost_states+1:n_states, :) = evecs_temp(:, 1:n_ghost_states, :)

      deallocate(evecs_temp)

   end subroutine filter_ghost_bands

   !> In order to avoid code duplication, we reshape the arrays of the second-variational eigenvalues and eigenstates
   !> to fit the shape of the first-variational ones. This allows us to use only one `filter_ghost_bands` routine.
   subroutine reshape_sv_arrays(evalsv, evecsv, n_states, ik, evalsv_reshape, evecsv_reshape)
      !> second variational eigenvalues 
      real(dp), intent(in) :: evalsv(:, :)
      !> second variational eigenvectors 
      complex(dp), intent(in) :: evecsv(:, :)
      !> number of states
      integer, intent(in) :: n_states
      !> k-point index
      integer, intent(in) :: ik
      !> reshaped second variational eigenvalues 
      real(dp), intent(out) :: evalsv_reshape(:, :)
      !> reshaped second variational eigenvectors 
      complex(dp), intent(out) :: evecsv_reshape(:, :, :)

      ! local variables 
      integer :: ist

      evalsv_reshape(1:n_states, 1) = evalsv(1:n_states, ik)
      do ist=1,n_states
        evecsv_reshape(1, ist, 1:n_states) = evecsv(1:n_states, ist)
      end do

   end subroutine reshape_sv_arrays

   !> Restores the original shape of the second-variational eigenvalue and eigenvector arrays.
   subroutine restore_shape_sv_arrays(evalsv_reshape, evecsv_reshape, n_states, ik, evalsv, evecsv)
      !> reshaped second variational eigenvalues 
      real(dp), intent(in) :: evalsv_reshape(:, :)
      !> reshaped second variational eigenvectors 
      complex(dp), intent(in) :: evecsv_reshape(:, :, :)
      !> reshaped number of states
      integer, intent(in) :: n_states
      !> k-point index
      integer, intent(in) :: ik
      !> second variational eigenvalues 
      real(dp), intent(inout) :: evalsv(:, :)
      !> second variational eigenvectors 
      complex(dp), intent(out) :: evecsv(:, :)

      ! local variables 
      integer :: ist

      evalsv(1:n_states, ik) = evalsv_reshape(1:n_states, 1)
      do ist=1,n_states
        evecsv(1:n_states, ist) = evecsv_reshape(1, ist, 1:n_states)
      end do

   end subroutine restore_shape_sv_arrays

   !> Ensures GhostBands structure is initialized with default values
   !> if not provided in the input file.
   subroutine set_default_ghost_band_parameters
      use modinput
      use inputdom
      
      if ( .not. (associated(input%groundstate%GhostBands))) then
         ! set the default values if element not present
         input%groundstate%GhostBands => getstructGhostBands (emptynode)
      end if
   end subroutine set_default_ghost_band_parameters

end module ghost_band_filter
