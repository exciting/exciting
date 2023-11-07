!> Module for setting up the transition energies from independent particle (IP) eigen energies.
module bse_diagonal
  use precision, only: dp
  use asserts, only: assert
  use modmpi, only: mpiinfo, distribute_loop, terminate_if_false
  use modinput, only: input_type
  use exciting_mpi, only: xmpi_bcast
  use multi_index_conversion, only: composite_index_to_indices
  use os_utils, only: path_exists
  use grid_utils, only: partial_grid

  implicit none

  private
  public :: setup_transition_energies

  

  contains

  !> Set up diagonal part of the Bethe Salpeter Hamiltonian. The diagonal part
  !> is defined as all possible differences between the unoccupied and occupied
  !> eigen energies of the independent particle problem at one k-point:
  !> \[
  !>   D_{i_o i_u \mathbf{k}} = \epsilon_{i_u \mathbf{k}} - \epsilon_{i_o \mathbf{k}},
  !> \]
  !> where \( i_o \) and \( i_u \) are the indices of the unoccupied bands respectively
  !> and \( \mathbf{k} \) is the corresponding k-point.
  subroutine setup_transition_energies(occupied_bands, unoccupied_bands, transition_energies)
    !> Eigen energies of the occupied bands
    real(dp), intent(in) :: occupied_bands(:, :)
    !> Eigen energies of the unoccupied bands
    real(dp), intent(in) :: unoccupied_bands(:, :)
    !> Transition energies \(D\).
    real(dp), allocatable, intent(out) :: transition_energies(:, :, :)

    integer :: n_o, n_u, n_k, i_o, i_u, i_k

    call assert(size(occupied_bands, dim=2) == size(unoccupied_bands, dim=2), &
            'Occupied and unoccupied eigen energies have not the same number of k-points.')


    n_o = size(occupied_bands, dim=1)
    n_u = size(unoccupied_bands, dim=1)
    n_k = size(occupied_bands, dim=2)

    allocate(transition_energies(n_u, n_o, n_k))
    do i_k=1, n_k
      do i_o=1, n_o
        do i_u=1, n_u
          transition_energies(i_u, i_o, i_k) = unoccupied_bands(i_u, i_k) - occupied_bands(i_o, i_k)
        end do
      end do
    end do
  end subroutine setup_transition_energies

end module bse_diagonal