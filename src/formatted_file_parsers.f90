!> Routines for reading objects needed for BSE calculation from file.
module formatted_file_parsers
  use precision, only: dp
  use constants, only: zone
  use modmpi, only: terminate_if_false
  use mod_hdf5
  use math_utils, only: boundary_mask
  use modinput, only: input_type

  implicit none
    
  private

  public :: read_eigen_energies, &
            read_grid_coordinates, &
            read_QP_energies

  contains


  !> Read eigen energies from file
  subroutine read_eigen_energies(file_name, E)
    !> Name of the file tp read from
    character(*), intent(in) :: file_name

    real(dp), intent(out), allocatable :: E(:, :)

    integer :: i, ik, is, N_states, n_k
    real(dp) :: line(3)
    character(256) :: charline

    open(unit=42, file=file_name)

    read(42,'(A)') charline
    read(charline(1:8),*) n_k

    read(42,'(A)') charline
    read(charline(1:8),*) N_states
            
    allocate(E(N_states, n_k))

    ! The first three lines are not needed.
    do i=1, 3
      read(42,*)
    end do
            
    do ik=1, n_k
      do is=1, N_states
        read(42,*) line
        E(is,ik) = line(2)
      end do
      if (ik < n_k) then  
        ! After each k-point section, four lines are not needed
        do i=1, 4
          read(42,*)
        end do
      end if
    end do
    
    close(42)
  end subroutine read_eigen_energies



  !> Read lattice coordinates of the grid points from exciting grid point file such as `KPOINTS.OUT` etc.
  subroutine read_grid_coordinates(file_name, coordinates) 
    !> File name
    character(*), intent(in) :: file_name
    !> Cartesian coordinates
    real(dp), intent(out), allocatable :: coordinates(:, :)

    integer :: i, N_points
    real(dp) :: line(9)
    character(256) :: charline

    open(unit=42, file=file_name)
    read(42, '(A)') charline
    read(charline(1 : 8), *) N_points

    allocate(coordinates(3, N_points))

    do i=1, N_points
      read(42, *) line
      coordinates(:, i) = line(5 : 7)
    end do

    close(42)
  end subroutine read_grid_coordinates


  !> Read QP energies from file. See the file header for the different content of the columns of the slices.
  function read_QP_energies(file_name, state_range, n_k) result(E_qp)
    !> File name
    character(*),  intent(in) :: file_name 
    !> Range of states, starting with the first, ending with the last.
    integer, intent(in) :: state_range(2)
    !> Number of \(\mathbf{k}\)-points 
    integer, intent(in) :: n_k

    real(dp), allocatable :: E_qp(:, :, :)

    real(dp) :: line(11)

    integer :: ik, is

    allocate(E_qp(state_range(1) : state_range(2), 11, n_k))

    open(42, file=file_name)

    read(42, *)
    read(42, *)

    do ik=1, n_k 
      do is=state_range(1), state_range(2) 
        read(42, *) line
        E_qp(is, :, ik) = line
      end do 

      read(42, *)
      read(42, *)
      read(42, *)
    end do
  end function read_QP_energies

end module formatted_file_parsers