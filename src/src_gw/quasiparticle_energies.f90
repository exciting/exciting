module quasiparticle_energies
#include "asserts.fpp"
  use mod_large_io, only: inquire_large, open_direct_unformatted_large
  use precision, only: dp, i32, long_int

  implicit none

  private

  character(len=*), parameter :: file_name_qp_energies = 'EVALQP'
  character(len=*), parameter :: extension_text_format ='.DAT'

  public :: checkpoint_matches_expected_qp_window
  public :: write_qp_energies_text_format

contains

!> Return whether an evGW0 checkpoint matches the expected k-point set and QP-band window.
logical function checkpoint_matches_expected_qp_window(file_name, kpoint_vectors, first_qp_band, last_qp_band) result(flag)
  !> Checkpoint file to inspect.
  character(len=*), intent(in) :: file_name
  !> Expected irreducible k-point vectors in lattice coordinates.
  real(dp), intent(in) :: kpoint_vectors(:, :)
  !> Expected first QP band.
  integer(i32), intent(in) :: first_qp_band
  !> Expected last QP band.
  integer(i32), intent(in) :: last_qp_band

  real(dp), allocatable :: eqp_in_file(:), eks_in_file(:)
  real(dp) :: efqp_in_file, efks_in_file, kpoint_vector_in_file(3)
  logical :: exists
  integer(i32) :: ib_in_file, ik, io_status, nb_in_file, nkpt_in_file, unit_number
  integer(long_int) :: record_length

  CALL_ASSERT(size(kpoint_vectors, 1) == 3, 'kpoint_vectors must have size 3 along the first dimension')

  flag = .false.
  inquire(file=trim(file_name), exist=exists)
  if (exists) then
    allocate(eqp_in_file(first_qp_band:last_qp_band), eks_in_file(first_qp_band:last_qp_band))
    call inquire_large(record_length, [nkpt_in_file, ib_in_file, nb_in_file], kpoint_vector_in_file, &
      eqp_in_file, eks_in_file, [efqp_in_file, efks_in_file])
    call open_direct_unformatted_large(unit_number, trim(file_name), "read", record_length, "old")
    flag = .true.
    ik = 1_i32
    do while (flag .and. ik <= size(kpoint_vectors, 2))
      read(unit_number, rec=ik, iostat=io_status) nkpt_in_file, ib_in_file, nb_in_file, &
        kpoint_vector_in_file, eqp_in_file, eks_in_file, efqp_in_file, efks_in_file

      if (io_status == 0) then
        flag = nkpt_in_file == size(kpoint_vectors, 2) &
          .and. ib_in_file == first_qp_band &
          .and. nb_in_file == last_qp_band &
          .and. maxval(abs(kpoint_vector_in_file - kpoint_vectors(:, ik))) <= 1.0e-6_dp
      else
        flag = .false.
      end if

      ik = ik + 1_i32
    end do

    close(unit_number)
    deallocate(eqp_in_file, eks_in_file)
  end if
end function checkpoint_matches_expected_qp_window

!> Write the QP energies to an output file with text format
subroutine write_qp_energies_text_format( kpt_indexes, kpt_coordinates_lattice, &
    kpt_weights, first_band, KS_eigenvalues, QP_eigenvalues, VXC_diag_elements, &
    sigma_x, sigma_c, Znk )
  

  !> Indexes of the k-points to be printed out
  integer(i32), intent(in) :: kpt_indexes(:)
  !> k-points coordinates in terms of the lattice vectors
  real(dp), intent(in) :: kpt_coordinates_lattice(:, :)
  !> k-points weights
  real(dp), intent(in) :: kpt_weights(:)
  !> Index of the 1st band
  integer(i32), intent(in) :: first_band
  !> KS eigenvalues (1st index: band, 2nd index: k-point)
  real(dp), intent(in) :: KS_eigenvalues(first_band:, :)
  !> QP eigenvalues (1st index: band, 2nd index: k-point)
  real(dp), intent(in) :: QP_eigenvalues(first_band:, :)
  !> Diagonal elements of the VXC matrix (1st index: band, 2nd index: k-point)
  real(dp), intent(in) :: VXC_diag_elements(first_band:, :)
  !> Exchange part of the self-energy
  complex(dp), intent(in) :: sigma_x(first_band:, :)
  !> Correlation part of the self-energy
  complex(dp), intent(in) :: sigma_c(first_band:, :)
  !> Renormalization factor
  real(dp), intent(in) :: Znk(first_band:, :)
  

  integer(i32) :: ie, ikp, fid, last_band
  real(dp) :: de, dx
  real(dp) :: ehf, eks, egw
  real(dp) :: vxc, sx, scr, sci, z

  last_band = ubound( KS_eigenvalues, 1 )
  CALL_ASSERT( size( kpt_coordinates_lattice, 1 ) == 3, 'kpt_coordinates_lattice must have size = 3 along 1st dim.' )
  CALL_ASSERT( size( kpt_indexes ) == size( kpt_coordinates_lattice, 2 ), 'kpt_indexes and kpt_coordinates_lattice have incompatible sizes' )

  
  open( newunit=fid, file=file_name_qp_energies//extension_text_format, action='WRITE', form='FORMATTED' )

  do ikp = 1, size( kpt_indexes )
    write( fid,'("k-point #",I6,":",4F12.6)' ) kpt_indexes(ikp), kpt_coordinates_lattice(:, ikp), kpt_weights(ikp)
    write( fid, '(A6,10(A12,5X))' ) 'state', 'E_KS[Ha]', 'E_HF[Ha]', 'E_GW[Ha]', 'Sx[Ha]', 'Re(Sc)[Ha]', 'Im(Sc)[Ha]', 'Vxc[Ha]', 'DE_HF[Ha]', 'DE_GW[Ha]', 'Znk'
    do ie = first_band, last_band
      eks = KS_eigenvalues(ie, ikp)
      egw = QP_eigenvalues(ie, ikp)
      de = egw-eks
      vxc = VXC_diag_elements(ie, ikp)
      sx = real( sigma_x(ie, ikp) )
      scr = real( sigma_c(ie, ikp) )
      sci = aimag( sigma_c(ie, ikp) )
      z = Znk(ie, ikp)
      dx = sx-vxc
      ehf = eks+dx
      write( fid, '(I4,2X,10(F16.8,1X))' ) ie, eks, ehf, egw, sx, scr, sci, vxc, dx, de, z
    end do ! ie
    write( fid, * )
  end do ! ikp

  close( fid )
  
  end subroutine
  

end module
