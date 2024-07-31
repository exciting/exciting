module quasiparticle_energies
  use asserts, only: assert
  use m_getunit, only: getunit
  use precision, only: dp, i32

  implicit none

  private

  character(len=*), parameter :: file_name_qp_energies = 'EVALQP'
  character(len=*), parameter :: extension_text_format ='.DAT'

  public :: write_qp_energies_text_format

contains

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
  call assert( size( kpt_coordinates_lattice, 1 ) == 3, 'kpt_coordinates_lattice must have size = 3 along 1st dim.' )
  call assert( size( kpt_indexes ) == size( kpt_coordinates_lattice, 2 ), 'kpt_indexes and kpt_coordinates_lattice have incompatible sizes' )

  
  call getunit(fid)
  open( fid, file=file_name_qp_energies//extension_text_format, action='WRITE', form='FORMATTED' )

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