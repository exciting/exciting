!
!
!
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine writedveff (iq, is, ia, ip, dveffmt, dveffir)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: iq
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ip
      Complex (8), Intent (In) :: dveffmt (lmmaxvr, nrcmtmax, natmtot0)
      Complex (8), Intent (In) :: dveffir (ngrtot0)
! local variables
      Integer :: js
      Character (256) :: fext
      Character (256) :: chdummy
      Call phfext (iq, is, ia, ip, 0, 1, chdummy, fext, chdummy)
      call write_potential_response( 'DVEFF'//trim(fext), nspecies, natoms(1:nspecies), nrcmt(1:nspecies), input%groundstate%lmaxvr, ngrid, &
        dveffmt, shape( dveffmt ), dveffir, shape( dveffir ) )
End Subroutine

subroutine write_potential_response( path, nspecies, natoms, nrmt, lmax, ngrid, dveffmt, shapemt, dveffir, shapeir )
  use precision, only: dp
#include "asserts.fpp"
  use modmpi, only: terminate_if_false
  use os_utils, only: path_exists
  use mod_misc, only: version
  character(*), intent(in) :: path
  integer, intent(in) :: nspecies
  integer, intent(in) :: natoms(nspecies)
  integer, intent(in) :: nrmt(nspecies)
  integer, intent(in) :: lmax
  integer, intent(in) :: ngrid(3)
  integer, intent(in) :: shapemt(3), shapeir(1)
  complex(dp), intent(in) :: dveffmt(shapemt(1), shapemt(2), shapemt(3))
  complex(dp), intent(in) :: dveffir(shapeir(1))

  integer :: nrmtmax, natmtot, lmmax, ngrtot
  integer :: un, ierr, is

  nrmtmax = maxval( nrmt(1:nspecies) )
  natmtot = sum( natoms(1:nspecies) )
  lmmax = (lmax + 1)**2
  ngrtot = product( ngrid )

  CALL_ASSERT( shapemt(1) == lmmax,  'Inconsistent 1st dimension of argument `dveffmt`.' )
  CALL_ASSERT( shapemt(2) == nrmtmax,  'Inconsistent 2nd dimension of argument `dveffmt`.' )
  CALL_ASSERT( shapemt(3) == natmtot,  'Inconsistent 3rd dimension of argument `dveffmt`.' )
  CALL_ASSERT( shapeir(1) == ngrtot,  'Inconsistent size of argument `dveffir`.' )
  
  open( newunit=un, file=trim( path ), action='write', form='unformatted', iostat=ierr )

  write( un ) version
  write( un ) nspecies
  write( un ) lmmax
  do is = 1, nspecies
    write( un ) natoms(is)
    write( un ) nrmt(is)
  end do
  write( un ) ngrid
  write( un ) dveffmt, dveffir

  close( un )
end subroutine write_potential_response
