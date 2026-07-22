!
!
!
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!
Subroutine readdveff (iq, is, ia, ip, dveffmt, dveffir)
      Use modmain
      Use modinput
      Implicit None
! arguments
      Integer, Intent (In) :: iq
      Integer, Intent (In) :: is
      Integer, Intent (In) :: ia
      Integer, Intent (In) :: ip
      Complex (8), Intent (Out) :: dveffmt (lmmaxvr, nrcmtmax, natmtot)
      Complex (8), Intent (Out) :: dveffir (ngrtot)
! local variables
      Character (256) :: fext
      Character (256) :: chdummy
      Call phfext (iq, is, ia, ip, 0, 1, chdummy, fext, chdummy)
      call read_potential_response( 'DVEFF'//trim(fext), nspecies, natoms(1:nspecies), nrcmt(1:nspecies), input%groundstate%lmaxvr, ngrid, &
        dveffmt, shape( dveffmt ), dveffir, shape( dveffir ) )
End Subroutine

subroutine read_potential_response( path, nspecies, natoms, nrmt, lmax, ngrid, dveffmt, shapemt, dveffir, shapeir )
  use, intrinsic :: iso_fortran_env, only: output_unit
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
  complex(dp), intent(out) :: dveffmt(shapemt(1), shapemt(2), shapemt(3))
  complex(dp), intent(out) :: dveffir(shapeir(1))

  integer :: nrmtmax, natmtot, lmmax, ngrtot
  integer :: un, ierr, is
  integer :: version_f(3), nspecies_f, lmmax_f, natoms_f, nrmt_f, ngrid_f(3)
  character(1024) :: errmsg

  nrmtmax = maxval( nrmt(1:nspecies) )
  natmtot = sum( natoms(1:nspecies) )
  lmmax = (lmax + 1)**2
  ngrtot = product( ngrid )

  CALL_ASSERT( shapemt(1) == lmmax,  'Inconsistent 1st dimension of argument `dveffmt`.' )
  CALL_ASSERT( shapemt(2) == nrmtmax,  'Inconsistent 2nd dimension of argument `dveffmt`.' )
  CALL_ASSERT( shapemt(3) == natmtot,  'Inconsistent 3rd dimension of argument `dveffmt`.' )
  CALL_ASSERT( shapeir(1) == ngrtot,  'Inconsistent size of argument `dveffir`.' )
  
  call terminate_if_false( path_exists( path, ierr ), '(read_potential_response) &
    Requested file path does not exist.' )

  open( newunit=un, file=trim( path ), action='read', form='unformatted', status='old', iostat=ierr )

  read( un ) version_f
  if( any( version_f /= version ) ) then
    call warning( 'Warning(read_potential_response):' )
    call warning( ' Different versions' )
    write( errmsg, '(" current	 : ", i3.3, ".", i3.3, ".", i3.3)' ) version
    call warning( errmsg )
    write( errmsg, '(" file   	 : ", i3.3, ".", i3.3, ".", i3.3)' ) version_f
    call warning( errmsg )
  end if

  read( un ) nspecies_f
  write( errmsg, '("Different values for variable `",a,"` in file ",a," and input argument.",a," &
    input: ",i6,a,"&
    file : ",i6)' ) 'nspecies', trim( path ), new_line( 'a' ), nspecies, new_line( 'a' ), nspecies_f
  call terminate_if_false( nspecies_f == nspecies, '(read_potential_response)'//new_line( 'a' )//trim( errmsg ) )

  read( un ) lmmax_f
  write( errmsg, '("Different values for variable `",a,"` in file ",a," and input argument.",a," &
    input: ",i6,a,"&
    file : ",i6)' ) 'lmmax', trim( path ), new_line( 'a' ), lmmax, new_line( 'a' ), lmmax_f
  call terminate_if_false( lmmax_f == lmmax, '(read_potential_response)'//new_line( 'a' )//trim( errmsg ) )

  do is = 1, nspecies
    read( un ) natoms_f
    write( errmsg, '("Different values for variable `",a,"` in file ",a," and input argument.",a," &
      input: ",i6,a,"&
      file : ",i6)' ) 'natoms', trim( path ), new_line( 'a' ), natoms(is), new_line( 'a' ), natoms_f
    call terminate_if_false( natoms_f == natoms(is), '(read_potential_response)'//new_line( 'a' )//trim( errmsg ) )

    read( un ) nrmt_f
    write( errmsg, '("Different values for variable `",a,"` in file ",a," and input argument.",a," &
      input: ",i6,a,"&
      file : ",i6)' ) 'nrmt', trim( path ), new_line( 'a' ), nrmt(is), new_line( 'a' ), nrmt_f
    call terminate_if_false( nrmt_f == nrmt(is), '(read_potential_response)'//new_line( 'a' )//trim( errmsg ) )
  end do

  read( un ) ngrid_f
  write( errmsg, '("Different values for variable `",a,"` in file ",a," and input argument.",a," &
    input: ",3i6,a,"&
    file : ",3i6)' ) 'ngrid', trim( path ), new_line( 'a' ), ngrid, new_line( 'a' ), ngrid_f
  call terminate_if_false( all( ngrid_f == ngrid ), '(read_potential_response)'//new_line( 'a' )//trim( errmsg ) )

  read( un ) dveffmt, dveffir

  close( un )
end subroutine read_potential_response
