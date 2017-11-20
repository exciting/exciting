! Copyright (C) 2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!
!BOP
! !ROUTINE: getematrad
! !INTERFACE:
Subroutine getematrad(iqr, iq)
! !USES:
  use modinput, only: input
  use modxs, only: riaa, riloa, rilolo, ngq, ematraddir
  use mod_APW_LO, only: lolmax, apwordmax, nlomax
  use mod_atoms, only: natmtot
  use m_genfilname
  use m_getunit
! !INPUT/OUTPUT PARAMETERS:
! IN:
!   integer(4), iqr : Index of non-reduced q-point
!   integer(4), iq  : Index of reduced q-point
!
! !DESCRIPTION:
!   The routine deallocates, reallocates and
!   reads the radial integral arryas form file. 
!   {\tt EMATRAD} is used as file basename.
!
! !REVISION HISTORY:
!   Added to documentation scheme. (Aurich 2016)
!   Changed formatting. (Aurich 2016)
!EOP
!BOC

  implicit none

  ! arguments
  integer(4), intent(in) :: iqr, iq

  ! local variables
  integer(4) :: lmax1, lmax2, lmax3, un
  character(256) :: fname

  lmax1 = max(input%xs%lmaxapwwf, lolmax)
  lmax2 = input%xs%lmaxemat

  ! lmax1 and lmax3 should be the same!
  lmax3 = lmax1

  if(allocated(riaa)) deallocate(riaa)
  allocate(riaa(0:lmax1, apwordmax, 0:lmax3, apwordmax, 0:lmax2, natmtot, ngq(iq)))

  if(allocated(riloa)) deallocate(riloa)
  allocate(riloa(nlomax, 0:lmax3, apwordmax, 0:lmax2, natmtot, ngq(iq)))

  if(allocated(rilolo)) deallocate(rilolo)
  allocate(rilolo(nlomax, nlomax, 0:lmax2, natmtot, ngq(iq)))

  call genfilname(basename=trim(adjustl(ematraddir))//'/'//'EMATRAD',&
    & iq=iqr, appfilext=.true., filnam=fname)

  call getunit(un)
  open(un, file=trim(fname), form='unformatted', action='read', status='old')
  read(un) riaa, riloa, rilolo
  close(un)
end subroutine getematrad
!EOC
