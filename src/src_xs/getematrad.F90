
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getematrad(iq)
  use modmain
  use modxs
  use m_genfilname
  use m_getunit
  implicit none
  ! arguments
  integer, intent(in) :: iq
  ! local variables
  integer :: lmax1,lmax2,lmax3,un
  character(256) :: fname
  lmax1=max(lmaxapwwf,lolmax)
  lmax2=lmaxemat
  ! lmax1 and lmax3 should be the same!
  lmax3=lmax1
  if (allocated(riaa)) deallocate(riaa)
  allocate(riaa(0:lmax1,apwordmax,0:lmax3,apwordmax,0:lmax2,natmtot, &
       ngq(iq)))
  if (allocated(riloa)) deallocate(riloa)
  allocate(riloa(nlomax,0:lmax3,apwordmax,0:lmax2,natmtot,ngq(iq)))
  if (allocated(rilolo)) deallocate(rilolo)
  allocate(rilolo(nlomax,nlomax,0:lmax2,natmtot,ngq(iq)))
  call genfilname(basename='EMATRAD',iq=iq,filnam=fname)
  call getunit(un)
  open(un,file=trim(fname),form='unformatted',action='read', &
       status='old')
  read(un) riaa,riloa,rilolo
  close(un)
end subroutine getematrad
