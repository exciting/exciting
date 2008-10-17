
! Copyright (C) 2004-2008 S. Sagmeister and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine tetcalccw
  use modmain
  use modxs
  use modmpi
  use modtetra
  use m_genfilname
  implicit none
  ! local variables
  character(*), parameter :: thisnam='tetcalccw'
  integer :: iq
  logical :: tet
  call init0
  ! initialise universal variables
  tet=tetradf
  tetradf=.true.
  call init1
  ! save Gamma-point variables
  call xssave0
  ! initialize q-point set
  call init2
  ! read Fermi energy
  call readfermi
  ! w-point interval for process
  if (tscreen) then
     nwdf=1
     call genparidxran('q',nqpt)
  else
     call genparidxran('w',nwdf)
  end if
  ! loop over q-points
  do iq=qpari,qparf
     ! call for q-point
     call tetcalccwq(iq)
     write(unitout,'(a,i8)') 'Info('//thisnam//'): weights for tetrahedron &
          &method finished for q-point:',iq
  end do
  ! synchronize
  call barrier
  if ((procs.gt.1).and.(rank.eq.0).and.(.not.tscreen)) call tetgather
  call barrier
  tetradf=tet
  write(unitout,'(a)') "Info("//trim(thisnam)//"): weights for tetrahedron &
       &method finished"
  call genfilname(setfilext=.true.)
end subroutine tetcalccw
