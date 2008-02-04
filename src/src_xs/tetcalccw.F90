
! Copyright (C) 2004-2007 S. Sagmeister and C. Ambrosch-Draxl.
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
  if (calledxs.eq.1) call init0
  ! initialise universal variables
  tet=tetra
  tetra=.true.
  call init1
  ! save Gamma-point variables
  call tdsave0
  ! initialize q-point set
  call init2xs
  ! read Fermi energy
  call readfermi
  ! w-point interval for process
  call genparidxran('w')
  ! loop over q-points
  do iq=1,nqpt
     ! call for q-point
     call tetcalccwq(iq)
     write(unitout,'(a,i8)') 'Info('//thisnam//'): weights for tetrahedron &
          &method finished for q-point:',iq
  end do
  ! synchronize
  call barrier
  if ((procs.gt.1).and.(rank.eq.0)) call tetgather
  call barrier
  tetra=tet
  write(unitout,'(a)') "Info("//trim(thisnam)//"): weights for tetrahedron &
       &method finished"
  call genfilname(setfilext=.true.)
end subroutine tetcalccw
