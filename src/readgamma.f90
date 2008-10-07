
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine readgamma(gq)
use modmain
implicit none
! arguments
real(8), intent(out) :: gq(3*natmtot,nqpt)
! local variables
integer iq,i
integer natmtot_,nqpt_,iq_,i_
real(8) vql_(3),vqc_(3)
! external functions
real(8) r3dist
external r3dist
open(50,file='GAMMAQ.OUT',action='READ',form='FORMATTED',status='OLD')
read(50,*)
read(50,*) natmtot_
if (natmtot.ne.natmtot_) then
  write(*,*)
  write(*,'("Error(readgamma): differing natmtot")')
  write(*,'(" current    : ",I4)') natmtot
  write(*,'(" GAMMAQ.OUT : ",I4)') natmtot_
  write(*,*)
  stop
end if
read(50,*) nqpt_
if (nqpt.ne.nqpt_) then
  write(*,*)
  write(*,'("Error(readgamma): differing nqpt")')
  write(*,'(" current    : ",I6)') nqpt
  write(*,'(" GAMMAQ.OUT : ",I6)') nqpt_
  write(*,*)
  stop
end if
read(50,*)
do iq=1,nqpt
  read(50,*) iq_
  if (iq.ne.iq_) then
    write(*,*)
    write(*,'("Error(readgamma): incorrect q-point index in GAMMAQ.OUT for&
     & q-point ",I6)') iq
    write(*,*)
    stop
  end if
  read(50,*) vql_
  if (r3dist(vql(:,iq),vql_).gt.epslat) then
    write(*,*)
    write(*,'("Error(readgamma): differing q-vectors in lattice coordinates for&
     & q-point ",I6)') iq
    write(*,'(" current    : ",3G18.10)') vql(:,iq)
    write(*,'(" GAMMAQ.OUT : ",3G18.10)') vql_
    write(*,*)
    stop
  end if
  read(50,*) vqc_
  if (r3dist(vqc(:,iq),vqc_).gt.epslat) then
    write(*,*)
    write(*,'("Error(readgamma): differing q-vectors in Cartesian coordinates&
     & for q-point ",I6)') iq
    write(*,'(" current    : ",3G18.10)') vqc(:,iq)
    write(*,'(" GAMMAQ.OUT : ",3G18.10)') vqc_
    write(*,*)
    stop
  end if
  do i=1,3*natmtot
    read(50,*) i_,gq(i,iq)
    if (i.ne.i_) then
      write(*,*)
      write(*,'("Error(readgamma): incorrect mode index in GAMMAQ.OUT for&
       & q-point ",I6)') iq
      write(*,*)
      stop
    end if
  end do
  read(50,*)
end do
close(50)
return
stop
end subroutine

