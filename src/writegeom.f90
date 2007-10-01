
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writegeom
! !INTERFACE:
subroutine writegeom(topt)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   topt : .true. if GEOMETRY-OPT.OUT is to be written
! !DESCRIPTION:
!   Outputs the lattice vectors and atomic positions to file, in a format
!   which may be then used directly in {\tt exciting.in}. If {\tt topt} is
!   {\tt .false.} then the file name is {\tt GEOMETRY.OUT}, otherwise it is
!   {\tt GEOMETRY\_OPT.OUT}.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: topt
! local variables
integer is,ia
real(8) v(3)
if (topt) then
  open(50,file='GEOMETRY_OPT'//trim(filext),action='WRITE',form='FORMATTED')
else
  open(50,file='GEOMETRY'//trim(filext),action='WRITE',form='FORMATTED')
end if
write(50,*)
write(50,'("scale")')
write(50,'(" 1.0")')
write(50,*)
write(50,'("scale1")')
write(50,'(" 1.0")')
write(50,*)
write(50,'("scale2")')
write(50,'(" 1.0")')
write(50,*)
write(50,'("scale3")')
write(50,'(" 1.0")')
write(50,*)
write(50,'("avec")')
write(50,'(3G18.10)') avec(:,1)
write(50,'(3G18.10)') avec(:,2)
write(50,'(3G18.10)') avec(:,3)
if (molecule) then
  write(50,*)
  write(50,'("molecule")')
  write(50,'(" ",L1)') molecule
end if
write(50,*)
write(50,'("atoms")')
write(50,'(I4,T40," : nspecies")') nspecies
do is=1,nspecies
  write(50,'(" ''",A,"''",T40," : spfname")') trim(spfname(is))
  write(50,'(I4,T40," : natoms; atpos, bfcmt below")') natoms(is)
  do ia=1,natoms(is)
    v(:)=atposl(:,ia,is)
! use Cartesian coordinates for the molecular case
    if (molecule) call r3mv(avec,v,v)
    write(50,'(3F14.8,"  ",3F12.8)') v(:),bfcmt(:,ia,is)
  end do
end do
close(50)
return
end subroutine
!EOC

